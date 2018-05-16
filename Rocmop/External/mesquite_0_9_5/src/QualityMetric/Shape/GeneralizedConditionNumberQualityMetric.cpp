/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   GeneralizedConditionNumberQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include "GeneralizedConditionNumberQualityMetric.hpp"
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

#include <math.h>
#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
   using namespace std;
#endif



using namespace Mesquite;

GeneralizedConditionNumberQualityMetric::GeneralizedConditionNumberQualityMetric()
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  set_name("Generalized Condition Number");
}
/*
#undef __FUNC__
#define __FUNC__ "GeneralizedConditionNumberQualityMetric::evaluate_element"

double GeneralizedConditionNumberQualityMetric::evaluate_element(PatchData *pd, int element_index,
                                                      MsqError &err)
{
  // THIS FUNCTION SHOULD EVENTUALLY BE MADE VERY EFFICIENT BY USING
  // RAW ARRAYS DIRECTLY.  RIGHT NOW IT IS A HACK TO PORT OPT-MS

  if ( pd->get_storage_mode() != PatchData::RAW_ARRAYS ) {
    err.set_msg("Need raw arrays in evaluate element function taking pd as argument\n");
    return(0.0);
  }
  if (element_index >= pd->get_element_array_size()) {
    err.set_msg("element index exceeds element array size\n");
    return(0.0);
  }

  Vector3D* coords = pd->get_coords_array(err);
  ConnectivityArrayT element_connectivity = (pd->get_connectivity_array(err))[element_index];

  switch (element_connectivity.entity_type) {
  case TRIANGLE:
    MsqNode node1(coords[element_connectivity.indices[0]][0], 
                  coords[element_connectivity.indices[0]][1], 
                  coords[element_connectivity.indices[0]][2]);
    MsqNode node2(coords[element_connectivity.indices[1]][0], 
                  coords[element_connectivity.indices[1]][1], 
                  coords[element_connectivity.indices[1]][2]);
    MsqNode node3(coords[element_connectivity.indices[2]][0], 
                  coords[element_connectivity.indices[2]][1], 
                  coords[element_connectivity.indices[2]][2]);
    MsqTri tri(&node1,&node2,&node3);
    return( this->evaluate_element(&tri, err) );
    break;
    //  default:
    //    err.set_msg("only supporting triangles in evaluate element for now\n");
    //    return(0.0);
  }
}
*/

bool GeneralizedConditionNumberQualityMetric::evaluate_element(PatchData &pd,
                                                      MsqMeshEntity *element,
                                                               double &fval,
                                                      MsqError &err)
{
  int num_sample_points;
  bool return_flag;
  vector<Vector3D> sample_points;
  ElementEvaluationMode eval_mode = get_element_evaluation_mode();
  element->get_sample_points(eval_mode,sample_points,err);  MSQ_ERRZERO(err);
  vector<Vector3D>::iterator iter=sample_points.begin();
    // loop over sample points
  Vector3D jacobian_vectors[3];
  short num_jacobian_vectors;
  int i=0;
  num_sample_points=sample_points.size();
  std::vector<double> metric_values(num_sample_points);
    //Vector3D* current_sample_point;
  for(i=0;i<num_sample_points;++i){
      // compute weighted jacobian
    element->compute_weighted_jacobian(pd, (*iter),
                                       jacobian_vectors,
                                       num_jacobian_vectors, err); MSQ_ERRZERO(err);
      // evaluate condition number at ith sample point
      //if 2 jacobian vectors (2D elem)
    
    return_flag=compute_condition_number(pd, element, jacobian_vectors,
                                         num_jacobian_vectors,
                                         metric_values[i],err);
    if(MSQ_CHKERR(err) || !return_flag){
      return false;
    }
    
    ++iter;
  }// end loop over sample points
  fval=average_metrics(&metric_values[0],num_sample_points,err); MSQ_ERRZERO(err);
  return true;
}

bool GeneralizedConditionNumberQualityMetric::evaluate_vertex(PatchData &/*pd*/,
                                                              MsqVertex* /*vertex*/,
                                                              double &fval,
                                                              MsqError &err)
{
  MSQ_SETERR(err)("Condition Number's evaluate_vertex is currently "
                  "being implemented", MsqError::NOT_IMPLEMENTED);
  fval=0.0;
  return false;
  
    /*
  fval=0.0;
  size_t this_vert = pd.get_vertex_index(vert);
  vector<size_t> adj_elems;
  pd.get_vertex_element_indices(this_vert, adg_elems, err);
  double num_elems = adj_elems.size();
  double *metric_values=new double[ num_elems ];
  int num_jacobian_vectors;
  Vector3D sample_point;
  Vecotr3D jacobian_vectors[3];
  
  int i;
  for ( i = 0; i<num_elems; i++){
    elems[i]->get_sample_point(vertex, sample_point, err); MSQ_CHKERR(err);
    elems[i]->compute_weighted_jacobian(&current_sample_point,
                                        jacobian_vectors,
                                        num_jacobian_vectors, err);
                                        MSQ_CHKERR(err);
    metric_values[i]=compute_condition_number(jacobian_vectors,
                                              num_jacobian_vectors, err);
                                              MSQ_CHKERR(err);
  }
  total_metric=average_metrics(metric_values, num_elems, err);
  MSQ_CHKERR(err);
  delete metric_values;
    */
}
  
