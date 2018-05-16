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
  \file   ASMQualityMetric.cpp
  \brief  QualityMetric class for ASM (Area Smoothness Metric)

  \author Michael Brewer
  \date   2002-06-9
*/
#include "ASMQualityMetric.hpp"
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
   using std::vector;
#endif

using namespace Mesquite;


ASMQualityMetric::ASMQualityMetric()
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err);
  MSQ_CHKERR(err);
  avgMethod=QualityMetric::MAXIMUM;
  feasible=0;
  set_name("Area Smoothness");
}

bool ASMQualityMetric::evaluate_element(PatchData &pd,
                                                    MsqMeshEntity *element,
                                                    double &fval,
                                                    MsqError &err)
{
  double temp_double;
  size_t elem_ind=pd.get_element_index(element);
  vector<size_t> adj_elems;
 
  MsqMeshEntity *elems = pd.get_element_array(err);
  switch(element->get_element_type()){
    case TRIANGLE:
    case QUADRILATERAL:
      pd.get_adjacent_entities_via_n_dim(1,elem_ind,adj_elems,err);  MSQ_ERRZERO(err);
      break;
    case TETRAHEDRON:
    case HEXAHEDRON:
      pd.get_adjacent_entities_via_n_dim(2,elem_ind,adj_elems,err);  MSQ_ERRZERO(err);
      break;
    default:
      MSQ_SETERR(err)("ASM quality metric not implemented for this "
                      "element type.", MsqError::NOT_IMPLEMENTED);
      return false;
  };
  int num_samp=adj_elems.size();
  if(num_samp < 1){
    fval=0.0;
  }
  else{
    double* met_vals = new double [num_samp];
    int i=0;
    switch(element->get_element_type()){
      case TRIANGLE:
      case QUADRILATERAL:
        temp_double=element->compute_unsigned_area(pd,err);  MSQ_ERRZERO(err);
          //PRINT_INFO("\nunsigned area = %f",temp_double);
        for(i=0;i<num_samp;++i){
          met_vals[i]=elems[adj_elems[i]].compute_unsigned_area(pd,err);  MSQ_ERRZERO(err);
            //PRINT_INFO("neighboring nunsigned area = %f",met_vals[i]);
          if((temp_double+met_vals[i])>MSQ_MIN){
            met_vals[i]=fabs((temp_double-met_vals[i])/(temp_double+
                                                        met_vals[i]));
          }
          else
            met_vals[i]=0.0;
        }
        break;                                
        
      case TETRAHEDRON:
      case HEXAHEDRON:
        temp_double=element->compute_unsigned_volume(pd,err);  MSQ_ERRZERO(err);
        for(i=0;i<num_samp;++i){
          met_vals[i]=elems[adj_elems[i]].compute_unsigned_volume(pd,err);  MSQ_ERRZERO(err);
          if((temp_double+met_vals[i])>MSQ_MIN){
            met_vals[i]=fabs((temp_double-met_vals[i])/(temp_double+
                                                        met_vals[i]));
          }
          else
            met_vals[i]=0.0;
        }
        break;
      default:
      MSQ_SETERR(err)("ASM quality metric not implemented for this "
                      "element type.", MsqError::NOT_IMPLEMENTED);
      return false;
    };
    fval=average_metrics(met_vals,num_samp,err);  MSQ_ERRZERO(err);
      //PRINT_INFO("\nRETURNING %f \n",fval);
    delete []met_vals;
  }
  return true;
}


