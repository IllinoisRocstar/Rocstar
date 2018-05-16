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
  \file   VertexConditionNumberQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include "VertexConditionNumberQualityMetric.hpp"
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

#include <math.h>
#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
   using std::vector;
#endif

using namespace Mesquite;

VertexConditionNumberQualityMetric::VertexConditionNumberQualityMetric()
{
  set_metric_type(VERTEX_BASED);
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  set_name("Vertex Condition Number");
}

bool VertexConditionNumberQualityMetric::evaluate_vertex(PatchData &pd,
                                                         MsqVertex* vert,
                                                         double &fval,
                                                         MsqError &err)
{
    //pd.generate_vertex_to_element_data();
  bool return_flag;
  fval=MSQ_MAX_CAP;
    //get the element array
  MsqMeshEntity* elems = pd.get_element_array(err);
    //convert the MsqVertex pointer into an index
  size_t this_vert = pd.get_vertex_index(vert);
    //get the vertex to element array and the offset array
  //const size_t* elem_offset = pd.get_vertex_to_elem_offset(err);  MSQ_ERRZERO(err);
  //const size_t* v_to_e_array = pd.get_vertex_to_elem_array(err);  MSQ_ERRZERO(err);
    //find the offset for this vertex
  //size_t this_offset = elem_offset[this_vert];
    //get the number of elements attached to this vertex (given by the
    //first entry in the vertex to element array)
  //size_t num_elems = v_to_e_array[this_offset];
    //PRINT_INFO("\nIN LOCAL SIZE CPP, num_elements = %i",num_elems);
    //if no elements, then return true
  size_t num_elems, *v_to_e_array;
  v_to_e_array = pd.get_vertex_element_adjacencies( this_vert, num_elems, err );
  MSQ_ERRZERO(err);  
    
  if(num_elems <= 0){
    return true;
  }
  
    //create an array to store the local metric values before averaging
    //Can we remove this dynamic allocatio?
  double* met_vals = new double[num_elems];
    //vector to hold the other verts which form a corner.
  vector<size_t> other_vertices;
  other_vertices.reserve(4);
  size_t i=0;
    //only 3 temp_vec will be sent to cond-num calculator, but the
    //additional vector3Ds may be needed during the calculations
  Vector3D temp_vec[6];
  MsqVertex *vertices=pd.get_vertex_array(err);
  //loop over the elements attached to this vertex
  for(i=0;i<num_elems;++i){
      //get the vertices connected to this vertex for this element
    elems[v_to_e_array[i]].get_connected_vertices(this_vert,
                                                                other_vertices,
                                                                err);  MSQ_ERRZERO(err);
      //switch over the element type of this element
    switch(elems[v_to_e_array[i]].get_element_type()){
    
      case TRIANGLE:
        temp_vec[0]=vertices[other_vertices[0]]-vertices[this_vert];
        temp_vec[2]=vertices[other_vertices[1]]-vertices[this_vert];
          //make relative to equilateral
        temp_vec[1]=((2*temp_vec[2])-temp_vec[0])*MSQ_SQRT_THREE_INV;
        return_flag=condition_number_2d(temp_vec,this_vert,pd,met_vals[i],err);  MSQ_ERRZERO(err);
        if(!return_flag)
          return return_flag;
        break;
      case QUADRILATERAL:
        temp_vec[0]=vertices[other_vertices[0]]-vertices[this_vert];
        temp_vec[1]=vertices[other_vertices[1]]-vertices[this_vert];
        return_flag=condition_number_2d(temp_vec,this_vert,pd,met_vals[i],err);  MSQ_ERRZERO(err);
        if(!return_flag)
          return return_flag;
        break;
      case TETRAHEDRON:
        temp_vec[0]=vertices[other_vertices[0]]-vertices[this_vert];
        temp_vec[3]=vertices[other_vertices[1]]-vertices[this_vert];
        temp_vec[4]=vertices[other_vertices[2]]-vertices[this_vert];
          //transform to equilateral tet
        temp_vec[1]=((2*temp_vec[3])-temp_vec[0])/MSQ_SQRT_THREE;
        temp_vec[2]=((3*temp_vec[4])-temp_vec[0]-temp_vec[3])/
          (MSQ_SQRT_THREE*MSQ_SQRT_TWO);
        return_flag=condition_number_3d(temp_vec,pd,met_vals[i],err);  MSQ_ERRZERO(err);
        if(!return_flag)
          return return_flag;
        break;
      case HEXAHEDRON:
        temp_vec[0]=vertices[other_vertices[0]]-vertices[this_vert];
        temp_vec[1]=vertices[other_vertices[1]]-vertices[this_vert];
        temp_vec[2]=vertices[other_vertices[2]]-vertices[this_vert];
        return_flag=condition_number_3d(temp_vec,pd,met_vals[i],err);  MSQ_ERRZERO(err);
        if(!return_flag)
          return return_flag;
        break;
      default:
        fval=MSQ_MAX_CAP;
        return false;
        
    }// end switch over element type
    other_vertices.clear();
  }//end loop over elements
  fval = average_metrics(met_vals, num_elems, err);  MSQ_ERRZERO(err);
  delete []met_vals;
  return true;
}


