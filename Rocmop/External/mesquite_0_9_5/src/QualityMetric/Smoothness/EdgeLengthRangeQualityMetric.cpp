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
/*! \file EdgeLengthRangeQualityMetric.cpp
  \author Michael Brewer
  \date 2002-05-14
  Evaluates the lengths of the edges attached to the given vertex.
  By default, the averaging method is set to SUM.
*/


#include "EdgeLengthRangeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MsqDebug.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
   using std::vector;
#endif


using namespace Mesquite;


/*!For the given vertex, vert, with connected edges of lengths l_j for
  j=1...k, the metric value is the average (where the default average
  type is SUM) of
        u_j = ( | l_j - lowVal | - (l_j - lowVal) )^2 +
              ( | highVal - l_j | - (highVal - l_j) )^2.
*/
bool EdgeLengthRangeQualityMetric::evaluate_vertex(PatchData &pd, MsqVertex* vert,
                                             double &fval, MsqError &err)
{
  fval=0.0;
  size_t this_vert = pd.get_vertex_index(vert);
  size_t other_vert;
  vector<size_t> adj_verts;
  Vector3D edg;
  pd.get_adjacent_vertex_indices(this_vert,adj_verts,err);  MSQ_ERRZERO(err);
  int num_sample_points=adj_verts.size();
  double *metric_values=new double[num_sample_points];
  MsqVertex* verts = pd.get_vertex_array(err);  MSQ_ERRZERO(err);
  int point_counter=0;
    //store the length of the edge, and the first and second component of
    //metric values, respectively.
  double temp_length=0.0;
  double temp_first=0.0;
  double temp_second=0.0;
    //PRINT_INFO("INSIDE ELR, vertex = %f,%f,%f\n",verts[this_vert][0],verts[this_vert][1],verts[this_vert][2]);
    //loop while there are still more adjacent vertices.
  while(!adj_verts.empty()){
    other_vert=adj_verts.back();
    adj_verts.pop_back();
    edg[0]=verts[this_vert][0]-verts[other_vert][0];
    edg[1]=verts[this_vert][1]-verts[other_vert][1];
    edg[2]=verts[this_vert][2]-verts[other_vert][2];
      //compute the edge length
    temp_length=edg.length();
      //get the first component
    temp_first = temp_length - lowVal;
    temp_first = fabs(temp_first) - (temp_first);
    temp_first*=temp_first;
      //get the second component
    temp_second = highVal - temp_length;
    temp_second = fabs(temp_second) - (temp_second);
    temp_second*=temp_second;
      //combine the two components
    metric_values[point_counter]=temp_first+temp_second;
      //increment the counter
    ++point_counter;
  }
    //average the metric values of the edges
  fval=average_metrics(metric_values,num_sample_points,err);  MSQ_ERRZERO(err);
    //clean up
  delete[] metric_values;
    //always return true because mesh is always valid wrt this metric.
  return true;
  
}

