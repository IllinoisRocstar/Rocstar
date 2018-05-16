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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file LocalSizeQualityMetric.hpp

Header file for the Mesquite::LocalSizeQualityMetric class

  \author Michael Brewer
  \date   April 9, 2003
 */


#ifndef LocalSizeQualityMetric_hpp
#define LocalSizeQualityMetric_hpp


#include "Mesquite.hpp"
#include "VolumeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"


namespace Mesquite
{
     /*! \class LocalSizeQualityMetric
       \brief Computes the local size metric for a given vertex.
       
        LocalSizeQualityMetric is a vertex based metric which computes
        the corner volume (or area) for the element corners attached
        to a given element.  Then these volumes (or areas) are averaged
        together.  The default averaging method is QualityMetric::RMS.
     */
   class MsqVertex;
   
   class LocalSizeQualityMetric : public VolumeQualityMetric
   {
  public:
        //Default constructor. 
      LocalSizeQualityMetric()
       {
         avgMethod=RMS;
         feasible=0;
         set_metric_type(QualityMetric::VERTEX_BASED);
         set_name("Local Size Quality Metric");
       }

       // virtual destructor ensures use of polymorphism during destruction
     virtual ~LocalSizeQualityMetric()
        {}
     
     
  protected:
       //!For the given vertex, vert, calculate the local size metric value.
     bool evaluate_vertex(PatchData &pd, MsqVertex *vert, double &fval,
                          MsqError &err);

  private:
       //!Calculate the area of the triangle formed by the three vertices.
     inline double compute_corner_area(PatchData &pd, size_t vert_1,
                                       size_t vert_2, size_t vert_3,
                                       MsqError &err);
     
       //!Calculate the volume of the tetrahedron formed by the four vertices.
     inline double compute_corner_volume(PatchData &pd, size_t vert_1,
                                         size_t vert_2, size_t vert_3,
                                         size_t vert_4, MsqError &err);
         
  };

   //!Calculate the area of the triangle formed by the three vertices.
   inline double LocalSizeQualityMetric::compute_corner_area(PatchData &pd,
                                                             size_t vert_1,
                                                             size_t vert_2,
                                                             size_t vert_3,
                                                             MsqError &err)
   {
     MsqVertex* verts = pd.get_vertex_array(err);
     Vector3D vec_1=verts[vert_2]-verts[vert_1];
     Vector3D vec_2=verts[vert_3]-verts[vert_1];
     Vector3D cross_vec=vec_1*vec_2;
     return (cross_vec.length()/2.0);
   }
   
   //!Calculate the volume of the tetrahedron formed by the four vertices.
   inline double LocalSizeQualityMetric::compute_corner_volume(PatchData &pd,
                                                               size_t vert_1,
                                                               size_t vert_2,
                                                               size_t vert_3,
                                                               size_t vert_4,
                                                               MsqError &err)
   {
     MsqVertex* verts = pd.get_vertex_array(err);
     Vector3D vec_1=verts[vert_2]-verts[vert_1];
     Vector3D vec_2=verts[vert_3]-verts[vert_1];
     Vector3D vec_3=verts[vert_4]-verts[vert_1];
     return fabs((vec_3%(vec_1*vec_2))/6.0);
     
   }  


} //namespace


#endif // LocalSizeQualityMetric_hpp


