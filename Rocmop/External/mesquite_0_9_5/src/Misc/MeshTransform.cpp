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
  \file   MeshTransform.cpp
  \brief  

  The MeshTransform Class is the base class for all the smoothing algorythms 

  \author Michael Brewer
  \date   2004-11-06
*/


#include "MeshTransform.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

namespace Mesquite {
/*! Constructor sets the matrix and the vector for the affine transformation.
  It also sets the patch type and culling method.  My default the patches
  are just created such that a patch is created around each vertex but
  no elements are added to the patch (each patch consists of a single vertex,
  ELEMENTS_ON_VERTEX_PATCH with a layer depth of 1).
*/
  MeshTransform::MeshTransform(Matrix3D &in_mat, Vector3D &in_vec, MsqError &err)
  {
    mMat = in_mat;
    mVec = in_vec;
    set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH ,err, 0, 1);
    no_culling_method();
  }
  void MeshTransform::set_patch_type(PatchData::PatchType patch_type, MsqError &err,
                                     int param1=0, int param2=0) {
    if(patch_type != PatchData::ELEMENTS_ON_VERTEX_PATCH ||
       param1 != 0){
      MSQ_SETERR(err)("Patch type must be ELEMENTS_ON_VERTEX_PATCH, depth 1.\n", MsqError::INVALID_ARG);
    }
    PatchDataUser::set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH ,err,
                                  0, param2);
  }
  
/*! \fn MeshTransform::loop_over_mesh(MeshSet &ms, MsqError &err)
  Actually apply the affine transformation
    \param const MeshSet &: this MeshSet is looped over. Only the
    mutable data members are changed (such as currentVertexInd).
  */
  double MeshTransform::loop_over_mesh(MeshSet &ms, MsqError &err)
  {
      //This shouldn't be an issue, but make sure the patch type is still
      //correct.
    if(get_patch_type() != PatchData::ELEMENTS_ON_VERTEX_PATCH){
      MSQ_SETERR(err)("Incorrect patch type being used.\n",
                      MsqError::INVALID_STATE );
      return 1.;
    }
      //get the first patch
    PatchData patch_data;
    bool next_patch=ms.get_next_patch(patch_data, this, err);
    assert(!err);
    if(err)
      return 1.0;
      //loop over the patches (ie, loop over the vertices.
    while( next_patch )
    {
        //make sure we have a vertex.
      if(patch_data.num_vertices()!=1){
        MSQ_SETERR(err)( "Incorrect patch depth being used.",
                         MsqError::INVALID_STATE ); 
        return 1.;
      }
      MsqVertex* vert_array = patch_data.get_vertex_array(err);
      if(err)
        return 1.0;
        //perform the affine transormation
      Vector3D temp_vec = vert_array[0];
      vert_array[0]=mMat*temp_vec+mVec;
        //update the vertex postion
      patch_data.update_mesh(err);
        //get the next patch
      next_patch =  ms.get_next_patch(patch_data, this, err);
      if(err)
        return 1.0;
    }
      //return
    return 0.;
  }
  
} // namespace Mesquite
