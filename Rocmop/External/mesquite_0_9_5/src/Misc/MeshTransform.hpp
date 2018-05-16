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
  \file   MeshTransform.hpp
  \brief  

  Class for performing an affine transformation on the mesh.

  \author Michael Brewer      
  \date   2004-11-06
*/

#ifndef Mesquite_MeshTransform_hpp 
#define Mesquite_MeshTransform_hpp


#include "Mesquite.hpp"
#include "QualityImprover.hpp"
#include "PatchData.hpp"
#include "ObjectiveFunction.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"

namespace Mesquite
{

  /*! \class MeshTransform
    Perform an Affine transformation on Mesh vertex positions.
    Essentially define the new vertex position, v_new, from the original
    vertex position, v_old, s.t.
    v_new = (mMat * v_old) + mVec,
    where mMat is a constant matrix and mVec is a constant vector.
   */  
  class MeshTransform : public PatchDataUser 
  {
  public:
    MeshTransform(Matrix3D &in_mat, Vector3D &in_vec, MsqError &err);

      // virtual destructor ensures use of polymorphism during destruction
    virtual ~MeshTransform() { };
      //virtual functions from PatchDataUser...
      //!Loop over the mesh and perform the affine transformation
    virtual double loop_over_mesh(MeshSet &ms, MsqError &err);
      //! Return the name of this PatchDataUser:  Mesh Transform
    virtual msq_std::string get_name() { return "Mesh Transform.";}
      //! Return the AlgorithmType of thsi PatchDataUser:  MESH_TRANSFORM
    virtual AlgorithmType get_algorithm_type() { return MESH_TRANSFORM; }
      //! MeshTransform is built to use ELEMETNS_ON_VERTEX_PATCH with a depth
      //! one.  We implement set_patch_type to ensure this is the patch that
      //! is used.
    virtual void set_patch_type(PatchData::PatchType patch_type, MsqError &err,
                                int param1, int param2);
    
  private:
    Matrix3D mMat;//!Matrix for the affine transformation
    Vector3D mVec;//!Vector for the affine transformation
  };

  
} // namespace
#endif // Mesquite_MeshTransform_hpp
