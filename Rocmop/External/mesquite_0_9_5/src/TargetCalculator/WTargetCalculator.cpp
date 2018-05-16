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
 
/*! \file WTargetCalculators.cpp

file for the Mesquite::WTargetCalculator class

  \author Thomas Leurent
  \date   2004-05-31
 */


#include "WTargetCalculator.hpp"
#include "PatchDataUser.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"

using namespace Mesquite;



/*! The type of targets computed by this function is selected by the constructor of
    the base classes. */
void WTargetCalculator::compute_target_matrices(PatchData &pd, MsqError &err)
{
  MSQ_FUNCTION_TIMER( "WTargetCalculator::compute_target_matrice" );

  size_t num_elements=pd.num_elements();
  
  PatchData ref_pd, *ref_pd_ptr;
  if ( refMesh != 0 ) {
    // If there is a reference mesh, gets a patch ref_pd equivalent to the patch pd of the main mesh.
    PatchDataParameters ref_pd_params(this->get_all_parameters());
    refMesh->get_next_patch(ref_pd, ref_pd_params, err); MSQ_ERRRTN(err);
    // Make sure topology of ref_pd and pd are equal
    assert( num_elements == ref_pd.num_elements() );
    size_t num_vertices=pd.num_vertices();
    assert( num_vertices == ref_pd.num_vertices() );
    ref_pd_ptr = &ref_pd;
  }
  else {
    // the reference patch is the same as the working patch if there is no reference mesh.
    ref_pd_ptr = &pd;
  }
  
  MsqMeshEntity* elems = pd.get_element_array(err); MSQ_ERRRTN(err);
  MsqMeshEntity* elems_ref = ref_pd_ptr->get_element_array(err); MSQ_ERRRTN(err);

  Matrix3D W_guides[MSQ_MAX_NUM_VERT_PER_ENT];
  TargetMatrix matrices[MSQ_MAX_NUM_VERT_PER_ENT];
  
  for (size_t i=0; i<num_elements; ++i) {
    int nve = elems[i].vertex_count();
    assert( nve = elems_ref[i].vertex_count() );

    compute_guide_matrices(guideMatrix, *ref_pd_ptr, i, W_guides, nve, err); MSQ_ERRRTN(err);
    for (int c = 0; c < nve; ++c)
      matrices[c] = W_guides[c];
    pd.targetMatrices.set_element_corner_tags( &pd, i, matrices, err ); MSQ_ERRRTN(err);
  }

  //if ( refMesh != 0 ) delete ref_pd;
}

