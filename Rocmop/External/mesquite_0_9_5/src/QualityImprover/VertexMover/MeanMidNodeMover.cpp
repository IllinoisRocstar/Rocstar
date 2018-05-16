/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#include "MeanMidNodeMover.hpp"
#include "MsqError.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "TopologyInfo.hpp"

namespace Mesquite {

MeanMidNodeMover::MeanMidNodeMover() 
{
  MsqError err;
  PatchDataUser::set_patch_type( PatchData::GLOBAL_PATCH, err );
  MSQ_CHKERR(err);
}

MeanMidNodeMover::~MeanMidNodeMover() {}

void MeanMidNodeMover::set_patch_type( PatchData::PatchType, MsqError &err, int , int )
{
  MSQ_SETERR(err)("Patch type cannot be changed for MeanMidNodeMover.",MsqError::INVALID_ARG);
  return;
}

msq_std::string MeanMidNodeMover::get_name()
{ return "MeanMidNodeMover"; }

PatchDataUser::AlgorithmType MeanMidNodeMover::get_algorithm_type()
{ return PatchDataUser::QUALITY_IMPROVER; }
  

double MeanMidNodeMover::loop_over_mesh( MeshSet& ms, MsqError& err )
{
  if (get_patch_type() != PatchData::GLOBAL_PATCH) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return 0;
  }
  
  if (!get_global_patch()) {
    MSQ_SETERR(err)("No global patch", MsqError::INVALID_STATE);
    return 0;
  }
  
  PatchData& pd = *get_global_patch();
  
    // For each higher-order node in the patch
  for (size_t vtx = pd.num_vertices(); vtx < pd.num_nodes(); ++vtx)
  {
      // Get an adjacent element
    size_t num_elems, *elem_list;
    elem_list = pd.get_vertex_element_adjacencies( vtx, num_elems, err ); MSQ_ERRZERO(err);
    if (num_elems < 1) // Mid-node without adjacent elements????
      continue;
    size_t element_index = elem_list[0];
    MsqMeshEntity& element = pd.element_by_index( element_index );
    
      // Find position of current vertex in element connectivity list
    unsigned elem_vert_index;
    size_t num_verts = element.node_count(); 
    size_t* vert_list = element.get_vertex_index_array();
    for (elem_vert_index = 0; elem_vert_index < num_verts; ++elem_vert_index)
      if (vert_list[elem_vert_index] == vtx)
        break;
    
      // Sanity checks 
      // Mesh is corrupt if vertex not in connectivity list
    if (elem_vert_index == num_verts)
    {
      MSQ_SETERR(err)("Inconsistent connectivity/adjacency data.", MsqError::INVALID_MESH );
      return 0;
    }
      // Must be appropriate element type (e.g. can't have higher-order
      // nodes on a polygon)
    EntityTopology topo = element.get_element_type();
    unsigned num_corners = TopologyInfo::corners( topo );
    if (!num_corners)
    {
      MSQ_SETERR(err)( MsqError::INVALID_MESH, 
                       "Element type %d cannot have higher-order nodes.",
                       (int)topo );
      return 0;
    }
      // The current vertex must be one of the higher-order nodes of the element
    if (elem_vert_index < num_corners)
    {
      MSQ_SETERR(err)( "Invalid mid-node flag for mesh (mixed connectivity?)", 
                       MsqError::INVALID_MESH );
      return 9;
    }
  
      // Get the element "side" that the node is a mid-node of.
    unsigned side, dimension;
    TopologyInfo::side_number( topo, element.node_count(), elem_vert_index,
                               side, dimension, err ); MSQ_ERRZERO(err);
      // Get the indices of the vertices defining the "side"
    unsigned side_size;
    const unsigned* side_indices = 
      TopologyInfo::side_vertices( topo, side, dimension, side_size, err ); MSQ_ERRZERO(err);
    if (!side_size)
    {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return 0;
    }
  
      // Calculate average position of side vertices
    const msq_stdc::size_t* conn_array = element.get_vertex_index_array();
    Vector3D pos(0,0,0);
    for (unsigned i = 0; i < side_size; ++i)
      pos += pd.vertex_by_index( conn_array[side_indices[i]] );
    pos /= side_size;

      // Set new vertex position
    pd.set_vertex_coordinates( pos, vtx, err ); MSQ_ERRZERO(err);
    pd.snap_vertex_to_domain( vtx, err );       MSQ_ERRZERO(err);
  }
  
  pd.update_mesh(err); MSQ_ERRZERO(err);
  
  return 0;
}
                       
  
  
} // namespace Mesquite

  
  
  
