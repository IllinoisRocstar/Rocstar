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
  \file   PatchData.cpp

  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-01-17
*/

#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "TargetMatrix.hpp"
#include "TargetCalculator.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <vector.h>
#  include <map.h>
#else
#  include <list>
#  include <vector>
#  include <map>
   using std::list;
   using std::map;
   using std::vector;
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <iostream.h>
#else
#  include <iostream>
   using std::ostream;
   using std::endl;
#endif

namespace Mesquite {

PatchData::PatchData()
  : targetMatrices( "MSQ_TARGET_MATRIX", Mesh::DOUBLE ),
    meshSet(0),
    domainSet(false),
    domainHint(NO_DOMAIN_HINT),
    mType(UNDEFINED_PATCH_TYPE),
    numCornerVertices(0),
    haveComputedInfos(0)
    
{
}


// Destructor
PatchData::~PatchData()
{
}


void PatchData::get_minmax_element_unsigned_area(double& min, double& max, MsqError &err)
{
  if (!have_computed_info(MAX_UNSIGNED_AREA) ||
      !have_computed_info(MIN_UNSIGNED_AREA))
  {
    max=0;
    min=MSQ_DBL_MAX;
    size_t count = num_elements();
    for (size_t i=0; i<count; ++i) {
      double vol;
      vol = elementArray[i].compute_unsigned_area(*this, err); MSQ_ERRRTN(err);
      if (vol > max)
        max = vol;
      if (vol < min)
        min = vol;
    }
    note_have_info(MAX_UNSIGNED_AREA);
    note_have_info(MIN_UNSIGNED_AREA);
    computedInfos[MAX_UNSIGNED_AREA] = max;
    computedInfos[MIN_UNSIGNED_AREA] = min;
  }
  else
  {
    max = computedInfos[MAX_UNSIGNED_AREA];
    min = computedInfos[MIN_UNSIGNED_AREA];
  }
    
  if (max <= 0 || min < 0 || min == MSQ_DBL_MAX)
    MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
}

double PatchData::get_barrier_delta(MsqError &err)
{
  double result;
  if (have_computed_info(MINMAX_SIGNED_DET3D))
  {
    result = computedInfos[MINMAX_SIGNED_DET3D];
  }
  else
  {
    double min= MSQ_DBL_MAX;
    double max=-MSQ_DBL_MAX;
    size_t count = num_elements();
    for (size_t i=0; i<count; ++i) {
      Matrix3D A[MSQ_MAX_NUM_VERT_PER_ENT];
      size_t nve = elementArray[i].vertex_count();
      elementArray[i].compute_corner_matrices(*this, A, nve, err);
      MSQ_ERRZERO(err);
      for (size_t j=0; j<nve; ++j) {
        min = det(A[j]) < min ? det(A[j]) : min;
        max = det(A[j]) > max ? det(A[j]) : max;
      }
    }

    if (max <= 0) {
      MSQ_SETERR(err)("Sigma_max is not positive.", MsqError::INVALID_MESH);
      return 0;
    }
      //We set delta to zero if everything in the initial mesh is valid.
      //  This causes metrics with a barrier between valid and inverted
      //  meshes to retain that barrier.  If there is a negative jacobian
      //  corner in the mesh, we set delta to a small fraction of the
      //  maximum jacobian in the mesh.
    result = (min<=MSQ_MIN) ? 0.001 * max : 0;
    computedInfos[MINMAX_SIGNED_DET3D] = result;
    note_have_info(MINMAX_SIGNED_DET3D);
  }
  
  return result;
}


double PatchData::get_average_Lambda_3d( MsqError &err)
{
  double avg;
  if (have_computed_info(AVERAGE_DET3D))
  {
    avg = computedInfos[AVERAGE_DET3D];
  }
  else 
  {
    avg =0.;
    int total_num_corners =0;
    Matrix3D A[MSQ_MAX_NUM_VERT_PER_ENT];
    for (size_t i=0; i<elementArray.size(); ++i) {
      int nve = elementArray[i].vertex_count();
      elementArray[i].compute_corner_matrices(*this, A, nve, err); 
      MSQ_ERRZERO(err);
      total_num_corners += nve;
      for (int c=0; c<nve; ++c) {
        avg += TargetCalculator::compute_Lambda(A[c], err); 
        MSQ_ERRZERO(err);
      }
    }

    avg = avg / total_num_corners;
    computedInfos[AVERAGE_DET3D] = avg;
    note_have_info(AVERAGE_DET3D);
  }
  return avg;
}



/*! \fn PatchData::reorder()
   Physically reorder the vertices and elements in the PatchData to improve
   the locality of reference.  This method implements a Reverse Breadth First 
   Search order starting with the vertex furthest from the origin.  Other
   orderings can also be implemented.
*/
void PatchData::reorder()
{
  const size_t numv = num_vertices();
  const size_t nume = num_elements();

  size_t *vtx;
  size_t *tmp;

  size_t *sta = new size_t[numv + 1];
  size_t *vte;
  size_t *ord = new size_t[numv];
  size_t *per = new size_t[numv];
  size_t *pel;
  size_t *que1 = new size_t[numv];
  size_t *que2 = new size_t[numv];

  //MsqVertex *v2a;
  //Mesh::VertexHandle *v2h;
  //MsqMeshEntity *e2a;
  //Mesh::ElementHandle *e2h;
    // Copy arrays so higher-order nodes get copied
  PatchDataMem<MsqVertex> v2a( vertexArray );
  PatchDataMem<Mesh::VertexHandle> v2h( vertexHandlesArray );
  PatchDataMem<MsqMeshEntity> e2a( elementArray.size() );
  PatchDataMem<Mesh::ElementHandle> e2h( elementHandlesArray.size() );

  double val, max;

  size_t toc;
  size_t vtc;
  size_t idx;
  size_t loc;
  size_t i, j;
  size_t q1l, q2l, q;
  size_t st, en;

  // Step -1: Clear any data that will be invalidated by this
  clear_tag_data();

  // Step 0:  Make sure patch data is valid.

  // Step 1:  Find the length of the element to vertex list for each 
  //          individual vertex.

  memset(sta, 0, (numv+1)*sizeof(size_t));
  for (i = 0; i < nume; ++i) {
    vtc = elementArray[i].vertex_count();
    vtx = elementArray[i].get_vertex_index_array();

    for (j = 0; j < vtc; ++j) {
      ++sta[vtx[j]];
    }
  }

  // Step 2:  Compute the offsets, total length of the element to vertex
  //          list, and allocate the data.

  toc = sta[0];
  sta[0] = 0;
  for (i = 1; i <= numv; ++i) {
    j = sta[i];
    sta[i] = toc;
    toc += j;
  }

  vte = new size_t[toc];

  // Step 3:  Finish constructing the vertex to element list.

  for (i = 0; i < nume; ++i) {
    vtc = elementArray[i].vertex_count();
    vtx = elementArray[i].get_vertex_index_array();

    for (j = 0; j < vtc; ++j) {
      vte[sta[vtx[j]]++] = i;
    }
  }

  for (i = numv; i > 0; --i) {
    sta[i] = sta[i-1];
  }
  sta[i] = 0;

  // Step 4:  Begin the reodering by computing the vertex furthest from the
  //          origin.

  max = -1.0;
  idx =  0;

  for (i = 0; i < numv; ++i) {
    val = vertexArray[i].length_squared();
    if (val > max) {
      max = val;
      idx = i+1;
    }
  }

  // Step 5:  Perform a breadth first search to find the ordering.

  memset(per, 0, numv*sizeof(size_t));

  loc = 0;
  while (idx > 0) {
    // The vertex referenced by idx has not been visited yet.  Insert it
    // into the queue for processing.
    --idx;

    q1l = 1;
    que1[0] = idx;
    per[idx] = 1;

    while (q1l) {
      q = 0;
      q2l = 0;

      while (q < q1l) {
        idx = que1[q++];
	ord[loc++] = idx;

        st = sta[idx];
        en = sta[idx+1];
        while (st < en) {
          vtc = elementArray[vte[st]].vertex_count();
          vtx = elementArray[vte[st]].get_vertex_index_array();
	  ++st;

          for (j = 0; j < vtc; ++j) {
            idx = vtx[j];
            if (!per[idx]) {
              que2[q2l++] = idx;
              per[idx] = 1;
            }
          }
        }
      }

      q1l = q2l;

      tmp  = que1;
      que1 = que2;
      que2 = tmp;
    }

    if (loc >= numv) {
      break;
    }

    // The mesh is not connected.  There is another piece with some vertices
    // remaining.  Repeat the breadth-first search algorithm on the new
    // mesh component.

    max = -1.0;
    idx =  0;
    for (i = 0; i < numv; ++i) {
      if (!per[i]) {
        val = vertexArray[i].length_squared();
        if (val > max) {
          max = val;
          idx = i+1;
        }
      }
    }
  }

  delete[] que1;
  delete[] que2;

  // Step 6:  Compute the permutation vectors

  pel = new size_t[nume];
  for (i = 0; i < nume; ++i) {
    pel[i] = nume;
  }

  toc = 0;
  for (i = 0; i < numv; ++i) {
    loc = ord[numv-1-i];

    per[loc] = i;

    st = sta[loc];
    en = sta[loc+1];
    while (st < en) {
      loc = vte[st++];
      if (nume == pel[loc]) {
        pel[loc] = toc++;
      }
    }
  }

  delete[] ord;
  delete[] vte;
  delete[] sta;

  // Step 7:  Permute the vertices

  //v2a = new MsqVertex[vertexArraySize];
  //v2h = new Mesh::VertexHandle[vertexArraySize];

  for (i = 0; i < numv; ++i) {
    loc = per[i];
    v2a[loc] = vertexArray[i];
    v2h[loc] = vertexHandlesArray[i];
  }
  
  if (!vertexNormals.empty() && domainHint != PLANAR_DOMAIN) {
    PatchDataMem<Vector3D> v2n(numv);
    for (i = 0; i < numv; ++i) {
      loc = per[i];
      v2n[loc] = vertexNormals[i];
    }
    vertexNormals = v2n;
  } 

  //delete[] vertexArray;
  //delete[] vertexHandlesArray;

  vertexArray = v2a;
  vertexHandlesArray = v2h;

  // Step 8: Permute the elements and vertex indices for the elements

  //e2a = new MsqMeshEntity[elemArraySize];
  //e2h = new Mesh::ElementHandle[elemArraySize];

  for (i = 0; i < nume; ++i) {
    vtc = elementArray[i].vertex_count();
    vtx = elementArray[i].get_vertex_index_array();

    for (j = 0; j < vtc; ++j) {
      vtx[j] = per[vtx[j]];
    }

    loc = pel[i];
    e2a[loc] = elementArray[i];
    e2h[loc] = elementHandlesArray[i];
  }

  //delete[] elementArray;
  //delete[] elementHandlesArray;

  elementArray = e2a;
  elementHandlesArray = e2h;

  // Step 9: Finish by deleting allocated memory

  delete[] per;
  delete[] pel;

  // Step 10: Recompute vertex to element mapping if it existed.
 
  if (vertAdjacencyOffsets.size()) {
    vertAdjacencyOffsets.clear();
    vertAdjacencyArray.clear();
    generate_vertex_to_element_data();
  }
  return;
}

/*! \fn PatchData::num_free_vertices()
   This function has to iterate through all the PatchData vertices to determine
   the number of free vertices. Use with care ! */
int PatchData::num_free_vertices(MsqError &/*err*/) const
{
  int num_free_vertices=0;
  
  size_t count = num_vertices();
  for (size_t i = 0; i < count; ++i)
    if (vertexArray[i].is_free_vertex())
      ++num_free_vertices;
  
  return num_free_vertices;
}

unsigned PatchData::num_free_nodes( MsqError& ) const
{
  unsigned result = 0;
  size_t count = num_nodes();
  for (size_t i = 0; i < count; ++i)
    if (vertexArray[i].is_free_vertex())
      ++result;
  return result;
}


// #undef __FUNC__
// #define __FUNC__ "PatchData::add_element"
// /*! \fn PatchData::add_element(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh, int* vertex_indices, EntityTopology topo,  MsqError &err)

// \param int* vertex_indices ... those indices corresponds to the indices of
// the element's vertices in the PatchData arrays -- see output
// of the add_vertex function.
// */
// int PatchData::add_element(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
//                            size_t* vertex_indices, EntityTopology topo,
//                            MsqError &err)
// {
//   int num_verts = MsqMeshEntity::vertex_count(topo);
//   if (!num_verts)
//     err.set_msg("Attempting to add unknown element type to PatchData.");
//   else if (numElements >= elemArraySize)
//     err.set_msg("No space available. Use reserve_element_capacity().");
//   else
//   {
//       // Set the element's type
//     elementArray[numElements].set_element_type(topo);
//     elementHandlesArray[numElements].mesh = mh;
//     elementHandlesArray[numElements].entity = eh;
//       // Go through each vertex
//     for (int n=0; n<num_verts; ++n)
//     {
//         // Make sure it's a valid index
//       if (vertex_indices[n]>=numVertices)
//         err.set_msg("invalid vertex indices");
//         // Set the element's vertex indices
//       elementArray[numElements].set_vertex_index(n, vertex_indices[n]);
//     }
//     return numElements++;
//   }
//   return -1;
// }

// #undef __FUNC__
// #define __FUNC__ "PatchData::add_triangle"
// /*! \fn PatchData::add_triangle(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh, size_t index_vtx1, size_t index_vtx2, size_t index_vtx3, MsqError &err)

// \brief adds a triangle element to the PatchData object.

// \param int index_vertex1 ... those 3 indices corresponds to the indices of
// the triangle's vertices in the PatchData arrays -- see output
// of the add_vertex function.
// */
// void PatchData::add_triangle(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
//                              size_t index_vtx1,
//                              size_t index_vtx2,
//                              size_t index_vtx3,
//                              MsqError &err)
// {
//     // make sure we've got space to add this element
//   if (elemArraySize == numElements)
//   {
//     err.set_msg("No more space in PatchData element array");
//     return;
//   }
  
//     // checks the indices are valid
//   if (index_vtx1>=numVertices || // index_vtx1<0 ||
//       index_vtx2>=numVertices || // index_vtx2<0 ||
//       index_vtx3>=numVertices // || index_vtx3<0
//       )
//     err.set_msg("invalid vertex indices");
  
//   elementHandlesArray[numElements].mesh = mh;
//   elementHandlesArray[numElements].entity = eh;
//   elementArray[numElements].set_element_type(TRIANGLE);
//   elementArray[numElements].set_vertex_index(0, index_vtx1);
//   elementArray[numElements].set_vertex_index(1, index_vtx2);
//   elementArray[numElements].set_vertex_index(2, index_vtx3);
//   ++numElements;
  
//   return;
// }


/*! \fn PatchData::move_free_vertices_constrained(Vector3D dk[], int nb_vtx, double step_size, MsqError &err)
   PatchData::move_free_vertices_constrained() moves the free vertices
   (see MsqVertex::is_free() ) as specified by the search direction (dk)
   and scale factor (step_size). After being moved, the vertices are
   projected onto the appropriate geometry.  Fixed vertices are not moved
   regardless of their corresponding dk direction.
   It is often useful to use the create_coords_momento() function before
   calling this function.
   Compile with -DMSQ_DBG3 to check that fixed vertices
   are not moved with that call.

   \param dk: must be a [nb_vtx] array of Vector3D that contains
   the direction in which to move each vertex. Fixed vertices moving
   direction should be zero, although fixed vertices will not be
   moved regardless of their corresponding dk value.
   \param nb_vtx is the number of vertices to move. must corresponds
   to the number of vertices in the PatchData.
   \param step_size will multiply the moving direction given in dk
   for each vertex.
  */
void PatchData::move_free_vertices_constrained(Vector3D dk[], size_t nb_vtx,
                                               double step_size, MsqError &err)
{
  if (nb_vtx != num_vertices())
  {
    MSQ_SETERR(err)("The directional vector must be of length numVertices.",
                    MsqError::INVALID_ARG);
    return;
  }
  
  MsqFreeVertexIndexIterator free_iter(this, err);
  free_iter.reset();
  while (free_iter.next())
  {
    vertexArray[free_iter.value()] += (step_size * dk[free_iter.value()]);
    snap_vertex_to_domain(free_iter.value(), err);
    MSQ_ERRRTN(err);
  }
  
    // Checks that moving direction is zero for fixed vertices.
  if (MSQ_DBG(3)) {
  for (size_t m=0; m<num_vertices(); ++m) {
    Vector3D zero_3d(0.,0.,0.);
    if (!vertexArray[m].is_free_vertex() 
     && dk[m] != zero_3d 
     && dk[m] != -zero_3d ) 
    {
      MSQ_DBGOUT(3) << "dk["<<m<<"]: " << dk[m] << endl;
      MSQ_DBGOUT(3) << "moving a fixed vertex." << endl;
    }
  }     
  }
}


/*! set_free_vertices_constrained is similar to 
PatchData::move_free_vertices_constrained() except the original vertex positions
are those stored in the PatchDataVerticesMemento instead of the actual vertex
position stored in the PatchData Vertex array.  The final location is stored
in the PatchData; the PatchDataVerticesMemento is unchanged.

   \param dk: must be a [nb_vtx] array of Vector3D that contains
   the direction in which to move each vertex. Fixed vertices moving
   direction should be zero, although fixed vertices will not be
   moved regardless of their corresponding dk value.
   \param nb_vtx is the number of vertices to move. must corresponds
   to the number of vertices in the PatchData.
   \param step_size will multiply the moving direction given in dk
   for each vertex.
  */
void PatchData::set_free_vertices_constrained(PatchDataVerticesMemento* memento,
                                              Vector3D dk[],
                                              size_t nb_vtx,
                                              double step_size,
                                              MsqError &err)
{
  if (nb_vtx != memento->numVertices)
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return;
  }
  
  size_t m=0;
  MsqFreeVertexIndexIterator free_iter(this, err);
  MSQ_ERRRTN(err);
  free_iter.reset();
  while (free_iter.next())
  {
    m=free_iter.value();
    vertexArray[m] = memento->vertices[m] + (step_size * dk[m]);
    snap_vertex_to_domain(m, err);
    MSQ_ERRRTN(err);
  }
  
    // Checks that moving direction is zero for fixed vertices.
  if (MSQ_DBG(3)) {
  for (m=0; m<num_vertices(); ++m)
  {
    Vector3D zero_3d(0.,0.,0.);
    if (   ! vertexArray[m].is_free_vertex()
           && ( dk[m] != zero_3d && dk[m] != -zero_3d)  ) 
    {
      MSQ_DBGOUT(3) << "dk["<<m<<"]: " << dk[m] << endl;
      MSQ_DBGOUT(3) <<"moving a fixed vertex." << endl;
    }
  }
  }
}


/*! Finds the maximum movement (in the distance norm) of the vertices in a
  patch.  The previous vertex positions are givena as a
  PatchDataVerticesMemento (memento).  The distance squared which each
  vertex has moved is then calculated, and the largest of those distances
  is returned.  This function only considers the movement of vertices
  that are currently 'free'.
  \param memento  a memento of this patch's vertex position at some
  (prior) time in the optimization.  
      */
double PatchData::get_max_vertex_movement_squared(PatchDataVerticesMemento*
                                                  memento,
                                                  MsqError &err)
{
  int m=0;
  Vector3D temp_vec;
  double temp_dist=0.0;
  double max_dist=0.0;
  MsqFreeVertexIndexIterator free_iter(this, err); MSQ_ERRZERO(err);
  free_iter.reset();
  while (free_iter.next())
  {
    m=free_iter.value();
    temp_vec=vertexArray[m] - memento->vertices[m];
    temp_dist=temp_vec.length_squared();
    if(temp_dist>max_dist)
    {
      max_dist=temp_dist;
    }
  }
  return max_dist;
}

/*!
 */
void PatchData::set_all_vertices_soft_fixed(MsqError &/*err*/)
{
  for(size_t i=0;i<num_vertices();++i)
    vertexArray[i].set_soft_fixed_flag();
}

/*!
 */
void PatchData::set_free_vertices_soft_fixed(MsqError &/*err*/)
{
  for(size_t i=0;i<num_vertices();++i){
    if(vertexArray[i].is_free_vertex())
      vertexArray[i].set_soft_fixed_flag();
  }
}

/*!
 */
void PatchData::set_all_vertices_soft_free(MsqError &/*err*/)
  {
    for(size_t i=0;i<num_vertices();++i)
      vertexArray[i].remove_soft_fixed_flag();
  }
  
/*! \fn PatchData::get_element_vertex_coordinates(size_t elem_index, vector<Vector3D> &coords, MsqError &err)

    \param elem_index The element index in the Patch
    \param coords This vector will have the coordinates appended to it.
    If necessary, make sure to clear the vector before calling the function.
  */
void PatchData::get_element_vertex_coordinates(
  size_t elem_index,
  vector<Vector3D> &coords,
  MsqError& /*err*/)
{
    // Check index
  if (elem_index >= num_elements())
    return;
  
    // Ask the element for its vertex indices
  const size_t *vertex_indices = elementArray[elem_index].get_vertex_index_array();
  
    // Get the coords for each indicated vertex
  size_t num_verts = elementArray[elem_index].vertex_count();
  coords.reserve(coords.size() + num_verts);
  for (size_t i = 0; i < num_verts; i++)
    coords.push_back(Vector3D(vertexArray[vertex_indices[i]]));
}

/*! This is an inefficient way of retrieving vertex_indices.
    Use PatchData::get_element_array followed by 
    MsqMeshEntity::get_vertex_index_array() if you don't need
    to fill an STL vector.
*/ 
void PatchData::get_element_vertex_indices(
  size_t elem_index,
  vector<size_t> &vertex_indices,
  MsqError& /*err*/)
{
    // Ask the element for its vertex indices
  elementArray[elem_index].get_vertex_indices(vertex_indices);
}


void PatchData::get_vertex_element_indices(size_t vertex_index,
                                           vector<size_t> &elem_indices,
                                           MsqError &err) 
{
  size_t count, *ptr;
  ptr = get_vertex_element_adjacencies( vertex_index, count, err );
  elem_indices.resize( count );
  memcpy( &elem_indices[0], ptr, count * sizeof(size_t) );
}

size_t* PatchData::get_vertex_element_adjacencies( size_t vertex_index,
                                                   size_t& array_len_out,
                                                   MsqError& err )
{
    // Make sure we've got the data
  if (vertAdjacencyArray.empty())
  {
    generate_vertex_to_element_data();
  }
  
  const size_t begin = vertAdjacencyOffsets[vertex_index];
  const size_t end = vertAdjacencyOffsets[vertex_index+1];
  array_len_out = end - begin;
  return &vertAdjacencyArray[begin];
}


/*!
    \brief This function fills a vector<size_t> with the indices
    to vertices connected to the given vertex by an edge.  If vert_indices
    is not initially empty, the function will not delete the current
    contents.  Instead, it will append the new indices at the end of
    the vector.

*/
void PatchData::get_adjacent_vertex_indices(size_t vertex_index,
                                            vector<size_t> &vert_indices,
                                            MsqError &err)
{
    //First get elems attached to vertex[vertex_index]
  vector<size_t> elem_indices;
  vector<size_t> temp_vert_indices;
  msq_std::vector<size_t>::iterator iter;
  size_t cur_vert;
  int found=0;
  get_vertex_element_indices(vertex_index, elem_indices,err);
  MSQ_ERRRTN(err);
  MsqMeshEntity* elems=get_element_array(err);
  MSQ_ERRRTN(err);
    //get nodes attached to vertex_index... with some duplication
  while(!elem_indices.empty()){
    elems[elem_indices.back()].get_connected_vertices(vertex_index, temp_vert_indices,err); 
    MSQ_ERRRTN(err);
    elem_indices.pop_back();
  }
    //eliminate duplication.
  while(!temp_vert_indices.empty()){
    cur_vert=temp_vert_indices.back();
    temp_vert_indices.pop_back();
    iter=vert_indices.begin();
    found=0;
    while(iter!=vert_indices.end() && !found){
      if(*(iter)==cur_vert)
        found=1;
      ++iter;
    }
    if(!found)
      vert_indices.push_back(cur_vert);
  }
}

/*! Fills a vector of indices into the entities array. The entities
    in the vector are connected the given entity (ent_ind) via an
    n-diminsional entity (where 'n' is a given integer).
    Thus, if n = 0, the entities must be connected via a vertex.
    If n = 1, the entities must be connected via an edge.
    If n = 2, the entities must be connected via a two-dimensional element.
    NOTE:  if n is 2 and the elements in the entity array are
    two-dimensional, no entities should meet this criterion.
    The adj_ents vector is cleared at the beginning of the call.

*/
void PatchData::get_adjacent_entities_via_n_dim(int n, size_t ent_ind,
                                                vector<size_t> &adj_ents,
                                                MsqError &err)
{
  //reset the vector
  adj_ents.clear();
    //vertices of this entity (given by ent_ind)
  vector<size_t> verts;
    //vector to store elements attached to the vertices in verts
  vector<size_t> elem_on_vert[MSQ_MAX_NUM_VERT_PER_ENT];
    //length of above vectos
  int length_elem_on_vert[MSQ_MAX_NUM_VERT_PER_ENT];
    //get verts on this element
  get_element_vertex_indices(ent_ind, verts, err);
  int num_vert=verts.size();
  int i=0;
  int j=0;
  for(i=0;i<num_vert;++i){
      //get elements on the vertices in verts and the number of vertices
    get_vertex_element_indices(verts[i],elem_on_vert[i],err);
    MSQ_ERRRTN(err);
    length_elem_on_vert[i]=elem_on_vert[i].size();
  }
    //this_ent is the index for an entity which is a candidate to be placed
    //into adj_ents
  size_t this_ent;
    //num of times this_ent has been found in the vectors of entity indices
  int counter=0;
  i = 0;
    //loop of each vert on ent_ind
  while(i<num_vert){
      //loop over each ent connected to vert
    j=0;
    while(j<length_elem_on_vert[i]){
        //get candidate element
      this_ent=elem_on_vert[i][j];
        //if we haven't already consider this ent
      if(this_ent!=ent_ind){
          //if this_ent occurred earlier we would have already considered it
          //so start at i and j+1
        int k1=i;
        int k2=j+1;
          //this_ent has occured once so far
        counter=1;
          //while k1 < num_vert
        while(k1<num_vert){
            //loop over entries in the elem on vert vector
          while(k2<length_elem_on_vert[k1]){
              //if it matches this_ent
            if(elem_on_vert[k1][k2]==this_ent){
                //mark it as 'seen', by making it the same as ent_ind
                //i.e., the entity  passed to us.
              elem_on_vert[k1][k2]=ent_ind;
              ++counter;
                //do not look at remaining elems in this vector
              k2+=length_elem_on_vert[k1];
            }
            else
              ++k2;
          }
          ++k1;
          k2=0;
          
        }
          //if this_ent occured enough times and isn't ent_ind
        if(counter>n && this_ent!=ent_ind){
          adj_ents.push_back(this_ent);
        }
      }
      ++j;
    }
    ++i;
  }
}

    
  


/*! \fn PatchData::update_mesh(MsqError &err)

    \brief This function copies to the TSTT mesh  the changes made to the
    free vertices / elements of the PatchData object.

*/
void PatchData::update_mesh(MsqError &err)
{
  if (!meshSet)
    return;

  meshSet->update_mesh(*this, err);
  MSQ_CHKERR(err);
}

void PatchData::generate_vertex_to_element_data()
{
  MSQ_FUNCTION_TIMER( "PatchData::generate_vertex_to_element_data" );
  
    // Skip if data already exists
  if (!vertAdjacencyArray.empty())
    return;
  
    // Allocate offset array
  vertAdjacencyOffsets.resize( num_nodes() + 1 );
  memset( &vertAdjacencyOffsets[0], 0, sizeof(size_t)*vertAdjacencyOffsets.size() );
  
    // Temporarily use offsets array to hold per-vertex element count
  PatchDataMem<MsqMeshEntity>::iterator elem_iter;
  const PatchDataMem<MsqMeshEntity>::iterator elem_end = elementArray.end();
  for (elem_iter = elementArray.begin(); elem_iter != elem_end; ++elem_iter)
  {
    size_t* conn_iter = elem_iter->get_vertex_index_array();
    const size_t* conn_end = conn_iter + elem_iter->node_count();
    for ( ; conn_iter != conn_end; ++conn_iter )
      ++vertAdjacencyOffsets[*conn_iter];
  }
  
    // Convert counts to end indices.
    // When done, vertAdjacencyOffsets will contain, for each vertex,
    // one more than the *last* index for that vertex's data in the
    // adjacency array.  This is *not* the final state for this data.
    // See comments for next loop.
  PatchDataMem<size_t>::iterator off_iter = vertAdjacencyOffsets.begin();
  const PatchDataMem<size_t>::iterator off_end = vertAdjacencyOffsets.end();
  size_t prev = *off_iter;
  ++off_iter;
  for ( ; off_iter != off_end; ++off_iter)
  {
    prev += *off_iter;
    *off_iter = prev;
  }
  
    // Allocate space for element numbers
  const size_t num_vert_uses = vertAdjacencyOffsets[num_nodes()-1];
  assert( num_vert_uses == elemConnectivityArray.size() );
  vertAdjacencyArray.resize( num_vert_uses );
  
    // Fill vertAdjacencyArray, using the indices in vertAdjacencyOffsets
    // as the location to insert the next element number in
    // vertAdjacencyArray.  When done, vertAdjacenyOffsets will contain
    // the start index for each vertex, rather than one past the last
    // index.
  for (size_t i = 0; i < elementArray.size(); ++i)
  {
    size_t* conn_iter = elementArray[i].get_vertex_index_array();
    const size_t* conn_end = conn_iter + elementArray[i].node_count();
    for ( ; conn_iter != conn_end; ++conn_iter )
    {
      const size_t array_index = --vertAdjacencyOffsets[*conn_iter];
      vertAdjacencyArray[array_index] = i;
    }
  }
  
    // Last entry should be number of vertex uses (one past the
    // last index of the last vertex.)
  vertAdjacencyOffsets[num_nodes()] = num_vert_uses;
}

void PatchData::get_subpatch(size_t center_vertex_index,
                             PatchData &subpatch,
                             MsqError &err)
{
    // Make sure we're in range
  if (center_vertex_index >= num_vertices())
  {
    MSQ_SETERR(err)("Invalid index for center vertex",MsqError::INVALID_ARG);
    return;
  }

    // Map vertex indices from this patch data to indices in new patch data
  msq_std::map<size_t, size_t> vertex_index_map;
  
  size_t num_elems;
  const size_t* vertex_adjacencies;
  vertex_adjacencies = get_vertex_element_adjacencies( center_vertex_index,
                                                       num_elems, err );
  MSQ_ERRRTN(err);
  
    // Loop through each element, populating vertex_index_map.
  size_t which_elem, which_vert, vertex_count = 0;
  size_t vertex_uses = 0;
  for (which_elem = 0; which_elem < num_elems; ++which_elem )
  {
    size_t elem_index = vertex_adjacencies[which_elem];
    MsqMeshEntity& elem = elementArray[elem_index];
    for (which_vert = elem.node_count(); which_vert--; )
    {
      size_t vert_index = elem.get_vertex_index(which_vert);
      if (vertex_index_map.find(vert_index) != vertex_index_map.end())
        vertex_index_map[vert_index] = vertex_count++;
    }
    vertex_uses += elem.node_count();
  }
  
    // Allocate storage in subpatch
  subpatch.vertexHandlesArray.resize( vertex_count );
  subpatch.elementArray.resize( num_elems );
  subpatch.elementHandlesArray.resize( num_elems );
  subpatch.elemConnectivityArray.resize( vertex_uses );
  
    // For now, put reverse of index map into handles array (such that
    // the handles array in the subpatch contains the index in this
    // patch for each vertex.)  When PatchData::initalize_data re-orders 
    // the vertices, we can then use this array to determine which vertices
    // where placed where.
  msq_std::map<size_t,size_t>::iterator iter;
  for (iter = vertex_index_map.begin(); iter != vertex_index_map.end(); ++iter)
    subpatch.vertexHandlesArray[iter->second] = (Mesh::VertexHandle)(iter->first);
  
    // Store elements with updated connectivity list in subpatch
  msq_std::vector<size_t> elem_conn_offsets(vertex_uses);
  size_t* elem_connectivity = &(subpatch.elemConnectivityArray[0]);
  size_t elem_conn_index = 0;
  for (which_elem = 0; which_elem < num_elems; ++which_elem)
  {
    size_t elem_index = vertex_adjacencies[which_elem];
    MsqMeshEntity& elem = elementArray[elem_index];
    elem_conn_offsets[which_elem] = elem_conn_index;
    for (which_vert = 0; which_vert < elem.node_count(); which_vert++ )
    {
      size_t vert_index = elem.get_vertex_index(which_vert);
      
        // Add this vertex to the new element's array
      elem_connectivity[elem_conn_index++] = vertex_index_map[vert_index];
    }
    subpatch.elementArray[which_elem].set_element_type( elem.get_element_type() );
  }
  
    // Re-order vertices and initialize other data in subpatch
  subpatch.initialize_data( &elem_conn_offsets[0], err ); MSQ_ERRRTN(err);
  
    // Copy vertex data into subpatch.  subpatch.vertexHandlesArray contains
    // the indices into this PatchData for each vertex
  subpatch.vertexArray.resize( vertex_count );
  for (which_vert = 0; which_vert < vertex_count; ++which_vert)
  {
    size_t vert_index = (size_t)(subpatch.vertexHandlesArray[which_vert]);
    subpatch.vertexHandlesArray[which_vert] = vertexHandlesArray[vert_index];
    subpatch.vertexArray[which_vert] = vertexArray[vert_index];
  }
}

//! Adjust the position of the specified vertex so that it
//! lies on its constraining domain.  The actual domain constraint
//! is managed by the TSTT mesh implementation
void PatchData::snap_vertex_to_domain(size_t vertex_index, MsqError &err)
{
  if (meshSet && meshSet->get_domain_constraint())
  {
    if (domainHint == SMOOTH_DOMAIN && !vertexNormals.empty())
    {
      meshSet->get_domain_constraint()->closest_point( 
                                          vertexHandlesArray[vertex_index],
                                          Vector3D(vertexArray[vertex_index]),
                                          vertexArray[vertex_index],
                                          vertexNormals[vertex_index],
                                          err ); MSQ_ERRRTN(err);
    }
    else
    {
      meshSet->get_domain_constraint()->snap_to(vertexHandlesArray[vertex_index],
                                                vertexArray[vertex_index]);
    }
  }
}


void PatchData::get_domain_normal_at_vertex(size_t vertex_index,
                                   bool normalize,
                                   Vector3D &surf_norm,
                                   MsqError &err) 
{
  if (domainHint != NO_DOMAIN_HINT)
  {
    if (vertexNormals.empty())
    {
      update_cached_normals( err ); MSQ_ERRRTN(err);
    }
    
    if (domainHint == PLANAR_DOMAIN)
      vertex_index = 0;
    surf_norm = vertexNormals[vertex_index];
  }
  else
  {
    MeshDomain* domain = meshSet ? meshSet->get_domain_constraint() : 0;
    if (!domain)
    {
      MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
      return;
    }
  
    surf_norm = vertexArray[vertex_index];
    domain->normal_at( vertexHandlesArray[vertex_index], surf_norm );
  }

  if (normalize)
    surf_norm.normalize();
}


void PatchData::update_cached_normals( MsqError& err )
{
  MeshDomain* domain = meshSet ? meshSet->get_domain_constraint() : 0;
  if (!domain)
  {
    MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
    vertexNormals.clear();
    return;
  }
  
  if (domainHint == PLANAR_DOMAIN)
  {
    vertexNormals.resize(1);
    vertexNormals[0] = vertex_by_index(0);
    domain->normal_at( elementHandlesArray[0], vertexNormals[0] );
  }
      
  else
  {
    vertexNormals.resize( num_vertices() );
    for ( unsigned i = 0; i < num_vertices(); ++i)
    {
      Vector3D& norm = vertexNormals[i];
      norm = vertexArray[i];
      domain->normal_at( vertexHandlesArray[i], norm );
      //norm.normalize();
      //if (!finite(norm.x()) || !finite(norm.y()) || !finite(norm.z()))
      //  norm.set(0,0,0);
    }
  }
}


void PatchData::get_domain_normal_at_element(size_t elem_index,
                                             Vector3D &surf_norm,
                                             MsqError &err) const
{
  if (meshSet && meshSet->get_domain_constraint())
  {
    elementArray[elem_index].get_centroid(surf_norm, *this, err); 
    MSQ_ERRRTN(err);
    meshSet->get_domain_constraint()->normal_at(
      elementHandlesArray[elem_index],
      surf_norm);
  }
  else
    MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
}
/*
void PatchData::get_domain_normal_at_corner( size_t elem_index,
                                             size_t elem_corner,
                                             Vector3D& normal_out,
                                             MsqError& err ) const
{
  size_t vert_index;
  
  if (meshSet && meshSet->get_domain_constraint())
  {
    vert_index = elementArray[elem_index].get_vertex_index_array()[elem_corner];
    normal_out = vertexArray[ vert_index ];
    meshSet->get_domain_constraint()->normal_at( elementHandlesArray[elem_index],
                                                 normal_out );
  }
  else
  {
    MSQ_SETERR(err)("No domain constraint set.", MsqError::INVALID_STATE );
  }
}
*/

void PatchData::get_domain_normals_at_corners( size_t elem_index,
                                               Vector3D normals_out[],
                                               MsqError& err ) 
{
  const MsqMeshEntity& elem = elementArray[elem_index];
  const unsigned count = elem.vertex_count();
  const size_t* const vertex_indices = elem.get_vertex_index_array();

  if (domainHint != NO_DOMAIN_HINT)
  {
    if (vertexNormals.empty())
    {
      update_cached_normals( err ); MSQ_ERRRTN(err);
    }
    
    if (domainHint == PLANAR_DOMAIN)
    {
      for (unsigned i = 0; i < count; ++i)
        normals_out[i] = vertexNormals[0];
    }
    else
    {
      for (unsigned i = 0; i < count; ++i)
        normals_out[i] = vertexNormals[vertex_indices[i]];
    }
  }
  else
  {
    MeshDomain* domain = meshSet ? meshSet->get_domain_constraint() : 0;
    if (!domain)
    {
      MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
      return;
    }
    
    for (unsigned i = 0; i < count; ++i)
      normals_out[i] = vertexArray[ vertex_indices[i] ];
    
    domain->normal_at( elementHandlesArray[elem_index], 
                       normals_out, count,
                       err ); MSQ_CHKERR(err);
  }
}  
  

void PatchData::set_mesh_set(MeshSet* ms)
{ meshSet = ms;
  if (ms->get_domain_constraint()!=NULL) domainSet = true; }

ostream& operator<<( ostream& stream, const PatchData& pd )
{
   size_t i;
   
   stream << "Vertices: " << endl;
   for (i = 0; i < pd.num_nodes(); ++i)
   {
      if (i == pd.num_vertices())
        stream << "Higher-Order Nodes: " << endl;
      
      stream << i << ". (" 
             << pd.vertexArray[i].x() << ","
             << pd.vertexArray[i].y() << ","
             << pd.vertexArray[i].z()
             << ") ";
      if (pd.vertexArray[i].is_flag_set( MsqVertex::MSQ_SOFT_FIXED ))
        stream << "S";
      if (pd.vertexArray[i].is_flag_set( MsqVertex::MSQ_HARD_FIXED ))
        stream << "H";
      if (pd.vertexArray[i].is_flag_set( MsqVertex::MSQ_COORDS_CHANGED ))
        stream << "C";
      
      if (pd.vertAdjacencyArray.size())
      {
        size_t j = pd.vertAdjacencyOffsets[i];
        size_t end = pd.vertAdjacencyOffsets[i+1];
        for ( ; j < end; ++j )
          stream << " " << pd.vertAdjacencyArray[j];
      }
      
      stream << endl;
   }
   
   stream << "Elements: " << endl;
   for (i = 0; i < pd.num_elements(); ++i)
   {
      stream << i << ". ";
      switch (pd.elementArray[i].get_element_type()) {
        case POLYGON:       stream << "Polygon";    break;
        case TRIANGLE:      stream << "Tri";        break;
        case QUADRILATERAL: stream << "Quad";       break;
        case POLYHEDRON:    stream << "Polyhedron"; break;
        case TETRAHEDRON:   stream << "Tet";        break;
        case HEXAHEDRON:    stream << "Hex";        break;
        case PRISM:         stream << "Wedge";      break;
        case PYRAMID:       stream << "Pyr";        break;
        default:            stream << "Unknown";    break;
      }
      stream << pd.elementArray[i].node_count() << ": ";
      for (size_t j = 0; j < pd.elementArray[i].node_count(); ++j)
        stream << pd.elementArray[i].get_vertex_index_array()[j] << " ";
      stream << endl;
   }
   stream << endl;
   
   stream << "MeshSet: " << (pd.meshSet?"yes":"no") << endl;
   stream << "domainSet: " << (pd.domainSet?"true":"false") << endl;
   stream << "mType: " << (pd.mType==PatchData::VERTICES_ON_VERTEX_PATCH?"vert-on-vert":
                           pd.mType==PatchData::ELEMENTS_ON_VERTEX_PATCH?"elem-on-vert":
                           pd.mType==PatchData::GLOBAL_PATCH?"global":"unknown") << endl;
   
   if (pd.haveComputedInfos)
   {
     stream << "ComputedInfos:" << endl;
     if (pd.have_computed_info(PatchData::MIN_UNSIGNED_AREA))
       stream << "\t MIN_UNSINGED_AREA = " << pd.computedInfos[PatchData::MIN_UNSIGNED_AREA] << endl;
     if (pd.have_computed_info(PatchData::MAX_UNSIGNED_AREA))
       stream << "\t MAX_UNSIGNED_AREA = " << pd.computedInfos[PatchData::MAX_UNSIGNED_AREA] << endl;
     if (pd.have_computed_info(PatchData::MIN_EDGE_LENGTH))
       stream << "\t MIN_EDGE_LENGTH = " << pd.computedInfos[PatchData::MIN_EDGE_LENGTH] << endl;
     if (pd.have_computed_info(PatchData::MAX_EDGE_LENGTH))
       stream << "\t MAX_EDGE_LENGTH = " << pd.computedInfos[PatchData::MAX_EDGE_LENGTH] << endl;
     if (pd.have_computed_info(PatchData::MINMAX_SIGNED_DET2D))
       stream << "\t MINMAX_SIGNED_DET2D = " << pd.computedInfos[PatchData::MINMAX_SIGNED_DET2D] << endl;
     if (pd.have_computed_info(PatchData::MINMAX_SIGNED_DET3D))
       stream << "\t MINMAX_SIGNED_DET3D = " << pd.computedInfos[PatchData::MINMAX_SIGNED_DET3D] << endl;
     if (pd.have_computed_info(PatchData::AVERAGE_DET3D))
       stream << "\t AVERAGE_DET3D = " << pd.computedInfos[PatchData::AVERAGE_DET3D] << endl;
  }
  
  return stream << endl;
}


void PatchData::initialize_data( size_t* elem_offset_array, MsqError& err )
{
  if (numCornerVertices)
  {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return;
  }
  
    // Clear out data specific to patch
  clear_tag_data();
  if (domainHint != PLANAR_DOMAIN)
    vertexNormals.clear();
  
    // Clear any vertex->element adjacency data.  It
    // is probably invalid, and certainly will be by the time
    // this function completes if the mesh contains higher-order
    // elements.
  vertAdjacencyArray.clear();
  vertAdjacencyOffsets.clear();
  
    // Initialize connectivity data in each element
  size_t i, j;
  bool higher_order = false;
  for (i = 0; i < elementArray.size(); ++i)
  {
    size_t start = elem_offset_array[i];
    size_t conn_len = elem_offset_array[i+1] - start;
    assert(conn_len > 0);
    elementArray[i].set_connectivity( &elemConnectivityArray[start], conn_len );
    if (conn_len != elementArray[i].vertex_count())
      higher_order = true;
  }
  
    // If no higher-order elements, then we're done
  if (!higher_order)
  {
      // All nodes are corner vertices
    numCornerVertices = vertexHandlesArray.size();
    return;
  }
  
    // Need to move higher-order nodes to end of list.
  
    // Use vertAdjacencyOffsets array as temporary storage.
  vertAdjacencyOffsets.resize( vertexHandlesArray.size() + 1 );
  size_t* vertex_index_map = &vertAdjacencyOffsets[0];
  
    // Note which are mid-nodes (1) and which are corner vertices (0)
  for (i = 0; i < elementArray.size(); ++i)
  {
    MsqMeshEntity& elem = elementArray[i];
    size_t* conn_array = elem.get_vertex_index_array();
    size_t num_vertices = elem.vertex_count();
    size_t num_nodes = elem.node_count();
    for (j = 0; j < num_vertices; ++j)
      vertex_index_map[ conn_array[j] ] = 0;
    for ( ; j < num_nodes; ++j)
      vertex_index_map[ conn_array[j] ] = 1;
  }
  
    // Shuffle nodes around such that all higher-order nodes
    // are at the end of the array.
    // Store new index for each node in index_map, replacing
    // flag value currently there for each node.
  i = 0;
  j = vertexHandlesArray.size() - 1;
  while (i < j)
  {
    if (!vertex_index_map[i])
    {
      // Is a corner vertex in first part of array, skip it.
      vertex_index_map[i] = i;
      ++i;
    }
    else if (vertex_index_map[j])
    {
      // Is a mid-node in latter part of array, skip it
      vertex_index_map[j] = j;
      --j;
    }
    else
    {
      // Swap mid-node from first part of array with 
      // corner vertex from latter part of array.
      msq_std::swap( vertexHandlesArray[i], vertexHandlesArray[j] );
      
      vertex_index_map[i] = j;
      vertex_index_map[j] = i;
      ++i;
      --j;
    }
  }
  
    // Finish up - set numCornerVertices to indicate 
    // where corner vertices end and mid-nodes begin in
    // vertexArray, and get the handle remaining vertex
    // if missed in the above loop.
  if (i > j)
  {
    numCornerVertices = i;
  }
  else // (i == j)
  {
    if (vertex_index_map[i])
      numCornerVertices = i;
    else
      numCornerVertices = i+1;
    vertex_index_map[i] = i;
  }
  
    // Update element connectivity data for new vertex indices
  for (i = 0; i < elemConnectivityArray.size(); ++i)
    elemConnectivityArray[i] = vertex_index_map[elemConnectivityArray[i]];
    
  vertAdjacencyOffsets.clear();
}

void PatchData::allocate_storage( size_t vertex_count,
                                  size_t element_count,
                                  size_t vertex_use_count,
                                  MsqError& err )
{
  vertexArray.resize( vertex_count );
  vertexHandlesArray.resize( vertex_count );
  elementArray.resize( element_count );
  elementHandlesArray.resize( element_count );
  elemConnectivityArray.resize( vertex_use_count );
  clear_tag_data();
  numCornerVertices = 0;
}

void PatchData::clear_tag_data()
{
  targetMatrices.clear();
}


size_t PatchData::num_corners() 
{
  size_t result = 0;
  for (unsigned i = 0; i < elementArray.size(); ++i)
    result += elementArray[i].vertex_count();
  return result;
}


} // namespace Mesquite
