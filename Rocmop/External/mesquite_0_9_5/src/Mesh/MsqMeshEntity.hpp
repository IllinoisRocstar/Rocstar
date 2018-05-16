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

/*! \file MsqMeshEntity.hpp

\author Darryl Melander
\author Thomas Leurent
\author Michael Brewer

*/

#ifndef MSQMESHENTITY_HPP
#define MSQMESHENTITY_HPP

#include "Mesquite.hpp"
#include "QualityMetric.hpp"
#include "TopologyInfo.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
#endif

#ifdef MSQ_USE_OLD_C_HEADERS
#  include <string.h>
#  include <assert.h>
#else
#  include <cstring>
#  include <cassert>
#endif


namespace Mesquite
{
  class PatchData;

  /*!
      \class MsqMeshEntity
      \brief MsqMeshEntity is the Mesquite object that stores information about
      the elements in the mesh.

      
  */
  class MsqMeshEntity
  {
  public:
    
    MsqMeshEntity()
      : mType(MIXED), vertexIndices(0), numVertexIndices(0)
      {}

      //! Returns element type
    inline EntityTopology get_element_type() const
      { return mType; }
    
      //! Returns the number of vertices in this element,
      //! based on its element type.
    inline msq_stdc::size_t vertex_count() const;
      //! Return number of nodes in element (number of corner
      //! vertices + number of higher-order nodes).
    inline msq_stdc::size_t node_count() const 
      { return numVertexIndices; }
    
      //! gets the vertices of the mesh entity
    void get_vertex_indices(msq_std::vector<msq_stdc::size_t> &vertex_list) const;
    void append_vertex_indices(msq_std::vector<msq_stdc::size_t> &vertex_list) const;
      //! gets the vertices of the mesh entity
    void get_node_indices(msq_std::vector<msq_stdc::size_t> &vertex_list) const;
    void append_node_indices(msq_std::vector<msq_stdc::size_t> &vertex_list) const;
    //! Very efficient retrieval of vertices indexes 
    //! (corresponding to the PatchData vertex array).
    inline const msq_stdc::size_t *get_vertex_index_array() const;
    inline msq_stdc::size_t *get_vertex_index_array();
    
      //! Sets element data
    void set_element_type(EntityTopology type)
      { mType = type; }

      //! Set connectivity data (vertex array) for element.
      //! MsqMeshEntity keeps the pointer to the passed array, it
      //! does not copy it.  The caller is still responsible for
      //! releasing the memory for the passed index array after
      //! the MsqMeshEntity is destroyed.  The intention is that
      //! this is a pointer to a portion of a larger connectivity array
      //! managed by the owning \ref PatchData.
    void set_connectivity( msq_stdc::size_t *indices, size_t num_vertices);

    msq_stdc::size_t get_vertex_index(msq_stdc::size_t vertex_in_element) const;
    
      //fills array of Vector3D's with the jacobian vectors and the 
      //number of jacobian vectors
    void compute_weighted_jacobian(PatchData &pd, Vector3D& sample_point,
                                   Vector3D jacobian_vectors[],
                                   short &num_jacobian_vectors,
                                   MsqError &err );
    
      //!Returns a list of sample points given an evaluationmode 
    void get_sample_points(QualityMetric::ElementEvaluationMode mode,
                           msq_std::vector<Vector3D> &coords,
                           MsqError &err);

    //! Returns the centroid of the element.
    void get_centroid(Vector3D& centroid, const PatchData &pd, MsqError &err) const;
    
      //!Fills a vector<size_t> with vertices connected to the given
      //!vertex through the edges of this MsqMeshEntity.
    void get_connected_vertices(msq_stdc::size_t vertex_index,
                                msq_std::vector<msq_stdc::size_t> &vert_indices,
                                MsqError &err);
    
      //!Computes the area of the element.
      //!The returned value is always non-negative.
    double compute_unsigned_area(PatchData &pd, MsqError &err );
    
      //!Computes the volume of the element.
      //!The returned value is always non-negative.
    double compute_unsigned_volume(PatchData &pd, MsqError &err );
    
      //!Computes the signed area of the element.
    double compute_signed_area(PatchData &pd, MsqError &err );

      //!Computes the signed volume of the element.
    double compute_signed_volume(PatchData &pd, MsqError &err );
    
      //! Uses a MeshDomain call-back function to compute the normal at the corner.
    //void compute_corner_normal( msq_stdc::size_t corner_pt, 
    //                            Vector3D &normal,
    //                            PatchData &pd, 
    //                            MsqError &err);
    void compute_corner_normals( Vector3D normals[], PatchData& pd, MsqError& err );

      //! Compute matrices which column are the vectors issued from a corner.
      //! Stores those corner matrices in the mTag data member.  
    void compute_corner_matrices(PatchData &pd, Matrix3D A[], int num_m3d, MsqError &err );

  private:
    static void get_linear_quad_jac(Vector3D *sp,
                                    Vector3D &coord0, Vector3D &coord1,
                                    Vector3D &coord2, Vector3D &coord3,
                                    Vector3D* jac);
    
    EntityTopology mType;
    /** Pointer to connectivity array.
     *  NOTE: The memory occupied by this array is assumed to be
     *        owned/managed by the owning PatchData.  Do not try
     *        to delete/resize/etc. this array directly!
     */
    size_t* vertexIndices;  
    size_t numVertexIndices;
    
      // output operator for debugging.
    friend msq_stdio::ostream& operator<<(msq_stdio::ostream &stream, 
                                          const MsqMeshEntity &entity);
    
  };
  
  
    // Returns the number of vertices in this type
    // of element, or 0 if a variable number.
  inline size_t MsqMeshEntity::vertex_count() const
  { return mType == POLYGON || mType == POLYHEDRON ? 
           node_count() : TopologyInfo::corners( mType ); }

  inline void MsqMeshEntity::set_connectivity( msq_stdc::size_t *indices,
                                               size_t num_vertices )
  {
    vertexIndices = indices;
    numVertexIndices = num_vertices;
  }
  
  inline const msq_stdc::size_t *MsqMeshEntity::get_vertex_index_array() const
  { return vertexIndices; }
  
  inline msq_stdc::size_t *MsqMeshEntity::get_vertex_index_array()
  { return vertexIndices; }
 
  inline msq_stdc::size_t 
  MsqMeshEntity::get_vertex_index(msq_stdc::size_t vertex_in_element) const
  {
      // Make sure we're in range
    assert(vertex_in_element < vertex_count());
      // Return the index
    return vertexIndices[vertex_in_element];
  }

  inline void MsqMeshEntity::get_linear_quad_jac(Vector3D *sp,
                                                 Vector3D &coord0,
                                                 Vector3D &coord1,
                                                 Vector3D &coord2,
                                                 Vector3D &coord3,
                                                 Vector3D* jac)
  {
    jac[0]=coord1-coord0+(*sp)[1]*(coord2+coord0-coord3-coord1);
    jac[1]=coord3-coord0+(*sp)[0]*(coord2+coord0-coord3-coord1);
  }
} //namespace


#endif // MsqMeshEntity_hpp
