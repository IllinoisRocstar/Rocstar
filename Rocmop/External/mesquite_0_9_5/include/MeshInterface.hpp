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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu
   
  ***************************************************************** */
/*!
  \file   MeshInterface.hpp
  \brief  This file contains the Mesquite mesh interface.
          Many users will want to implement a concrete class derived from
          the MeshInterface class to access their mesh.


  \author Darryl Melander
  \author Thomas Leurent
  \date   2003-04-17
*/
#ifndef MESQUITE_INTERFACE_HPP
#define MESQUITE_INTERFACE_HPP

#include "Mesquite.hpp"
#include "TopologyInfo.hpp"

#ifdef MSQ_USE_OLD_C_HEADERS
#  include <stddef.h>
#else
#  include <cstddef>
#endif

#include <string>

namespace Mesquite
{
  class EntityIterator;
  class MsqError;
  class MsqVertex;
  class Vector3D;
  typedef EntityIterator VertexIterator;
  typedef EntityIterator ElementIterator;
    
    /** Type used to refer to a tag defintion */
  typedef void* TagHandle;

  inline size_t vertices_in_topology(EntityTopology);

  /*! \class Mesh
     \brief  A Mesquite::Mesh is a collection of mesh elements which are
     composed of mesh vertices.  Intermediate objects are not accessible
     through this interface (where intermediate objects include things
     like the faces of a hex, or an element's edges).
  */
  class Mesh
  {
  public:
//************ Type Definitions **************
      //! Opaque EntityHandle type and tag type.
    typedef void* EntityHandle;
    
    
      // We typedef specific types of EntityHandles just
      // to make it clear what kind of entity is to be
      // returned or used as a function parameter, but the
      // different handle types are not actually distinct.
    typedef EntityHandle VertexHandle;
    typedef EntityHandle ElementHandle;

//************ Operations on entire mesh ****************
      //! Returns whether this mesh lies in a 2D or 3D coordinate system.
    virtual int get_geometric_dimension(MsqError &err) = 0;
    
    /** \brief get sizes for calling \ref get_all_mesh
     *
     * Get counts of entities in mesh.
     *
     *\param vertex_count  - Number of vertices connected to active mesh
     *\param element_count - Number of elements in active mesh
     *\param vertex_use_count - Number of vertex uses (sum of the length
     *                          of the connectivity list for all elements
     *                          in active.)
     */
    virtual void get_all_sizes( size_t& vertex_count,
                                size_t& element_count,
                                size_t& vertex_use_count,
                                MsqError& err ) = 0;
    
    /** \brief Get entities and connectivity 
     *
     * Get vertex handles, element handles, and connectivty
     * for active mesh.  Use \ref get_all_sizes to determine
     * required array sizes.
     *
     *\param vert_array        Array to store vertex handles in
     *\param vert_len          Length of \ref vert_array
     *\param elem_array        Array to store element handles in
     *\param elem_len          Length of \ref elem_array
     *\param elem_conn_offsets Offsets into \ref elem_conn_indices at
     *                         which the connectivity data for each
     *                         element begins.  
     *\param offset_len        Length of \ref elem_conn_offsets.  Should
     *                         be \ref elem_len + 1.
     *\param elem_conn_indices Indices into \ref vert_array
     *\param index_len         Length of \ref elem_conn_indices.
     */
    virtual void get_all_mesh( VertexHandle*  vert_array, size_t vert_len,
                               ElementHandle* elem_array, size_t elem_len,
                               size_t* elem_conn_offsets, size_t offset_len,
                               size_t* elem_conn_indices, size_t index_len,
                               MsqError& err ) = 0;
    
      //! Returns a pointer to an iterator that iterates over the
      //! set of all vertices in this mesh.  The calling code should
      //! delete the returned iterator when it is finished with it.
      //! If vertices are added or removed from the Mesh after obtaining
      //! an iterator, the behavior of that iterator is undefined.
    virtual VertexIterator* vertex_iterator(MsqError &err) = 0;
    
      //! Returns a pointer to an iterator that iterates over the
      //! set of all top-level elements in this mesh.  The calling code should
      //! delete the returned iterator when it is finished with it.
      //! If elements are added or removed from the Mesh after obtaining
      //! an iterator, the behavior of that iterator is undefined.
    virtual ElementIterator* element_iterator(MsqError &err) = 0;

//************ Vertex Properties ********************
      //! Returns true or false, indicating whether the vertex
      //! is allowed to be repositioned.  True indicates that the vertex
      //! is fixed and cannot be moved.  Note that this is a read-only
      //! property; this flag can't be modified by users of the
      //! Mesquite::Mesh interface.
    virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err) = 0;

      //! Returns true or false, indicating whether the vertex
      //! is on the boundary.  Boundary nodes may be treated as
      //! a special case by some algorithms or culling methods.
      //! Note that this is a read-only
      //! property; this flag can't be modified by users of the
      //! Mesquite::Mesh interface.
    virtual void vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
                                  size_t num_vtx, MsqError &err) = 0;
    
      //! Get/set location of a vertex
    virtual void vertices_get_coordinates(const VertexHandle vert_array[],
                                        MsqVertex* coordinates,
                                        size_t num_vtx,
                                        MsqError &err) = 0;
    virtual void vertex_set_coordinates(VertexHandle vertex,
                                        const Vector3D &coordinates,
                                        MsqError &err) = 0;
    
      //! Each vertex has a byte-sized flag that can be used to store
      //! flags.  This byte's value is neither set nor used by the mesh
      //! implementation.  It is intended to be used by Mesquite algorithms.
      //! Until a vertex's byte has been explicitly set, its value is 0.
    virtual void vertex_set_byte (VertexHandle vertex,
                                  unsigned char byte, MsqError &err) = 0;
    virtual void vertices_set_byte (VertexHandle *vert_array,
                                    unsigned char *byte_array,
                                    size_t array_size, MsqError &err) = 0;
    
      //! Retrieve the byte value for the specified vertex or vertices.
      //! The byte value is 0 if it has not yet been set via one of the
      //! *_set_byte() functions.
    virtual void vertex_get_byte(VertexHandle vertex,
                                 unsigned char *byte, MsqError &err) = 0;
    virtual void vertices_get_byte(VertexHandle *vertex,
                                   unsigned char *byte_array,
                                   size_t array_size, MsqError &err) = 0;
    
//**************** Vertex Topology *****************    
      //! Gets the number of elements attached to this vertex.
      //! Useful to determine how large the "elem_array" parameter
      //! of the vertex_get_attached_elements() function must be.
    virtual size_t vertex_get_attached_element_count(VertexHandle vertex,
                                                     MsqError &err) = 0;
    
      //! Gets the elements attached to this vertex.
    virtual void vertex_get_attached_elements(VertexHandle vertex,
                                              ElementHandle* elem_array,
                                              size_t sizeof_elem_array,
                                              MsqError &err) = 0;
    
//*************** Element Topology *************
    
      //! Gets the number of vertices in this element.
      //! This data can also be found by querying the
      //! element's topology and getting the number
      //! of vertices per element for that topology type.
    virtual size_t element_get_attached_vertex_count(ElementHandle elem,
                                                     MsqError &err) = 0;
    virtual size_t get_vertex_use_count( ElementHandle* handle_array,
                                         size_t num_handles,
                                         MsqError& err ) = 0;

    /*! \brief  
  Returns the vertices that are part of the topological definition of each
  element in the "elem_handles" array.  
  
  When this function is called, the
  following must be true:
   -# "elem_handles" points at an array of "num_elems" element handles.
   -# "vert_handles" points at an array of size "sizeof_vert_handles"
   -# "csr_data" points at an array of size "sizeof_csr_data"
   -# "csr_offsets" points at an array of size "num_elems+1"
      
  When this function returns, adjacency information will be stored
  in csr format:
    -# "vert_handles" stores handles to all vertices found in one
       or more of the elements.  Each vertex appears only
       once in "vert_handles", even if it is in multiple elements.
    -# "sizeof_vert_handles" is set to the number of vertex
       handles placed into "vert_handles".
    -# "sizeof_csr_data" is set to the total number of vertex uses (for
       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
       the two triangles share some vertices).
    -# "csr_offsets" is filled such that csr_offset[i] indicates the location
       of entity i's first adjacency in "csr_data".  The number of vertices
       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
       reason, csr_offsets[num_elems] is set to the new value of
       "sizeof_csr_data".
    -# "csr_data" stores integer offsets which give the location of
       each adjacency in the "vert_handles" array.

  As an example of how to use this data, you can get the handle of the first
  vertex in element #3 like this:
   \code VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ] \endcode

  and the second vertex of element #3 like this:
   \code VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ] \endcode
*/
    virtual void elements_get_attached_vertices(ElementHandle *elem_handles,
                                                size_t num_elems,
                                                VertexHandle *vert_handles,
                                                size_t &sizeof_vert_handles,
                                                size_t *csr_data,
                                                size_t &sizeof_csr_data,
                                                size_t *csr_offsets,
                                                MsqError &err) = 0;
    
    
      //! Returns the topology of the given entity.
    virtual EntityTopology element_get_topology(ElementHandle entity_handle,
                                              MsqError &err) = 0;
      //! Returns the topologies of the given entities.  The "entity_topologies"
      //! array must be at least "num_elements" in size.
    virtual void elements_get_topologies(ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err) = 0;

    
//***************  Tags  ***********
    
      /** The type of a tag */
    enum TagType { BYTE, BOOL, INT, DOUBLE, HANDLE };

      /** \brief Create a tag
       *
       * Create a user-defined data type that can be attached
       * to any element or vertex in the mesh.  For an opaque or
       * undefined type, use type=BYTE and length=sizeof(..).
       *
       * \param tag_name  A unique name for the data object
       * \param type      The type of the data
       * \param length    Number of values per entity (1->scalar, >1 ->vector)
       * \param default_value Default value to assign to all entities - may be NULL
       * \return - Handle for tag definition 
       */
    virtual TagHandle tag_create( const msq_std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err) = 0;
    
      /** \brief Remove a tag and all corresponding data
       *
       * Delete a tag.
       */
    virtual void tag_delete( TagHandle handle, MsqError& err ) = 0;
    
    
      /** \brief Get handle for existing tag, by name. 
        *
        * Check for the existance of a tag given it's name and
        * if it exists return a handle for it.  If the specified
        * tag does not exist, zero should be returned WITHOUT 
        * flagging an error.
        */
    virtual TagHandle tag_get( const msq_std::string& name, 
                               MsqError& err ) = 0;
     
      /** \brief Get properites of tag
       *
       * Get data type and number of values per entity for tag.
       * \param handle     Tag to get properties of.
       * \param name_out   Passed back tag name.
       * \param type_out   Passed back tag type.
       * \param length_out Passed back number of values per entity.
       */
    virtual void tag_properties( TagHandle handle,
                                 msq_std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err ) = 0;
    
      /** \brief Set tag values on elements
       * 
       * Set the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err ) = 0;

      /** \brief Set tag values on vertices
       * 
       * Set the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of node_array
       * \param node_array Array of vertices for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err ) = 0;
    
    
      /** \brief Get tag values on elements
       * 
       * Get the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err ) = 0;
    
      /** \brief Get tag values on vertices.
       * 
       * Get the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of vertices for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err ) = 0;

    
    
//**************** Memory Management ****************
      //! Tells the mesh that the client is finished with a given
      //! entity handle.  
    virtual void release_entity_handles(EntityHandle *handle_array,
                                        size_t num_handles, MsqError &err) = 0;
    
      //! Instead of deleting a Mesh when you think you are done,
      //! call release().  In simple cases, the implementation could
      //! just call the destructor.  More sophisticated implementations
      //! may want to keep the Mesh object to live longer than Mesquite
      //! is using it.
    virtual void release() = 0;
    
  protected:
      //! Don't allow a Mesh to be deleted directly.
    virtual ~Mesh()
      {}
  };
  
  /*! \class EntityIterator 
  \brief   Iterates through a set of entities.  An EntityIterator 
  is typically obtained via Mesh::vertex_iterator() or
  Mesh::element_iterator().  An iterator obtained in this
  way iterates over the set of all vertices/elements in
  the Mesh from which the iterator was obtained. */
  class EntityIterator
  {
  public:
    virtual ~EntityIterator()
      {}

      //! Moves the iterator back to the first
      //! entity in the list.
    virtual void restart() = 0;
    
      //! *iterator.  Return the handle currently
      //! being pointed at by the iterator.
    virtual Mesh::EntityHandle operator*() const = 0;
    
      //! ++iterator
    virtual void operator++() = 0;
    
      //! Returns false until the iterator has
      //! been advanced PAST the last entity.
      //! Once is_at_end() returns true, *iterator
      //! returns 0.
    virtual bool is_at_end() const = 0;
  };
  

    /** A hint on the characteristics of the domain
      * that Mesquite may use to determine what, if
      * any, scheme to use to cache characteristics 
      * of the geometric domain.
      */ 
  enum DomainHint {
    NO_DOMAIN_HINT,     //< Do not make any assumptions about the domain.
    PLANAR_DOMAIN,      //< Domain is planar
    LOCAL_PALANAR,      //< Domain is close to planar 
    SMOOTH_DOMAIN       //< Domain is C1-continuous
 };
  /*! \class MeshDomain
      The MeshDomain class provides geometrical information concerning the Mesh.
      It is called during surface meshes optimization to figure out the surface normal,
      how to snap vertices back to the surface, etc... . 
    */
  class MeshDomain
  {
  public:
    virtual ~MeshDomain()
      {}
      
      //! Give a hint about the nature of the domain for
      //! better performance.  For implementations, if
      //! unsure, return NO_DOMAIN_HINT.  SMOOTH_DOMAIN is
      //! a good default choice if the domain is a single
      //! geometric surface.
    virtual DomainHint hint() const = 0;

      //! Modifies "coordinate" so that it lies on the
      //! domain to which "entity_handle" is constrained.
      //! The handle determines the domain.  The coordinate
      //! is the proposed new position on that domain.
    virtual void snap_to(Mesh::EntityHandle entity_handle,
                         Vector3D &coordinate) const = 0;
    
      //! Returns the normal of the domain to which
      //! "entity_handle" is constrained.  For non-planar surfaces,
      //! the normal is calculated at the point on the domain that
      //! is closest to the passed in value of "coordinate".  If the
      //! domain does not have a normal, or the normal cannot
      //! be determined, "coordinate" is set to (0,0,0).  Otherwise,
      //! "coordinate" is set to the domain's normal at the
      //! appropriate point.
      //! In summary, the handle determines the domain.  The coordinate
      //! determines the point of interest on that domain.
      //!
      //! User should see also PatchData::get_domain_normal_at_vertex and
      //! PatchData::get_domain_normal_at_element .
    virtual void normal_at(Mesh::EntityHandle entity_handle,
                           Vector3D &coordinate) const = 0;
                           
      /**\brief evaluate surface normals
       *
       * Returns normals for a domain.
       *
       *\param entity_handle The domain evaluated is the one in which
       *                     this mesh entity is constrained.
       *\param coordinates   As input, a list of positions at which to
       *                     evaluate the domain.  As output, the resulting
       *                     domain normals.
       *\param count         The length of the coordinates array.
       */
    virtual void normal_at( Mesh::EntityHandle handle,
                            Vector3D coordinates[],
                            unsigned count,
                            MsqError& err ) const = 0;
                            
      /**\brief evaluate closest point and normal
       *
       * Given a position in space, return the closest 
       * position in the domain and the domain normal
       * at that point.
       *
       *\param entity_handle Evaluate the subset of the domain contianing
       *                     this entity
       *\param position      Input position for which to evaluate
       *\param closest       Closest position in the domain.
       *\param normal        Domain normal at the location of 'closest'
       */
    virtual void closest_point( Mesh::EntityHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const = 0;
       
  };
}

inline size_t Mesquite::vertices_in_topology(Mesquite::EntityTopology topo)
{
  return TopologyInfo::corners( topo );
}

#endif
