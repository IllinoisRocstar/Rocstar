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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov  ,
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */
/*!
  \file   MeshImpl.hpp
  \brief  

  \author Darryl Melander
  \author Thomas Leurent
  \author Jason Kraftcheck
  \date   2003-04-17
*/

#ifndef MESQUITE_MESH_IMPL_HPP
#define MESQUITE_MESH_IMPL_HPP

#include "MeshInterface.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#include <map.h>
#else
#include <map>
#endif

namespace Mesquite
{
  class FileTokenizer;
  class MeshImplData;
  class MeshImplTags;
  struct TagDescription;
  
  /*!  \class MeshImpl

  \brief MeshImpl is a Mesquite implementation of the Mesh interface. 
  Applications can also provide their own implementation of the interface.
    
  MeshImpl can read in mesh files in VTK format and ExodusII format. 
  */
  class MeshImpl : public Mesquite::Mesh
  {
  public:
//********* Functions that are NOT inherited ************
    MeshImpl();
    virtual ~MeshImpl();
    void read_vtk(const char* in_filename, Mesquite::MsqError &err);
    void write_vtk(const char* out_filename, Mesquite::MsqError &err);
    void read_exodus(const char* in_filename, Mesquite::MsqError &err);
    void write_exodus(const char* out_filename, Mesquite::MsqError &err);
//********* Functions that ARE inherited ************
      // Returns whether this mesh lies in a 2D or 3D coordinate system.
    virtual int get_geometric_dimension(MsqError &err) ;
    
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
                                MsqError& err );
    
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
     *                         be t\ref elem_len + 1.
     *\param elem_conn_indices Indices into \ref vert_array
     *\param index_len         Length of \ref elem_conn_indices.
     */
    virtual void get_all_mesh( VertexHandle*  vert_array, size_t vert_len,
                               ElementHandle* elem_array, size_t elem_len,
                               size_t* elem_conn_offsets, size_t offset_len,
                               size_t* elem_conn_indices, size_t index_len,
                               MsqError& err );
    
    /** Get sum of number of vertices in each element */
    virtual size_t get_vertex_use_count( ElementHandle* elem_array,
                                         size_t elem_array_length,
                                         MsqError& err );
    
      // Returns a pointer to an iterator that iterates over the
      // set of all vertices in this mesh.  The calling code should
      // delete the returned iterator when it is finished with it.
      // If vertices are added or removed from the Mesh after obtaining
      // an iterator, the behavior of that iterator is undefined.
    virtual VertexIterator* vertex_iterator(MsqError &err);
    
      // Returns a pointer to an iterator that iterates over the
      // set of all top-level elements in this mesh.  The calling code should
      // delete the returned iterator when it is finished with it.
      // If elements are added or removed from the Mesh after obtaining
      // an iterator, the behavior of that iterator is undefined.
    virtual ElementIterator* element_iterator(MsqError &err);

//************ Vertex Properties ********************
      // Returns true or false, indicating whether the vertex
      // is allowed to be repositioned.  True indicates that the vertex
      // is fixed and cannot be moved.  Note that this is a read-only
      // property; this flag can't be modified by users of the
      // Mesquite::Mesh interface.
    virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err);

      // Returns true or false, indicating whether the vertex
      // is on the boundary.  Boundary nodes may be treated as
      // a special case by some algorithms or culling methods.
      // Note that this is a read-only
      // property; this flag can't be modified by users of the
      // Mesquite::Mesh interface.
    virtual void vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
                                  size_t num_vtx, MsqError &err);
    
      // Get/set location of a vertex
    virtual void vertices_get_coordinates(const Mesh::VertexHandle vert_array[],
                                          Mesquite::MsqVertex* coordinates,
                                          size_t num_vtx,
                                          MsqError &err);
    virtual void vertex_set_coordinates(VertexHandle vertex,
                                        const Vector3D &coordinates,
                                        MsqError &err);
    
      // Each vertex has a byte-sized flag that can be used to store
      // flags.  This byte's value is neither set nor used by the mesh
      // implementation.  It is intended to be used by Mesquite algorithms.
      // Until a vertex's byte has been explicitly set, its value is 0.
    virtual void vertex_set_byte (VertexHandle vertex,
                                  unsigned char byte,
                                  MsqError &err);
    virtual void vertices_set_byte (VertexHandle *vert_array,
                                    unsigned char *byte_array,
                                    size_t array_size,
                                    MsqError &err);
    
      // Retrieve the byte value for the specified vertex or vertices.
      // The byte value is 0 if it has not yet been set via one of the
      // *_set_byte() functions.
    virtual void vertex_get_byte(VertexHandle vertex,
                                 unsigned char *byte,
                                 MsqError &err);
    virtual void vertices_get_byte(VertexHandle *vertex,
                                   unsigned char *byte_array,
                                   size_t array_size,
                                   MsqError &err);
    
//**************** Vertex Topology *****************    
      // Gets the number of elements attached to this vertex.
      // Useful to determine how large the "elem_array" parameter
      // of the vertex_get_attached_elements() function must be.
    virtual size_t vertex_get_attached_element_count(VertexHandle vertex,
                                                     MsqError &err);
    
      // Gets the elements attached to this vertex.
    virtual void vertex_get_attached_elements(VertexHandle vertex,
                                              ElementHandle* elem_array,
                                              size_t sizeof_elem_array,
                                              MsqError &err);
    
    
//*************** Element Topology *************
    
      // Gets the number of vertices in this element.
      // This data can also be found by querying the
      // element's topology and getting the number
      // of vertices per element for that topology type.
    virtual size_t element_get_attached_vertex_count(ElementHandle elem,
                                                     MsqError &err);
    
// Returns the vertices that are part of the topological definition of each
// element in the "elem_handles" array.  When this function is called, the
// following must be true:
//   a) "elem_handles" points at an array of "num_elems" element handles.
//   b) "vert_handles" points at an array of size "sizeof_vert_handles"
//   c) "csr_data" points at an array of size "sizeof_csr_data"
//   d) "csr_offsets" points at an array of size "num_elems+1"
//      
// When this function returns, adjacency information will be stored
// in csr format:
//    a) "vert_handles" stores handles to all vertices found in one
//       or more of the elements.  Each vertex appears only
//       once in "vert_handles", even if it is in multiple elements.
//    b) "sizeof_vert_handles" is set to the number of vertex
//       handles placed into "vert_handles".
//    c) "sizeof_csr_data" is set to the total number of vertex uses (for
//       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
//       the two triangles share some vertices).
//    c) "csr_offsets" is filled such that csr_offset[i] indicates the location
//       of entity i's first adjacency in "csr_data".  The number of vertices
//       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
//       reason, csr_offsets[num_elems] is set to the new value of
//       "sizeof_csr_data".
//    d) "csr_data" stores integer offsets which give the location of
//       each adjacency in the "vert_handles" array.
//
// As an example of how to use this data, you can get the handle of the first
// vertex in element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ]
//
// and the second vertex of element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ]
// 
    virtual void elements_get_attached_vertices(ElementHandle *elem_handles,
                                                size_t num_elems,
                                                VertexHandle *vert_handles,
                                                size_t &sizeof_vert_handles,
                                                size_t *csr_data,
                                                size_t &sizeof_csr_data,
                                                size_t *csr_offsets,
                                                MsqError &err);
    
    void element_get_connectivity( ElementHandle element,
                                   VertexHandle* vert_handles,
                                   size_t sizeof_vert_handles, 
                                   MsqError& err );
    
      // Returns the topology of the given entity.
    virtual EntityTopology element_get_topology(ElementHandle entity_handle,
                                                MsqError &err);
      // Returns the topologies of the given entities.  The "entity_topologies"
      // array must be at least "num_elements" in size.
    virtual void elements_get_topologies(ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements,
                                         MsqError &err);

    
//*************** Tags  ***********

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
                                  MsqError &err);
     
      /** \brief Remove a tag and all corresponding data
       *
       * Delete a tag.
       */
    virtual void tag_delete( TagHandle handle, MsqError& err );
    
    
      /** \brief Get handle for existing tag, by name. */
    virtual TagHandle tag_get( const msq_std::string& name, 
                               MsqError& err );
     
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
                                 MsqError& err );
    
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
                                       MsqError& err );

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
                                       MsqError& err );
    
    
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
                                       MsqError& err );
    
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
                                       MsqError& err );

//**************** Memory Management ****************
      // Tells the mesh that the client is finished with a given
      // entity handle.  
    virtual void release_entity_handles(EntityHandle *handle_array,
                                        size_t num_handles,
                                        MsqError &err);
    
      // Instead of deleting a Mesh when you think you are done,
      // call release().  In simple cases, the implementation could
      // just call the destructor.  More sophisticated implementations
      // may want to keep the Mesh object to live longer than Mesquite
      // is using it.
    virtual void release();
    
      // Remove all data
    void clear();
    
  private:
    
    /** Coordinate values per vertex */
    int numCoords;
    
    MeshImplData* myMesh;
    MeshImplTags* myTags;


//**************** VTK Parsing ****************

      /** Read a data block from the file */
    void vtk_read_dataset( FileTokenizer& file, MsqError& err );
    
      /** Read structured point mesh */
    void vtk_read_structured_points( FileTokenizer& file, MsqError& err );
      /** Read structured grid mesh */
    void vtk_read_structured_grid  ( FileTokenizer& file, MsqError& err );
      /** Read rectilinear grid structured mesh */
    void vtk_read_rectilinear_grid ( FileTokenizer& file, MsqError& err );
      /** Read polydata mesh */
    void vtk_read_polydata         ( FileTokenizer& file, MsqError& err );
      /** Read unstructured mesh */
    void vtk_read_unstructured_grid( FileTokenizer& file, MsqError& err );
      /** Read file-level field data */
    void vtk_read_field            ( FileTokenizer& file, MsqError& err );
    
      /** Helper function for \ref vtk_read_polydata - reads polygon subsection */
    void vtk_read_polygons( FileTokenizer& file, MsqError& err );
      /** Helper function for readers of structured mesh - create elements */
    void vtk_create_structured_elems( const long* dims, MsqError& err );
    
      /** Read attribute data for vertices */
    void vtk_read_point_data( FileTokenizer& file, MsqError& err );
      /** Read attribute data for elements */
    void vtk_read_cell_data ( FileTokenizer& file, MsqError& err );
      /** Read actual data for both \ref vtk_read_point_data and \ref vtk_read_cell_data 
       *  Initializes all fields of passed TagDescription */
    void* vtk_read_attrib_data( FileTokenizer& file, 
                                long num_data_to_read, 
                                TagDescription& tag_out,
                                MsqError& err );
      /** Read a 2-D array of data of the specified type from the file 
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_typed_data( FileTokenizer& file, int type,
                               size_t per_elem, size_t num_elem,
                               TagDescription& tag_out,
                               MsqError& err );
                             
      /** Read scalar attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_scalar_attrib ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read color attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_color_attrib  ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read vector or normal attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_vector_attrib ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read texture attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_texture_attrib( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read tensor (3x3 matrix) data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_tensor_attrib ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );

      /** Write tag data to VTK attributes */
    void vtk_write_attrib_data( msq_stdio::ostream& file,
                                const TagDescription& desc,
                                const void* data, size_t count,
                                MsqError& err ) const;

//**************** End VTK Parsing ****************
  };
}

#endif
