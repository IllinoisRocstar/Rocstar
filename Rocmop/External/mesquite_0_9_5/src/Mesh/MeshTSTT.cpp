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
  \file   MeshTSTT.cpp
  \brief  


  \author Jason Kraftcheck
  \date   2004-10-26
*/

#include <sidl_cxx.hh>
#include "TSTTB.hh"
#include "TSTTM.hh"
#include "MeshTSTT.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "MeshInterface.hpp"
#include "TSTTUtil.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
# include <map.h>
# include <set.h>
#else
# include <map>
# include <set>
#endif

#define OPAQUE_TYPE_OPAQUE_PADDED 0
#define OPAQUE_TYPE_OPAQUE_PACKED 1
#define OPAQUE_TYPE_CHAR          2
#define OPAQUE_TYPE_UCHAR         3
#define OPAQUE_TYPE_BYTE          OPAQUE_TYPE_UCHAR
#define TSTT_OPAQUE_TAG_TYPE OPAQUE_TYPE_CHAR


namespace Mesquite
{



/*************************************************************************
 *                          Iterator Definitions
 *
 * This file contains three alternate iterator implementations.  After
 * some performance testing, two should be removed.
 ************************************************************************/

/** \brief Wrapper around single-entity TSTT interator */
class TSTTIterator : public EntityIterator
{
  private:
    TSTTM::Entity tsttMesh;        /**< TSTT mesh interface */
    void* tsttIter;                /**< TSTT iterator handle */
    bool notAtEnd;                 /**< Flag to mark if end of iterator had been reached */
    void* entityHandle;            /**< Current TSTT entity handle */

      /** 
       * Advance the iterator to the next entity.  Updates \ref entityHandle and 
       * \ref notAtEnd
       */
    inline void get_next_entity()
      { notAtEnd = tsttMesh.getNextEntIter( tsttIter, entityHandle ); }

  public:

      /**
       *\param mesh    The TSTT mesh interface instance
       *\param meshset The meshset to iterate over
       *\param type    TSTT type in \ref meshset to iterate over
       *\param topo    TSTT topology in \ref meshset to iterate over
       */
    TSTTIterator( TSTTM::Entity& mesh, 
                  void* meshset, 
                  TSTTM::EntityType type,
                  TSTTM::EntityTopology topo,
                  MsqError& err );
      /**\brief reset iterator */
    virtual void restart();
      /**\brief get current entity handle */
    virtual Mesh::EntityHandle operator*() const;
      /**\brief check if any remaining entity handles */
    virtual bool is_at_end() const;
      /**\biref step */
    virtual void operator++();
    
    virtual ~TSTTIterator();
};



TSTTIterator::TSTTIterator( TSTTM::Entity& mesh, void* meshset, 
                            TSTTM::EntityType type, TSTTM::EntityTopology topo,
                            MsqError& err )
  : tsttMesh(mesh), tsttIter(0), notAtEnd(true), entityHandle(0)
{
  try {
    tsttMesh.initEntIter( meshset, type, topo, tsttIter );
    get_next_entity();
  } 
  catch (TSTTB::Error& tstt_err) {
    MSQ_SETERR(err)( process_tstt_error( tstt_err ), MsqError::INTERNAL_ERROR );
  }
}

TSTTIterator::~TSTTIterator()
{
  if (tsttIter)
    tsttMesh.endEntIter( tsttIter );
}

void TSTTIterator::restart()
{
  try {
    tsttMesh.resetEntIter( tsttIter );
    get_next_entity();
  }
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}

Mesh::EntityHandle TSTTIterator::operator*() const
{
  return entityHandle;
}

bool TSTTIterator::is_at_end() const
{
  return !notAtEnd;
}

void TSTTIterator::operator++() 
{
  try {
    get_next_entity();
  }
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}
       

/** \brief Iterate over a sidl array of TSTT entity handles */
class SIDLIterator : public EntityIterator
{
  private:
    sidl::array<void*>::const_iterator iter;
    const sidl::array<void*> array;
  
  public:

      /**\param array Array to iterate over */
    SIDLIterator( const sidl::array<void*>& a )
      : iter( a.begin() ), array( a )            {}

      /**\brief reset iterator */
    virtual void restart()                       { iter = array.begin(); }
    
      /**\brief get current entity handle */
    virtual Mesh::EntityHandle operator*() const { return *iter; }
    
      /**\brief check if any remaining entity handles */
    virtual bool is_at_end() const               { return iter == array.end(); }
    
      /**\biref step */
    virtual void operator++()                    { ++iter; }
};

/** \brief TSTT iterator using array-iterator interface and buffer of handles */
class TSTTArrIter : public EntityIterator
{
  private:
  
    sidl::array<void*> handleArray;
    void* tsttIter;
    int index, count;
    bool notAtEnd;
    TSTTM::Arr& myMesh;
    
    inline void get_next_array() 
    {
      index = count = 0;
      notAtEnd = myMesh.getEntArrNextIter( tsttIter, handleArray, count ) && count;
    }
    
  public:
  
    TSTTArrIter( TSTTM::Arr& mesh, 
                 void* meshset, 
                 TSTTM::EntityType type,
                 TSTTM::EntityTopology topo,
                 unsigned buffer_count = 1024 );
    
      /**\brief reset iterator */
    virtual void restart();
      /**\brief get current entity handle */
    virtual Mesh::EntityHandle operator*() const;
      /**\brief check if any remaining entity handles */
    virtual bool is_at_end() const;
      /**\biref step */
    virtual void operator++();
    
    virtual ~TSTTArrIter();
};

TSTTArrIter::TSTTArrIter( TSTTM::Arr& mesh, 
                          void* meshset, 
                          TSTTM::EntityType type,
                          TSTTM::EntityTopology topo,
                          unsigned buffer_count )
                        : handleArray( alloc_sidl_vector<void*>(buffer_count) ),
                          tsttIter(0),
                          index(0),
                          count(0),
                          notAtEnd(false),
                          myMesh( mesh )
{
  mesh.initEntArrIter( meshset, type, topo, buffer_count, tsttIter );
  get_next_array();
}

TSTTArrIter::~TSTTArrIter() 
{
  try {
    if (tsttIter)
      myMesh.endEntArrIter( tsttIter );
  }
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}

void TSTTArrIter::restart()
{
  try {
    myMesh.resetEntArrIter( tsttIter );
    get_next_array();
  }
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}

Mesh::EntityHandle TSTTArrIter::operator*() const
{
  return handleArray.get(index);
}

bool TSTTArrIter::is_at_end() const
{
  return index == count && !notAtEnd;
}

void TSTTArrIter::operator++()
{
  try {
    ++index;
    if (index == count && notAtEnd)
      get_next_array();
  }                        
  catch (TSTTB::Error& tstt_err) {
    process_tstt_error( tstt_err );
    throw;
  }
}
     
      



/*************************************************************************
 *                            Mesh Interface 
 ************************************************************************/

  /** \brief Implementation of MeshTSTT
   */
  class MeshTSTTImpl : public MeshTSTT
  {
  public:

    MeshTSTTImpl(TSTTM::Mesh& tstt_mesh, MsqError &err);
    virtual ~MeshTSTTImpl();
    
      /** \brief set mesh to be smoothed.
       *
       * Set the mesh which Mesquite is to smooth.  Optionally
       * specify fixed vertices.
       * NOTE: If an active set is not specified, the default
       *       is to use the global set (the ENTIRE mesh.)
       *
       *\param element_set TSTT entity set handle for set containing
       *                  mesh elements and vertices for which quality 
       *                  is to be improved.
       */
    virtual void set_active_set( void* element_set, MsqError& );
    

      /**\brief Get dimension of vertex coordinates (2D vs. 3D). */
    virtual int get_geometric_dimension(Mesquite::MsqError &/*err*/);
    
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
     *                         be \ref elem_len + 1.
     *\param elem_conn_indices Indices into \ref vert_array
     *\param index_len         Length of \ref elem_conn_indices.
     */
    virtual void get_all_mesh( VertexHandle*  vert_array, size_t vert_len,
                               ElementHandle* elem_array, size_t elem_len,
                               size_t* elem_conn_offsets, size_t offset_len,
                               size_t* elem_conn_indices, size_t index_len,
                               MsqError& err );
    
      /**\brief Create iterator for vertices in active set */
    virtual VertexIterator* vertex_iterator(MsqError &err);
    
      /**\brief Create iterator for elements in active set */
    virtual ElementIterator* element_iterator(MsqError &err);

      /**\brief Query "fixed" flag for a vertex */
    virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err);

      /**\brief Query "boundary" flag for an array of vertices */
    virtual void vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
                                 size_t num_vtx, MsqError &err);
    
      /**\brief Get vertex coordinates */
    virtual void vertices_get_coordinates(const VertexHandle vert_array[],
                                  MsqVertex* coordinates,
                                  size_t num_vtx, MsqError &err);
      /**\brief Set vertex coordinates */
    virtual void vertex_set_coordinates(VertexHandle vertex,
                                 const Vector3D &coordinates, MsqError &err);
    
      /**\brief Set vertex mark */
    virtual void vertex_set_byte (VertexHandle vertex,
                            unsigned char byte, MsqError &err);
      /**\brief Set vertex mark */
    virtual void vertices_set_byte (VertexHandle *vert_array,
                              unsigned char *byte_array,
                              size_t array_size, MsqError &err);
    
      /**\brief Get vertex mark */
    virtual void vertex_get_byte(VertexHandle vertex,
                                 unsigned char *byte, MsqError &err);
      /**\brief Get vertex mark */
    virtual void vertices_get_byte(VertexHandle *vert_array,
                                   unsigned char *byte_array,
                                   size_t array_size, MsqError &err);
    
      /**\brief Get vertex adjacencies */
    virtual size_t vertex_get_attached_element_count(VertexHandle vertex, MsqError &err);
    
      /**\brief Get vertex adjacencies */
    virtual void vertex_get_attached_elements(VertexHandle vertex,
                                              ElementHandle* elem_array,
                                              size_t sizeof_elem_array,
                                              MsqError &err);
    
      /**\biref Get length of connectivity list */
    virtual size_t element_get_attached_vertex_count(ElementHandle elem,
                                                  MsqError &err);
    
    virtual size_t get_vertex_use_count( ElementHandle* array,
                                         size_t length, 
                                         MsqError& err );
    
/**\brief Get element connectivity in overly-complex CSR rep.
 *
 * Returns the vertices that are part of the topological definition of each
 * element in the "elem_handles" array.  When this function is called, the
 * following must be true:
 *   a) "elem_handles" points at an array of "num_elems" element handles.
 *   b) "vert_handles" points at an array of size "sizeof_vert_handles"
 *   c) "csr_data" points at an array of size "sizeof_csr_data"
 *   d) "csr_offsets" points at an array of size "num_elems+1"
 *      
 * When this function returns, adjacency information will be stored
 * in csr format:
 *    a) "vert_handles" stores handles to all vertices found in one
 *       or more of the elements.  Each vertex appears only
 *       once in "vert_handles", even if it is in multiple elements.
 *    b) "sizeof_vert_handles" is set to the number of vertex
 *       handles placed into "vert_handles".
 *    c) "sizeof_csr_data" is set to the total number of vertex uses (for
 *       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
 *       the two triangles share some vertices).
 *    c) "csr_offsets" is filled such that csr_offset[i] indicates the location
 *       of entity i's first adjacency in "csr_data".  The number of vertices
 *       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
 *       reason, csr_offsets[num_elems] is set to the new value of
 *       "sizeof_csr_data".
 *    d) "csr_data" stores integer offsets which give the location of
 *       each adjacency in the "vert_handles" array.
 *
 * As an example of how to use this data, you can get the handle of the first
 * vertex in element #3 like this:
 *   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ]
 *
 * and the second vertex of element #3 like this:
 *   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ]
 */ 
    virtual void elements_get_attached_vertices(ElementHandle *elem_handles,
                                                size_t num_elems,
                                                VertexHandle *vert_handles,
                                                size_t &sizeof_vert_handles,
                                                size_t *csr_data,
                                                size_t &sizeof_csr_data,
                                                size_t *csr_offsets,
                                                MsqError &err);
    
  
      /**\brief Return topology type enum for specified element */
    virtual EntityTopology element_get_topology(ElementHandle entity_handle,
                                           MsqError &err);
      /**\brief Return topology type enum for an array of elements */
    virtual void elements_get_topologies(ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err);
    
//**************** Memory Management ****************
      /**\brief no-op */ 
    virtual void release_entity_handles(EntityHandle *handle_array,
                                        size_t num_handles, MsqError &err);
    
      // Instead of deleting a Mesh when you think you are done,
      // call release().  In simple cases, the implementation could
      // just call the destructor.  More sophisticated implementations
      // may want to keep the Mesh object to live longer than Mesquite
      // is using it.
    virtual void release();

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
                                       
  protected:
        
    void set_int_tag( void* tag, void* meshset, int value, MsqError& err );

    /** Populate \ref inputElements from \ref elemetnSet */
    void popupate_input_elements( ) throw( TSTTB::Error );

  private:
      /** \brief Set tag values */
    void tag_set_data ( TagHandle handle,
                        size_t num_elems,
                        const EntityHandle* handle_array,
                        const void* tag_data,
                        MsqError& err );
    
      /** \brief Get tag values */
    void tag_get_data( TagHandle handle,
                       size_t num_elems,
                       const EntityHandle* handle_array,
                       void* tag_data,
                       MsqError& err );
    
    /** TSTT basic mesh interface instance */
    TSTTM::Mesh meshIFace;
    /** TSTT interface for per-entity queries */
    TSTTM::Entity entIFace;
    TSTTB::EntTag tagIFace;
    /** TSTT interface for multi-entity (array) queries */
    TSTTM::Arr arrIFace;
    TSTTB::ArrTag arrTagIFace;
    /** TSTT interface for modifying mesh */
    TSTTM::Modify modIFace;
    /** TSTT interface for entity set operations */
    TSTTB::EntSet setIFace;
    
    /** Have mesh */
    bool haveMesh;
    /** TSTTM entity set handle for elements to improve */
    void* elementSet;
    /** TSTTM entity set handle for nodes to move */
    void* nodeSet;
    /** std::set containing elements in elementSet, used
     *  to constrain vertex->element adjaceny queries to
     *  only those elements that are in the input element set.
     */
    msq_std::set<void*> inputElements;
    
    /** The type of elements contained in the input element set.
     * Should be one of:
     * - TSTTM::EntityType_REGION    - volume elements
     * - TSTTM::EntityType_FACE      - face/2d elements
     * - TSTTM::EntityType_ALL_TYPES - mixed volume and face elements
     */
    TSTTM::EntityType inputSetType;
    
    /** Handle for tag used to hold vertex byte */
    TagHandle byteTag; 
    /** Tag was created in constructor */
    bool createdByteTag;
    /** Handle for tag used to hold vertex-fixed flag */
    TagHandle fixedTag;
    /** Fixed tag was created in constructor */
    bool createdFixedTag;
    /** Handle for the tag used internally to remove duplicates from lists */
//    TagHandle vertexIndexTag;
    /** vertexIndexTag was created in constructor */
//    bool createdVertexIndexTag;
    
    /** Map TSTTM::EntityTopology to Mesquite::EntityTopology */
    EntityTopology topologyMap[TSTTM::EntityTopology_ALL_TOPOLOGIES+1];
    
    /** Cached result for vertex->element query */
    sidl::array<void*> vertexAdjElements;
    /** Number of valid entries \ref vertexAdjElements */
    int vertexAdjElementSize;
    /** Vertex for which \ref vertexAdjElements is cached */
    void* cachedAdjVertex;
    /** Get elements adjacent to vertex and store in \ref vertexAdjElements */
    void cache_adjacent_elements( void* vertex_handle, MsqError& err );
 };

/*************************************************************************
 *                          Mesh Definition
 ************************************************************************/

MeshTSTT* MeshTSTT::create( TSTTM::Mesh& mesh, void* meshset, MsqError& err )
{
  MeshTSTT* result = new MeshTSTTImpl( mesh, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  result->set_active_set( meshset, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  return result;
}

MeshTSTT* MeshTSTT::create( TSTTM::Mesh& mesh, MsqError& err )
{
  MeshTSTT* result = new MeshTSTTImpl( mesh, err );
  if (MSQ_CHKERR(err))
  {
    delete result;
    return 0;
  }
  return result;
}

MeshTSTT::~MeshTSTT() {}


MeshTSTTImpl::MeshTSTTImpl(TSTTM::Mesh& tstt_mesh, Mesquite::MsqError& err) 
  : meshIFace(tstt_mesh), 
    elementSet(0), nodeSet(0), 
    inputSetType( TSTTM::EntityType_ALL_TYPES ),
    byteTag(0), createdByteTag(false),
    fixedTag(0), createdFixedTag(false),
    vertexAdjElementSize(0), cachedAdjVertex(0)
//    ,vertexIndexTag(0), createdVertexIndexTag(false)
{
    // Initialize topology map 
  
  const size_t mapsize = sizeof(topologyMap) / sizeof(Mesquite::EntityTopology);
  if (mapsize < TSTTM::EntityTopology_ALL_TOPOLOGIES)
  {
    MSQ_SETERR(err)("MeshTSTT needs to be updated for new TSTT element topologies.",
                    MsqError::INTERNAL_ERROR);
    return;
  }
  
  for (size_t i = 0; i <= TSTTM::EntityTopology_ALL_TOPOLOGIES; ++i)
    topologyMap[i] = Mesquite::MIXED;
  
  topologyMap[TSTTM::EntityTopology_TRIANGLE     ] = Mesquite::TRIANGLE;
  topologyMap[TSTTM::EntityTopology_QUADRILATERAL] = Mesquite::QUADRILATERAL;
  topologyMap[TSTTM::EntityTopology_TETRAHEDRON  ] = Mesquite::TETRAHEDRON;
  topologyMap[TSTTM::EntityTopology_HEXAHEDRON   ] = Mesquite::HEXAHEDRON;
  topologyMap[TSTTM::EntityTopology_PRISM        ] = Mesquite::PRISM;
  topologyMap[TSTTM::EntityTopology_PYRAMID      ] = Mesquite::PYRAMID;
  
  try {
      // Attempt to cast to single-entity query interface
    entIFace  = tstt_mesh;
    if (!entIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::EntTag",
                       MsqError::INVALID_STATE );
      return;
    } 
      // Attempt cast to multi-entity query interface
    arrIFace  = tstt_mesh;
    if (!arrIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::ArrTag",
                       MsqError::INVALID_STATE );
      return;
    }
      // Attempt cast to modify interface
    modIFace  = tstt_mesh;
    if (!arrIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::Modify",
                       MsqError::INVALID_STATE );
      return;
    }
      // Attempt cast to entity set interface
    setIFace  = tstt_mesh;
    if (!arrIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTM::EntSet",
                       MsqError::INVALID_STATE );
      return;
    }
      // Get tag interface
    tagIFace = tstt_mesh;
    if (!tagIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTB::EntTag",
                       MsqError::INVALID_STATE );
      return;
    }
      // Get tag interface
    arrTagIFace = tstt_mesh;
    if (!arrTagIFace)
    {
      MSQ_SETERR(err)( "TSTTM::Mesh does not implement TSTTB::ArrTag",
                       MsqError::INVALID_STATE );
      return;
    }
    
      // Get tag for fixed flag
    try {
      fixedTag = tagIFace.getTagHandle( VERTEX_FIXED_TAG_NAME );

        // Double-check types incase tag already existed
      if (tagIFace.getTagSize(fixedTag) != sizeof(int) || 
          tagIFace.getTagType(fixedTag) != TSTTB::TagValueType_INTEGER)
      {
        MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                         "Tag \"%s\" exists with invalid type/size", 
                         VERTEX_FIXED_TAG_NAME );
        return;
      }
    } catch (...) { 
      fixedTag = 0;
    }
    
      // Get/create tag for vertex byte
    try {
      byteTag = tagIFace.getTagHandle( VERTEX_BYTE_TAG_NAME );
    } catch (...) {}

    if (!byteTag) {
      tagIFace.createTag( VERTEX_BYTE_TAG_NAME, sizeof(int), TSTTB::TagValueType_INTEGER, byteTag );
      createdByteTag = true;
    } 
      // Double-check types incase tag already existed
    if (tagIFace.getTagSize(byteTag) != sizeof(int) || 
        tagIFace.getTagType(byteTag) != TSTTB::TagValueType_INTEGER)
    {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                       "Tag \"%s\" exists with invalid type/size", 
                       VERTEX_BYTE_TAG_NAME );
      return;
    }
    
      // Clear vertex byte tag
    set_int_tag( byteTag, nodeSet, 0, err );
    if (err)
      return;
    
      // Get/create tag for vertex index
/*
    //try {
    //  vertexIndexTag = tagIFace.getTagHandle( VERTEX_INDEX_TAG_NAME );
    //} catch (...) {}

    if (!vertexIndexTag) {
      tagIFace.createTag( VERTEX_INDEX_TAG_NAME, sizeof(int), TSTTB::TagValueType_INTEGER, byteTag );
      createdVertexIndexTag = true;
    } 
      // Double-check types incase tag already existed
    if (tagIFace.getTagSize(vertexIndexTag) != sizeof(int) || 
        tagIFace.getTagType(vertexIndexTag) != TSTTB::TagValueType_INTEGER)
    {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
                       "Tag \"%s\" exists with invalid type/size", 
                       VERTEX_INDEX_TAG_NAME );
      return;
    }
*/
  }
  catch (TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    try {
      if (createdFixedTag)
        tagIFace.destroyTag( fixedTag, false );
      if (createdByteTag)
        tagIFace.destroyTag( byteTag, false );
//      if (createdVertexIndexTag)
//        tagIFace.destroyTag( vertexIndexTag, false );
    } 
    catch (...) {}
  }
  catch(...) {
    MSQ_SETERR(err)("Uknown exception",MsqError::INTERNAL_ERROR);
  }
}

Mesquite::MeshTSTTImpl::~MeshTSTTImpl() 
{
  try {
    if (elementSet) 
      setIFace.destroyEntSet( elementSet );
    if (nodeSet)
      setIFace.destroyEntSet( nodeSet );
      
    if (createdFixedTag)
      tagIFace.destroyTag( fixedTag, false );
    if (createdByteTag)
      tagIFace.destroyTag( byteTag, false );
//    if (createdVertexIndexTag)
//      tagIFace.destroyTag( vertexIndexTag, false );
  } 
  catch (TSTTB::Error& tstt_err ) {
    process_tstt_error(tstt_err);
    throw;
  }
}


void MeshTSTTImpl::set_int_tag( void* tag,
                                void* elem_set, 
                                int value, 
                                MsqError& err )
{
  const unsigned BUFFER_COUNT = 1024;
  
  sidl::array<int> value_array( alloc_sidl_vector<int>( BUFFER_COUNT, value ) );
  
  sidl::array<void*> handle_array( alloc_sidl_vector<void*>( BUFFER_COUNT ) );
  int count = BUFFER_COUNT;
  void* iter = 0;
  bool more;
  
  try {
    
    arrIFace.initEntArrIter( elem_set, 
                             TSTTM::EntityType_VERTEX, 
                             TSTTM::EntityTopology_POINT, 
                             BUFFER_COUNT, iter );
                             
    do {
      more = arrIFace.getEntArrNextIter( iter, handle_array, count );
      if (count > 0)
        arrTagIFace.setIntArrData( handle_array, count, tag, value_array, count );
    } while (more);
    
    arrIFace.endEntArrIter( iter );
  }
  
  catch (TSTTB::Error &tstt_err) {
    if( iter ) try {  arrIFace.endEntArrIter( iter ); } catch (...) {}
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


void MeshTSTTImpl::set_active_set( void* elem_set, MsqError& err )
{
  const int ELEM_BUFFER_SIZE = 1024;
  const int NODE_BUFFER_SIZE = 27 * ELEM_BUFFER_SIZE; 
  sidl::array<void*> elements( alloc_sidl_vector<void*>( ELEM_BUFFER_SIZE ) );
  sidl::array<void*>    nodes( alloc_sidl_vector<void*>( NODE_BUFFER_SIZE ) );
  sidl::array<int>    offsets( alloc_sidl_vector<int  >( ELEM_BUFFER_SIZE ) );
  void* iter = 0;
 
  try {
    
      // Create new sets for nodes and elements
    if (elementSet)
    {
      setIFace.destroyEntSet( elementSet );
      elementSet = 0;
    }
    if (nodeSet)
    {
      setIFace.destroyEntSet( nodeSet );
      nodeSet = 0;
    }
    setIFace.createEntSet( false, elementSet );
    setIFace.createEntSet( false, nodeSet );
    
      // Iterate over set twice, once for FACEs and once for REGIONs
    bool have_faces = false, have_regions = false;
    for (int i = 0; i < 2; ++i)
    {
      TSTTM::EntityType type = i ? TSTTM::EntityType_REGION : 
                                   TSTTM::EntityType_FACE;
      bool& have_some = i ? have_regions : have_faces;
                                   
      arrIFace.initEntArrIter( elem_set, 
                               type, 
                               TSTTM::EntityTopology_ALL_TOPOLOGIES,
                               ELEM_BUFFER_SIZE, 
                               iter );
      
      int count = 0;
      bool more = false;
      do {
          // Add elements to element set
        more = arrIFace.getEntArrNextIter( iter, elements, count );
        if (!count) break;
        setIFace.addEntArrToSet( elements, count, elementSet );
        
          // Add nodes to node set
        int num_nodes, num_offsets;
        arrIFace.getEntArrAdj( elements, count, TSTTM::EntityType_VERTEX,
                               nodes, num_nodes, offsets, num_offsets );
        setIFace.addEntArrToSet( nodes, num_nodes, nodeSet );
        
        have_some = true;
      } while (more);
        
      arrIFace.endEntArrIter( iter );
      iter = 0;
      
    } // for (type)
    
    if (!have_faces)
      inputSetType = TSTTM::EntityType_REGION;
    else if (!have_regions)
      inputSetType = TSTTM::EntityType_FACE;
    else
      inputSetType = TSTTM::EntityType_ALL_TYPES;
    
    //set_fixed_tag( nodeSet, false, err ); MSQ_ERRRTN(err);
  }
  catch (TSTTB::Error &tstt_err) {
      // If an error occured, try to clean up anything created above

    try {
      if (iter) {
        arrIFace.endEntArrIter( iter );
        iter = 0;
      }
    }
    catch( ... ) {}
    
    try {
      if (elementSet) {
        setIFace.destroyEntSet( elementSet );
        elementSet = 0;
      }
    }
    catch( ... ) {}
    
    try {
      if (nodeSet) {
        setIFace.destroyEntSet( nodeSet );
        nodeSet = 0;
      }
    }
    catch( ... ) {}
    
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
  
    // clear cached data
  inputElements.clear();
  vertexAdjElementSize = 0;
  cachedAdjVertex = 0;
}

void MeshTSTTImpl::popupate_input_elements( ) throw( TSTTB::Error )
{
  const int ELEM_BUFFER_SIZE = 1024;
  void* handle_array[ELEM_BUFFER_SIZE];
  sidl::array<void*> handles = convert_to_sidl_vector( handle_array, ELEM_BUFFER_SIZE );
  void* iter;
 
  inputElements.clear();
 
  arrIFace.initEntArrIter( elementSet, 
                           TSTTM::EntityType_ALL_TYPES, 
                           TSTTM::EntityTopology_ALL_TOPOLOGIES,
                           ELEM_BUFFER_SIZE, 
                           iter );
      
  int count = 0;
  bool more = false;
  do {
      // Add elements to element set
    more = arrIFace.getEntArrNextIter( iter, handles, count );
    if (!count) break;
    setIFace.addEntArrToSet( handles, count, elementSet );

    void** const end_ptr = handle_array + count;
    for (void** ptr = handle_array; ptr != end_ptr; ++ptr)
      inputElements.insert( *ptr );

  } while (more);

  arrIFace.endEntArrIter( iter );
} 

  

// Returns whether this mesh lies in a 2D or 3D coordinate system.
int MeshTSTTImpl::get_geometric_dimension(Mesquite::MsqError &err)
{
  try {
    return meshIFace.getGeometricDim();
  }
  catch (TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}
    
    
// Returns a pointer to an iterator that iterates over the
// set of all vertices in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If vertices are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
VertexIterator* MeshTSTTImpl::vertex_iterator(MsqError& err)
{
  try {
    return new TSTTArrIter( arrIFace, 
                            nodeSet, 
                            TSTTM::EntityType_ALL_TYPES,
                            TSTTM::EntityTopology_ALL_TOPOLOGIES );
  }
  catch (TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}
    
// Returns a pointer to an iterator that iterates over the
// set of all top-level elements in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If elements are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
ElementIterator* MeshTSTTImpl::element_iterator(MsqError &err)
{
  try {
    return new TSTTArrIter( arrIFace, 
                            elementSet, 
                            TSTTM::EntityType_ALL_TYPES, 
                            TSTTM::EntityTopology_ALL_TOPOLOGIES );
  }
  catch (TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

//************ Vertex Properties ********************
// Returns true or false, indicating whether the vertex
// is allowed to be repositioned.  True indicates that the vertex
// is fixed and cannot be moved.  Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
bool MeshTSTTImpl::vertex_is_fixed( VertexHandle vertex, MsqError &err )
{
    // If mesh does not contain a fixed tag, assume no vertices are fixed
  if (!fixedTag) 
    return false;
    
  try {
    return (bool)tagIFace.getIntData( vertex, fixedTag );
  }
  catch( TSTTB::Error& tstt_err ) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return true;
  }
}

// Returns true or false, indicating whether the vertex
// is on the boundary.  Boundary nodes may be treated as
// a special case by some algorithms or culling methods.
// Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
void MeshTSTTImpl::vertices_are_on_boundary(
  VertexHandle vert_array[], 
  bool bool_array[],
  size_t num_vtx, MsqError &err)
{
    // If mesh does not contain a fixed tag, assume no vertices are fixed
  if (!fixedTag) {
    memset( bool_array, 0, num_vtx * sizeof(bool) );
    return;
  }
  
    // Get per-vertex flags from fixedTag
  try {
    sidl::array<void*> vert_wrapper( convert_to_sidl_vector( vert_array, num_vtx ) );
    sidl::array<int> bools = alloc_sidl_vector<int>(num_vtx);
    
    int num_bools = num_vtx;
    arrTagIFace.getIntArrData( vert_wrapper, num_vtx, fixedTag, bools, num_bools );
    if (num_vtx != (unsigned)num_bools)
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);

    for (size_t i = 0; i < num_vtx; ++i)
      bool_array[i] = bools.get(i);
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

// Get vertex coordinates 
void MeshTSTTImpl::vertices_get_coordinates(
  const Mesquite::Mesh::VertexHandle vert_array[],
  MsqVertex* coordinates, 
  size_t num_vtx, 
  MsqError &err)
{
  double* dbl_array = 0;
  
  try {
    int dim = meshIFace.getGeometricDim();
    int dbl_len = dim * num_vtx;

    sidl::array<void*> vertex_wrapper( 
      convert_to_sidl_vector( const_cast<void**>(vert_array), num_vtx ) );
    dbl_array = new double[dbl_len];
    sidl::array<double> dbl_wrapper( convert_to_sidl_vector( dbl_array, dbl_len ));
    
    TSTTM::StorageOrder order = TSTTM::StorageOrder_UNDETERMINED;
    meshIFace.getVtxArrCoords( vertex_wrapper, num_vtx, order, dbl_wrapper, dbl_len );
    if ((unsigned)dbl_len != dim*num_vtx) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      delete [] dbl_array;
      return;
    }
    
    if (dim == 2)
    {
      if (order == TSTTM::StorageOrder_INTERLEAVED)
      {
        double* iter = dbl_array;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].x(*iter); ++iter;
          coordinates[i].y(*iter); ++iter;
          coordinates[i].z(0);
        }
      }
      else
      {
        double *xiter = dbl_array;
        double *yiter = dbl_array + num_vtx;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].x(*xiter); ++xiter;
          coordinates[i].y(*yiter); ++yiter;
          coordinates[i].z(0);
        }
      }
    }
    else if(dim == 3)
    {
      if (order == TSTTM::StorageOrder_INTERLEAVED)
      {
        double* iter = dbl_array;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].set(iter);
          iter += 3;
        }
      }
      else
      {
        double *xiter = dbl_array;
        double *yiter = dbl_array + num_vtx;
        double *ziter = dbl_array + 2*num_vtx;
        for (size_t i = 0; i < num_vtx; ++i)
        {
          coordinates[i].x(*xiter); ++xiter;
          coordinates[i].y(*yiter); ++yiter;
          coordinates[i].z(*ziter); ++ziter;
        }
      }
    }
    else
    {
      MSQ_SETERR(err)(MsqError::INVALID_STATE, "TSTTB database dimension = %d", dim);
      delete [] dbl_array;
      return;
    }
      
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
  
  delete [] dbl_array;
}

void MeshTSTTImpl::vertex_set_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  const Vector3D &coordinates, MsqError &err)
{
  try {
    sidl::array<double> coords( 
      convert_to_sidl_vector( const_cast<double*>(coordinates.to_array()), 3 ) );
    modIFace.setVtxCoords( vertex, coords, 1 );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

// Each vertex has a byte-sized flag that can be used to store
// flags.  This byte's value is neither set nor used by the mesh
// implementation.  It is intended to be used by Mesquite algorithms.
// Until a vertex's byte has been explicitly set, its value is 0.
void MeshTSTTImpl::vertex_set_byte (
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char byte, MsqError &err)
{
  try {
    int value = byte;
    tagIFace.setIntData( vertex, byteTag, value );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::vertices_set_byte (
  VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
    // TSTT implementations seem to be inconsistant with
    // regard to setting opaque tags.  Set one at a time
    // to be safe.  This would be much easier if Mesquite
    // used a TSTT-defined type for the data, rather than
    // a single byte.
  try {
    sidl::array<void*> handles( convert_to_sidl_vector( vert_array, array_size ));
    sidl::array<int> data(alloc_sidl_vector<int>(array_size));
    for (size_t i = 0; i < array_size; ++i)
      data.set( i, byte_array[i] );
    arrTagIFace.setIntArrData( handles, array_size, byteTag, data, array_size );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

// Retrieve the byte value for the specified vertex or vertices.
// The byte value is 0 if it has not yet been set via one of the
// *_set_byte() functions.
void MeshTSTTImpl::vertex_get_byte(
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char *byte, MsqError &err)
{
  try {
    *byte = (unsigned char)tagIFace.getIntData( vertex, byteTag );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::vertices_get_byte(
  VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
    // TSTT implementations seem to be inconsistant with
    // regard to setting opaque tags.  Set one at a time
    // to be safe.  This would be much easier if Mesquite
    // used a TSTT-defined type for the data, rather than
    // a single byte.
  try {
    sidl::array<void*> handles( convert_to_sidl_vector( vert_array, array_size ));
    sidl::array<int> data( alloc_sidl_vector<int>(array_size) );
    int32_t junk;
    arrTagIFace.getIntArrData( handles, array_size, byteTag, data, junk );

    for (size_t i = 0; i< array_size; ++i )
      byte_array[i] = (unsigned char)data.get(i);
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


//**************** Vertex Topology *****************

void MeshTSTTImpl::cache_adjacent_elements( void* vtx, MsqError& err )
{
    // If already have data cached for this vertex, just return
  if (vtx == cachedAdjVertex)
    return;
  
    // make sure to clear this so that if we fail in the middle,
    // we don't end up with invlaid cache data for the prev vertex.
  cachedAdjVertex = 0;
  vertexAdjElementSize = 0;
  
  try {
      // First try using current array size
    bool success;
    try {
      entIFace.getEntAdj( vtx, inputSetType, vertexAdjElements, vertexAdjElementSize );
      success = true;
    } 
    catch( ... ) {
      success = false;
    }
  
    // If failed, try again and let implementation allocate storage
    if (!success) {
      vertexAdjElements = sidl::array<void*>();
      entIFace.getEntAdj( vtx, inputSetType, vertexAdjElements, vertexAdjElementSize );
    }
    
      // Store which vertex we have cached adj data for
    cachedAdjVertex = vtx;

      // We need to use inputElements.  Fill it if it hasn't
      // been filled yet.
    if (inputElements.empty())
      popupate_input_elements();
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
  
    // Remove from list any elements not in the input set
  void** write_ptr = convert_from_sidl_vector( vertexAdjElements );
  void** const end_ptr = write_ptr + vertexAdjElementSize;
  for (void** read_ptr = write_ptr; read_ptr != end_ptr; ++read_ptr)
  {
    if (inputElements.find( *read_ptr ) != inputElements.end())
    {
      *write_ptr = *read_ptr;
      ++write_ptr;
    }
  }
  
  vertexAdjElementSize = write_ptr - convert_from_sidl_vector( vertexAdjElements );
}
    
    


// Gets the number of elements attached to this vertex.
// Useful to determine how large the "elem_array" parameter
// of the vertex_get_attached_elements() function must be.
size_t MeshTSTTImpl::vertex_get_attached_element_count(
  VertexHandle vertex, MsqError &err)
{
  cache_adjacent_elements( vertex, err ); MSQ_ERRZERO(err);
  return vertexAdjElementSize;
}

// Gets the elements attached to this vertex.
void MeshTSTTImpl::vertex_get_attached_elements(
  VertexHandle vertex,
  ElementHandle* elem_array,
  size_t sizeof_elem_array, MsqError &err)
{
  cache_adjacent_elements( vertex, err ); MSQ_ERRRTN(err);
  if (sizeof_elem_array < (unsigned)vertexAdjElementSize) {
    MSQ_SETERR(err)("Insufficient space in array.", MsqError::INVALID_ARG);
    return;
  }
  
  assert( sizeof(ElementHandle) == sizeof(void*) );
  void** array = convert_from_sidl_vector( vertexAdjElements );
  memcpy( elem_array, array, vertexAdjElementSize * sizeof(void*) );
}


// Gets the number of vertices in this element.
// This data can also be found by querying the
// element's topology and getting the number
// of vertices per element for that topology type.
size_t MeshTSTTImpl::element_get_attached_vertex_count(
  ElementHandle elem,
  MsqError &err)
{
  try {
    sidl::array<void*> junk;
    int result = 0;
    entIFace.getEntAdj( elem, TSTTM::EntityType_VERTEX, junk, result );
    return result;
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

size_t MeshTSTTImpl::get_vertex_use_count( ElementHandle* array,
                                           size_t length, 
                                           MsqError& err )
{
  try {
    sidl::array<void*> handles( convert_to_sidl_vector( array, length ) ); 
    sidl::array<void*> adj;
    sidl::array<int> offsets;
    int32_t adj_size, off_size;
    arrIFace.getEntArrAdj( handles, length, TSTTM::EntityType_VERTEX, adj, adj_size, offsets, off_size );
    return adj_size;
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

/** Get connectivity
 *\param elements - Array of length num_elems containing elements
 *                  handles of elements for which connectivity is to
 *                  be queried.
 *\param vertices - Array of vertex handles in connectivity list.
 *\param offsets  - Indices into \ref indices array, one per element.
 *\param indices  - Indices into \ref vertex_handles
 */
void MeshTSTTImpl::elements_get_attached_vertices(
  ElementHandle *elements,
  size_t num_elems,
  VertexHandle *vertices,
  size_t &vert_len,
  size_t *indices,
  size_t &indices_len,
  size_t *offsets,
  Mesquite::MsqError &err)
{
  if (num_elems == 0)
    return;

  try {
      // Constuct arguments to TSTT
    int vert_count, off_count;
    sidl::array<void*> elem_arr( convert_to_sidl_vector( elements, num_elems ) );
    sidl::array<void*> vert_arr( convert_to_sidl_vector( (void**)indices, indices_len ) );
    sidl::array<int> off_arr;
    if (sizeof(size_t) == sizeof(int))
      off_arr = convert_to_sidl_vector( (int*)offsets, num_elems + 1 );
    else
      off_arr = alloc_sidl_vector<int>( num_elems + 1 );
    
      // Query TSTT for element connectivity
    arrIFace.getEntArrAdj( elem_arr, num_elems,
                           TSTTM::EntityType_VERTEX,
                           vert_arr, vert_count,
                           off_arr,  off_count );
  
    indices_len = vert_count;
      // If couldn't do array borrow, copy offsets into output array
    if (sizeof(size_t) != sizeof(int))
      copy_from_sidl( off_arr, offsets );
      
      // Mesquite expects index array length as last value in offset array.
      // If TSTT mesh doesn't return it, then add it.
    if ((unsigned)off_count == num_elems)
      offsets[num_elems] = vert_count;
    else
      assert( (unsigned)off_count == num_elems+1 );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
    
    // Construct unique list of vertex handles
  vert_len = 0;
  msq_std::map<VertexHandle,size_t> index_map;
  const size_t* end = indices + indices_len;
  /*
  for (size_t* iter = indices; iter != end; ++iter)
  {
    VertexHandle vertex = *(VertexHandle*)iter;
    msq_std::map<VertexHandle,size_t>::iterator pos = index_map.find( vertex );
    size_t index;
    if (pos == index_map.end())
      index_map[vertex] = index = vert_len++;
    else
      index = pos->second;
    
    vertices[index] = vertex;
    *iter = index;
  }
  */
    // Rather than the above simple code, do the following more complicated 
    // code so the resulting list of vertices are sorted (easier to debug)
  for (size_t* iter = indices; iter != end; ++iter)
  {
    VertexHandle vertex = *(VertexHandle*)iter;
    msq_std::map<VertexHandle,size_t>::iterator pos = index_map.find( vertex );
    if (pos == index_map.end())
      index_map[vertex] = 0;
  }
  for (msq_std::map<VertexHandle,size_t>::iterator m_iter = index_map.begin();
       m_iter != index_map.end(); ++m_iter)
  {
    vertices[vert_len] = m_iter->first;
    m_iter->second = vert_len++;
  }
  for (size_t* iter = indices; iter != end; ++iter)
  {
    VertexHandle vertex = *(VertexHandle*)iter;
    msq_std::map<VertexHandle,size_t>::iterator pos = index_map.find( vertex );
    *iter = pos->second;
  }      
  
  
}

void MeshTSTTImpl::get_all_sizes( size_t& vertex_count,
                                  size_t& element_count,
                                  size_t& vertex_use_count,
                                  MsqError& err )
{
  try {
    
      // Get number of vertices
    vertex_count   = meshIFace.getNumOfType( nodeSet,    TSTTM::EntityType_VERTEX );
    
      // Get number of elements
    element_count  = meshIFace.getNumOfType( elementSet, TSTTM::EntityType_FACE );
    element_count += meshIFace.getNumOfType( elementSet, TSTTM::EntityType_REGION );
    
      // Get number of vertex uses
    sidl::array<void*> handles1, handles2;
    sidl::array<int> offsets1, offsets2;
    sidl::array<int> flags1, flags2;
    int num_handles, num_offsets, num_flags;
      // Face elements
    meshIFace.getAdjEntities( elementSet, 
                              TSTTM::EntityType_FACE,
                              TSTTM::EntityTopology_ALL_TOPOLOGIES,
                              TSTTM::EntityType_VERTEX,
                              handles1, num_handles,
                              offsets1, num_offsets,
                              flags1  , num_flags );
    vertex_use_count = num_handles;
      // Region elements
    meshIFace.getAdjEntities( elementSet, 
                              TSTTM::EntityType_REGION,
                              TSTTM::EntityTopology_ALL_TOPOLOGIES,
                              TSTTM::EntityType_VERTEX,
                              handles2, num_handles,
                              offsets2, num_offsets,
                              flags2  , num_flags );
    vertex_use_count += num_handles;
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

void MeshTSTTImpl::get_all_mesh( VertexHandle* vert_array,  size_t vert_len,
                                 ElementHandle* elem_array, size_t elem_len,
                                 size_t* offset_array,      size_t offset_len,
                                 size_t* conn_array,        size_t conn_len ,
                                 MsqError& err )
{
  int num_elem;
  try {
    sidl::array<void*> elements( convert_to_sidl_vector( elem_array, elem_len ) );
    meshIFace.getEntities( elementSet, 
                           TSTTM::EntityType_ALL_TYPES,
                           TSTTM::EntityTopology_ALL_TOPOLOGIES,
                           elements,
                           num_elem );
    if ((unsigned)num_elem > elem_len)
    {
      MSQ_SETERR(err)("Insufficient space in array", MsqError::OUT_OF_MEMORY);
      return;
    }
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
  
  if (offset_len <= (unsigned)num_elem)
  {
    MSQ_SETERR(err)("Insufficient space in array", MsqError::OUT_OF_MEMORY);
    return;
  }
  
  elements_get_attached_vertices( elem_array, num_elem,
                                  vert_array, vert_len,
                                  conn_array, conn_len,
                                  offset_array, err );
}
      

// Returns the topology of the given entity.
EntityTopology MeshTSTTImpl::element_get_topology(
  ElementHandle element, MsqError &err)
{
  try {
    TSTTM::EntityTopology topo = entIFace.getEntTopo( element );
    return topologyMap[topo];
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return MIXED;
  }
}

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
void MeshTSTTImpl::elements_get_topologies(
  ElementHandle *element_handle_array,
  EntityTopology *element_topologies,
  size_t num_elements, MsqError &err)
{
  try {
    
    sidl::array<TSTTM::EntityTopology> topologies( alloc_sidl_vector<TSTTM::EntityTopology>( num_elements ) );
    sidl::array<void*> handles( convert_to_sidl_vector( element_handle_array, num_elements ) );
    int topo_len_out;
    arrIFace.getEntArrTopo( handles, num_elements, topologies, topo_len_out );
    if ((size_t)topo_len_out != num_elements) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return;
    }
    
    for (unsigned i = 0; i < num_elements; ++i)
      element_topologies[i] = topologyMap[topologies.get(i)];
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void MeshTSTTImpl::release_entity_handles(
  Mesquite::Mesh::EntityHandle */*handle_array*/,
  size_t /*num_handles*/, MsqError &/*err*/)
{
    // Do nothing
}

// Instead of deleting a Mesh when you think you are done,
// call release().  In simple cases, the implementation could
// just call the destructor.  More sophisticated implementations
// may want to keep the Mesh object to live longer than Mesquite
// is using it.
void MeshTSTTImpl::release()
{
}

//**************** Tags ****************
TagHandle MeshTSTTImpl::tag_create( const msq_std::string& name, 
                                    TagType type, unsigned length,
                                    const void* default_val,
                                    MsqError& err )
{
  TSTTB::TagValueType tstt_type;
  size_t size = 0;
  switch (type) {
    case Mesquite::Mesh::BYTE:   size = sizeof(char  ); tstt_type = TSTTB::TagValueType_OPAQUE;        break;
    case Mesquite::Mesh::INT:    size = sizeof(int   ); tstt_type = TSTTB::TagValueType_INTEGER;       break;
    case Mesquite::Mesh::DOUBLE: size = sizeof(double); tstt_type = TSTTB::TagValueType_DOUBLE;        break;
    case Mesquite::Mesh::HANDLE: size = sizeof(void* ); tstt_type = TSTTB::TagValueType_ENTITY_HANDLE; break;
    default:
      MSQ_SETERR(err)("Invalid tag type", MsqError::INVALID_ARG );
      return 0;
  }
  size *= length;
  
  try {
    void* handle = 0;
    tagIFace.createTag( name, size, tstt_type, handle );
    return handle;
  } 
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return 0;
  }
}

void MeshTSTTImpl::tag_delete( TagHandle handle, MsqError& err )
{
  try {
    tagIFace.destroyTag( handle, true );
  } 
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

TagHandle MeshTSTTImpl::tag_get( const msq_std::string& name, MsqError& err )
{
  try {
    return tagIFace.getTagHandle( name );
  } 
  catch(::TSTTB::Error &tstt_err) {
    if (tstt_err.getErrorType() != TSTTB::ErrorType_TAG_NOT_FOUND) {
      MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    }
    return 0;
  }
}

void MeshTSTTImpl::tag_properties( TagHandle handle,
                                   msq_std::string& name,
                                   TagType& type_out,
                                   unsigned& length_out,
                                   MsqError& err )
{
  size_t size;
  TSTTB::TagValueType type;
  try {
    name = tagIFace.getTagName( handle );
    size = tagIFace.getTagSize( handle );
    type = tagIFace.getTagType( handle );
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
    return;
  }
  
  size_t tsize;
  switch (type) {
    case TSTTB::TagValueType_OPAQUE       : tsize = sizeof(char  ); type_out = Mesquite::Mesh::BYTE  ; break;
    case TSTTB::TagValueType_INTEGER      : tsize = sizeof(int   ); type_out = Mesquite::Mesh::INT   ; break;
    case TSTTB::TagValueType_DOUBLE       : tsize = sizeof(double); type_out = Mesquite::Mesh::DOUBLE; break;
    case TSTTB::TagValueType_ENTITY_HANDLE: tsize = sizeof(void* ); type_out = Mesquite::Mesh::HANDLE; break;
    default:
      MSQ_SETERR(err)("Unsupported TSTT tag type", MsqError::NOT_IMPLEMENTED );
      return;
  }
  
  if (size % tsize != 0)
  {
    MSQ_SETERR(err)("Invalid TSTT tag size", MsqError::INTERNAL_ERROR );
    return;
  }
  
  length_out = size / tsize;
}

void MeshTSTTImpl::tag_set_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         const void* data,
                                         MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}

void MeshTSTTImpl::tag_set_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        const void* data,
                                        MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}
    
void MeshTSTTImpl::tag_set_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 const void* data,
                                 MsqError& err )
{
  try {
    size_t len, size = tagIFace.getTagSize( tag );
    int count;
    sidl::array<void*> handles( convert_to_sidl_vector( (void**)array, num_elems ) );
    switch (tagIFace.getTagType( tag ))
    {
      case TSTTB::TagValueType_ENTITY_HANDLE:
      {
        len = size / sizeof(void*);
        sidl::array<void*> sdata( convert_to_sidl_vector( (void**)data, len*num_elems ));
        arrTagIFace.setEHArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_DOUBLE:
      {
        len = size / sizeof(double);
        sidl::array<double> sdata( convert_to_sidl_vector( (double*)data, len*num_elems ));
        arrTagIFace.setDblArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_INTEGER:
      {
        len = size / sizeof(int);
        sidl::array<int> sdata( convert_to_sidl_vector( (int*)data, len*num_elems ));
        arrTagIFace.setIntArrData( handles, num_elems, tag, sdata, count );
      }
      break;

      default:
      {
#if TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_OPAQUE_PACKED
        len = num_elems / sizeof(void*);
        if (num_elems % sizeof(void*))
          ++len;
        sidl::array<void*> sdata( convert_to_sidl_vector( (void**)data, len ) );
#elif TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_OPAQUE_PADDED
        assert( size <= sizeof(void*) );
        sidl::array<void*> sdata( alloc_sidl_vector<void*>(num_elems) );
        const char* ptr = reinterpret_cast<const char*>(data);
        for (size_t i = 0; i < num_elems; ++i)
          sdata.set( i, reinterpret_cast<void*>(*(ptr + i*size)) );
#elif TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_CHAR
        len = size * num_elems;
        sidl::array<char> sdata( convert_to_sidl_vector( (char*)data, len ) );
#elif TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_UCHAR \
   || TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_BYTE
        len = size * num_elems;
        sidl::array<unsigned char> sdata( alloc_sidl_vector( (unsigned char*)data, len ) );
#else
#error
#endif
        arrTagIFace.setArrData( handles, num_elems, tag, sdata, count, size );
      }
    }
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}


void MeshTSTTImpl::tag_get_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         void* data,
                                         MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}

void MeshTSTTImpl::tag_get_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        void* data,
                                        MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}
    
void MeshTSTTImpl::tag_get_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 void* data,
                                 MsqError& err )
{
  try {
    size_t len, size = tagIFace.getTagSize( tag );
    int32_t lower = 0, upper = num_elems - 1, stride = 1, count = num_elems;
    sidl::array<void*> handles;
    handles.borrow( const_cast<void**>(array), 1, &lower, &upper, &stride );
    switch (tagIFace.getTagType( tag ))
    {
      case TSTTB::TagValueType_ENTITY_HANDLE:
      {
        len = size / sizeof(void*);
        sidl::array<void*> sdata( convert_to_sidl_vector( (void**)data, len*num_elems ));
        arrTagIFace.getEHArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_DOUBLE:
      {
        len = size / sizeof(double);
        sidl::array<double> sdata( convert_to_sidl_vector( (double*)data, len*num_elems ));
        arrTagIFace.getDblArrData( handles, num_elems, tag, sdata, count );
      }
      break;
      
      case TSTTB::TagValueType_INTEGER:
      {
        len = size / sizeof(int);
        sidl::array<int> sdata( convert_to_sidl_vector( (int*)data, len*num_elems ));
        arrTagIFace.getIntArrData( handles, num_elems, tag, sdata, count );
      }
      break;

      default:
      {
#if TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_OPAQUE_PACKED
        len = num_elems / sizeof(void*);
        if (num_elems % sizeof(void*))
          ++len;
        sidl::array<void*> sdata( convert_to_sidl_vector( (void**)data, len ) );
#elif TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_OPAQUE_PADDED
        assert( size <= sizeof(void*) );
        sidl::array<void*> sdata( alloc_sidl_vector<void*>(num_elems) );
#elif TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_CHAR
        len = size * num_elems;
        sidl::array<char> sdata( convert_to_sidl_vector( (char*)data, len ) );
#elif TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_UCHAR \
   || TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_BYTE
        len = size * num_elems;
        sidl::array<unsigned char> sdata( alloc_sidl_vector( (unsigned char*)data, len ) );
#else
#error
#endif
        int32_t junk;
        arrTagIFace.getArrData( handles, num_elems, tag, sdata, count, junk );
        
#if TSTT_OPAQUE_TAG_TYPE == OPAQUE_TYPE_OPAQUE_PADDED
        const char* ptr = reinterpret_cast<const char*>(data);
        for (size_t i = 0; i < num_elems; ++i)
          sdata.set( i, reinterpret_cast<void*>(*(ptr + i*size)) );
#endif
      }
    }
  }
  catch(::TSTTB::Error &tstt_err) {
    MSQ_SETERR(err)( process_tstt_error(tstt_err), MsqError::INTERNAL_ERROR );
  }
}

} // namespace Mesquite
