/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
// $Id: MesqPane_95.h,v 1.5 2008/12/06 08:45:23 mtcampbe Exp $
/** \file MesqPane.h
 *   MesqPane.h
 *  
 *  An implementation of Mesquite's mesh class as defined in 
 *  MeshInterface.hpp and MeshInterface.cpp.
 *
 *  Many of the comments in this code come from MeshInterface.hpp,
 *  which is Copyright 2004 Sandia Corporation and Argonne National
 *  Laboratory and dsitributed under the GNU Lesser General Public License.
 *
 *  Please see that file for more information.
 *
 */

#ifndef MESQPANE_H
#define MESQPANE_H

#ifndef USE_STD_INCLUDES
#define USE_STD_INCLUDES 1
#endif

#ifndef USE_C_PREFIX_INCLUDES 
#define USE_C_PREFIX_INCLUDES 1
#endif

// MESQUITE INCLUDES
#include "Mesquite.hpp"
#include "TopologyInfo.hpp"

// ROCMOP INCLUDES
#include "mopbasic.h"

// ROCCOM INCLUDES
#include "com_devel.hpp"
#include "Pane.hpp"
#include "Pane_boundary.h"
#include "Dual_connectivity.h"
#include "com_devel.hpp"


#ifdef MSQ_USE_OLD_C_HEADERS
#  include <stddef.h>
#else
#  include <cstddef>
#endif

#include <string>

#define MESQUITE_BEGIN_NAMESPACE namespace Mesquite{
#define MESQUITE_END_NAMESPACE }

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

//using namespace Mesquite;
using namespace std;
using std::string;
USE_MOP_NAMESPACE
MESQUITE_BEGIN_NAMESPACE

class EntityIterator;
class MsqError;
class MsqVertex;
class Vector3D;
typedef EntityIterator VertexIterator;
typedef EntityIterator ElementIterator;


/** Type used to refer to a tag defintion */

class MesqPane : public Mesquite::Mesh {

  inline size_t vertices_in_topology(EntityTopology);

public:
  //************ Type Definitions **************
    //! Opaque EntityHandle type and tag type.
  typedef void* EntityHandle;
  typedef void* TagHandle;        

  // We typedef specific types of EntityHandles just
  // to make it clear what kind of entity is to be
  // returned or used as a function parameter, but the
  // different handle types are not actually distinct.
  typedef EntityHandle VertexHandle;
  typedef EntityHandle ElementHandle;

  //! Invert Tetrahedrons
  /**
   *  Assumes that the connectivities of the Pane are all unstructured
   *  tetrahedral, and inverts the tetrahedrons by swapping the 2nd
   *  and 4th nodes in each element.
   */
  void invert();

  //! Construct from a pane
  /**
   * Build a MesqPane for the given Pane either on all nodes
   * and elements, or just on real nodes and elements.
   *
   * \param p the Pane whose mesh we would like to smooth.
   * \param with_ghost use ghost nodes and cells?
   */
  MesqPane(COM::Pane* p, bool with_ghost = true) 
    : _pane(p), _verb(0), _with_ghost(with_ghost){
    init();
  }
  
  //! Set the MesqPane verbose level (int, >= 0)
  void set_verb(int verb){
    _verb = verb;
  }

  //! Destructor.
  virtual ~MesqPane();

  //! Initialize the MesqPane
  /**
   *  Builds the dual connectivity, determines border nodes,
   *  and resizes data structures.
   */
  void init();

  //************ Operations on entire mesh ****************
    //! Returns whether this mesh lies in a 2D or 3D coordinate system.
  virtual int get_geometric_dimension(MsqError &err) {
    if(_verb >1) 
      std::cout << "MesqPane::get_geometric_dimension\n";
    return 3;  
  }
    
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
			      MsqError& err ){
    vertex_count = _with_ghost ? _pane->size_of_nodes() : 
      _pane->size_of_real_nodes();
    element_count = _with_ghost ? _pane->size_of_elements() : 
      _pane->size_of_real_elements();
    std::vector<const COM::Connectivity*> conns;
    _pane->connectivities(conns);
    vertex_use_count = 0;
    for(int i=0, ni = conns.size(); i<ni; ++i){
      vertex_use_count += conns[i]->size_of_nodes_pe() *
	(_with_ghost ? conns[i]->size_of_elements() :
	 conns[i]->size_of_real_elements());
    }
  }
    
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
    
  //! Returns a pointer to an iterator that iterates over the
  //! set of all vertices in this mesh.  The calling code should
  //! delete the returned iterator when it is finished with it.
  //! If vertices are added or removed from the Mesh after obtaining
  //! an iterator, the behavior of that iterator is undefined.
  virtual VertexIterator* vertex_iterator(MsqError &err);
    
  //! Returns a pointer to an iterator that iterates over the
  //! set of all top-level elements in this mesh.  The calling code should
  //! delete the returned iterator when it is finished with it.
  //! If elements are added or removed from the Mesh after obtaining
  //! an iterator, the behavior of that iterator is undefined.
  virtual ElementIterator* element_iterator(MsqError &err);

  //************ Vertex Properties ********************
    //! Returns true or false, indicating whether the vertex
    //! is allowed to be repositioned.  True indicates that the vertex
    //! is fixed and cannot be moved.  Note that this is a read-only
    //! property; this flag can't be modified by users of the
    //! Mesquite::Mesh interface.
  virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err);

  //! Returns true or false, indicating whether the vertex
  //! is on the boundary.  Boundary nodes may be treated as
  //! a special case by some algorithms or culling methods.
  //! Note that this is a read-only
  //! property; this flag can't be modified by users of the
  //! Mesquite::Mesh interface.
  virtual void vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
					size_t num_vtx, MsqError &err);
    
  //! Get/set location of a vertex
  virtual void vertices_get_coordinates(const VertexHandle vert_array[],
					MsqVertex* coordinates,
					size_t num_vtx,
					MsqError &err);

  virtual void vertex_set_coordinates(VertexHandle vertex,
				      const Vector3D &coordinates,
				      MsqError &err);
    
  //! Each vertex has a byte-sized flag that can be used to store
  //! flags.  This byte's value is neither set nor used by the mesh
  //! implementation.  It is intended to be used by Mesquite algorithms.
  //! Until a vertex's byte has been explicitly set, its value is 0.
  virtual void vertex_set_byte (VertexHandle vertex,
				unsigned char byte, MsqError &err);

  virtual void vertices_set_byte (VertexHandle *vert_array,
				  unsigned char *byte_array,
				  size_t array_size, MsqError &err);
    
  //! Retrieve the byte value for the specified vertex or vertices.
  //! The byte value is 0 if it has not yet been set via one of the
  //! *_set_byte() functions.
  virtual void vertex_get_byte(VertexHandle vertex,
			       unsigned char *byte, MsqError &err);
  virtual void vertices_get_byte(VertexHandle *vertex,
				 unsigned char *byte_array,
				 size_t array_size, MsqError &err);
    
  //**************** Vertex Topology *****************    
    //! Gets the number of elements attached to this vertex.
    //! Useful to determine how large the "elem_array" parameter
    //! of the vertex_get_attached_elements() function must be.
  virtual size_t vertex_get_attached_element_count(VertexHandle vertex,
						   MsqError &err);
    
  //! Gets the elements attached to this vertex.
  virtual void vertex_get_attached_elements(VertexHandle vertex,
					    ElementHandle* elem_array,
					    size_t sizeof_elem_array,
					    MsqError &err);
    
  //*************** Element Topology *************
    
    //! Gets the number of vertices in this element.
    //! This data can also be found by querying the
    //! element's topology and getting the number
    //! of vertices per element for that topology type.

  virtual size_t element_get_attached_vertex_count(ElementHandle elem,
						   MsqError &err);

  virtual size_t get_vertex_use_count( ElementHandle* handle_array,
				       size_t num_handles,
				       MsqError& err );

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
					      MsqError &err);
    
  //! Returns the topology of the given entity.
  virtual EntityTopology element_get_topology(ElementHandle entity_handle,
					      MsqError &err);

  //! Returns the topologies of the given entities.  The "entity_topologies"
  //! array must be at least "num_elements" in size.
  virtual void elements_get_topologies(ElementHandle *element_handle_array,
				       EntityTopology *element_topologies,
				       size_t num_elements, MsqError &err);

    
  //***************  Tags  ***********

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
    
    
  /** \brief Get handle for existing tag, by name. 
   *
   * Check for the existance of a tag given it's name and
   * if it exists return a handle for it.  If the specified
   * tag does not exist, zero should be returned WITHOUT 
   * flagging an error.
   */
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
    //! Tells the mesh that the client is finished with a given
    //! entity handle.  
  virtual void release_entity_handles(EntityHandle *handle_array,
				      size_t num_handles, MsqError &err);
    
  //! Instead of deleting a Mesh when you think you are done,
  //! call release().  In simple cases, the implementation could
  //! just call the destructor.  More sophisticated implementations
  //! may want to keep the Mesh object to live longer than Mesquite
  //! is using it.
  virtual void release();

protected:
  COM::Pane* _pane;
  MAP::Pane_dual_connectivity* _dc;
  std::vector<bool> _is_border;
  std::vector<unsigned char> _vertexBytes;
    
  int  _verb; // 1 for function start/end
  // 2 for all information
  // 0 = nothing
    
  bool _with_ghost; // If true(default), then use ghost information.
  
  struct tagStruct {
    unsigned size; // tag size: # of values per entity
    TagType  type; // tag type
    void *  edata; // elemental data
    void *  ndata; // nodal data
    msq_std::string name;
  };

#if 0
  tagStruct::tagStruct()
    : size(0), type(0), edata(NULL), ndata(NULL), dval(NULL)
  {
    ;
  }
#endif

  std::map<msq_std::string,tagStruct> s_to_t;

};
  

#if 0
inline size_t Mesquite::vertices_in_topology(Mesquite::EntityTopology topo)
{
  return TopologyInfo::corners( topo );
}
#endif

//! A class for iterating through a Panes vertices or elements.
class MyEntityIterator : public Mesquite::EntityIterator{
public:
  
  //! Constructor
  MyEntityIterator(uint size){ _cur = 1; _size = size;}

  //! Destructor
  virtual ~MyEntityIterator()
  {_cur = 0;
    _size = 0;
  }
  
  //! Moves the iterator back to the first entity in the list.
  virtual void restart(){
    _cur = 1;
  };
  
  //! Return the handle currently being pointed to by the iterator.
  virtual Mesh::EntityHandle operator*() const {
    if (_cur > _size)
      return 0;
    else
      return (char*)NULL+_cur;
  }
  
  //! Move to the next entity.
  virtual void operator++(){
    _cur ++;
  }
  
  //! Check if we are at the end of the entites
  /**
   * Returns false until the iterator has been advanced PAST the last 
   * entity once is_at_end() returns true, *iterator returns 0.
   */
  virtual bool is_at_end() const{
    if (_cur <= _size) return 0;
    else return 1;
  }
private:
  uint  _cur; // The current node's position.
  uint  _size;// The size of the element array
};

MESQUITE_END_NAMESPACE

#endif






