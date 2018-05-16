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
// $Id: MesqPane.h,v 1.3 2008/12/06 08:45:24 mtcampbe Exp $
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

#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "MesquiteError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorithms
#include "MeanRatioQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;
#include "MeshInterface.hpp"
#include "MesquiteError.hpp"
#include "MsqVertex.hpp"

#include "Pane.h"
#include "Pane_boundary.h"
#include "Dual_connectivity.h"
#include "mopbasic.h"
#include "roccom_devel.h"

#define MESQUITE_BEGIN_NAMESPACE namespace Mesquite{
#define MESQUITE_END_NAMESPACE }

using namespace std;
USE_MOP_NAMESPACE
MESQUITE_BEGIN_NAMESPACE

using std::string;

class EntityIterator;
typedef EntityIterator VertexIterator;
typedef EntityIterator ElementIterator;

//! A class enabling Mesquite calls on Rocmop panes.
/**
 * Mesquite's Mesh base class defines an interface which allows
 * Mesquite to interact with non-native mesh types. MesqPane
 * is an implementation of this class which enables Mesquite
 * to use a Roccom Pane as an input mesh.
 */
class MesqPane : public Mesquite::Mesh {
  
  inline size_t vertices_in_topology(EntityTopology);
  
public:

  /** \name Constructors and Destructors
   *\{
   */
  
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

  //! Destructor.
  virtual ~MesqPane();

  //\}
  // end of Constructors and Destructors

  /** \name Miscellaneous functions
   *\{
   */

  //! Invert Tetrahedrons
  /**
   *  Assumes that the connectivities of the Pane are all unstructured
   *  tetrahedral, and inverts the tetrahedrons by swapping the 2nd
   *  and 4th nodes in each element.
   */
  void invert();
  
  //! Initialize the MesqPane
  /**
   *  Builds the dual connectivity, determines border nodes,
   *  and resizes data structures.
   */
  void init();

  //! Set the MesqPane verbose level (int, >= 0)
  void set_verb(int verb){
    _verb = verb;
  }

  //\} End of Miscellaneous functions
  
  
  /** \name Inherited Operations on Entire Mesh
   * \{
   */

  //! Returns whether this mesh lies in a 2D or 3D coordinate system.
  virtual int get_geometric_dimension(MsqError &err) const {
    if(_verb >1) std::cout << "MesqPane::get_geometric_dimension" << std::endl;
    return 3;  
  }
  
  //! Returns the number of nodes
  virtual size_t get_total_vertex_count(MsqError &err) const {
    if(_verb >1) std::cout << "MesqPane::get_total_vertex_count" << std::endl;
    return _with_ghost ? _pane->size_of_nodes() : _pane->size_of_real_nodes();  
  }

  //! Returns the number of elements
  virtual size_t get_total_element_count(MsqError &err) const {
    if(_verb >1) std::cout << "MesqPane::get_total_element_count" << std::endl;
    return _with_ghost ? _pane->size_of_elements() : _pane->size_of_real_elements();  }

  
  //! Fills array with handles to all vertices.
  /**
   * \param vert_array a pointer to the vertex array to be filled.
   * \param array_size Must be at least the number of vertices.
   * If less than the mesh number of vertices, causes a run time
   * error.
   */
  virtual void get_all_vertices(VertexHandle *vert_array,
				size_t array_size, MsqError &err);

  //! Fills array with handles to all elements in the mesh.
  /**
   * \param elem_array a pointer to the element handle array.
   * \param array_size Must be at least the number of elements.
   * If less than the mesh number of elements, causes a run time
   * error.
   */
  virtual void get_all_elements(ElementHandle *elem_array,
				size_t array_size, MsqError &err);
  

  //! Returns a pointer to a vertex iterator over this Pane
  /** 
   * The calling code should
   * delete the returned iterator when it is finished with it.
   * If vertices are added or removed from the Mesh after obtaining
   * an iterator, the behavior of that iterator is undefined.
   */
  virtual VertexIterator* vertex_iterator(MsqError &err);
  
  //! Returns a pointer to an element iterator over this Pane
  /** 
   * The calling code should
   * delete the returned iterator when it is finished with it.
   * If elements are added or removed from the Mesh after obtaining
   * an iterator, the behavior of that iterator is undefined.
   */
  virtual ElementIterator* element_iterator(MsqError &err);

  //\}
  // end Inherited Operations on Entire Mesh
  
  /** \name Inherited vertex Properties
   * \{
   */
  
  //! Tells wheter the vertex is allowed to be repositioned.
  /** 
   * True indicates that the vertex is fixed and cannot be moved.  
   * This flag is not modifiable by users of the Mesquite::Mesh interface.
   */
  virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err);
  
  //! Tells whether the vertex  is on the boundary.  
  /** 
   * Boundary nodes may be treated as a special case by some algorithms 
   * or culling methods.
   * This flag is not modifiable by users of the Mesquite::Mesh interface.
   */
  virtual void vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
					size_t num_vtx, MsqError &err);
  
  //! Get the location of a vertex
  /**
   * \param vert_array the array of vertex handles for this Pane.
   * \param coordinates Mesquite's array of vertex coordinates.
   * \param num_vtx number of vertices to process.
   */
  virtual void vertices_get_coordinates(VertexHandle vert_array[],
                                        Mesquite::MsqVertex* const &coordinates,
                                        const size_t &num_vtx,
                                        MsqError &err);

  //! Set the Location of a vertex
  /**
   * \param vertex the vertex in question.
   * \param coordinates Mesquite's coordinates for this vertex.
   */
  virtual void vertex_set_coordinates(VertexHandle vertex,
				      const Vector3D &coordinates,
				      MsqError &err);
  
  //! Set a vertex flag
  /**
   * Each vertex has a byte-sized flag that can be used to store
   * flags.  This byte's value is neither set nor used by the Pane.  
   * It is intended to be used by Mesquite algorithms.
   * Until a vertex's byte has been explicitly set, its value is 0.
   *
   * \param vertex a handle to the vertex in question.
   * \param byte the new flag value.
   */
  virtual void vertex_set_byte (VertexHandle vertex,
				unsigned char byte, MsqError &err);

  //! Set mutliple vertex flags
  /**
   * Each vertex has a byte-sized flag that can be used to store
   * flags.  This byte's value is neither set nor used by the Pane.  
   * It is intended to be used by Mesquite algorithms.
   * Until a vertex's byte has been explicitly set, its value is 0.
   *
   * \param vert_array handles to the vertices in question.
   * \param byte the new flag value.
   */
  virtual void vertices_set_byte (VertexHandle *vert_array,
				  unsigned char *byte_array,
				  size_t array_size, MsqError &err);
  
  //! Retrieve the byte value for the specified vertex.
  /**
   * The byte value is 0 if it has not yet been set via one of the
   * _set_byte() functions.
   */
  virtual void vertex_get_byte(VertexHandle vertex,
			       unsigned char *byte, MsqError &err);

  //! Retrieve the byte value for the specified vertices.
  /**
   * The byte value is 0 if it has not yet been set via one of the
   * _set_byte() functions.
   */
  virtual void vertices_get_byte(VertexHandle *vertex,
				 unsigned char *byte_array,
				 size_t array_size, MsqError &err);

  //\}
  // end Inherited Vertex Properties
  
  /** \name Inherited Vertex Topology
   * \{
   */
  
  //! Gets the number of elements attached to this vertex.
  //! Useful to determine how large the "elem_array" parameter
  //! of the vertex_get_attached_elements() function must be.
  virtual size_t vertex_get_attached_element_count(VertexHandle vertex,
						   MsqError &err) const;
  
  //! Gets the elements attached to this vertex.
  virtual void vertex_get_attached_elements(VertexHandle vertex,
					    ElementHandle* elem_array,
					    size_t sizeof_elem_array,
					    MsqError &err);
  
  //\}
  // end Inherited Vertex Topology

  /** \name Inherited Element Topology
   * \{
   */

  //! Gets the number of vertices in this element.
  virtual size_t element_get_attached_vertex_count(ElementHandle elem,
						   MsqError &err) const;
  
  //! Returns vertices which are part of each element in the "elem_handles" array
  /**
   * When this function is called, the
   * following must be true:
   * -# "elem_handles" points at an array of "num_elems" element handles.
   * -# "vert_handles" points at an array of size "sizeof_vert_handles"
   *  -# "csr_data" points at an array of size "sizeof_csr_data"
   * -# "csr_offsets" points at an array of size "num_elems+1"
   * 
   * When this function returns, adjacency information will be stored
   * in csr format:
   * -# "vert_handles" stores handles to all vertices found in one
   * or more of the elements.  Each vertex appears only
   * once in "vert_handles", even if it is in multiple elements.
   * -# "sizeof_vert_handles" is set to the number of vertex
   * handles placed into "vert_handles".
   * -# "sizeof_csr_data" is set to the total number of vertex uses (for
   * example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
   * the two triangles share some vertices).
   * -# "csr_offsets" is filled such that csr_offset[i] indicates the location
   * of entity i's first adjacency in "csr_data".  The number of vertices
   * in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
   *  reason, csr_offsets[num_elems] is set to the new value of
   * "sizeof_csr_data".
   * -# "csr_data" stores integer offsets which give the location of
   * each adjacency in the "vert_handles" array.
   * 
   * As an example of how to use this data, you can get the handle of the first
   * vertex in element #3 like this:
   * \code VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ] \endcode
   * 
   * and the second vertex of element #3 like this:
   * \code VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ] \endcode
  */
  virtual void elements_get_attached_vertices(ElementHandle *elem_handles,
					      size_t num_elems,
					      VertexHandle *vert_handles,
					      size_t &sizeof_vert_handles,
					      size_t *csr_data,
					      size_t &sizeof_csr_data,
					      size_t *csr_offsets,
					      MsqError &err);
  
  //! Return each vertex's global index
  /**
  * Identifies the vertices attached to the elements by returning
  * each vertex's global index.  The vertex's global index indicates
  * where that vertex can be found in the array returned by get_all_vertices()
  * \param elems Array of element handles.
  * \param num_elems number of elements in the array elems
  * \param index_array Array containing the indexes of the elements vertices.
  *                    Indexes can be repeated.
  * \param index_array_size indicates the size of index_array
  * \param offsets Has length num_elems+1.
  *                An array of offsets into the index_array that indicates where
  *                the indexes corresponding to an element start. First entry is 0.
  * For example element, to get the vertex indices of elems[3], we look at the entries
  * index_array[offsets[3]] to index_array[offsets[4]] .
  */
  virtual void elements_get_attached_vertex_indices(ElementHandle element[],
						    size_t num_elems,
						    size_t index_array[],
						    size_t array_size,
						    size_t* offsets,
						    MsqError &err);
  
  //! Returns the topology of the given entity.
  virtual EntityTopology element_get_topology(ElementHandle entity_handle,
                                              MsqError &err) const;

  //! Returns the topologies of the given entities.  
  /**
   *The "entity_topologies" array must be at least "num_elements" in size.
   */
  virtual void elements_get_topologies(ElementHandle *element_handle_array,
				       EntityTopology *element_topologies,
				       size_t num_elements, MsqError &err);
  
  //\}
  // end Inherited Element Topology


  /** \name Inherited Tags
   * \{
   */

  //! Add a new tag with given name and size to the dense tag map.
  virtual void element_tag_create(const string tag_name,
				  int tag_size,
				  TagHandle& tag_handle,
				  MsqError &err);
    
  //! Returns the tag handle corresponding to the string.
   virtual void* tag_get_handle(const string tag_name, MsqError &err);
  
  
  //! Sets the Tag data pointers for an array of entities.
  /**
  * The memory pointed to is the reponsibility of the function caller, i.e. either the
  * memory was obtained from entities_get_tag_data or the function caller allocated the memory
  * and will keep it allocated as long as the tag pointer is set to that memory.
  */
  virtual void elements_set_tag_data(const size_t num_elements, 
				     const TagHandle tag_handle,
				     TagDataPt const tag_data_array,
				     const int& tag_size,
				     MsqError &err);
  
  //! Gets the Tag data for an array of entities.
  /**
   * The implementation handles the pointers to the tag data, not the memory pointed to.
   * This means that the implementation is not responsible for making sure the
   * memory pointed to is valid.
   */
  virtual void elements_get_tag_data(const size_t num_elements,
				     const TagHandle tag_handle,
				     TagDataPt &tag_data_array,
				     int& tag_size,
				     MsqError &err);
  //\}
  // end Inherited Tags
  
  /** \name Inherited Memory Management
   * \{
   */

  //! Tells the mesh that the client is finished with a given entity handle.  
  virtual void release_entity_handles(EntityHandle *handle_array,
				      size_t num_handles, MsqError &err);
  
  //! Release a MesqPane from use.
  /**
   * Instead of deleting a Mesh when you think you are done,
   * call release().  In simple cases, the implementation could
   * just call the destructor.  More sophisticated implementations
   * may want to keep the Mesh object to live longer than Mesquite
   * is using it.
   */
  virtual void release();
  
  //\}
  // end of Inherited Memory Management.


protected:
  COM::Pane* _pane;
  MAP::Pane_dual_connectivity* _dc;
  std::vector<bool> _is_border;
  std::vector<unsigned char> _vertexBytes;

  int  _verb; // 1 for function start/end
               // 2 for all information
               // 0 = nothing

  bool _with_ghost; // If true(default), then use ghost information.


  //! A tag structure required by Mesquite
  struct tag {
    void* pt; // points to the beginning of the dense tag array
    int size; // size is the increment to use to go from one tag to the next
  };
  std::map<std::string, MesqPane::tag> denseTags; 
};

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
// end namespace MESQUITE
#endif






