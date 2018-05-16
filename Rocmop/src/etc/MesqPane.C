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
/* ****
   MesqPane.C

*** */

#ifndef MESQPANE_C
#define MESQPANE_C

#include "MesqPane.h"
#include "Element_accessors.h"

using namespace MOP;
using namespace std;
namespace Mesquite{

  MesqPane::~MesqPane(){
     _pane = NULL;
     if(_dc)
       delete _dc;
     _dc = NULL;
  }

  void MesqPane::invert(){
    std::vector<COM::Connectivity*> conn;
    _pane->connectivities(conn);

    int temp = 0;
    int *cptr = NULL;
    uint offset = 0;
    for(uint i=0; i<conn.size(); ++i){
      if(_verb>1){
	cout << "Connectivity " << i <<" has name " << conn[i]->name() << endl;
	cout << "offset = " << offset << endl;
      }
      cptr = conn[i]->pointer();
      uint nj = _with_ghost ? (uint)conn[i]->size_of_elements() : 
	(uint)conn[i]->size_of_real_elements();
      for(uint j = 0; j< nj; ++j){
	temp = cptr[offset+4*j];
	cptr[offset+4*j] = cptr[offset+4*j+2];
	cptr[offset+4*j+2] = temp;
      }
    }
  }

  void MesqPane::init(){
    if(_verb)
      cout << "MesqPane::init" << endl;
    _dc = new MAP::Pane_dual_connectivity(_pane, _with_ghost);
    if(_verb>1)
      cout << "Building pb" << endl;
    MAP::Pane_boundary pb(_pane);
    if(_verb>1)
      cout << "Finished building pb" << endl;
    std::vector<bool> _is_isolated;
    _is_isolated.clear();
    if(_verb>1)
      cout << "Determining border nodes" << endl;
    if(_with_ghost)
      pb.determine_border_nodes(_is_border, _is_isolated, NULL, 1);
    else
      pb.determine_border_nodes(_is_border, _is_isolated, NULL, 0);
    if(_verb>1)
      cout << "Finished determining border nodes" << endl;
    int siz = _is_border.size();
    for(int i =0; i < siz; ++i){
      if(_is_border[i] && _verb >1)
	cout << " node " << i+1 << " is on the border" << endl;
      if(_is_isolated[i] && _verb >1)
	cout << " node " << i+1 << " is isolated" << endl;
      if(!(_is_border[i] || _is_isolated[i]) && _verb >1)
	cout << " node " << i+1 << " is free" << endl;
      _is_border[i] = (_is_border[i] || _is_isolated[i]);
      
    }
    _vertexBytes.resize(_with_ghost ? _pane->size_of_nodes() :
			_pane->size_of_real_nodes());
    if(_verb>1){
      cout << "  size_of_elements       = " << _pane->size_of_elements() << endl;
      cout << "  size_of_ghost_elements = " << _pane->size_of_ghost_elements() << endl;
      cout << "  size_of_real_elements  = " << _pane->size_of_real_elements() << endl;
      cout << "  size_of_nodes          = " << _pane->size_of_nodes() << endl;
      cout << "  size_of_ghost_nodes    = " << _pane->size_of_ghost_nodes() << endl;
      cout << "  size_of_real_nodes     = " << _pane->size_of_real_nodes() << endl;
    }
    if(_verb>0)
      cout << "MesqPane::init done" << endl;
  }
    
  void MesqPane::get_all_vertices(VertexHandle *vert_array,
				   size_t array_size, MsqError &err){
    if(_verb)
      cout << "MesqPane::get_all_vertices" << endl;
    int pane_size = _with_ghost ? _pane->size_of_nodes() : 
      _pane->size_of_real_nodes();
    COM_assertion_msg( (array_size >= (size_t)pane_size),
		       "Array_Size Must Be At Least Number of Nodes");
    for(int i = 1; i <= pane_size; ++i){
      vert_array[i-1] = ((char*)NULL)+i;
    }
    if(_verb)
      cout << "MesqPane::get_all_vertices done" << endl;
  }

  void MesqPane::get_all_elements(ElementHandle *elem_array,
				   size_t array_size, MsqError &err){
    if(_verb)
      cout << "MesqPane::get_all_elements" << endl;
    int pane_size = _with_ghost ? _pane->size_of_elements() :
      _pane->size_of_real_elements();
    COM_assertion_msg( (array_size >= (size_t)pane_size),
		       "Array_Size Must Be At Least Number of Nodes");
    for(int i = 1; i <= pane_size; ++i){
      elem_array[i-1] = ((char*)NULL)+i;
    }
    if(_verb)
      cout << "MesqPane::get_all_elements done" << endl;
  }

  VertexIterator* MesqPane::vertex_iterator(MsqError &err){
    if(_verb>1)
      cout << "MesqPane::vertex_iterator" << endl;
    return new MyEntityIterator(_with_ghost ? _pane->size_of_nodes() :
				_pane->size_of_real_nodes());
  }

  ElementIterator* MesqPane::element_iterator(MsqError &err){
    if(_verb>1)
      cout << "MesqPane::element_iterator" << endl;
    return new MyEntityIterator(_with_ghost ? _pane->size_of_elements() :
				_pane->size_of_real_elements());
  }

  bool MesqPane::vertex_is_fixed(VertexHandle vertex, MsqError & err){
    if(_verb>1)
      cout << "MesqPane::vertex_is_fixed" << endl;
    uint node_id = (char*)vertex - (char*)NULL;
    COM_assertion_msg(node_id <= _is_border.size(),
		      "Invalid vertex for MesqPane::vertex_is_fixed()");
    return _is_border[node_id-1];
    if(_verb>1)
      cout << "MesqPane::vertex_is_fixed finished" << endl;
  }

  void MesqPane::vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
 					  size_t num_vtx, MsqError &err){
    if(_verb)
      cout << "MesqPane::vertices_are_on_boundary" << endl;
    for (size_t i = 0; i < num_vtx; ++i){
      uint node_id = (((char*)vert_array[i]) - (char*)NULL);
      if (_is_border[node_id-1]){
	on_bnd[i] = true;
	if(_verb>1)
	  std::cout << "Node " << node_id << " is on the boundary." << endl;
      }
      else{
	on_bnd[i] = false;
      }
    }
    if(_verb)
      cout << "MesqPane::vertices_are_on_boundary finished" << endl;
  }

  void MesqPane::vertices_get_coordinates(VertexHandle vert_array[],
					  Mesquite::MsqVertex* const &coordinates,
					  const size_t &num_vtx,
					  MsqError &err){
    if(_verb)
      cout << "MesqPane::vertices_get_coordinates" << endl;
    const double *xptr = _pane->coordinates(), *yptr=NULL, *zptr=NULL;
    const COM::Attribute *nc = _pane->attribute(COM::COM_NC);
    int stride = nc->stride();
    if(!stride){
      const COM::Attribute *xc = _pane->attribute(COM::COM_NC1);
      const COM::Attribute *yc = _pane->attribute(COM::COM_NC2);
      const COM::Attribute *zc = _pane->attribute(COM::COM_NC3);
      int stride1 = xc->stride();
      int stride2 = yc->stride(); 
      int stride3 = zc->stride();
      xptr = _pane->x_coordinates();
      yptr = _pane->y_coordinates();
      zptr = _pane->z_coordinates();
      for (size_t i = 0; i < num_vtx; ++i){
	int offset = (char*)vert_array[i]-(char*)NULL -1;
	coordinates[i].set(*(xptr+offset*stride1),
			   *(yptr+offset*stride2),
			   *(zptr+offset*stride3));
      }
    }
    else if (stride==1){
      int nn = _with_ghost ? _pane->size_of_nodes() : 
	_pane->size_of_real_nodes();
      xptr = _pane->x_coordinates();
      for (size_t i = 0; i < num_vtx; ++i){
	int offset = (char*)vert_array[i]-(char*)NULL -1;
	coordinates[i].set(*(xptr+offset),
			   *(xptr+offset+nn),
			   *(xptr+offset+2*nn));
      }      
    }
    else if (stride>=3){
      xptr = _pane->x_coordinates();
      for (size_t i = 0; i < num_vtx; ++i){
	int offset = (char*)vert_array[i]-(char*)NULL -1;
	coordinates[i].set(xptr + offset * stride);
	if(_verb>1){
	  cout << "Coords of node " << offset+1 << " = ["
	       << coordinates[i][0] << ","
	       << coordinates[i][1] << ","
	       << coordinates[i][2] << "]" << endl;
	}
      }      
    }
    else{
      COM_assertion_msg( 0,
			 "Invalid stride length for nodal coordinates.");
    }
    if(_verb)
      cout << "MesqPane::vertices_get_coordinates finished" << endl;
  }
  
  void MesqPane::vertex_set_coordinates(VertexHandle vertex,
					const Vector3D &coordinates,
					MsqError &err){
    int node_id = ((char*)vertex-(char*)NULL);
    if(_verb > 1)
      cout << "MesqPane::vertex_set_coordinates, node = " << node_id << endl;
    int offset = node_id - 1;
    double * xptr = _pane->coordinates();
    if(xptr){
      xptr += offset *3;
      if(_verb>1){
	cout << "Coordinate: [" << *xptr << " , "
	     << *(xptr+1) << " , " << *(xptr+2) << "], ";      
      }
      coordinates.get_coordinates(xptr);
      if(_verb>1){
	cout << "-> [" << *xptr << " , "
	     << *(xptr+1) << " , " << *(xptr+2) << "]" << endl;
      }
    }
    else{
      coordinates.get_coordinates(*(_pane->x_coordinates()+offset),
				  *(_pane->y_coordinates()+offset),
				  *(_pane->z_coordinates()+offset));
    }
    if(_verb)
      cout << "MesqPane::vertex_set_coordinates finished" << endl;
  }

  void MesqPane::vertex_set_byte (VertexHandle vertex,
				unsigned char byte, MsqError &err){
    if(_verb)
      cout << "MesqPane::vertex_set_byte" << endl;
    int offset = ((char*)vertex - (char*)NULL) - 1;
    _vertexBytes[offset] = byte;    
    if(_verb)
      cout << "MesqPane::vertex_set_byte finished" << endl;  
}

  void MesqPane::vertices_set_byte (VertexHandle *vert_array,
				  unsigned char *byte_array,
				  size_t array_size, MsqError &err){
    if(_verb)
      cout << "MesqPane::vertices_set_byte" << endl;
    int offset =0;
    for(size_t i = 0; i < array_size; ++i){ 
      offset = ((char*)vert_array[i] - (char*)NULL) - 1;
      _vertexBytes[offset] = byte_array[i];
    }
    if(_verb)
      cout << "MesqPane::vertices_set_byte" << endl;
  }

  void MesqPane::vertex_get_byte(VertexHandle vertex,
			       unsigned char *byte, MsqError &err){
    if(_verb)
      cout << "MesqPane::vertex_get_byte" << endl;
    int offset = ((char*)vertex - (char*)NULL) - 1;
    *byte = _vertexBytes[offset];
  }

  void MesqPane::vertices_get_byte(VertexHandle *vertex,
				 unsigned char *byte_array,
				 size_t array_size, MsqError &err){
    if(_verb)
      cout << "MesqPane::vertices_get_byte" << endl;
    int offset =0;
    for(size_t i = 0; i < array_size; ++i){ 
      offset = ((char*)vertex[i] - (char*)NULL) - 1;
      byte_array[i] = _vertexBytes[offset];
    }
    if(_verb)
      cout << "MesqPane::vertices_get_byte finished" << endl;
  }

  size_t MesqPane::vertex_get_attached_element_count(VertexHandle vertex,
						     MsqError &err) const{
    if(_verb)
      cout << "MesqPane::vertex_get_attached_element_count" << endl;
    int node_id = ((char*)vertex - (char*)NULL);		    
    std::vector<int> elist;
    _dc->incident_elements(node_id,elist);
    return elist.size();
    if(_verb)
      cout << "MesqPane::vertex_get_attached_element_count finished" << endl;
  }

  void MesqPane::vertex_get_attached_elements(VertexHandle vertex,
					      ElementHandle* elem_array,
					      size_t sizeof_elem_array,
					      MsqError &err){
    if(_verb)
      cout << "MesqPane::vertex_get_attached_elements" << endl;
    int node_id = ((char*)vertex - (char*)NULL);		    
    std::vector<int> elist;
    _dc->incident_elements(node_id,elist);    
    if(sizeof_elem_array > elist.size())
      sizeof_elem_array = elist.size();
    for(uint i = 0; i < sizeof_elem_array; i++){
      elem_array[i] = ((char*)NULL+elist[i]);
    }
    if(_verb)
      cout << "MesqPane::vertex_get_attached_elements finished" << endl;
  }

  size_t MesqPane::element_get_attached_vertex_count(ElementHandle elem,
						   MsqError &err) const{
    if(_verb)
      cout << "MesqPane::element_get_attached_vertex_count" << endl;
    int element_id = ((char*)elem - (char*)NULL);		        
    const COM::Connectivity *con = _pane->connectivity(element_id);
    return con->size_of_nodes_pe();
  }

  void MesqPane::elements_get_attached_vertices(ElementHandle *elem_handles,
					      size_t num_elems,
					      VertexHandle *vert_handles,
					      size_t &sizeof_vert_handles,
					      size_t *csr_data,
					      size_t &sizeof_csr_data,
					      size_t *csr_offsets,
					      MsqError &err){
    if(_verb)
      cout << "MesqPane::elements_get_attached_vertices" << endl;
    int vert_count = 0;
    // maps vertex id's to position in vert_handles
    std::map<int,int> itop_map;
    std::map<int,int>::iterator pos;
    if (num_elems == 0)
      return;
    csr_offsets[0] = 0;
    // loop through elements
    for (size_t i = 0; i < num_elems; ++i){
      int elem_id = ((char*)elem_handles[i] - (char*)NULL);
      std::cout << "\n  element " << elem_id << " has nodes : ";
      std::vector<int> elist;
      _dc->incident_elements(elem_id,elist);
      size_t nodes_in_elem = elist.size();
      csr_offsets[i+1] = csr_offsets[i] + nodes_in_elem;
      // Check for space in csr_data
      COM_assertion_msg( (sizeof_csr_data >= csr_offsets[i+1]),
			 "Not enough space in arg. csr_data");	

      // Loop through vertices of current element
      for( uint j = 0 ; j < nodes_in_elem; ++j){
	int node_id = ((char*)elist[j] - (char*)NULL);		        	
	std::cout << node_id << " ";
	pos = itop_map.find(node_id);
	// current vertex isn't in vert_handles, add it
	if(pos == itop_map.end()){
	  itop_map.insert(std::map<int,int>::value_type(node_id,num_elems));
	  vert_handles[vert_count] = ((char*)(NULL)+node_id);
	  csr_data[csr_offsets[i]+j] = vert_count;
	  ++vert_count;  
	}
	// add current vertex to csr_data
	csr_data[csr_offsets[i] + j] = pos->second;
      }
      elist.clear();
    }
    sizeof_csr_data = csr_offsets[num_elems];
    sizeof_vert_handles = vert_count;
    if(_verb)
      cout << "MesqPane::elements_get_attached_vertices" << endl;
  }

  void MesqPane::elements_get_attached_vertex_indices(ElementHandle element[],
						      size_t num_elems,
						      size_t index_array[],
						      size_t array_size,
						      size_t* offsets,
						      MsqError &err){
    if(_verb)
      cout << "MesqPane::elements_get_attached_vertex_indices" << endl;

    int elem_id;
    offsets[0] = 0;

    Element_node_enumerator *ene;
    for(uint i = 0; i < num_elems; ++i){

      elem_id = ((char*)element[i] - (char*)NULL);
      if(_verb>1){
	cout << "    element " << elem_id << " contains nodes ";
      }
	
      ene = new Element_node_enumerator(_pane,elem_id);
      int nodes_in_elem = ene->size_of_nodes();
      offsets[i+1] = offsets[i] + nodes_in_elem;
      // loop through vertices in element
      std::vector<int> nodes;
      ene->get_nodes(nodes);
      
      for(int j = 0; j < nodes_in_elem; ++j){
	//index_array[offsets[i] + j] = con_array[first_offset+nodes_in_elem*(elem_id-1)+j] -1;
	index_array[offsets[i]+j] = nodes[j]-1;
	if(_verb>1){
	  cout << nodes[j]
	     << " ";
	}
	delete ene;
	ene = NULL;
      }
      if(_verb>1)
	cout << endl;
    }
    if(_verb>1)
      cout << "  MesqPane::e_g_a_v_i loop finished  " << endl;
  }

  EntityTopology MesqPane::element_get_topology(ElementHandle entity_handle,
						MsqError &err) const{
    if(_verb)
      cout << "MesqPane::element_get_topology" << endl;
    int elem_id = ((char*)entity_handle-(char*)NULL);

    const COM::Connectivity *con = _pane->connectivity(elem_id);
    COM_assertion_msg(con!=NULL, "con is NULL");
    switch(con->element_type()){
    case COM::Connectivity::TRI3 :
      return TRIANGLE;
      break;
    case COM::Connectivity::QUAD4 :
      return QUADRILATERAL;
      break;
    case COM::Connectivity::TET4 :
      return TETRAHEDRON;
      break;
    case COM::Connectivity::HEX8 :
      return HEXAHEDRON;
      break;
    default:
      COM_assertion_msg( 0,
			 "Element type not supported by Mesquite.");
      break;
    }
    // Control should never reach here...
    // this 'return' gets rid of compiler warnings
    return TRIANGLE;
  }

  void MesqPane::elements_get_topologies(ElementHandle *element_handle_array,
					 EntityTopology *element_topologies,
					 size_t num_elements, MsqError &err){
    if(_verb >0)
      cout << "MesqPane::element_get_topologies" << endl;
    for (size_t i = 0; i < num_elements; i++){
      element_topologies[i] = element_get_topology(element_handle_array[i],err);
    }
    if(_verb >0)
      cout << "MesqPane::element_get_topologies done" << endl;
 }

  void MesqPane::element_tag_create(const string tag_name,
				    int tag_size,
				    TagHandle& tag_handle,
				    MsqError &err)
  {
    if(_verb)
      cout << "MesqPane::element_tag_create" << endl;
    MesqPane::tag new_tag;
    new_tag.pt = 0;
    new_tag.size = tag_size;
    denseTags[tag_name] = new_tag;
    tag_handle = tag_get_handle(tag_name, err);
    MSQ_CHKERR(err);
    if(_verb)
      cout << "MesqPane::element_tag_create finished" << endl;  }
    
  //! Returns the tag handle corresponding to the string.
  //! Sets an error if the string is not recognised
  void* MesqPane::tag_get_handle(const string tag_name, MsqError &err){
    if(_verb)
      cout << "MesqPane::tag_get_handle" << endl;
    return (void*) &(denseTags[tag_name]);
    if(_verb)
      cout << "MesqPane::tag_get_handle finished" << endl;
  }
  
  
  //! Sets the Tag data pointers for an array of entities.
  //! The memory pointed to is the reponsibility of the function caller, i.e. either the
  //! memory was obtained from entities_get_tag_data or the function caller allocated the memory
  //! and will keep it allocated as long as the tag pointer is set to that memory.
  void MesqPane::elements_set_tag_data(const size_t num_elements, 
				       const TagHandle tag_handle,
				       TagDataPt const tag_data_array,
				       const int& tag_size,
				       MsqError &err)
  {
    if(_verb)
      cout << "MesqPane::element_set_tag_data" << endl;
    if (_with_ghost && (num_elements != _pane->size_of_elements())) {
      err.set_msg("Incorrect num_elements. Must be equal to the total number of elements "
		  "since only dense tags are supported.");
      return;
    }
    else if(num_elements != _pane->size_of_real_elements()){
      err.set_msg("Incorrect num_elements. Must be equal to the number of real elements "
		  "since only dense tags are supported.");
      return;
    }
    else if (((tag*)tag_handle)->size != tag_size) {
      err.set_msg("tag_size does not correspong to actual tag size.");
      return;
    }
    else
      ((tag*)tag_handle)->pt = tag_data_array;
    if(_verb)
      cout << "MesqPane::element_set_tag_data" << endl;
  }
  
  //! Gets the Tag data for an array of entities.
  //! The implementation handles the pointers to the tag data, not the memory pointed to.
  //! This means that the implementation is not responsible for making sure the
  //! memory pointed to is valid
  void MesqPane::elements_get_tag_data(const size_t num_elements,
				       const TagHandle tag_handle,
				       TagDataPt &tag_data_array,
				       int& tag_size,
				       MsqError &err){
    if(_verb)
      cout << "MesqPane::elements_get_tag_data" << endl;
    if (_with_ghost && (num_elements != _pane->size_of_elements())) {
      err.set_msg("Incorrect num_elements. Must be equal to the total number of elements "
		  "since only dense tags are supported.");
      return;    }
    else if (num_elements != _pane->size_of_real_elements()) {
      err.set_msg("Incorrect num_elements. Must be equal to the number of real elements "
		  "since only dense tags are supported.");
      return;
    }
    else if (((tag*)tag_handle)->size != tag_size) {
      err.set_msg("tag_size does not correspong to actual tag size.");
      return;
    }
    else
      tag_data_array = ((tag*)tag_handle)->pt;
    if(_verb)
      cout << "MesqPane::elements_get_tag_data finished" << endl;
  }  

  void MesqPane::release_entity_handles(EntityHandle *handle_array,
					size_t num_handles, MsqError &err){

    if(_verb)
      cout << "MesqPane::release_entity_handles" << endl;
  }
  
  void MesqPane::release()
  {
    if(_verb)
      cout << "MesqPane::release()" << endl;
  }
  
} // end namespace Mesquite
#endif






