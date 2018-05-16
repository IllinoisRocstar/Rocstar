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

#include "MsqVertex.hpp"
#include "MesqPane_1_1.h"
#include "Element_accessors.h"

using namespace MOP;
using namespace std;
namespace Mesquite{

  MesqPane::~MesqPane(){
    _pane = NULL;
    if(_dc)
      delete _dc;
    _dc = NULL;
    std::map<string,tagStruct>::iterator pos;
    for(pos = s_to_t.begin(); pos!= s_to_t.end(); ++pos){
      switch ((pos->second).type){
      case BYTE : {
	if(pos->second.edata)
	  delete ((std::vector<char>*)pos->second.edata);
	if(pos->second.ndata)
	  delete ((std::vector<char>*)pos->second.ndata);
	break;
      }
      case BOOL :{
	if(pos->second.edata)
	  delete ((std::vector<char>*)pos->second.edata);
	if(pos->second.ndata)
	  delete ((std::vector<char>*)pos->second.ndata);
	break;
      }
      case INT :{
	if(pos->second.edata)
	  delete ((std::vector<char>*)pos->second.edata);
	if(pos->second.ndata)
	  delete ((std::vector<char>*)pos->second.ndata);
	break;
      }
      case DOUBLE :{
	if(pos->second.edata)
	  delete ((std::vector<char>*)pos->second.edata);
	if(pos->second.ndata)
	  delete ((std::vector<char>*)pos->second.ndata);
	break;
      }
      case HANDLE : {
	if(pos->second.edata)
	  delete ((std::vector<char>*)pos->second.edata);
	if(pos->second.ndata)
	  delete ((std::vector<char>*)pos->second.ndata);
	break;
      }
      }
    }
  }
  
  void MesqPane::invert(int conn_type){
    std::vector<COM::Connectivity*> conn;
    _pane->connectivities(conn);

    int temp = 0;
    int *cptr = NULL;
    uint offset = 0;
    for(uint i=0; i<conn.size(); ++i){
      if(conn[i]->element_type() == conn_type &&
	 conn_type == COM::Connectivity::TET4){
	cptr = conn[i]->pointer();
	uint nj = _with_ghost ? (uint)conn[i]->size_of_elements() : 
	  (uint)conn[i]->size_of_real_elements();
	for(uint j = 0; j< nj; ++j){
	  temp = cptr[offset+4*j];
	  cptr[offset+4*j] = cptr[offset+4*j+2];
	  cptr[offset+4*j+2] = temp;
	}
      }
      else if(conn[i]->element_type() == conn_type &&
	      conn_type == COM::Connectivity::HEX8){
	cptr = conn[i]->pointer();
	uint nj = _with_ghost ? (uint)conn[i]->size_of_elements() : 
	  (uint)conn[i]->size_of_real_elements();
	for(uint j = 0; j< nj; ++j){
	  temp = cptr[offset+8*j+5];
	  cptr[offset+8*j+5] = cptr[offset+8*j+7];
	  cptr[offset+8*j+7] = temp;
	  temp = cptr[offset+8*j+1];
	  cptr[offset+8*j+1] = cptr[offset+8*j+3];
	  cptr[offset+8*j+3] = temp;
	}	
      }
    }
  }

  void MesqPane::init(){
    if(_verb)
      cout << "MOP> MesqPane::init" << endl;
    _dc = new MAP::Pane_dual_connectivity(_pane, _with_ghost);
    if(_verb>1)
      cout << "MOP> Building pb" << endl;
    MAP::Pane_boundary pb(_pane);
    if(_verb>1)
      cout << "MOP> Finished building pb" << endl;
    std::vector<bool> _is_isolated;
    _is_isolated.clear();
    if(_verb>1)
      cout << "MOP> Determining border nodes" << endl;
    if(_with_ghost)
      pb.determine_border_nodes(_is_border, _is_isolated, NULL, 1);
    else
      pb.determine_border_nodes(_is_border, _is_isolated, NULL, 0);
    if(_verb>1)
      cout << "MOP> Finished determining border nodes" << endl;
    int siz = _is_border.size();
    for(int i =0; i < siz; ++i){
      if(_is_border[i] && _verb >1)
	cout << "MOP> node " << i+1 << " is on the border" << endl;
      if(_is_isolated[i] && _verb >1)
	cout << "MOP> node " << i+1 << " is isolated" << endl;
      if(!(_is_border[i] || _is_isolated[i]) && _verb >1)
	cout << "MOP> node " << i+1 << " is free" << endl;
      _is_border[i] = (_is_border[i] || _is_isolated[i]);
      
    }
    _vertexBytes.resize(_with_ghost ? _pane->size_of_nodes() :
			_pane->size_of_real_nodes());
    if(_verb>1){
      cout << "MOP>  size_of_elements       = " << _pane->size_of_elements() << endl;
      cout << "MOP>  size_of_ghost_elements = " << _pane->size_of_ghost_elements() << endl;
      cout << "MOP>  size_of_real_elements  = " << _pane->size_of_real_elements() << endl;
      cout << "MOP>  size_of_nodes          = " << _pane->size_of_nodes() << endl;
      cout << "MOP>  size_of_ghost_nodes    = " << _pane->size_of_ghost_nodes() << endl;
      cout << "MOP>  size_of_real_nodes     = " << _pane->size_of_real_nodes() << endl;
    }
    if(_verb>0)
      cout << "MOP> MesqPane::init done" << endl;
  }

  // NEW
  void MesqPane::get_all_elements( msq_std::vector<ElementHandle>& elements,
				   MsqError& err ){
    elements.clear();
    int e_size = _with_ghost ? _pane->size_of_elements() :
      _pane->size_of_real_elements();
    elements.resize(e_size);
    for(int i=1; i<= e_size; ++i){
      elements[i-1] = ((char*)NULL)+i;
    }
  }
  
  // NEW
  void MesqPane::get_all_vertices( msq_std::vector<VertexHandle>& vertices,
                                   MsqError& err ){
    vertices.clear();
    int v_size = _with_ghost ? _pane->size_of_nodes() :
      _pane->size_of_real_nodes();
    vertices.resize(v_size);
    for(int i=1; i<= v_size; ++i){
      vertices[i-1] = ((char*)NULL)+i;
    }
  }


  void MesqPane::get_all_mesh( VertexHandle*  vert_array, size_t vert_len,
			       ElementHandle* elem_array, size_t elem_len,
			       size_t* elem_conn_offsets, size_t offset_len,
			       size_t* elem_conn_indices, size_t index_len,
			       MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::get_all_mesh" << endl;
    int pane_size = _with_ghost ? _pane->size_of_nodes() : 
      _pane->size_of_real_nodes();

    COM_assertion_msg( (vert_len >= (size_t)pane_size),
		       "Vert_Len Must Be At Least Number of Nodes");

    for(int i = 1; i <= pane_size; ++i)
      vert_array[i-1] = ((char*)NULL)+i;
    
    pane_size = _with_ghost ? _pane->size_of_elements() :
      _pane->size_of_real_elements();

    COM_assertion_msg( (elem_len >= (size_t)pane_size),
		       "Array_Size Must Be At Least Number of Nodes");

    elem_conn_offsets[0] = 0;
    Element_node_enumerator ene(_pane,1);
    std::vector<int> nodes;

    for(int i = 1; i <= pane_size; ++i, ene.next()){
      
      elem_array[i-1] = ((char*)NULL)+i;

      elem_conn_offsets[i] = elem_conn_offsets[i-1]+ ene.size_of_nodes();

      ene.get_nodes(nodes);
      for(int j=0, nj=ene.size_of_nodes(); j<nj; ++j)
	elem_conn_indices[elem_conn_offsets[i-1]+j] = nodes[j]-1;
    }
  }
    
  VertexIterator* MesqPane::vertex_iterator(MsqError &err){
    if(_verb>1)
      cout << "MOP> MesqPane::vertex_iterator" << endl;
    return new MyEntityIterator(_with_ghost ? _pane->size_of_nodes() :
				_pane->size_of_real_nodes());
  }

  ElementIterator* MesqPane::element_iterator(MsqError &err){
    if(_verb>1)
      cout << "MOP> MesqPane::element_iterator" << endl;
    return new MyEntityIterator(_with_ghost ? _pane->size_of_elements() :
				_pane->size_of_real_elements());
  }

  void MesqPane::vertices_get_fixed_flag( const VertexHandle vert_array[], 
					  bool fixed_flag_array[],
					  size_t num_vtx, 
					  MsqError &err ){
    for(size_t i =0; i< num_vtx; ++i){
      uint node_id = (char*)vert_array[i] - (char*)NULL;
      fixed_flag_array[i] = _is_border[node_id-1];
    }
  }

  bool MesqPane::vertex_is_fixed(VertexHandle vertex, MsqError & err){
    if(_verb>1)
      cout << "MOP> MesqPane::vertex_is_fixed" << endl;
    uint node_id = (char*)vertex - (char*)NULL;
    COM_assertion_msg(node_id <= _is_border.size(),
		      "Invalid vertex for MesqPane::vertex_is_fixed()");
    return _is_border[node_id-1];
    if(_verb>1)
      cout << "MOP> MesqPane::vertex_is_fixed finished" << endl;
  }

  void MesqPane::vertices_are_on_boundary(VertexHandle vert_array[], bool on_bnd[],
 					  size_t num_vtx, MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertices_are_on_boundary" << endl;
    for (size_t i = 0; i < num_vtx; ++i){
      uint node_id = (((char*)vert_array[i]) - (char*)NULL);
      if (_is_border[node_id-1]){
	on_bnd[i] = true;
	if(_verb>1)
	  std::cout << "MOP> Node " << node_id << " is on the boundary." << endl;
      }
      else{
	on_bnd[i] = false;
      }
    }
    if(_verb)
      cout << "MOP> MesqPane::vertices_are_on_boundary finished" << endl;
  }

  void MesqPane::vertices_get_coordinates(const VertexHandle vert_array[],
					  MsqVertex* coordinates,
					  size_t num_vtx,
					  MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertices_get_coordinates" << endl;
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
	  cout << "MOP> Coords of node " << offset+1 << " = ["
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
      cout << "MOP> MesqPane::vertices_get_coordinates finished" << endl;
  }
  
  void MesqPane::vertex_set_coordinates(VertexHandle vertex,
					const Vector3D &coordinates,
					MsqError &err){
    int node_id = ((char*)vertex-(char*)NULL);
    if(_verb)
      cout << "MOP> MesqPane::vertex_set_coordinates, node = " << node_id << endl;
    int offset = node_id - 1;
    double * xptr = _pane->coordinates();
    if(xptr){
      xptr += offset *3;
      if(_verb>2){
	cout << "MOP> Coordinate: [" << *xptr << " , "
	     << *(xptr+1) << " , " << *(xptr+2) << "], ";      
      }
      coordinates.get_coordinates(xptr);
      if(_verb>2){
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
      cout << "MOP> MesqPane::vertex_set_coordinates finished" << endl;
  }

  void MesqPane::vertex_set_byte (VertexHandle vertex,
				  unsigned char byte, MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertex_set_byte" << endl;
    int offset = ((char*)vertex - (char*)NULL) - 1;
    _vertexBytes[offset] = byte;    
    if(_verb)
      cout << "MOP> MesqPane::vertex_set_byte finished" << endl;  
  }

  void MesqPane::vertices_set_byte (VertexHandle *vert_array,
				    unsigned char *byte_array,
				    size_t array_size, MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertices_set_byte" << endl;
    int offset =0;
    for(size_t i = 0; i < array_size; ++i){ 
      offset = ((char*)vert_array[i] - (char*)NULL) - 1;
      _vertexBytes[offset] = byte_array[i];
    }
    if(_verb)
      cout << "MOP> MesqPane::vertices_set_byte" << endl;
  }

  // NEW
  void MesqPane::vertices_set_byte( const VertexHandle *vert_array,
				    const unsigned char *byte_array,
				    size_t array_size, 
				    MsqError &err ){
    int offset =0;
    for(size_t i = 0; i < array_size; ++i){ 
      offset = ((char*)vert_array[i] - (char*)NULL) - 1;
      _vertexBytes[offset] = byte_array[i];
    }
  }

  void MesqPane::vertex_get_byte(VertexHandle vertex,
				 unsigned char *byte, MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertex_get_byte" << endl;
    int offset = ((char*)vertex - (char*)NULL) - 1;
    *byte = _vertexBytes[offset];
  }

  void MesqPane::vertices_get_byte(VertexHandle *vertex,
				   unsigned char *byte_array,
				   size_t array_size, MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertices_get_byte" << endl;
    int offset =0;
    for(size_t i = 0; i < array_size; ++i){ 
      offset = ((char*)vertex[i] - (char*)NULL) - 1;
      byte_array[i] = _vertexBytes[offset];
    }
    if(_verb)
      cout << "MOP> MesqPane::vertices_get_byte finished" << endl;
  }

  size_t MesqPane::vertex_get_attached_element_count(VertexHandle vertex,
						     MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertex_get_attached_element_count" << endl;
    int node_id = ((char*)vertex - (char*)NULL);		    
    std::vector<int> elist;
    _dc->incident_elements(node_id,elist);
    return elist.size();
    if(_verb)
      cout << "MOP> MesqPane::vertex_get_attached_element_count finished" << endl;
  }

  void MesqPane::vertex_get_attached_elements(VertexHandle vertex,
					      ElementHandle* elem_array,
					      size_t sizeof_elem_array,
					      MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::vertex_get_attached_elements" << endl;
    int node_id = ((char*)vertex - (char*)NULL);		    
    std::vector<int> elist;
    _dc->incident_elements(node_id,elist);    
    if(sizeof_elem_array > elist.size())
      elist.resize(sizeof_elem_array);
    for(uint i = 0; i < sizeof_elem_array; i++){
      elem_array[i] = ((char*)NULL+elist[i]);
    }
    if(_verb)
      cout << "MOP> MesqPane::vertex_get_attached_elements finished" << endl;
  }

  void MesqPane::
  vertices_get_attached_elements(const VertexHandle* vertex_array,
				 size_t num_vertex,
				 msq_std::vector<ElementHandle>& elements,
				 msq_std::vector<size_t>& offsets,
				 MsqError& err ){

    // Guess 6 elements to each node for initial memory allocation
    // to avoid frequent memory allocation and swapping.
    elements.reserve(6*num_vertex);
    offsets.resize(num_vertex+1);
    if(num_vertex > 0)
      offsets[0] = 0;
    std::vector<int> elist;
    elist.clear();
    
    for(size_t i=0; i<num_vertex; ++i){
      offsets[i] += elist.size();
      offsets[i+1] = offsets[i];
      int node_id = ((char*)vertex_array[i] - (char*)NULL);      
      _dc->incident_elements(node_id,elist);
      for(uint i=0, ni=elist.size(); i<ni; ++i)
	elements.push_back((char*)NULL + elist[i]);
    }
    offsets[num_vertex] = offsets[num_vertex-1] + elist.size();

    // trim excess capacity
    vector<ElementHandle>(elements).swap(elements);
  }
  
  size_t MesqPane::element_get_attached_vertex_count(ElementHandle elem,
						     MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::element_get_attached_vertex_count" << endl;
    int element_id = ((char*)elem - (char*)NULL);		        
    const COM::Connectivity *con = _pane->connectivity(element_id);
    return con->size_of_nodes_pe();
  }


  size_t MesqPane::get_vertex_use_count( ElementHandle* handle_array,
					 size_t num_handles,
					 MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::get_vertex_use_count" << endl;
    size_t count = 0;
    for(uint i=0; i< num_handles; ++i){
      int e_id = ((char*)handle_array[i] - (char*)NULL);		        
      Element_node_enumerator ene(_pane,e_id);
      count += (size_t)ene.size_of_nodes();
    }
    return count;
  }

  void MesqPane::elements_get_attached_vertices(const ElementHandle *elem_handles,
						size_t num_elems,
						msq_std::vector<VertexHandle>& vert_handles,
						msq_std::vector<size_t>& offsets, 
						MsqError &err){

    //    std::cout << "elements_get_attached_vertices\n";

    offsets.clear();
    vert_handles.clear();
    offsets.resize(num_elems+1);
    offsets[0] = 0;
    std::vector<int> nodes;
    nodes.clear();

    // guess 4 vertices per element
    vert_handles.reserve(4*num_elems);

    //    std::cout << "vert_handles = \n";

    for( size_t i=0; i< num_elems; ++i){
      offsets[i] += nodes.size();
      offsets[i+1] = offsets[i];
      nodes.clear();

      int elem_id = ((char*)elem_handles[i] - (char*)NULL);
      Element_node_enumerator ene(_pane,elem_id);
      ene.get_nodes(nodes);

      for(uint j=0, nj= nodes.size(); j<nj; ++j){
	vert_handles.push_back((char*)NULL + nodes[j]);
	//	std::cout << nodes[j] << " ";
      }
      //      std::cout << "\n";
    }
    offsets[num_elems] += nodes.size();
    msq_std::vector<VertexHandle>(vert_handles).swap(vert_handles);
    //    std::cout << "offsets = \n";
    for(uint i=0; i < offsets.size(); ++i){
      //      std::cout << offsets[i] << " ";
    }
    //    std::cout << "\n";
    
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
      cout << "MOP> MesqPane::elements_get_attached_vertices" << endl;
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
      cout << "MOP> MesqPane::elements_get_attached_vertices" << endl;
  }

  EntityTopology MesqPane::element_get_topology(ElementHandle entity_handle,
						MsqError &err){
    if(_verb)
      cout << "MOP> MesqPane::element_get_topology" << endl;
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
    case COM::Connectivity::PYRIMID5 :
      return PYRAMID;
      break;
    case COM::Connectivity::PRISM6 :
      return PRISM;
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
      cout << "MOP> MesqPane::element_get_topologies" << endl;
    for (size_t i = 0; i < num_elements; i++){
      element_topologies[i] = element_get_topology(element_handle_array[i],err);
    }
    if(_verb >0)
      cout << "MOP> MesqPane::element_get_topologies done" << endl;
  }

   void MesqPane::elements_get_topologies(const ElementHandle *element_handle_array,
					 EntityTopology *element_topologies,
					 size_t num_elements, MsqError &err){
    if(_verb >0)
      cout << "MOP> MesqPane::element_get_topologies" << endl;
    for (size_t i = 0; i < num_elements; i++){
      element_topologies[i] = element_get_topology(element_handle_array[i],err);
    }
    if(_verb >0)
      cout << "MOP> MesqPane::element_get_topologies done" << endl;
  }
 
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
  void* MesqPane::tag_create( const msq_std::string& tag_n,
			      TagType type, unsigned length,
			      const void* default_value,
			      MsqError &err){

    std::string tag_name(*(const_cast<std::string *>(&tag_n)));
    if(_verb)
      cout << "MOP> MesqPane::tag_create" << endl;

    tagStruct new_tag = {length, type, NULL, 
			     NULL, tag_name};
    int esize = 1+(_with_ghost ? _pane->size_of_elements() : 
		 _pane->size_of_real_elements());
    int nsize = 1+(_with_ghost ? _pane->size_of_nodes() :
		 _pane->size_of_real_nodes());
    switch (type){
    case BYTE :{

      std::vector<char>* temp = new std::vector<char>;
      temp->reserve(length*esize);
      temp->resize(length*esize);
      new_tag.edata = (void *)temp;

      temp = new std::vector<char>;
      temp->reserve(nsize*length);
      temp->resize(nsize*length);
      new_tag.ndata = (void *)temp;

      for(uint i=0; i<length*esize; ++i)
	(*static_cast<std::vector<char>*>(new_tag.edata))[i] 
	  = ((char *)default_value)[i];
      for(uint i=0; i<length*nsize; ++i)      
	(*((std::vector<char>*)new_tag.ndata))[i] 
	  = ((char *)default_value)[i];
      break;
    }
    case BOOL :{

      std::vector<bool>* temp = new std::vector<bool>;
      temp->reserve(esize*length);
      temp->resize(esize*length);
      new_tag.edata = (void *) temp;

      temp = new std::vector<bool>;
      temp->reserve(nsize*length);
      temp->resize(nsize*length);
      new_tag.ndata = (void *) temp;

      for(uint i=0; i<esize*length; ++i)
	(*((std::vector<bool>*)new_tag.edata))[i] 
	  = ((bool *)default_value)[i];
      for(uint i=0; i<nsize*length; ++i)
	(*((std::vector<bool>*)new_tag.ndata))[i] 
	  = ((bool *)default_value)[i];
      break;
    }
    case INT :{

      std::vector<int>* temp = new std::vector<int>;
      temp->reserve(esize*length);
      temp->resize(esize*length);
      new_tag.edata = (void *) temp;

      temp = new std::vector<int>;
      temp->reserve(nsize*length);
      temp->resize(nsize*length);     
      new_tag.ndata = (void *) temp;

      for(uint i=0; i<esize*length; ++i)
	(*((std::vector<int>*)new_tag.edata))[i] 
	  = ((int *)default_value)[i];
      for(uint i=0; i<nsize*length; ++i)
	(*((std::vector<int>*)new_tag.ndata))[i] 
	  = ((int *)default_value)[i];
      break;
    }
    case DOUBLE :{

      std::vector<double>* temp = new std::vector<double>;
      temp->reserve(esize*length);
      temp->resize(esize*length);
      new_tag.edata = (void *) temp;

      temp = new std::vector<double>;
      temp->reserve(nsize*length);
      temp->resize(nsize*length);
      new_tag.ndata = (void *) temp;

      for(uint i=0; i<esize*length; ++i)
	(*((std::vector<double>*)new_tag.edata))[i] 
	  = ((double *)default_value)[i];
      for(uint i=0; i<nsize*length; ++i)
	(*((std::vector<double>*)new_tag.ndata))[i] 
	  = ((double *)default_value)[i];
      break;
    }
    case HANDLE : {

      std::vector<void*>* temp = new std::vector<void*>;
      temp->reserve(esize*length);
      temp->resize(esize*length);
      new_tag.edata = (void *) temp;

      temp = new std::vector<void*>;
      temp->reserve(nsize*length);
      temp->resize(nsize*length);      
      new_tag.ndata = (void *) temp;

      for(uint i=0; i<esize*length; ++i)
	(*((std::vector<void*>*)new_tag.edata))[i] 
	  = ((void* *)default_value)[i];
      for(uint i=0; i<nsize*length; ++i)
	(*((std::vector<void*>*)new_tag.ndata))[i] 
	  = ((void* *)default_value)[i];
      break;
    }
    default :
      COM_assertion_msg (0,
			 "Unrecognized TagType.");
      break;
    }
    s_to_t.insert(std::make_pair<msq_std::string,tagStruct>(tag_name,new_tag));
    std::map<msq_std::string,tagStruct>::iterator pos;
    pos = s_to_t.find(tag_name);
    COM_assertion_msg(pos != s_to_t.end(),
		      "Unable to create tag for Mesquite.");
    return &(pos->second);
  }
      
  /** \brief Remove a tag and all corresponding data
   *
   * Delete a tag.
   */
  void MesqPane::tag_delete( TagHandle handle, MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::tag_delete" << endl;
    tagStruct* tag = (tagStruct*)handle;
    switch ( (*tag).type){
    case BYTE : {
      std::vector<char>* v = (std::vector<char>*)(*tag).edata;
      if (v!= NULL)
	delete v;
      v = (std::vector<char>*)(*tag).ndata;
      if (v!= NULL)
	delete v;
      break;
    }
    case BOOL :{
      std::vector<bool>* v = (std::vector<bool>*) (*tag).edata;
      if (v!= NULL)
	delete v;
      v = (std::vector<bool>*) (*tag).ndata;
      if (v!= NULL)
	delete v;
      break;
    }
    case INT :{
      std::vector<int>* v = (std::vector<int>*) (*tag).edata;
      if (v!= NULL)
	delete v;
      v = (std::vector<int>*) (*tag).ndata;
      if (v!= NULL)
	delete v;
      break;
    }
    case DOUBLE :{
      std::vector<double>* v = (std::vector<double>*) (*tag).edata;
      if (v!= NULL)
	delete v;
      v = (std::vector<double>*) (*tag).ndata;
      if (v!= NULL)
	delete v;
      break;
    }
    case HANDLE :{
      std::vector<void*>* v = (std::vector<void*>*) (*tag).edata;
      if (v!= NULL)
	delete v;
      v = (std::vector<void*>*) (*tag).ndata;
      if (v!= NULL)
	delete v;
      break;
    }
    default:
      COM_assertion_msg (0,
			 "Unrecognized TagType.");
      break;
    }
    const msq_std::string n= (*tag).name;
    std::map<msq_std::string,tagStruct>::iterator pos;
    pos = s_to_t.find(n);
    assert(pos != s_to_t.end());
    s_to_t.erase(n);
  }

    
    
  /** \brief Get handle for existing tag, by name. 
   *
   * Check for the existance of a tag given it's name and
   * if it exists return a handle for it.  If the specified
   * tag does not exist, zero should be returned WITHOUT 
   * flagging an error.
   */
  TagHandle MesqPane::tag_get( const msq_std::string& name, 
			       MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::tag_get" << endl;
    std::map<msq_std::string,tagStruct>::iterator pos;
    pos = s_to_t.find(name);
    if(pos == s_to_t.end())
      return NULL;
    else
      return (TagHandle)(&(*pos));
  }
     
  /** \brief Get properites of tag
   *
   * Get data type and number of values per entity for tag.
   * \param handle     Tag to get properties of.
   * \param name_out   Passed back tag name.
   * \param type_out   Passed back tag type.
   * \param length_out Passed back number of values per entity.
   */
  void MesqPane::tag_properties( TagHandle handle,
				 msq_std::string& name_out,
				 TagType& type_out,
				 unsigned& length_out,
				 MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::tag_properties" << endl;
    tagStruct* tag = (tagStruct*)handle;
    name_out = (*tag).name;
    type_out = (*tag).type;
    length_out = (*tag).size;
  }
    
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
  void MesqPane::tag_set_element_data( TagHandle handle,
				       size_t num_elems,
				       const ElementHandle* elem_array,
				       const void* tag_data,
				       MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::tag_set_element_data" << endl;
    tagStruct* tag = (tagStruct*)handle;
    int ne = 1 + (_with_ghost ? _pane->size_of_elements() :
		  _pane->size_of_real_elements());
    switch ((*tag).type){
    case BYTE : {
      std::vector<char>* v = (std::vector<char>*)(*tag).edata;
      uint tsize = ne*(*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(int i=0; i < ne; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i*(*tag).size+j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((char*)elem_array[i]-(char*)NULL); 
	for(uint j=0; j< (*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((char*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case BOOL :{
      std::vector<bool>* v = (std::vector<bool>*)(*tag).edata;
      uint tsize = ne*(*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(int i=0; i < ne; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i*(*tag).size+j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((bool*)elem_array[i]-(bool*)NULL); 
	for(uint j=0; j< (*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((bool*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case INT :{
      std::vector<int>* v = (std::vector<int>*)(*tag).edata;
      uint tsize = ne*(*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(int i=0; i < ne; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i*(*tag).size+j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((int*)elem_array[i]-(int*)NULL); 
	for(uint j=0; j< (*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((int*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case DOUBLE :{
      std::vector<double>* v = (std::vector<double>*)(*tag).edata;
      uint tsize = ne*(*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(int i=0; i < ne; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i*(*tag).size+j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((double*)elem_array[i]-(double*)NULL); 
	for(uint j=0; j< (*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((double*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case HANDLE : {
      std::vector<void*>* v = (std::vector<void*>*)(*tag).edata;
      uint tsize = ne*(*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(int i=0; i < ne; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i*(*tag).size+j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((void**)elem_array[i]-(void**)NULL); 
	for(uint j=0; j< (*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((void**)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    default :
      COM_assertion_msg(0,
		        "TagHandle with unrecgonized TagType received.");
    }
  }


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
  void MesqPane::tag_set_vertex_data ( TagHandle handle,
				       size_t num_elems,
				       const VertexHandle* node_array,
				       const void* tag_data,
				       MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::tag_set_vertex_data" << endl;
    tagStruct* tag = (tagStruct*)handle;
    uint nn = 1 + (_with_ghost ? _pane->size_of_nodes() :
		   _pane->size_of_real_nodes());
    switch ((*tag).type){
    case BYTE : {
      std::vector<char>* v = (std::vector<char>*)(*tag).ndata;
      uint tsize = (nn)* (*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(uint i=0; i < nn; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i * (*tag).size +j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((char*)node_array[i]-(char*)NULL); 
	for(uint j=0; j<(*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((char*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case BOOL :{
      std::vector<bool>* v = (std::vector<bool>*)(*tag).ndata;
      uint tsize = (nn)* (*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(uint i=0; i < nn; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i * (*tag).size +j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((bool*)node_array[i]-(bool*)NULL); 
	for(uint j=0; j<(*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((bool*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case INT :{
      std::vector<int>* v = (std::vector<int>*)(*tag).ndata;
      uint tsize = (nn)* (*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(uint i=0; i < nn; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i * (*tag).size +j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((int*)node_array[i]-(int*)NULL); 
	for(uint j=0; j<(*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((int*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case DOUBLE :{
      std::vector<double>* v = (std::vector<double>*)(*tag).ndata;
      uint tsize = (nn)* (*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(uint i=0; i < nn; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i * (*tag).size +j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((double*)node_array[i]-(double*)NULL); 
	for(uint j=0; j<(*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((double*)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    case HANDLE : {
      std::vector<void*>* v = (std::vector<void*>*)(*tag).ndata;
      uint tsize = (nn)* (*tag).size;
      if (v->size() != tsize){
	// This tag isn't used yet... should only have default value
	COM_assertion_msg(v->size() == (*tag).size,
			  "Invalid TagHandle received from Mesquite");
	v->reserve(tsize);
	v->resize(tsize);
	for(uint i=0; i < nn; ++i){
	  for(uint j=0; j< (*tag).size; ++j)
	    (*v)[i * (*tag).size +j]=(*v)[j];
	}
      }
      for(uint i=0; i<num_elems; ++i){
	int eid = ((void**)node_array[i]-(void**)NULL); 
	for(uint j=0; j<(*tag).size; ++j)
	  (*v)[eid * (*tag).size +j]=((void**)tag_data)[i*(*tag).size+j];
      }
      break;
    }
    }
  }
    
    
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
  void MesqPane::tag_get_element_data( TagHandle handle,
				       size_t num_elems,
				       const ElementHandle* elem_array,
				       void* tag_data,
				       MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::tag_get_element_data" << endl;
    tagStruct* tag = (tagStruct*)handle;
    switch ((*tag).type){
    case BYTE : {
      for(uint i=0; i<num_elems; ++i){
	int eid = ((char*)elem_array[i]-(char*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((char*)tag_data)[i * (*tag).size +j] 
	    = ((char*)((*tag).edata))[eid*(*tag).size+j];
      }
      break;
    }
    case BOOL :{
      for(uint i=0; i<num_elems; ++i){
	int eid = ((bool*)elem_array[i]-(bool*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((bool*)tag_data)[i * (*tag).size +j] 
	    = ((bool*)((*tag).edata))[eid*(*tag).size+j];
      }
      break;
    }
    case INT :{
      for(uint i=0; i<num_elems; ++i){
	int eid = ((int*)elem_array[i]-(int*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((int*)tag_data)[i * (*tag).size +j] 
	    = ((int*)((*tag).edata))[eid*(*tag).size+j];
      }
      break;
    }
    case DOUBLE :{
      for(uint i=0; i<num_elems; ++i){
	int eid = ((double*)elem_array[i]-(double*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((double*)tag_data)[i * (*tag).size +j] 
	    = ((double*)((*tag).edata))[eid*(*tag).size+j];
      }
      break;
    }
    case HANDLE : {
      for(uint i=0; i<num_elems; ++i){
	int eid = ((void**)elem_array[i]-(void**)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((void**)tag_data)[i * (*tag).size +j] 
	    = ((void**)((*tag).edata))[eid*(*tag).size+j];
      }
      break;
    }
    }
  }
    
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
  void MesqPane::tag_get_vertex_data ( TagHandle handle,
				       size_t num_elems,
				       const VertexHandle* node_array,
				       void* tag_data,
				       MsqError& err ){
    if(_verb)
      cout << "MOP> MesqPane::tag_get_vertex_data" << endl;
    tagStruct* tag = (tagStruct*)handle;
    switch ((*tag).type){
    case BYTE : {
      for(uint i=0; i<num_elems; ++i){
	int vid = ((char*)node_array[i]-(char*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((char*)tag_data)[i * (*tag).size +j] 
	    = ((char*)((*tag).ndata))[vid*(*tag).size+j];
      }
      break;
    }
    case BOOL :{
      for(uint i=0; i<num_elems; ++i){
	int vid = ((bool*)node_array[i]-(bool*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((bool*)tag_data)[i * (*tag).size +j] 
	    = ((bool*)((*tag).ndata))[vid*(*tag).size+j];
      }
      break;
    }
    case INT :{
      for(uint i=0; i<num_elems; ++i){
	int vid = ((int*)node_array[i]-(int*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((int*)tag_data)[i * (*tag).size +j] 
	    = ((int*)((*tag).ndata))[vid*(*tag).size+j];
      }
      break;
    }
    case DOUBLE :{
      for(uint i=0; i<num_elems; ++i){
	int vid = ((double*)node_array[i]-(double*)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((double*)tag_data)[i * (*tag).size +j] 
	    = ((double*)((*tag).ndata))[vid*(*tag).size+j];
      }
      break;
    }
    case HANDLE : {
      for(uint i=0; i<num_elems; ++i){
	int vid = ((void**)node_array[i]-(void**)NULL); 	
	for(uint j=0; j<(*tag).size; ++j)
	  ((void**)tag_data)[i * (*tag).size +j] 
	    = ((void**)((*tag).ndata))[vid*(*tag).size+j];
      }
      break;
    }
    }
  }
  
  
  void MesqPane::release_entity_handles(EntityHandle *handle_array,
					size_t num_handles, MsqError &err){

    if(_verb)
      cout << "MOP> MesqPane::release_entity_handles" << endl;
  }

  void MesqPane::release_entity_handles(const EntityHandle *handle_array,
					size_t num_handles, MsqError &err){
    
    if(_verb)
      cout << "MOP> MesqPane::release_entity_handles" << endl;
  }

  void MesqPane::release()
  {
    if(_verb)
      cout << "MOP> MesqPane::release()" << endl;
  }

  
} // end namespace Mesquite
#endif






