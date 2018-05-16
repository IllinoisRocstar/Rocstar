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
#ifdef MESQUITE
#define USE_STD_INCLUDES 1
#define USE_C_PREFIX_INCLUDES 1
#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
//#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "ShapeImprovementWrapper.hpp"

// algorithms
#include "MeanRatioFunctions.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MesqPane_95.h"

using namespace Mesquite;
#endif

#include "Rocmop_2.h"
#include "roccom.h"
#include "Pane.h"
#include "Rocblas.h"
#include "Rocmap.h"
#include "Geometric_Metrics_3.h"
#include "Pane_connectivity.h"
#include "Pane_boundary.h"
#include "Pane_communicator.h"

MOP_BEGIN_NAMESPACE

using MAP::Rocmap;
/*
//FIXME, much of the information gathering can be moved to 
// the smoother specific initialization routine.
void Rocmop::smooth_vol_mesq_ng(double initial_quality){

  if(_verb) 
    std::cout << "      Entering Rocmop::smooth_mesquite_ng" << std::endl;

  const std::string surf_attr("is_surface");
  COM::Attribute* w_is_surface = _buf_window->attribute(surf_attr);
  COM_assertion_msg( w_is_surface, "Unexpected NULL pointer");

  // First Perform Element-based Laplacian smoothing on pane boundary volume nodes
  // implemented as follows:

  // * this part should go into smoother_specific_init
  // 0 Declare and resize variables and buffers
  // 1 Loop through panes
  //   a. Find the set of shared nodes.
  //   b. Remove surface nodes from above set, "leaving submerged boundary nodes"
  //   c. Loop through all elements
  //      if(element contains submersed boundary nodes)
  //          increment submersed_boundary nodes' adj. element count
  //          add element center to each submersed_boundary nodes' new position


  // 2 Sum reduction on nodal positions and adj. element counts
  // 3 Divide new nodal positions by element counts

  // 0 Declare and resize data structures.
  if(_verb > 2) std::cout << "  Declaring Variables\n";

  // Allocate buffer space for new nodal positions
  COM::Attribute *w_new_coords = 
    _buf_window->new_attribute("new_coords", 'n', COM_DOUBLE, 3, "");
  _buf_window->resize_array(w_new_coords, 0);

  // Allocate buffer space for adjacent element counts
  COM::Attribute *w_adj_elem_cnts =
    _buf_window->new_attribute("adj_elem_cnts", 'n', COM_DOUBLE, 1, "");
  _buf_window->resize_array(w_adj_elem_cnts, 0);
  _buf_window->init_done();

  // Allocate buffer space for safe distance
  COM::Attribute *w_safe_dist =
    _buf_window->new_attribute("safe_dist", 'n', COM_DOUBLE, 1, "");
  _buf_window->resize_array(w_safe_dist, 0);
  _buf_window->init_done();

  // Allocate buffer space for backup coords
  COM::Attribute *w_backup =
    _buf_window->new_attribute("backup", 'n', COM_DOUBLE, 3, "");
  _buf_window->resize_array(w_backup, 0);
  _buf_window->init_done();

  // Allocate buffer space for reset flags
  COM::Attribute *w_reset =
    _buf_window->new_attribute("reset", 'n', COM_INT, 1, "");
  _buf_window->resize_array(w_reset, 0);
  _buf_window->init_done();

  // Initialize buffer spaces to zeroes.
  double dzero = 0.;
  double dlarge =99999999;
  int izero = 0;
  Rocblas::copy_scalar( &dzero, w_new_coords);
  Rocblas::copy_scalar( &dzero, w_adj_elem_cnts);
  Rocblas::copy_scalar( &dlarge, w_safe_dist);
  Rocblas::copy_scalar( &izero, w_reset);

  //backup nodal coords.
  Rocblas::copy( _buf_window->attribute(COM::COM_NC), w_backup);

  // Get worst qualities, pre smoothing.

  // Get a list of local panes
  std::vector<COM::Pane*> allpanes;
  _buf_window->panes(allpanes);
  int total_npanes = (int)allpanes.size();

  // Allocate space for the set of submerged_boundary.
  std::vector<std::vector<bool> > is_sub_bnd; // Is submerged boundary node?
  std::vector<std::vector<bool> > is_pan_bnd; // Is pane boundary node?
  is_sub_bnd.resize(total_npanes);
  is_pan_bnd.resize(total_npanes);

  // Get attribute ids
  int new_coords_id = w_new_coords->id();
  int adj_elem_cnts_id = w_adj_elem_cnts->id();
  int is_surface_id = w_is_surface->id();
  int safe_dist_id = w_safe_dist->id();
  int reset_id = w_reset->id();
  int backup_id = w_backup->id();

  // get the pconn offset
  int pconn_offset = MAP::Pane_connectivity::pconn_offset();

  // 1 Loop through panes

  if(_verb > 2) std::cout << "  Step 1\n";
  for(int i=0; i < total_npanes; ++i){
    int size_of_real_nodes = allpanes[i]->size_of_real_nodes();
    int size_of_real_elems = allpanes[i]->size_of_real_elements();
    is_sub_bnd[i].resize(size_of_real_nodes,0);
    is_pan_bnd[i].resize(size_of_real_nodes,0);

    std::vector<bool> is_isolated; // is a node isolated?
    MAP::Pane_boundary pb (allpanes[i]);
    pb.determine_border_nodes(is_pan_bnd[i], is_isolated);

    //   a. Find the set of shared nodes.

    // get pane level pointers
    COM::Attribute *p_is_surface = allpanes[i]->attribute(is_surface_id);
    int *is_surface_ptr = (int*)p_is_surface->pointer();

    double * adj_elem_cnts_ptr = reinterpret_cast<double*>
      (allpanes[i]->attribute(adj_elem_cnts_id)->pointer());

    double * safe_dist_ptr = 
      (double*)(allpanes[i]->attribute(safe_dist_id)->pointer());

    const Vector_3<double> *pnts = reinterpret_cast<Vector_3<double>*>
      (allpanes[i]->attribute(COM_NC)->pointer());

    Vector_3<double> *new_coords_ptr = reinterpret_cast<Vector_3<double>*>
      (allpanes[i]->attribute(new_coords_id)->pointer());
    
    // Obtain the pane connectivity of the local pane
    const COM::Attribute *pconn = allpanes[i]->attribute(COM::COM_PCONN);
    const int *vs = (const int*)pconn->pointer() + pconn_offset;
    int vs_size = pconn->size_of_real_items() - pconn_offset;

    // Determine the number of communicating panes for shared nodes.
    int count=0;
    for (int j=0, nj=vs_size; j<nj; j+=vs[j+1]+2) {
      if (_buf_window->owner_rank( vs[j]) >=0) ++count;
    }

    int index = 0;
    // Loop through communicating panes for shared nodes.
    for ( int j=0; j<count; ++j, index+=vs[index+1]+2) {
      // We skip the panes that are not in the current window 
      while ( _buf_window->owner_rank(vs[index])<0) {
	index+=vs[index+1]+2;
	COM_assertion_msg( index<=vs_size, "Invalid communication map");
      }	
      // Mark node as shared
      for(int k=0; k<vs[index+1]; ++k){
	is_sub_bnd[i][vs[index+2+k]-1]=1;
      }
    }
    
    //   b. Remove surface nodes from above set, leaving submerged boundary nodes
    for (int j =0; j < size_of_real_nodes; ++j){
      if(is_surface_ptr[j])
	is_sub_bnd[i][j] = 0;
    }

    //   c. Loop through all elements
    //      if(element contains submersed boundary nodes)
    //          increment submersed_boundary nodes' adj. element count
    //          add element center to each submersed_boundary nodes' new position
    Element_node_enumerator ene(allpanes[i],1);
    for(int j =0; j < size_of_real_elems; ++j, ene.next()){
      // Collect element indices of submerged boundary nodes.
      int nn = ene.size_of_nodes(), flag = 0;
      for(int k =0; k < nn; ++k){
	if(is_sub_bnd[i][ene[k]-1]){
	  flag = 1;
	  adj_elem_cnts_ptr[ene[k]-1] += 1.0;
	}
      }
      // If the element contains submerged boundary nodes
      if(flag!=0){
	// Calculate the element's center
	Vector_3<double> elem_cntr_pos(0.0,0.0,0.0);
	for(int k =0; k < nn; ++k){
	  elem_cntr_pos += pnts[ene[k]-1];
	}	       
	elem_cntr_pos /= (double)nn;
	// Add the center position to the accumulating new positions
	// of the constituent submerged boundary nodes.
	// And update safe distance if needed	
	for(int k =0; k < nn; ++k){
	  if(is_sub_bnd[i][ene[k]-1]){
	    // initialize safe distance with distance to center
	    double local_safe_dist = (pnts[ene[k]-1] - elem_cntr_pos).norm();
	    int long_src = -1;
	    double long_dist = -1.0;
	    // replace safe distance with .5 shortest adj. edge length
	    // if less than distance to center.
	    int safe_source = -1;
	    for(int ii=1; ii<nn; ++ii){
	      double cur_edge_dist =
		.5*(pnts[ene[k]-1] - pnts[ene[(k+ii)%nn]-1]).norm();
	      if(cur_edge_dist < local_safe_dist){
		safe_source = ii;
		local_safe_dist = cur_edge_dist;
	      }
	      if(cur_edge_dist > long_dist){
		long_dist = cur_edge_dist;
		long_src = ii;
	      }
	    }
	    new_coords_ptr[ene[k]-1] += elem_cntr_pos;
#if 0
	    if(safe_source == -1)
	      new_coords_ptr[ene[k]-1] += elem_cntr_pos;
	    else {
	      Vector_3<double> direction = 
		(pnts[ene[(k+long_src)%nn]-1] - pnts[ene[k]-1]);
	      direction.normalize();
	      new_coords_ptr[ene[k]-1] += pnts[ene[k]-1];
	      new_coords_ptr[ene[k]-1] += direction*local_safe_dist;	      
	      //new_coords_ptr[ene[k]-1] -= pnts[ene[(k+safe_source)%nn]-1];
	    }
#endif
	    if(local_safe_dist < safe_dist_ptr[ene[k]-1])
	      safe_dist_ptr[ene[k]-1] = local_safe_dist;
	  }
	}
      }
    }
  }

  // 2 Sum reduction on nodal positions and adj. element counts
  if(_verb > 1) std::cout << "  Step 2\n";
  reduce_sum_on_shared_nodes(w_new_coords);
  reduce_sum_on_shared_nodes(w_adj_elem_cnts);
  if(COMMPI_Initialized()){
    MAP::Pane_communicator pc(_buf_window, _buf_window->get_communicator());
    pc.init(w_safe_dist);
    pc.begin_update_shared_nodes();
    pc.reduce_on_shared_nodes(MPI_MIN);
    pc.end_update_shared_nodes();
  }

  // 3 Divide new nodal positions by element counts
  if(_verb > 1) std::cout << "  Step 3\n";
  Rocblas::div(w_new_coords, w_adj_elem_cnts, w_new_coords);

  //std::string msg("\n  I submerged quality = ");
  //print_mquality(msg, is_sub_bnd);
  //msg = "  I pane boundary quality = ";
  //print_mquality(msg, is_pan_bnd);

  std::vector<std::vector<bool> > elem_to_check;
  mark_elems_from_nodes(is_sub_bnd, elem_to_check);
  
  // get min,max angles for later use
  double max_angle = 0.0;
  double min_angle = 180.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk = elem_to_check[i].size(); k<nk; ++k){
      if(elem_to_check[i][k]){
	Element_node_enumerator ene(allpanes[i],k+1);
	Angle_Metric_3 am;
	am.initialize(ene);
	am.compute(angles);
	if(angles[1]>max_angle)
	  max_angle = angles[1];
	if(angles[0]<min_angle)
	  min_angle = angles[0];
      }
    }
  }
  int rank =0;
  double temp = min_angle;
  if(COMMPI_Initialized()){
    agree_double(max_angle,MPI_MAX);
    agree_double(min_angle,MPI_MIN);
    int ierr = MPI_Comm_rank( _buf_window->get_communicator()
			      , &rank); 
    assert( ierr == 0);
  }  

  if(_verb > 1) std::cout << "  Step 4\n";
  // 4 Copy new nodal coords into mesh for submerged boundary nodes

  for(int i=0; i < total_npanes; ++i){

    // get pane level attributes
    Vector_3<double>* new_coords_ptr = 
      reinterpret_cast<Vector_3<double>* >(allpanes[i]->
					   attribute(new_coords_id)->
					   pointer());
    Vector_3<double>* coords_ptr = 
      reinterpret_cast<Vector_3<double>* >(allpanes[i]->
					   attribute(COM::COM_NC)
					   ->pointer());
    
    double* safe_dist_ptr = 
      (double*)(allpanes[i]->attribute(safe_dist_id)->pointer());

    for (int j = 0, nj = allpanes[i]->size_of_real_nodes(); j<nj; ++j){
      if(is_sub_bnd[i][j]){
	Vector_3<double> direction = 
	  new_coords_ptr[j] - coords_ptr[j];
	double dist = direction.norm();
	if (dist > safe_dist_ptr[j]){
	  double alpha = .9*safe_dist_ptr[j]/dist;
	  coords_ptr[j] = 
	    alpha*new_coords_ptr[j] +(1.0-alpha)*coords_ptr[j];
	}
	else{
	  coords_ptr[j] = .5*new_coords_ptr[j] +.5*coords_ptr[j];
	}
      }
    }
  }
  //msg = ("\n  L submerged quality = ");
  //print_mquality(msg, is_sub_bnd);
  //msg = "  L pane boundary quality = ";
  //print_mquality(msg, is_pan_bnd);
  //msg = "  L overall quality = ";
  //print_quality(msg);

  // Decide which positions to reset
  for(int i=0; i < total_npanes; ++i){
    
    int* reset_ptr =
      (int*)(allpanes[i]->attribute(reset_id)->pointer());
    
    for(int k =0,nk = elem_to_check[i].size(); k<nk; ++k){
      if(elem_to_check[i][k]){
	double angles[] = {0.0, 0.0};
	
	Element_node_enumerator ene(allpanes[i],k+1);
	Angle_Metric_3 am;
	am.initialize(ene);
	am.compute(angles);
	if(angles[1]>max_angle ||
	   angles[0]<min_angle){
	  std::vector<int> nodes;
	  ene.get_nodes(nodes);
	  for(int j =0, nj=nodes.size(); j<nj; ++j){
	    if (is_sub_bnd[i][nodes[j]-1])
	      reset_ptr[nodes[j]-1] = 1;
	  }
	}
      }
    }
  }

  // All nodes being moved are shared, so this should be valid
  Rocmap::reduce_maxabs_on_shared_nodes(w_reset);

  for(int i=0; i < total_npanes; ++i){
    
    // get pane level attributes
    Vector_3<double>* backup_ptr = 
      reinterpret_cast<Vector_3<double>* >(allpanes[i]->
					   attribute(backup_id)->
					   pointer());
    Vector_3<double>* coords_ptr = 
      reinterpret_cast<Vector_3<double>* >(allpanes[i]->
					   attribute(COM::COM_NC)
					   ->pointer());  
    
    int* reset_ptr =
      (int*)(allpanes[i]->attribute(reset_id)->pointer());
    
    for(int j=0, nj = allpanes[i]->size_of_real_nodes(); j<nj; ++j){
      if(reset_ptr[j])
	coords_ptr[j] = backup_ptr[j];
    }
  }

  //msg = ("\n  R submerged quality = ");
  //print_mquality(msg, is_sub_bnd);
  //msg = "  R pane boundary quality = ";
  //print_mquality(msg, is_pan_bnd);
  //msg = "  R overall quality = ";
  //print_quality(msg);

  // Smooth using mesquite.
  smooth_mesquite(allpanes,0);

  //msg = "\n  M submerged quality = ";
  //print_mquality(msg, is_sub_bnd);
  // msg = "  M pane boundary quality = ";
  //print_mquality(msg, is_pan_bnd);
  //msg = "  M overall quality = ";
  //print_quality(msg);
    
  //Deallocate buffer space
  _buf_window->delete_attribute( w_adj_elem_cnts->name());
  _buf_window->delete_attribute( w_new_coords->name());
  _buf_window->init_done();

  if(_verb > 1) 
    std::cout << "      Exiting Rocmop::smooth_mesquite_ng" << std::endl;
}
*/

//! A structure used to represent element faces.
struct Four_tuple {
  Four_tuple( int a, int b, int c, int d) 
    :  _a(a), _b(b), _c(c), _d(d) {}
  int &operator[](int i) { return (&_a)[i]; }
  const int &operator[](int i) const { return (&_a)[i]; }
  bool operator<( const Four_tuple &x) const {
    return _a<x._a || (_a==x._a && (_b<x._b || _b==x._b && 
				    (_c<x._c || _c==x._c && _d<x._d)));
  }
private:
  int _a, _b, _c, _d;
};

typedef std::map< Four_tuple, MAP::Facet_ID>  Corners2Face_Map;
typedef std::map< int, int> Id_Map;

void Rocmop::determine_physical_border(COM::Attribute* w_is_surface){
  COM_assertion_msg( w_is_surface, "Unexpected NULL pointer");
  COM_assertion_msg( COM_compatible_types( w_is_surface->data_type(), COM_INT),
		     "Surface-list must have integer type");
  COM_assertion_msg( w_is_surface->is_nodal() == 1,
		     "Surface-list must be nodal");
  COM_assertion_msg( w_is_surface->size_of_components() == 1,
		     "Surface-list must have a single component");
  int w_is_surface_id = w_is_surface->id();

  // This function is implemented as follows:
  // 0 Declare and resize data structures.
  // 
  // Loop through panes
  //   1 Get the pane boundary nodes and faces using Rocmap::Pane_boundary
  //   2 Identify adjacent panes, and shared nodes.  Also, create an id mappings
  //     from a node's id to its index in the shared node array and vice versa.
  //   3 Identify "possibly shared faces", those faces whose nodes are all
  //     shared with the same adjacent pane. Make a list of these faces.
  // 4 Communicate the size of the local face list to send to / recieve from
  //   adjacent panes.
  // 5 Create buffers for the incoming face list and communicate faces lists.
  // 6 For every face received, remove any matching faces in the list of 
  //   boundary faces.
  // 7 Mark all nodes still contained in boundary faces as physical boundary
  //   nodes.


  // 0 Declare and resize data structures.

  // vector corresponds to adjacent panes (reused amongst local panes)
  std::vector<bool> is_border; // is a node a border node?
  std::vector<bool> is_shared; // is a node a shared node?
  std::vector<bool> is_isolated; // is a node isolated?
  std::vector<MAP::Facet_ID> border_facet_list;

  // Space in which false pconns are built. (reused amongs local panes)
  // outside vectors correspond to adjacent panes
  // inside vectors correspond to connectivity information
  //std::vector<std::vector<int> > false_pconn;

  // Per-local-pane data structures
  // outside vectors correspond to the local panes
  // inside vectors correspond to adjacent panes
  std::vector<std::set<MAP::Facet_ID> > border_facets; // List of all border faces.
  std::vector<std::vector<Corners2Face_Map> > maybe_shared; // Possibly shared facets
  std::vector<std::vector<Id_Map > > lid2ind; // map from local id to pconn index
  std::vector<std::vector<Id_Map > > ind2lid; // map from pconn index to local id
  std::vector<std::vector<std::set<int> > > adj_pane_set; // nodes shared with adjacent panes
  std::vector<std::vector<int> > adj_pane_id; // ids of adjacent panes
  std::vector<std::vector<int> > adj_pane_recv; // amount of data coming from adjacent panes.

  // The amount of data to be sent by each pane
  std::vector<int> send_sizes;

  std::set<MAP::Facet_ID>::iterator bf_it, bf_it2;
  Id_Map::iterator idm_it,idm_it2; // iterator for lid2ind or ind2lid
  Corners2Face_Map::iterator ms_it,ms_it2;
  std::set<int >::iterator aps_it, aps_it2; //iterator for adj_pane_map

  // get the pconn offset
  int pconn_offset = MAP::Pane_connectivity::pconn_offset();

  // Resize per-local-pane data structures
  std::vector<COM::Pane*> allpanes;
  COM::Window * wrk_window = w_is_surface->window();
  wrk_window->panes(allpanes);
  int total_npanes = (int)allpanes.size(); 
  border_facets.resize(total_npanes);
  maybe_shared.resize(total_npanes);
  lid2ind.resize(total_npanes);
  ind2lid.resize(total_npanes);
  adj_pane_set.resize(total_npanes);
  adj_pane_id.resize(total_npanes);
  adj_pane_recv.resize(total_npanes);
  send_sizes.resize(total_npanes);

  // Register an attribute for our fake pconn.
  COM::Attribute * false_pconn = wrk_window->new_attribute("false_pconn", 'p', COM_INT, 1, "");

  // Register an attribute for sending and receiving sizes
  // also used for sending buffer
  COM::Attribute *com_buff = wrk_window->new_attribute("com_buff", 'p', COM_INT, 1,"");
  int w_com_buff_id = com_buff->id();

  // Loop through panes
  for(int i=0; i<total_npanes; ++i){
    int size_of_real_nodes = allpanes[i]->size_of_real_nodes();
    is_border.clear();
    is_shared.clear();
    is_isolated.clear();
    is_border.resize(size_of_real_nodes,0);
    is_shared.resize(size_of_real_nodes,0);
    is_isolated.resize(size_of_real_nodes,0);
    send_sizes[i] = 0;
    
    //   1 Get the pane boundary nodes and faces using Rocmap::Pane_boundary
    // Determine the border nodes.
    MAP::Pane_boundary pb(allpanes[i]);
    pb.determine_border_nodes(is_border, is_isolated, &border_facet_list);

    // put the list of border facets into a set for faster searching
    for(int j =0, nj=border_facet_list.size(); j<nj; ++j)
      border_facets[i].insert(border_facet_list[j]);
    
    //   2 Identify adjacent panes, and shared nodes.  Also, create an id mappings
    //     from a node's id to its index in the shared node array and vice versa.

    // Obtain the pane connectivity of the local pane
    const COM::Attribute *pconn = allpanes[i]->attribute(COM::COM_PCONN);
    const int *vs = (const int*)pconn->pointer() + pconn_offset;
    int vs_size = pconn->size_of_real_items() - pconn_offset;

    // Determine the number of communicating panes for shared nodes.
    int count=0;
    for (int j=0, nj=vs_size; j<nj; j+=vs[j+1]+2) {
      if (wrk_window->owner_rank( vs[j]) >=0) ++count;
    }

    // Resize per-adjacent-pane data structures
    //false_pconn.resize(count);
    maybe_shared[i].resize(count);
    lid2ind[i].resize(count);
    ind2lid[i].resize(count);
    adj_pane_set[i].resize(count);
    adj_pane_id[i].resize(count);
    adj_pane_recv[i].resize(count);
    
    int index = 0;
    // Loop through communicating panes for shared nodes.
    for ( int j=0; j<count; ++j, index+=vs[index+1]+2) {
      // We skip the panes that are not in the current window 
      while ( wrk_window->owner_rank(vs[index])<0) {
	index+=vs[index+1]+2;
	COM_assertion_msg( index<=vs_size, "Invalid communication map");
      }	
      adj_pane_id[i][j] = vs[index];
      // Add this pane to current shared node's list, mark node as shared and
      // keep track of the mapping between array index and local node id
      for(int k=0; k<vs[index+1]; ++k){
	is_shared[vs[index+2+k]-1]=1;
	adj_pane_set[i][j].insert(vs[index+2+k]);
	lid2ind[i][j].insert(std::make_pair(vs[index+2+k],k));
	ind2lid[i][j].insert(std::make_pair(k,vs[index+2+k]));
      }
    } 
#if 0

    std::cout << "Pane " << allpanes[i]->id() << " lid2ind = \n";
    for(int j=0; j < count; ++j){
      std::cout << " to Pane " << adj_pane_id[i][j] << " =\n";
      for(int k =0; k < ind2lid[i][j].size();++k){
	std::cout << "(" << k << "," << ind2lid[i][j][k] << ") ";
      }
      std::cout << "\n";
    }
    std::cout << "Pane " << i << "\n"
	      << " shared nodes = \n";

    for(int j=0; j<is_shared.size(); ++j){
      if(is_shared[j])
	std::cout << j+1 << " ";
    }
    std::cout << "\n\n nodes shared with othe panes";

    aps_it = adj_pane_set[i][0].begin();
    for(; aps_it != adj_pane_set[i][0].end(); ++aps_it){
      std::cout << *aps_it << " ";
    }

#endif

    //   3 Identify "possibly shared faces", those faces whose nodes are all
    //     shared with the same adjacent pane. Make a list of these faces.

    bf_it2 = border_facets[i].end();
    bf_it = border_facets[i].begin();
    //for(int j=0, nj=border_facets[i].size(); j<nj; ++j){
    for(; bf_it != bf_it2; ++bf_it){
      Element_node_enumerator ene(allpanes[i], (*bf_it).eid());
      Facet_node_enumerator fne (&ene, (*bf_it).lid());
      //if all the nodes are shared
      if( is_shared[fne[0]-1] &&
	  is_shared[fne[1]-1] &&
	  is_shared[fne[2]-1] &&
	  fne.size_of_edges()>3?is_shared[fne[3]-1]:1
	  ){
	// then see if they are all shared by the same panes
	Four_tuple ns( fne[0], fne[1], fne[2], fne.size_of_edges()>3?fne[3]:-1);
        for(int k=0; k<count; ++k){
	  aps_it2 = adj_pane_set[i][k].end();
	  if(aps_it2 != adj_pane_set[i][k].find(fne[0]) &&
	     aps_it2 != adj_pane_set[i][k].find(fne[1]) &&
	     aps_it2 != adj_pane_set[i][k].find(fne[2]) &&
	     (fne.size_of_edges()>3?(adj_pane_set[i][k].find(fne[3]) != aps_it2):1)){
	    // and if they are
	    // then this facet is possibly shared with the adjacent node
	    Four_tuple ns( fne[0], fne[1], fne[2], fne.size_of_edges()>3?fne[3]:-1);
	    std::sort(&ns[0],&ns[4]);
	    maybe_shared[i][k].insert(std::make_pair(ns,(*bf_it)));
	  }
	}
      }
    }
    
#if 0
    // print out the list of possibly shared nodes:
    std::cout << "Pane " << allpanes[i]->id() << " possibly shared panes\n";
    for(int j = 0; j<count; ++j){
      std::cout << "  faces possibly shared with Pane " << adj_pane_id[i][j] << "\n";
      ms_it2 = maybe_shared[i][j].end();
      for(ms_it = maybe_shared[i][j].begin(); ms_it != ms_it2; ++ms_it){
	std::cout << "(" << (ms_it->first)[0] << " " 
		  << (ms_it->first)[1] << " "
		  << (ms_it->first)[2] << " "
		  << (ms_it->first)[3] << ")  ";	  
      }
      std::cout << "\n";
    }
#endif

    wrk_window->set_size("com_buff", allpanes[i]->id(), count, 0);
    int* my_buff;
    wrk_window->resize_array("com_buff", (const int)allpanes[i]->id(),
			     reinterpret_cast<void**>(&my_buff));

    // find the size of the pconn to create
    // 3*(pconn_offset)           // number of communicating blocks
    // + 3                        // shared node block
    // + 6*count                  // (pane id,node count,array-id for each adj. pane)*2blocks

    int my_size = 3*pconn_offset + 3 + 6*count;
    int* my_pconn;
    wrk_window->set_size("false_pconn", allpanes[i]->id(), my_size, 
			 my_size-(3+pconn_offset));
    wrk_window->resize_array("false_pconn", (const int)allpanes[i]->id(),
			     reinterpret_cast<void**>(&my_pconn));
    // Add shared node data
    if(pconn_offset){my_pconn[0]=1;++my_pconn;}
    my_pconn[0]=1; my_pconn[1]=1; my_pconn[2]=1;
    my_pconn += 3;

    // filling in send buffer and real nodes to send info.
    if(pconn_offset){my_pconn[0]=count;++my_pconn;}
    for(int j=0; j<count; ++j){
      // print communicating pane id
      my_pconn[0] = adj_pane_id[i][j];
      // print number of ints to communicate
      my_pconn[1] = 1;
      // print index-id of information
      my_pconn[2] = j+1;
      my_pconn+=3;
      // store size of information to send
      my_buff[j] = 4*maybe_shared[i][j].size();
      send_sizes[i] += my_buff[j];
    }

    // add ghost nodes to receive info.
    if(pconn_offset){my_pconn[0]=count;++my_pconn;}
    for(int j=0; j<count; ++j){
      my_pconn[0] = adj_pane_id[i][j];
      my_pconn[1] = 1;
      my_pconn[2] = j+1;
      my_pconn+=3;
    }
    
#if 0
    for ( int j=0; j<count; ++j, index+=vs[index+1]+2) {
      // We skip the panes that are not in the current window 
      while ( wrk_window->owner_rank(vs[index])<0)
	index+=vs[index+1]+2;
      // List the communicating pane, 1 item to send, and the index of that item in a 
      // local panel array
      my_pconn[3*j] = vs[index];
      my_pconn[3*j+1] = 1;
      my_pconn[3*j+2] = j+1;

      my_pconn[3*(count+j)+pconn_offset] = vs[index];
      my_pconn[3*(count+j)+1+pconn_offset] = 1;
      my_pconn[3*(count+j)+2+pconn_offset] = j+1;

      my_pconn[3*(2*count+j)+2*pconn_offset] = vs[index];
      my_pconn[3*(2*count+j)+1+2*pconn_offset] = 1;
      my_pconn[3*(2*count+j)+2+2*pconn_offset] = j+1;


      // store the size of the buffer to send across
      my_buff[j] = 4*maybe_shared[i][j].size();
      send_sizes[i] += my_buff[j];
    }
#endif
#if 0
    std::cout << "pconn_offset = " << pconn_offset << "\n";
    std::cout << "total size = " << my_size << "\n";
    std::cout << "ghost sized = " << my_size - pconn_offset << "\n\n";
#endif
  }
  wrk_window->init_done();

#if 0
  //print out pconn, for debugging.
  for(int i=0; i<total_npanes; ++i){
    std::cout << "FALSE pconn for pane " << allpanes[i]->id() << std::endl;
    COM::Attribute *my_pconn = allpanes[i]->attribute(w_false_pconn_id);
    int* ptr = (int*)my_pconn->pointer();
    for(int k=0,nk=my_pconn->size_of_items();k<nk;++k){
      std::cout << ptr[k] << " ";
    }
    std::cout << "\n\n";

    // print out size of buffer being sent across
    std::cout << "size of face list\n";
    COM::Attribute *my_bsize = allpanes[i]->attribute(w_com_buff_id);
    ptr = (int*)my_bsize->pointer();
    for(int j =0; j<my_bsize->size_of_items(); ++j)
      std::cout << ptr[j] << " ";
    std::cout << "\n\n";
  }
#endif

  // 4 Communicate the size of the local face list to send to / recieve from
  //   adjacent panes.Pane_connectivity::
  Rocmap::update_ghosts(com_buff,
			false_pconn);

  //print out pconn, for debugging.
#if 0
  for(int i=0; i<total_npanes; ++i){
    COM::Attribute *my_pconn = allpanes[i]->attribute(w_false_pconn_id);
    int* ptr = (int*)my_pconn->pointer();
    for(int k=0,nk=my_pconn->size_of_items();k<nk;++k){
      std::cout << ptr[k] << " ";
    }
    std::cout << "\n\n";
    // print out possibly shared faces
    std::cout << "possibly shared faces\n";
    for(ms_it = maybe_shared[i][0].begin(), ms_it2 = maybe_shared[i][0].end(); 
	ms_it != ms_it2; ++ms_it){
      std::cout << "(" 
		<< (ms_it->first)[0] << " "
		<< (ms_it->first)[1] << " "
		<< (ms_it->first)[2] << " "
		<< (ms_it->first)[3] << ") ";
      }
    std::cout << "\n\n";
    // print out size of buffer being sent across
    std::cout << "size of communicated face list\n";
    COM::Attribute *my_bsize = allpanes[i]->attribute(w_com_buff_id);
    ptr = (int*)my_bsize->pointer();
    for(int j =0; j <adj_pane_id[i].size(); ++j)
      std::cout << ptr[j] << " ";
    std::cout << "\n\n";
  }
#endif

  // 5 Create buffers for the incoming face list and communicate faces lists.
  // loop through panes
  for(int i=0; i<total_npanes; ++i){

    // get pane level attributes and pointers
    COM::Attribute *p_com_buff = allpanes[i]->attribute(w_com_buff_id);
    int *com_buff_ptr = (int*)p_com_buff->pointer();

    // find the number of adjacent panes
    int count = adj_pane_id[i].size();

    // find the size of data to be received and sent.
    int recv_size = 0;
    // add up the sizes of all the buffers being sent
    for(int j=0; j< count; ++j){
      recv_size += com_buff_ptr[j];
      adj_pane_recv[i][j] = com_buff_ptr[j];
    }
    int send_size = send_sizes[i];

    // find the size of the pconn to create
    // 3*(pconn_offset)           // number of communicating blocks
    // + 3                        // shared node block
    // + 4*count                  // (pane id and node count for each adj. pane)*2blocks
    // + recv_size+send_size      // number of nodes in facet lists
    int pconn_size = 3*(pconn_offset)+3+4*count+recv_size+send_size;

#if 0

    std::cout << " SIZES\nrecv_size = " << recv_size << " send_size = " << send_size
	      << " pconn size = " << pconn_size << "\n"
	      << " 3*pconn_offset = " << 3*pconn_offset
	      << " 4*count = " << 4*count << "\n\n";

#endif

    wrk_window->set_size("com_buff", allpanes[i]->id(),std::max(send_size,recv_size),0);
    wrk_window->resize_array("com_buff", (const int)allpanes[i]->id(),
			     reinterpret_cast<void**>(&com_buff_ptr));    
    int* pconn_ptr;
    wrk_window->set_size("false_pconn", allpanes[i]->id(), pconn_size, 
			 pconn_size - (3+pconn_offset));
    wrk_window->resize_array("false_pconn", (const int)allpanes[i]->id(),
			     reinterpret_cast<void**>(&pconn_ptr));

    // add shared node data
    if(pconn_offset){pconn_ptr[0]=1;++pconn_ptr;}
    pconn_ptr[0]=1; pconn_ptr[1]=1; pconn_ptr[2]=1;
    pconn_ptr += 3;
    
    // now fill in the send buffer and real nodes to send info. in the pconn
    if(pconn_offset){pconn_ptr[0]=count;++pconn_ptr;}
    int com_ind = 0; // index of next place in com_buff
    for(int j =0; j<count; ++j){
      pconn_ptr[0]=adj_pane_id[i][j];
      pconn_ptr[1]=4*maybe_shared[i][j].size();
      pconn_ptr+=2;
      ms_it2 = maybe_shared[i][j].end();
      for(ms_it = maybe_shared[i][j].begin(); ms_it != ms_it2; 
	  ++ms_it, com_ind+=4, pconn_ptr+=4){
	// Convert nodal id's to pconn based index id's and place in the comm buffer
	const Four_tuple* ft = &(ms_it->first);	
	if ((*ft)[0]==-1)
	  com_buff_ptr[com_ind] = -1;
	else
	  idm_it = lid2ind[i][j].find((*ft)[0]);
	idm_it = lid2ind[i][j].find((*ft)[1]);
	com_buff_ptr[com_ind+1] = idm_it->second;
	idm_it = lid2ind[i][j].find((*ft)[2]);
	com_buff_ptr[com_ind+2] = idm_it->second;
	idm_it = lid2ind[i][j].find((*ft)[3]);
	com_buff_ptr[com_ind+3] = idm_it->second;
	std::sort(&com_buff_ptr[com_ind],&com_buff_ptr[com_ind+4]);
	// Also fill in the pconn
	pconn_ptr[0] = com_ind+1;
	pconn_ptr[1] = com_ind+2;
	pconn_ptr[2] = com_ind+3;
	pconn_ptr[3] = com_ind+4;
      }
    }

    // now fill in the ghost nodes to receive info. in the pconn
    if(pconn_offset){pconn_ptr[0]=count;++pconn_ptr;}
    com_ind = 1;
    for(int j=0; j<count; ++j){
      pconn_ptr[0]=adj_pane_id[i][j];
      pconn_ptr[1]=adj_pane_recv[i][j];
      for(int k=0; k<adj_pane_recv[i][j]; ++k){
	pconn_ptr[k+2] = k+com_ind;	
      }   
      com_ind += adj_pane_recv[i][j];
      pconn_ptr += 2+adj_pane_recv[i][j];
    }
  }

#if 0
  //print out pconn, for debugging.
  for(int i=0; i<total_npanes; ++i){
    std::cout << "Facet list FALSE pconn for pane " << allpanes[i]->id() << std::endl;
    COM::Attribute *my_pconn = allpanes[i]->attribute(w_false_pconn_id);
    int* ptr = (int*)my_pconn->pointer();
    for(int k=0,nk=my_pconn->size_of_items();k<nk;++k){
      std::cout << ptr[k] << " ";
    }
    std::cout << "\n\n";

    std::cout << "Facet list buffer for pane " << allpanes[i]->id() << std::endl;
    COM::Attribute *my_buff = allpanes[i]->attribute(w_com_buff_id);
    ptr = (int*)my_buff->pointer();
    for(int k=0,nk=my_buff->size_of_items();k<nk;++k){
      std::cout << ptr[k] << " ";
    }
    std::cout << "\n\n";
  }
#endif

  // update
  Rocmap::update_ghosts(com_buff,
			false_pconn);

#if 0
  //print out pconn, for debugging.
  for(int i=0; i<total_npanes; ++i){
    std::cout << "Updated Facet list buffer for pane " << allpanes[i]->id() << std::endl;
    COM::Attribute *my_buff = allpanes[i]->attribute(w_com_buff_id);
    int* ptr = (int*)my_buff->pointer();
    for(int k=0,nk=my_buff->size_of_items();k<nk;++k){
      std::cout << ptr[k] << " ";
    }
    std::cout << "\n\n";
  }
#endif
  // 6 For every face received, remove any matching faces in the list of 
  //   boundary faces.
  // loop through all panes
  for(int i =0; i<total_npanes; ++i){
    COM::Attribute *my_buff = allpanes[i]->attribute(w_com_buff_id);
    int* buff_ptr = (int*)my_buff->pointer();

    int count = adj_pane_recv[i].size();
    // loop through adjacent panes
    //std::cout << "Pane " << allpanes[i]->id() << "\n";
    for(int j=0;j<count; buff_ptr += adj_pane_recv[i][j], ++j){

      //std::cout << "Converting node array\n";
      // convert all node array indices to local id
      //std::cout << adj_pane_recv[i][j] << " items being converted\n";
      for(int k=0, nk =adj_pane_recv[i][j]; k<nk; ++k){
	if(buff_ptr[k]!=-1){
	  buff_ptr[k] = ind2lid[i][j][buff_ptr[k]];
	}
      }
      //      std::cout << "Checking for matching facets with pane" << adj_pane_id[i][j] << "\n";
      // check for matching facets
      ms_it2 = maybe_shared[i][j].end();
      bf_it2 = border_facets[i].end();
      for(int k=0, nk=adj_pane_recv[i][j]; k<nk; k+=4){
	Four_tuple ns(buff_ptr[k],buff_ptr[k+1],buff_ptr[k+2],buff_ptr[k+3]);
	std::sort(&ns[0],&ns[4]);
	//std::cout << " (" << ns[0] << " " << ns[1] << " " << ns[2] << " " << ns[3] << ") ";
	ms_it = maybe_shared[i][j].find(ns);
	if(ms_it != ms_it2){
	  //std::cout << "found ";
	  Element_node_enumerator ene(allpanes[i], (ms_it->second).eid());
	  Facet_node_enumerator fne (&ene, (ms_it->second).lid());
	  //std::cout << "(" << fne[0] << " " << fne[1]
	  //	    << " " << fne[2];
	  //	  if(fne.size_of_edges()>3)
	  //  std::cout << " " << fne[3];
	  //std::cout << ")\n";
	  MAP::Facet_ID fid = ms_it->second;
	  bf_it = border_facets[i].find(fid);
	  if(bf_it != bf_it2)
	    border_facets[i].erase(fid);
	}
	else
	  ;
	  //std::cout << "not found\n";
      }
    }
  }

  //  std::cout << "On to step 7 \n";

  // 7 Mark all nodes still contained in boundary faces as physical boundary
  //   nodes.  
  for(int i =0; i < total_npanes; ++i){
    //    std::cout << "Pane " << i << " detected surface faces\n";
    bf_it2 = border_facets[i].end();
    for(bf_it = border_facets[i].begin(); bf_it!= bf_it2; ++bf_it){
      Element_node_enumerator ene(allpanes[i], (*bf_it).eid());
      Facet_node_enumerator fne (&ene, (*bf_it).lid());
      
      //std::cout << "(" << fne[0] << " " << fne[1]
      //		<< " " << fne[2];
      //  if(fne.size_of_edges()>3)
      //	std::cout << " " << fne[3];
      //std::cout << ")";
    }
    //std::cout << "\n\n";
  }

  for(int i = 0; i < total_npanes; ++i){
    COM::Attribute *p_is_surface = allpanes[i]->attribute(w_is_surface_id);
    int* surf_ptr = (int*)p_is_surface->pointer();
    // initialize surface to 0s
    for(int j=0, nj= p_is_surface->size_of_items();j<nj; ++j){
      surf_ptr[j] = 0;
    }
    bf_it2 = border_facets[i].end();
    //std::cout << "Pane " << allpanes[i]->id() << " marking surface nodes\n";
    for(bf_it = border_facets[i].begin(); bf_it!= bf_it2; ++bf_it){
      Element_node_enumerator ene(allpanes[i], (*bf_it).eid());
      Facet_node_enumerator fne (&ene, (*bf_it).lid());
      surf_ptr[fne[0]-1] = 1;
      surf_ptr[fne[1]-1] = 1;
      surf_ptr[fne[2]-1] = 1;
      //std::cout << fne[0] << " " << fne[1] << " " << fne[2] << "\n";
      if(fne.size_of_edges()>3)
	surf_ptr[fne[3]-1] = 1;
    }
    //    std::cout << "\n\n";
  }

#if 0
  for(int i = 0; i < total_npanes; ++i){
    COM::Attribute *p_pconn = allpanes[i]->attribute(COM::COM_PCONN);
    //    std::cout << "PANE " << allpanes[i]->id() << "'s real pconn size = "
    //	      << p_pconn->size_of_real_items() << " and ghost size = "
    //	      << p_pconn->size_of_items()-p_pconn->size_of_real_items() << "\n\n";

    //int* ptr = (int*)p_pconn->pointer();
    //for(int j=0, nj=p_pconn->size_of_items(); j<nj; ++j){
      //      std::cout << ptr[j] << " ";
    //}
    //    std::cout << "\n\n";
  }
#endif

  // update
  Rocmap::reduce_maxabs_on_shared_nodes(w_is_surface);
  Rocmap::update_ghosts(w_is_surface);

  // delete buffer and false pconn attributes
  wrk_window->delete_attribute(com_buff->name());
  wrk_window->delete_attribute(false_pconn->name());  

}

MOP_END_NAMESPACE






