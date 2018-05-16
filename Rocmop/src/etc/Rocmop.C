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
// $Id: Rocmop.C,v 1.3 2008/12/06 08:45:25 mtcampbe Exp $

#ifdef MESQUITE
#define USE_STD_INCLUDES 1
#define USE_C_PREFIX_INCLUDES 1
#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MesquiteError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "ShapeImprovementWrapper.hpp"

// algorithms
#include "MeanRatioQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MsqMessage.hpp"
#include "MesqPane.h"

using namespace Mesquite;
#endif

#include "Rocblas.h"
#include "Rocmop.h"
#include "Rocmap.h"
#include "roccom.h"
#include "Pane_communicator.h"
#include "Pane_connectivity.h"
#include "Geometric_Metrics_3.h"
#include "PN_patch.h"
#include "geometry.h"
#include "mapbasic.h"

// for debugging
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
// end for debugging

MOP_BEGIN_NAMESPACE

using MAP::Pane_connectivity;
using MAP::Pane_communicator;
using MAP::Rocmap;
using std::cout;
using std::endl;

Rocmop::~Rocmop() { 
  if(_verb)
    cout << "Entering Rocmop::~Rocmop\n";

  if (_wrk_window) {
    delete _wrk_window; _wrk_window = NULL;
  }

  if(_verb > 1)
    cout << "Exiting Rocmop::~Rocmop\n";
}

void Rocmop::load( const std::string &mname) {

  Rocmop *mop = new Rocmop();

  COM_new_window( mname.c_str());

  std::string glb=mname+".global";

  COM_new_attribute( glb.c_str(), 'w', COM_VOID, 1, "");
  COM_set_object( glb.c_str(), 0, mop);

  COM_Type types[3];
  types[0] = COM_RAWDATA; types[1] = COM_METADATA;
  types[2] = COM_METADATA;

  COM_set_member_function( (mname+".smooth").c_str(), 
			   (Member_func_ptr)(&Rocmop::smooth),
			   glb.c_str(), "biB", types);

  COM_set_member_function( (mname+".smooth_in_place").c_str(), 
			   (Member_func_ptr)(&Rocmop::smooth_in_place),
			   glb.c_str(), "bb", types);

  types[1] = COM_STRING; types[2] = COM_VOID;
  COM_set_member_function( (mname+".set_value").c_str(), 
			   (Member_func_ptr)(&Rocmop::set_value),
			   glb.c_str(), "bii", types);    
  
  COM_window_init_done( mname.c_str());

}

void Rocmop::unload( const std::string &mname) {

  Rocmop *mop;
  std::string glb=mname+".global";

  COM_get_object( glb.c_str(), 0, &mop);
  delete mop;

  COM_delete_window( mname.c_str());

}

void Rocmop::smooth(const COM::Attribute *pmesh,
		    COM::Attribute *disp){
  COM_assertion_msg( validate_object()==0, "Invalid object");
  if(_verb)
    cout << "Entering Rocmop::smooth\n";

  int pmesh_id = pmesh->id();
  COM_assertion_msg( pmesh &&
		     (pmesh_id==COM::COM_PMESH || pmesh_id==COM::COM_MESH),
		     "Input to Rocmop::smooth must be a pmesh.");
  _is_pmesh = (pmesh_id == COM_PMESH) ? true : false;

  if(_wrk_window) delete _wrk_window;
  _usr_window = pmesh->window();

  // Create a buffer window by inheriting(clone) from the user's mesh.
  if (_verb >3) 
    std::cout << "  Creating(clone) buffer window." << std::endl;
  std::string buf_name(_usr_window->name()+"-Rocmopbuf");
  _wrk_window = new COM::Window(buf_name, 
				_usr_window->get_communicator());
  _wrk_window->inherit( const_cast<COM::Attribute*>(pmesh), "", 
			COM::Pane::INHERIT_CLONE, true, NULL, 0);
  _wrk_window->init_done();

  std::vector<const Pane*> allpanes;
  _wrk_window->panes(allpanes);

  double pre_worst = 180.0, post_worst = 0.0;

  // Check initial quality if we are in lazy or monotone mode
  if(_monotone || _lazy){
    pre_worst = check_all_elem_quality(allpanes);
    // get the worst pre smoothing quality across all panes
    agree_double(pre_worst, MPI_MAX);
  }

  bool to_smooth = (!_lazy || (pre_worst > 180.*_tol));
    
  int degraded = 0;
  if (to_smooth){
    perform_smoothing(pre_worst);
    
    // Check final mesh quality if a tolerance is set, or we're in monotonic mode
    if(_monotone || (_tol!=0.0)){
      post_worst = check_all_elem_quality(allpanes);
      agree_double(post_worst, MPI_MAX);
    }
    
    // If in monotone mode, see if quality degraded on any pane
    if(_monotone && (pre_worst > post_worst)){
      cerr << "Warning, smoothing would have degraded mesh quality \n" 
	   << "so the mesh has not been modified.\n";
      degraded = 1;
    }
    
    // Warn the user if could not reach a specified minimum mesh quality.
    if(_tol && (post_worst < _tol)){
      cerr << "Warning, post-smoothing mesh quality is below specified minimum " 
	   << "quality \n";
    }
  }
  
  const COM::Attribute *old_nc = pmesh->window()->attribute( COM::COM_NC);
  const COM::Attribute *new_nc = (degraded || !to_smooth) ? 
    old_nc :_wrk_window->attribute( COM::COM_NC);
  Rocblas::sub (new_nc, old_nc, disp);

  if (_verb > 1) 
    std::cout << "Exiting rocmop::smooth" << std::endl;
}

void Rocmop::smooth_in_place(COM::Attribute *pmesh){
  COM_assertion_msg( validate_object()==0, "Invalid object");
  if (_verb)
    std::cout << "Entering rocmop::smooth_in_place" << std::endl;

  int pmesh_id = pmesh->id();

  COM_assertion_msg( pmesh && 
		     (pmesh_id==COM::COM_PMESH || pmesh_id==COM::COM_MESH) , 
		     "Input to Rocmop::smooth_in_place must be a mesh or pmesh");

  if(_wrk_window) 
    delete _wrk_window;
  _usr_window = pmesh->window();

  // FIXME
  // Create a buffer window by inheriting(use) from the user's mesh.
  // can't use the _usr_window directly because of problems with
  // Element_node_enumerator
  _wrk_window = new COM::Window(_usr_window->name()+"-Rocmopbuf", 
				_usr_window->get_communicator());
  _wrk_window->inherit( const_cast<COM::Attribute*>(pmesh), "", 
			false, true, NULL, 0);
  _wrk_window->init_done();

  std::vector<const Pane*> allpanes;
  _wrk_window->panes(allpanes);

  double pre_worst = 180.0, post_worst = 0.0;

  // Check initial quality if we are in lazy or monotone mode
  if(_monotone || _lazy) {
    pre_worst = check_all_elem_quality(allpanes);
    // get the worst pre smoothing quality on all panes
    agree_double(pre_worst, MPI_MIN);
  }

  bool to_smooth = (!_lazy || (pre_worst > 180.*_tol));

  if (to_smooth){
    perform_smoothing(pre_worst);
    
    // Check final mesh quality if a tolerance is set, or we're in monotonic mode
    if(_monotone || (_tol!=0.0)){
      post_worst = check_all_elem_quality(allpanes);
      agree_double(post_worst, MPI_MIN);
    }
    
    // If in monotone mode, see if quality degraded on any pane
    if(_monotone && (pre_worst > post_worst))
      cerr << "Warning, mesh quality degraded during smoothing in place." << endl;

    // If a quality tolerance is set, and we didn't meet it, then warn
    if(_tol && (post_worst < _tol)){
      cerr << "Warning, post-smoothing mesh quality is below specified minimum " 
	   << "quality \n";
    }
  }

  if (_wrk_window){ delete _wrk_window;_wrk_window = NULL;}

  if (_verb > 1)
    std::cout << "Exiting rocmop::smooth_in_place" << std::endl;
}

void Rocmop::perform_smoothing(double pre_quality){
  if(_verb) 
    std::cout << "  Entering Rocmop::perform_smoothing" << std::endl;

  COM_assertion_msg(_wrk_window, "Unexpected NULL pointer encountered.");

  // Perform initialization specific to the smoother.
  smoother_specific_init();

  // Select iterative or noniterative option depending on the smoother.
  if (0) // No noniterative smoothers currently
    perform_noniterative_smoothing();
  else if( _method < SMOOTH_NONE)
    perform_iterative_smoothing(pre_quality);
  else
    COM_assertion_msg(0, "No valid smoothing method selected");

  if(_verb) 
    std::cout << "  Exiting Rocmop::perform_smoothing" << std::endl;
}

void Rocmop::perform_iterative_smoothing(double pre_quality){
  if(_verb) 
    std::cout << "    Entering Rocmop::perform_iterative_smoothing" << std::endl;

  COM_assertion_msg(_wrk_window, "Unexpected NULL pointer encountered.");

  // bool_iter is true until the max number of iterations (_niter)
  // is reached or the minimum quality (_tol) is reached.
  int to_iter = true; 
  std::vector<const Pane*> allpanes;
  _wrk_window->panes(allpanes);
  double cur_qual = pre_quality;
  
  for(int i = 0; 
      ( (i < _niter) && (to_iter) ); 
      ++i, agree_int(to_iter, MPI_MIN)){

    if(_verb > 2)
      std::cout << "      Smoothing iteration " << i << std::endl;

#ifdef MESQUITE
    if(_method==SMOOTH_VOL_MESQ_WG){
      //std::string msg("\nQuality prior to smoothing = ");
      //print_quality(msg);
      smooth_vol_mesq_wg(cur_qual);
    }
    //    else if(_method == SMOOTH_LAPLACE)
    //  smooth_laplace();
    else if(_method == SMOOTH_VOL_MESQ_NG){
      //std::string msg("Quality prior to smoothing = ");
      //print_quality(msg);
      smooth_vol_mesq_ng(cur_qual);
    }
#else
    if((_method==SMOOTH_VOL_MESQ_WG) || (_method==SMOOTH_VOL_MESQ_NG))
      COM_assertion_msg(0,"Rocmop not compiled with MESQUITE");
#endif
    else if(_method==SMOOTH_SURF_MEDIAL)
      smooth_surf_medial();
    else COM_assertion_msg(0, "No valid iterative smoothing method selected");  

    // If a non zero quality tolerance is set, then determine if quality
    // is low enough to warrant smoothing.
    if( _tol != 0.0 ){
      cur_qual = check_all_elem_quality(allpanes);
      if(cur_qual >= _tol)
	to_iter = false;
    }
  }

  if(_verb > 1) 
    std::cout << "    Exiting Rocmop::perform_iterative_smoothing" << std::endl;
}

void Rocmop::perform_noniterative_smoothing(){
  if(_verb) 
    std::cout << "    Entering Rocmop::perform_noniterative_smoothing" << std::endl;

  if(_niter != 1){
    std::cerr << "Although the maximum number of iterations > 1 is selected,\n"
	      << "the smoothing method is noniterative, and will run once.\n";
  }

  if (0)
    ;// Currently, no noniterative methods exist.
  else COM_assertion_msg(0, "No valid noniterative smoothing method selected");

  if(_verb > 1) 
    std::cout << "    Exiting Rocmop::perform_noniterative_smoothing" << std::endl;
}

// Perform smoother specific initializing, for example initializing the
// Window_manifold_2 for a surface mesh, adding smoother specific attributes
// to the window, etc.
void Rocmop::smoother_specific_init(){
  if(_verb) 
    std::cout << "    Entering Rocmop::smoother_specific_init" << std::endl;

  // get rid of any old data
  if(_wm){ delete _wm; _wm = NULL; }

  // perform initialization common to surface smoothers
  if(_method==SMOOTH_SURF_MEDIAL){
    // Initialize the Window_manifold
    if(_wm == NULL)
      Rocsurf::initialize(_wrk_window->attribute(_is_pmesh ? COM_PMESH : COM_MESH));    

  }
  
  // perform smoother specific initialization 
  switch (_method){
    
#ifdef MESQUITE
  case SMOOTH_VOL_MESQ_WG: {
    // Obtain a list of elements containing shared nodes for each pane.
    // Don't bother if _ctol == 0 or _ncycle ==1
    //    if(!(_ctol==0.) && !(_ncycle==1))
    determine_shared_border();
    if(_invert)
      invert_tets();
    break;
  }
  case SMOOTH_VOL_MESQ_NG: {

    // Check to see if the physical surface boundary exists
    const std::string surf_attr("is_surface");
    const COM::Attribute *w_is_surface = _usr_window->attribute(surf_attr);
    if(w_is_surface){
      
      COM_assertion_msg( COM_compatible_types( w_is_surface->data_type(), COM_INT),
			 "Surface-list must have integer type");
      COM_assertion_msg( w_is_surface->is_nodal() == 1,
			 "Surface-list must be nodal");
      COM_assertion_msg( w_is_surface->size_of_components() == 1,
			 "Surface-list must have a single component");
      COM_assertion_msg( w_is_surface->initialized() == 1,
			 "Surface-list must be initialized");      

      // Clone the attribute
      COM::Attribute * new_attr = 
	_wrk_window->inherit( const_cast<COM::Attribute *>(w_is_surface), 
			      surf_attr, COM::Pane::INHERIT_CLONE, true, NULL, 0);
    }
    // else, detect the physical boundary ourselves
    else{
      COM::Attribute* w_surf_attr =
	_wrk_window->new_attribute( "is_surface", 'n', COM_INT, 1, "");
      _wrk_window->resize_array( w_surf_attr, 0);
      
      determine_physical_border(w_surf_attr);
    }
    _wrk_window->init_done();
    
    if(_invert)
      invert_tets();
    break;
  }
#else
  case SMOOTH_VOL_MESQ_WG:
  case SMOOTH_VOL_MESQ_NG:
    COM_assertion_msg(0, "Not compiled with MESQUITE");
    break;
#endif

  case SMOOTH_SURF_MEDIAL: {

    // Extend buffer window
    COM::Attribute* w_disps =
      _wrk_window->new_attribute( "disps", 'n', COM_DOUBLE, 3, "");
    _wrk_window->resize_array( w_disps, 0);

    COM::Attribute* w_facenormals = 
      _wrk_window->new_attribute( "facenormals", 'e', COM_DOUBLE, 3, "");
    _wrk_window->resize_array( w_facenormals, 0);
    
    COM::Attribute* w_facecenters = 
      _wrk_window->new_attribute( "facecenters", 'e', COM_DOUBLE, 3, "");
    _wrk_window->resize_array( w_facecenters, 0);

    COM::Attribute* w_eigvalues = 
    _wrk_window->new_attribute( "lambda", 'n', COM_DOUBLE, 3, "");
    _wrk_window->resize_array( w_eigvalues, 0);

    COM::Attribute* w_vnormals = 
    _wrk_window->new_attribute( "vnormals", 'n', COM_DOUBLE, 3, "");
    _wrk_window->resize_array( w_vnormals, 0);

    COM::Attribute* w_awnormals = 
    _wrk_window->new_attribute( "awnormals", 'n', COM_DOUBLE, 3, "");
    _wrk_window->resize_array( w_awnormals, 0);

    COM::Attribute* w_uwnormals = 
    _wrk_window->new_attribute( "uwnormals", 'n', COM_DOUBLE, 3, "");
    _wrk_window->resize_array( w_uwnormals, 0);

    COM::Attribute* w_eigvecs = 
    _wrk_window->new_attribute( "eigvecs", 'n', COM_DOUBLE, 9, "");
    _wrk_window->resize_array( w_eigvecs, 0);

    COM::Attribute* w_tangranks = 
    _wrk_window->new_attribute( "tangranks", 'n', COM_INT, 1, "");
    _wrk_window->resize_array( w_tangranks, 0);

    COM::Attribute* w_cntnranks = 
    _wrk_window->new_attribute( "cntnranks", 'n', COM_INT, 1, "");
    _wrk_window->resize_array( w_cntnranks, 0);
    
    COM::Attribute* w_cntnvecs =
    _wrk_window->new_attribute( "cntnvecs", 'n', COM_DOUBLE, 6, "");
    _wrk_window->resize_array( w_cntnvecs, 0);

    COM::Attribute* w_scales = 
    _wrk_window->new_attribute( "scales", 'n', COM_DOUBLE, 1, "");
    _wrk_window->resize_array( w_scales, 0);
    
    COM::Attribute* w_weights = 
    _wrk_window->new_attribute( "weights", 'n', COM_DOUBLE, 1, "");
    _wrk_window->resize_array( w_weights, 0);

    COM::Attribute* w_weights2 = 
    _wrk_window->new_attribute( "weights2", 'n', COM_DOUBLE, 1, "");
    _wrk_window->resize_array( w_weights2, 0);

    COM::Attribute* w_barycrds = 
    _wrk_window->new_attribute( "barycrds", 'n', COM_DOUBLE, 2, "");
    _wrk_window->resize_array( w_barycrds, 0);

    COM::Attribute* w_PNelemids = 
    _wrk_window->new_attribute( "PNelemids", 'n', COM_INT, 1, "");
    _wrk_window->resize_array( w_PNelemids, 0);

    // Extend the buffer window to hold local contributions to new placement
    // and the number of contributing faces for Laplacian smoothing.
    COM::Attribute * w_pnt_contrib = 
      _wrk_window->new_attribute("pnt_contrib", 'n', COM_DOUBLE, 3, "");
    _wrk_window->resize_array(w_pnt_contrib, 0);

    COM::Attribute * w_disp_count =
      _wrk_window->new_attribute("disp_count", 'n', COM_DOUBLE, 1, "");
    _wrk_window->resize_array(w_disp_count, 0);

    _wrk_window->init_done();

    break;
  }
    //  case SMOOTH_LAPLACE : {
    //break;
    //}
  default:
    COM_assertion_msg(0, "Can't initialize for invalid smoother.");
    break;
  }

  COM_assertion_msg(_wrk_window, "Unexpected NULL pointer encountered.");

  if(_verb > 1) 
    std::cout << "    Exiting Rocmop::smoother_specific_init" << std::endl;
}

void Rocmop::set_value(const char* opt, const void* value)
{
  if(_verb) 
    std::cout << "Entering Rocmop::set_value" << std::endl;

  COM_assertion_msg( validate_object()==0, "Invalid object");

  COM_assertion_msg( opt && value,
		      "Rocmop::set_value does not except NULL parameters");
  std::string option;
  if(opt) option = opt;
  if ( option == "method") {
    COM_assertion_msg( *((int*)value) <= SMOOTH_NONE && *((int*)value)>=0 
                      ,"Illegal value for 'method' option");
    _method = *((int*)value);
  }
  else if ( option == "verbose"){ 
    _verb = *((int*)value); }
  else if ( option == "lazy"){ 
    _lazy = *((int*)value); }
  else if ( option == "tol"){ 
    COM_assertion_msg( *((float*)value) <= 180. && *((float*)value)>=0. 
                      ,"Illegal value for 'method' option");
    _tol = *((float*)value); }
  else if ( option == "niter"){ 
    _niter = *((int*)value); }
  else if ( option == "ctol"){ 
    COM_assertion_msg( *((float*)value) <= 1. && *((float*)value)>=0. 
                      ,"Illegal value for 'method' option");
    _ctol = *((float*)value); }
  else if ( option == "ncycle"){
    _ncycle = *((int*)value); }
  else if ( option == "inverted"){
    _invert = *((int*)value); 
  }
  else COM_assertion_msg( false, "Unknown option");

  if(_verb > 1) 
    std::cout << "Exiting Rocmop::set_value" << std::endl;
}

void Rocmop::smooth_mesquite(std::vector<COM::Pane*> &allpanes,
			     int ghost_level){

  Mesquite::MsqError err;
  MesqPane *mp;
  ShapeImprovementWrapper mesh_quality_algorithm;

  int total_npanes = (int)allpanes.size();
  bool wg = (ghost_level == 0) ? false : true;
  for(int j=0; j < total_npanes; ++j){
    MeshSet* mesh_set1 = new MeshSet;
    mp = new MesqPane(allpanes[j], wg);
    if(_verb > 4) 
      mp->set_verb(_verb - 4);
    mesh_set1->add_mesh(mp, err); 
    MSQ_CHKERR(err);
    mesh_quality_algorithm.run_instructions(*mesh_set1, err); 
    MSQ_CHKERR(err);
    delete mesh_set1;
    if(mp)
      delete mp;
    mp = NULL;
  }
}

void Rocmop::reduce_sum_on_shared_nodes(COM::Attribute *att){
  Pane_communicator pc(att->window(), att->window()->get_communicator());
  pc.init(att);
  pc.begin_update_shared_nodes();
  pc.reduce_on_shared_nodes(MPI_SUM);
  pc.end_update_shared_nodes();
}

void Rocmop::determine_pane_border(){
  if(_verb) 
    std::cout << "Entering Rocmop::determine_pane_border" << std::endl;
  
  std::vector<const COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();

  _is_pane_bnd_node.resize(local_npanes);
  _is_pane_bnd_elem.resize(local_npanes);

  for(int i=0; i< local_npanes; ++i){
    int size_of_real_nodes = allpanes[i]->size_of_real_nodes();
    int size_of_real_elems = allpanes[i]->size_of_real_elements();
    _is_pane_bnd_node[i].resize(size_of_real_nodes,0);

    std::vector<bool> is_isolated; // is a node isolated?
    MAP::Pane_boundary pb (allpanes[i]);
    pb.determine_border_nodes(_is_pane_bnd_node[i], is_isolated);
  }

  mark_elems_from_nodes(_is_pane_bnd_node,_is_pane_bnd_elem);
  
  if(_verb > 1) 
    std::cout << "Exiting Rocmop::determine_pane_border" << std::endl;    
}


void Rocmop::determine_shared_border(){
  if(_verb) 
    std::cout << "Entering Rocmop::determine_shared_nodes" << std::endl;

  std::vector<const COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();
    
  _is_shared_node.resize(local_npanes);
  
  //First, get the list of shared nodes.
  for(int i=0; i < (int)(local_npanes); ++i){
    // Obtain the pane connectivity of the local pane.
    const COM::Attribute *pconn = allpanes[i]->attribute(COM::COM_PCONN);
    // Use the pconn offset
    const int *vs = (const int*)pconn->pointer()+Pane_connectivity::pconn_offset();
    int vs_size=pconn->size_of_real_items()-Pane_connectivity::pconn_offset();    
    _is_shared_node[i].resize(allpanes[i]->size_of_real_nodes(),0);
    
    // Determine the number of communicating panes for shared nodes.
    int count=0;
    for (int j=0, nj=vs_size; j<nj; j+=vs[j+1]+2) {
      if (_wrk_window->owner_rank( vs[j]) >=0) ++count;
    }
    
    int index = 0;
    // Loop through communicating panes for shared nodes.
    for ( int j=0; j<count; ++j, index+=vs[index+1]+2) {
      // We skip the panes that are not in the current window 
      while ( _wrk_window->owner_rank(vs[index])<0) {
	index+=vs[index+1]+2;
	COM_assertion_msg( index<=vs_size, "Invalid communication map");
      }	
      // Add shared nodes to the list
      for(int k=0; k<vs[index+1]; ++k){
	_is_shared_node[i][vs[index+2+k]-1] = 1;
      }
    }
  }

  mark_elems_from_nodes(_is_shared_node,_is_shared_elem);

  if(_verb > 1) 
    std::cout << "Exiting Rocmop::determine_shared_nodes" << std::endl;
}

void Rocmop::determine_physical_border(){
  if(_verb) 
    std::cout << "Entering Rocmop::determine_physical_border()" << std::endl;

  const std::string surf_attr("is_surface");
  COM::Attribute* w_is_surface = _wrk_window->attribute(surf_attr);
  COM_assertion_msg( w_is_surface, "Unexpected NULL pointer");
  int is_surface_id = w_is_surface->id();

  std::vector<const COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();
    
  _is_pane_bnd_node.resize(local_npanes);

  for(int i=0; i < local_npanes; ++i){
    _is_pane_bnd_node[i].resize(allpanes[i]->size_of_real_nodes());

    // get pane level pointer to physical border property.
    const COM::Attribute *p_is_surface = allpanes[i]->attribute(is_surface_id);
    int *is_surface_ptr = (int*)p_is_surface->pointer();
    
    // loop through real nodes
    for(int j=0; j< allpanes[i]->size_of_real_nodes(); ++j){
      if (is_surface_ptr[j])
	_is_pane_bnd_node[i][j] = true;
    }
  }

  mark_elems_from_nodes(_is_pane_bnd_node,_is_pane_bnd_elem);

  if(_verb > 1)
    std::cout << "Exiting Rocmop::determine_physical_border()" << std::endl;
}

void Rocmop::mark_elems_from_nodes(std::vector<std::vector<bool> > &marked_nodes,
				   std::vector<std::vector<bool> > &marked_elems){
  if(_verb) 
    std::cout << "Entering Rocmop::mark_elems_from_nodes" << std::endl;

  std::vector<const COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();

  marked_elems.clear();
  marked_elems.resize(local_npanes);
  
  //Loop through panes
  for(int i=0; i < (int)(local_npanes); ++i){

    marked_elems[i].clear();
    marked_elems[i].resize(allpanes[i]->size_of_real_elements(),false);

    // Loop through real elements.
    // Mark for quality check if they contain shared nodes.
    int s_real_elems = allpanes[i]->size_of_real_elements();
    std::vector<int> nodes;
    for(int j=1; j<= s_real_elems; ++j){
      Element_node_enumerator ene(allpanes[i],j);
      ene.get_nodes(nodes);
      for(int k=0, nk=nodes.size(); k<nk; ++k){
	if (marked_nodes[i][nodes[k]-1])
	  marked_elems[i][j-1] = true;
      }
    }
  }

  if(_verb > 1) 
    std::cout << "Exiting Rocmop::mark_elems_from_nodes" << std::endl;  
}

void Rocmop::invert_tets(){
  std::vector<Pane*> allpanes;
  _wrk_window->panes(allpanes);
  for(int i=0, ni = allpanes.size(); i<ni; ++i){
    MesqPane* mp = new MesqPane(allpanes[i]);
    mp->invert();
    if(mp)
      delete mp;
    mp = NULL;
  }  
}

double Rocmop::check_marked_elem_quality(std::vector<std::vector<bool> > &marked_elems,
					 std::vector<COM::Pane*> &allpanes){
  if(_verb) 
    std::cout << "Entering Rocmop::check_marked_elems" << std::endl;

  double worst_angle = 0.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk=allpanes[i]->size_of_real_elements(); k<nk; ++k){
      if(marked_elems[i][k]){
	Element_node_enumerator ene(allpanes[i],k+1);
	Angle_Metric_3 am;	
	am.initialize(ene);
	am.compute(angles);
	if(angles[1]>worst_angle)
	  worst_angle = angles[1];
      }
    }
  }
  return worst_angle;

  if(_verb > 1) 
    std::cout << "Exiting Rocmop::check_marked_elems" << std::endl;
}

double Rocmop::check_all_elem_quality(std::vector<const COM::Pane*> &allpanes,
				      bool with_ghosts){
  if(_verb) 
    std::cout << "Exiting Rocmop::check_all_elem_quality" << std::endl;

  int rank =0;
  int ierr = MPI_Comm_rank( _wrk_window->get_communicator(),
			    &rank); assert( ierr == 0);
  double worst_angle = 0.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    int nk=allpanes[i]->size_of_real_elements();
    if(with_ghosts)
      nk = allpanes[i]->size_of_elements();
    for(int k =0; k<nk; ++k){
      Element_node_enumerator ene(allpanes[i],k+1);
      Angle_Metric_3 am;
      am.initialize(ene);
      am.compute(angles);
      if(angles[1]>worst_angle)
	worst_angle = angles[1];
    }
  }
  return worst_angle;

  if(_verb > 1) 
    std::cout << "Exiting Rocmop::check_all_elem_quality" << std::endl;
}

double Rocmop::check_all_elem_quality(std::vector<COM::Pane*> &allpanes,
				      bool with_ghosts){
  if(_verb) 
    std::cout << "Entering Rocmop::check_all_elem_quality" << std::endl;

  int rank =0;
  int ierr = MPI_Comm_rank( _wrk_window->get_communicator(),
			    &rank); assert( ierr == 0);

  double worst_angle = 0.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    int nk=allpanes[i]->size_of_real_elements();
    if(with_ghosts){
      nk = allpanes[i]->size_of_elements();
    }
    for(int k =0; k<nk; ++k){
      Element_node_enumerator ene(allpanes[i],k+1);
      Angle_Metric_3 am;
      am.initialize(ene);
      am.compute(angles);

      if(angles[1]>worst_angle)
	worst_angle = angles[1];
    }
  }
  return worst_angle;

  if(_verb > 1) 
    std::cout << "Exiting Rocmop::check_all_elem_quality" << std::endl;
}

void Rocmop::print_quality(std::string &s){

  string outstr("angles.txt");
  ofstream file (outstr.c_str(), std::ios::app);
  COM_assertion_msg( file, "File failed to open\n");
  
  std::vector<Pane*> allpanes;
  _wrk_window->panes(allpanes);

  double max_angle = 0.0;
  double min_angle = 180.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk=allpanes[i]->size_of_real_elements(); k<nk; ++k){
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
  int rank =0;
  if(COMMPI_Initialized()){
    agree_double(max_angle,MPI_MAX);
    agree_double(min_angle,MPI_MIN);
    int ierr = MPI_Comm_rank( _wrk_window->get_communicator(),
			      &rank); assert( ierr == 0);
  }
  if (rank==0){
    file << std::left << std::setw(30) << s << std::setw(0)
	 << "(" << min_angle << " , " << max_angle << ")\n";
  }
}

void Rocmop::print_mquality(std::string &s,
			    std::vector<std::vector<bool> > &to_check){

  std::vector<std::vector<bool> > elem_to_check;
  mark_elems_from_nodes(to_check, elem_to_check);

  std::vector<Pane*> allpanes;
  _wrk_window->panes(allpanes);
  
  double max_angle = 0.0;
  double min_angle = 180.0;
  double angles[] = {0.0, 0.0};
  int id = -1;
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk = elem_to_check[i].size(); k<nk; ++k){
      if(elem_to_check[i][k]){
	Element_node_enumerator ene(allpanes[i],k+1);
	Angle_Metric_3 am;
	am.initialize(ene);
	am.compute(angles);
	if(angles[1]>max_angle)
	  max_angle = angles[1];
	if(angles[0]<min_angle){
	  id = k;
	  min_angle = angles[0];
	}
      }
    }
  }
  int rank =0;
  double temp = min_angle;
  if(COMMPI_Initialized()){
    agree_double(max_angle,MPI_MAX);
    agree_double(min_angle,MPI_MIN);
    int ierr = MPI_Comm_rank( _wrk_window->get_communicator()
			      , &rank); 
    assert( ierr == 0);
  }  

  if (rank==0){
    string outstr("angles.txt");
    ofstream file (outstr.c_str(), std::ios::app);
    COM_assertion_msg( file, "File failed to open\n");

    file << std::left << std::setw(30) << s << std::setw(0)
	 << "(" << min_angle << " , " << max_angle << ")";

    file.close();
  }

  if(COMMPI_Initialized())
    int ierr = MPI_Barrier(_wrk_window->get_communicator());

  if(min_angle == temp){
    
    string outstr("angles.txt");
    ofstream file (outstr.c_str(), std::ios::app);
    COM_assertion_msg( file, "File failed to open\n");
    
    file << "  worst = (" << rank << " , " << id << ")\n";
    file.close();
  }

  if(COMMPI_Initialized())
    int ierr = MPI_Barrier(_wrk_window->get_communicator());
}

extern "C" void Rocmop_load_module( const char *mname) 
{ Rocmop::load( mname); }

extern "C" void Rocmop_unload_module( const char *mname) 
{ Rocmop::unload( mname); }

// Fortran bindings
extern "C" void COM_F_FUNC2(rocmop_load_module, ROCMOP_LOAD_MODULE)( const char *mname, long int length) 
{ Rocmop::load( std::string(mname, length)); }

extern "C" void COM_F_FUNC2(rocmop_unload_module, ROCMOP_UNLOAD_MODULE)( const char *mname, long int length) 
{ Rocmop::unload( std::string(mname, length)); }

MOP_END_NAMESPACE






