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
// $Id: Rocprop.C,v 1.27 2008/12/06 08:45:28 mtcampbe Exp $

#include "Rocprop.h"
#include "Rocblas.h"

#include "MarkerParticles_3.h"
#include "FaceOffset_3.h"

PROP_BEGIN_NAMESPACE

using namespace SURF;

Rocprop::~Rocprop() {
  // Delete Window_manifold_2, Window, and Propagation_3 objects.
  if ( _wm) { delete _wm; _wm = NULL; }
  if ( _buf){ delete _buf; _buf = NULL; }
  
  if ( _prop) { delete _prop; _prop = NULL; }
  if ( _rem && _rem_owner) { delete _rem; _rem = NULL; }
}

void Rocprop::initialize( const COM::DataItem *pmesh, Rocsurf *rsurf)
{
  COM_assertion_msg( validate_object()==0, "Invalid object");
  COM_assertion_msg( pmesh, "Mesh must be present");

  _win = const_cast<COM::Window*>(pmesh->window());
  _parent = rsurf;

  MPI_Comm comm = _win->get_communicator();
  if ( COMMPI_Initialized())  _rank = COMMPI_Comm_rank(comm);
  else _rank = 0;

  if ( _parent==NULL) {
    // Initialize the superclass
    Rocsurf::initialize( pmesh);
  }
  else {   
    // If has a parent, then use the Window_manifold_2 of parent
    // and clear up this object's Window_manifold_2 if created.
    if ( _wm) { delete _wm; _wm = NULL; }
    COM_assertion_msg(rsurf->manifold(), "Manifold of parent must exist");
  }

  // Delete propagation object
  if ( _prop) { delete _prop; _prop = NULL; }
  
  // Create a buffer window to clone the coordinates but use nc and pconn.
  if ( _buf) delete _buf;
  _buf = new COM::Window( _win->name()+"-propbuffer", comm);
  _buf->inherit( _win->dataitem( COM::COM_CONN), "", 
		 COM::Pane::INHERIT_USE, true, NULL, 0);
  _buf->inherit( _win->dataitem( COM::COM_NC), "", 
		 COM::Pane::INHERIT_USE, true, NULL, 0);
  _buf->inherit( _win->dataitem( COM::COM_NC), "oldnc", 
		 COM::Pane::INHERIT_CLONE, true, NULL, 0);
  _buf->inherit( _win->dataitem( COM::COM_PCONN), "", 
		 COM::Pane::INHERIT_USE, true, NULL, 0);

  _buf->init_done();
}

void Rocprop::perturb_mesh( COM::DataItem *pmesh, const double &alpha) {
  COM_assertion_msg( validate_object()==0, "Invalid object");

  // If Rocprop is not yet initialized, call the initialization routine.
  SURF::Window_manifold_2 *wm = manifold();
  if ( wm == NULL) 
  { initialize( pmesh, _parent); wm = manifold(); }
  
  wm->perturb_mesh( alpha);
}

void Rocprop::set_constraints( const COM::DataItem *cnstr_types) { 
  COM_assertion_msg( validate_object()==0, "Invalid object");
  _cnstr_types = cnstr_types; 
}

void Rocprop::set_bounds( const COM::DataItem *bnd) {
  COM_assertion_msg( validate_object()==0, "Invalid object");
  _cnstr_bound = bnd; 
}

void Rocprop::propagate( const COM::DataItem *pmesh,
			 COM::DataItem *spds_io,
			 const double *dt, 
			 COM::DataItem *du,
			 double *dt_elapsed,
			 int *code)
{
  COM_assertion_msg(false,"Made it here.");
  COM_assertion_msg( validate_object()==0, "Invalid object");

  // If Rocprop is not yet initialized, call the initialization routine.
  SURF::Window_manifold_2 *wm = manifold();
  if ( wm == NULL) 
  { initialize( pmesh, _parent); wm = manifold(); }
  
  COM_assertion_msg( pmesh->window() == _win,
		     "Rocprop was initialized for a different window");

  // If the propagation object is not yet created, create a default one.
  if ( _prop == NULL) {
    if ( _prop_method == PROP_FO)
      _prop = new FaceOffset_3( wm, _buf);
    else
      _prop = new MarkerParticles_3( wm, _buf);
  }

  // Save the coordinates into buffers
  COM::DataItem *nc = _buf->dataitem( COM::COM_NC);
  COM::DataItem *oldnc = _buf->dataitem( "oldnc");
  Rocblas::copy( nc, oldnc);

  _prop->set_constraints( _cnstr_types);

  // std::cout << "made it HERE" << std::endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // std::cout << "made it HERE" << std::endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // std::cout << "made it HERE" << std::endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // std::cout << "made it HERE" << std::endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // std::cout << "made it HERE" << std::endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // std::cout << "made it HERE" << std::endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // exit(1);

  // Initialize constraints
  _prop->set_bounds( _cnstr_bound);
  _prop->set_verbose(_verb);
  
  // std::cout << "made it HERE" << std::endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // exit(1);



  // Set tolerance for eigenvalues.
  if ( _prop_method == PROP_FO) {
    if ( _eig_thres>=0)
      ((FaceOffset_3*)_prop)->set_eigen_threshold( _eig_thres);
    if ( _courant>=0)
      ((FaceOffset_3*)_prop)->set_courant_constant( _courant);
    if ( _fangle_strong>=0)
      ((FaceOffset_3*)_prop)->set_fangle_strong( _fangle_strong);
    if ( _fangle_weak>=0)
      ((FaceOffset_3*)_prop)->set_fangle_weak( _fangle_weak);
    if ( _fangle_turn>=0)
      ((FaceOffset_3*)_prop)->set_fangle_turn( _fangle_turn);
    if ( _wf_expn>=0)
      ((FaceOffset_3*)_prop)->set_wavefrontal_expansion( _wf_expn);
    if ( _nrm_dfsn>=0)
      ((FaceOffset_3*)_prop)->set_normal_diffusion( _nrm_dfsn);
    if ( _feature_layer>=0)
      ((FaceOffset_3*)_prop)->set_feature_layer( _feature_layer);
    if ( _wght_scheme>=0)
      ((FaceOffset_3*)_prop)->set_weighting_scheme( _wght_scheme);
    if ( _smoother>=0)
      ((FaceOffset_3*)_prop)->set_smoother( _smoother);
    if ( _conserv>=0)
      ((FaceOffset_3*)_prop)->set_conserve( _conserv);
  }

  COM::DataItem *spd = spds_io ? 
    _buf->inherit( spds_io, "rpr_spd_buf", 
		   COM::Pane::INHERIT_USE, true, NULL, 0) : NULL;
  _buf->init_done( false);

  // Subcycle to take small time steps
  double t = 0, t_rem = *dt;
  while ( t_rem > 0 || *dt==0) {
    // Time stepping
    int smoothed=(_rediter>0);
    double max_dt = _prop->time_stepping( spd, t_rem, du, &smoothed);

    if ( _verb && _rank==0)
      std::cout << "Rocprop: Subcycling with time step " << max_dt << std::endl;

    t = *dt - (t_rem-max_dt);
    t_rem = *dt - t;

    if ( t_rem > 0 && max_dt < *dt * _time_lb) {
      // If code is present, then set it to nonzero value 
      // if time step is large.
      if ( code) {
	*code = -1;
	std::cerr << "Rocprop: Given time step could not be reached" << std::endl;
	break; 
      }
      else if ( dt_elapsed==NULL) {
	std::cerr << "Rocprop: Time step is smaller than the lower bound "
		  << *dt*_time_lb << '=' << *dt << '*' << _time_lb 
		  << "\nMesh is too distorted. Stopping..." << std::endl;
	COM_assertion( max_dt >= *dt * _time_lb); abort();
      }
    }
    if ( code) *code = 0;

    // Add the dispments to the current coordinates
    Rocblas::add( nc, du, nc);

    // Reset the speed function
    if ( _cnstr_bound && spd && spd->is_elemental())
      _prop->bound_facial_speed( spd);

    // If propagation method is face offsetting, then perform mesh smoothing
    // by face offsetting with time step equal to 0.
    if ( _prop_method == PROP_FO) {
      for ( int i=smoothed; i<_rediter; ++i) {
	if ( _verb && _rank==0)
	  std::cout << "Rocprop: Perform additional mesh smoothing" << std::endl;
	_prop->time_stepping( NULL, 0, du, NULL);
	Rocblas::add( nc, du, nc);
      }
    }

    // If user request for elapsed time, then do not subcycle.
    if ( dt_elapsed) { *dt_elapsed = t; break; }
    else if ( *dt==0) break;
  }

  // Compute displacements and recover the original coordinates
  Rocblas::sub( nc, oldnc, du);
  Rocblas::copy( oldnc, nc);

  // Deallocate spds_buf
  if ( spd) _buf->delete_dataitem( spd->name());
  _buf->init_done( false);
}

void Rocprop::set_option( const char *opt, const char *val) 
{
  COM_assertion_msg( validate_object()==0, "Invalid object");

  std::string option;
  if ( opt) option = opt;

  if ( option == "method") {
    std::string value(val);
    int old_method = _prop_method;

    if ( value == "fo") // Face offsetting
      _prop_method = PROP_FO; 
    else
      _prop_method = PROP_MP; // Marker particle

    // Remove old object if method is changed.
    if ( _prop_method != old_method && _prop) 
    { delete _prop; _prop=NULL; }
  }
  else if ( option == "eigthres") {
    _eig_thres = std::atof( val);
  }
  else if ( option == "smoother") {
    std::string value(val);
    if ( value == "none")
      _smoother = SMOOTHER_NONE;
    else if ( value == "centroid")
      _smoother = SMOOTHER_ANGLE_WEIGHTED_CENTROID;
    else if ( value == "anisotropic")
      _smoother = SMOOTHER_ANISOTROPIC;
    else {
      COM_assertion( value == "Laplacian" || value == "laplacian");
      _smoother = SMOOTHER_LAPLACIAN;
    }
  }
  else if ( option == "rediter") {
    _rediter = std::atoi( val);
  }
  else if ( option == "conserv") {
    _conserv = std::atoi( val);
  }
  else if ( option == "courant") {
    _courant = std::atof( val);
  }
  else if ( option == "wavefrontal") {
    _wf_expn = std::atoi( val);
  }
  else if ( option == "weight") {
    _wght_scheme = std::atoi( val);
  }
  else if ( option == "normaldif") {
    _nrm_dfsn = std::atoi( val);
  }
  else if ( option == "feature_layer") {
    _feature_layer = std::atoi( val);
  }
  else if ( option == "sangle" || option == "fangle_strong") {
    _fangle_strong = std::atof( val);
  }
  else if ( option == "fangle" || option=="fangle_weak") {
    _fangle_weak = std::atof( val);
  }
  else if ( option == "tangle" || option == "fangle_tun") {
    _fangle_turn = std::atof( val);
  }
  else if ( option == "verbose")
    _verb = atoi( val);
  else
    COM_assertion_msg( false, (std::string("Unknown option ")+option).c_str());
}

void Rocprop::
remesh_serial( COM::DataItem *mesh_out, double *lave, double *fangle) {
  COM_assertion( !COMMPI_Initialized() || COMMPI_Comm_rank(MPI_COMM_WORLD)==0);
  
  if ( _rem) {
    double l = lave?*lave:0.;
    if ( l==0) compute_edge_lengths( &l, NULL, NULL);

    _rem->remesh_serial( _wm, mesh_out, l, fangle?*fangle:0.);
  }
  else {
    if ( !COMMPI_Initialized() || COMMPI_Comm_rank( MPI_COMM_WORLD)==0)
      std::cerr << "No remesher was registered" << std::endl;
  }
}

void Rocprop::load( const std::string &mname) {
  Rocprop *rp = new Rocprop();

  COM_new_window( mname.c_str());

  std::string glb=mname+".global";

  COM_new_dataitem( glb.c_str(), 'w', COM_VOID, 1, "");
  COM_set_object( glb.c_str(), 0, rp);

  COM_Type types[7];
  types[0] = types[2] = COM_RAWDATA; types[1] = COM_METADATA;
  COM_set_member_function( (mname+".initialize").c_str(), 
			   (Member_func_ptr)(&Rocprop::initialize), 
			   glb.c_str(), "biB", types);
    
  COM_set_member_function( (mname+".set_constraints").c_str(),
			   (Member_func_ptr)(&Rocprop::set_constraints), 
			   glb.c_str(), "bi", types);
  
  COM_set_member_function( (mname+".set_bounds").c_str(),
			   (Member_func_ptr)(&Rocprop::set_bounds),
			   glb.c_str(), "bi", types);
  
  types[2] = COM_DOUBLE; 
  COM_set_member_function( (mname+".perturb_mesh").c_str(),
			   (Member_func_ptr)(&Rocprop::perturb_mesh), 
			   glb.c_str(), "bbi", types);
    
  types[2] = types[4] = COM_METADATA; types[3] = types[5] = COM_DOUBLE; 
  types[6]=COM_INT; 
  COM_set_member_function( (mname+".propagate").c_str(),
			   (Member_func_ptr)(&Rocprop::propagate), 
			   glb.c_str(), "bibioOO", types);
    
  types[1] = types[2] = COM_STRING;
  COM_set_member_function( (mname+".set_option").c_str(),
			   (Member_func_ptr)(&Rocprop::set_option), 
			   glb.c_str(), "bii", types);

  types[1] = COM_METADATA; types[2] = types[3] = COM_DOUBLE; 
  COM_set_member_function( (mname+".remesh_serial").c_str(),
			   (Member_func_ptr)(&Rocprop::remesh_serial), 
			   glb.c_str(), "boII", types);

  types[1] = COM_VOID; types[2] = COM_INT;
  COM_set_member_function( (mname+".set_remesher").c_str(),
			   (Member_func_ptr)(&Rocprop::set_remesher), 
			   glb.c_str(), "biI", types);

  // Register inherited member functions from Rocsurf.
  std::string glb_surf=mname+".global_surf";
  Rocsurf *surf = rp;
  COM_new_dataitem( glb_surf.c_str(), 'w', COM_VOID, 1, "");
  COM_set_object( glb_surf.c_str(), 0, surf);
  
  types[1] = types[2] = types[3] = COM_DOUBLE;
  COM_set_member_function( (mname+".compute_edge_lengths").c_str(),
			   (Member_func_ptr)(&Rocsurf::compute_edge_lengths),
			   glb_surf.c_str(), "boOO", types);

  COM_window_init_done( mname.c_str());
}

void Rocprop::unload( const std::string &mname) {
  Rocprop *rp;
  std::string glb=mname+".global";

  COM_get_object( glb.c_str(), 0, &rp);
  
  COM_assertion_msg( rp->validate_object()==0, "Invalid object");

  delete rp;
  COM_delete_window( mname.c_str());
}

// C/C++ binding
extern "C" void Rocprop_load_module( const char *name)
{ Rocprop::load( std::string(name)); }
extern "C" void Rocprop_unload_module( const char *name) 
{ Rocprop::unload( std::string(name)); }

// Fortran binding
extern "C" void rocprop_load_module( const char *name, long int length)
{ Rocprop::load( std::string(name, length)); }
extern "C" void rocprop_unload_module( const char *name, long int length) 
{ Rocprop::unload( std::string(name, length)); }

extern "C" void ROCPROP_LOAD_MODULE( const char *name, long int length)
{ Rocprop::load( std::string(name, length)); }
extern "C" void ROCPROP_UNLOAD_MODULE( const char *name, long int length) 
{ Rocprop::unload( std::string(name, length)); }

extern "C" void rocprop_load_module_( const char *name, long int length)
{ Rocprop::load( std::string(name, length)); }
extern "C" void rocprop_unload_module_( const char *name, long int length) 
{ Rocprop::unload( std::string(name, length)); }

extern "C" void ROCPROP_LOAD_MODULE_( const char *name, long int length)
{ Rocprop::load( std::string(name, length)); }
extern "C" void ROCPROP_UNLOAD_MODULE_( const char *name, long int length) 
{ Rocprop::unload( std::string(name, length)); }

PROP_END_NAMESPACE






