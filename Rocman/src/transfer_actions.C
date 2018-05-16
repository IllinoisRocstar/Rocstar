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
// $Id: transfer_actions.C,v 1.59 2008/12/06 08:45:22 mtcampbe Exp $

#include "basic_actions.h"
#include "transfer_actions.h"
#include "FluidAgent.h"
#include "SolidAgent.h"
#include "BurnAgent.h"

COM_EXTERN_MODULE(SurfUtil);
COM_EXTERN_MODULE(SurfX);

static inline void load_rocsurf()
{
  if (COM_get_window_handle("SURF") <= 0) {
    COM_LOAD_MODULE_STATIC_DYNAMIC( SurfUtil, "SURF");
  }
}

extern void compute_overlay( FluidAgent *fagent, SolidAgent *sagent, double t);

void _load_rocface(FluidAgent *fagent, SolidAgent *sagent, const RocmanControl_parameters *param)
{
  int rank = fagent->get_comm_rank();

    // INIT_ROCFACE of rocman.f90
  if (COM_get_window_handle( "RFC") <=0 ) {
    if(rank == 0)
      MAN_DEBUG(3, ("Rocstar: load module RocFace.\n"));
    COM_LOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");

    int RFC_setv = COM_get_function_handle( "RFC.set_verbose");
    COM_call_function( RFC_setv, &param->rfc_verb);

#if 0
      // if remeshed
    MAN_DEBUG(3, ("Rocstar: remeshed: %d.\n", param->remeshed));
    if (param->remeshed) {
      double t = fagent->get_coupling()->get_control_param()->current_time;
      compute_overlay( fagent, sagent, t);
    }
#endif

    // load overlay
    MPI_Comm comm = fagent->get_communicator();

    int RFC_read = COM_get_function_handle( "RFC.read_overlay");
    int fluid_mesh = COM_get_dataitem_handle_const( fagent->fluidBufNG+".mesh");
    int solid_mesh = COM_get_dataitem_handle_const( sagent->solidBuf+".mesh");
    std::string rfc_dir = "Rocman/"+fagent->get_modinstance_name()+sagent->get_modinstance_name()+"/";
    //std::string rfc_dir = "Rocman/RocfloRocsolid/";
    std::string fluid_dir = rfc_dir+"ifluid";
    std::string solid_dir = rfc_dir+"isolid";

    if(rank == 0)
      MAN_DEBUG(2, ("Rocstar: read RocFace overlay mesh from %s.\n", rfc_dir.c_str()));

    //COM_call_function( RFC_read, &fluid_mesh, &solid_mesh, &comm, fluid_dir.c_str(), solid_dir.c_str(), "HDF");
    COM_call_function( RFC_read, &fluid_mesh, &solid_mesh, &comm, fluid_dir.c_str(), solid_dir.c_str(), "CGNS");

    int RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
    int RFC_interpolate = COM_get_function_handle("RFC.interpolate");
    COM_set_profiling_barrier( RFC_transfer, comm);
    COM_set_profiling_barrier( RFC_interpolate, comm);
  }
}

void InterMeshTransfer::load_rocface(const RocmanControl_parameters *param)
{
    // INIT_ROCFACE of rocman.f90
  _load_rocface(fagent, sagent, param);
}

InterMeshTransfer::InterMeshTransfer( FluidAgent *fag, SolidAgent *sag, char *name):
                     Action(0, (const char **)NULL, NULL, NULL, name), 
		     fagent(fag), sagent(sag)
{
}


// Transfer load from fluid to solid
// Arguments: tf_FF (IN), mdot (IN), rb (IN), 
// tf_SF (OUT)
LoadTransfer_FS::LoadTransfer_FS( FluidAgent *fag, SolidAgent *sag, 
                                  const std::string f_ts, 
                                  const std::string s_ts, 
                                  const std::string s_pf) : 
  InterMeshTransfer(fag, sag, (char *)"LoadTransfer_FS"),
    f_ts_str(f_ts), s_ts_str(s_ts), s_pf_str(s_pf)
{
  int io[] = {IN, IN, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = f_ts;
  atts[1] = s_ts;
  atts[2] = s_pf;
  set_attr(3, atts);

  traction_mode = fagent->get_coupling()->get_rocmancontrol_param()->traction_mode;

  // create_dataitem
  // for SolidAgent ??????????????????
//  if (traction_mode == NO_SHEER && size_ts == 3) {
    sagent->register_new_dataitem( sagent->solidBufBase, ".pf", 'e', COM_DOUBLE, 1, "Pa");
//  }

  // for FluidAgent
  std::string::size_type pos = f_ts_str.find( ".");
  COM_assertion_msg(pos!=std::string::npos, "LoadTransfer_FS failed!");
  std::string f = f_ts_str.substr( 0, pos);
  std::string ts = f_ts_str.substr( pos, f_ts_str.size());
  if (traction_mode == NO_SHEER)  {
    fagent->register_clone_dataitem( 0, f, ts, fagent->get_surface_window(), ".pf");
  }
  else
    fagent->register_clone_dataitem( 0, f, ts, fagent->get_surface_window(), ".tf");
}

// Create buffer data
void LoadTransfer_FS::init( double t) {
  // find out size_ts ???
  int dummy;
  std::string unit;
  char loc;
  COM_get_dataitem( sagent->solidBufBase+".ts", &loc, &dummy, &size_ts, &unit);
  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions must be vectors!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());

  // for SolidAgent
  f_ts_hdl = get_dataitem_handle(0);
  s_ts_hdl = get_dataitem_handle(1);
  if (traction_mode == NO_SHEER && size_ts == 3)
    s_pf_hdl = get_dataitem_handle(2);
  else 
    s_pf_hdl = -1;

  f_tf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufNG+".tf");
  f_pf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufNG+".pf");

  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");

  load_rocsurf();
  SURF_compute_face_normals = COM_get_function_handle( "SURF.compute_element_normals");

  MAN_DEBUG(3, ("LoadTransfer_FS::init() called - traction_mode: %d size_ts: %d.\n", traction_mode, size_ts));
}

void LoadTransfer_FS::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: LoadTransfer_FS::run() with t:%e dt:%e.\n", t, dt));
  // ts = tf + ( mdot*Vs

  // part 1  (POST_UPDATE_FLUID in fluid_agent.f90)
  if (traction_mode != NO_SHEER) {
    COM_call_function( RocBlas::copy, &f_tf_hdl, &f_ts_hdl);
  }
  else {   // ts is a scalar
    COM_call_function( RocBlas::copy, &f_pf_hdl, &f_ts_hdl);
  }

  // part 2  (solid_agent.f90 INIT_INBUFF_SOLID())
  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_pf_hdl);
    COM_call_function( SURF_compute_face_normals, &s_ts_hdl);
    COM_call_function( RocBlas::mul, &s_ts_hdl, &s_pf_hdl, &s_ts_hdl);
    COM_call_function( RocBlas::neg, &s_ts_hdl, &s_ts_hdl);
  }
  else if ( traction_mode == NO_SHEER && size_ts == 1) {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  }
  else if (size_ts == 1) {
    COM_assertion_msg(0, "ERROR: NOT IMPLEMENTED!");
  }
  else {
    COM_assertion_msg(0, "ERROR: NOT IMPLEMENTED!");
  }
}

// Transfer load from fluid to solid with Burn
// Arguments: tf_FF (IN), mdot (IN), rb (IN), 
// tf_SF (OUT)
LoadTransfer_FSc_ALE::LoadTransfer_FSc_ALE( FluidAgent *fag, SolidAgent *sag, 
                                  BurnAgent *bag,
                                  const std::string f_pf,
                                  const std::string fb_mdot, 
                                  const std::string b_rb, 
                                  const std::string s_ts, 
                                  const std::string s_pf) : 
  InterMeshTransfer(fag, sag, (char *)"LoadTransfer_FSc_ALE"),
    bagent(bag) 
{
  int io[] = {IN, IN, IN, IN, OUT};
  set_io( 5, io); 
  
  std::string atts[5];
  atts[0] = f_pf;
  atts[1] = fb_mdot;
  atts[2] = b_rb;
  atts[3] = s_ts;
  atts[4] = s_pf;
  set_attr(5, atts);

  traction_mode = fagent->get_coupling()->get_rocmancontrol_param()->traction_mode;

  // create_dataitem
  // for SolidAgent ??????????????????
//  if (traction_mode == NO_SHEER && size_ts == 3) {
    sagent->register_new_dataitem( sagent->solidBufBase, ".pf", 'e', COM_DOUBLE, 1, "Pa");
//  }

  // for FluidAgent
  if (traction_mode == NO_SHEER)  {
    fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".ts", fagent->get_surface_window(), ".pf");
  }
  else
    fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".ts", fagent->get_surface_window(), ".tf");

//  fagent->register_new_dataitem( fagent->fluidBufB, ".mdot_tmp", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
}

// Create buffer data
void LoadTransfer_FSc_ALE::init( double t) {
  // find out size_ts 
  int dummy;
  std::string unit;
  char loc;
  COM_get_dataitem( sagent->solidBufBase+".ts", &loc, &dummy, &size_ts, &unit);
  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions must be vectors!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  //f_pf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufNG+".pf");
  f_pf_hdl = get_dataitem_handle_const( 0);
  fb_mdot_hdl = get_dataitem_handle( 1);
  b_rb_hdl = get_dataitem_handle( 2);
  s_ts_hdl = get_dataitem_handle( 3);
  if (traction_mode == NO_SHEER && size_ts == 3)
    s_pf_hdl = get_dataitem_handle( 4);
  else 
    s_pf_hdl = -1;

  f_ts_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".ts");
  f_tf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufNG+".tf");
  fb_ts_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".ts");
  fb_pf_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".pf");

  fb_rhof_alp_hdl = COM_get_dataitem_handle ( fagent->fluidBufB+".rhof_alp");
  fb_mdot_tmp_hdl = COM_get_dataitem_handle ( fagent->fluidBufB+".mdot_tmp");

  fb_nf_alp_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".nf_alp");
  fb_tf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufB+".tf");

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");

  load_rocsurf();
  SURF_compute_face_normals = COM_get_function_handle( "SURF.compute_element_normals");

  MAN_DEBUG(3, ("LoadTransfer_FSc_ALE::init() called - traction_mode: %d size_ts: %d.\n", traction_mode, size_ts));
}

void LoadTransfer_FSc_ALE::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: LoadTransfer_FSc_ALE::run() with t:%e dt:%e.\n", t, dt));
  // ts = tf + ( mdot*Vs

#if 1
  // part 1  (POST_UPDATE_FLUID in fluid_agent.f90)
  if (traction_mode != NO_SHEER) {
    COM_call_function( RocBlas::copy, &f_tf_hdl, &f_ts_hdl);
                                                                                
      // Compute tfmts = t_s-t_f (see developers guide).
      // Here fb_mdot_tmp is used as buffer space.
    COM_call_function( RocBlas::div, &fb_mdot_hdl, &fb_rhof_alp_hdl, &fb_mdot_tmp_hdl);
                                                                                
    if ( sagent->withALE) {
      COM_call_function( RocBlas::sub, &fb_mdot_tmp_hdl, &b_rb_hdl, &fb_mdot_tmp_hdl);
    }      
                                                                                
    COM_call_function( RocBlas::mul, &fb_mdot_tmp_hdl, &fb_mdot_hdl, &fb_mdot_tmp_hdl);
    COM_call_function( RocBlas::mul, &fb_nf_alp_hdl, &fb_mdot_tmp_hdl, &fb_ts_hdl);
                                                                                
       // Compute ts = tf - (tf-ts)
    COM_call_function( RocBlas::sub, &fb_tf_hdl, &fb_ts_hdl, &fb_ts_hdl);
  }
  else {   // ts is a scalar

    COM_call_function( RocBlas::copy, &f_pf_hdl, &f_ts_hdl);

    //    Compute tsmpf = p_f-t_s (see developers guide).
    //    Here fb_mdot_tmp is used as buffer space.
    COM_call_function( RocBlas::div, &fb_mdot_hdl, &fb_rhof_alp_hdl, &fb_mdot_tmp_hdl);

    if ( sagent->withALE) {
        COM_call_function( RocBlas::sub, &fb_mdot_tmp_hdl, &b_rb_hdl, &fb_mdot_tmp_hdl);
    }

    COM_call_function( RocBlas::mul, &fb_mdot_tmp_hdl, &fb_mdot_hdl, &fb_ts_hdl);

     //    Compute ts = pf - (pf-ts)
    COM_call_function( RocBlas::sub, &fb_pf_hdl, &fb_ts_hdl, &fb_ts_hdl);

    //debug_print(fagent->fluidBufB+".mdot", 102, 0, "LOADTRANSFER");
    //debug_print(fagent->fluidBufB+".ts", 102, 0, "LOADTRANSFER");

     // Subtract from P_ambient if not zeros 
    double P_ambient = fagent->get_coupling()->get_rocmancontrol_param()->P_ambient;
    if ( P_ambient != 0.0) {
       COM_call_function( RocBlas::sub_scalar, &f_ts_hdl, &P_ambient, &f_ts_hdl);
    }
  }

#endif

  // part 2  (solid_agent.f90 INIT_INBUFF_SOLID())
  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_pf_hdl);
    COM_call_function( SURF_compute_face_normals, &s_ts_hdl);
    COM_call_function( RocBlas::mul, &s_ts_hdl, &s_pf_hdl, &s_ts_hdl);
    COM_call_function( RocBlas::neg, &s_ts_hdl, &s_ts_hdl);
  }
  else if ( traction_mode == NO_SHEER && size_ts == 1) {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_ts_hdl);
    //debug_print(sagent->solidBuf+".ts", 102, 0, fagent->get_communicator(), "LOADTRANSFER");
  }
  else if (size_ts == 1) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions must be vectors!");
  }
  else {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  }
}

// Transfer load from fluid to solid with Burn  (part 2 only)
// Arguments: f_ts (IN)
// s_ts (OUT) tf_SF (OUT)
LoadTransferOnly_FSc_ALE::LoadTransferOnly_FSc_ALE( 
                                  FluidAgent *fag, SolidAgent *sag, 
                                  BurnAgent *bag,
                                  const std::string f_ts,
                                  const std::string s_ts, 
                                  const std::string s_pf) : 
  InterMeshTransfer(fag, sag, (char *)"LoadTransferOnly_FSc_ALE"),
    bagent(bag) 
{
  int io[] = {IN, OUT, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = f_ts;
  atts[1] = s_ts;
  atts[2] = s_pf;
  set_attr(3, atts);

  traction_mode = fagent->get_coupling()->get_rocmancontrol_param()->traction_mode;

  // create_dataitem
  // for SolidAgent ??????????????????
//  if (traction_mode == NO_SHEER && size_ts == 3) {
    sagent->register_new_dataitem( sagent->solidBufBase, ".pf", 'e', COM_DOUBLE, 1, "Pa");
//  }

  // for FluidAgent
  if (traction_mode == NO_SHEER)  {
    fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".ts", fagent->get_surface_window(), ".pf");
  }
  else
    fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".ts", fagent->get_surface_window(), ".tf");

//  fagent->register_new_dataitem( fagent->fluidBufB, ".mdot_tmp", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
}

// Create buffer data
void LoadTransferOnly_FSc_ALE::init( double t) {
  // find out size_ts 
  int dummy;
  std::string unit;
  char loc;
  COM_get_dataitem( sagent->solidBufBase+".ts", &loc, &dummy, &size_ts, &unit);
  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions must be vectors!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  //f_pf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufNG+".pf");
  f_ts_hdl = get_dataitem_handle_const( 0);
  s_ts_hdl = get_dataitem_handle( 1);
  if (traction_mode == NO_SHEER && size_ts == 3)
    s_pf_hdl = get_dataitem_handle( 2);
  else 
    s_pf_hdl = -1;

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");

  load_rocsurf();
  SURF_compute_face_normals = COM_get_function_handle( "SURF.compute_element_normals");

  MAN_DEBUG(3, ("LoadTransferOnly_FSc_ALE::init() called - traction_mode: %d size_ts: %d.\n", traction_mode, size_ts));
}

void LoadTransferOnly_FSc_ALE::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: LoadTransferOnly_FSc_ALE::run() with t:%e dt:%e.\n", t, dt));

  // part 2  (solid_agent.f90 INIT_INBUFF_SOLID())
  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_pf_hdl);
    COM_call_function( SURF_compute_face_normals, &s_ts_hdl);
    COM_call_function( RocBlas::mul, &s_ts_hdl, &s_pf_hdl, &s_ts_hdl);
    COM_call_function( RocBlas::neg, &s_ts_hdl, &s_ts_hdl);
  }
  else if ( traction_mode == NO_SHEER && size_ts == 1) {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  }
  else if (size_ts == 1) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions must be vectors!");
  }
  else {
    COM_call_function( RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  }
}

/*
// Arguments: tf_FF (IN), mdot (IN), rb (IN), 
// tf_SF (OUT)
LoadTransfer_FSc_ALE::LoadTransfer_FSc_ALE( FluidAgent *fag, SolidAgent *sag, const char *at[], int *is, void *p) : 
  InterMeshTransfer(fag, sag, "LoadTransfer_FS") {
  int io[] = {IN, OUT};
  set_io( 2, io); 
}

// Create buffer data
void LoadTransfer_FSc_ALE::init( double t) {
}

void LoadTransfer_FSc_ALE::run( double t, double dt, double alpha) {
  // ts = tf + ( mdot*Vs
}
*/


// solid_agent.f90 POST_UPDATE_SOLID
GetDeformedMesh::GetDeformedMesh(FluidAgent *fag, SolidAgent *sag, const std::string s_x, const std::string s_uhat, const std::string s_y):
  InterMeshTransfer(fag, sag, (char *)"GetDeformedMesh"),
    s_x_str(s_x), s_uhat_str(s_uhat), s_y_str(s_y)
{
  int io[] = {IN, IN, OUT};
  set_io( 3, io); 
  
  std::string atts[4];
  atts[0] = s_x;
  atts[1] = s_uhat;
  atts[2] = s_y;
  set_attr(3, atts);

    // register dataitems
/*
  std::string::size_type pos = s_x_str.find( ".");
  COM_assertion_msg(pos!=std::string::npos, "GetDeformedMesh::create_dataitem failed!");
  std::string s = s_x_str.substr( 0, pos);
  std::string x = s_x_str.substr( pos, s_x_str.size());
  sagent->register_use_dataitem( s, x, sagent->solidBufBase, ".nc");

  pos = s_y_str.find( ".");
  COM_assertion_msg(pos!=std::string::npos, "GetDeformedMesh::create_dataitem failed!");
  s = s_y_str.substr( 0, pos);
  std::string y = s_y_str.substr( pos, s_y_str.size());
  sagent->register_clone_dataitem( 0, s, y, sagent->get_surface_window(), ".nc");
*/
}

void GetDeformedMesh::init( double t) {
  s_x_hdl = get_dataitem_handle(0);
  s_uhat_hdl = get_dataitem_handle(1);
  s_y_hdl = get_dataitem_handle(2);

#if 0
  // POST_INIT_SOLID
  COM_call_function( RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
    // initial_start  ?????
  if (t == 0.0) {
//    double alpha = 0.0;
//     ????????????????? strange bug called in bcinitscheduler::init()
//    sagent->obtain_bc(&alpha);
  }
#endif
}

// POST_UPDATE_SOLID
void GetDeformedMesh::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: calling GetDeformedMesh::run() with t:%e dt:%e alpha:%e.\n", t, dt, alpha));
  COM_call_function( RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
}

// solid_agent.f90 POST_UPDATE_SOLID
GetDeformedMesh_ALE::GetDeformedMesh_ALE(FluidAgent *fag, SolidAgent *sag, const std::string s_x, const std::string s_uhat, const std::string s_y, double z):
  InterMeshTransfer(fag, sag, (char *)"GetDeformedMesh_ALE"),
    s_x_str(s_x), s_uhat_str(s_uhat), s_y_str(s_y), zoom(z)
{
  int io[] = {IN, IN, OUT};
  set_io( 3, io); 
  
  std::string atts[4];
  atts[0] = s_x;
  atts[1] = s_uhat;
  atts[2] = s_y;
  set_attr(3, atts);

/*
    // register dataitems
  std::string::size_type pos = s_x_str.find( ".");
  COM_assertion_msg(pos!=std::string::npos, "GetDeformedMesh_ALE::create_dataitem failed!");
  std::string s = s_x_str.substr( 0, pos);
  std::string x = s_x_str.substr( pos, s_x_str.size());
  sagent->register_use_dataitem( s, x, sagent->solidBufBase, ".nc");

    // s_y = solidBuf+".nc"
  pos = s_y_str.find( ".");
  COM_assertion_msg(pos!=std::string::npos, "GetDeformedMesh_ALE::create_dataitem failed!");
  s = s_y_str.substr( 0, pos);
  std::string y = s_y_str.substr( pos, s_y_str.size());
  sagent->register_clone_dataitem( 0, s, y, sagent->get_surface_window(), y);
*/
}

void GetDeformedMesh_ALE::init( double t) {
  s_x_hdl = get_dataitem_handle(0);
  s_uhat_hdl = get_dataitem_handle(1);
  s_y_hdl = get_dataitem_handle(2);

  withALE = sagent->withALE;

  s_rb_hdl = COM_get_dataitem_handle( sagent->solidBuf+".rb");
  s_mdot_hdl = COM_get_dataitem_handle( sagent->solidBuf+".mdot");
  s_areas_hdl = COM_get_dataitem_handle( sagent->solidBuf+".areas");
  s_rhos_hdl = COM_get_dataitem_handle( sagent->solidBuf+".rhos");

  load_rocsurf();
  SURF_compute_bounded_volumes = COM_get_function_handle( "SURF.compute_bounded_volumes");
  SURF_compute_face_areas = COM_get_function_handle( "SURF.compute_element_areas");

#if 0
  // POST_INIT_SOLID
  COM_call_function( RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
    // initial_start  ?????
  if (t == 0.0) {
//    double alpha = 0.0;
//     ????????????????? strange bug called in bcinitscheduler::init()
//    sagent->obtain_bc(&alpha);
  }
#endif
}

// POST_UPDATE_SOLID
void GetDeformedMesh_ALE::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: calling GetDeformedMesh_ALE::run() with t:%e dt:%e alpha:%e withALE=%d zoom=%e.\n", t, dt, alpha, withALE, zoom));

  double mdot_min;

  if ( withALE && zoom>0 && !fagent->get_coupling()->get_rocmancontrol_param()->PROP_fom) {
      // compute volumes only for faces with nonzero burning rate
    COM_call_function( RocBlas::copy, &s_rb_hdl, &s_mdot_hdl);

    //debug_print(sagent->solidBuf+".nc", 202, 1);
    //debug_print(sagent->solidBuf+".x", 202, 1);

      // Compute mass injection using s_x and s_y as buffers.
      //  So far s_y stores the undeformed configuration s_x.
    int one = 1;
    COM_call_function( SURF_compute_bounded_volumes, &s_x_hdl, &s_y_hdl, &s_mdot_hdl, &one);

    //debug_print(sagent->solidBuf+".nc", 202, 1);
    //debug_print(attr[2], 202, 1);

#if 0
      //  Check whether the burning rates are nonnegative.
    int comm = MPI_COMM_SELF;
    COM_call_function( RocBlas::min_scalar_MPI, &s_mdot_hdl, &mdot_min, &comm);

    if ( mdot_min < 0) {
          printf("Rocstar ERROR: Negative mdot found %e. Aborting...\n", mdot_min);
          MPI_Abort( MPI_COMM_WORLD, -1);
    }
#else
    double zero = 0.0;
    COM_call_function( RocBlas::maxof_scalar, &s_mdot_hdl, &zero, &s_mdot_hdl);
#endif

    // Compute deformed configuration
    COM_call_function( RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
    COM_call_function( SURF_compute_face_areas, &s_areas_hdl);

    // Compute mdot=dV*rhos/area/dt
    COM_call_function( RocBlas::mul, &s_mdot_hdl, &s_rhos_hdl, &s_mdot_hdl);
    COM_call_function( RocBlas::div, &s_mdot_hdl, &s_areas_hdl, &s_mdot_hdl);
    COM_call_function( RocBlas::div_scalar, &s_mdot_hdl, &dt, &s_mdot_hdl);
  }
  else {
      // Compute deformed configuration
    COM_call_function( RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
      // Compute mdot as rb*rhos
    if ( withALE)
        COM_call_function( RocBlas::mul, &s_rb_hdl, &s_rhos_hdl, &s_mdot_hdl);
  }
}


// fluid_agent.f90  INIT_INBUFF_FLUID
MeshMotionTransfer_SF::MeshMotionTransfer_SF(FluidAgent *fag, SolidAgent *sag, 
    const std::string s_u, const std::string f_total_disp, 
    const std::string f_vm):
  InterMeshTransfer(fag, sag, (char *)"MeshMotionTransfer_SF")
{
  // isolid_i+".nc" (in), isolid_i+".u" (in), fagent->ifluid_i+".vm" (out)
  int io[] = {IN, IN, INOUT, OUT};
  set_io( 4, io); 
  
  std::string atts[4];
  atts[0] = sagent->solidBufBase+".nc";
  atts[1] = s_u;
  atts[2] = f_total_disp;
  atts[3] = f_vm;
  set_attr(4, atts);

    // Change to transfer total_displacement from solid to fluids, and obtain
    // incremental displacement as (nc_t0+total_disp)-nc_tn. This is prone
    // to cancellation errors but avoids accumulation of any errors.
  fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".total_disp", fagent->get_surface_window(), ".du_alp");
  fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".nc_t0", fagent->get_surface_window(), ".nc");
    // nc_tmp is used to compare against solid nodal coordinates.
  //fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".nc_tmp", fagent->get_surface_window(), ".nc");
  //fagent->register_new_dataitem( fagent->ifluid_i, ".sq_dist", 'n', COM_DOUBLE, 1, "m");
}

void MeshMotionTransfer_SF::init( double t) {
  // skip 0
  s_u_hdl = get_dataitem_handle(1);
  f_total_disp_hdl = get_dataitem_handle(2);
  f_vm_hdl = get_dataitem_handle(3);

  f_nc_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".nc");
  f_nc_t0_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".nc_t0");

  // POST_INIT_FLUID in fluid_agent.f90
    // initial_start  ?????
//  if (t == 0.0) {
    // TODO ???
//    COM_call_function( RocBlas::copy, &f_nc_hdl, &f_nc_t0_hdl);
      // UPDATE_INBUFF_GM_FLUID     ???????????????
//    double zero = 0.0;
//    fagent->obtain_gm(&zero);
//  }

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
}

// INIT_INBUFF_FLUID
void MeshMotionTransfer_SF::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: calling MeshMotionTransfer_SF::run() with t:%e dt:%e alpha:%e.\n", t, dt, alpha));

    // b. Interpolates total displacement from solid nodes to fluid nodes
  //debug_print(sagent->solidBuf+".u", 109, 0, sagent->get_communicator(), "U");
  COM_call_function( RFC_interpolate, &s_u_hdl, &f_total_disp_hdl);

    // Compute incremental mesh volocity
  if (!fagent->get_coupling()->initial_start()) {
    // Compute incremental mesh volocity  
    COM_call_function( RocBlas::add, &f_nc_t0_hdl, &f_total_disp_hdl, &f_vm_hdl);
    COM_call_function( RocBlas::sub, &f_vm_hdl, &f_nc_hdl, &f_vm_hdl);
    COM_call_function( RocBlas::div_scalar, &f_vm_hdl, &dt, &f_vm_hdl);
  }
}


// transfer vs from solid nodes to fluid faces using Rocface
DeformationVelTransfer_SF::DeformationVelTransfer_SF(FluidAgent *fag, 
             SolidAgent *sag, const std::string s_vs, const std::string f_vs):
  InterMeshTransfer(fag, sag, (char *)"DeformationVelTransfer_SF"),
             s_vs_str(s_vs), f_vs_str(f_vs)
{
  // added a builtin in (ifluid_i+".vm")
  // ifluid_i+".vm" (in),  isolid_i+".vs", ifluid_i+".vs"
  int io[] = {IN, IN, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = fagent->ifluid_i+".vm";
  atts[1] = s_vs;
  atts[2] = f_vs;
  set_attr(3, atts);
}

void DeformationVelTransfer_SF::init( double t) {
  s_vs_hdl = get_dataitem_handle(1);
  f_vs_hdl = get_dataitem_handle(2);

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());

  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void DeformationVelTransfer_SF::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: DeformationVelTransfer_SF::run() with t:%e dt:%e alpha:%e.\n", t, dt, alpha));

    // c. Transfer vs from solid nodes to fluid faces using Rocface
  COM_call_function( RFC_transfer, &s_vs_hdl, &f_vs_hdl);
}


MeshMotionTransferISS::MeshMotionTransferISS(FluidAgent *fag, SolidAgent *sag, const std::string s_u, const std::string s_vs, const std::string f_vm):
  InterMeshTransfer(fag, sag, (char *)"MeshMotionTransferISS"),
    s_u_str(s_u), s_vs_str(s_vs), f_vm_str(f_vm)
{
  // isolid_i+".u" (in), isolid_i+".vs" (in), ifluid_i+".vm" (out)
  int io[] = {IN, IN, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = s_u;
  atts[1] = s_vs;
  atts[2] = f_vm;
  set_attr(3, atts);

  std::string::size_type pos = f_vm_str.find( ".");
  COM_assertion_msg(pos!=std::string::npos, "MeshMotionTransferISS failed!");
  std::string f = f_vm_str.substr( 0, pos);
  std::string vm = f_vm_str.substr( pos, f_vm_str.size());
  fagent->register_new_dataitem( f, vm, 'n', COM_DOUBLE, 3, "m/s");
    // used internally
  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".u", fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".uold", fagent->fluidBufNG, ".vm");
   // avoid vs because vs is 'e'
  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".fvs", fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".vsold", fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".utmp", fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".vtmp", fagent->fluidBufNG, ".vm");
}

void MeshMotionTransferISS::init( double t) {
  s_u_hdl = get_dataitem_handle(0);
  s_vs_hdl = get_dataitem_handle(1);
  f_vm_hdl = get_dataitem_handle(2);

  f_u_hdl = COM_get_dataitem_handle(fagent->fluidBufNG+".u");
  f_vs_hdl = COM_get_dataitem_handle(fagent->fluidBufNG+".fvs");
  f_uold_hdl = COM_get_dataitem_handle(fagent->fluidBufNG+".uold");
  f_vsold_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".vsold");
  f_utmp_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".utmp");
  f_vtmp_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".vtmp");

    // set zero
  double zero = 0.0;
  COM_call_function( RocBlas::copy_scalar, &zero, &f_uold_hdl);
  COM_call_function( RocBlas::copy_scalar, &zero, &f_vsold_hdl);

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
}

// INIT_INBUFF_FLUID
void MeshMotionTransferISS::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("Rocstar: calling MeshMotionTransferISS::run() with t:%e dt:%e.\n", t, dt));
  // fomular:  Vm = (un-un-1)/deltaT + (vsn - vsn-1) /2

  COM_call_function( RFC_interpolate, &s_u_hdl, &f_u_hdl);
  COM_call_function( RFC_interpolate, &s_vs_hdl, &f_vs_hdl);

  COM_call_function( RocBlas::sub, &f_u_hdl, &f_uold_hdl, &f_utmp_hdl);
  
  COM_call_function( RocBlas::div_scalar, &f_utmp_hdl, &dt, &f_utmp_hdl);

  COM_call_function( RocBlas::sub, &f_vs_hdl, &f_vsold_hdl, &f_vtmp_hdl);

  double half = 0.5;
  COM_call_function( RocBlas::mul_scalar, &f_vtmp_hdl, &half, &f_vtmp_hdl);

  COM_call_function( RocBlas::add, &f_utmp_hdl, &f_vtmp_hdl, &f_vm_hdl);

  // save copies
  COM_call_function( RocBlas::copy, &f_u_hdl, &f_uold_hdl);
  COM_call_function( RocBlas::copy, &f_vs_hdl, &f_vsold_hdl);
}

// burn_agent.f90  INIT_INBUFF_BURN
// Transfer rhos from solid faces to fluid faces
TransferSolidDensity::TransferSolidDensity(FluidAgent *fag, SolidAgent *sag,
                           const std::string s_rhos, const std::string f_rhos):
  InterMeshTransfer(fag, sag, (char *)"TransferSolidDensity"),
                s_rhos_str(s_rhos), f_rhos_str(f_rhos)
{
  int io[] = {IN, OUT};
  set_io( 2, io); 
  
  std::string atts[2];
  atts[0] = s_rhos;
  atts[1] = f_rhos;
  set_attr(2, atts);

  fagent->register_new_dataitem(fagent->fluidBufNG, ".rhos", 'e', COM_DOUBLE, 1, "kg/(m^3)");

/*
  sagent->register_use_dataitem( sagent->solidBuf, ".x", sagent->solidBufBase, ".nc");
  sagent->register_clone_dataitem( 0, sagent->solidBuf, ".nc", sagent->get_surface_window(), ".nc");
*/
}

void TransferSolidDensity::init( double t) {
  MAN_DEBUG(3, ("Rocstar: TransferSolidDensity::init() with t:%e.\n", t));
  s_rhos_hdl = get_dataitem_handle(0);
  f_rhos_hdl = get_dataitem_handle(1);

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());

  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void TransferSolidDensity::run( double t, double dt, double alpha) {
  int rhos_mode = sagent->rhos_mode;

  MAN_DEBUG(3, ("Rocstar: TransferSolidDensity::run() with rhs_mode=%d t:%e dt:%e alpha:%e.\n", rhos_mode, t, dt, alpha));

  if (rhos_mode == 1)   // Constant density
    COM_call_function( RocBlas::copy, &s_rhos_hdl, &f_rhos_hdl);
  else if (fagent->get_coupling()->initial_start())
    COM_call_function( RFC_transfer, &s_rhos_hdl, &f_rhos_hdl);
  else if (rhos_mode == 3)    // Varying density
    COM_call_function( RFC_transfer, &s_rhos_hdl, &f_rhos_hdl);
}


// solid_agent.f90  INIT_INBUFF_SOLID
// Transfer rb (burn rate) from fluid faces to solid faces
TransferBurnRate_FS_ALE::TransferBurnRate_FS_ALE(FluidAgent *fag, 
                SolidAgent *sag,
                const std::string b_rb, const std::string s_rb):
  InterMeshTransfer(fag, sag, (char *)"TransferBurnRate_FS_ALE"),
                b_rb_str(b_rb), s_rb_str(s_rb)
{
  int io[] = {IN, OUT};
  set_io( 2, io); 
  
  std::string atts[2];
  atts[0] = b_rb;
  atts[1] = s_rb;
  set_attr(2, atts);
}

void TransferBurnRate_FS_ALE::init( double t) {
  if (sagent->withALE == 0) return;

  b_rb_hdl = get_dataitem_handle(0);
  s_rb_hdl = get_dataitem_handle(1);
  int PROP_fom = fagent->get_coupling()->get_rocmancontrol_param()->PROP_fom;
  if ( PROP_fom) {
     p_rb_hdl = COM_get_dataitem_handle( sagent->propBufAll+".rb");
  }
  else {
     p_rb_hdl = COM_get_dataitem_handle( sagent->propBuf+".rb");
  }
  f_mdot_tmp_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".mdot_tmp");
  fb_mdot_tmp_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".mdot_tmp");

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void TransferBurnRate_FS_ALE::run( double t, double dt, double alpha) {
  if (sagent->withALE == 0) return;

  MAN_DEBUG(3, ("Rocstar: TransferBurnRate_FS_ALE::run() with t:%e dt:%e alpha:%e.\n", t, dt, alpha));

    // Transfer rb from fluid faces to solid faces
  double zero = 0.0;
  COM_call_function( RocBlas::copy_scalar, &zero, &p_rb_hdl);

    // Here we use f_mdot_tmp as buffer space for holding f_rb
  COM_call_function( RocBlas::copy_scalar, &zero, &f_mdot_tmp_hdl);

  COM_call_function( RocBlas::copy, &b_rb_hdl, &fb_mdot_tmp_hdl);

  COM_call_function( RFC_transfer, &f_mdot_tmp_hdl, &s_rb_hdl);
}

// INIT_INBUFF_FLUID
// Depending on whether solid has ALE or not, we need to compute mdot
// differently. If solid has ALE, to ensure mass conservation, we compute
// mdot using a geometric construction on the solid side and transfer
// mdot to fluids. If solid has no ALE, mass conservation is not possible
// and the above approach will generate zero mdot, so we have to compute
// mdot as rhos*rb in this case.
MassTransfer_SF_ALE::MassTransfer_SF_ALE(FluidAgent *fag, SolidAgent *sag, 
                 BurnAgent *bag, const std::string f_mdot): 
  InterMeshTransfer(fag, sag, (char *)"MassTransfer_SF_ALE"),
                 bagent(bag)
{
  int io[] = {OUT};
  set_io( 1, io); 
  
  std::string atts[1];
  atts[0] = f_mdot;
  set_attr(1, atts);

  fagent->register_new_dataitem(fagent->ifluid_i, ".rhos", 'e', COM_DOUBLE, 1, "kg/(m^3)");
  fagent->register_new_dataitem(fagent->ifluid_i, ".mdot", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  bagent->register_use_dataitem( bagent->iburn_all, ".rhos", bagent->parentWin, ".rhos");
}

void MassTransfer_SF_ALE::init( double t) {
  //f_mdot_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".mdot");
  f_mdot_hdl = get_dataitem_handle( 0);

  b_rhos_hdl = COM_get_dataitem_handle( bagent->iburn_ng+".rhos");
  b_rb_hdl = COM_get_dataitem_handle( bagent->iburn_ng+".rb");
  fb_mdot_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".mdot");

    // ALE
  f_total_disp_hdl = COM_get_dataitem_handle(fagent->fluidBufNG+".total_disp");
  f_nc_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".nc");
  f_nc_t0_hdl = COM_get_dataitem_handle( fagent->fluidBufNG+".nc_t0");
  s_mdot_hdl = COM_get_dataitem_handle( sagent->solidBuf+".mdot");

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void MassTransfer_SF_ALE::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("[%d] Rocstar: MassTransfer_SF_ALE::run() with t:%e dt:%e.\n", fagent->get_comm_rank(), t, dt));

  if (sagent->withALE) {
    COM_call_function( RocBlas::add,  &f_nc_t0_hdl, &f_total_disp_hdl, &f_nc_hdl);
    COM_call_function( RFC_transfer, &s_mdot_hdl, &f_mdot_hdl);
  }
  else {
    double zero = 0.0;
    COM_call_function( RocBlas::copy_scalar, &zero, &f_mdot_hdl);
    COM_call_function( RocBlas::mul, &b_rhos_hdl, &b_rb_hdl, &fb_mdot_hdl);
  }
}


// trasnfer temperature from solid to fluid
// s_Ts in solid is 'n'  f_Tf in fluid is 'e' 
TemperatureTransfer_SF::TemperatureTransfer_SF(SolidAgent *sag, FluidAgent *fag,
                 const std::string s_Ts, 
		 const std::string fb_Tflm, const std::string fn_Tb):
  InterMeshTransfer(fag, sag, (char *)"TemperatureTransfer_SF")
{
  int io[] = {IN, OUT, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = s_Ts;
  atts[1] = fb_Tflm;
  atts[2] = fn_Tb;
  set_attr(3, atts);

  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".ts", fagent->get_surface_window(), ".pf", 0);

    // Tb is rocman internal buffer
  std::string::size_type pos = fn_Tb.find( ".");
  COM_assertion_msg(pos!=std::string::npos, "LoadTransfer_FS failed!");
  std::string fn = fn_Tb.substr( 0, pos);
  std::string Tb = fn_Tb.substr( pos, fn_Tb.size());
  fagent->register_clone_dataitem( 0, fn, Tb, fagent->get_surface_window(), ".Tb_alp", 0);
  fagent->register_clone_dataitem( 0, fn, Tb+"_old", fagent->get_surface_window(), ".Tb_alp", 0);
}

void TemperatureTransfer_SF::init( double t) {
  s_Ts_hdl = get_dataitem_handle( 0);
  fb_Tflm_hdl = get_dataitem_handle( 1);
  fn_Tb_hdl = get_dataitem_handle( 2);

  f_Ts_hdl = get_dataitem_handle( fagent->fluidBufNG+".ts");

  bcflag_hdl = get_dataitem_handle( fagent->fluidBufNG+".bcflag");

  load_rocface( fagent->get_coupling()->get_rocmancontrol_param());
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
}

void TemperatureTransfer_SF::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("[%d] Rocstar: TemperatureTransfer_SF::run() with t:%e dt:%e.\n", fagent->get_comm_rank(), t, dt));

  //debug_print(sagent->solidBuf+".Ts", 102, 0);
  COM_call_function( RFC_transfer, &s_Ts_hdl, &f_Ts_hdl);

    // extract burn and non-burn values
  COM_copy_dataitem( fb_Tflm_hdl, f_Ts_hdl, 1, bcflag_hdl, 1);
  COM_copy_dataitem( fn_Tb_hdl, f_Ts_hdl, 1, bcflag_hdl, 0);
  //debug_print(fagent->fluidBufNG+".Tf", 2, 0);
}


// Transfer normal heat flux from fluid to solid
// Arguments: qr_FF (IN), qc_FF (IN), qs_SF (OUT)
// qev: is face centered ('e')
HeatTransfer_FS::HeatTransfer_FS(FluidAgent *fag, SolidAgent *sag,
		 BurnAgent *bag,
                 const std::string f_qc, const std::string f_qr, 
                 const std::string b_qev, const std::string s_qs): 
  InterMeshTransfer(fag, sag, (char *)"HeatTransfer_FS"), bagent(bag)
{
  int io[] = {IN, IN, IN, OUT};
  set_io( 4, io); 
  
  std::string atts[4];
  atts[0] = f_qc;
  atts[1] = f_qr;
  atts[2] = b_qev;
  atts[3] = s_qs;
  set_attr(4, atts);

    // temporary buffer for solid
  sagent->register_clone_dataitem( 0, sagent->solidBuf, ".qc_tmp", sagent->get_surface_window(), ".qs");
  sagent->register_clone_dataitem( 0, sagent->solidBuf, ".qr_tmp", sagent->get_surface_window(), ".qs");
  sagent->register_clone_dataitem( 0, sagent->solidBuf, ".qev", sagent->get_surface_window(), ".qs", 0);
  fagent->register_clone_dataitem( 0, fagent->fluidBufNG, ".qev", fagent->get_surface_window(), ".qc", 0);
}

void HeatTransfer_FS::init( double t) {
  f_qc_hdl = get_dataitem_handle( 0);
  int with_qr = COM_get_dataitem_handle( attr[1]) > 0;
  if (with_qr)
    f_qr_hdl = get_dataitem_handle( 1);
  else 
    f_qr_hdl = -1;
  b_qev_hdl = get_dataitem_handle( attr[2]);
  s_qs_hdl = get_dataitem_handle( 3);

  s_qc_hdl = get_dataitem_handle( sagent->solidBuf+".qc_tmp");
  s_qr_hdl = get_dataitem_handle( sagent->solidBuf+".qr_tmp");
  f_qev_hdl = get_dataitem_handle( fagent->fluidBufNG+".qev");   // burn
  s_qev_hdl = get_dataitem_handle( sagent->solidBuf+".qev");   // burn

     // TODO:  set qev_flag to signal rocburn
  /*
   COM_new_dataitem( attr.c_str(),'w',COM_INT,1,"");
  */
  int one = 1;
  int qev_flag_hdl = get_dataitem_handle( bagent->get_surface_window()+".qev_flag");   // burn
  //COM_call_function( RocBlas::copy_scalar, &one, &qev_flag_hdl);
  int *vm =NULL;
  int strid, cap;
  COM_get_array((bagent->get_surface_window()+".qev_flag").c_str(), 0, &vm, &strid, &cap);
  COM_assertion_msg(vm!=NULL, "Error: qev_flag!");
  *vm = 1;
    
  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");

  MAN_DEBUG(3, ("[%d] Rocstar: HeatTransfer_FS::init() with t:%e with_qr=%d.\n", fagent->get_comm_rank(), t, with_qr));
}

// qc is 'e', qs is 'e'
void HeatTransfer_FS::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("[%d] Rocstar: HeatTransfer_FS::run() with t:%e dt:%e.\n", fagent->get_comm_rank(), t, dt));

  int size_ts = sagent->size_ts;
  int traction_mode = sagent->traction_mode;

  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_assertion_msg(0, "Not implemented!");
  }
  else if ( traction_mode == NO_SHEER && size_ts == 1) {
    //debug_print(fagent->get_surface_window()+".qc", 2, 0);
    //debug_print(fagent->get_surface_window()+".qr", 2, 0);
    if (f_qr_hdl != -1) {
      COM_call_function( RFC_transfer, &f_qc_hdl, &s_qc_hdl);
      COM_call_function( RFC_transfer, &f_qr_hdl, &s_qr_hdl);
      COM_call_function( RocBlas::add, &s_qc_hdl, &s_qr_hdl, &s_qs_hdl);
    }
    else
      COM_call_function( RFC_transfer, &f_qc_hdl, &s_qs_hdl);
    //debug_print(sagent->solidBuf+".qs", 101, 0, bagent->get_communicator(), "s_qs");
      // add qsource
    if (b_qev_hdl == -1) {
      std::cerr << "Rocstar WARNING: No qev is obtained! " << std::endl;
      // only on burn surface
      COM_call_function( RocBlas::neg, &s_qs_hdl, &s_qs_hdl);
    }
    else {
      //debug_print(bagent->iburn_ng+".qev", 102, 0, bagent->get_communicator(), "b_qev");
        // copy qev from burn to fluid, transfer to solid
      COM_copy_dataitem( f_qev_hdl, b_qev_hdl, 0);
      COM_call_function( RFC_transfer, &f_qev_hdl, &s_qev_hdl);
        // only on burn surface
      COM_call_function( RocBlas::sub, &s_qev_hdl, &s_qs_hdl, &s_qs_hdl);
    }
  }
  else if (size_ts == 1) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions must be vectors!");
  }
  else {
    COM_call_function( RFC_transfer, &f_qc_hdl, &s_qc_hdl);
    COM_call_function( RFC_transfer, &f_qr_hdl, &s_qr_hdl);
    COM_call_function( RocBlas::add, &s_qc_hdl, &s_qr_hdl, &s_qs_hdl);
  }
}

// for remesh, reinitialize nc_t0
// called only when -remeshed in init_scheduler
RemeshInit::RemeshInit(FluidAgent *fag, SolidAgent *sag, 
    const std::string s_u, const std::string f_total_disp, 
    const std::string f_nc, const std::string f_nc_t0):
  InterMeshTransfer(fag, sag, (char *)"RemeshInit")
{
  int io[] = {IN, IN, INOUT, OUT};
  set_io( 4, io); 
  
  std::string atts[4];
  atts[0] = s_u;
  atts[1] = f_total_disp;
  atts[2] = f_nc;
  atts[3] = f_nc_t0;
  set_attr(4, atts);

  fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".total_disp", fagent->get_surface_window(), ".du_alp");
  fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".nc_t0", fagent->get_surface_window(), ".nc");
}

void RemeshInit::init( double t) {
  s_u_hdl = get_dataitem_handle(0);
  f_total_disp_hdl = get_dataitem_handle(1);
  f_nc_hdl = get_dataitem_handle(2);
  f_nc_t0_hdl = get_dataitem_handle(3);

  load_rocface(fagent->get_coupling()->get_rocmancontrol_param());
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
}

void RemeshInit::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: calling RemeshInit::run() with t:%e dt:%e alpha:%e.\n", t, dt, alpha));

  COM_call_function( RFC_interpolate, &s_u_hdl, &f_total_disp_hdl);

  COM_call_function( RocBlas::sub, &f_nc_hdl, &f_total_disp_hdl, &f_nc_t0_hdl);
}







