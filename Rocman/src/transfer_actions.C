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

#include "transfer_actions.h"

#include <iostream>

#include "BurnAgent.h"
#include "FluidAgent.h"
#include "RocBlas-SIM.h"
#include "RocstarCoupling.h"
#include "SolidAgent.h"
#include "rocman.h"

COM_EXTERN_MODULE(SurfUtil)
COM_EXTERN_MODULE(SurfX)

static inline void load_rocsurf() {
  if (COM_get_window_handle("SURF") <= 0) {
    COM_LOAD_MODULE_STATIC_DYNAMIC(SurfUtil, "SURF");
  }
}

extern void compute_overlay(FluidAgent *fagent, SolidAgent *sagent, double t);

void _load_rocface(FluidAgent *fagent, SolidAgent *sagent, int rfc_verb) {
  int rank = fagent->get_comm_rank();

  // INIT_ROCFACE of rocman.f90
  if (COM_get_window_handle("RFC") <= 0) {
    if (rank == 0)
      MAN_DEBUG(3, ("Rocstar: load module RocFace.\n"));
    COM_LOAD_MODULE_STATIC_DYNAMIC(SurfX, "RFC");

    int RFC_setv = COM_get_function_handle("RFC.set_verbose");
    COM_call_function(RFC_setv, &rfc_verb);

#if 0
    // if remeshed
    MAN_DEBUG(3, ("Rocstar: remeshed: %d.\n", param->remeshed));
    if (param->remeshed) {
      double t = fagent->get_coupling()->get_control_param()->current_time;
      compute_overlay(fagent, sagent, t);
    }
#endif

    // load overlay
    MPI_Comm comm = fagent->get_communicator();

    int RFC_read = COM_get_function_handle("RFC.read_overlay");
    int fluid_mesh =
        COM_get_dataitem_handle_const(fagent->get_surf_win_i() + ".mesh");
    int solid_mesh = COM_get_dataitem_handle_const(sagent->solidBuf + ".mesh");
    /*
    int fluid_mesh =
        COM_get_dataitem_handle_const(fagent->fluidBufNG + ".mesh");
    int solid_mesh = COM_get_dataitem_handle_const(sagent->solidBuf + ".mesh");
    */
    std::string rfc_dir = "Rocman/" + fagent->get_module_wname() +
                          sagent->get_module_wname() + "/";
    // std::string rfc_dir = "Rocman/RocfloRocsolid/";
    std::string fluid_dir = rfc_dir + "ifluid";
    std::string solid_dir = rfc_dir + "isolid";

    if (rank == 0)
      MAN_DEBUG(2, ("Rocstar: read RocFace overlay mesh from %s.\n",
                    rfc_dir.c_str()));

    /*
    sagent->write_data_files(0.0, "isolid_i", "isolid_i.all");
    */
    /*
    fagent->write_data_files(0.0, "ifluid_i", "ifluid_i.all");
    fagent->write_data_files(0.0, fagent->fluidBufNG,
                             fagent->fluidBufNG + ".all");
    */

    /*
    COM_call_function(RFC_read, &fluid_mesh, &solid_mesh, &comm,
                      fluid_dir.c_str(), solid_dir.c_str(), "HDF");
    */
    /*
    COM_call_function(RFC_read, &fluid_mesh, &solid_mesh, &comm,
                      fluid_dir.c_str(), solid_dir.c_str(), "CGNS");
    */
    COM_call_function(RFC_read, &solid_mesh, &fluid_mesh, &comm,
                      solid_dir.c_str(), fluid_dir.c_str(), "CGNS");

    int RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
    int RFC_interpolate = COM_get_function_handle("RFC.interpolate");
    COM_set_profiling_barrier(RFC_transfer, comm);
    COM_set_profiling_barrier(RFC_interpolate, comm);
  }
}

InterMeshTransfer::InterMeshTransfer(ActionDataList adl, FluidAgent *fag,
                                     SolidAgent *sag, const std::string &name)
    : Action(std::move(adl), name), fagent(fag), sagent(sag) {}

void InterMeshTransfer::load_rocface(int rfc_verb) {
  // INIT_ROCFACE of rocman.f90
  _load_rocface(fagent, sagent, rfc_verb);
}

// Transfer load from fluid to solid
// Arguments: tf_FF (IN), mdot (IN), rb (IN),
// tf_SF (OUT)
LoadTransfer_FS::LoadTransfer_FS(FluidAgent *fag, SolidAgent *sag,
                                 const std::string &f_ts,
                                 const std::string &s_ts,
                                 const std::string &s_pf)
    : InterMeshTransfer({{f_ts, 0, OUT}, {s_ts, 0, OUT}, {s_pf, 0, OUT}}, fag,
                        sag, "LoadTransfer_FS") {
  traction_mode =
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->traction_mode;

  // create_dataitem
  // for SolidAgent ??????????????????
  // if (traction_mode == NO_SHEER && size_ts == 3) {
  sagent->register_new_dataitem(sagent->solidBufBase, ".pf", 'e', COM_DOUBLE, 1,
                                "Pa");
  // }

  // for FluidAgent
  std::string::size_type pos = f_ts.find('.');
  COM_assertion_msg(pos != std::string::npos, "LoadTransfer_FS failed!");
  std::string f = f_ts.substr(0, pos);
  std::string ts = f_ts.substr(pos, f_ts.size());
  if (traction_mode == NO_SHEER) {
    fagent->register_clone_dataitem(false, f, ts, fagent->get_surf_win(),
                                    ".pf");
  } else {
    fagent->register_clone_dataitem(false, f, ts, fagent->get_surf_win(),
                                    ".tf");
  }
}

// Create buffer data
void LoadTransfer_FS::init(double t) {
  // find out size_ts ???
  int dummy;
  std::string unit;
  char loc;
  COM_get_dataitem(sagent->solidBufBase + ".ts", &loc, &dummy, &size_ts, &unit);

  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions "
                         "must be vectors!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);

  // for SolidAgent
  f_ts_hdl = get_dataitem_handle(0);
  s_ts_hdl = get_dataitem_handle(1);
  if (traction_mode == NO_SHEER && size_ts == 3)
    s_pf_hdl = get_dataitem_handle(2);

  f_tf_hdl = COM_get_dataitem_handle_const(fagent->get_surf_win_i() + ".tf");
  f_pf_hdl = COM_get_dataitem_handle_const(fagent->get_surf_win_i() + ".pf");
  /*
  f_tf_hdl = COM_get_dataitem_handle_const(fagent->fluidBufNG + ".tf");
  f_pf_hdl = COM_get_dataitem_handle_const(fagent->fluidBufNG + ".pf");
  */

  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");

  load_rocsurf();
  SURF_compute_face_normals =
      COM_get_function_handle("SURF.compute_element_normals");

  MAN_DEBUG(
      3, ("LoadTransfer_FS::init() called - traction_mode: %d size_ts: %d.\n",
          traction_mode, size_ts));
}

void LoadTransfer_FS::run(double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: LoadTransfer_FS::run() with t:%e dt:%e.\n", t, dt));
  // ts = tf + ( mdot*Vs )

  // part 1  (POST_UPDATE_FLUID in fluid_agent.f90)
  if (traction_mode != NO_SHEER) {
    COM_call_function(RocBlas::copy, &f_tf_hdl, &f_ts_hdl);
  } else { // ts is a scalar
    COM_call_function(RocBlas::copy, &f_pf_hdl, &f_ts_hdl);
  }

  // part 2  (solid_agent.f90 INIT_INBUFF_SOLID())
  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_pf_hdl);
    COM_call_function(SURF_compute_face_normals, &s_ts_hdl);
    COM_call_function(RocBlas::mul, &s_ts_hdl, &s_pf_hdl, &s_ts_hdl);
    COM_call_function(RocBlas::neg, &s_ts_hdl, &s_ts_hdl);
  } else if (traction_mode == NO_SHEER && size_ts == 1) {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  } else if (size_ts == 1) {
    COM_assertion_msg(false, "ERROR: NOT IMPLEMENTED!");
  } else {
    COM_assertion_msg(false, "ERROR: NOT IMPLEMENTED!");
  }
}

// Transfer load from fluid to solid with Burn
// Arguments: tf_FF (IN), mdot (IN), rb (IN),
// tf_SF (OUT)
LoadTransfer_FSc_ALE::LoadTransfer_FSc_ALE(
    FluidAgent *fag, SolidAgent *sag, BurnAgent *bag, const std::string &f_pf,
    const std::string &fb_mdot, const std::string &b_rb,
    const std::string &s_ts, const std::string &s_pf)
    : InterMeshTransfer({{f_pf, 0, IN},
                         {fb_mdot, 0, IN},
                         {b_rb, 0, IN},
                         {s_ts, 0, IN},
                         {s_pf, 0, OUT}},
                        fag, sag, "LoadTransfer_FSc_ALE"),
      bagent(bag) {
  traction_mode =
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->traction_mode;

  // create_dataitem
  // for SolidAgent ??????????????????
  // if (traction_mode == NO_SHEER && size_ts == 3) {
  sagent->register_new_dataitem(sagent->solidBufBase, ".pf", 'e', COM_DOUBLE, 1,
                                "Pa");
  // }

  // for FluidAgent
  if (traction_mode == NO_SHEER) {
    fagent->register_clone_dataitem(false, fagent->get_surf_win_i(), ".ts",
                                    fagent->get_surf_win(), ".pf");
  } else
    fagent->register_clone_dataitem(false, fagent->get_surf_win_i(), ".ts",
                                    fagent->get_surf_win(), ".tf");

  /*
  fagent->register_new_dataitem(fagent->fluidBufB, ".mdot_tmp", 'e', COM_DOUBLE,
                                1, "kg/(m^2 s)");
  */
}

// Create buffer data
void LoadTransfer_FSc_ALE::init(double t) {
  // find out size_ts
  int dummy;
  std::string unit;
  char loc;
  COM_get_dataitem(sagent->solidBufBase + ".ts", &loc, &dummy, &size_ts, &unit);
  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions "
                         "must be vectors!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // f_pf_hdl = COM_get_dataitem_handle_const(fagent->fluidBufNG + ".pf");
  f_pf_hdl = get_dataitem_handle(0);
  fb_mdot_hdl = get_dataitem_handle(1);
  b_rb_hdl = get_dataitem_handle(2);
  s_ts_hdl = get_dataitem_handle(3);
  if (traction_mode == NO_SHEER && size_ts == 3)
    s_pf_hdl = get_dataitem_handle(4);

  f_ts_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".ts");
  f_tf_hdl = COM_get_dataitem_handle_const(fagent->fluidBufNG + ".tf");
  fb_ts_hdl = COM_get_dataitem_handle(fagent->fluidBufB + ".ts");
  fb_pf_hdl = COM_get_dataitem_handle(fagent->fluidBufB + ".pf");

  fb_rhof_alp_hdl = COM_get_dataitem_handle(fagent->fluidBufB + ".rhof_alp");
  fb_mdot_tmp_hdl = COM_get_dataitem_handle(fagent->fluidBufB + ".mdot_tmp");

  fb_nf_alp_hdl = COM_get_dataitem_handle(fagent->fluidBufB + ".nf_alp");
  fb_tf_hdl = COM_get_dataitem_handle_const(fagent->fluidBufB + ".tf");

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");

  load_rocsurf();
  SURF_compute_face_normals =
      COM_get_function_handle("SURF.compute_element_normals");

  MAN_DEBUG(
      3,
      ("LoadTransfer_FSc_ALE::init() called - traction_mode: %d size_ts: %d.\n",
       traction_mode, size_ts));
}

void LoadTransfer_FSc_ALE::run(double t, double dt, double alpha) {
  MAN_DEBUG(3,
            ("Rocstar: LoadTransfer_FSc_ALE::run() with t:%e dt:%e.\n", t, dt));
  // ts = tf + ( mdot*Vs )

#if 1
  // part 1  (POST_UPDATE_FLUID in fluid_agent.f90)
  if (traction_mode != NO_SHEER) {
    COM_call_function(RocBlas::copy, &f_tf_hdl, &f_ts_hdl);

    // Compute tfmts = t_s-t_f (see developers guide).
    // Here fb_mdot_tmp is used as buffer space.
    COM_call_function(RocBlas::div, &fb_mdot_hdl, &fb_rhof_alp_hdl,
                      &fb_mdot_tmp_hdl);

    if (sagent->with_ALE) {
      COM_call_function(RocBlas::sub, &fb_mdot_tmp_hdl, &b_rb_hdl,
                        &fb_mdot_tmp_hdl);
    }

    COM_call_function(RocBlas::mul, &fb_mdot_tmp_hdl, &fb_mdot_hdl,
                      &fb_mdot_tmp_hdl);
    COM_call_function(RocBlas::mul, &fb_nf_alp_hdl, &fb_mdot_tmp_hdl,
                      &fb_ts_hdl);

    // Compute ts = tf - (tf-ts)
    COM_call_function(RocBlas::sub, &fb_tf_hdl, &fb_ts_hdl, &fb_ts_hdl);
  } else { // ts is a scalar
    COM_call_function(RocBlas::copy, &f_pf_hdl, &f_ts_hdl);

    // Compute tsmpf = p_f-t_s (see developers guide).
    // Here fb_mdot_tmp is used as buffer space.
    COM_call_function(RocBlas::div, &fb_mdot_hdl, &fb_rhof_alp_hdl,
                      &fb_mdot_tmp_hdl);

    if (sagent->with_ALE) {
      COM_call_function(RocBlas::sub, &fb_mdot_tmp_hdl, &b_rb_hdl,
                        &fb_mdot_tmp_hdl);
    }

    COM_call_function(RocBlas::mul, &fb_mdot_tmp_hdl, &fb_mdot_hdl, &fb_ts_hdl);

    // Compute ts = pf - (pf-ts)
    COM_call_function(RocBlas::sub, &fb_pf_hdl, &fb_ts_hdl, &fb_ts_hdl);

    // debug_print(fagent->fluidBufB + ".mdot", 102, 0, "LOADTRANSFER");
    // debug_print(fagent->fluidBufB + ".ts", 102, 0, "LOADTRANSFER");

    // Subtract from P_ambient if not zeros
    double P_ambient =
        fagent->get_rocstar_coupling()->get_rocmancontrol_param()->P_ambient;
    if (P_ambient != 0.0) {
      COM_call_function(RocBlas::sub_scalar, &f_ts_hdl, &P_ambient, &f_ts_hdl);
    }
  }

#endif

  // part 2  (solid_agent.f90 INIT_INBUFF_SOLID())
  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_pf_hdl);
    COM_call_function(SURF_compute_face_normals, &s_ts_hdl);
    COM_call_function(RocBlas::mul, &s_ts_hdl, &s_pf_hdl, &s_ts_hdl);
    COM_call_function(RocBlas::neg, &s_ts_hdl, &s_ts_hdl);
  } else if (traction_mode == NO_SHEER && size_ts == 1) {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_ts_hdl);
    // debug_print(sagent->solidBuf+".ts", 102, 0, fagent->get_communicator(),
    // "LOADTRANSFER");
  } else if (size_ts == 1) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions "
                         "must be vectors!");
  } else {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  }
}

// Transfer load from fluid to solid with Burn  (part 2 only)
// Arguments: f_ts (IN)
// s_ts (OUT) tf_SF (OUT)
LoadTransferOnly_FSc_ALE::LoadTransferOnly_FSc_ALE(
    FluidAgent *fag, SolidAgent *sag, BurnAgent *bag, const std::string &f_ts,
    const std::string &s_ts, const std::string &s_pf)
    : InterMeshTransfer({{f_ts, 0, IN}, {s_ts, 0, OUT}, {s_pf, 0, OUT}}, fag,
                        sag, "LoadTransferOnly_FSc_ALE"),
      bagent(bag) {
  traction_mode =
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->traction_mode;

  // create_dataitem
  // for SolidAgent ??????????????????
  // if (traction_mode == NO_SHEER && size_ts == 3) {
  sagent->register_new_dataitem(sagent->solidBufBase, ".pf", 'e', COM_DOUBLE, 1,
                                "Pa");
  // }

  // for FluidAgent
  if (traction_mode == NO_SHEER) {
    fagent->register_clone_dataitem(false, fagent->get_surf_win_i(), ".ts",
                                    fagent->get_surf_win(), ".pf");
  } else {
    fagent->register_clone_dataitem(false, fagent->get_surf_win_i(), ".ts",
                                    fagent->get_surf_win(), ".tf");
  }

  /*
  fagent->register_new_dataitem(fagent->fluidBufB, ".mdot_tmp", 'e', COM_DOUBLE,
                                1, "kg/(m^2 s)");
  */
}

// Create buffer data
void LoadTransferOnly_FSc_ALE::init(double t) {
  // find out size_ts
  int dummy;
  std::string unit;
  char loc;
  COM_get_dataitem(sagent->solidBufBase + ".ts", &loc, &dummy, &size_ts, &unit);
  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions "
                         "must be vectors!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // f_pf_hdl = COM_get_dataitem_handle_const(fagent->fluidBufNG + ".pf");
  f_ts_hdl = get_dataitem_handle(0);
  s_ts_hdl = get_dataitem_handle(1);
  if (traction_mode == NO_SHEER && size_ts == 3)
    s_pf_hdl = get_dataitem_handle(2);

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");

  load_rocsurf();
  SURF_compute_face_normals =
      COM_get_function_handle("SURF.compute_element_normals");

  MAN_DEBUG(3, ("LoadTransferOnly_FSc_ALE::init() called - traction_mode: %d "
                "size_ts: %d.\n",
                traction_mode, size_ts));
}

void LoadTransferOnly_FSc_ALE::run(double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: LoadTransferOnly_FSc_ALE::run() with t:%e dt:%e.\n",
                t, dt));

  // part 2  (solid_agent.f90 INIT_INBUFF_SOLID())
  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_pf_hdl);
    COM_call_function(SURF_compute_face_normals, &s_ts_hdl);
    COM_call_function(RocBlas::mul, &s_ts_hdl, &s_pf_hdl, &s_ts_hdl);
    COM_call_function(RocBlas::neg, &s_ts_hdl, &s_ts_hdl);
  } else if (traction_mode == NO_SHEER && size_ts == 1) {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  } else if (size_ts == 1) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions "
                         "must be vectors!");
  } else {
    COM_call_function(RFC_transfer, &f_ts_hdl, &s_ts_hdl);
  }
}

/*
// Arguments: tf_FF (IN), mdot (IN), rb (IN),
// tf_SF (OUT)
LoadTransfer_FSc_ALE::LoadTransfer_FSc_ALE(FluidAgent *fag, SolidAgent *sag,
                                           const char *at[], int *is, void *p)
    : InterMeshTransfer(fag, sag, "LoadTransfer_FS") {
  int io[] = {IN, OUT};
  set_io(2, io);
}

// Create buffer data
void LoadTransfer_FSc_ALE::init(double t) {}

void LoadTransfer_FSc_ALE::run(double t, double dt, double alpha) {
  // ts = tf + ( mdot*Vs
}
*/

// solid_agent.f90 POST_UPDATE_SOLID
GetDeformedMesh::GetDeformedMesh(FluidAgent *fag, SolidAgent *sag,
                                 const std::string &s_x,
                                 const std::string &s_uhat,
                                 const std::string &s_y)
    : InterMeshTransfer({{s_x, 0, IN}, {s_uhat, 0, IN}, {s_y, 0, OUT}}, fag,
                        sag, "GetDeformedMesh") {
  /*
  // register dataitems
  std::string::size_type pos = s_x_str.find(".");
  COM_assertion_msg(pos != std::string::npos,
                    "GetDeformedMesh::create_dataitem failed!");
  std::string s = s_x_str.substr(0, pos);
  std::string x = s_x_str.substr(pos, s_x_str.size());
  sagent->register_use_dataitem(s, x, sagent->solidBufBase, ".nc");

  pos = s_y_str.find(".");
  COM_assertion_msg(pos != std::string::npos,
                    "GetDeformedMesh::create_dataitem failed!");
  s = s_y_str.substr(0, pos);
  std::string y = s_y_str.substr(pos, s_y_str.size());
  sagent->register_clone_dataitem(false, s, y, sagent->get_surf_win(), ".nc");
  */
}

void GetDeformedMesh::init(double t) {
  s_x_hdl = get_dataitem_handle(0);
  s_uhat_hdl = get_dataitem_handle(1);
  s_y_hdl = get_dataitem_handle(2);
}

// POST_UPDATE_SOLID
void GetDeformedMesh::run(double t, double dt, double alpha) {
  MAN_DEBUG(
      3, ("Rocstar: calling GetDeformedMesh::run() with t:%e dt:%e alpha:%e.\n",
          t, dt, alpha));
  COM_call_function(RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
}

// solid_agent.f90 POST_UPDATE_SOLID
GetDeformedMesh_ALE::GetDeformedMesh_ALE(FluidAgent *fag, SolidAgent *sag,
                                         const std::string &s_x,
                                         const std::string &s_uhat,
                                         const std::string &s_y, double z)
    : InterMeshTransfer({{s_x, 0, IN}, {s_uhat, 0, IN}, {s_y, 0, OUT}}, fag,
                        sag, "GetDeformedMesh_ALE"),
      zoom(z) {
  /*
  // register dataitems
  std::string::size_type pos = s_x.find('.');
  COM_assertion_msg(pos != std::string::npos,
                    "GetDeformedMesh_ALE::create_dataitem failed!");
  std::string s = s_x.substr(0, pos);
  std::string x = s_x.substr(pos, s_x.size());
  sagent->register_use_dataitem(s, x, sagent->solidBufBase, ".nc");

  // s_y = solidBuf + ".nc";
  pos = s_y.find('.');
  COM_assertion_msg(pos != std::string::npos,
                    "GetDeformedMesh_ALE::create_dataitem failed!");
  s = s_y.substr(0, pos);
  std::string y = s_y.substr(pos, s_y.size());
  sagent->register_clone_dataitem(false, s, y, sagent->get_surf_win(), y);
  */
}

void GetDeformedMesh_ALE::init(double t) {
  s_x_hdl = get_dataitem_handle(0);
  s_uhat_hdl = get_dataitem_handle(1);
  s_y_hdl = get_dataitem_handle(2);

  with_ALE = sagent->with_ALE;

  s_rb_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".rb");
  s_mdot_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".mdot");
  s_areas_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".areas");
  s_rhos_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".rhos");

  load_rocsurf();
  SURF_compute_bounded_volumes =
      COM_get_function_handle("SURF.compute_bounded_volumes");
  SURF_compute_face_areas =
      COM_get_function_handle("SURF.compute_element_areas");
}

// POST_UPDATE_SOLID
void GetDeformedMesh_ALE::run(double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: calling GetDeformedMesh_ALE::run() with t:%e dt:%e "
                "alpha:%e with_ALE=%d zoom=%e.\n",
                t, dt, alpha, with_ALE, zoom));

  if (with_ALE && zoom > 0 &&
      !fagent->get_rocstar_coupling()->get_rocmancontrol_param()->PROP_fom) {
    // compute volumes only for faces with nonzero burning rate
    COM_call_function(RocBlas::copy, &s_rb_hdl, &s_mdot_hdl);

    // debug_print(sagent->solidBuf + ".nc", 202, 1);
    // debug_print(sagent->solidBuf + ".x", 202, 1);

    // Compute mass injection using s_x and s_y as buffers.
    // So far s_y stores the undeformed configuration s_x.
    int one = 1;
    COM_call_function(SURF_compute_bounded_volumes, &s_x_hdl, &s_y_hdl,
                      &s_mdot_hdl, &one);

    // debug_print(sagent->solidBuf + ".nc", 202, 1);
    // debug_print(attr[2], 202, 1);

#if ENABLE_DEBUG
    // Check whether the burning rates are nonnegative.
    double mdot_min{0.0};
    MPI_Comm comm = MPI_COMM_SELF;
    COM_call_function(RocBlas::min_scalar_MPI, &s_mdot_hdl, &mdot_min, &comm);

    if (mdot_min < 0) {
      printf("Rocstar ERROR: Negative mdot found %e. Aborting...\n", mdot_min);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
#endif

    double zero = 0.0;
    COM_call_function(RocBlas::maxof_scalar, &s_mdot_hdl, &zero, &s_mdot_hdl);

    // Compute deformed configuration
    COM_call_function(RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
    COM_call_function(SURF_compute_face_areas, &s_areas_hdl);

    // Compute mdot=dV*rhos/area/dt
    COM_call_function(RocBlas::mul, &s_mdot_hdl, &s_rhos_hdl, &s_mdot_hdl);
    COM_call_function(RocBlas::div, &s_mdot_hdl, &s_areas_hdl, &s_mdot_hdl);
    COM_call_function(RocBlas::div_scalar, &s_mdot_hdl, &dt, &s_mdot_hdl);
  } else {
    // Compute deformed configuration
    COM_call_function(RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);
    // Compute mdot as rb*rhos
    if (with_ALE)
      COM_call_function(RocBlas::mul, &s_rb_hdl, &s_rhos_hdl, &s_mdot_hdl);
  }
}

// fluid_agent.f90  INIT_INBUFF_FLUID
MeshMotionTransfer_SF::MeshMotionTransfer_SF(FluidAgent *fag, SolidAgent *sag,
                                             const std::string &s_u,
                                             const std::string &f_total_disp,
                                             const std::string &f_vm)
    : InterMeshTransfer({{sag->solidBufBase + ".nc", 0, IN},
                         {s_u, 0, IN},
                         {f_total_disp, 0, OUT},
                         {f_vm, 0, OUT}},
                        fag, sag, "MeshMotionTransfer_SF") {
  // isolid_i+".nc" (in), isolid_i+".u" (in), fagent->ifluid_i+".vm" (out)

  // Change to transfer total_displacement from solid to fluids, and obtain
  // incremental displacement as (nc_t0+total_disp)-nc_tn. This is prone
  // to cancellation errors but avoids accumulation of any errors.
  fagent->register_clone_dataitem(false, fagent->get_surf_win_i(),
                                  ".total_disp", fagent->get_surf_win(),
                                  ".du_alp");
  fagent->register_clone_dataitem(false, fagent->get_surf_win_i(), ".nc_t0",
                                  fagent->get_surf_win(), ".nc");

  // nc_tmp is used to compare against solid nodal coordinates.
  fagent->register_clone_dataitem(false, fagent->get_surf_win_i(), ".nc_tmp",
                                  fagent->get_surf_win(), ".nc");
  fagent->register_new_dataitem(fagent->get_surf_win_i(), ".sq_dist", 'n',
                                COM_DOUBLE, 1, "m");
}

void MeshMotionTransfer_SF::init(double t) {
  // skip 0
  s_u_hdl = get_dataitem_handle(1);
  f_total_disp_hdl = get_dataitem_handle(2);
  f_vm_hdl = get_dataitem_handle(3);

  f_nc_hdl = COM_get_dataitem_handle(fagent->get_surf_win_i() + ".nc");
  f_nc_t0_hdl = COM_get_dataitem_handle(fagent->get_surf_win_i() + ".nc_t0");
  /*
  f_nc_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".nc");
  f_nc_t0_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".nc_t0");
  */

  /*
  // POST_INIT_FLUID in fluid_agent.f90
  // initial_start  ?????
  if (t == 0.0) {
    // TODO ???
    COM_call_function(RocBlas::copy, &f_nc_hdl, &f_nc_t0_hdl);
    // UPDATE_INBUFF_GM_FLUID     ???????????????
    double zero = 0.0;
    fagent->obtain_gm(&zero);
  }
  */

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
}

// INIT_INBUFF_FLUID
void MeshMotionTransfer_SF::run(double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: calling MeshMotionTransfer_SF::run() with t:%e dt:%e "
                "alpha:%e.\n",
                t, dt, alpha));

  // b. Interpolates total displacement from solid nodes to fluid nodes
  /*
  debug_print(sagent->solidBuf + ".u", 109, 0, sagent->get_communicator(), "U");
  */
  COM_call_function(RFC_interpolate, &s_u_hdl, &f_total_disp_hdl);

  // Compute incremental mesh velocity
  if (!fagent->get_coupling()->is_initial_start()) {
    // Compute incremental mesh velocity
    COM_call_function(RocBlas::add, &f_nc_t0_hdl, &f_total_disp_hdl, &f_vm_hdl);
    COM_call_function(RocBlas::sub, &f_vm_hdl, &f_nc_hdl, &f_vm_hdl);
    COM_call_function(RocBlas::div_scalar, &f_vm_hdl, &dt, &f_vm_hdl);
  }
}

// transfer vs from solid nodes to fluid faces using Rocface
DeformationVelTransfer_SF::DeformationVelTransfer_SF(FluidAgent *fag,
                                                     SolidAgent *sag,
                                                     const std::string &s_vs,
                                                     const std::string &f_vs)
    : InterMeshTransfer({{fag->get_surf_win_i() + ".vm", 0, IN},
                         {s_vs, 0, IN},
                         {f_vs, 0, OUT}},
                        fag, sag, "DeformationVelTransfer_SF") {
  // added a builtin in (ifluid_i+".vm")
  // ifluid_i+".vm" (in),  isolid_i+".vs", ifluid_i+".vs"
}

void DeformationVelTransfer_SF::init(double t) {
  s_vs_hdl = get_dataitem_handle(1);
  f_vs_hdl = get_dataitem_handle(2);

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);

  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void DeformationVelTransfer_SF::run(double t, double dt, double alpha) {
  MAN_DEBUG(
      3,
      ("Rocstar: DeformationVelTransfer_SF::run() with t:%e dt:%e alpha:%e.\n",
       t, dt, alpha));

  // c. Transfer vs from solid nodes to fluid faces using Rocface
  COM_call_function(RFC_transfer, &s_vs_hdl, &f_vs_hdl);
}

MeshMotionTransferISS::MeshMotionTransferISS(FluidAgent *fag, SolidAgent *sag,
                                             const std::string &s_u,
                                             const std::string &s_vs,
                                             const std::string &f_vm)
    : InterMeshTransfer({{s_u, 0, IN}, {s_vs, 0, IN}, {f_vm, 0, OUT}}, fag, sag,
                        "MeshMotionTransferISS") {
  // isolid_i+".u" (in), isolid_i+".vs" (in), ifluid_i+".vm" (out)
  std::string::size_type pos = f_vm.find('.');
  COM_assertion_msg(pos != std::string::npos, "MeshMotionTransferISS failed!");
  std::string f = f_vm.substr(0, pos);
  std::string vm = f_vm.substr(pos, f_vm.size());
  fagent->register_new_dataitem(f, vm, 'n', COM_DOUBLE, 3, "m/s");
  // used internally
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".u",
                                  fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".uold",
                                  fagent->fluidBufNG, ".vm");
  // avoid vs because vs is 'e'
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".fvs",
                                  fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".vsold",
                                  fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".utmp",
                                  fagent->fluidBufNG, ".vm");
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".vtmp",
                                  fagent->fluidBufNG, ".vm");
}

void MeshMotionTransferISS::init(double t) {
  s_u_hdl = get_dataitem_handle(0);
  s_vs_hdl = get_dataitem_handle(1);
  f_vm_hdl = get_dataitem_handle(2);

  f_u_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".u");
  f_vs_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".fvs");
  f_uold_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".uold");
  f_vsold_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".vsold");
  f_utmp_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".utmp");
  f_vtmp_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".vtmp");

  // set zero
  double zero = 0.0;
  COM_call_function(RocBlas::copy_scalar, &zero, &f_uold_hdl);
  COM_call_function(RocBlas::copy_scalar, &zero, &f_vsold_hdl);

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
}

// INIT_INBUFF_FLUID
void MeshMotionTransferISS::run(double t, double dt, double alpha) {
  MAN_DEBUG(3,
            ("Rocstar: calling MeshMotionTransferISS::run() with t:%e dt:%e.\n",
             t, dt));
  // formula:  Vm = (un-un-1)/deltaT + (vsn - vsn-1) /2

  COM_call_function(RFC_interpolate, &s_u_hdl, &f_u_hdl);
  COM_call_function(RFC_interpolate, &s_vs_hdl, &f_vs_hdl);

  COM_call_function(RocBlas::sub, &f_u_hdl, &f_uold_hdl, &f_utmp_hdl);

  COM_call_function(RocBlas::div_scalar, &f_utmp_hdl, &dt, &f_utmp_hdl);

  COM_call_function(RocBlas::sub, &f_vs_hdl, &f_vsold_hdl, &f_vtmp_hdl);

  double half = 0.5;
  COM_call_function(RocBlas::mul_scalar, &f_vtmp_hdl, &half, &f_vtmp_hdl);

  COM_call_function(RocBlas::add, &f_utmp_hdl, &f_vtmp_hdl, &f_vm_hdl);

  // save copies
  COM_call_function(RocBlas::copy, &f_u_hdl, &f_uold_hdl);
  COM_call_function(RocBlas::copy, &f_vs_hdl, &f_vsold_hdl);
}

// burn_agent.f90  INIT_INBUFF_BURN
// Transfer rhos from solid faces to fluid faces
TransferSolidDensity::TransferSolidDensity(FluidAgent *fag, SolidAgent *sag,
                                           const std::string &s_rhos,
                                           const std::string &f_rhos)
    : InterMeshTransfer({{s_rhos, 0, IN}, {f_rhos, 0, OUT}}, fag, sag,
                        "TransferSolidDensity") {
  fagent->register_new_dataitem(fagent->fluidBufNG, ".rhos", 'e', COM_DOUBLE, 1,
                                "kg/(m^3)");
  /*
  sagent->register_use_dataitem(sagent->solidBuf, ".x", sagent->solidBufBase,
                                ".nc");
  sagent->register_clone_dataitem(false, sagent->solidBuf, ".nc",
                                  sagent->get_surf_win(), ".nc");
  */
}

void TransferSolidDensity::init(double t) {
  MAN_DEBUG(3, ("Rocstar: TransferSolidDensity::init() with t:%e.\n", t));
  s_rhos_hdl = get_dataitem_handle(0);
  f_rhos_hdl = get_dataitem_handle(1);

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);

  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void TransferSolidDensity::run(double t, double dt, double alpha) {
  int rhos_mode = sagent->rhos_mode;

  MAN_DEBUG(3, ("Rocstar: TransferSolidDensity::run() with rhs_mode=%d t:%e "
                "dt:%e alpha:%e.\n",
                rhos_mode, t, dt, alpha));

  if (rhos_mode == 1) { // Constant density
    COM_call_function(RocBlas::copy, &s_rhos_hdl, &f_rhos_hdl);
  } else if (fagent->get_coupling()->is_initial_start()) {
    COM_call_function(RFC_transfer, &s_rhos_hdl, &f_rhos_hdl);
  } else if (rhos_mode == 3) { // Varying density
    COM_call_function(RFC_transfer, &s_rhos_hdl, &f_rhos_hdl);
  }
}

// solid_agent.f90  INIT_INBUFF_SOLID
// Transfer rb (burn rate) from fluid faces to solid faces
TransferBurnRate_FS_ALE::TransferBurnRate_FS_ALE(FluidAgent *fag,
                                                 SolidAgent *sag,
                                                 const std::string &b_rb,
                                                 const std::string &s_rb)
    : InterMeshTransfer({{b_rb, 0, IN}, {s_rb, 0, OUT}}, fag, sag,
                        "TransferBurnRate_FS_ALE") {}

void TransferBurnRate_FS_ALE::init(double t) {
  if (!sagent->with_ALE)
    return;

  b_rb_hdl = get_dataitem_handle(0);
  s_rb_hdl = get_dataitem_handle(1);
  int PROP_fom =
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->PROP_fom;
  if (PROP_fom) {
    p_rb_hdl = COM_get_dataitem_handle(sagent->propBufAll + ".rb");
  } else {
    p_rb_hdl = COM_get_dataitem_handle(sagent->propBuf + ".rb");
  }
  f_mdot_tmp_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".mdot_tmp");
  fb_mdot_tmp_hdl = COM_get_dataitem_handle(fagent->fluidBufB + ".mdot_tmp");

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void TransferBurnRate_FS_ALE::run(double t, double dt, double alpha) {
  if (!sagent->with_ALE)
    return;

  MAN_DEBUG(
      3, ("Rocstar: TransferBurnRate_FS_ALE::run() with t:%e dt:%e alpha:%e.\n",
          t, dt, alpha));

  // Transfer rb from fluid faces to solid faces
  double zero = 0.0;
  COM_call_function(RocBlas::copy_scalar, &zero, &p_rb_hdl);

  // Here we use f_mdot_tmp as buffer space for holding f_rb
  COM_call_function(RocBlas::copy_scalar, &zero, &f_mdot_tmp_hdl);

  COM_call_function(RocBlas::copy, &b_rb_hdl, &fb_mdot_tmp_hdl);

  COM_call_function(RFC_transfer, &f_mdot_tmp_hdl, &s_rb_hdl);
}

// INIT_INBUFF_FLUID
// Depending on whether solid has ALE or not, we need to compute mdot
// differently. If solid has ALE, to ensure mass conservation, we compute
// mdot using a geometric construction on the solid side and transfer
// mdot to fluids. If solid has no ALE, mass conservation is not possible
// and the above approach will generate zero mdot, so we have to compute
// mdot as rhos*rb in this case.
MassTransfer_SF_ALE::MassTransfer_SF_ALE(FluidAgent *fag, SolidAgent *sag,
                                         BurnAgent *bag,
                                         const std::string &f_mdot)
    : InterMeshTransfer({{f_mdot, 0, OUT}}, fag, sag, "MassTransfer_SF_ALE"),
      bagent(bag) {
  fagent->register_new_dataitem(fagent->get_surf_win_i(), ".rhos", 'e',
                                COM_DOUBLE, 1, "kg/(m^3)");
  fagent->register_new_dataitem(fagent->get_surf_win_i(), ".mdot", 'e',
                                COM_DOUBLE, 1, "kg/(m^2 s)");
  bagent->register_use_dataitem(bagent->get_surf_win(), ".rhos",
                                bagent->parentWin, ".rhos");
}

void MassTransfer_SF_ALE::init(double t) {
  // f_mdot_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".mdot");
  f_mdot_hdl = get_dataitem_handle(0);

  b_rhos_hdl = COM_get_dataitem_handle(bagent->iburn_ng + ".rhos");
  b_rb_hdl = COM_get_dataitem_handle(bagent->iburn_ng + ".rb");
  fb_mdot_hdl = COM_get_dataitem_handle(fagent->fluidBufB + ".mdot");

  // ALE
  f_total_disp_hdl =
      COM_get_dataitem_handle(fagent->fluidBufNG + ".total_disp");
  f_nc_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".nc");
  f_nc_t0_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".nc_t0");
  s_mdot_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".mdot");

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void MassTransfer_SF_ALE::run(double t, double dt, double alpha) {
  MAN_DEBUG(3, ("[%d] Rocstar: MassTransfer_SF_ALE::run() with t:%e dt:%e.\n",
                fagent->get_comm_rank(), t, dt));

  if (sagent->with_ALE) {
    COM_call_function(RocBlas::add, &f_nc_t0_hdl, &f_total_disp_hdl, &f_nc_hdl);
    COM_call_function(RFC_transfer, &s_mdot_hdl, &f_mdot_hdl);
  } else {
    double zero = 0.0;
    COM_call_function(RocBlas::copy_scalar, &zero, &f_mdot_hdl);
    COM_call_function(RocBlas::mul, &b_rhos_hdl, &b_rb_hdl, &fb_mdot_hdl);
  }
}

// transfer temperature from solid to fluid
// s_Ts in solid is 'n'  f_Tf in fluid is 'e'
TemperatureTransfer_SF::TemperatureTransfer_SF(SolidAgent *sag, FluidAgent *fag,
                                               const std::string &s_Ts,
                                               const std::string &fb_Tflm,
                                               const std::string &fn_Tb)
    : InterMeshTransfer({{s_Ts, 0, IN}, {fb_Tflm, 0, OUT}, {fn_Tb, 0, OUT}},
                        fag, sag, "TemperatureTransfer_SF") {
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".ts",
                                  fagent->get_surf_win(), ".pf", 0);

  // Tb is rocman internal buffer
  std::string::size_type pos = fn_Tb.find('.');
  COM_assertion_msg(pos != std::string::npos, "LoadTransfer_FS failed!");
  std::string fn = fn_Tb.substr(0, pos);
  std::string Tb = fn_Tb.substr(pos, fn_Tb.size());
  fagent->register_clone_dataitem(false, fn, Tb, fagent->get_surf_win(),
                                  ".Tb_alp", 0);
  fagent->register_clone_dataitem(false, fn, Tb + "_old",
                                  fagent->get_surf_win(), ".Tb_alp", 0);
}

void TemperatureTransfer_SF::init(double t) {
  s_Ts_hdl = get_dataitem_handle(0);
  fb_Tflm_hdl = get_dataitem_handle(1);
  fn_Tb_hdl = get_dataitem_handle(2);

  f_Ts_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".ts");

  bcflag_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".bcflag");

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
}

void TemperatureTransfer_SF::run(double t, double dt, double alpha) {
  MAN_DEBUG(3,
            ("[%d] Rocstar: TemperatureTransfer_SF::run() with t:%e dt:%e.\n",
             fagent->get_comm_rank(), t, dt));

  // debug_print(sagent->solidBuf + ".Ts", 102, 0);
  COM_call_function(RFC_transfer, &s_Ts_hdl, &f_Ts_hdl);

  // extract burn and non-burn values
  COM_copy_dataitem(fb_Tflm_hdl, f_Ts_hdl, 1, bcflag_hdl, 1);
  COM_copy_dataitem(fn_Tb_hdl, f_Ts_hdl, 1, bcflag_hdl, 0);
  // debug_print(fagent->fluidBufNG + ".Tf", 2, 0);
}

// Transfer normal heat flux from fluid to solid
// Arguments: qr_FF (IN), qc_FF (IN), qs_SF (OUT)
// qev: is face centered ('e')
HeatTransfer_FS::HeatTransfer_FS(FluidAgent *fag, SolidAgent *sag,
                                 BurnAgent *bag, const std::string &f_qc,
                                 const std::string &f_qr,
                                 const std::string &b_qev,
                                 const std::string &s_qs)
    : InterMeshTransfer(
          {{f_qc, 0, IN}, {f_qr, 0, IN}, {b_qev, 0, IN}, {s_qs, 0, OUT}}, fag,
          sag, "HeatTransfer_FS"),
      bagent(bag) {
  // temporary buffer for solid
  sagent->register_clone_dataitem(false, sagent->solidBuf, ".qc_tmp",
                                  sagent->get_surf_win(), ".qs");
  sagent->register_clone_dataitem(false, sagent->solidBuf, ".qr_tmp",
                                  sagent->get_surf_win(), ".qs");
  sagent->register_clone_dataitem(false, sagent->solidBuf, ".qev",
                                  sagent->get_surf_win(), ".qs", 0);
  fagent->register_clone_dataitem(false, fagent->fluidBufNG, ".qev",
                                  fagent->get_surf_win(), ".qc", 0);
}

void HeatTransfer_FS::init(double t) {
  f_qc_hdl = get_dataitem_handle(0);
  int with_qr = get_dataitem_handle(1) > 0;
  if (with_qr)
    f_qr_hdl = get_dataitem_handle(1);
  b_qev_hdl = get_dataitem_handle(2);
  s_qs_hdl = get_dataitem_handle(3);

  s_qc_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".qc_tmp");
  s_qr_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".qr_tmp");
  f_qev_hdl = COM_get_dataitem_handle(fagent->fluidBufNG + ".qev"); // burn
  s_qev_hdl = COM_get_dataitem_handle(sagent->solidBuf + ".qev");   // burn

  // TODO:  set qev_flag to signal rocburn
  /*
  COM_new_dataitem( attr.c_str(),'w',COM_INT,1,"");
  // int one = 1;
  // int qev_flag_hdl = COM_get_dataitem_handle(bagent->get_surf_win() +
  // ".qev_flag");   // burn COM_call_function( RocBlas::copy_scalar, &one,
  // &qev_flag_hdl);
  */
  int *vm = nullptr;
  int strid, cap;
  COM_get_array((bagent->get_surf_win() + ".qev_flag").c_str(), 0, &vm, &strid,
                &cap);
  COM_assertion_msg(vm != nullptr, "Error: qev_flag!");
  *vm = 1;

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");

  MAN_DEBUG(3, ("[%d] Rocstar: HeatTransfer_FS::init() with t:%e with_qr=%d.\n",
                fagent->get_comm_rank(), t, with_qr));
}

// qc is 'e', qs is 'e'
void HeatTransfer_FS::run(double t, double dt, double alpha) {
  MAN_DEBUG(3, ("[%d] Rocstar: HeatTransfer_FS::run() with t:%e dt:%e.\n",
                fagent->get_comm_rank(), t, dt));

  int size_ts = sagent->size_ts;
  int traction_mode = sagent->traction_mode;

  if (traction_mode == NO_SHEER && size_ts == 3) {
    COM_assertion_msg(0, "Not implemented!");
  } else if (traction_mode == NO_SHEER && size_ts == 1) {
    // debug_print(fagent->get_surf_win() + ".qc", 2, 0);
    // debug_print(fagent->get_surf_win() + ".qr", 2, 0);
    if (f_qr_hdl != -1) {
      COM_call_function(RFC_transfer, &f_qc_hdl, &s_qc_hdl);
      COM_call_function(RFC_transfer, &f_qr_hdl, &s_qr_hdl);
      COM_call_function(RocBlas::add, &s_qc_hdl, &s_qr_hdl, &s_qs_hdl);
    } else {
      COM_call_function(RFC_transfer, &f_qc_hdl, &s_qs_hdl);
    }
    /*
    debug_print(sagent->solidBuf + ".qs", 101, 0, bagent->get_communicator(),
                "s_qs");
    */
    // add qsource
    if (b_qev_hdl == -1) {
      std::cerr << "Rocstar WARNING: No qev is obtained! " << std::endl;
      // only on burn surface
      COM_call_function(RocBlas::neg, &s_qs_hdl, &s_qs_hdl);
    } else {
      /*
      debug_print(bagent->iburn_ng + ".qev", 102, 0, bagent->get_communicator(),
                  "b_qev");
      */
      // copy qev from burn to fluid, transfer to solid
      COM_copy_dataitem(f_qev_hdl, b_qev_hdl, 0);
      COM_call_function(RFC_transfer, &f_qev_hdl, &s_qev_hdl);
      // only on burn surface
      COM_call_function(RocBlas::sub, &s_qev_hdl, &s_qs_hdl, &s_qs_hdl);
    }
  } else if (size_ts == 1) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions "
                         "must be vectors!");
  } else {
    COM_call_function(RFC_transfer, &f_qc_hdl, &s_qc_hdl);
    COM_call_function(RFC_transfer, &f_qr_hdl, &s_qr_hdl);
    COM_call_function(RocBlas::add, &s_qc_hdl, &s_qr_hdl, &s_qs_hdl);
  }
}

// for remesh, reinitialize nc_t0
// called only when -remeshed in init_scheduler
RemeshInit::RemeshInit(FluidAgent *fag, SolidAgent *sag, const std::string &s_u,
                       const std::string &f_total_disp, const std::string &f_nc,
                       const std::string &f_nc_t0)
    : InterMeshTransfer({{s_u, 0, IN},
                         {f_total_disp, 0, IN},
                         {f_nc, 0, IN},
                         {f_nc_t0, 0, OUT}},
                        fag, sag, "RemeshInit") {
  fagent->register_clone_dataitem(false, fagent->get_surf_win_i(),
                                  ".total_disp", fagent->get_surf_win(),
                                  ".du_alp");
  fagent->register_clone_dataitem(false, fagent->get_surf_win_i(), ".nc_t0",
                                  fagent->get_surf_win(), ".nc");
}

void RemeshInit::init(double t) {
  s_u_hdl = get_dataitem_handle(0);
  f_total_disp_hdl = get_dataitem_handle(1);
  f_nc_hdl = get_dataitem_handle(2);
  f_nc_t0_hdl = get_dataitem_handle(3);

  load_rocface(
      fagent->get_rocstar_coupling()->get_rocmancontrol_param()->rfc_verb);
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
}

void RemeshInit::run(double t, double dt, double alpha) {
  MAN_DEBUG(3,
            ("Rocstar: calling RemeshInit::run() with t:%e dt:%e alpha:%e.\n",
             t, dt, alpha));

  COM_call_function(RFC_interpolate, &s_u_hdl, &f_total_disp_hdl);

  COM_call_function(RocBlas::sub, &f_nc_hdl, &f_total_disp_hdl, &f_nc_t0_hdl);
}
