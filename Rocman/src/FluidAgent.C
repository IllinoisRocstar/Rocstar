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
// $Id: FluidAgent.C,v 1.69 2009/09/14 14:18:27 mtcampbe Exp $

#include <utility>

#include "FluidAgent.h"

#include "Control_parameters.h"
#include "RocBlas.h"
#include "RocstarCoupling.h"
#include "rocman.h"

#ifdef STATIC_LINK
extern "C" void COM_F_FUNC2(rocflo_load_module,
                            ROCFLO_LOAD_MODULE)(const char *, int);
extern "C" void COM_F_FUNC2(rocflo_unload_module,
                            ROCFLO_UNLOAD_MODULE)(const char *, int);
extern "C" void COM_F_FUNC2(rocflu_load_module,
                            ROCFLU_LOAD_MODULE)(const char *, int);
extern "C" void COM_F_FUNC2(rocflu_unload_module,
                            ROCFLU_UNLOAD_MODULE)(const char *, int);
#endif

FluidAgent::FluidAgent(RocstarCoupling *cp, std::string mod, std::string obj,
                       MPI_Comm com, bool withSolid)
    : RocstarAgent(cp, std::move(mod), std::move(obj), "FluidAgent", "fluid",
                   "ifluid", com, true),
      with_plag(false), with_solid(withSolid) {
  fluid_plag = "fluid_plag";

  // SimIN windows, Input buffer windows
  fluidPlagIn = "fluidPlagIn";
  fluidVPIn = "fluidVolIn fluidPlagIn";

  ifluid_i = surf_i; // FluidBuf

  propBufAll = "FluidPropAll";
  fluidBufNG = "FluidBufNG"; // fluidBufNG
  propBuf = "FluidProp";
  fluidBufB = "FluidBufB";
  fluidBufNB = "FluidBufNB";

  fluidBufPC = "FluidBufPC";   // surface data to be backed up
  fluidBufBak = "fluidBufBak"; // back up PC iterations
  fluidBufPRE = "FluidBufPRE"; //  for convergence check

  fluidVolBak = "FluidVolBak";   // backup of volume data
  fluidPlagBak = "FluidPlagBak"; // backup of volume data

  if (!with_solid)
    dobackup = false;
}

/* AEG: Functionality moved to parent class. The STATIC_LINK flag is now in
 * IMPACT. Since load_module() is called during constructing of the parent, it
 * cannot be virtual and defined in the derived. *//*
void FluidAgent::load_module() {
#ifdef STATIC_LINK
  MAN_DEBUG(3, ("[%d] Rocstar: FluidAgent::load_module %s %s.\n", comm_rank,
                module_lname.c_str(), module_wname.c_str()));
#ifdef RFLO
  if (module_lname == "Rocflo")
    COM_F_FUNC2(rocflo_load_module, ROCFLO_LOAD_MODULE)
    (module_wname.c_str(), module_wname.length());
#elif defined(RFLU)
  if (module_lname == "Rocflu")
    COM_F_FUNC2(rocflu_load_module, ROCFLU_LOAD_MODULE)
    (module_wname.c_str(), module_wname.length());
#else
  COM_assertion_msg(0, "Unknown FluidAgent mod!");
#endif
#endif
}

void FluidAgent::unload_module() {
#ifdef STATIC_LINK
  MAN_DEBUG(3,
            ("Rocstar: FluidAgent::unload_module %s.\n", module_lname.c_str()));
#ifdef RFLO
  if (module_lname == "Rocflo")
    COM_F_FUNC2(rocflo_unload_module, ROCFLO_UNLOAD_MODULE)
    (module_wname.c_str(), module_wname.length());
#else
#ifdef RFLU
  if (module_lname == "Rocflu")
    COM_F_FUNC2(rocflu_unload_module, ROCFLU_UNLOAD_MODULE)
    (module_wname.c_str(), module_wname.length());
#endif
#endif
#endif
}
*/

/* AEG: Functionality moved to parent class. However, FluidAgent was jury-rigged
 * to support Rocflu/o's PLAG. This support is now disabled. The additional code
 * snippets are kept below. *//*
void FluidAgent::input(double t) {
#ifndef NATIVE_MP_IO
  if (read_by_control_file(t, fluid_plag, fluidPlagIn) == -1) {
    COM_new_window(fluidPlagIn);
    COM_window_init_done(fluidPlagIn);
  }
#else
  COM_new_window(fluidPlagIn);
  COM_window_init_done(fluidPlagIn);
#endif
}

void FluidAgent::init_module(double t) {
#ifndef NATIVE_MP_IO
  COM_delete_window(fluidPlagIn);
#endif
  // INITIALIZE_XXX common portion

  // Reassign vol_window
  std::string::size_type pos = vol_window.find(' ');
  if (pos != std::string::npos) {
    plag_window = vol_window.substr(pos + 1, vol_window.size());
    vol_window = vol_window.substr(0, pos);
  }
  with_plag = !plag_window.empty();
#ifdef NATIVE_MP_IO
  with_plag = false;
#endif
}

void FluidAgent::output_restart_files(double t) {
  // Write out Plag window
  if (with_plag) {
    write_data_files(t, fluid_plag, plag_window + ".all");
    write_control_file(t, fluid_plag, plag_window);
  }
}
*/

// called in init callback (ic_handle)
void FluidAgent::create_buffer_all() {
  std::string bcflag = surf_window + ".bcflag";

  /**
   * propBufAll is a USE of surf_window mesh with NEW data:
   *   .vm, .rb, .positions, .cflag
   * If .cnstr_type exists on surf_window, then also USE:
   *   .bflag, .cnstr_type, .bcflag,
   * and make a CLONE of .cnstr_type into .cnstr_type_bak
   */
  COM_new_window(propBufAll);
  COM_use_dataitem(propBufAll, surf_window + ".mesh");
  // COM_use_dataitem(propBufAll, surf_bcflag);
  // Mesh motion velocities
  COM_new_dataitem(propBufAll + ".vm", 'n', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem(propBufAll + ".rb", 'e', COM_DOUBLE, 1, "m/s");
  COM_new_dataitem(propBufAll + ".positions", 'n', COM_DOUBLE, 3, "m");
  COM_new_dataitem(propBufAll + ".cflag", 'n', COM_INTEGER, 1, "");
  COM_resize_array(propBufAll + ".data");
  if (COM_get_dataitem_handle(surf_window + ".cnstr_type") > 0) {
    COM_use_dataitem(propBufAll, surf_window + ".bflag");
    COM_use_dataitem(propBufAll, surf_window + ".cnstr_type");
    COM_use_dataitem(propBufAll, bcflag);
    // a temporary buffer to store ".cnstr_type" of fluid at restart
    COM_clone_dataitem(propBufAll + ".cnstr_type_bak",
                       surf_window + ".cnstr_type");
  }
  create_registered_window_dataitems(propBufAll);
  COM_window_init_done(propBufAll);

  // create a window for surface propagation if run in fluid-alone
  if (!with_solid) {
    /**
     * propBuf is a USE of surf_window mesh without ghosts restricted to burning
     * Then it will USE propBufAll for pconn and data
     */
    COM_new_window(propBuf);
    COM_use_dataitem(propBuf, surf_window + ".mesh", 0, bcflag, 1);
    COM_use_dataitem(propBuf, propBufAll + ".pconn");
    COM_use_dataitem(propBuf, propBufAll + ".data");
    create_registered_window_dataitems(propBuf);
    COM_window_init_done(propBuf);
  }

  /*
  // update
  split_surface_window(surf_all, surf_i, surf_nb, surf_b, surf_ni);

  COM_delete_window(surf_nb);
  COM_delete_window(surf_b);
  COM_delete_window(surf_ni);
  */

  // surf_i is the fluidBuf in old rocman
  // Precondition: The window fluidSurf should have been defined and contain
  // the following variables: pf, qc, qr, rhof_alp, nf_alp, tf, Tf, mdot_alp,
  // Tflm_alp, du_alp, and rhofvf_alp.
  /**
   * surf_i is a USE of surf_window mesh restricted to interacting
   * NEW data:
   *   .rhos, .mdot, .mdot_old, .mdot_grad, .mdot_tmp, .vs, .vs_alp, .vs_old
   */
  COM_new_window(surf_i);
  COM_use_dataitem(surf_i, surf_window + ".mesh", 1, bcflag, 0);
  COM_use_dataitem(surf_i, surf_window + ".mesh", 1, bcflag, 1);

  COM_new_dataitem(surf_i + ".rhos", 'e', COM_DOUBLE, 1, "kg/(m^3)");
  COM_new_dataitem(surf_i + ".mdot", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  COM_new_dataitem(surf_i + ".mdot_old", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  COM_new_dataitem(surf_i + ".mdot_grad", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  COM_new_dataitem(surf_i + ".mdot_tmp", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");

  COM_new_dataitem(surf_i + ".vs", 'e', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem(surf_i + ".vs_alp", 'e', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem(surf_i + ".vs_old", 'e', COM_DOUBLE, 3, "m/s");

  COM_resize_array(surf_i + ".data");

  // for compute distances
  if (with_solid) {
    COM_clone_dataitem(surf_i + ".nc_tmp", surf_window + ".nc");
    COM_new_dataitem(surf_i + ".sq_dist", 'n', COM_DOUBLE, 1, "m");
    COM_resize_array(surf_i + ".sq_dist");
  }

  create_registered_window_dataitems(surf_i);
  COM_window_init_done(surf_i);

  /**
   * fluidBuffNG is a USE of surf_i mesh without ghosts
   * It USE data from surf_window, surf_i, and propBufAll
   */
  COM_new_window(fluidBufNG);
  COM_use_dataitem(fluidBufNG, surf_i + ".mesh", 0);
  COM_use_dataitem(fluidBufNG, surf_window + ".data");
  COM_use_dataitem(fluidBufNG, surf_i + ".data");
  COM_use_dataitem(fluidBufNG, propBufAll + ".data");
  create_registered_window_dataitems(fluidBufNG);
  COM_window_init_done(fluidBufNG);

  /**
   * fluidBufB is a USE of surf_window mesh without ghosts restricted to burning
   * It USE data from surf_window, surf_i, and propBufAll
   */
  COM_new_window(fluidBufB);
  COM_use_dataitem(fluidBufB, surf_window + ".mesh", 0, bcflag, 1);
  COM_use_dataitem(fluidBufB, surf_window + ".data");
  COM_use_dataitem(fluidBufB, surf_i + ".data");
  COM_use_dataitem(fluidBufB, propBufAll + ".data");
  create_registered_window_dataitems(fluidBufB);
  COM_window_init_done(fluidBufB);

  /**
   * fluidBufNB is a USE of surf_window mesh without ghosts restricted to
   * non-burning
   * It USE data from surf_window, surf_i, and propBufAll
   */
  COM_new_window(fluidBufNB);
  COM_use_dataitem(fluidBufNB, surf_window + ".mesh", 0, bcflag, 0);
  COM_use_dataitem(fluidBufNB, surf_window + ".data");
  COM_use_dataitem(fluidBufNB, surf_i + ".data");
  COM_use_dataitem(fluidBufNB, propBufAll + ".data");
  create_registered_window_dataitems(fluidBufNB);
  COM_window_init_done(fluidBufNB);

  // Create windows for output of burning patches, including ghost nodes/cells
  // (variables include all of fluidSurf, and mdot, mdot_old, vm, vs, vs_old,
  //  ts, nc_t0, sq_dist)
  /**
   * surf_b is a USE of surf_window mesh and data restricted to burning
   * It USE pconn and vm from propBufAll
   */
  COM_new_window(surf_b);
  COM_use_dataitem(surf_b, surf_window + ".mesh", 1, bcflag, 1);
  COM_use_dataitem(surf_b, surf_window + ".data");

  COM_use_dataitem(surf_b, propBufAll + ".pconn");
  COM_use_dataitem(surf_b, propBufAll + ".vm");

  if (with_solid) {
    COM_use_dataitem(surf_b, surf_i + ".vs");
    COM_use_dataitem(surf_b, surf_i + ".vs_old");
    if (COM_get_dataitem_handle(surf_i + ".ts") > 0) {
      COM_use_dataitem(surf_b, surf_i + ".ts");
      COM_use_dataitem(surf_b, surf_i + ".nc_t0");
      COM_use_dataitem(surf_b, surf_i + ".sq_dist");
    }
  }

  COM_use_dataitem(surf_b, surf_i + ".mdot");
  COM_use_dataitem(surf_b, surf_i + ".mdot_old");

  COM_window_init_done(surf_b);

  // Create windows for output of nonburning patches, including ghosts.
  // (variables include all of fluidSurf, and vm, vm_old, vs, vs_old, ts, nc_t0,
  //  sq_dist)
  /**
   * surf_nb is a USE of surf_window mesh and data restricted to non-burning
   * It USE pconn and vm from propBufAll
   */
  COM_new_window(surf_nb);
  COM_use_dataitem(surf_nb, surf_window + ".mesh", 1, bcflag, 0);
  COM_use_dataitem(surf_nb, surf_window + ".data");

  COM_use_dataitem(surf_nb, propBufAll + ".pconn");
  COM_use_dataitem(surf_nb, propBufAll + ".vm");

  if (with_solid) {
    COM_use_dataitem(surf_nb, surf_i + ".vs");
    COM_use_dataitem(surf_nb, surf_i + ".vs_old");
    if (COM_get_dataitem_handle(surf_i + ".ts") > 0) {
      COM_use_dataitem(surf_nb, surf_i + ".ts");
      COM_use_dataitem(surf_nb, surf_i + ".nc_t0");
      COM_use_dataitem(surf_nb, surf_i + ".sq_dist");
    }
  }
  COM_window_init_done(surf_nb);

  // Create windows for output of nonburning patches, including ghosts.
  // (variables include all of fluidSurf)
  /**
   * surf_ni is a USE of surf_window mesh and data restricted to
   * non-interacting
   * It USE pconn and vm from propBufAll
   */
  COM_new_window(surf_ni);
  COM_use_dataitem(surf_ni, surf_window + ".mesh", 1, bcflag, 2);
  COM_use_dataitem(surf_ni, surf_window + ".data");

  COM_use_dataitem(surf_ni, propBufAll + ".pconn");
  COM_use_dataitem(surf_ni, propBufAll + ".vm");

  COM_window_init_done(surf_ni);

  // setup for PC iterations
  int maxPredCorr = get_coupling()->get_max_ipc();
  if (maxPredCorr > 1) {
    // Create a window to encapsulate the surface data to be backed up
    COM_new_window(fluidBufPC);
    COM_use_dataitem(fluidBufPC, surf_i + ".mesh");
    COM_use_dataitem(fluidBufPC, surf_i + ".nc");
    COM_window_init_done(fluidBufPC);

    // Create a window to store backed-up surface data
    COM_new_window(fluidBufBak);
    COM_use_dataitem(fluidBufBak, fluidBufPC + ".mesh");
    COM_clone_dataitem(fluidBufBak, fluidBufPC + ".nc");
    COM_window_init_done(fluidBufBak);

    // Create a window for backing up volume data
    COM_new_window(fluidVolBak);
    COM_use_dataitem(fluidVolBak, vol_window + ".mesh");
    COM_clone_dataitem(fluidVolBak, vol_window + ".nc");
    COM_clone_dataitem(fluidVolBak, vol_window + ".data");
    COM_window_init_done(fluidVolBak);

    // Create a window for backing up PLAG data
    if (with_plag) {
      COM_new_window(fluidPlagBak);
      COM_use_dataitem(fluidPlagBak, plag_window + ".mesh");
      COM_clone_dataitem(fluidPlagBak, plag_window + ".nc");
      COM_clone_dataitem(fluidPlagBak, plag_window + ".data");
      COM_window_init_done(fluidPlagBak);
    }

    // Create a window for convergence check
    COM_new_window(fluidBufPRE);
    COM_use_dataitem(fluidBufPRE, fluidBufNG + ".mesh");
    COM_clone_dataitem(fluidBufPRE, fluidBufNG + ".vm");
    COM_clone_dataitem(fluidBufPRE, fluidBufNG + ".vs");
    COM_clone_dataitem(fluidBufPRE, fluidBufNG + ".mdot");
    if (with_solid)
      COM_clone_dataitem(fluidBufPRE, fluidBufNG + ".ts");
    COM_window_init_done(fluidBufPRE);

    // Initialize the DataItem handles to be stored/restored
    pc_hdls.push_back({COM_get_dataitem_handle(fluidBufPC + ".all"),
                       COM_get_dataitem_handle(fluidBufBak + ".all")});

    pc_hdls.push_back({COM_get_dataitem_handle(vol_window + ".all"),
                       COM_get_dataitem_handle(fluidVolBak + ".all")});

    if (with_plag) {
      pc_hdls.push_back({COM_get_dataitem_handle(plag_window + ".all"),
                         COM_get_dataitem_handle(fluidPlagBak + ".all")});
    }

    f_mdot_hdl = COM_get_dataitem_handle(fluidBufNG + ".mdot");
    f_mdot_pre_hdl = COM_get_dataitem_handle(fluidBufPRE + ".mdot");
    f_vm_hdl = COM_get_dataitem_handle(fluidBufNG + ".vm");
    f_vm_pre_hdl = COM_get_dataitem_handle(fluidBufPRE + ".vm");
    if (with_solid) {
      f_ts_hdl = COM_get_dataitem_handle(fluidBufNG + ".ts");
      f_ts_pre_hdl = COM_get_dataitem_handle(fluidBufPRE + ".ts");
    }

    tolerMass = get_rocstar_coupling()->get_control_param()->tolerMass;
    tolerVelo = get_rocstar_coupling()->get_control_param()->tolerVelo;
    if (with_solid)
      tolerTrac = get_rocstar_coupling()->get_control_param()->tolerTrac;
  }

  // for compute_distances
  nc_hdl = COM_get_dataitem_handle(fluidBufNG + ".nc");
  nc_tmp_hdl = COM_get_dataitem_handle(fluidBufNG + ".nc_tmp");
  sq_dist_hdl = COM_get_dataitem_handle(fluidBufNG + ".sq_dist");

  if (!get_coupling()->is_initial_start()) {
    int cnstr_type_hdl = COM_get_dataitem_handle(surf_window + ".cnstr_type");
    int cnstr_type_bak_hdl =
        COM_get_dataitem_handle(propBufAll + ".cnstr_type_bak");
    if (cnstr_type_hdl > 0) {
      // backup ".cnstr_type" of fluid and restore it later
      // this will discard the values in HDF restart files!
      COM_call_function(RocBlas::copy, &cnstr_type_hdl, &cnstr_type_bak_hdl);
    }
    read_restart_data();
    if (cnstr_type_hdl > 0) {
      COM_call_function(RocBlas::copy, &cnstr_type_bak_hdl, &cnstr_type_hdl);
    }
  }
}

void FluidAgent::read_restart_data() {
  // Obtain data for fluidBuf
  int atts_hdl = COM_get_dataitem_handle_const(surf_window_in + ".data");
  int buf_hdl = COM_get_dataitem_handle(surf_i + ".data");
  COM_call_function(obtain_attr_handle, &atts_hdl, &buf_hdl);

  // Obtain data for propBufAll
  int propBufAll_hdl = COM_get_dataitem_handle(propBufAll + ".data");
  COM_call_function(obtain_attr_handle, &atts_hdl, &propBufAll_hdl);

  // Obtain pconn for surface propagation
  int pconn_hdl_const =
      COM_get_dataitem_handle_const(surf_window_in + ".pconn");
  int pconn_hdl = COM_get_dataitem_handle(surf_window_in + ".pconn");
  COM_call_function(obtain_attr_handle, &pconn_hdl_const, &pconn_hdl);
  COM_clone_dataitem(propBufAll, surf_window_in + ".pconn");
}

/*
// Visualization window buffers
static const char *ifluid_vis = "ifluid_vis";
static const char *fluid_vis = "fluid_vis";
static const char *fluid_plag_vis = "fluid_plag_vis";

static const char *ifluid_b_vis = "ifluid_b_vis";
static const char *ifluid_nb_vis = "ifluid_nb_vis";
static const char *ifluid_ni_vis = "ifluid_ni_vis";
*/

void FluidAgent::output_visualization_files(double t) {
  // TODO: Define visualization sub-windows
}

void FluidAgent::finalize_windows() {
  if (get_coupling()->get_max_ipc() > 1) {
    COM_delete_window(fluidBufPC);
    COM_delete_window(fluidBufBak);
    COM_delete_window(fluidVolBak);
    if (with_plag)
      COM_delete_window(fluidPlagBak);
    COM_delete_window(fluidBufPRE);
  }

  COM_delete_window(propBufAll);
  COM_delete_window(fluidBufB);
  COM_delete_window(fluidBufNG);

  if (!with_solid)
    COM_delete_window(propBuf);
}

void FluidAgent::init_convergence(int iPredCorr) {
  MAN_DEBUG(2, ("Rocstar: FluidAgent::init_convergence at %d.\n", iPredCorr));

  // Copy current solution to pre for convergence check
  COM_call_function(RocBlas::copy, &f_mdot_hdl, &f_mdot_pre_hdl);
  COM_call_function(RocBlas::copy, &f_vm_hdl, &f_vm_pre_hdl);
  if (with_solid)
    COM_call_function(RocBlas::copy, &f_ts_hdl, &f_ts_pre_hdl);

  // STORE_SOLUTIONS
  RocstarAgent::init_convergence(iPredCorr);
}

bool FluidAgent::check_convergence() const {
  MAN_DEBUG(3, ("Rocstar: FluidAgent::check_convergence .\n"));

  bool b1, b2, b3 = true;
  b1 = check_convergence_helper(f_vm_hdl, f_vm_pre_hdl, tolerVelo, "vm");
  b2 = check_convergence_helper(f_mdot_hdl, f_mdot_pre_hdl, tolerMass, "mdot");
  if (with_solid)
    b3 = check_convergence_helper(f_ts_hdl, f_ts_pre_hdl, tolerTrac, "ts");

  return b1 && b2 && b3;
}

int FluidAgent::compute_integrals() {
  if (compute_integrals_handle > 0) {
    MAN_DEBUG(3, ("Rocstar: FluidAgent::compute_integrals.\n"));
    COM_call_function(compute_integrals_handle, integrals);
  }
  return 1;
}
