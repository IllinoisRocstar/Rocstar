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
// $Id: SolidAgent.C,v 1.53 2009/05/12 20:17:48 mtcampbe Exp $

#include <utility>

#include "SolidAgent.h"

#include "RocstarCoupling.h"
#include "rocman.h"

#ifdef STATIC_LINK
extern "C" void COM_F_FUNC2(rocsolid_load_module,
                            ROCSOLID_LOAD_MODULE)(const char *, int);
extern "C" void COM_F_FUNC2(rocsolid_unload_module,
                            ROCSOLID_UNLOAD_MODULE)(const char *, int);
#ifdef ROCFRAC3
extern "C" void COM_F_FUNC2(rocfrac3_load_module,
                            ROCFRAC3_LOAD_MODULE)(const char *, int);
extern "C" void COM_F_FUNC2(rocfrac3_unload_module,
                            ROCFRAC3_UNLOAD_MODULE)(const char *, int);
#else
extern "C" void COM_F_FUNC2(rocfrac_load_module,
                            ROCFRAC_LOAD_MODULE)(const char *, int);
extern "C" void COM_F_FUNC2(rocfrac_unload_module,
                            ROCFRAC_UNLOAD_MODULE)(const char *, int);
#endif
#endif

// const char *SolidAgent::window_name = "SolidAgent";

SolidAgent::SolidAgent(RocstarCoupling *coup, std::string mod, std::string obj,
                       MPI_Comm com, bool withFluid)
    : RocstarAgent(coup, std::move(mod), std::move(obj), "SolidAgent", "solid",
                   "isolid", com, false),
      with_fluid(withFluid) {
  with_ALE = false;

  solidBufBase = "SolidBufBase";
  solidBuf = "SolidBuf";       // for registering dataitems
  propBufAll = "SolidPropAll"; // for surface propagation (containing all panes)
  propBuf = "SolidProp";

  solidBufBak = "SolidBufBak";
  solidVolBak = "SolidVolBak";
}

/* AEG: Functionality moved to parent class. The STATIC_LINK flag is now in
 * IMPACT. Since load_module() is called during constructing of the parent, it
 * cannot be virtual and defined in the derived. *//*
void SolidAgent::load_module() {
#ifdef STATIC_LINK
  MAN_DEBUG(3, ("[%d] Rocstar: SolidAgent::load_module %s %s.\n", comm_rank,
                module_lname.c_str(), module_wname.c_str()));
  if (module_lname == "Rocsolid")
    COM_F_FUNC2(rocsolid_load_module, ROCSOLID_LOAD_MODULE)
  (module_wname.c_str(), module_wname.length());
#ifdef ROCFRAC3
  else if (module_lname == "Rocfrac3")
      COM_F_FUNC2(rocfrac3_load_module, ROCFRAC3_LOAD_MODULE)(
          module_wname.c_str(), module_wname.length());
#else
  else if (module_lname == "Rocfrac")
      COM_F_FUNC2(rocfrac_load_module, ROCFRAC_LOAD_MODULE)(
          module_wname.c_str(), module_wname.length());
#endif
  else COM_assertion_msg(0, "Unknown SolidAgent mod!");
#endif
}

void SolidAgent::unload_module() {
#ifdef STATIC_LINK
  MAN_DEBUG(3, ("[%d] Rocstar: SolidAgent::unload_module %s.\n", comm_rank,
                module_lname.c_str()));
  if (module_lname == "Rocsolid")
    COM_F_FUNC2(rocsolid_unload_module, ROCSOLID_UNLOAD_MODULE)
    (module_wname.c_str(), module_wname.length());
#ifdef ROCFRAC3
  else if (module_lname == "Rocfrac3")
    COM_F_FUNC2(rocfrac3_unload_module, ROCFRAC3_UNLOAD_MODULE)
    (module_wname.c_str(), module_wname.length());
#else
  else if (module_lname == "Rocfrac")
    COM_F_FUNC2(rocfrac_unload_module, ROCFRAC_UNLOAD_MODULE)
    (module_wname.c_str(), module_wname.length());
#endif
#endif
}
*/

// called from callback
void SolidAgent::create_buffer_all() {
  int dummy = COM_get_dataitem_handle_const(surf_window + ".vbar_alp");
  with_ALE = dummy > 0;
  if (comm_rank == 0)
    std::cout << "Rocstar: *** Solid with ALE is " << with_ALE << std::endl;

  if (with_ALE) {
    COM_new_window(propBufAll);
    COM_use_dataitem(propBufAll + ".mesh", surf_window + ".mesh");
    COM_use_dataitem(propBufAll + ".vbar_alp", surf_window + ".vbar_alp");
    COM_clone_dataitem(propBufAll + ".vbar", surf_window + ".vbar_alp");
    COM_clone_dataitem(propBufAll + ".vbar_old", surf_window + ".vbar_alp");
    COM_clone_dataitem(propBufAll + ".vbar_grad", surf_window + ".vbar_alp");
    COM_new_dataitem(propBufAll + ".rb", 'e', COM_DOUBLE, 1, "m/s");
    COM_new_dataitem(propBufAll + ".positions", 'n', COM_DOUBLE, 3, "m");
    COM_new_dataitem(propBufAll + ".constrained", 'n', COM_INTEGER, 1, "");
    COM_resize_array(propBufAll + ".positions");
    COM_resize_array(propBufAll + ".constrained");
    COM_resize_array(propBufAll + ".rb");
    create_registered_window_dataitems(propBufAll);
    COM_window_init_done(propBufAll);

    COM_new_window(propBuf);
    std::string bc_str = surf_window + ".bcflag";
    COM_use_dataitem(propBuf + ".mesh", surf_window + ".mesh", 1, bc_str, 1);
    COM_use_dataitem(propBuf + ".pconn", propBufAll + ".pconn");
    if (COM_get_dataitem_handle(surf_window + ".cnstr_type") > 0) {
      COM_use_dataitem(propBuf, surf_window + ".cnstr_type");
    }
    COM_use_dataitem(propBuf + ".vbar_alp", propBufAll + ".vbar_alp");
    COM_use_dataitem(propBuf + ".vbar", propBufAll + ".vbar");
    COM_use_dataitem(propBuf + ".vbar_old", propBufAll + ".vbar_old");
    COM_use_dataitem(propBuf + ".vbar_grad", propBufAll + ".vbar_grad");
    COM_use_dataitem(propBuf + ".rb", propBufAll + ".rb");
    create_registered_window_dataitems(propBuf);
    COM_window_init_done(propBuf);
  }

  COM_new_window(solidBufBase);
  std::string bcflag = surf_window + ".bcflag";
  COM_use_dataitem(solidBufBase + ".mesh", surf_window + ".mesh", 1, bcflag, 0);
  COM_use_dataitem(solidBufBase + ".mesh", surf_window + ".mesh", 1, bcflag, 1);
  if (with_ALE) {
    COM_new_dataitem(solidBufBase + ".areas", 'e', COM_DOUBLE, 1, "m^2");
    COM_resize_array(solidBufBase + ".areas");
    COM_new_dataitem(solidBufBase + ".mdot", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
    COM_resize_array(solidBufBase + ".mdot");
  }
  COM_clone_dataitem(solidBufBase + ".ts", surf_window + ".ts_alp");
  COM_clone_dataitem(solidBufBase + ".ts_old", surf_window + ".ts_alp");
  COM_clone_dataitem(solidBufBase + ".ts_grad", surf_window + ".ts_alp");
  create_registered_window_dataitems(solidBufBase);
  COM_window_init_done(solidBufBase);

  COM_new_window(solidBuf);
  COM_use_dataitem(solidBuf, solidBufBase + ".all");
  COM_use_dataitem(solidBuf, surf_window + ".data");
  if (with_ALE)
    COM_use_dataitem(solidBuf, propBufAll + ".data");
  if (with_fluid) { // coupled
    COM_use_dataitem(solidBuf + ".x", solidBufBase + ".nc");
    COM_clone_dataitem(solidBuf + ".nc", surf_window + ".nc");
  }
  create_registered_window_dataitems(solidBuf);
  COM_window_init_done(solidBuf);

  /*
  split_surface_window(surf_all, surf_i, surf_nb, surf_b, surf_ni);

  COM_delete_window(surf_nb);
  COM_delete_window(surf_b);
  COM_delete_window(surf_ni);
  */

  // Create a window for output solid non-burning patch data
  COM_new_window(surf_nb);
  COM_use_dataitem(surf_nb + ".mesh", surf_window + ".mesh", 1, bcflag, 0);
  COM_use_dataitem(surf_nb, surf_window + ".data");
  if (with_ALE) {
    COM_use_dataitem(surf_nb + ".pconn", propBufAll + ".pconn");
    COM_use_dataitem(surf_nb, propBufAll + ".vbar");
    COM_use_dataitem(surf_nb, propBufAll + ".vbar_old");
    COM_use_dataitem(surf_nb, propBufAll + ".rb");
  }
  COM_use_dataitem(surf_nb, solidBufBase + ".ts");
  COM_use_dataitem(surf_nb, solidBufBase + ".ts_old");
  COM_window_init_done(surf_nb);

  // Create a window for output solid interface data
  COM_new_window(surf_b);
  COM_use_dataitem(surf_b + ".mesh", surf_window + ".mesh", 1, bcflag, 1);
  COM_use_dataitem(surf_b, surf_window + ".data");
  if (with_ALE) {
    COM_use_dataitem(surf_b + ".pconn", propBufAll + ".pconn");
    COM_use_dataitem(surf_b, propBufAll + ".vbar");
    COM_use_dataitem(surf_b, propBufAll + ".vbar_old");
    COM_use_dataitem(surf_b, propBufAll + ".rb");
  }
  COM_use_dataitem(surf_b, solidBufBase + ".ts");
  COM_use_dataitem(surf_b, solidBufBase + ".ts_old");
  COM_window_init_done(surf_b);

  // Create a window for non-solid/fluid interface
  COM_new_window(surf_ni);
  COM_use_dataitem(surf_ni + ".mesh", surf_window + ".mesh", 1, bcflag, 2);
  COM_use_dataitem(surf_ni, surf_window + ".data");
  if (with_ALE) {
    COM_use_dataitem(surf_ni + ".pconn", propBufAll + ".pconn");
    COM_use_dataitem(surf_ni, propBufAll + ".vbar");
    COM_use_dataitem(surf_ni, propBufAll + ".vbar_old");
    COM_use_dataitem(surf_ni, propBufAll + ".rb");
  }
  COM_window_init_done(surf_ni);

  // setup for pc iterations
  if (get_coupling()->get_max_ipc() > 1) {
    // Create window for backing surface data
    COM_new_window(solidBufBak);
    COM_use_dataitem(solidBufBak + ".mesh", solidBuf + ".mesh");
    COM_clone_dataitem(solidBufBak + ".x", solidBuf + ".x");
    COM_clone_dataitem(solidBufBak + ".y", solidBuf + ".nc");
    COM_window_init_done(solidBufBak);

    // Create window for backing up volume data
    COM_new_window(solidVolBak);
    COM_use_dataitem(solidVolBak + ".mesh", vol_window + ".mesh");
    COM_clone_dataitem(solidVolBak, vol_window + ".nc");
    COM_clone_dataitem(solidVolBak, vol_window + ".data");
    COM_window_init_done(solidVolBak);

    // Initlaize the dataitem handles to be stored/restored
    pc_hdls.push_back({COM_get_dataitem_handle(vol_window + ".all"),
                       COM_get_dataitem_handle(solidVolBak + ".all")});
    pc_hdls.push_back({COM_get_dataitem_handle(solidBuf + ".x"),
                       COM_get_dataitem_handle(solidBufBak + ".x")});
    pc_hdls.push_back({COM_get_dataitem_handle(solidBuf + ".nc"),
                       COM_get_dataitem_handle(solidBufBak + ".y")});
  }

  // setup variables
  traction_mode =
      get_rocstar_coupling()->get_rocmancontrol_param()->traction_mode;

  std::string unit;
  char loc;
  COM_get_dataitem(surf_all + ".rhos", &loc, &dummy, &dummy, &unit);
  if (loc == 'w' || loc == 'p')
    rhos_mode = 1;
  else if (loc == 'n' || loc == 'e') {
    rhos_mode = 2;
  } else {
    COM_abort_msg(EXIT_FAILURE,
                  "Rocstar: Error: Unknown type of location of rhos: " +
                      std::to_string(loc));
  }

  COM_get_dataitem(solidBufBase + ".ts", &loc, &dummy, &size_ts, &unit);
  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_abort_msg(EXIT_FAILURE, "If traction mode is with sheer, then solid "
                                "tractions must be vectors!");
  }

  // for compute_distances
  y_hdl = COM_get_dataitem_handle(solidBuf + ".nc");

  if (!get_coupling()->is_initial_start())
    read_restart_data();
}

void SolidAgent::read_restart_data() {
  int atts_hdl = COM_get_dataitem_handle_const(surf_window_in + ".data");
  int buf_hdl = COM_get_dataitem_handle(solidBufBase + ".data");
  COM_call_function(obtain_attr_handle, &atts_hdl, &buf_hdl);

  if (with_ALE) {
    // Obtain data for surface propagation if ALE is enabled
    int propBufAll_hdl = COM_get_dataitem_handle(propBufAll + ".data");
    COM_call_function(obtain_attr_handle, &atts_hdl, &propBufAll_hdl);

    // Obtain pconn for surface propagation
    int pconn_hdl_const =
        COM_get_dataitem_handle_const(surf_window_in + ".pconn");
    int pconn_hdl = COM_get_dataitem_handle(surf_window_in + ".pconn");
    COM_call_function(obtain_attr_handle, &pconn_hdl_const, &pconn_hdl);
    COM_clone_dataitem(propBufAll + ".pconn", surf_window_in + ".pconn");
  }
}

/*
// Visualization window buffers
static const char *isolid_vis = "isolid_vis";
static const char *solid_vis = "solid_vis";
static const char *solid_plag_vis = "solid_plag_vis";

static const char *isolid_b_vis = "isolid_b_vis";
static const char *isolid_nb_vis = "isolid_nb_vis";
static const char *isolid_ni_vis = "isolid_ni_vis";
*/

void SolidAgent::output_visualization_files(double t) {
  // TODO: Define visualization sub-windows
}

void SolidAgent::finalize_windows() {
  if (get_coupling()->get_max_ipc() > 1) {
    COM_delete_window(solidBufBak);
    COM_delete_window(solidVolBak);
  }

  COM_delete_window(solidBufBase);
  COM_delete_window(solidBuf);

  if (with_ALE) {
    COM_delete_window(propBufAll);
    COM_delete_window(propBuf);
  }
}

int SolidAgent::compute_integrals() {
  if (compute_integrals_handle > 0) {
    MAN_DEBUG(3, ("Rocstar: SolidAgent::compute_integrals.\n"));
    COM_call_function(compute_integrals_handle, integrals);
  }
  return 1;
}
