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
// $Id: BurnAgent.C,v 1.36 2010/02/18 21:47:40 juzhang Exp $

#include "BurnAgent.h"

#include <utility>

#include "Control_parameters.h"
#include "RocstarCoupling.h"
#include "rocman.h"

#ifdef STATIC_LINK
extern "C" void COM_F_FUNC2(rocburn_load_module,
                            ROCBURN_LOAD_MODULE)(const char *, int);
extern "C" void COM_F_FUNC2(rocburn_unload_module,
                            ROCBURN_UNLOAD_MODULE)(const char *, int);
#endif

BurnAgent::BurnAgent(RocstarCoupling *coup, std::string mod, std::string obj,
                     MPI_Comm com, std::string parentwin)
    : RocstarAgent(coup, std::move(mod), std::move(obj), "BurnAgent", "burn",
                   "iburn", com, false),
      parentWin(std::move(parentwin)) {
  ignmodel = false;
  std::string ModName(coup->get_control_param()->burn_module);
  if (ModName == "RocburnPY")
    ignmodel = true;

  iburn_ng = "iburn_ng"; // burnBuf

  burnIntBak = "BurnIntBak";
}

/* AEG: Functionality moved to parent class. The STATIC_LINK flag is now in
 * IMPACT. Since load_module() is called during constructing of the parent, it
 * cannot be virtual and defined in the derived. *//*
void BurnAgent::load_module() {
#ifdef STATIC_LINK
  MAN_DEBUG(3, ("[%d] Rocstar: BurnAgent::load_module %s %s.\n", comm_rank,
                module_lname.c_str(), module_wname.c_str()));
  COM_assertion_msg(module_lname == "Rocburn", "Unknown BurnAgent mod!");
  COM_F_FUNC2(rocburn_load_module, ROCBURN_LOAD_MODULE)
  (module_wname.c_str(), module_wname.length());
#endif
}

void BurnAgent::unload_module() {
#ifdef STATIC_LINK
  MAN_DEBUG(3, ("Rocstar: BurnAgent::unload_module %s %s.\n",
                module_lname.c_str(), module_wname.c_str()));
  COM_assertion_msg(module_lname == "Rocburn", "Unknown BurnAgent mod!");
  COM_F_FUNC2(rocburn_unload_module, ROCBURN_UNLOAD_MODULE)
  (module_wname.c_str(), module_wname.length());
#endif
}
*/

/* AEG: Functionality moved to input_fallback_surf() and input_fallback_vol().
void BurnAgent::init_module(double t) {
  if (get_coupling()->is_initial_start() ||
      COM_get_dataitem_handle_const(surf_window_in + ".bflag") <= 0) {
    COM_use_dataitem(surf_window_in + ".mesh", parentWin + ".mesh", 0);
    COM_use_dataitem(surf_window_in + ".bflag", parentWin + ".bflag", 0);
    COM_window_init_done(surf_window_in);
  }

  // COM_window_init_done(surf_window_in);

  if (get_coupling()->is_initial_start()) {
    // Change burnVolIN to use parent windows' mesh without ghost
    COM_use_dataitem(vol_window_in + ".mesh", parentWin + ".mesh", 0);
    COM_window_init_done(vol_window_in);
  }
}
*/

void BurnAgent::parse_ic_options(void *option_data) {
  tbl_flag = *(int *)option_data;
  MAN_DEBUG(3, ("BurnAgent: tbl_flag = %d.\n", tbl_flag));
}

void BurnAgent::input_fallback_surf(const std::string &surface_window_in) {
  COM_use_dataitem(surface_window_in + ".mesh", parentWin + ".mesh", 0);
  COM_use_dataitem(surface_window_in + ".bflag", parentWin + ".bflag", 0);
}
void BurnAgent::input_fallback_vol(const std::string &volume_window_in) {
  COM_use_dataitem(volume_window_in + ".mesh", parentWin + ".mesh", 0);
}

void BurnAgent::create_buffer_all() {
  // Precondition: The window burnWin should have been defined and contain
  // the following variables: pf_alp, qc_alp, qr_alp, rhos_alp, Tf_alp,
  // rb, and Tflm. It should have uses the dataitems rb,
  // pf, qc, qr, rhos_alp, Tf, and Tflm_alp from parent window.

  bool with_qc = COM_get_dataitem_handle_const(surf_window + ".qc_alp") > 0;
  bool with_qr = COM_get_dataitem_handle_const(surf_window + ".qr_alp") > 0;
  bool with_Tf = COM_get_dataitem_handle_const(surf_window + ".Tf_alp") > 0;
  bool with_Tv = COM_get_dataitem_handle_const(surf_window + ".Tv_alp") > 0;
  bool with_dn = COM_get_dataitem_handle_const(surf_window + ".dn_alp") > 0;
  bool with_rhos = COM_get_dataitem_handle_const(surf_window + ".rhos_alp") > 0;

  MAN_DEBUG(
      3,
      ("with_qc=%d with_qr=%d with_Tf=%d with_rhos=%d with_Tv=%d with_dn=%d\n",
       with_qc, with_qr, with_Tf, with_rhos, with_Tv, with_dn));

  /*
   * surf_all is a
   */
  COM_use_dataitem(surf_all, parentWin + ".bcflag");
  COM_use_dataitem(surf_all, parentWin + ".Tflm_alp");
  COM_clone_dataitem(surf_all + ".Tflm_old", surf_window + ".Tflm");
  COM_use_dataitem(surf_all, parentWin + ".pf");
  COM_clone_dataitem(surf_all + ".pf_old", surf_window + ".pf_alp");
  if (with_qc) {
    COM_use_dataitem(surf_all, parentWin + ".qc");
    COM_clone_dataitem(surf_all + ".qc_old", surf_window + ".qc_alp");
    COM_clone_dataitem(surf_all + ".qc_grad", surf_window + ".qc_alp");
  }
  if (with_qr) {
    COM_use_dataitem(surf_all, parentWin + ".qr");
    COM_clone_dataitem(surf_all + ".qr_old", surf_window + ".qr_alp");
    COM_clone_dataitem(surf_all + ".qr_grad", surf_window + ".qr_alp");
  }
  if (with_Tf) {
    COM_use_dataitem(surf_all, parentWin + ".Tf");
    COM_clone_dataitem(surf_all + ".Tf_old", surf_window + ".Tf_alp");
  }
  if (with_Tv) {
    COM_use_dataitem(surf_all, parentWin + ".Tv");
    COM_clone_dataitem(surf_all + ".Tv_old", surf_window + ".Tv_alp");
  }
  if (with_dn) {
    COM_use_dataitem(surf_all, parentWin + ".dn");
    COM_clone_dataitem(surf_all + ".dn_old", surf_window + ".dn_alp");
  }
  if (with_rhos) {
    COM_use_dataitem(surf_all, parentWin + ".rhos");
    COM_clone_dataitem(surf_all + ".rhos_old", surf_window + ".rhos_alp");
  }
  COM_clone_dataitem(surf_all + ".rb_alp", surf_window + ".rb");
  COM_clone_dataitem(surf_all + ".rb_old", surf_window + ".rb");
  COM_clone_dataitem(surf_all + ".rb_grad", surf_window + ".rb");
  create_registered_window_dataitems(surf_all);
  COM_window_init_done(surf_all);

  // no ghost
  /**
   * iburn_ng is a USE of all of surf_window and the data of surf_all
   */
  COM_new_window(iburn_ng);
  COM_use_dataitem(iburn_ng + ".all", surf_window + ".all", 0);
  COM_use_dataitem(iburn_ng, surf_all + ".data");
  create_registered_window_dataitems(iburn_ng);
  COM_window_init_done(iburn_ng);

  // Create a window for backing up the internal data of Rocburn
  // if predictor-corrector iteration is on.
  if (get_coupling()->get_max_ipc() > 1) {
    // Create window for backing up Rocburn's internal data
    COM_new_window(burnIntBak);
    COM_use_dataitem(burnIntBak + ".mesh", vol_window + ".mesh");
    COM_clone_dataitem(burnIntBak, vol_window + ".data");
    COM_window_init_done(burnIntBak);

    pc_hdls.push_back({COM_get_dataitem_handle(vol_window + ".data"),
                       COM_get_dataitem_handle(burnIntBak + ".data")});
  }

  if (!get_coupling()->is_initial_start())
    read_restart_data();
}

void BurnAgent::read_restart_data() {
  int atts_hdl = COM_get_dataitem_handle_const(surf_window_in + ".data");
  int buf_hdl = COM_get_dataitem_handle(surf_all + ".data");
  COM_call_function(obtain_attr_handle, &atts_hdl, &buf_hdl);
}

/*
// Visualization window buffers
static const char *iburn_vis = "iburn_vis";
static const char *burn_vis = "burn_vis";
static const char *burn_plag_vis = "burn_plag_vis";

static const char *iburn_b_vis = "iburn_b_vis";
static const char *iburn_nb_vis = "iburn_nb_vis";
static const char *iburn_ni_vis = "iburn_ni_vis";
*/

void BurnAgent::output_visualization_files(double t) {
  // TODO: Define visualization sub-windows
}

void BurnAgent::finalize_windows() {
  if (get_coupling()->get_max_ipc() > 1)
    COM_delete_window(burnIntBak);

  COM_delete_window(iburn_ng);
}
