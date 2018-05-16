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

#include "sys/types.h"
#include "sys/stat.h"

#include "rocman.h"
#include "BurnAgent.h"
#include "Interpolate.h"

#ifdef STATIC_LINK
extern "C" void COM_F_FUNC2(rocburn_load_module, ROCBURN_LOAD_MODULE)( const char *, int);
extern "C" void COM_F_FUNC2(rocburn_unload_module, ROCBURN_UNLOAD_MODULE)( const char *, int);
#endif

const char *BurnAgent::window_name = "BurnAgent";

BurnAgent::BurnAgent(Coupling *coup, std::string mod, std::string obj, MPI_Comm com, const std::string parentwin) :
  Agent(coup, mod, obj, "Burn-Agent", com, false), parentWin(parentwin)
{
  load_module();

  ignmodel = false;
  std::string ModName(coup->get_control_param()->burn_module);
  if(ModName == "RocburnPY")
    ignmodel = true;

  iburn_all = "iburn_all";   //  for registering dataitems

  iburn_ng = "iburn_ng";    //    burnBuf

  // Material names in files.
  iburn = "iburn";		//  for input
  burn = "burn";

  // SimIN windows
  burnSurfIN = "BurnSurfIN";
  burnVolIN  = "BurnVolIN";

  burnIntBak = "BurnIntBak";

  burnBufOUT = "BurnBufOUT";

  tmp_window = iburn_all;
}

void BurnAgent::load_module()
{
  MAN_DEBUG(3, ("[%d] Rocstar: BurnAgent::load_module %s %s.\n", comm_rank, rocmod_name.c_str(), mod_instance.c_str()));
#ifdef STATIC_LINK
  COM_assertion_msg(rocmod_name == "Rocburn", "Unknown BurnAgent mod!");
  COM_F_FUNC2( rocburn_load_module, ROCBURN_LOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#else
  COM_load_module(rocmod_name.c_str(), mod_instance.c_str());
#endif

  init_function_handles();	// defined in Agent
}

void BurnAgent::unload_module()
{
  MAN_DEBUG(3, ("Rocstar: BurnAgent::unload_module %s %s.\n", rocmod_name.c_str(), mod_instance.c_str()));
#ifdef STATIC_LINK
  COM_assertion_msg(rocmod_name == "Rocburn", "Unknown BurnAgent mod!");
  COM_F_FUNC2( rocburn_unload_module, ROCBURN_UNLOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#else
  if (get_coupling()->in_restart())
    COM_close_module(rocmod_name.c_str(), mod_instance.c_str());
  else
    COM_unload_module(rocmod_name.c_str(), mod_instance.c_str());
#endif
}

void BurnAgent::input( double t) {
  int status;
  status = read_by_control_file( t, iburn, burnSurfIN);
  if (status == -1) {
    COM_new_window(burnSurfIN);
    COM_window_init_done(burnSurfIN);
  }

  status = read_by_control_file( t, burn, burnVolIN);
  if (status == -1) {
    COM_new_window(burnVolIN);
    COM_window_init_done(burnVolIN);
  }
}

void BurnAgent::init_module( double t, double dt) {
  MAN_DEBUG(3, ("BurnAgent::init_module called t:%f.\n", t));

  Agent::init_module(t, dt);

  int initial_start = get_coupling()->initial_start();

  if (initial_start || 
           COM_get_dataitem_handle_const( burnSurfIN+".bflag") <= 0) 
  {
    COM_use_dataitem( burnSurfIN+".mesh", parentWin+".mesh", 0);
    COM_use_dataitem( burnSurfIN+".bflag", parentWin+".bflag", 0);
//    COM_window_init_done( burnSurfIN);
  }

  COM_window_init_done( burnSurfIN);

  if (initial_start)
  {
    // Change burnVolIN to use parent windows' mesh without ghost
    COM_use_dataitem( burnVolIN+".mesh", parentWin+".mesh", 0);
    COM_window_init_done( burnVolIN);
  }

  // Call initialization routine of physics module
  // rocman.f90:INITIALIZE()
  // call to INIT_WRAPPER() in rocburn_2D.f90
  COM_call_function( init_handle, &t, &communicator, &ic_handle, 
		     burnSurfIN.c_str(), burnVolIN.c_str(), 
                     &obtain_attr_handle);

  tbl_flag = *(int *)option_data;
  MAN_DEBUG(3, ("BurnAgent: tbl_flag = %d.\n", tbl_flag));

  // Delete input buffer windows
  COM_delete_window( burnSurfIN);
  COM_delete_window( burnVolIN);
}

void BurnAgent::create_buffer_all()
{
  Agent::create_buffer_all();

  // Precondition: The window burnWin should have been defined and contain
  // the following variables: pf_alp, qc_alp, qr_alp, rhos_alp, Tf_alp,
  // rb, and Tflm. It should have uses the dataitems rb,
  // pf, qc, qr, rhos_alp, Tf, and Tflm_alp from parent window.

  int with_qc = COM_get_dataitem_handle_const( surf_window+".qc_alp") > 0;
  int with_qr = COM_get_dataitem_handle_const( surf_window+".qr_alp") > 0;
  int with_Tf = COM_get_dataitem_handle_const( surf_window+".Tf_alp") > 0;
  int with_Tv = COM_get_dataitem_handle_const( surf_window+".Tv_alp") > 0;
  int with_dn = COM_get_dataitem_handle_const( surf_window+".dn_alp") > 0;
  int with_rhos = COM_get_dataitem_handle_const( surf_window+".rhos_alp") > 0;

  MAN_DEBUG(3, ("with_qc=%d with_qr=%d with_Tf=%d with_rhos=%d\n with_Tv=%d\n with_dn=%d\n", 
		with_qc, with_qr, with_Tf,with_rhos,with_Tv,with_dn));

  COM_use_dataitem( iburn_all, parentWin+".Tflm_alp");
  COM_clone_dataitem( iburn_all+".Tflm_old", surf_window+".Tflm");

  COM_use_dataitem( iburn_all, parentWin+".pf");
  COM_use_dataitem( iburn_all, parentWin+".rhos");
  if (with_qc) COM_use_dataitem( iburn_all, parentWin+".qc");
  if (with_qr) COM_use_dataitem( iburn_all, parentWin+".qr");
  if (with_Tf) COM_use_dataitem( iburn_all, parentWin+".Tf");
  if (with_Tv) COM_use_dataitem( iburn_all, parentWin+".Tv");
  if (with_dn) COM_use_dataitem( iburn_all, parentWin+".dn");

  COM_clone_dataitem( iburn_all+".rb_alp", surf_window+".rb");
  COM_clone_dataitem( iburn_all+".rb_old", surf_window+".rb");
  COM_clone_dataitem( iburn_all+".rb_grad", surf_window+".rb");

    // Create back-up data fields if predictor-corrector iteration is on
  COM_clone_dataitem( iburn_all+".pf_old", surf_window+".pf_alp");
  if ( with_qc ) {
       COM_clone_dataitem( iburn_all+".qc_old", surf_window+".qc_alp");
       COM_clone_dataitem( iburn_all+".qc_grad", surf_window+".qc_alp");
  }
  if ( with_qr) {
       COM_clone_dataitem( iburn_all+".qr_old", surf_window+".qr_alp");
       COM_clone_dataitem( iburn_all+".qr_grad", surf_window+".qr_alp");
  }
  if ( with_rhos) 
       COM_clone_dataitem( iburn_all+".rhos_old", surf_window+".rhos_alp");
  if (with_Tf) 
       COM_clone_dataitem( iburn_all+".Tf_old", surf_window+".Tf_alp");
  if (with_Tv) 
       COM_clone_dataitem( iburn_all+".Tv_old", surf_window+".Tv_alp");
  if (with_dn) 
       COM_clone_dataitem( iburn_all+".dn_old", surf_window+".dn_alp");

  create_registered_window_dataitems( iburn_all);
  COM_window_init_done( iburn_all);

   // no ghost
  COM_new_window( iburn_ng);
  COM_use_dataitem( iburn_ng+".all", surf_window+".all", 0);
  COM_use_dataitem( iburn_ng, iburn_all+".data");
  create_registered_window_dataitems( iburn_ng);
  COM_window_init_done( iburn_ng);

  // Create windows for writing interface data
  COM_new_window( burnBufOUT);
  COM_use_dataitem( burnBufOUT, surf_window+".all");

  COM_use_dataitem( burnBufOUT, iburn_all+".rb_old");
  COM_use_dataitem( burnBufOUT, iburn_all+".pf_old");
  COM_use_dataitem( burnBufOUT, iburn_all+".Tflm_old");
  if ( with_qc)
         COM_use_dataitem( burnBufOUT, iburn_all+".qc_old");
  if ( with_qr) 
         COM_use_dataitem( burnBufOUT, iburn_all+".qr_old");
  if (with_Tf) 
         COM_use_dataitem( burnBufOUT, iburn_all+".Tf_old");
  if (with_Tv) 
         COM_use_dataitem( burnBufOUT, iburn_all+".Tv_old");
  if (with_dn) 
         COM_use_dataitem( burnBufOUT, iburn_all+".dn_old");
  COM_window_init_done( burnBufOUT);

    // Create a window for backing up the internal data of Rocburn
    //    if predictor-corrector iteration is on.
  int maxPredCorr = get_coupling()->get_max_ipc();
  if ( maxPredCorr>1) {
       // Create window for backing up Rocburn's internal data
     COM_new_window( burnIntBak);
     COM_use_dataitem( burnIntBak+".mesh", vol_window+".mesh");
     COM_clone_dataitem( burnIntBak, vol_window+".data");

     COM_window_init_done( burnIntBak);

     pc_hdls[0][0] = COM_get_dataitem_handle( vol_window+".data");
     pc_hdls[1][0] = COM_get_dataitem_handle( burnIntBak+".data");
     pc_count = 1;
  }

  if (!get_coupling()->initial_start()) read_restart_data();
}

void BurnAgent::read_restart_data()
{
  int atts_hdl = COM_get_dataitem_handle_const(burnSurfIN+".data");
  int buf_hdl = COM_get_dataitem_handle( iburn_all+".data");
  COM_call_function( obtain_attr_handle, &atts_hdl, &buf_hdl);
}

void BurnAgent::output_restart_files( double t) {
  static const std::string iburn_prefix = "iburn_all";

  // Write out surface sub-windows
  //write_data_files( t, iburn_prefix, iburn_all+".all");
  write_data_files( t, iburn_prefix, burnBufOUT+".all");
  write_control_file( t, "iburn*", surf_window.c_str());

  // Write out volume window
  write_data_files( t, burn, (vol_window+".all").c_str());
  write_control_file( t, burn, vol_window.c_str());
}

// Visualization window buffers
static const char *iburn_vis = "iburn_vis";
static const char *burn_vis = "burn_vis";
static const char *burn_plag_vis = "burn_plag_vis";

static const char *iburn_b_vis = "iburn_b_vis";
static const char *iburn_nb_vis = "iburn_nb_vis";
static const char *iburn_ni_vis = "iburn_ni_vis";


void BurnAgent::output_visualization_files( double t) {
  // TODO: Define visualization sub-windows
}

void BurnAgent::finalize()
{
  COM_delete_window( iburn_ng);
  COM_delete_window( iburn_all);
  int maxPredCorr = get_coupling()->get_max_ipc();
  if ( maxPredCorr>1) 
       COM_delete_window( burnIntBak);

  Agent::finalize();
}








