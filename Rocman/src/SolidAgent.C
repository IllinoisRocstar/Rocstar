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

#include "sys/types.h"
#include "sys/stat.h"

#include "rocman.h"
#include "SolidAgent.h"
#include "Interpolate.h"

#ifdef STATIC_LINK
extern "C" void COM_F_FUNC2(rocsolid_load_module, ROCSOLID_LOAD_MODULE)( const char *, int);
extern "C" void COM_F_FUNC2(rocsolid_unload_module, ROCSOLID_UNLOAD_MODULE)( const char *, int);
#ifdef ROCFRAC3
extern "C" void COM_F_FUNC2(rocfrac3_load_module, ROCFRAC3_LOAD_MODULE)( const char *, int);
extern "C" void COM_F_FUNC2(rocfrac3_unload_module, ROCFRAC3_UNLOAD_MODULE)( const char *, int);
#else
extern "C" void COM_F_FUNC2(rocfrac_load_module, ROCFRAC_LOAD_MODULE)( const char *, int);
extern "C" void COM_F_FUNC2(rocfrac_unload_module, ROCFRAC_UNLOAD_MODULE)( const char *, int);
#endif
#endif

//const char *SolidAgent::window_name = "SolidAgent";

static std::string agentCount = "0";       // only support up to 10 solid agents

SolidAgent::SolidAgent(Coupling *coup, std::string mod, std::string obj, MPI_Comm com, int withFluid) :
  Agent(coup, mod, obj, "Solid", com, false), with_fluid(withFluid)
{
  load_module();

  withALE = 0;

    // TODO: need to personalize these variables
  if (agentCount.length()==1) agentCount[0] ++;
  isolid_i = "isolid";  isolid_i += agentCount;    // "SolidBuf"
  isolid_all = "isolid_all"; isolid_all += agentCount; //for registering dataitems

  solidBufBase = "SolidBufBase"; solidBufBase += agentCount;  
  solidBuf = "SolidBuf"; solidBuf += agentCount; //  for registering dataitems
  propBufAll = "SolidPropAll"; // for surface propagation (containing all panes)
  propBufAll += agentCount;
  propBuf = "SolidProp"; propBuf += agentCount;

  // Surface windows for Rocout
  isolid_b = "isolid_b";	// buring
  isolid_b += agentCount;
  isolid_nb = "isolid_nb";     // non-buring
  isolid_nb += agentCount;
  isolid_ni = "isolid_ni";     // noninteracting
  isolid_ni += agentCount;

  // Material names in files.
  isolid = "isolid";		//  for input
  solid = "solid";

  // SimIN windows
  solidSurfIN = "SolidSurfIN"; solidSurfIN += agentCount;
  solidVolIN = "SolidVolIN";   solidVolIN += agentCount;

  solidBufBak    = "SolidBufBak"; solidBufBak += agentCount;
  solidVolBak    = "SolidVolBak"; solidVolBak += agentCount;

  tmp_window = isolid_all;
}

void SolidAgent::load_module()
{
  MAN_DEBUG(3, ("[%d] Rocstar: SolidAgent::load_module %s %s.\n", comm_rank, rocmod_name.c_str(), mod_instance.c_str()));

#ifdef STATIC_LINK
  if (rocmod_name == "Rocsolid")
    COM_F_FUNC2( rocsolid_load_module, ROCSOLID_LOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#ifdef ROCFRAC3
  else if (rocmod_name == "Rocfrac3")
    COM_F_FUNC2( rocfrac3_load_module, ROCFRAC3_LOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#else
  else if (rocmod_name == "Rocfrac")
    COM_F_FUNC2( rocfrac_load_module, ROCFRAC_LOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#endif
  else
    COM_assertion_msg(0, "Unknown SolidAgent mod!");
#else
    // dynamic loading
  COM_assertion_msg(rocmod_name == "Rocsolid" || rocmod_name == "Rocfrac" 
                    || rocmod_name == "Rocfrac3", 
               (std::string("Unknown SolidAgent module:")+rocmod_name).c_str());
  COM_load_module(rocmod_name.c_str(), mod_instance.c_str());
#endif

  init_function_handles();      // defined in Agent
}

void SolidAgent::unload_module()
{
  MAN_DEBUG(3, ("[%d] Rocstar: SolidAgent::unload_module %s.\n", comm_rank, rocmod_name.c_str()));
#ifdef STATIC_LINK
  if (rocmod_name == "Rocsolid")
    COM_F_FUNC2( rocsolid_unload_module, ROCSOLID_UNLOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#ifdef ROCFRAC3
  else if (rocmod_name == "Rocfrac3")
    COM_F_FUNC2( rocfrac3_unload_module, ROCFRAC3_UNLOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#else
  else if (rocmod_name == "Rocfrac")
    COM_F_FUNC2( rocfrac_unload_module, ROCFRAC_UNLOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#endif
#else
  if (get_coupling()->in_restart())
    COM_close_module(rocmod_name.c_str(), mod_instance.c_str());
  else
    COM_unload_module(rocmod_name.c_str(), mod_instance.c_str());
#endif
}

void SolidAgent::input( double t) {
  read_by_control_file( t, isolid, solidSurfIN);
  read_by_control_file( t, solid, solidVolIN);
}

void SolidAgent::init_module( double t, double dt) {
  MAN_DEBUG(3, ("Rocstar: SolidAgent::init_module t=%e dt=%e.\n", t, dt));

  Agent::init_module(t, dt);

  // Call initialization routine of physics module
  // rocman.f90:INITIALIZE()
  COM_call_function( init_handle, &t, &communicator, &ic_handle, 
		     solidSurfIN.c_str(), solidVolIN.c_str(), 
                     &obtain_attr_handle);

  // Delete input buffer windows
  COM_delete_window( solidSurfIN);
  COM_delete_window( solidVolIN);

  // Split surface window for output
  // ifluid is fluidBuf in INITIALIZE_FLUID() ???
  //split_surface_window( isolid_all, isolid_i, isolid_nb, isolid_b, isolid_ni);
}

// called from callback
void SolidAgent::create_buffer_all()
{

  int dummy = COM_get_dataitem_handle_const( surf_window+".vbar_alp");
  withALE = (dummy > 0);

  if (comm_rank == 0)
   MAN_DEBUG(2, ("Rocstar: *** Solid with ALE is %d \n", withALE));

  Agent::create_buffer_all();
  split_surface_window( isolid_all, isolid_i, isolid_nb, isolid_b, isolid_ni);

  if (withALE) {
    COM_new_window( propBufAll);
    COM_use_dataitem( propBufAll+".mesh", surf_window+".mesh");
    COM_use_dataitem( propBufAll+".vbar_alp", surf_window+".vbar_alp");
    COM_clone_dataitem( propBufAll+".vbar", surf_window+".vbar_alp");
    COM_clone_dataitem( propBufAll+".vbar_old", surf_window+".vbar_alp");
    COM_clone_dataitem( propBufAll+".vbar_grad", surf_window+".vbar_alp");
    COM_new_dataitem( propBufAll+".rb", 'e', COM_DOUBLE, 1, "m/s");
    COM_new_dataitem( propBufAll+".positions",'n',COM_DOUBLE,3,"m");
    COM_new_dataitem( propBufAll+".constrained",'n',COM_INTEGER,1,"");
    COM_resize_array( propBufAll+".positions" );
    COM_resize_array( propBufAll+".constrained");
    COM_resize_array( propBufAll+".rb" );
    create_registered_window_dataitems( propBufAll);
    COM_window_init_done( propBufAll);

    COM_new_window( propBuf);
    std::string bc_str = surf_window+".bcflag";
    COM_use_dataitem( propBuf+".mesh", surf_window+".mesh", 1, bc_str.c_str(), 1);
    COM_use_dataitem( propBuf+".pconn", propBufAll+".pconn");
    if ( COM_get_dataitem_handle( surf_window+".cnstr_type") >0) {
       COM_use_dataitem( propBuf, surf_window+".cnstr_type");
    }
    COM_use_dataitem( propBuf+".vbar_alp", propBufAll+".vbar_alp");
    COM_use_dataitem( propBuf+".vbar", propBufAll+".vbar");
    COM_use_dataitem( propBuf+".vbar_old", propBufAll+".vbar_old");
    COM_use_dataitem( propBuf+".vbar_grad", propBufAll+".vbar_grad");
    COM_use_dataitem( propBuf+".rb", propBufAll+".rb");
    create_registered_window_dataitems( propBuf);
    COM_window_init_done( propBuf);
  }

  COM_new_window( solidBufBase);
  std::string bcflag = surf_window+".bcflag";
  COM_use_dataitem( solidBufBase+".mesh", surf_window+".mesh", 1, bcflag.c_str(), 0);
  COM_use_dataitem( solidBufBase+".mesh", surf_window+".mesh", 1, bcflag.c_str(), 1);
  if (withALE) {
    COM_new_dataitem(solidBufBase+".areas", 'e', COM_DOUBLE, 1, "m^2");
    COM_resize_array(solidBufBase+".areas");
    COM_new_dataitem(solidBufBase+".mdot", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
    COM_resize_array(solidBufBase+".mdot");
  }
  COM_clone_dataitem( solidBufBase+".ts", surf_window+".ts_alp");
  COM_clone_dataitem( solidBufBase+".ts_old", surf_window+".ts_alp");
  COM_clone_dataitem( solidBufBase+".ts_grad", surf_window+".ts_alp");
  create_registered_window_dataitems( solidBufBase);
  COM_window_init_done( solidBufBase);

  COM_new_window( solidBuf);
  COM_use_dataitem( solidBuf, solidBufBase+".all");
  COM_use_dataitem( solidBuf, surf_window+".data");
  if (withALE) 
    COM_use_dataitem( solidBuf, propBufAll+".data");
  if (with_fluid == 1) {    // coupled
    COM_use_dataitem( solidBuf+".x", solidBufBase+".nc");
    COM_clone_dataitem( solidBuf+".nc", surf_window+".nc");
  }
  create_registered_window_dataitems( solidBuf);
  COM_window_init_done( solidBuf);

  split_surface_window( isolid_all, isolid_i, isolid_nb, isolid_b, isolid_ni);

  COM_delete_window( isolid_b);
  COM_delete_window( isolid_nb);
  COM_delete_window( isolid_ni);

  // Create a window for output solid non-burning patch data
  COM_new_window( isolid_nb);
  COM_use_dataitem( isolid_nb+".mesh", surf_window+".mesh", 1, bcflag.c_str(), 0);
  COM_use_dataitem( isolid_nb, surf_window+".data");

  if ( withALE) {
    COM_use_dataitem( isolid_nb+".pconn", propBufAll+".pconn");
    COM_use_dataitem( isolid_nb, propBufAll+".vbar");
    COM_use_dataitem( isolid_nb, propBufAll+".vbar_old");
    COM_use_dataitem( isolid_nb, propBufAll+".rb");
  }

  COM_use_dataitem( isolid_nb, solidBufBase+".ts");
  COM_use_dataitem( isolid_nb, solidBufBase+".ts_old");
  COM_window_init_done( isolid_nb);

  // Create a window for output solid interface data
  COM_new_window( isolid_b);
  COM_use_dataitem( isolid_b+".mesh", surf_window+".mesh", 1, bcflag.c_str(), 1);
  COM_use_dataitem( isolid_b, surf_window+".data");

  if ( withALE) {
       COM_use_dataitem( isolid_b+".pconn", propBufAll+".pconn");
       COM_use_dataitem( isolid_b, propBufAll+".vbar");
       COM_use_dataitem( isolid_b, propBufAll+".vbar_old");
       COM_use_dataitem( isolid_b, propBufAll+".rb");
  }

  COM_use_dataitem( isolid_b, solidBufBase+".ts");
  COM_use_dataitem( isolid_b, solidBufBase+".ts_old");
  COM_window_init_done( isolid_b);

  //   Create a window for non-solid/fluid interface
  COM_new_window( isolid_ni);
  COM_use_dataitem( isolid_ni+".mesh", surf_window+".mesh", 1, bcflag.c_str(), 2);
  COM_use_dataitem( isolid_ni, surf_window+".data");
  if ( withALE) {
       COM_use_dataitem( isolid_ni+".pconn", propBufAll+".pconn");
       COM_use_dataitem( isolid_ni, propBufAll+".vbar");
       COM_use_dataitem( isolid_ni, propBufAll+".vbar_old");
       COM_use_dataitem( isolid_ni, propBufAll+".rb");
  }
  COM_window_init_done( isolid_ni);

    //  setup for pc iterations
  int maxPredCorr = get_coupling()->get_max_ipc();
  if ( maxPredCorr>1) {
       // Create window for backing surface data
       COM_new_window( solidBufBak);
       COM_use_dataitem( solidBufBak+".mesh", solidBuf+".mesh");
       COM_clone_dataitem( solidBufBak+".x", solidBuf+".x");
       COM_clone_dataitem( solidBufBak+".y", solidBuf+".nc");
       COM_window_init_done( solidBufBak);

       // Create window for backing up volume data
       COM_new_window( solidVolBak);
       COM_use_dataitem( solidVolBak+".mesh", vol_window+".mesh");
       COM_clone_dataitem( solidVolBak, vol_window+".nc");
       COM_clone_dataitem( solidVolBak, vol_window+".data");
       COM_window_init_done( solidVolBak);

       // Initlaize the dataitem handles to be stored/restored
       pc_hdls[0][0] = COM_get_dataitem_handle( vol_window+".all");
       pc_hdls[1][0] = COM_get_dataitem_handle( solidVolBak+".all");

       pc_hdls[0][1] = COM_get_dataitem_handle( solidBuf+".x");
       pc_hdls[1][1] = COM_get_dataitem_handle( solidBufBak+".x");
       pc_hdls[0][2] = COM_get_dataitem_handle( solidBuf+".nc");
       pc_hdls[1][2] = COM_get_dataitem_handle( solidBufBak+".y");
       pc_count = 3;
  }

  // update
  split_surface_window( isolid_all, isolid_i, isolid_nb, isolid_b, isolid_ni);

    // setup variables
  traction_mode = get_coupling()->get_rocmancontrol_param()->traction_mode;

  std::string unit;
  char loc;
  COM_get_dataitem( isolid_all+".rhos", &loc, &dummy, &dummy, &unit);
  if (loc=='w' || loc=='p')
    rhos_mode = 1;
  else {
    if (loc != 'n' && loc != 'e') {
      printf("Rocstar: Error: Unknown type of location of rohs: %c\n", loc);
      exit(1);
    }
    rhos_mode = 2;
  }

  COM_get_dataitem( solidBufBase+".ts", &loc, &dummy, &size_ts, &unit);
  if (size_ts == 1 && traction_mode != NO_SHEER) {
    COM_assertion_msg(0, "If traction mode is with sheer, then solid tractions must be vectors!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // for compute_distances
  y_hdl = COM_get_dataitem_handle( solidBuf+".nc");

  if (!get_coupling()->initial_start()) read_restart_data();
}

void SolidAgent::read_restart_data()
{
  int atts_hdl = COM_get_dataitem_handle_const( solidSurfIN+".data");
  int buf_hdl = COM_get_dataitem_handle( solidBufBase+".data");
  COM_call_function( obtain_attr_handle, &atts_hdl, &buf_hdl);

  if ( withALE) {
          // Obtain data for surface propagation if ALE is enabled
        int propBufAll_hdl = COM_get_dataitem_handle( propBufAll+".data");
        COM_call_function( obtain_attr_handle, &atts_hdl, &propBufAll_hdl);

          // Obtain pconn for surface propagation
        int pconn_hdl_const = COM_get_dataitem_handle_const( solidSurfIN+".pconn");
        int pconn_hdl = COM_get_dataitem_handle( solidSurfIN+".pconn");
        COM_call_function( obtain_attr_handle, &pconn_hdl_const, &pconn_hdl);
        COM_clone_dataitem( propBufAll+".pconn", solidSurfIN+".pconn");
  }
}

void SolidAgent::output_restart_files( double t) {
  static const std::string isolid_prefix = "isolid_all";

  // Write out surface sub-windows
#if 1
  write_data_files( t, isolid_b, isolid_b+".all");
  write_data_files( t, isolid_nb, isolid_nb+".all");
  write_data_files( t, isolid_ni, isolid_ni+".all");
#else
  write_data_files( t, isolid_prefix, isolid_all+".all");
#endif
  write_control_file( t, "isolid*", surf_window.c_str());

  // Write out volume window
  write_data_files( t, solid, (vol_window+".all").c_str());
  write_control_file( t, solid, vol_window.c_str());
}

// Visualization window buffers
static const char *isolid_vis = "isolid_vis";
static const char *solid_vis = "solid_vis";
static const char *solid_plag_vis = "solid_plag_vis";

static const char *isolid_b_vis = "isolid_b_vis";
static const char *isolid_nb_vis = "isolid_nb_vis";
static const char *isolid_ni_vis = "isolid_ni_vis";


void SolidAgent::output_visualization_files( double t) {
  // TODO: Define visualization sub-windows
}

void SolidAgent::finalize()
{
  int maxPredCorr = get_coupling()->get_max_ipc();
  if ( maxPredCorr>1) {
       COM_delete_window( solidBufBak);
       COM_delete_window( solidVolBak);
  }
                                                                                
  COM_delete_window( solidBufBase);
  COM_delete_window( isolid_nb);
  COM_delete_window( isolid_b);
  COM_delete_window( isolid_ni);
  COM_delete_window( solidBuf);
  COM_delete_window( isolid_i);
  COM_delete_window( isolid_all);
  if ( withALE) {
      COM_delete_window( propBufAll);
      COM_delete_window( propBuf);
  }

  Agent::finalize();
}

int SolidAgent::compute_integrals()
{
  if (compute_integrals_handle > 0) {
    MAN_DEBUG(3, ("Rocstar: SolidAgent::compute_integrals.\n"));
    COM_call_function(compute_integrals_handle, integrals);
  }
  return 1;
}






