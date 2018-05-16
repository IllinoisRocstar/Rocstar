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

#include "sys/types.h"
#include "sys/stat.h"

#include "rocman.h"
#include "FluidAgent.h"
#include "Interpolate.h"

#ifdef STATIC_LINK
extern "C" void COM_F_FUNC2(rocflo_load_module, ROCFLO_LOAD_MODULE)( const char *, int);
extern "C" void COM_F_FUNC2(rocflo_unload_module, ROCFLO_UNLOAD_MODULE)( const char *, int);
extern "C" void COM_F_FUNC2(rocflu_load_module, ROCFLU_LOAD_MODULE)( const char *, int);
extern "C" void COM_F_FUNC2(rocflu_unload_module, ROCFLU_UNLOAD_MODULE)( const char *, int);
#endif

//const char *FluidAgent::window_name = "FluidAgent";

FluidAgent::FluidAgent(Coupling *cp, std::string mod, std::string obj, MPI_Comm com, int withSolid) :
    Agent(cp, mod, obj, "Fluid", com, true), 
    with_plag(false), with_solid(withSolid)
{
  load_module();

  ifluid = "ifluid";	       // for input
  fluid = "fluid";	    
  fluid_plag = "fluid_plag";

  // SimIN windows, Input buffer windows
  fluidSurfIn = "fluidSurfIn";
  fluidVolIn = "fluidVolIn";
  fluidPlagIn = "fluidPlagIn";
  fluidVPIn = "fluidVolIn fluidPlagIn";

  ifluid_all = "ifluid_all";	// tmp buffer
  ifluid_i = "ifluid";	       // FluidBuf

    // Surface windows for Rocout, Surface window buffers
  ifluid_b = "ifluid_b";		// fluidBufBOUT
  ifluid_nb = "ifluid_nb";		// fluidBufNBOUT
  ifluid_ni = "ifluid_ni";		// fluidBufNIOUT

  propBufAll = "FluidPropAll";
  fluidBufNG = "FluidBufNG";		// fluidBufNG
  propBuf = "FluidProp";
  fluidBufB      = "FluidBufB";
  fluidBufNB      = "FluidBufNB";

  fluidBufPC     = "FluidBufPC";        // surface data to be backed up
  fluidBufBak     = "fluidBufBak";      // back up PC iterations
  fluidBufPRE    = "FluidBufPRE";       //  for convergence check

  fluidVolBak    = "FluidVolBak";       // backup of volume data
  fluidPlagBak    = "FluidPlagBak";     // backup of volume data

  tmp_window = ifluid_all;

  if (!with_solid) dobackup = 0;
}

void FluidAgent::load_module() 
{
  MAN_DEBUG(3, ("[%d] Rocstar: FluidAgent::load_module %s %s.\n", comm_rank, rocmod_name.c_str(), mod_instance.c_str()));

#ifdef STATIC_LINK
# ifdef RFLO
  if (rocmod_name == "Rocflo")
    COM_F_FUNC2( rocflo_load_module, ROCFLO_LOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
# elif defined(RFLU)
  if (rocmod_name == "Rocflu")
    COM_F_FUNC2( rocflu_load_module, ROCFLU_LOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
# else
    COM_assertion_msg(0, "Unknown FluidAgent mod!");
# endif
#else     // dynamic loading
  COM_assertion_msg(rocmod_name == "Rocflo" || rocmod_name == "Rocflu", 
               (std::string("Unknown FluidAgent module:")+rocmod_name).c_str());
  COM_load_module(rocmod_name.c_str(), mod_instance.c_str());
#endif

  init_function_handles();	// defined in Agent
}

void FluidAgent::unload_module() 
{
  MAN_DEBUG(3, ("Rocstar: FluidAgent::unload_module %s.\n", rocmod_name.c_str()));
#ifdef STATIC_LINK
# ifdef RFLO
  if (rocmod_name == "Rocflo")
    COM_F_FUNC2( rocflo_unload_module, ROCFLO_UNLOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
# else
#  ifdef RFLU
  if (rocmod_name == "Rocflu")
    COM_F_FUNC2( rocflu_unload_module, ROCFLU_UNLOAD_MODULE)( mod_instance.c_str(), mod_instance.length());
#  endif
# endif
#else
    // in restarting, close_module does not dlclose the shared lib
  if (get_coupling()->in_restart())
    COM_close_module(rocmod_name.c_str(), mod_instance.c_str());
  else
    COM_unload_module(rocmod_name.c_str(), mod_instance.c_str());
#endif
}

void FluidAgent::input( double t) {
  
  MPI_Barrier(communicator);
  //  if(!comm_rank)
  //    std::cout << "Fluid agent reading surface input." << std::endl;
  read_by_control_file( t, ifluid, fluidSurfIn);
  //  std::cout << comm_rank << ": done reading surface." << std::endl;
  MPI_Barrier(communicator);
  //  if(!comm_rank)
  //    std::cout << "Fluid agent done reading surface input, reading volume" << std::endl;
  read_by_control_file( t, fluid, fluidVolIn);
  //  std::cout << comm_rank << ": done reading volume." << std::endl;
  MPI_Barrier(communicator);
  //  if(!comm_rank)
  //    std::cout << "Fluid agent done reading volume" << std::endl;

#ifndef NATIVE_MP_IO
  if (read_by_control_file( t, fluid_plag, fluidPlagIn)==-1) {
    COM_new_window(fluidPlagIn);
    COM_window_init_done(fluidPlagIn);
  }
#else
  COM_new_window(fluidPlagIn);
  COM_window_init_done(fluidPlagIn);
#endif

  MPI_Barrier(communicator);
  //  std::cout << "All processors made it past reading" << std::endl;
}

// called in Coupling::init()
void FluidAgent::init_module( double t, double dt) {
  MAN_DEBUG(3, ("Rocstar: FluidAgent::init_module t=%e dt=%e.\n", t, dt));

  Agent::init_module(t, dt);

  // Call initialization routine of physics module
  // was in Rocman.f90
  COM_call_function( init_handle, &t, &communicator, &ic_handle, 
		     fluidSurfIn.c_str(), fluidVPIn.c_str(), 
                     &obtain_attr_handle);

  // Delete input buffer windows
  COM_delete_window( fluidSurfIn);
  COM_delete_window( fluidVolIn);
#ifndef NATIVE_MP_IO
  COM_delete_window( fluidPlagIn);
#endif
  // INITIALIZE_XXX common portion

  // Reassign vol_window
  std::string::size_type pos = vol_window.find( " ");
  if ( pos != std::string::npos) {
    plag_window = vol_window.substr( pos+1, vol_window.size());
    vol_window = vol_window.substr( 0, pos);
  }
  with_plag = !plag_window.empty();
#ifdef NATIVE_MP_IO
  with_plag = false;
#endif
  MAN_DEBUG(3, ("Rocstar: vol_window = %s with_plag = %d.\n", vol_window.c_str(), with_plag));

  // Split surface window for output
  // ifluid is fluidBuf in INITIALIZE_FLUID() ???
  //split_surface_window( ifluid_all, ifluid_i, ifluid_nb, ifluid_b, ifluid_ni);
/*
  std::string bcflag = surf_window+".bcflag";
  COM_use_dataitem( ifluid_i, surf_window+".mesh", 1, bcflag.c_str(), 0);
  COM_use_dataitem( ifluid_i, surf_window+".mesh", 1, bcflag.c_str(), 1);
  COM_window_init_done( ifluid_i);
*/
}

// called in init callback (ic_handle)
void FluidAgent::create_buffer_all()
{
  Agent::create_buffer_all();

  COM_new_window( propBufAll);
  COM_use_dataitem( propBufAll, surf_window+".mesh",  1);
  //  COM_use_dataitem( propBufAll, surf_window+".bcflag");
    // Mesh motion velocities
  COM_new_dataitem( propBufAll+".vm", 'n', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem( propBufAll+".rb", 'e', COM_DOUBLE, 1, "m/s");
  COM_new_dataitem( propBufAll+".positions",'n',COM_DOUBLE,3,"m");
  COM_new_dataitem( propBufAll+".cflag",'n',COM_INTEGER,1,"");
  COM_resize_array(propBufAll+".data");
  if ( COM_get_dataitem_handle( surf_window+".cnstr_type") >0) {
    COM_use_dataitem( propBufAll, surf_window+".bflag", 1);
    COM_use_dataitem( propBufAll, surf_window+".cnstr_type");
    COM_use_dataitem( propBufAll, surf_window+".bcflag");
     // a temporary buffer to store ".cnstr_type" of fluid at restart
    COM_clone_dataitem( propBufAll+".cnstr_type_bak", surf_window+".cnstr_type");
  }
  create_registered_window_dataitems( propBufAll);
  COM_window_init_done( propBufAll);

  std::string bcflag = surf_window+".bcflag";
    // create a window for surface propagation if run in fluid-alone
  if ( !with_solid)  {
    COM_new_window( propBuf);
    COM_use_dataitem( propBuf, surf_window+".mesh", 0, bcflag.c_str(), 1);
    COM_use_dataitem( propBuf, propBufAll+".pconn");
    COM_use_dataitem( propBuf, propBufAll+".data");
    create_registered_window_dataitems(propBuf);
    COM_window_init_done( propBuf);
  }

  // update
  //split_surface_window( ifluid_all, ifluid_i, ifluid_nb, ifluid_b, ifluid_ni);

  // ifluid_i is the fluidBuf in old rocman
  // Precondition: The window fluidSurf should have been defined and contain
  // the following variables: pf, qc, qr, rhof_alp, nf_alp, tf, Tf, mdot_alp,
  // Tflm_alp, du_alp, and rhofvf_alp.

  COM_new_window( ifluid_i);
  COM_use_dataitem( ifluid_i, surf_window+".mesh", 1, bcflag.c_str(), 0);
  COM_use_dataitem( ifluid_i, surf_window+".mesh", 1, bcflag.c_str(), 1);

  COM_new_dataitem(ifluid_i+".rhos", 'e', COM_DOUBLE, 1, "kg/(m^3)");
  COM_new_dataitem(ifluid_i+".mdot", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  COM_new_dataitem(ifluid_i+".mdot_old", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  COM_new_dataitem(ifluid_i+".mdot_grad", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  COM_new_dataitem(ifluid_i+".mdot_tmp", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  COM_resize_array(ifluid_i+".rhos");
  COM_resize_array(ifluid_i+".mdot");
  COM_resize_array(ifluid_i+".mdot_old");
  COM_resize_array(ifluid_i+".mdot_grad");
  COM_resize_array(ifluid_i+".mdot_tmp");

  COM_new_dataitem( ifluid_i+".vs", 'e', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem( ifluid_i+".vs_alp", 'e', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem( ifluid_i+".vs_old", 'e', COM_DOUBLE, 3, "m/s");
  COM_resize_array( ifluid_i+".vs");
  COM_resize_array( ifluid_i+".vs_alp");
  COM_resize_array( ifluid_i+".vs_old");

    // for compute distances
  if ( with_solid) {
    COM_clone_dataitem( ifluid_i+".nc_tmp", surf_window+".nc");
    COM_new_dataitem( ifluid_i+".sq_dist", 'n', COM_DOUBLE, 1, "m");
    COM_resize_array( ifluid_i+".sq_dist");
  }

  create_registered_window_dataitems(ifluid_i);
//  COM_resize_array( ifluid_i+".data");
  COM_window_init_done( ifluid_i);

  COM_new_window( fluidBufNG);
  COM_use_dataitem( fluidBufNG, ifluid_i+".mesh", 0, "", 0);
  COM_use_dataitem( fluidBufNG, surf_window+".data");
  COM_use_dataitem( fluidBufNG, ifluid_i+".data");
  COM_use_dataitem( fluidBufNG, propBufAll+".data");
  create_registered_window_dataitems(fluidBufNG);
  COM_window_init_done( fluidBufNG);

  COM_new_window( fluidBufB);
  COM_use_dataitem( fluidBufB, surf_window+".mesh", 0, bcflag.c_str(), 1);
  COM_use_dataitem( fluidBufB, surf_window+".data");
  COM_use_dataitem( fluidBufB, ifluid_i+".data");
  COM_use_dataitem( fluidBufB, propBufAll+".data");
  create_registered_window_dataitems(fluidBufB);
  COM_window_init_done( fluidBufB);

  COM_new_window( fluidBufNB);
  COM_use_dataitem( fluidBufNB, surf_window+".mesh", 0, bcflag.c_str(), 0);
  COM_use_dataitem( fluidBufNB, surf_window+".data");
  COM_use_dataitem( fluidBufNB, ifluid_i+".data");
  COM_use_dataitem( fluidBufNB, propBufAll+".data");
  create_registered_window_dataitems(fluidBufNB);
  COM_window_init_done( fluidBufNB);

  //COM_delete_window( ifluid_b);
  //COM_delete_window( ifluid_nb);
  //COM_delete_window( ifluid_ni);
  
  //   Create windows for output of burning patches, including ghost nodes/cells
  //   (variables include all of fluidSurf, and mdot, mdot_old, vm,
  //    vs, vs_old, ts, nc_t0, sq_dist)
  COM_new_window( ifluid_b);
  COM_use_dataitem( ifluid_b, surf_window+".mesh", 1, bcflag.c_str(), 1);
  COM_use_dataitem( ifluid_b, surf_window+".data");

  COM_use_dataitem( ifluid_b, propBufAll+".pconn");
  COM_use_dataitem( ifluid_b, propBufAll+".vm");

  if ( with_solid) {
       COM_use_dataitem( ifluid_b, ifluid_i+".vs");
       COM_use_dataitem( ifluid_b, ifluid_i+".vs_old");
       if (COM_get_dataitem_handle( ifluid_i+".ts")>0) {
         COM_use_dataitem( ifluid_b, ifluid_i+".ts");
         COM_use_dataitem( ifluid_b, ifluid_i+".nc_t0");
         COM_use_dataitem( ifluid_b, ifluid_i+".sq_dist");
       }
  } 

  COM_use_dataitem( ifluid_b, ifluid_i+".mdot");
  COM_use_dataitem( ifluid_b, ifluid_i+".mdot_old");

  COM_window_init_done( ifluid_b);

    //  Create windows for output of nonburning patches, including ghosts.
    // (variables include all of fluidSurf, and vm, vm_old, vs, vs_old,
    //  ts, nc_t0, sq_dist)
  COM_new_window( ifluid_nb);
  COM_use_dataitem( ifluid_nb, surf_window+".mesh", 1, bcflag.c_str(),0);
  COM_use_dataitem( ifluid_nb, surf_window+".data");

  COM_use_dataitem( ifluid_nb, propBufAll+".pconn");
  COM_use_dataitem( ifluid_nb, propBufAll+".vm");

  if ( with_solid) {
       COM_use_dataitem( ifluid_nb, ifluid_i+".vs");
       COM_use_dataitem( ifluid_nb, ifluid_i+".vs_old");
       if (COM_get_dataitem_handle( ifluid_i+".ts")>0) {
         COM_use_dataitem( ifluid_nb, ifluid_i+".ts");
         COM_use_dataitem( ifluid_nb, ifluid_i+".nc_t0");
         COM_use_dataitem( ifluid_nb, ifluid_i+".sq_dist");
       }
  }
  COM_window_init_done( ifluid_nb);

  //   Create windows for output of nonburning patches, including ghosts.
  //   (variables include all of fluidSurf)
  COM_new_window( ifluid_ni);
  COM_use_dataitem( ifluid_ni, surf_window+".mesh", 1, bcflag.c_str(),2);
  COM_use_dataitem( ifluid_ni, surf_window+".data");

  COM_use_dataitem( ifluid_ni, propBufAll+".pconn");
  COM_use_dataitem( ifluid_ni, propBufAll+".vm");

  COM_window_init_done( ifluid_ni);

    // setup for pc iterations
  int maxPredCorr = get_coupling()->get_max_ipc();
  if ( maxPredCorr>1) {
         // Create a window to encapsulate the surface data to be backed up
       COM_new_window( fluidBufPC);
       COM_use_dataitem( fluidBufPC, ifluid_i+".mesh");
       COM_use_dataitem( fluidBufPC, ifluid_i+".nc");
       COM_window_init_done( fluidBufPC);
                                                                                
         // Create a window to store backed-up surface data
       COM_new_window( fluidBufBak);
       COM_use_dataitem( fluidBufBak, fluidBufPC+".mesh");
       COM_clone_dataitem( fluidBufBak, fluidBufPC+".nc");
       COM_window_init_done( fluidBufBak);
                                                                                
         // Create window for backing up volume data
       COM_new_window( fluidVolBak);
       COM_use_dataitem( fluidVolBak, vol_window+".mesh");
       COM_clone_dataitem( fluidVolBak, vol_window+".nc");
       COM_clone_dataitem( fluidVolBak, vol_window+".data");
       COM_window_init_done( fluidVolBak);
  
         // Create window for backing up Plag data
       if ( with_plag) {
          COM_new_window( fluidPlagBak);
          COM_use_dataitem( fluidPlagBak, plag_window+".mesh");
          COM_clone_dataitem( fluidPlagBak, plag_window+".nc");
          COM_clone_dataitem( fluidPlagBak, plag_window+".data");
          COM_window_init_done( fluidPlagBak);
       }

         // Create window for convergence check
       COM_new_window( fluidBufPRE);
       COM_use_dataitem( fluidBufPRE, fluidBufNG+".mesh");
       COM_clone_dataitem( fluidBufPRE, fluidBufNG+".vm");
       COM_clone_dataitem( fluidBufPRE, fluidBufNG+".vs");
       COM_clone_dataitem( fluidBufPRE, fluidBufNG+".mdot");
       COM_clone_dataitem( fluidBufPRE, fluidBufNG+".ts");
       COM_window_init_done( fluidBufPRE);

        // Initlaize the dataitem handles to be stored/restored
       pc_hdls[0][0] = COM_get_dataitem_handle( fluidBufPC+".all");
       pc_hdls[1][0] = COM_get_dataitem_handle( fluidBufBak+".all");
                                                                                
       pc_hdls[0][1] = COM_get_dataitem_handle( vol_window+".all");
       pc_hdls[1][1] = COM_get_dataitem_handle( fluidVolBak+".all");
                                                                                
       pc_count = 2;
       if ( with_plag) {
          pc_hdls[0][2] = COM_get_dataitem_handle( plag_window+".all");
          pc_hdls[1][2] = COM_get_dataitem_handle( fluidPlagBak+".all");
          pc_count ++;
       }
       f_mdot_hdl = COM_get_dataitem_handle( fluidBufNG+".mdot");
       f_mdot_pre_hdl = COM_get_dataitem_handle( fluidBufPRE+".mdot");
       f_ts_hdl = COM_get_dataitem_handle( fluidBufNG+".ts");
       f_ts_pre_hdl = COM_get_dataitem_handle( fluidBufPRE+".ts");
       f_vm_hdl = COM_get_dataitem_handle( fluidBufNG+".vm");
       f_vm_pre_hdl = COM_get_dataitem_handle( fluidBufPRE+".vm");
  }

  // Split surface window for output
  // ifluid is fluidBuf in INITIALIZE_FLUID() ???
  // update
//  split_surface_window( ifluid_all, ifluid_i, ifluid_nb, ifluid_b, ifluid_ni);

  // for compute_distances
  nc_hdl = COM_get_dataitem_handle( fluidBufNG+".nc");
  nc_tmp_hdl = COM_get_dataitem_handle( fluidBufNG+".nc_tmp");
  sq_dist_hdl = COM_get_dataitem_handle( fluidBufNG+".sq_dist");

  if (!get_coupling()->initial_start()) {
       int cnstr_type_hdl = COM_get_dataitem_handle( surf_window+".cnstr_type");
       int cnstr_type_bak_hdl = COM_get_dataitem_handle( propBufAll+".cnstr_type_bak");
       if (cnstr_type_hdl > 0) {
           // backup ".cnstr_type" of fluid and restore it later
           // this will discard the values in HDF restart files!
         COM_call_function( RocBlas::copy, &cnstr_type_hdl, &cnstr_type_bak_hdl);
       }
       read_restart_data();
       if (cnstr_type_hdl > 0) {
         COM_call_function( RocBlas::copy, &cnstr_type_bak_hdl, &cnstr_type_hdl);
       }
  }
}

void FluidAgent::read_restart_data()
{
    // Obtain data for fluidBuf
  int atts_hdl = COM_get_dataitem_handle_const(fluidSurfIn+".data");
  int buf_hdl = COM_get_dataitem_handle( ifluid_i+".data");
  COM_call_function( obtain_attr_handle, &atts_hdl, &buf_hdl);

    // Obtain data for propBufAll
  int propBufAll_hdl = COM_get_dataitem_handle( propBufAll+".data");
  COM_call_function( obtain_attr_handle, &atts_hdl, &propBufAll_hdl);

    // Obtain pconn for surface propagation
  int pconn_hdl_const = COM_get_dataitem_handle_const( fluidSurfIn+".pconn");
  int pconn_hdl = COM_get_dataitem_handle( fluidSurfIn+".pconn");
  COM_call_function( obtain_attr_handle, &pconn_hdl_const, &pconn_hdl);
  COM_clone_dataitem( propBufAll, fluidSurfIn+".pconn");
}

void FluidAgent::output_restart_files( double t) {
  static const std::string ifluid_prefix = "ifluid_all";

  if (pre_hdf_handle != -1)
    COM_call_function( pre_hdf_handle);

  // Write out surface sub-windows
#if 1
  write_data_files( t, ifluid_b, ifluid_b+".all");
  write_data_files( t, ifluid_nb, ifluid_nb+".all");
  write_data_files( t, ifluid_ni, ifluid_ni+".all");
#else
  write_data_files( t, ifluid_prefix, ifluid_all+".all");
#endif
  write_control_file( t, "ifluid*", surf_window.c_str());

  // Write out volume window
  write_data_files( t, fluid, (vol_window+".all").c_str());
  write_control_file( t, fluid, vol_window.c_str());

  // Write out Plag window
  if ( with_plag) {
    write_data_files( t, fluid_plag, (plag_window+".all").c_str());
    write_control_file( t, fluid_plag, plag_window.c_str());
  }

  if (post_hdf_handle != -1)
    COM_call_function( post_hdf_handle);

}

// Visualization window buffers
static const char *ifluid_vis = "ifluid_vis";
static const char *fluid_vis = "fluid_vis";
static const char *fluid_plag_vis = "fluid_plag_vis";

static const char *ifluid_b_vis = "ifluid_b_vis";
static const char *ifluid_nb_vis = "ifluid_nb_vis";
static const char *ifluid_ni_vis = "ifluid_ni_vis";


void FluidAgent::output_visualization_files( double t) {
  // TODO: Define visualization sub-windows
}

void FluidAgent::finalize()
{
  if ( !with_solid) 
       COM_delete_window( propBuf);
  COM_delete_window( propBufAll);
  COM_delete_window( ifluid_all);
  COM_delete_window( ifluid_i);
  COM_delete_window( fluidBufB);
  COM_delete_window( fluidBufNG);
  COM_delete_window( ifluid_nb);
  COM_delete_window( ifluid_b);
  COM_delete_window( ifluid_ni);
                                                                                
  int maxPredCorr = get_coupling()->get_max_ipc();
  if ( maxPredCorr>1) {
       COM_delete_window( fluidBufPC);
       COM_delete_window( fluidBufBak);
       COM_delete_window( fluidVolBak);
       if ( with_plag) COM_delete_window( fluidPlagBak);
       COM_delete_window( fluidBufPRE);
  }

  Agent::finalize();
}

void FluidAgent::init_convergence( int iPredCorr) 
{
  MAN_DEBUG(2, ("Rocstar: FluidAgent::init_convergence at %d.\n", iPredCorr));

    // Copy current solution to pre for convergence check
  COM_call_function( RocBlas::copy, &f_mdot_hdl, &f_mdot_pre_hdl);

  COM_call_function( RocBlas::copy, &f_ts_hdl, &f_ts_pre_hdl);

  COM_call_function( RocBlas::copy, &f_vm_hdl, &f_vm_pre_hdl);

    // STORE_SOLUTIONS
  Agent::init_convergence( iPredCorr);
}

int FluidAgent::check_convergence( double tolerMass, double tolerTract, double tolerVelo) 
{
  MAN_DEBUG(3, ("Rocstar: FluidAgent::check_convergence .\n"));

  int result = 0;
  if (!check_convergence_help(f_vm_hdl, f_vm_pre_hdl, tolerVelo, "vm")) return result;
  if (!check_convergence_help(f_ts_hdl, f_ts_pre_hdl, tolerTract, "ts")) return result;
  if (!check_convergence_help(f_mdot_hdl, f_mdot_pre_hdl, tolerMass, "mdot")) return result;
  return 1;
}

int FluidAgent::compute_integrals()
{
  if (compute_integrals_handle > 0) {
    MAN_DEBUG(3, ("Rocstar: FluidAgent::compute_integrals.\n"));
    COM_call_function(compute_integrals_handle, integrals);
  }
  return 1;
}








