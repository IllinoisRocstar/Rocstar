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
// $Id: Surfdiver.C,v 1.11 2009/10/08 15:36:00 mtcampbe Exp $

#include <iostream>
#include <cstring>
#include <cstdlib>

#include "com.h"
#include "basic_actions.h"
#include "FluidAgent.h"
#include "SolidAgent.h"

using namespace std;

extern void _load_rocface(FluidAgent *fagent, SolidAgent *sagent, const RocmanControl_parameters *param);

void read_file( const char *fname, const string &wname, double alpha)
{
  char *lastdot=strrchr( const_cast<char *>(fname), '.');

  COM_new_window( wname.c_str(), MPI_COMM_SELF);

  // Read in HDF files or a SimIN control file
  if(man_verbose > 2)
    std::cout << "Reading file " << fname << " " << wname << std::endl;

  // Read in HDF format
  //COM_UNLOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");
  //COM_LOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");
    
  int IN_read;
  // Read in HDF format using SimIN::read_window or ::read_by_control_file 
  if ( strcmp( lastdot, ".hdf")==0)
    IN_read = COM_get_function_handle( "IN.read_window");
  else
    IN_read = COM_get_function_handle( "IN.read_by_control_file");

  // Pass MPI_COMM_NULL to SimIN so that the rank becomes a wildcard.
  MPI_Comm comm_null = MPI_COMM_NULL;
  std::string bufwin("bufwin");
  COM_call_function( IN_read, fname, bufwin.c_str(), &comm_null);

  int IN_obtain = COM_get_function_handle( "IN.obtain_dataitem");

  // Check whether bcflag exists. If so, retain only the panes with flag<=1.
  int bcflag = COM_get_dataitem_handle((bufwin+".bcflag").c_str());
  if (bcflag > 0) {
    // Read in bcflags.
    COM_call_function( IN_obtain, &bcflag, &bcflag);
    
    // Obtain the IDs of the panes of the window
    int npanes, *pane_ids;
    COM_get_panes( bufwin.c_str(), &npanes, &pane_ids);
    
    // Loop through the panes to remove those with bcflag >1.
    for ( int i=0; i<npanes; ++i) {
      int *flag;
      COM_get_array( (bufwin+".bcflag").c_str(), pane_ids[i], &flag);
      if ( flag==NULL || *flag>1)
	COM_delete_pane( bufwin.c_str(), pane_ids[i]);
    }

    // remove buffers.
    COM_free_buffer( &pane_ids);
  }

  // Remove all dataitems except for the mesh
  COM_delete_dataitem(  (bufwin+".data").c_str());

  // Read in the mesh.
  int buf_mesh = COM_get_dataitem_handle((bufwin+".mesh").c_str());
  COM_call_function( IN_obtain, &buf_mesh, &buf_mesh);
  //COM_UNLOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");
 
  if(man_verbose > 2) 
    std::cout << "Obtained window " << wname << " from file " << fname << std::endl;

  // Change the memory layout to contiguous.
  COM_clone_dataitem( (wname+".mesh").c_str(), (bufwin+".mesh").c_str(), 0);
  COM_delete_window( bufwin.c_str());
}

SurfDiver::SurfDiver(FluidAgent *fag, SolidAgent *sag):
  Action( 0, (const char**)NULL, NULL, NULL, (char *)"SurfDiver"),
  fagent(fag), sagent(sag)
{
  outdir = "Rocman/"+fagent->get_rocmod_name()+sagent->get_rocmod_name()+"/";
}

void SurfDiver::init(double t) 
{
  fluid_mesh = COM_get_dataitem_handle_const( fagent->fluidBufNG+".mesh");
  solid_mesh = COM_get_dataitem_handle_const( sagent->solidBuf+".mesh");

  fluid_mesh_str = fagent->get_rocmod_name()+".mesh";
  solid_mesh_str = sagent->get_rocmod_name()+".mesh";

  _load_rocface(fagent, sagent, fagent->get_coupling()->get_rocmancontrol_param());
  RFC_readcntr = COM_get_function_handle( "RFC.read_control_file");
  RFC_overlay = COM_get_function_handle( "RFC.overlay");
  RFC_write = COM_get_function_handle( "RFC.write_overlay");
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
  RFC_interpolate = COM_get_function_handle("RFC.interpolate");
  RFC_read = COM_get_function_handle( "RFC.read_overlay");
}

void SurfDiver::run( double t, double dt, double alpha) {
  MAN_DEBUG(1, ("[%d] Rocstar: SurfDiver::run() with t:%e dt:%e.\n", fagent->get_comm_rank(), t, dt));

  MPI_Comm comm = fagent->get_communicator();

  //  dump meshes
  MAN_DEBUG(1, ("[%d] Rocstar: SurfDiver::run() dumping output files for time %e.\n", fagent->get_comm_rank(), t));
  fagent->output_restart_files( t);
  sagent->output_restart_files( t);

  MPI_Barrier(comm);

    // run sequentially
  if (fagent->get_comm_rank() == 0) {

    MPI_Comm oldcomm = COM_get_default_communicator();
    COM_set_default_communicator( MPI_COMM_NULL);
      // need to turn off profiling as it may hang for npes > 1
    COM_set_profiling(0);

    //  read meshes 

    std::string time_str;
    fagent->get_time_string(t, time_str);
    std::string fluid_file = fagent->get_rocmod_name()+"/Rocout/ifluid_in_"+time_str+".txt";
    std::string solid_file = sagent->get_rocmod_name()+"/Rocout/isolid_in_"+time_str+".txt";

    std::string fluid_wname = outdir+"ifluid";
    std::string solid_wname = outdir+"isolid";

    read_file( fluid_file.c_str(), fluid_wname.c_str(), 1.);
    COM_window_init_done( fluid_wname);

    read_file( solid_file.c_str(), solid_wname.c_str(), 1.);
    COM_window_init_done( solid_wname);

    int fluid_mesh1 = COM_get_dataitem_handle( (fluid_wname+".mesh").c_str());
    int solid_mesh1 = COM_get_dataitem_handle( (solid_wname+".mesh").c_str());

    const char *format = "HDF";
  
    // mesh overlay
    MAN_DEBUG(2,("Starting mesh overlay..."));
    COM_call_function( RFC_overlay, &fluid_mesh1, &solid_mesh1);

    // output overlay mesh
    COM_call_function( RFC_write, &fluid_mesh1, &solid_mesh1, 
		     fluid_wname.c_str(), solid_wname.c_str(), format);

    COM_set_default_communicator( oldcomm);
    COM_set_profiling(1);
  }   // end of PE 0

  MPI_Barrier(comm);

  // reload overlay
  std::string fluid_dir = outdir+"ifluid";
  std::string solid_dir = outdir+"isolid";

  MAN_DEBUG(2,("Reload partitioned mesh overlay... "));
  COM_call_function( RFC_read, &fluid_mesh, &solid_mesh, &comm, fluid_dir.c_str(), solid_dir.c_str(), "HDF");
}

//
void compute_overlay( FluidAgent *fagent, SolidAgent *sagent, double t) {
  MAN_DEBUG(1, ("[%d] Rocstar: compute_overlay with t:%e .\n", fagent->get_comm_rank(), t));

  MPI_Comm comm = fagent->get_communicator();

  MPI_Barrier(comm);

    // run sequentially
  if (fagent->get_comm_rank() == 0) {

    MPI_Comm oldcomm = COM_get_default_communicator();
    COM_set_default_communicator( MPI_COMM_NULL);
      // need to turn off profiling as it may hang for npes > 1
    COM_set_profiling(0);

    //  read meshes 

    std::string time_str;
    fagent->get_time_string(t, time_str);
    std::string fluid_file = fagent->get_rocmod_name()+"/Rocout/ifluid_"+time_str+".txt";
    std::string solid_file = sagent->get_rocmod_name()+"/Rocout/isolid_"+time_str+".txt";

    std::string outdir = "Rocman/"+fagent->get_rocmod_name()+sagent->get_rocmod_name()+"/";
    std::string fluid_wname = outdir+"ifluid";
    std::string solid_wname = outdir+"isolid";

      // load two meshes
    read_file( fluid_file.c_str(), fluid_wname.c_str(), 1.);
    COM_window_init_done( fluid_wname);

    read_file( solid_file.c_str(), solid_wname.c_str(), 1.);
    COM_window_init_done( solid_wname);

    int fluid_mesh1 = COM_get_dataitem_handle( (fluid_wname+".mesh").c_str());
    int solid_mesh1 = COM_get_dataitem_handle( (solid_wname+".mesh").c_str());

    const char *format = "HDF";
  
    // call Rocblas to get deformed data HERE
    int s_x_hdl = COM_get_dataitem_handle( sagent->solidBuf + ".x");
    int s_uhat_hdl = COM_get_dataitem_handle( sagent->solidBuf + ".uhat");
    int s_y_hdl = COM_get_dataitem_handle( sagent->solidBuf + ".nc");

      // get deformed
    COM_call_function( RocBlas::add, &s_x_hdl, &s_uhat_hdl, &s_y_hdl);

    int RFC_overlay = COM_get_function_handle( "RFC.overlay");
    int RFC_write = COM_get_function_handle( "RFC.write_overlay");

    // mesh overlay
    MAN_DEBUG(2,("Starting mesh overlay... "));
    COM_call_function( RFC_overlay, &fluid_mesh1, &solid_mesh1);

    // output overlay mesh
    COM_call_function( RFC_write, &fluid_mesh1, &solid_mesh1, 
		     fluid_wname.c_str(), solid_wname.c_str(), format);

    COM_set_default_communicator( oldcomm);
    COM_set_profiling(1);
  }   // end of PE 0

  MPI_Barrier(comm);

}

void SurfDiverAfterRemeshing::run( double t, double dt, double alpha) 
{
  compute_overlay( fagent, sagent, t);
}







