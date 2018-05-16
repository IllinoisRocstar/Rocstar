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
// $Id: mptest.C,v 1.4 2008/12/06 08:45:28 mtcampbe Exp $

/** \file mptest.C Testing marker-particle method in parallel */

#include "com.h"
#include "IM_Reader.h"
#include <iostream>

using namespace std;

COM_EXTERN_MODULE( SimOut);
COM_EXTERN_MODULE( Simpal);
COM_EXTERN_MODULE( Rocprop);

int main(int argc, char *argv[]) {
  COM_init( &argc, &argv);

  if ( argc < 3) {
    std::cout << "Usage:\n\tTo run in serial: " << argv[0] 
	      << " <surffile> <hdffile> " << endl;
    std::cout << "\n\tTo run in parallel: <mpirun-command> " << argv[0] 
	      << " -com-mpi <Rocin control file> <hdfoutput-prefix> " << endl;
    exit(-1);
  }

  std::cout << "Reading surface mesh file \"" << argv[1] << '"' << endl;

  std::string fname(argv[1]), wname;
  string::size_type n0 = fname.find_last_of( "/");

  if ( n0 != std::string::npos) 
    fname = fname.substr( n0+1, fname.size());

  string::size_type ni;
  ni = fname.find_first_of( ".:-*[]?\\\"\'0123456789");
  COM_assertion_msg(ni, "File name must start with a letter");

  if ( ni == std::string::npos) {
    wname = fname;
  }
  else {
    while (fname[ni-1]=='_') --ni; // Remove the '_' at the end.
    wname = fname.substr( 0, ni);
  }

  std::cout << "Creating window \"" << wname << '"' << endl;

  // Read in IM/HDF format
  int err = IM_Reader().read_winmesh( argv[1], wname); 
  COM_assertion( err>=0);

  // Allocate memory for normals
  const string rb = wname+".rb";
  COM_new_dataitem(rb.c_str(), 'e', COM_DOUBLE, 1, "");
  COM_resize_array( rb.c_str());

  // Allocate an element and a nodal dataitem for testing elements_to_nodes
  const string disps = wname+".disps";
  COM_new_dataitem(disps.c_str(), 'e', COM_DOUBLE, 3, "m");
  COM_resize_array( disps.c_str());
  
  COM_window_init_done( wname.c_str());

  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocout, "OUT");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocprop, "PROP");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocblas, "BLAS");

  int PROP_init = COM_get_function_handle( "PROP.initialize");
  int mesh_hdl = COM_get_dataitem_handle( (wname+".mesh").c_str());

  COM_call_function( PROP_init, &mesh_hdl);

  int PROP_prop = COM_get_function_handle( "PROP.propagate");

  int rb_hdl = COM_get_dataitem_handle( rb.c_str());
  int disps_hdl = COM_get_dataitem_handle( disps.c_str());

  double spd = 1.e-1;
  int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
  COM_call_function( BLAS_copy_scalar, &spd, &rb_hdl);

  double dt = 1.e-6;
  COM_call_function( PROP_prop, &mesh_hdl, &rb_hdl, &dt, &disps_hdl);

  std::cout << "Output window into file..." << endl;

  // Output normals
  int OUT_set = COM_get_function_handle( "OUT.set_option");
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");

  COM_call_function( OUT_set, "mode", "w");
  int all_hdl = COM_get_dataitem_handle( (wname+".all").c_str());
  COM_call_function( OUT_write, argv[2], &all_hdl, 
		     (char*)wname.c_str(), "000");
  
  COM_finalize();
}






