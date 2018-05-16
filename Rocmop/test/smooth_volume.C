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
// $Id: smooth_volume.C,v 1.8 2008/12/06 08:45:25 mtcampbe Exp $

#include "com.h"
#include "mapbasic.h"
#include "Rocblas.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>

using namespace std;

// Necessary for handling modules in static or dynamic fashion
COM_EXTERN_MODULE( SimIN);
COM_EXTERN_MODULE( SimOUT);
COM_EXTERN_MODULE( Rocmop);
COM_EXTERN_MODULE( SurfMap);
COM_EXTERN_MODULE( Simpal);

// Get the MPI rank of the process
int get_comm_rank( MPI_Comm comm) {
  int rank;
  int ierr = MPI_Comm_rank( comm, &rank); assert( ierr == 0);
  return rank;
}

// Get the size of the MPI communicator
int get_comm_size( MPI_Comm comm) {
  int size;
  int ierr = MPI_Comm_size( comm, &size); assert( ierr == 0);
  return size;
}


int main(int argc, char *argv[]) {
  MPI_Init( &argc, &argv);

  if ( argc < 5) {
    std::cout << "Usage: " << argv[0] 
	      << " <Rocin input file><Output file name><Niter><Invert Tets>" << endl;
    exit(-1);
  }
  string inputfname(argv[1]);
  string outputfname(argv[2]);
  int niter = atoi(argv[3]);
  int invert_tets = atoi(argv[4]);
  
  std::cout << "niter = " << niter << "\ninvert_tets = " << invert_tets << "\n\n";

  MPI_Comm comm = MPI_COMM_WORLD;
  int size_of_p = get_comm_size(comm);
  int myrank = get_comm_rank(comm);
    
  COM_init( &argc, &argv);

  // LOAD MODULES
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimOUT, "OUT");
  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocmop, "MOP");
  COM_LOAD_MODULE_STATIC_DYNAMIC( SurfMap, "MAP");
  COM_LOAD_MODULE_STATIC_DYNAMIC( Simpal, "BLAS");

  // GET FUNCTION HANDLES
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int OUT_add = COM_get_function_handle( "OUT.add_dataitem");
  int IN_read = COM_get_function_handle( "IN.read_by_control_file");
  int MOP_init = COM_get_function_handle( "MOP.initialize");
  int MOP_det = COM_get_function_handle( "MOP.determine_physical_border");
  int MOP_set = COM_get_function_handle ( "MOP.set_value");
  int MOP_smth = COM_get_function_handle( "MOP.smooth");
  int MOP_asp = COM_get_function_handle( "MOP.add_aspect_ratios");
  int BLAS_add = COM_get_function_handle("BLAS.add");

  // SET OUTPUT LEVELS
  COM_set_verbose(11);
  COM_set_profiling(1);

  // BUILD FILE and WINDOW NAMES
  ostringstream rankstream;
  rankstream << myrank+1;
  //rankstream << myrank;
  const string wname("smoothed_mesh");

  // READ IN FILES
  COM_call_function( IN_read,
		     inputfname.c_str(), 
		     wname.c_str(),
		     &comm);

  COM_new_dataitem((wname+".disp").c_str(),'n',COM_DOUBLE,3,"m");
  COM_resize_array((wname+".disp").c_str(),0,NULL,3);
  COM_window_init_done(wname.c_str());

  // SMOOTH THE WINDOW
  std::cout << "Smoothing the window\n";
  ostringstream MyStream;
  int method = 0;
  COM_call_function(MOP_set, "method", &method);
  int verb = 2;
  COM_call_function(MOP_set, "verbose", &verb);
  COM_call_function(MOP_set, "inverted", &invert_tets);

  int VOL_pmesh = COM_get_dataitem_handle((wname+".pmesh").c_str());
  int VOL_disp = COM_get_dataitem_handle((wname+".disp").c_str());
  int VOL_nc = COM_get_dataitem_handle((wname+".nc").c_str());
  int VOL_all = COM_get_dataitem_handle((wname+".all").c_str());
  string MyString = "00";

  niter++;
  for (int i = 1; i < niter; i++){
    COM_call_function( MOP_smth, &VOL_pmesh, &VOL_disp);
    string::size_type len;
    len = MyString.length();
    MyString.erase(0,len+1);
    MyStream.str("");
    MyStream << "00" << i;
    MyString = MyStream.str();

    // ADD DISPS ONTO NC
    COM_call_function( BLAS_add,&VOL_disp,&VOL_nc,&VOL_nc);  

    // CALCULATE ASPECT RATIOS
    //COM_call_function(MOP_asp, &VOL_all);
  }

  // PRINT SMOOTHED WINDOW
  COM_call_function( OUT_write, 
		     outputfname.c_str(),   //filename_pre 
		     &VOL_all,           // attribute
		     wname.c_str(),      // material
		     MyString.c_str()    // time
                   ); 

  COM_print_profile("", "");

  COM_finalize();
  MPI_Finalize();
}






