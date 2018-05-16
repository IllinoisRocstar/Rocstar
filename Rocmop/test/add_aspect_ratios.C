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
// $Id: add_aspect_ratios.C,v 1.3 2008/12/06 08:45:25 mtcampbe Exp $

#include "com.h"
#include "mapbasic.h"
#include "Rocblas.h"
#include <iostream>
#include <sstream>
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

  if ( argc < 3) {
    std::cout << "Usage: " << argv[0] 
	      << " <Rocin input file> <material name>" << endl;
    exit(-1);
  }
  string inputfname(argv[1]);
  const string matname(argv[2]);

  MPI_Comm comm = MPI_COMM_WORLD;
    
  COM_init( &argc, &argv);

  // LOAD MODULES
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimOUT, "OUT");
  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocmop, "MOP");

  // GET FUNCTION HANDLES
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int IN_read   = COM_get_function_handle( "IN.read_by_control_file");
  int MOP_asp   = COM_get_function_handle( "MOP.add_aspect_ratios");

  // SET OUTPUT LEVELS
  COM_set_verbose(11);
  COM_set_profiling(1);

  const string newfname("with_aspect_ratios");

  // READ IN FILES
  COM_call_function( IN_read,
		     inputfname.c_str(), 
		     matname.c_str(),
		     &comm);

  int MAT_all = COM_get_dataitem_handle((matname+".all").c_str());
  string MyString = "001";

  // CALCULATE ASPECT RATIOS
  COM_call_function(MOP_asp, &MAT_all);

 // PRINT SMOOTHED WINDOW
  COM_call_function( OUT_write, 
		     newfname.c_str(),   //filename_pre 
 		     &MAT_all,           // attribute
		     matname.c_str(),      // material
		     "001"               // time
		     ); 

  COM_print_profile("", "");

  COM_finalize();
  MPI_Finalize();
}






