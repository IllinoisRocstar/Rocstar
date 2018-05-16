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
// $Id: detectfeas.C,v 1.4 2008/12/06 08:45:28 mtcampbe Exp $

#include "com.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <sstream>
#include "com_assertion.h"
#include "IM_Reader.h"


using namespace std;

COM_EXTERN_MODULE( Rocprop);
COM_EXTERN_MODULE( Rocout);

void load_modules() {
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocprop, "PROP");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocout, "OUT");
}

static int rank = 0;

void print_usage( int argc, char *argv[]) {
  if ( argc < 2) {
    cout << "Usage: " << argv[0] << " <surffile> [<controlfile>]" << std::endl;
    
    exit(-1);
  }
}

struct Control_parameter {
  Control_parameter() {}
  
  string fangle;
  string sangle;
  string format;
};

void read_control_file( const char *fname, Control_parameter &cp) {
  /* Sample control file:
   * fangle: 15          # weak feature edge angle: between 0 and 180
   * sangle: 60          # strong feature edge angle: between 0 and 180
   * format: CGNS        # file format: CGNS or HDF
   */
  ifstream is(fname); COM_assertion_msg( is, "File does not exist");

  while ( !is.eof()) {
    char buf[255];
    is.get( buf, 255, ':');
    if ( buf[0] == '\0') { is.getline( buf, 255); continue; }

    istringstream istr(buf); 
    string keyword; istr >> keyword;
    is.getline( buf, 255, ':'); 

    if ( keyword == "fangle")
      is >> cp.fangle;
    else if ( keyword == "sangle")
      is >> cp.sangle;
    else if ( keyword == "format")
      is >> cp.format;
    else
      std::cerr << "Unknow keyword " << keyword << std::endl;
    is.getline( buf, 255); 
  }
}

void init_parameters( const Control_parameter &cntr_param) {
  int PROP_set_option = COM_get_function_handle( "PROP.set_option");

  if ( !cntr_param.fangle.empty()) {
    COM_call_function( PROP_set_option, "fangle", cntr_param.fangle.c_str());
    if ( rank==0) std::cout << "Set weak-feature angle to " << cntr_param.fangle << std::endl;
  }
  else
    COM_call_function( PROP_set_option, "fangle", "15");

  if ( !cntr_param.sangle.empty()) {
    COM_call_function( PROP_set_option, "sangle", cntr_param.sangle.c_str());
    if ( rank==0) std::cout << "Set strong-feature angle to " << cntr_param.sangle << std::endl;
  }
  else
    COM_call_function( PROP_set_option, "sangle", "60");

  COM_call_function( PROP_set_option, "rediter", "0");
  COM_call_function( PROP_set_option, "weight", "3"); // Angle-based
}

// Read in a surface pmesh, and return its window name.
std::string read_in_mesh ( const char *fname) {
  if ( rank==0) cout << "Reading surface mesh file \"" << fname << '"' << endl;

  std::string fname_str(fname);

  std::string::size_type pos = fname_str.find_first_of( ".");
  const string wname = fname_str.substr( 0, pos);

  if ( rank==0) cout << "Creating window \"" << wname << '"' << endl;

  IM_Reader im_reader;
  int npanes = im_reader.read_winmesh( fname, wname, false); 
  COM_assertion_msg( npanes>=0, "Failed to read in mesh file. File empty?");

  return wname;
}

void init_attributes( const string &wname, 
		      const Control_parameter &cntr_param) {
  COM_new_dataitem((wname+".disps").c_str(), 'n', COM_DOUBLE, 3, "m");

  // Ridges of each pane.
  COM_new_dataitem((wname+".ridges").c_str(), 'p', COM_INT, 2, "");
  COM_set_size((wname+".ridges").c_str(), 0, 0);

  COM_new_dataitem((wname+".tangranks").c_str(), 'n', COM_INT, 1, "");
  // COM_new_dataitem((wname+".weaks").c_str(), 'n', COM_INT, 1, "");

  COM_resize_array( (wname+".data").c_str());
  COM_window_init_done( wname.c_str());
}

void output_solution( const string &wname, const char *timelevel, 
		      const string &format) {
  static int OUT_write = 0, hdl;

  if ( OUT_write==0) {
    OUT_write = COM_get_function_handle( "OUT.write_dataitem");
    hdl = COM_get_dataitem_handle( (wname+".all").c_str());

    int OUT_set = COM_get_function_handle( "OUT.set_option");
#if USE_CGNS
    if ( format != "HDF" && format != "hdf")
      COM_call_function( OUT_set, "format", "CGNS");
    else
#endif
      COM_call_function( OUT_set, "format", "HDF");
  }

  std::string fname = wname+"_"+timelevel;

  if ( !COMMPI_Initialized()) {
#if USE_CGNS
    if ( format != "HDF" && format != "hdf")
      fname.append(".cgns");
    else
#endif
      fname.append(".hdf");
  }
  else fname.append("_");

  COM_call_function( OUT_write, (char*)fname.c_str(),
		     &hdl, (char*)wname.c_str(), timelevel);
}

int main(int argc, char *argv[]) {
  COM_init( &argc, &argv);
  load_modules();
  print_usage( argc, argv);

  if ( COMMPI_Initialized()) rank = COMMPI_Comm_rank( MPI_COMM_WORLD);

  // Read in mesh file.
  string wname = read_in_mesh( argv[1]);

  // Read in control parameters
  Control_parameter cntr_param;
  if ( argc>2) read_control_file( argv[2], cntr_param);

  // Initialize the attributes.
  init_attributes( wname, cntr_param);
  
  // Attributes
  int pmesh = COM_get_dataitem_handle( (wname+".pmesh").c_str());
  int disps = COM_get_dataitem_handle( (wname+".disps").c_str());

  // Funcitions
  int PROP_init = COM_get_function_handle( "PROP.initialize");
  int PROP_propagate = COM_get_function_handle( "PROP.propagate");
  init_parameters( cntr_param);

  if ( rank==0) cout << "Initializing data structures..." << endl;
  // Initialize control parameters
  COM_call_function( PROP_init, &pmesh);

  if ( rank==0) cout << "Detecting features..." << endl;
  double zero=0.;
  // Perform one step of smoothing to detect features.
  COM_call_function( PROP_propagate, &pmesh, &disps, &zero, &disps, &zero);

  if ( rank==0) cout << "Writing out solutions..." << endl;
  output_solution( wname, "00000", cntr_param.format);

  COM_finalize();
}






