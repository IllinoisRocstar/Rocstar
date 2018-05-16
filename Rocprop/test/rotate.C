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
// $Id: rotate.C,v 1.8 2008/12/06 08:45:28 mtcampbe Exp $

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

COM_EXTERN_MODULE( Simpal);
COM_EXTERN_MODULE( SurfMap);
COM_EXTERN_MODULE( Rocprop);
COM_EXTERN_MODULE( SurfUtil);
COM_EXTERN_MODULE( SimOut);

void load_modules() {
  COM_LOAD_MODULE_STATIC_DYNAMIC(Simpal, "BLAS");
  COM_LOAD_MODULE_STATIC_DYNAMIC(SurfMap, "MAP");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocprop, "PROP");
  COM_LOAD_MODULE_STATIC_DYNAMIC(SurfUtil, "SURF");
  COM_LOAD_MODULE_STATIC_DYNAMIC(SimOut, "OUT");
}

static int rank = 0;

void print_usage( int argc, char *argv[]) {
  if ( argc <= 2) {
    cout << "Usage: " << argv[0] << " <surffile> <controlfile>" << std::endl;
    
    exit(-1);
  }
}

struct Control_parameter {
  Control_parameter() : perturb(0), speed(0), sploc('n'), timestep(0), 
			steps(0), interval(0) {}
  
  string method;
  string wavefrontal;
  string normaldif;
  string eigthres;
  string courant;
  string fangle;
  string smoother;
  string rediter;

  double perturb;
  double speed;
  char   sploc;
  double timestep;
  int    steps;
  int    interval;
  string verbose;
};

void read_control_file( const char *fname, Control_parameter &cp) {
  /* Sample control file:
   * method: fo          # method: fo and mp
   * wavefrontal: 1        # wavefrontal condition
   * normaldif: 1        # normal diffusion
   * eigthres: 1.e-4     # threshold for null space: 0..1 (1.e-4 default)
   * courant: 0.5        # courant constant
   * fangle: 35          # feature edge angle: between 0 and 180
   * smoother: angle     # type of mesh-smoothing algorithm
   * rediter: 1          # Number of iterations for vertex redistribution
   * speed: 0.1          # Speed
   * sploc: e            # location of speed: n or e
   * timestep: 0.001     # time step
   * steps: 100          # number of time steps
   * interval: 10        # output intervals
   * verbose: 1          # verbose level
   */
  ifstream is(fname); COM_assertion_msg( is, "File does not exist");

  while ( !is.eof()) {
    char buf[255];
    is.get( buf, 255, ':');
    if ( buf[0] == '\0') { is.getline( buf, 255); continue; }

    istringstream istr(buf); 
    string keyword; istr >> keyword;
    is.getline( buf, 255, ':'); 

    if ( keyword == "method")
      is >> cp.method;
    else if ( keyword == "wavefrontal")
      is >> cp.wavefrontal;
    else if ( keyword == "normaldif")
      is >> cp.normaldif;
    else if ( keyword == "eigthres")
      is >> cp.eigthres;
    else if ( keyword == "courant")
      is >> cp.courant;
    else if ( keyword == "fangle")
      is >> cp.fangle;
    else if ( keyword == "smoother")
      is >> cp.smoother;
    else if ( keyword == "rediter")
      is >> cp.rediter;
    else if ( keyword == "perturb")
      is >> cp.perturb;
    else if ( keyword == "speed")
      is >> cp.speed;
    else if ( keyword == "sploc")
      is >> cp.sploc;
    else if ( keyword == "timestep")
      is >> cp.timestep;
    else if ( keyword == "steps")
      is >> cp.steps;
    else if ( keyword == "interval")
      is >> cp.interval;
    else if ( keyword == "verbose")
      is >> cp.verbose;
    else
      std::cerr << "Unknow keyword " << keyword << std::endl;
    is.getline( buf, 255); 
  }

  if ( rank==0) std::cout << " speed is " << cp.speed << std::endl;
}

void init_parameters( const Control_parameter &cntr_param) {
  int PROP_set_option = COM_get_function_handle( "PROP.set_option");

  if ( !cntr_param.method.empty()) {
    COM_call_function( PROP_set_option, "method", cntr_param.method.c_str());

    if ( rank==0) std::cout << "Set propagation method to " << cntr_param.method << std::endl;
  }

  if ( !cntr_param.wavefrontal.empty()) {
    COM_call_function( PROP_set_option, "wavefrontal", cntr_param.wavefrontal.c_str());
    if ( rank==0) std::cout << "Set wavefrontal to " << cntr_param.wavefrontal << std::endl;
  }

  if ( !cntr_param.normaldif.empty()) {
    COM_call_function( PROP_set_option, "normaldif", cntr_param.normaldif.c_str());
    if ( rank==0) std::cout << "Set normaldif to " << cntr_param.normaldif << std::endl;
  }

  if ( !cntr_param.eigthres.empty()) {
    COM_call_function( PROP_set_option, "eigthres", cntr_param.eigthres.c_str());
    if ( rank==0) std::cout << "Set eigthres to " << cntr_param.eigthres << std::endl;
  }

  if ( !cntr_param.courant.empty()) {
    COM_call_function( PROP_set_option, "courant", cntr_param.courant.c_str());
    if ( rank==0) std::cout << "Set courant constant to " << cntr_param.courant << std::endl;
  }

  if ( !cntr_param.fangle.empty()) {
    COM_call_function( PROP_set_option, "fangle", cntr_param.fangle.c_str());
    if ( rank==0) std::cout << "Set feature angle to " << cntr_param.fangle << std::endl;
  }

  if ( !cntr_param.smoother.empty()) {
    COM_call_function( PROP_set_option, "smoother", cntr_param.smoother.c_str());
    if ( rank==0) std::cout << "Set smoother to " << cntr_param.smoother << std::endl;
  }
  
  if ( !cntr_param.rediter.empty()) {
    COM_call_function( PROP_set_option, "rediter", cntr_param.rediter.c_str());
    if ( rank==0) std::cout << "Set rediter to " << cntr_param.rediter << std::endl;
  }

  if ( !cntr_param.verbose.empty()) {
    COM_call_function( PROP_set_option, "verbose", cntr_param.verbose.c_str());

    if ( rank==0) std::cout << "Set verbose level to " << cntr_param.verbose << std::endl;
  }

}

// Read in a surface pmesh, and return its window name.
std::string read_in_mesh ( const char *fname) {
  if ( rank==0) cout << "Reading surface mesh file \"" << fname << '"' << endl;

  std::string fname_str(fname);

  std::string::size_type pos = fname_str.find_first_of( ".");
  const string wname = fname_str.substr( 0, pos);

  if ( rank==0) cout << "Creating window \"" << wname << '"' << endl;

  IM_Reader im_reader;
  int npanes = im_reader.read_winmesh( fname, wname); 
  COM_assertion_msg( npanes>=0, "Failed to read in mesh file. File empty?");

  return wname;
}

// Initialize the constraints for ACM Rocfrac mesh.
void init_constraints_acmfrac( const string &wname) {
  static int BLAS_copy_scalar=0, cnstr_handle, MAP_maxabs;
  
  if ( BLAS_copy_scalar==0) {
    COM_get_function_handle( "BLAS.copy_scalar");
    cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types").c_str());
    MAP_maxabs = COM_get_function_handle( "MAP.reduce_maxabs_on_shared_nodes");
  }

  int zero = 0;
  COM_call_function( BLAS_copy_scalar, &zero, &cnstr_handle);

  // Set for pane 2 separately
  int *cnstr_types;
  COM_get_array( (wname+".cnstr_types").c_str(), 2, (void**)&cnstr_types);
  
  MAP::Vector_3<double> *coors;
  COM_get_array( (wname+".nc").c_str(), 2, (void**)&coors);

  int nnodes;
  COM_get_size( (wname+".nc").c_str(), 2, &nnodes);

  // Loop through all the points
  for ( int j=0; j<nnodes; ++j) {
    cnstr_types[j] = 'x'; // Allow the nodes to move only along x dir.
  }

  COM_call_function( MAP_maxabs, &cnstr_handle);

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}

// Initialize the constraints for ACM Rocflu mesh.
void init_constraints_acmflu( const string &wname) {
  static int BLAS_copy_scalar=0, cnstr_handle, MAP_maxabs;
  
  if ( BLAS_copy_scalar==0) {
    COM_get_function_handle( "BLAS.copy_scalar");
    cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types").c_str());
    MAP_maxabs = COM_get_function_handle( "MAP.reduce_maxabs_on_shared_nodes");
  }

  int zero = 0;
  COM_call_function( BLAS_copy_scalar, &zero, &cnstr_handle);

  // Set for pane 3 separately
  int *cnstr_types;
  COM_get_array( (wname+".cnstr_types").c_str(), 3, (void**)&cnstr_types);
  
  MAP::Vector_3<double> *coors;
  COM_get_array( (wname+".nc").c_str(), 3, (void**)&coors);

  int nnodes;
  COM_get_size( (wname+".nc").c_str(), 3, &nnodes);

  // Loop through all the points
  for ( int j=0; j<nnodes; ++j) {
    if ( coors[j].x() <= 1.e-10 || coors[j].x() > 0.047624)
      cnstr_types[j] = 2; 	// Fix the nodes at two ends.
    else 
      cnstr_types[j] = 'x'; // Allow the nodes to move only along x dir.
  }

  COM_call_function( MAP_maxabs, &cnstr_handle);

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}

// Initialize the constraints for Star slice.
void init_constraints_starslice( const string &wname) {
  static int BLAS_copy_scalar=0, cnstr_handle, MAP_maxabs;
  
  if ( BLAS_copy_scalar==0) {
    COM_get_function_handle( "BLAS.copy_scalar");
    cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types").c_str());
    MAP_maxabs = COM_get_function_handle( "MAP.reduce_maxabs_on_shared_nodes");
  }
  int zero = 0;
  COM_call_function( BLAS_copy_scalar, &zero, &cnstr_handle);

  // Set for pane 2 separately
  int *cnstr_types;
  COM_get_array( (wname+".cnstr_types").c_str(), 2, (void**)&cnstr_types);
  
  MAP::Vector_3<double> *coors;
  COM_get_array( (wname+".nc").c_str(), 2, (void**)&coors);

  int nnodes;
  COM_get_size( (wname+".nc").c_str(), 2, &nnodes);

  // Loop through all the points
  for ( int j=0; j<nnodes; ++j) {
    if ( coors[j].y() > -20.054 || coors[j].x() < -20.5538)
      cnstr_types[j] = -'y';  // Move within xz plane.
    else
      cnstr_types[j] = 2;     // Fix other points
  }

  COM_call_function( MAP_maxabs, &cnstr_handle);

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}

void init_attributes( const string &wname, 
		      const Control_parameter &cntr_param) {
  COM_new_dataitem((wname+".flag").c_str(), 'p', COM_DOUBLE, 1, "");
  COM_set_size( (wname+".flag").c_str(), 0, 1);

  COM_new_dataitem((wname+".disps_total").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".cnstr_types").c_str(), 'n', COM_INT, 1, "");
  
  COM_new_dataitem((wname+".spds").c_str(), 'n', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem((wname+".disps").c_str(), 'n', COM_DOUBLE, 3, "m");

  COM_new_dataitem((wname+".disps_novis").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".facenormals").c_str(), 'e', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".faceheights").c_str(), 'e', COM_DOUBLE, 1, "");

  // Dataitem for storing the number of eigenvalues for each node.
  COM_new_dataitem((wname+".lambdas").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".eigvecs").c_str(), 'n', COM_DOUBLE, 9, "");
  COM_new_dataitem((wname+".tangranks").c_str(), 'n', COM_CHAR, 1, "");
  COM_new_dataitem((wname+".scales").c_str(), 'n', COM_DOUBLE, 1, "");

  COM_resize_array( (wname+".data").c_str());
  COM_window_init_done( wname.c_str());
}

inline double square(double x) { return x*x; }

// Initialize speed for n points.
// Coordinates and velocities are stored in two nx3 arrays, respectively.
void init_speed_arrays( const double *coors, double *velocities, int n) {
  const double pi = 3.14159265358979;
  for ( int i=0; i<n; ++i) {
#if 1
    velocities[3*i] = pi*(-coors[3*i+1]);
    velocities[3*i+1] = pi*(coors[3*i]-1.5);
    velocities[3*i+2] = 0;
#else
    double xs[]={pi*coors[3*i], pi*coors[3*i+1], pi*coors[3*i+2]};
    velocities[3*i] = 2*square(sin(xs[0]))* sin(2*xs[1])*sin(2*xs[2]);
    velocities[3*i+1] = -sin(2*xs[0])* square(sin(xs[1]))* sin(2*xs[2]);
    velocities[3*i+2] = -sin(2*xs[0])* sin(2*xs[1])* square(sin(xs[2]));
#endif
  }
}

void init_speed( const string &wname) {
  // Set speed for different panes.
  int *pane_ids, npanes;
  COM_get_panes( wname.c_str(), &npanes, &pane_ids);
  
  for ( int i=0; i<npanes; ++i) {
    double *coor;
    double *spds;
    COM_get_array( (wname+".nc").c_str(), pane_ids[i], &coor);
    COM_get_array( (wname+".spds").c_str(), pane_ids[i], &spds);

    int n;
    COM_get_size( (wname+".nc").c_str(), pane_ids[i], &n);
    init_speed_arrays( coor, spds, n);
  }

  COM_free_buffer( &pane_ids);
}

void output_solution( const string &wname, const char *timelevel) {
  static int OUT_write = 0, hdl;

  if ( OUT_write==0) {
    OUT_write = COM_get_function_handle( "OUT.write_dataitem");
    hdl = COM_get_dataitem_handle( (wname+".all").c_str());

    int OUT_set = COM_get_function_handle( "OUT.set_option");
    COM_call_function( OUT_set, "format", "HDF");
  }

  std::string fname = wname+"_"+timelevel;

  if ( !COMMPI_Initialized()) fname.append(".hdf");
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
  read_control_file( argc>2?argv[2]:"control.txt", cntr_param);
  
  // Initialize the attributes and output the initial solution.
  init_attributes( wname, cntr_param);

  // Attributes
  int pmesh = COM_get_dataitem_handle( (wname+".pmesh").c_str());
  int nc = COM_get_dataitem_handle( (wname+".nc").c_str());
  int spds = COM_get_dataitem_handle( (wname+".spds").c_str());
  int disps = COM_get_dataitem_handle( (wname+".disps").c_str());
  int disps_total = COM_get_dataitem_handle( (wname+".disps_total").c_str());

  // Funcitions
  int PROP_init = COM_get_function_handle( "PROP.initialize");
  int PROP_perturb = COM_get_function_handle( "PROP.perturb_mesh");
  int PROP_propagate = COM_get_function_handle( "PROP.propagate");
  int BLAS_add = COM_get_function_handle( "BLAS.add");

  if ( rank==0) cout << "Propagating interface..." << endl;
  // Initialize control parameters
  COM_call_function( PROP_init, &pmesh);
  init_parameters( cntr_param);

  if ( wname.substr(0,6) == "acmflu") 
    init_constraints_acmflu( wname);
  else if ( wname.substr(0,3) == "acm") 
    init_constraints_acmfrac( wname);
  else if ( wname.substr(0,4) == "star") 
    ; // init_constraints_starslice( wname);

  if ( cntr_param.perturb > 0)
    COM_call_function( PROP_perturb, &pmesh, &cntr_param.perturb);

  output_solution( wname, "00000");

  COM_set_profiling(1);
  COM_set_profiling_barrier( PROP_propagate, MPI_COMM_WORLD);

  char fname[40];
  std::sprintf(fname, "timedata_%03d.txt", rank);

  for ( int i=1; i<=cntr_param.steps; ++i) {
    if ( rank==0) cout << "Step " << i << endl;
    
    init_speed( wname);
    COM_call_function( PROP_propagate, &pmesh, &spds, 
		       &cntr_param.timestep, &disps);

    COM_call_function( BLAS_add, &nc, &disps, &nc);
    COM_call_function( BLAS_add, &disps_total, &disps, &disps_total);

    if ( i%cntr_param.interval == 0) {
      char   steps[10];
      std::sprintf( steps, "%05d", i);
      output_solution( wname, steps);

      COM_print_profile( fname, "Proptest");
    }
  }

  COM_finalize();
}






