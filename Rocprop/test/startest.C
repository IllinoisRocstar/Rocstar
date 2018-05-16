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
// $Id: startest.C,v 1.7 2008/12/06 08:45:28 mtcampbe Exp $

#include "roccom.h"
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
#include "roccom_assertion.h"

#include "../Rocsurf/test/IM_Reader.h"
#include "Remesher_Simmetrix.h"

using namespace std;

COM_EXTERN_MODULE( Rocblas);
COM_EXTERN_MODULE( Rocsurf);
COM_EXTERN_MODULE( Rocprop);
COM_EXTERN_MODULE( Rocout);

void load_modules() {
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocblas, "BLAS");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocprop, "PROP");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocsurf, "SURF");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocout, "OUT");
}

static int rank = 0;

void print_usage( int argc, char *argv[]) {
  if ( argc <= 2) {
    cout << "Usage: " << argv[0] << " <surffile> <controlfile>" << std::endl;
    
    exit(-1);
  }
}

struct Control_parameter {
  Control_parameter() : speed(0), sploc('n'), timestep(0), 
			steps(0), interval(0), remesh_interval(0) {}
  
  string method;
  string wavefrontal;
  string normaldif;
  string eigthres;
  string courant;
  string fangle;
  double speed;
  char   sploc;
  double timestep;
  int    steps;
  int    interval;
  int    remesh_interval;
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
   * speed: 0.1          # Speed
   * sploc: e            # location of speed: n or e
   * timestep: 0.001     # time step
   * steps: 100          # number of time steps
   * interval: 10        # output intervals
   * remesh_interval: 50 # remesh frequency
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
    else if ( keyword == "remesh_interval")
      is >> cp.remesh_interval;
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
  int npanes = im_reader.read_winmesh( fname, wname, false); 
  COM_assertion_msg( npanes>=0, "Failed to read in mesh file. File empty?");

  return wname;
}

void init_attributes( const string &wname, 
		      const Control_parameter &cntr_param) {
  COM_new_attribute((wname+".cnstr_types_facial").c_str(), 'e', COM_INT, 1, "");
  COM_new_attribute((wname+".spds").c_str(), cntr_param.sploc, 
		    COM_DOUBLE, 1, "m/s");
  COM_new_attribute((wname+".disps").c_str(), 'n', COM_DOUBLE, 3, "m");

  COM_new_attribute((wname+".facenormals").c_str(), 'e', COM_DOUBLE, 3, "");
  COM_new_attribute((wname+".facecenters").c_str(), 'e', COM_DOUBLE, 3, "");

  // Attribute for storing the number of eigenvalues for each node.
  COM_new_attribute((wname+".tangranks").c_str(), 'n', COM_CHAR, 1, "");

  COM_resize_array( (wname+".data").c_str());
  COM_window_init_done( wname.c_str());

  int spds = COM_get_attribute_handle( (wname+".spds").c_str());
  int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
  COM_call_function( BLAS_copy_scalar, &cntr_param.speed, &spds);
}

void output_solution( const string &wname, const string &material,
		      const char *timelevel, const char *attr=NULL) {
  int OUT_write = COM_get_function_handle( "OUT.write_attribute");
  int hdl;

  if (attr==NULL) 
    hdl = COM_get_attribute_handle( (wname+".all").c_str());
  else
    hdl = COM_get_attribute_handle( (wname+"."+attr).c_str());

  int OUT_set = COM_get_function_handle( "OUT.set_option");
  COM_call_function( OUT_set, "format", "HDF");

  std::string fname = material;
  if ( attr) fname = fname + "-"+attr;
  fname = fname +"_"+timelevel;

  if ( !COMMPI_Initialized()) fname.append(".hdf");
  else fname.append("_");

  COM_call_function( OUT_write, (char*)fname.c_str(),
		     &hdl, (char*)material.c_str(), timelevel);
}

template <class T>
T square(T t) { return t*t; }

// Initialize the constraints for Star slice.
void init_constraints_starslice( const string &wname) {
  int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
  int cnstr_handle = COM_get_attribute_handle( (wname+".cnstr_types_facial").c_str());
  int zero = 0;
  COM_call_function( BLAS_copy_scalar, &zero, &cnstr_handle);

  int SURF_normals = COM_get_function_handle( "SURF.compute_element_normals");
  int nrms_handle = COM_get_attribute_handle( (wname+".facenormals").c_str());
  COM_call_function( SURF_normals, &nrms_handle);

  int SURF_centers = COM_get_function_handle( "SURF.interpolate_to_centers");
  int nc_handle = COM_get_attribute_handle( (wname+".nc").c_str());
  int centers_handle = COM_get_attribute_handle( (wname+".facecenters").c_str());
  COM_call_function( SURF_centers, &nc_handle, &centers_handle);

  // Set speed for different panes.
  int *pane_ids, npanes;
  COM_get_panes( wname.c_str(), &npanes, &pane_ids);

  for ( int i=0; i<npanes; ++i) {
    int pid = pane_ids[i];
    MAP::Vector_3<double> *nrms, *centers;
    COM_get_array( (wname+".facenormals").c_str(), pid, (void**)&nrms);
    COM_get_array( (wname+".facecenters").c_str(), pid, (void**)&centers);

    // Set for pane 2 separately
    int *cnstr_types;
    COM_get_array( (wname+".cnstr_types_facial").c_str(), 
		   pid, (void**)&cnstr_types);
    
    int nelems;
    COM_get_size( (wname+".facenormals").c_str(), pid, &nelems);

    double eps = 1.e-2;
    // Loop through all faces
    for ( int j=0; j<nelems; ++j) {
      if ( std::abs(std::abs(nrms[j][0])-1)<eps && centers[j][0] < 0.001)
	cnstr_types[j] = 't';  // Move within yz plane.
      else if ( std::abs(std::abs(nrms[j][1])-1)<eps && 
		(centers[j][1]>-20.054 || centers[j][1] < -20.5538))
	cnstr_types[j] = 't';  // Move within xz plane.
      else if ( square(centers[j][0])+square(centers[j][2])>=3.9)
	cnstr_types[j] = 2;
    }   
  }
  COM_free_buffer( &pane_ids);

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}

void remesh( const string &wname_old, const string &wname,
	     const Control_parameter &cntr_param, 
	     double len, bool preserve_feature) {
  COM_new_window( wname.c_str());
  COM_window_init_done( wname.c_str());

  static PROP::Remesher_Simmetrix *rem=NULL;
  if ( rem==NULL) {
    rem = new PROP::Remesher_Simmetrix();
    int PROP_set_remesher = COM_get_function_handle( "PROP.set_remesher");
    int owner = 1;
    COM_call_function( PROP_set_remesher, rem, &owner);
  }
  
  int PROP_remesh = COM_get_function_handle( "PROP.remesh_serial");
  int newmesh_hdl = COM_get_attribute_handle( (wname+".mesh").c_str());

  double fangle = (preserve_feature&&!cntr_param.fangle.empty()) 
    ? std::atof(cntr_param.fangle.c_str()) : 0.;
  COM_call_function( PROP_remesh, &newmesh_hdl, &len, &fangle);

  // Initialize the attributes and output the initial solution.
  init_attributes( wname, cntr_param);
  init_constraints_starslice( wname);

  int PROP_init = COM_get_function_handle( "PROP.initialize");
  int pmesh = COM_get_attribute_handle( (wname+".pmesh").c_str());
  COM_call_function( PROP_init, &pmesh);

  // Delete the old window
  COM_delete_window( wname_old.c_str());
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
  int pmesh = COM_get_attribute_handle( (wname+".pmesh").c_str());
  int nc = COM_get_attribute_handle( (wname+".nc").c_str());
  int spds = COM_get_attribute_handle( (wname+".spds").c_str());
  int disps = COM_get_attribute_handle( (wname+".disps").c_str());

  // Funcitions
  int PROP_init = COM_get_function_handle( "PROP.initialize");
  int PROP_propagate = COM_get_function_handle( "PROP.propagate");
  int BLAS_add = COM_get_function_handle( "BLAS.add");

  if ( rank==0) cout << "Propagating interface..." << endl;
  // Initialize control parameters
  COM_call_function( PROP_init, &pmesh);
  init_parameters( cntr_param);

  // Compute average edge length
  double len;
  int PROP_length = COM_get_function_handle( "PROP.compute_edge_lengths");
  COM_call_function( PROP_length, &len);

  init_constraints_starslice( wname);

  COM_set_profiling(1);
  COM_set_profiling_barrier( PROP_propagate, MPI_COMM_WORLD);

  char fname[40];
  std::sprintf(fname, "timedata_%03d.txt", rank);

  std::string wname_rem = wname;
  for ( int i=1; i<=cntr_param.steps; ++i) {
    if ( rank==0) cout << "Step " << i << endl;

    for ( int j=0;j<2;++j) {
      // Perform remeshing
      if (  j || cntr_param.remesh_interval>0 && 
	   (i-1)%cntr_param.remesh_interval==0) {
        std::string wname_old = wname_rem;
	char buf[40];
	std::sprintf(buf, "%s%d", wname.c_str(), i);

	wname_rem = buf;
	remesh( wname_old, wname_rem, cntr_param, len, j==0);

	pmesh = COM_get_attribute_handle( (wname_rem+".pmesh").c_str());
	nc = COM_get_attribute_handle( (wname_rem+".nc").c_str());
	spds = COM_get_attribute_handle( (wname_rem+".spds").c_str());
	disps = COM_get_attribute_handle( (wname_rem+".disps").c_str());

	// Write out the new mesh
	char   steps[10];
	if ( cntr_param.steps>=1.e6)
	  std::sprintf( steps, "%07d", i);
	else if ( cntr_param.steps>=1.e5)
	  std::sprintf( steps, "%06d", i);
	else
	  std::sprintf( steps, "%05d", i);

	output_solution( wname_rem, wname, steps, "mesh");
      }

      int code;
      COM_call_function( PROP_propagate, &pmesh, &spds, 
			 &cntr_param.timestep, &disps, &code);

      if ( code) {
	COM_assertion_msg( j==0, "failed");
	continue;
      }
      else
	break;
    }

    COM_call_function( BLAS_add, &nc, &disps, &nc);

    if ( i%cntr_param.interval == 0) {
      char   steps[10];
      if ( cntr_param.steps>=1.e6)
	std::sprintf( steps, "%07d", i);
      else if ( cntr_param.steps>=1.e5)
	std::sprintf( steps, "%06d", i);
      else
	std::sprintf( steps, "%05d", i);

      output_solution( wname_rem, wname, steps);

      COM_print_profile( fname, "Proptest");
    }
  }

  COM_finalize();
}






