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
// $Id: proptest.C,v 1.33 2008/12/06 08:45:28 mtcampbe Exp $

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
#include "Rocblas.h"
#include "Rocsurf.h"

#include "IM_Reader.h"

#ifdef MESH_ADAPT
#include "AdaptCOMWindow.h"
#endif


using namespace std;

COM_EXTERN_MODULE( Rocblas);
COM_EXTERN_MODULE( Rocmap);
COM_EXTERN_MODULE( Rocprop);
COM_EXTERN_MODULE( Rocsurf);
COM_EXTERN_MODULE( Rocout);

void load_modules() {
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocblas, "BLAS");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocmap, "MAP");
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
  Control_parameter() : perturb(0), speed(0), sploc('n'), timestep(0), 
			steps(0), start(0), interval(0), remesh_interval(10), 
			subadapt(0.5), adapt_iter(2), do_redist(1), do_collapse(1), 
			do_flip(1), do_split(1), collapse_ratio(0), 
			split_angle(0), refine(1), edge_len(0) {}
  
  string method;
  string wavefrontal;
  string normaldif;
  string eigthres;
  string courant;
  string fangle;
  string sangle;
  string smoother;
  string rediter;
  string feature_layer;
  string weight;
  string format;

  double perturb;
  double speed;
  char   sploc;
  double timestep;
  int    steps;
  int    start;
  int    interval;
  int    remesh_interval;
  string verbose;

  // The following options are for mesh adaptation
  double subadapt;
  int    adapt_iter;
  int    do_redist;
  int    do_collapse;
  int    do_flip;
  int    do_split;
  double collapse_ratio;
  double split_angle;
  double refine;
  double edge_len;
};

void read_control_file( const char *fname, Control_parameter &cp) {
  /* Sample control file:
   * method: fo          # method: fo and mp
   * wavefrontal: 1      # wavefrontal expansion
   * normaldif: 1        # normal diffusion
   * eigthres: 1.e-4     # threshold for null space: 0..1 (1.e-4 default)
   * courant: 0.5        # courant constant
   * fangle: 25          # feature edge angle: between 0 and 180
   * sangle: 65          # feature edge angle: between 0 and 180
   * smoother: angle     # type of mesh-smoothing algorithm
   * rediter: 1          # Number of iterations for vertex redistribution
   * speed: 0.1          # Speed
   * sploc: e            # location of speed: n or e
   * timestep: 0.001     # time step
   * steps: 100          # number of time steps
   * interval: 10        # output intervals
   * remesh_interval: 10 # remesh intervals
   * weight: 3           # angle weighted
   * verbose: 1          # verbose level
   *
   * adapt_iter: 2
   * do_flip: 1
   * do_split: 1
   * do_redist: 1
   * do_collapse: 1
   * collapse_ratio: 0.1
   * split_angle: 2.618  # (150 degrees)
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
    else if ( keyword == "sangle")
      is >> cp.sangle;
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
    else if ( keyword == "start")
      is >> cp.start;
    else if ( keyword == "interval")
      is >> cp.interval;
    else if ( keyword == "remesh_interval")
      is >> cp.remesh_interval;
    else if ( keyword == "weight")
      is >> cp.weight;
    else if ( keyword == "format")
      is >> cp.format;
    else if ( keyword == "verbose")
      is >> cp.verbose;
    else if ( keyword == "verbose")
      is >> cp.verbose;
    else if ( keyword == "subadapt")
      is >> cp.subadapt;
    else if ( keyword == "adapt_iter")
      is >> cp.adapt_iter;
    else if ( keyword == "do_redist")
      is >> cp.do_redist;
    else if ( keyword == "do_contract" || keyword == "do_collapse")
      is >> cp.do_collapse;
    else if ( keyword == "do_flip")
      is >> cp.do_flip;
    else if ( keyword == "do_split")
      is >> cp.do_split;
    else if ( keyword == "feature_layer")
      is >> cp.feature_layer;
    else if ( keyword == "refine")
      is >> cp.refine;
    else if ( keyword == "edge_len")
      is >> cp.edge_len;
    else if ( keyword == "collapse_ratio")
      is >> cp.collapse_ratio;
    else if ( keyword == "split_angle")
      is >> cp.split_angle;
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
    if ( rank==0) std::cout << "Set weak-feature angle to " << cntr_param.fangle << std::endl;
  }
  if ( !cntr_param.sangle.empty()) {
    COM_call_function( PROP_set_option, "sangle", cntr_param.sangle.c_str());
    if ( rank==0) std::cout << "Set strong-feature angle to " << cntr_param.sangle << std::endl;
  }

  if ( !cntr_param.smoother.empty()) {
    COM_call_function( PROP_set_option, "smoother", cntr_param.smoother.c_str());
    if ( rank==0) std::cout << "Set smoother to " << cntr_param.smoother << std::endl;
  }

  if ( !cntr_param.rediter.empty()) {
    COM_call_function( PROP_set_option, "rediter", cntr_param.rediter.c_str());
    if ( rank==0) std::cout << "Set rediter to " << cntr_param.rediter << std::endl;
  }

  if ( !cntr_param.feature_layer.empty()) {
    COM_call_function( PROP_set_option, "feature_layer", 
		       cntr_param.feature_layer.c_str());
    if ( rank==0) 
      std::cout << "Set feature_layer to " << cntr_param.feature_layer << std::endl;
  }

  if ( !cntr_param.verbose.empty()) {
    COM_call_function( PROP_set_option, "verbose", cntr_param.verbose.c_str());

    if ( rank==0) std::cout << "Set verbose level to " << cntr_param.verbose << std::endl;
  }

  if ( !cntr_param.weight.empty()) {
    COM_call_function( PROP_set_option, "weight", cntr_param.weight.c_str());
    if ( rank==0) std::cout << "Set weight to " << cntr_param.weight << std::endl;
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

void rescale_object( std::string &wname, double alpha, 
		     const SURF::Vector_3<double> &origin) {
  // Rescale and translate the object.
  COM::Window *win=COM_get_com()->get_window_object(wname);
  COM::DataItem *nc=win->dataitem("nc");

  // Rescale the radius to 0.15
  Rocblas::mul_scalar( nc, &alpha, nc);

  // Translate the origin to (0.5,0.75, 0.5);
  COM::DataItem *x=win->dataitem("1-nc");
  COM::DataItem *y=win->dataitem("2-nc");
  COM::DataItem *z=win->dataitem("3-nc");

  Rocblas::add_scalar( x, &origin[0], x);
  Rocblas::add_scalar( y, &origin[1], y);
  Rocblas::add_scalar( z, &origin[2], z);
}

// Initialize the constraints for ACM Rocflu mesh.
void init_constraints_acmflu( const string &wname) {

  int SURF_centers = COM_get_function_handle( "SURF.interpolate_to_centers");
  int nc_handle = COM_get_dataitem_handle( (wname+".nc").c_str());
  int centers_handle = COM_get_dataitem_handle( (wname+".facecenters").c_str());
  COM_call_function( SURF_centers, &nc_handle, &centers_handle);

  // Set speed for different panes.
  int *pane_ids, npanes;
  COM_get_panes( wname.c_str(), &npanes, &pane_ids);
  
  for ( int i=0; i<npanes; ++i) {
    bool not_interacting;
    if ( COM_get_dataitem_handle((wname+".bcflag").c_str())>0) {
      int *bcflag;
      COM_get_array( (wname+".bcflag").c_str(), pane_ids[i], &bcflag);
      
      not_interacting = (*bcflag == 2);
    }
    else
      not_interacting = (pane_ids[i] == 3);
    
    if ( !not_interacting) continue;

    int pid = pane_ids[i], nelems, ngelems;
    COM_get_size( (wname+".cnstr_types_facial").c_str(), pid, &nelems, &ngelems);
    int *cnstr_types;
    COM_get_array( (wname+".cnstr_types_facial").c_str(), pid, (void**)&cnstr_types);    
    double *spds;
    COM_get_array( (wname+".spds").c_str(), pid, (void**)&spds);

    MAP::Vector_3<double> *centers; 
    COM_get_array( (wname+".facecenters").c_str(), pid, (void**)&centers);
   
    // Loop through all faces
    for ( int j=0; j<nelems-ngelems; ++j) {
      spds[j] = 0; 	// Set speed to zero.

      if ( centers[j].x() <= 1.e-10 || centers[j].x() > 0.047624)
	cnstr_types[j] = 2; 	// Fix the nodes at two ends.
      else 
	cnstr_types[j] = 'x'; // Allow nodes to move only along x dir.
    }
  }
  
  COM_free_buffer( &pane_ids);
  
  int cnstr_handle = COM_get_dataitem_handle((wname+".cnstr_types_facial").c_str());
  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}


// Initialize the constraints for ACM Rocfrac mesh.
void init_constraints_acmfrac( const string &wname) {
  // Set constraint for pane 2 so that nodes can move only along x-direction
  int *cnstr_types;
  COM_get_array( (wname+".cnstr_types_panel").c_str(), 2, (void**)&cnstr_types); 
  *cnstr_types = 'x';

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  int cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types_panel").c_str());

  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}

template <class T>
T square(T t) { return t*t; }

// Initialize the constraints for Star slice.
void init_constraints_starslice( const string &wname, 
				 const Control_parameter &cntr_param) {
  int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
  int cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types_facial").c_str());
  int zero = 0;
  COM_call_function( BLAS_copy_scalar, &zero, &cnstr_handle);

  int SURF_normals = COM_get_function_handle( "SURF.compute_element_normals");
  int nrms_handle = COM_get_dataitem_handle( (wname+".facenormals").c_str());
  COM_call_function( SURF_normals, &nrms_handle);

  int SURF_centers = COM_get_function_handle( "SURF.interpolate_to_centers");
  int nc_handle = COM_get_dataitem_handle( (wname+".nc").c_str());
  int centers_handle = COM_get_dataitem_handle( (wname+".facecenters").c_str());
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
    int *cnstr_types, nelems;
    COM_get_array( (wname+".cnstr_types_facial").c_str(), 
		   pid, (void**)&cnstr_types);    
    COM_get_size( (wname+".facenormals").c_str(), pid, &nelems);
    double *spds;
    COM_get_array( (wname+".spds").c_str(), pid, (void**)&spds);

    double eps = 1.e-5;
    // Loop through all faces
    for ( int j=0; j<nelems; ++j) {
      if ( 1-std::abs(nrms[j][2])<eps && 
	   ( centers[j][2] < eps || (1-centers[j][2]) < eps)) {
	cnstr_types[j] = 't';  // Move within xy plane.
	spds[j] = 0;
      }
      else if ( std::abs(nrms[j][2]) < 0.5 && 
		square(centers[j][0])+square(centers[j][1])>=3.999) {
	cnstr_types[j] = 2;
	spds[j] = 0;
      }
      else 
	spds[j] = cntr_param.speed;
    }
  }
  COM_free_buffer( &pane_ids);

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}

// Initialize the constraints for Star slice.
void init_constraints_staraft( const string &wname, 
			       const Control_parameter &cntr_param) {

  int *pane_ids, npanes;
  COM_get_panes( wname.c_str(), &npanes, &pane_ids);
  int cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types_facial").c_str());

  // Set constraints and speeds
  if ( cntr_param.start == 0) {
    int cnstr_handle_pn = COM_get_dataitem_handle((wname+".cnstr_type").c_str());
    if ( cnstr_handle_pn>0) {
      std::cout << "Obtaining constraints from input file" << std::endl;
      int BLAS_copy = COM_get_function_handle( "BLAS.copy");

      COM_call_function( BLAS_copy, &cnstr_handle_pn, &cnstr_handle);

      int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
      double zero = 0;
      int spd_elms = COM_get_dataitem_handle( (wname+".flag_spd").c_str());
      COM_call_function( BLAS_copy_scalar, &zero, &spd_elms);
    }
    else { 
      std::cout << "Did not find constraints from input file. "
		<< "Setting constraints from geometry." << std::endl;

      // Compute constraints from geometry
      int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
      int zero = 0;
      COM_call_function( BLAS_copy_scalar, &zero, &cnstr_handle); 

      int SURF_normals = COM_get_function_handle( "SURF.compute_element_normals");
      int nrms_handle = COM_get_dataitem_handle( (wname+".facenormals").c_str());
      COM_call_function( SURF_normals, &nrms_handle);

      int SURF_centers = COM_get_function_handle( "SURF.interpolate_to_centers");
      int nc_handle = COM_get_dataitem_handle( (wname+".nc").c_str());
      int centers_handle = COM_get_dataitem_handle( (wname+".facecenters").c_str());
      COM_call_function( SURF_centers, &nc_handle, &centers_handle);

      for ( int i=0; i<npanes; ++i) {
	int pid = pane_ids[i];
	MAP::Vector_3<double> *nrms, *centers; 
	COM_get_array( (wname+".facenormals").c_str(), pid, (void**)&nrms);
	COM_get_array( (wname+".facecenters").c_str(), pid, (void**)&centers);
      
	int *cnstr_types, nelems, ngelems;
	COM_get_array( (wname+".cnstr_types_facial").c_str(),
		       pid, (void**)&cnstr_types);
	COM_get_size( (wname+".facenormals").c_str(), pid, &nelems, &ngelems);

	double eps = 1.e-6;
	// Loop through all faces
	for ( int j=0; j<nelems-ngelems; ++j) {
	  if ( centers[j][0] < eps || centers[j][0]+eps > 1.76314)
	    cnstr_types[j] = 2;
	  else if ( centers[j][0] - 0.004826 < eps)
	    cnstr_types[j] = -'x';
	  else if ( centers[j][0] >= 1.6962 && std::abs(nrms[j][0])<0.1)
	    cnstr_types[j] = 'x';
	}
      }
    }
  }

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);

  // Set speeds based on constraints
  for ( int i=0; i<npanes; ++i) {
    int pid = pane_ids[i];

#if 1
    double *bnd;
    // Set bounds
    COM_set_size( (wname+".cnstr_bound").c_str(), pid, 1);
    COM_resize_array( (wname+".cnstr_bound").c_str(), pid, (void**)&bnd);

    bnd[0] = -'x'; bnd[1] = 0; bnd[2] = 0.0635;
    bnd[3] = bnd[4] = bnd[5] = 0; // Center
    bnd[6] = 0; bnd[7] = bnd[8] = -0.07;       // Lower bound of bounding box
    bnd[9] = 1.7; bnd[10] = bnd[11] = 0.07; // Upper bound of bounding box
#endif

    // Set speeds
    int *cnstr_types, nelems, ngelems;
    double *spds;

    COM_get_array( (wname+".cnstr_types_facial").c_str(),
		   pid, (void**)&cnstr_types);
    COM_get_size( (wname+".facenormals").c_str(), pid, &nelems, &ngelems);
    COM_get_array( (wname+".spds").c_str(), pid, (void**)&spds);

    double *flags;
    COM_get_array( (wname+".flag_spd").c_str(), pid, (void**)&flags);

    // Loop through all faces
    for ( int j=0; j<nelems-ngelems; ++j) {
      if ( cnstr_types[j]) {
	flags[j] = spds[j] = 0;
      }
      else {
	spds[j] = cntr_param.speed;
	flags[j] = 1.;
      }
    }
  }

#if 1
  int PROP_set_bounds = COM_get_function_handle( "PROP.set_bounds");
  int bnd_handle = COM_get_dataitem_handle( (wname+".cnstr_bound").c_str());
  COM_call_function( PROP_set_bounds, &bnd_handle);
#endif

  COM_free_buffer( &pane_ids);  
}

// Initialize the constraints for Star slice.
void init_constraints_myacm( const string &wname, 
			     const Control_parameter &cntr_param) {
  int cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types_facial").c_str());

  int SURF_normals = COM_get_function_handle( "SURF.compute_element_normals");
  int nrms_handle = COM_get_dataitem_handle( (wname+".facenormals").c_str());
  COM_call_function( SURF_normals, &nrms_handle);

  int SURF_centers = COM_get_function_handle( "SURF.interpolate_to_centers");
  int nc_handle = COM_get_dataitem_handle( (wname+".nc").c_str());
  int centers_handle = COM_get_dataitem_handle( (wname+".facecenters").c_str());
  COM_call_function( SURF_centers, &nc_handle, &centers_handle);

  // Set speed for different panes.
  int *pane_ids, npanes;
  COM_get_panes( wname.c_str(), &npanes, &pane_ids);

  double zero = 0;
  int spd_elms = COM_get_dataitem_handle( (wname+".flag_spd").c_str());
  int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
  COM_call_function( BLAS_copy_scalar, &zero, &spd_elms);
      
  for ( int i=0; i<npanes; ++i) {
    int pid = pane_ids[i];
    MAP::Vector_3<double> *nrms, *centers;
    COM_get_array( (wname+".facenormals").c_str(), pid, (void**)&nrms);
    COM_get_array( (wname+".facecenters").c_str(), pid, (void**)&centers);

    int *cnstr_types;
    COM_get_array( (wname+".cnstr_types_facial").c_str(), 
		   pid, (void**)&cnstr_types);
    double *spds;
    COM_get_array( (wname+".spds").c_str(), pid, (void**)&spds);
    double *flags;
    COM_get_array( (wname+".flag_spd").c_str(), pid, (void**)&flags);
    
    int nelems, ngelems;
    COM_get_size( (wname+".facenormals").c_str(), pid, &nelems, &ngelems);

    double eps = 5.e-2;
    if ( cntr_param.start == 0) {
      // Loop through all faces
      for ( int j=0; j<nelems-ngelems; ++j) {
	if ( std::abs(nrms[j][2])<eps && 
	     square(centers[j][0])+square(centers[j][1]) >= square(0.0065)) {
	  cnstr_types[j] = 'z';
	  flags[j] = spds[j]= 0;
	}
	else {
	  cnstr_types[j] = 0;
	  flags[j] = 1.; spds[j]= cntr_param.speed;
	}
      }
    }
    else {
      for ( int j=0; j<nelems-ngelems; ++j) {
	if ( cnstr_types[j] == 0) {
	  flags[j] = 1.; spds[j]= cntr_param.speed;
	}
	else {
	  flags[j] = spds[j]= 0;
	}
      }
    }
  }
  COM_free_buffer( &pane_ids);

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}


// Initialize the constraints for Lab scale rocket
void init_constraints_labscale( const string &wname, 
				const Control_parameter &cntr_param) {
  int BLAS_copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
  int cnstr_handle = COM_get_dataitem_handle( (wname+".cnstr_types_facial").c_str());
  int zero = 0;
  COM_call_function( BLAS_copy_scalar, &zero, &cnstr_handle);

  int SURF_normals = COM_get_function_handle( "SURF.compute_element_normals");
  int nrms_handle = COM_get_dataitem_handle( (wname+".facenormals").c_str());
  COM_call_function( SURF_normals, &nrms_handle);

  int SURF_centers = COM_get_function_handle( "SURF.interpolate_to_centers");
  int nc_handle = COM_get_dataitem_handle( (wname+".nc").c_str());
  int centers_handle = COM_get_dataitem_handle( (wname+".facecenters").c_str());
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
    
    int nelems, ngelems;
    COM_get_size( (wname+".facenormals").c_str(), pid, &nelems, &ngelems);

    double eps = 1.e-2;
    // Loop through all faces
    for ( int j=0; j<nelems-ngelems; ++j) {
      if ( std::abs(std::abs(nrms[j][0])-1)<eps && centers[j][0] < 0.153)
	cnstr_types[j] = 't';  // Move within yz plane.
    }
  }
  COM_free_buffer( &pane_ids);

  int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
  COM_call_function( PROP_set_cnstr, &cnstr_handle);
}

// Initialize the constraints
void init_constraints( const string &wname, 
		       const Control_parameter &cntr_param) {

  if ( wname.substr(0,5) == "acmfl") 
    init_constraints_acmflu( wname);
  else if ( wname.substr(0,9) == "starslice") 
    init_constraints_starslice( wname, cntr_param);
  else if ( wname.substr(0,7) == "staraft") 
    init_constraints_staraft( wname, cntr_param);
  else if ( wname.substr(0,8) == "labscale") 
    init_constraints_labscale( wname, cntr_param);
  else if ( wname.substr(0,3) == "acm") 
    init_constraints_myacm( wname, cntr_param);
  else {
    int cnstr_handle = COM_get_dataitem_handle((wname+".cnstr_type").c_str());
    if ( cnstr_handle>=0) {
      int PROP_set_cnstr = COM_get_function_handle( "PROP.set_constraints");
      COM_call_function( PROP_set_cnstr, &cnstr_handle);
    }
  }
}

void init_attributes( const string &wname, 
		      const Control_parameter &cntr_param) {
  COM_new_dataitem((wname+".flag").c_str(), 'p', COM_DOUBLE, 1, "");
  COM_set_size( (wname+".flag").c_str(), 0, 1);

  COM_new_dataitem((wname+".flag_spd").c_str(), 'e', COM_DOUBLE, 1, "");
  
  COM_new_dataitem((wname+".disps_total").c_str(), 'n', COM_DOUBLE, 3, "");
  // COM_new_dataitem((wname+".lambdas").c_str(), 'n', COM_DOUBLE, 3, "");
  // COM_new_dataitem((wname+".eigvecs").c_str(), 'n', COM_DOUBLE, 9, "");
  COM_new_dataitem((wname+".cnstr_types_nodal").c_str(), 'n', COM_INT, 1, "");
  COM_new_dataitem((wname+".cnstr_types_facial").c_str(), 'e', COM_INT, 1, "");
  COM_new_dataitem((wname+".cnstr_types_panel").c_str(), 'p', COM_INT, 1, "");
  COM_set_size((wname+".cnstr_types_panel").c_str(), 0, 1);

#if 1
  COM_new_dataitem((wname+".cnstr_bound").c_str(), 'p', COM_DOUBLE, 15, ""); 
  COM_set_size( (wname+".cnstr_bound").c_str(), 0, 1);
#endif

  COM_new_dataitem((wname+".spds").c_str(), cntr_param.sploc, 
		    COM_DOUBLE, 1, "m/s");
  COM_new_dataitem((wname+".disps").c_str(), 'n', COM_DOUBLE, 3, "m");

  COM_new_dataitem((wname+".faceareas").c_str(), 'e', COM_DOUBLE, 1, ""); 
  COM_new_dataitem((wname+".facenormals").c_str(), 'e', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".facecenters").c_str(), 'e', COM_DOUBLE, 3, "");

  // DataItem for storing the number of eigenvalues for each node.
  COM_new_dataitem((wname+".tangranks").c_str(), 'n', COM_INT, 1, "");

  // Ridges of each pane.
  COM_new_dataitem((wname+".ridges").c_str(), 'p', COM_INT, 2, "");
  COM_set_size((wname+".ridges").c_str(), 0, 0);

  COM_resize_array( (wname+".data").c_str());
  COM_window_init_done( wname.c_str());

  // Set speed for different panes.
  int *pane_ids, npanes;
  COM_get_panes( wname.c_str(), &npanes, &pane_ids);
  
  for ( int i=0; i<npanes; ++i) {
    double *flag;
    COM_get_array( (wname+".flag").c_str(), pane_ids[i], &flag);
    // Set constraints for nonburning panes and speeds for burning panes.
    if ( wname.substr(0,5) == "acmfl" ) {
      if ( COM_get_dataitem_handle((wname+".bcflag").c_str())>0) {
	int *bcflag;
	COM_get_array( (wname+".bcflag").c_str(), pane_ids[i], &bcflag);
	
	*flag = ( *bcflag == 1) ? cntr_param.speed : 0.;
      }
      else 
	*flag = (pane_ids[i]<=2)?cntr_param.speed:0.;
    }
    else if ( wname.substr(0,4) == "star") 
      *flag = (i%2)?0.:cntr_param.speed; 
    else if ( COM_get_dataitem_handle((wname+".bcflag").c_str())>0) {
      int *bcflag;
      COM_get_array( (wname+".bcflag").c_str(), pane_ids[i], &bcflag);
	
      *flag = ( *bcflag == 1) ? cntr_param.speed : 0.;
    }
    else
      *flag = cntr_param.speed; 

    std::cout << "Setting speed for pane " << pane_ids[i] 
	      << " to " << *flag << std::endl;
  }

  COM_free_buffer( &pane_ids);

  int spds = COM_get_dataitem_handle( (wname+".spds").c_str());
  int spd_pane = COM_get_dataitem_handle( (wname+".flag").c_str());

  int BLAS_copy = COM_get_function_handle( "BLAS.copy");
  COM_call_function( BLAS_copy, &spd_pane, &spds);

  int BLAS_div_scalar = COM_get_function_handle( "BLAS.div_scalar");
  COM_call_function( BLAS_div_scalar, &spd_pane, &cntr_param.speed, &spd_pane);

  int spd_elms = COM_get_dataitem_handle( (wname+".flag_spd").c_str());
  COM_call_function( BLAS_copy, &spd_pane, &spd_elms);
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

double compute_area( const string &wname) {
  static int SURF_area = 0, BLAS_mul, BLAS_sum_scalar, area_hdl, flag_hdl;

  if ( SURF_area==0) {
    SURF_area = COM_get_function_handle( "SURF.compute_element_areas");
    BLAS_mul = COM_get_function_handle( "BLAS.mul");
    BLAS_sum_scalar = COM_get_function_handle( "BLAS.sum_scalar_MPI");

    area_hdl = COM_get_dataitem_handle( (wname+".faceareas").c_str());
    flag_hdl = COM_get_dataitem_handle( (wname+".flag_spd").c_str());
  }

  COM_call_function( SURF_area, &area_hdl);
  COM_call_function( BLAS_mul, &area_hdl, &flag_hdl, &area_hdl);

  double sum;
  MPI_Comm comm=MPI_COMM_WORLD;
  COM_call_function( BLAS_sum_scalar, &area_hdl, &sum, &comm);
  return sum;
}

double compute_volume( const string &wname) {
  static int SURF_vol = 0, mesh_hdl;

  if ( SURF_vol==0) {
    SURF_vol = COM_get_function_handle( "SURF.compute_signed_volumes");

    mesh_hdl = COM_get_dataitem_handle( (wname+".mesh").c_str());
  }

  double vol;
  COM_call_function( SURF_vol, &mesh_hdl, &vol);

  return vol;
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
  
  // Initialize the dataitems and output the initial solution.
  init_attributes( wname, cntr_param);

  // DataItems
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

  init_constraints( wname, cntr_param);

  if ( cntr_param.perturb > 0)
    COM_call_function( PROP_perturb, &pmesh, &cntr_param.perturb);

  if ( cntr_param.start==0) {
    double zero=0.;
    // Run one step to detect features.
    COM_call_function( PROP_propagate, &pmesh, &spds, &zero, &disps, &zero);
  }
  output_solution( wname, "00000", cntr_param.format);

  COM_set_profiling(1);
  COM_set_profiling_barrier( PROP_propagate, MPI_COMM_WORLD);

  std::ofstream areas, vols;
  if ( rank==0) { 
    areas.open((wname+"_areas.m").c_str());
    vols.open((wname+"_vols.m").c_str());

    areas.precision(10); vols.precision(10);  cout.precision(10);
  }

  double a=compute_area( wname), v=compute_volume(wname);

  if ( rank==0) {
    std::cout << "Burning area is " << a 
	      << " and total volume is " << v << std::endl;
    areas << wname << "_as = [0 " << a << ";..." << std::endl;
    vols << wname << "_vs = [0 " << v << ";..." << std::endl;
  }

  char fname[40];
  std::sprintf(fname, "timedata_%03d.txt", rank);

#ifdef MESH_ADAPT
  AdaptCOMWindow acw( wname.c_str());
  acw.set_option( "adapt_iter", 1); // cntr_param.adapt_iter);
  acw.set_option( "redist_iter", 0); // cntr_param.do_redist);

  acw.set_option( "do_collapse", cntr_param.do_collapse);
  acw.set_option( "do_flip", cntr_param.do_flip);
  acw.set_option( "do_split", cntr_param.do_split);
  acw.set_option( "refine", cntr_param.refine);
  if ( cntr_param.edge_len>0)
    acw.set_option( "edge_len", cntr_param.edge_len);

  if ( !cntr_param.fangle.empty())
    acw.set_option( "fangle_weak", 
		    std::atof(cntr_param.fangle.c_str())/180*M_PI);
  if ( !cntr_param.sangle.empty())
    acw.set_option( "fangle_strong", 
		    std::atof(cntr_param.sangle.c_str())/180*M_PI);
  if ( !cntr_param.eigthres.empty())
    acw.set_option( "eig_thres", std::atof(cntr_param.eigthres.c_str()));
  if ( cntr_param.collapse_ratio>0)
    acw.set_option( "collapse_ratio", cntr_param.collapse_ratio);
  if ( cntr_param.split_angle > 0)
    acw.set_option( "split_angle", cntr_param.split_angle);

  acw.set_option( "anisotropic", cntr_param.smoother == "anisotropic");
  if ( !cntr_param.feature_layer.empty())
    acw.set_option( "feature_layer", std::atof(cntr_param.feature_layer.c_str()));
#endif

  double t=cntr_param.start*cntr_param.timestep; 
  for ( int i=1+cntr_param.start; i<=cntr_param.steps; ++i, t+=cntr_param.timestep) {
    cntr_param.start = i; // Save current time step
    if ( rank==0) cout << "Step " << i << endl;

    double local_t = 0, t_rem = cntr_param.timestep, t_elapsed = 0;
    while ( t_rem > 0) {
      // If speed is 0, then do smoothing only.
      if ( cntr_param.speed==0) t_rem=0;

#ifdef MESH_ADAPT
      if ( i>1 && cntr_param.adapt_iter>0 && 
	   ( local_t!=0 && t_elapsed < cntr_param.subadapt*cntr_param.timestep ||
	     local_t==0 && (i-1)%cntr_param.remesh_interval==0)) {
      
        for (int j=0; j<cntr_param.adapt_iter; ++j) {
	  if ( !acw.adapt()) break;

	  // Reset the constraints
	  init_constraints( wname, cntr_param);

	  // If the window has changed, reinitialize Rocprop.
	  COM_call_function( PROP_init, &pmesh);

	  for ( int k=0; k<std::max(1,cntr_param.do_redist); ++k) {
	    // Perform mesh smoothing
	    double zero=0.;
	    std::cout << "Perform anisotropic smoothing" << std::endl;
	    COM_call_function( PROP_propagate, &pmesh, &spds, 
			       &zero, &disps, &t_elapsed);
	    COM_call_function( BLAS_add, &nc, &disps, &nc);
	    COM_call_function( BLAS_add, &disps_total, &disps, &disps_total);
	  }
	}
      }
#endif

      COM_call_function( PROP_propagate, &pmesh, &spds, 
			 &t_rem, &disps, &t_elapsed);
            
      COM_call_function( BLAS_add, &nc, &disps, &nc);
      COM_call_function( BLAS_add, &disps_total, &disps, &disps_total);

      local_t = cntr_param.timestep - (t_rem-t_elapsed);
      t_rem = cntr_param.timestep - local_t;

      if ( t_rem>0 && t_elapsed*1.e3<cntr_param.timestep) {
	std::cout << "Time step too small. Stopping..." << std::endl;
	exit(-1);
      }
    }

    a=compute_area( wname); v=compute_volume(wname);
    if ( rank==0) {
      std::cout << "Burning area is " << a 
		<< " and total volume is " << v << std::endl;
      double t = i*cntr_param.timestep;
      areas << '\t' << t << '\t' << a << ";..." << std::endl;
      vols << '\t' << t << '\t' << v << ";..." << std::endl;
    }

    if ( i%cntr_param.interval == 0) {
      char   steps[10];
      if ( cntr_param.steps>=1.e6)
	std::sprintf( steps, "%07d", i);
      else if ( cntr_param.steps>=1.e5)
	std::sprintf( steps, "%06d", i);
      else
	std::sprintf( steps, "%05d", i);

      output_solution( wname, steps, cntr_param.format);

      COM_print_profile( fname, "Proptest");
    }
  }

  areas << "];" << std::endl;
  vols << "];" << std::endl;

  COM_finalize();
}






