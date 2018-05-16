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
// $Id: advectest.C,v 1.6 2008/12/06 08:45:28 mtcampbe Exp $

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
#include "PointPropagate.h"

#ifdef MESH_ADAPT
#include "AdaptCOMWindow.h"
#endif


using namespace std;
using namespace PROP;

COM_EXTERN_MODULE( Simpal);
COM_EXTERN_MODULE( SurfMap);
COM_EXTERN_MODULE( Rocprop);
COM_EXTERN_MODULE( SurfUtil);
COM_EXTERN_MODULE( SimOut);

void load_modules() {
  COM_LOAD_MODULE_STATIC_DYNAMIC(Simpal,   "BLAS");
  COM_LOAD_MODULE_STATIC_DYNAMIC(SurfMap,  "MAP");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocprop,  "PROP");
  COM_LOAD_MODULE_STATIC_DYNAMIC(SurfUtil, "SURF");
  COM_LOAD_MODULE_STATIC_DYNAMIC(SimOut,   "OUT");
}

static int rank = 0;

void print_usage( int argc, char *argv[]) {
  if ( argc <= 2) {
    cout << "Usage: " << argv[0] << " <surffile> <controlfile>" << std::endl;
    
    exit(-1);
  }
}

struct Control_parameter {
  Control_parameter() : speed(0), timestep(0), steps(0), start(0), interval(0), 
			remesh_interval(10), test("rotate"), 
			adapt_iter(2), do_redist(0), do_collapse(1), do_flip(1), 
			do_split(1), collapse_ratio(0), split_angle(0), refine(1) {}
  
  string method;
  string wavefrontal;
  string normaldif;
  string eigthres;
  string courant;
  string fangle;
  string sangle;
  string smoother;
  string rediter;
  string weight;
  string conserv;

  double speed;
  double timestep;
  int    steps;
  int    start;
  int    interval;
  int    remesh_interval;

  string test;
  string verbose;

  // The following options are for mesh adaptation
  int    adapt_iter;
  int    do_redist;
  int    do_collapse;
  int    do_flip;
  int    do_split;
  double collapse_ratio;
  double split_angle;
  double refine;
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
   * conserv: 0          # Whether or not to conserve volume
   * speed: 0.1          # Speed
   * timestep: 0.001     # time step
   * steps: 100          # number of time steps
   * interval: 10        # output intervals
   * remesh_interval: 10 # remesh intervals
   * weight: 3           # angle weighted
   * verbose: 1          # verbose level
   * test: rotate        # test case
   *
   * adapt_iter: 2
   * do_flip: 1
   * do_split: 1
   * do_redist: 0
   * do_collapse: 1
   * collapse_ratio: 0.2
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
    else if ( keyword == "conserv")
      is >> cp.conserv;
    else if ( keyword == "speed")
      is >> cp.speed;
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
    else if ( keyword == "test")
      is >> cp.test;
    else if ( keyword == "verbose")
      is >> cp.verbose;
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
    else if ( keyword == "refine")
      is >> cp.refine;
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

  if ( !cntr_param.conserv.empty()) {
    COM_call_function( PROP_set_option, "conserv", cntr_param.conserv.c_str());
    if ( rank==0) std::cout << "Set conserv to " << cntr_param.conserv << std::endl;
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

template <class T>
T square(T t) { return t*t; }


void init_attributes( const string &wname, 
		      const Control_parameter &cntr_param) {

  COM_new_dataitem((wname+".disps_total").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".cnstr_types_nodal").c_str(), 'n', COM_INT, 1, "");
  COM_new_dataitem((wname+".cnstr_types_facial").c_str(), 'e', COM_INT, 1, "");
  COM_new_dataitem((wname+".cnstr_types_panel").c_str(), 'p', COM_INT, 1, "");
  COM_set_size((wname+".cnstr_types_panel").c_str(), 0, 1);
  
  COM_new_dataitem((wname+".qvels").c_str(), 'e', COM_DOUBLE, 9, "m/s");
  COM_new_dataitem((wname+".vvels").c_str(), 'n', COM_DOUBLE, 3, "m/s");
  COM_new_dataitem((wname+".disps").c_str(), 'n', COM_DOUBLE, 3, "m");
  COM_new_dataitem((wname+".faceareas").c_str(), 'e', COM_DOUBLE, 1, ""); 

  // Ridges of each pane.
  COM_new_dataitem((wname+".ridges").c_str(), 'p', COM_INT, 2, "");
  COM_set_size((wname+".ridges").c_str(), 0, 0);

  // COM_new_dataitem((wname+".disps_novis").c_str(), 'n', COM_DOUBLE, 3, "");
  // COM_new_dataitem((wname+".facenormals").c_str(), 'e', COM_DOUBLE, 3, "");
  // COM_new_dataitem((wname+".facecenters").c_str(), 'e', COM_DOUBLE, 3, "");
  // COM_new_dataitem((wname+".faceheights").c_str(), 'e', COM_DOUBLE, 1, "");

  // Dataitem for storing the number of eigenvalues for each node.
  // COM_new_dataitem((wname+".lambdas").c_str(), 'n', COM_DOUBLE, 3, "");
  // COM_new_dataitem((wname+".eigvecs").c_str(), 'n', COM_DOUBLE, 9, "");
  COM_new_dataitem((wname+".tangranks").c_str(), 'n', COM_CHAR, 1, "");
  COM_new_dataitem((wname+".scales").c_str(), 'n', COM_DOUBLE, 1, "");

  COM_resize_array( (wname+".data").c_str());
  COM_window_init_done( wname.c_str());
}

// return the volume of the domain
double output_solution( const string &wname, const char *timelevel, 
			double ref=0.) {
  static int OUT_write = 0, hdl;

  if ( OUT_write==0) {
    OUT_write = COM_get_function_handle( "OUT.write_dataitem");
    hdl = COM_get_dataitem_handle( (wname+".all").c_str());

    int OUT_set = COM_get_function_handle( "OUT.set_option");
    COM_call_function( OUT_set, "format", "HDF");
  }

  if ( timelevel) {
    std::string fname = wname+"_"+timelevel;

    if ( !COMMPI_Initialized()) fname.append(".hdf");
    else fname.append("_");
  
    COM_call_function( OUT_write, (char*)fname.c_str(),
		       &hdl, (char*)wname.c_str(), timelevel);
  }

  static COM::DataItem *mesh=NULL;
  if ( mesh==NULL) 
    mesh = COM_get_com()->get_window_object( wname)->dataitem( "mesh");

  double vol;
  SURF::Rocsurf::compute_signed_volumes( mesh, &vol);
  
  std::cout << "Volume of object is " << vol << std::endl;
  if ( ref!=0) {
    std::cout << "Relative error in volume is " 
	      << std::abs((vol-ref)/ref) << std::endl;
  }

  return vol;
}

double compute_area( const string &wname) {
  static int SURF_area = 0, BLAS_mul, BLAS_sum_scalar, area_hdl;

  if ( SURF_area==0) {
    SURF_area = COM_get_function_handle( "SURF.compute_element_areas");
    BLAS_mul = COM_get_function_handle( "BLAS.mul");
    BLAS_sum_scalar = COM_get_function_handle( "BLAS.sum_scalar_MPI");

    area_hdl = COM_get_dataitem_handle( (wname+".faceareas").c_str());
  }

  COM_call_function( SURF_area, &area_hdl);

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
  
  // Initialize the attributes and output the initial solution.
  init_attributes( wname, cntr_param);

  // Attributes
  int pmesh = COM_get_dataitem_handle( (wname+".pmesh").c_str());
  int nc = COM_get_dataitem_handle( (wname+".nc").c_str());
  int vvels = COM_get_dataitem_handle( (wname+".vvels").c_str());
  int qvels = COM_get_dataitem_handle( (wname+".qvels").c_str());
  int disps = COM_get_dataitem_handle( (wname+".disps").c_str());
  int disps_total = COM_get_dataitem_handle( (wname+".disps_total").c_str());

  // Funcitions
  int PROP_init = COM_get_function_handle( "PROP.initialize");
  int PROP_propagate = COM_get_function_handle( "PROP.propagate");
  int BLAS_add = COM_get_function_handle( "BLAS.add");

  if ( rank==0) cout << "Propagating interface..." << endl;
  // Initialize control parameters
  COM_call_function( PROP_init, &pmesh);
  init_parameters( cntr_param);

  Translate_speed trans_spd( cntr_param.speed);
  Rotate_speed    rotate_spd;
  Vortex_flow     vortex_spd(cntr_param.speed);
  LeVeque_flow    leveque_spd(cntr_param.speed);

  bool  use_quadpoints = false;
  Speed *spd=NULL;
  std::cerr << "Initializing for test case " << cntr_param.test << std::endl;
  if ( cntr_param.test == "trans") {
    spd = &trans_spd;
  }
  else if ( cntr_param.test == "rotate") {
    // use_quadpoints = true;
    spd = &rotate_spd;
    rescale_object( wname, 1, SURF::Vector_3<double>(5./3., 0, 0));
  }
  else if ( cntr_param.test == "vortexcube") {
    spd = &vortex_spd;
    if ( cntr_param.speed>0 && cntr_param.start==0)
      rescale_object( wname, 0.3, SURF::Vector_3<double>(0.5, 0.75, 0.5));
  }
  else if ( cntr_param.test == "vortex") {
    spd = &vortex_spd;
    if ( cntr_param.speed>0 && cntr_param.start==0)
      rescale_object( wname, 0.15, SURF::Vector_3<double>(0.5, 0.75, 0.5));
  }
  else if ( cntr_param.test == "leveque") {
    spd = &leveque_spd;
    if ( cntr_param.speed>0 && cntr_param.start==0)
      rescale_object( wname, 0.15, SURF::Vector_3<double>(0.35, 0.35, 0.35));
  }
  else {
    COM_assertion_msg(false, "Unkonwn test case");
  }

  if ( cntr_param.start==0) {
    double zero=0.;
    // Perform one step of smoothing to detect features.
    COM_call_function( PROP_propagate, &pmesh, &vvels, &zero, &disps, &zero);
  }

  double vol = output_solution( wname, cntr_param.start==0?"00000":NULL);

  COM_set_profiling(1);
  COM_set_profiling_barrier( PROP_propagate, MPI_COMM_WORLD);

  std::ofstream areas, vols;
  areas.precision(10); vols.precision(10);
  cout.precision(10);

  if ( rank==0) areas.open((wname+"_areas.m").c_str());
  if ( rank==0) vols.open((wname+"_vols.m").c_str());

  double a=compute_area( wname), v=compute_volume(wname);

  if ( rank==0) {
    std::cout << "Area is " << a << " and total volume is " << v << std::endl;
    areas << wname << "_as = [0 " << a << ";..." << std::endl;
    vols << wname << "_vs = [0 " << v << ";..." << std::endl;
  }

  char fname[40];
  std::sprintf(fname, "timedata_%03d.txt", rank);

#ifdef MESH_ADAPT
  AdaptCOMWindow acw( wname.c_str());

  acw.set_option( "adapt_iter", cntr_param.adapt_iter);
  acw.set_option( "redist_iter", cntr_param.do_redist);

  acw.set_option( "do_collapse", cntr_param.do_collapse);
  acw.set_option( "do_flip", cntr_param.do_flip);
  acw.set_option( "do_split", cntr_param.do_split);
  acw.set_option( "refine", cntr_param.refine);

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
#endif

  double t=cntr_param.start*cntr_param.timestep; 
  for ( int i=1+cntr_param.start; i<=cntr_param.steps; ++i, t+=cntr_param.timestep) {
    if ( rank==0) cout << "Step " << i << endl;

    double local_t = 0, t_rem = cntr_param.timestep, t_elapsed;
    while ( t_rem > 0) {
      if ( use_quadpoints) {
	// Initialize the velocity at quadrature points.
	PointPropagate::propagate_faces( wname, *spd, t,
					 cntr_param.timestep, "qvels");
      }
      else {
	// Initialize the velocity at vertices.
	PointPropagate::propagate_nodes( wname, *spd, t,
					 cntr_param.timestep, "vvels");
      }
      
      // If speed is 0, then do smoothing only.
      if ( cntr_param.speed==0) t_rem=0;

      COM_call_function( PROP_propagate, &pmesh, use_quadpoints?&qvels:&vvels, 
			 &t_rem, &disps, &t_elapsed);
      
      
      COM_call_function( BLAS_add, &nc, &disps, &nc);
      COM_call_function( BLAS_add, &disps_total, &disps, &disps_total);

      local_t = cntr_param.timestep - (t_rem-t_elapsed);
      if ( local_t*1.e3<cntr_param.timestep) {
	std::cout << "Time step too small. Stopping..." << std::endl;
	exit(-1);
      }

#ifdef MESH_ADAPT
      if ( (t_elapsed<0.1*t_rem || t_elapsed==t_rem &&
	    i%cntr_param.remesh_interval==0) && cntr_param.adapt_iter>0) {

	for (int j=-1; j<cntr_param.adapt_iter; ++j) {
	  if ( j>=0) { // Note: we perform a step of smoothing first.
	    if ( acw.adapt()==0) break;
	    // If the window has changed, reinitialize Rocprop.
	    COM_call_function( PROP_init, &pmesh);
	  }

	  for ( int k=0; k<cntr_param.do_redist; ++k) {
	    // Perform mesh smoothing
	    double zero=0.;
	    std::cout << "Perform anisotropic smoothing" << std::endl;
	    COM_call_function( PROP_propagate, &pmesh, &vvels, 
			       &zero, &disps, &t_elapsed);
	    COM_call_function( BLAS_add, &nc, &disps, &nc);
	    COM_call_function( BLAS_add, &disps_total, &disps, &disps_total);
	  }
	}
      }
#endif

      t_rem = cntr_param.timestep - local_t;

    }

    a=compute_area( wname); v=compute_volume(wname);
    if ( rank==0) {
      std::cout << "Area is " << a << " and total volume is " << v << std::endl;
      double t = i*cntr_param.timestep;
      areas << '\t' << t << '\t';
      areas.precision(16); areas << a << ";..." << std::endl;
      vols << '\t' << t << '\t';
      vols.precision(16);  vols << v << ";..." << std::endl;
    }

    if ( i%cntr_param.interval == 0) {
      char   steps[10];
      if ( cntr_param.steps>=1.e6)
	std::sprintf( steps, "%07d", i);
      else if ( cntr_param.steps>=1.e5)
	std::sprintf( steps, "%06d", i);
      else
	std::sprintf( steps, "%05d", i);

      output_solution( wname, steps, vol);

      COM_print_profile( fname, "Proptest");
    }
  }

  areas << "];" << std::endl;
  vols << "];" << std::endl;

  COM_finalize();
}






