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
// $Id: compare_normals.C,v 1.8 2008/12/06 08:45:28 mtcampbe Exp $

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

typedef MAP::Vector_3<double> Vector_3;

COM_EXTERN_MODULE( Simpal);
COM_EXTERN_MODULE( SurfMap);
COM_EXTERN_MODULE( SurfUtil);
COM_EXTERN_MODULE( Rocprop);
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
			steps(0), interval(0), func(1) {}
  string method;
  string wavefrontal;
  string normaldif;
  string eigthres;
  string courant;
  string fangle;
  double perturb;
  double speed;
  char   sploc;
  double timestep;
  int    steps;
  int    interval;

  int    func;
  string verbose;
};

void read_control_file( const char *fname, Control_parameter &cp) {
  /* Sample control file:
   * method: fo          # method: fo and mp
   * wavefrontal: 1      # wavefrontal condition
   * normaldif: 1        # normal diffusion
   * eigthres: 1.e-4     # threshold for null space: 0..1 (1.e-4 default)
   * courant: 0.5        # courant constant
   * fangle: 35          # feature edge angle: between 0 and 180
   * speed: 0.1          # Speed
   * sploc: e            # location of speed: n or e
   * timestep: 0.001     # time step
   * steps: 100          # number of time steps
   * interval: 10        # output intervals
   * function: 1         # function
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
    else if ( keyword == "function")
      is >> cp.func;
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
  int npanes = im_reader.read_winmesh( fname, wname); 
  COM_assertion_msg( npanes>=0, "Failed to read in mesh file. File empty?");

  return wname;
}

void init_attributes( const string &wname) {
  COM_new_dataitem((wname+".facenormals").c_str(), 'e', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".facecenters").c_str(), 'e', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".fnormals_ana").c_str(), 'e', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".fnormals_cr").c_str(), 'e', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".ferrs").c_str(), 'e', COM_DOUBLE, 1, "");
  COM_new_dataitem((wname+".verrs").c_str(), 'n', COM_DOUBLE, 1, "");

  COM_new_dataitem((wname+".vnormals_ana").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".vnormals_nw").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".vnormals_area").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".vnormals_angle").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".vnormals_sphere").c_str(), 'n', COM_DOUBLE, 3, "");

  COM_new_dataitem((wname+".vnormals_eig").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".vnormals_mq").c_str(), 'n', COM_DOUBLE, 3, "");

  // Dataitem for storing the number of eigenvalues for each node.
  COM_new_dataitem((wname+".lambdas").c_str(), 'n', COM_DOUBLE, 3, "");
  COM_new_dataitem((wname+".eigvecs").c_str(), 'n', COM_DOUBLE, 9, "");
  COM_new_dataitem((wname+".tangranks").c_str(), 'n', COM_CHAR, 1, "");
  COM_new_dataitem((wname+".scales").c_str(), 'n', COM_DOUBLE, 1, "");

  COM_new_dataitem((wname+".spds").c_str(), 'e', COM_DOUBLE, 1, "m/s");
  COM_new_dataitem((wname+".disps").c_str(), 'n', COM_DOUBLE, 3, "m");

  COM_resize_array( (wname+".data").c_str());
  COM_window_init_done( wname.c_str());
}

void func1( Vector_3 &pnt, Vector_3 &nrm) {
  // z=sqrt(2-x^2-y^2)
  pnt[2] = std::sqrt(2-pnt[0]*pnt[0] - pnt[1]*pnt[1]);
  if ( pnt[2]==0) 
    nrm = Vector_3( -pnt[0], -pnt[1], 0);
  else
    nrm = Vector_3::cross_product( Vector_3(0, 1, -pnt[1]/pnt[2]),
				   Vector_3(1, 0, -pnt[0]/pnt[2]));
  nrm.normalize();
}

const double pi = 3.14159265358979;

void func2( Vector_3 &pnt, Vector_3 &nrm) {
  // z = x^2 + 0.25*sin( 2pi*y)
  pnt[2] = pnt[0]*pnt[0] + 0.25*std::sin(2*pi*pnt[1]);
  nrm = Vector_3::cross_product( Vector_3(0, 1, 0.5*pi*std::cos(2*pi*pnt[1])),
				 Vector_3(1, 0, 2*pnt[0]));
  nrm.normalize();
}

void func3( Vector_3 &pnt, Vector_3 &nrm) {
  // z = 0.25*sin( pi*x) + 0.25*cos( pi*y)
  pnt[2] = 0.25*(std::sin(pi*pnt[0]) + std::cos(pi*pnt[1]));
  nrm = Vector_3::cross_product( Vector_3(0, 1, -0.25*pi*std::sin(pi*pnt[1])),
				 Vector_3(1, 0, 0.25*pi*std::cos(pi*pnt[0])));
  nrm.normalize();
}

void init_normals( const string &wname, int func) {
  int *pane_ids, npanes;
  COM_get_panes( wname.c_str(), &npanes, &pane_ids);

  // Compute vertex normals
  for ( int i=0; i<npanes; ++i) {
    int stride, nn;
    Vector_3 *coors, *nrms;
    COM_get_array( (wname+".nc").c_str(), pane_ids[i], 
		   &(void*&)coors, &stride, &nn);
    COM_get_array( (wname+".vnormals_ana").c_str(), pane_ids[i], &(void*&)nrms);

    for ( int j=0; j<nn; ++j) {
      switch (func) {
      case 1: { func1( coors[j], nrms[j]); break; }
      case 2: { func2( coors[j], nrms[j]); break; }
      case 3: { func3( coors[j], nrms[j]); break; }
      }
    }
  }

  // Interpolate coordinates to face centers
  int SURF_n2f = COM_get_function_handle( "SURF.interpolate_to_centers");
  int nc = COM_get_dataitem_handle( (wname+".nc").c_str());
  int cnts = COM_get_dataitem_handle( (wname+".facecenters").c_str());
  COM_call_function( SURF_n2f, &nc, &cnts);

  // Compute face normals
  for ( int i=0; i<npanes; ++i) {
    int stride, ne;
    Vector_3 *coors, *nrms;
    COM_get_array( (wname+".facecenters").c_str(), pane_ids[i], 
		   &(void*&)coors, &stride, &ne);
    COM_get_array( (wname+".fnormals_ana").c_str(), pane_ids[i], &(void*&)nrms);

    for ( int j=0; j<ne; ++j) {
      switch (func) {
      case 1: { func1( coors[j], nrms[j]); break; }
      case 2: { func2( coors[j], nrms[j]); break; }
      case 3: { func3( coors[j], nrms[j]); break; }
      }
    }
  }

  COM_free_buffer( &pane_ids);
}

// RMS and L-inf error in face normals
std::pair<double,double> 
compute_fnormal_error( const int ref_nrms_hdl,
		       const int cur_nrms_hdl, 
		       const int buf_hdl, const char *scheme=NULL) {
  static int SURF_normals=0, SURF_integrate=0, SURF_comparea=0;
  static int BLAS_dot=0, BLAS_mul_scalar=0, BLAS_sub_scalar=0;
  static int BLAS_sum_scalar_MPI=0, BLAS_max_scalar_MPI=0;

  if ( SURF_normals==0) { // Load the functions
    SURF_normals = COM_get_function_handle( "SURF.compute_element_normals");
    SURF_integrate = COM_get_function_handle( "SURF.integrate");
    SURF_comparea = COM_get_function_handle( "SURF.compute_element_areas");
  
    BLAS_dot = COM_get_function_handle( "BLAS.dot");
    BLAS_mul_scalar = COM_get_function_handle( "BLAS.mul_scalar");
    BLAS_sub_scalar = COM_get_function_handle( "BLAS.sub_scalar");
    BLAS_sum_scalar_MPI = COM_get_function_handle( "BLAS.sum_scalar_MPI");
    BLAS_max_scalar_MPI = COM_get_function_handle( "BLAS.max_scalar_MPI");
  }
  
  // If two arguments are the same, then initialize ref_nrms
  if ( ref_nrms_hdl == cur_nrms_hdl) {
    COM_call_function( SURF_normals, &ref_nrms_hdl);
    return std::pair<double,double>(0.,0.);
  }
  else {
    MPI_Comm comm = MPI_COMM_WORLD;
  
    // Compute root-mean-square as \sqrt(\int_{2-ref_nrms*cur_nrm))
    double err_rms, err_max, area, minus_two=-2; 
    // Compute area
    COM_call_function( SURF_comparea, &buf_hdl);
    COM_call_function( BLAS_sum_scalar_MPI, &buf_hdl, &area, &comm);
    
    COM_call_function( BLAS_dot, &cur_nrms_hdl, &ref_nrms_hdl, &buf_hdl);
    COM_call_function( BLAS_mul_scalar, &buf_hdl, &minus_two, &buf_hdl);
    COM_call_function( BLAS_sub_scalar, &buf_hdl, &minus_two, &buf_hdl);

    COM_call_function( BLAS_max_scalar_MPI, &buf_hdl, &err_max, &comm);
    COM_call_function( SURF_integrate, &buf_hdl, &err_rms);
    
    std::pair<double,double> errs(err_rms/area, err_max);
    if ( scheme) {
      std::cout << "Root-mean-square error in " << scheme
		<< " is " << std::sqrt(errs.first) << std::endl;
      std::cout << "Max square-root of error in " << scheme
		<< " is " << std::sqrt(errs.second) << std::endl;
    }
    return errs;
  }
}

// RMS and L-inf error in vertex normals
std::pair<double,double> 
compute_vnormal_error( const int ref_nrms_hdl,
		       const int cur_nrms_hdl, 
		       const int vbuf_hdl,
		       const int fbuf_hdl, const char *scheme=NULL) {
  static int SURF_normals=0, SURF_integrate=0, SURF_v2f=0, SURF_comparea=0;
  static int BLAS_dot=0, BLAS_mul_scalar=0, BLAS_sub_scalar=0;
  static int BLAS_sum_scalar_MPI=0, BLAS_max_scalar_MPI=0;

  if ( SURF_normals==0) { // Load the functions
    SURF_normals = COM_get_function_handle( "SURF.compute_element_normals");
    SURF_v2f = COM_get_function_handle( "SURF.interpolate_to_centers");
    SURF_integrate = COM_get_function_handle( "SURF.integrate");
    SURF_comparea = COM_get_function_handle( "SURF.compute_element_areas");
  
    BLAS_dot = COM_get_function_handle( "BLAS.dot");
    BLAS_mul_scalar = COM_get_function_handle( "BLAS.mul_scalar");
    BLAS_sub_scalar = COM_get_function_handle( "BLAS.sub_scalar");
    BLAS_sum_scalar_MPI = COM_get_function_handle( "BLAS.sum_scalar_MPI");
    BLAS_max_scalar_MPI = COM_get_function_handle( "BLAS.max_scalar_MPI");
  }
  
  // If two arguments are the same, then initialize ref_nrms
  if ( ref_nrms_hdl == cur_nrms_hdl) {
    COM_call_function( SURF_normals, &ref_nrms_hdl);
    return std::pair<double,double>(0.,0.);
  }
  else {
    MPI_Comm comm = MPI_COMM_WORLD;
  
    // Compute root-mean-square as \sqrt(\int_{2-ref_nrms*cur_nrm))
    double err_rms, err_max, area, minus_two=-2; 
    // Compute area
    COM_call_function( SURF_comparea, &fbuf_hdl);
    COM_call_function( BLAS_sum_scalar_MPI, &fbuf_hdl, &area, &comm);
    
    COM_call_function( BLAS_dot, &cur_nrms_hdl, &ref_nrms_hdl, &vbuf_hdl);
    COM_call_function( BLAS_mul_scalar, &vbuf_hdl, &minus_two, &vbuf_hdl);
    COM_call_function( BLAS_sub_scalar, &vbuf_hdl, &minus_two, &vbuf_hdl);

    COM_call_function( BLAS_max_scalar_MPI, &vbuf_hdl, &err_max, &comm);

    COM_call_function( SURF_v2f, &vbuf_hdl, &fbuf_hdl);
    COM_call_function( SURF_integrate, &fbuf_hdl, &err_rms);
    
    std::pair<double,double> errs(err_rms/area, err_max);
    if ( scheme) {
      std::cout << "Root-mean-square error in " << scheme
		<< " is " << std::sqrt(errs.first) << std::endl;
      std::cout << "Max square-root of error in " << scheme
		<< " is " << std::sqrt(errs.second) << std::endl;
    }
    return errs;
  }
}

// Compute weighted schemes
void compute_weighted_normals( const string &wname) {
  int SURF_compnrm = COM_get_function_handle( "SURF.compute_normals");
  int mesh_handle = COM_get_dataitem_handle( (wname+".mesh").c_str());

  int nrm_ana = COM_get_dataitem_handle( (wname+".vnormals_ana").c_str());
  int vbuf_hdl = COM_get_dataitem_handle( (wname+".verrs").c_str());
  int fbuf_hdl = COM_get_dataitem_handle( (wname+".ferrs").c_str());

  int scheme = SURF::E2N_ONE;
  int nrm_nw = COM_get_dataitem_handle( (wname+".vnormals_nw").c_str());
  COM_call_function( SURF_compnrm, &mesh_handle, &nrm_nw, &scheme);

  compute_vnormal_error( nrm_ana, nrm_nw, vbuf_hdl, fbuf_hdl, "no-weight");

  scheme = SURF::E2N_AREA;
  int nrm_area = COM_get_dataitem_handle( (wname+".vnormals_area").c_str());
  COM_call_function( SURF_compnrm, &mesh_handle, &nrm_area, &scheme);
  compute_vnormal_error( nrm_ana, nrm_area, vbuf_hdl, fbuf_hdl, "area-weighted");

  scheme = SURF::E2N_ANGLE;
  int nrm_angle = COM_get_dataitem_handle( (wname+".vnormals_angle").c_str());
  COM_call_function( SURF_compnrm, &mesh_handle, &nrm_angle, &scheme);
  compute_vnormal_error( nrm_ana, nrm_angle, vbuf_hdl, fbuf_hdl, "angle-weighted");

  scheme = SURF::E2N_SPHERE;
  int nrm_sphere = COM_get_dataitem_handle( (wname+".vnormals_sphere").c_str());
  COM_call_function( SURF_compnrm, &mesh_handle, &nrm_sphere, &scheme);
  compute_vnormal_error( nrm_ana, nrm_sphere, vbuf_hdl, fbuf_hdl, "sphere-weighted");
}

// Compute weighted schemes
void compute_quadric_normals( const string &wname) {
  int PROP_set  = COM_get_function_handle( "PROP.set_option");
  int PROP_propagate = COM_get_function_handle( "PROP.propagate");
  int BLAS_add = COM_get_function_handle( "BLAS.add");
  int nc = COM_get_dataitem_handle( (wname+".nc").c_str());
  int pmesh = COM_get_dataitem_handle( (wname+".pmesh").c_str());
  int spds = COM_get_dataitem_handle( (wname+".spds").c_str());
  int disps = COM_get_dataitem_handle( (wname+".disps").c_str());

  double timestep = 1.e-4;

  int nrm_ana = COM_get_dataitem_handle( (wname+".vnormals_ana").c_str());
  int vbuf_hdl = COM_get_dataitem_handle( (wname+".verrs").c_str());
  int fbuf_hdl = COM_get_dataitem_handle( (wname+".ferrs").c_str());
  int nrm_eig = COM_get_dataitem_handle( (wname+".vnormals_eig").c_str());

  COM_call_function( PROP_set, "rediter", "0");

#if 0
  COM_call_function( PROP_set, "reorthog", "0");

  COM_call_function( PROP_set, "weight", "1");
  COM_call_function( PROP_propagate, &pmesh, &spds, &timestep, &disps);
  compute_vnormal_error( nrm_ana, nrm_eig, vbuf_hdl, fbuf_hdl, "eigen-noweight");

  COM_call_function( PROP_set, "weight", "2");
  COM_call_function( PROP_propagate, &pmesh, &spds, &timestep, &disps);
  compute_vnormal_error( nrm_ana, nrm_eig, vbuf_hdl, fbuf_hdl, "eigen-area");

  COM_call_function( PROP_set, "weight", "3");
  COM_call_function( PROP_propagate, &pmesh, &spds, &timestep, &disps);
  compute_vnormal_error( nrm_ana, nrm_eig, vbuf_hdl, fbuf_hdl, "eigen-angle");

  COM_call_function( PROP_set, "weight", "4");
  COM_call_function( PROP_propagate, &pmesh, &spds, &timestep, &disps);
  compute_vnormal_error( nrm_ana, nrm_eig, vbuf_hdl, fbuf_hdl, "eigen-sphere");
#endif

#if 1
  COM_call_function( PROP_set, "reorthog", "1");
  COM_call_function( PROP_set, "weight", "3");
  COM_call_function( PROP_propagate, &pmesh, &spds, &timestep, &disps);
  COM_call_function( BLAS_add, &nc, &disps, &nc);

  compute_vnormal_error( nrm_ana, nrm_eig, vbuf_hdl, fbuf_hdl, "mq-angle");
  int nrm_mq = COM_get_dataitem_handle( (wname+".vnormals_mq").c_str());
  int BLAS_copy = COM_get_function_handle( "BLAS.copy");
  COM_call_function( BLAS_copy, &nrm_eig, &nrm_mq);
#endif

//   COM_call_function( PROP_set, "reorthog", "0");
//   COM_call_function( PROP_set, "weight", "2");
//   COM_call_function( PROP_propagate, &pmesh, &spds, &timestep, &disps);
//   compute_vnormal_error( nrm_ana, nrm_eig, vbuf_hdl, fbuf_hdl, "eigen-area");
}

void output_solution( const string &wname, const char *timelevel) {
  static int OUT_write = 0, hdl;

  if ( OUT_write==0) {
    OUT_write = COM_get_function_handle( "OUT.write_dataitem");
    hdl = COM_get_dataitem_handle( (wname+".all").c_str());
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
  
  // Initialize the attributes
  init_attributes( wname);
  init_normals( wname, cntr_param.func);

  if ( cntr_param.perturb > 0) {
    int PROP_perturb = COM_get_function_handle( "PROP.perturb_mesh");
    int pmesh = COM_get_dataitem_handle( (wname+".pmesh").c_str());
    COM_call_function( PROP_perturb, &pmesh, &cntr_param.perturb);
  }
  output_solution( wname, "00000");

//   compute_weighted_normals( wname);

  // Initialize control parameters
  init_parameters( cntr_param);
  for ( int i=0; i<cntr_param.steps; ++i) {
    char str[6];
    compute_quadric_normals(wname);
    std::sprintf( str, "%05d", i+1);
    output_solution( wname, str);
  }

  COM_print_profile( "", "Compute-normals");

  COM_finalize();
}






