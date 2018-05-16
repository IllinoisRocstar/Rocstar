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
// $Id: surfdiver.C,v 1.13 2009/10/08 20:35:10 mtcampbe Exp $

#include <iostream>
#include <cstring>
#include <cstdlib>

#include "roccom.h"

COM_EXTERN_MODULE( Rocface);
COM_EXTERN_MODULE( Rocout);
COM_EXTERN_MODULE( Rocblas);
COM_EXTERN_MODULE( Rocsurf);
COM_EXTERN_MODULE( Rocin);

using namespace std;

void read_file( const char *fname, const string &wname, double alpha) {
  char *lastdot=strrchr( const_cast<char *>(fname), '.');
  COM_new_window( wname.c_str());
  // Read in HDF files or a Rocin control file
  std::cout << "Reading file " << fname << "..." << std::endl;

  // Read in HDF format
  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocin, "IN");
  int IN_read;

  // Read in HDF format using Rocin::read_window or ::read_by_control_file 
  bool processed = false;
  if(!(lastdot == NULL)){
    if ( strcmp( lastdot, ".hdf")==0){
      IN_read = COM_get_function_handle( "IN.read_window");
      processed = true;
    }
  }
  if(!processed){
    IN_read = COM_get_function_handle( "IN.read_by_control_file");
  }

  // Pass MPI_COMM_NULL to Rocin so that the rank becomes a wildcard.
  MPI_Comm comm_null = MPI_COMM_NULL;
  std::string bufwin("bufwin");
  COM_call_function( IN_read, fname, bufwin.c_str(), &comm_null);

  int IN_obtain = COM_get_function_handle( "IN.obtain_attribute");

  // MS
  /*
  int npn, *pnids;
  COM_get_panes( bufwin.c_str(), &npn, &pnids);
  std::cout << "There exisits "<< npn 
            << " panes in this window." << std::endl;
  for ( int i=0; i<npn; ++i) {
     std::cout << "Pane ID = " << pnids[i] << std::endl; 
  }
  */
  // MS End

  // Check whether bcflag exists. If so, retain only the panes with flag<=1.
  int bcflag = COM_get_attribute_handle((bufwin+".bcflag").c_str());
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
      //std::cout << "Existing pane ID = " << pane_ids[i] << std::endl;
      if ( flag==NULL || *flag>1) {
	COM_delete_pane( bufwin.c_str(), pane_ids[i]);
        //std::cout << "Removing..." << std::endl;
      }
    }

    // remove buffers.
    COM_free_buffer( &pane_ids);
  }

  // This is NOT correct for problems with regression.
  int disp_hndl = COM_get_attribute_handle((bufwin+".uhat").c_str());
  if(disp_hndl > 0){
    std::cout << "Applying total displacements..." << std::endl;
    COM_LOAD_MODULE_STATIC_DYNAMIC( Rocblas, "BLAS");
    int add = COM_get_function_handle( "BLAS.add");
    COM_call_function(IN_obtain,&disp_hndl,&disp_hndl);
    int nc_hndl = COM_get_attribute_handle( bufwin + ".nc");
    COM_call_function( add, &disp_hndl, &nc_hndl, &nc_hndl);
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocblas,"BLAS");
  }

  // Read in the mesh.
  int buf_mesh = COM_get_attribute_handle((bufwin+".mesh").c_str());
  COM_call_function( IN_obtain, &buf_mesh, &buf_mesh);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC( Rocin, "IN");

  // Remove all attributes except for the mesh
  COM_delete_attribute(  (bufwin+".atts").c_str());

  
  std::cout << "Obtained window " << wname
	    << " from file " << fname << std::endl;

  // Change the memory layout to contiguous.
  COM_clone_attribute( (wname+".mesh").c_str(), (bufwin+".mesh").c_str(), 0);
  COM_delete_window( bufwin.c_str());
}

int main(int argc, char *argv[]) {
  COM_init( &argc, &argv);
  
  if ( argc < 3) {
    std::cerr << "Usage: " << argv[0]
	      << " <HDF|RocinControlFile1> <HDF|RocinControlFile2> [<out_prefix>] [<RocfaceControlFile>]\n\n"
              << "\t<HDF|RocinControl File1> specifies the files for the first window\n"
              << "\t<HDF|RocinControl File2> specifies the files for the second window\n"
              << "\t<out_prefix> specifies a prefix for output files. \n\t\tDefault is the current directory\n"
	      << "\t<RocfaceControlFile> specifies a file name for Rocface control parameters. \n"
              << "\nExample:\t" 
              << argv[0] << " \"ifluid*.hdf\" \"isolid*.hdf\" " << "\n\t\t"
              << argv[0] << " Rocflo/Rocin/\"ifluid*.txt\" Rocfrac/Rocin/\"isolid*.txt\" Rocman/RocfloRocfrac/" << "\n\t\t"
              << std::endl;
    exit(-1);
  }

  COM_set_profiling( 1);

  string     fnames[2] = {string(argv[1]), string(argv[2])};
  string     pre = (argc>3)?argv[3]:"";
  // Append '/' to pre if not there
  if ( !pre.empty() && pre[pre.size()-1] != '/') pre.append("/");
  
  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocface, "RFC");

  int RFC_readcntr = COM_get_function_handle( "RFC.read_control_file");
  int RFC_overlay = COM_get_function_handle( "RFC.overlay");
  int RFC_write = COM_get_function_handle( "RFC.write_overlay");

  if ( argc>4) {
    std::cout << "Reading Rocface control file..." << std::endl;
    COM_call_function( RFC_readcntr, argv[4]);
    std::cout << "Finished reading Rocface control file." << std::endl;
  }

  string     wnames[2];
  for ( int k=0; k<2; ++k) {
    // Discard the directory name and suffix to obtain a window name.
    string::size_type n0 = fnames[k].find_last_of( "/");

    std::string fname;
    if ( n0 == std::string::npos) fname=fnames[k]; 
    else fname = fnames[k].substr( n0+1, fnames[k].size());

    string::size_type ni;
    ni = fname.find_first_of( ".:_-*[]?\\\"\'0123456789");
    COM_assertion_msg(ni, "File name must start with a letter");

    if ( ni == std::string::npos) {
      wnames[k] = fname;
      fnames[k].append(".hdf"); // Append the .hdf suffix to the file name.
    }
    else {
      if ( fname[ni] == '_' && (fname[ni+1] == 's' || fname[ni+1] == 'f'))
	ni += 2;
      wnames[k] = pre+fname.substr( 0, ni);
    }
    COM_assertion_msg( k==0 || wnames[0]!=wnames[1],
		       "Two input files must have different alphabetic prefix");

    read_file( fnames[k].c_str(), wnames[k], 1.);

    COM_window_init_done( wnames[k].c_str());
  }

  int tri1_mesh = COM_get_attribute_handle( (wnames[0]+".mesh").c_str());
  int tri2_mesh = COM_get_attribute_handle( (wnames[1]+".mesh").c_str());

  //const char *format = "HDF";
  const char *format = "CGNS";
  
  std::cout << "Starting mesh overlay..." << std::endl;
  COM_call_function( RFC_overlay, &tri1_mesh, &tri2_mesh);
  COM_call_function( RFC_write, &tri1_mesh, &tri2_mesh, 
		     wnames[0].c_str(), wnames[1].c_str(), format);

  COM_print_profile( "", "");

  COM_finalize();
}






