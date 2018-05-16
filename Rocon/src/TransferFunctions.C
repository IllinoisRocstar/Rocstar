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
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <map>
#include <sstream>
#include <iomanip>
#include <vector>

#include <unistd.h>

//#ifdef _TRAIL_MPI_
#include "mpi.h"
//#else
//typedef int MPI_Comm;
//#endif

#include "GEM.H"
#include "TRAIL.H"
#include "TRAIL_Flu.H"
#include "TRAIL_Remesh.H"

//#ifdef _ROCSTAR_X_
#include "com.h"
//#endif

COM_EXTERN_MODULE( SurfX);
COM_EXTERN_MODULE( SimOUT);
COM_EXTERN_MODULE( SimIN);
COM_EXTERN_MODULE( Simpal);
COM_EXTERN_MODULE( SurfUtil);

using namespace std;



void 
read_file( const string &fname, 
	   const string &wname, 
	   vector<int>  &bcflags,
	   MPI_Comm     comm,
	   bool         apply_disp,
	   bool         all,
	   bool         with_ghost) 
{
  //  char *lastdot=strrchr( fname, '.');

  COM_new_window( wname.c_str());
  // Read in HDF files or a Rocin control file
  std::cout << "Reading file " << fname << "..." << std::endl;

  // Read in HDF format
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimIN, "TRAILIN");
    
  int IN_read;
  // Read in HDF format using SimIN::read_window or ::read_by_control_file 
  IN_read = COM_get_function_handle( "TRAILIN.read_by_control_file");

  // Pass MPI_COMM_NULL to SimIN so that the rank becomes a wildcard.
  //  MPI_Comm comm_null = MPI_COMM_NULL;
  std::string bufwin("bufwin");
  COM_call_function( IN_read, fname.c_str(), bufwin.c_str(), &comm);
  int IN_obtain = COM_get_function_handle( "TRAILIN.obtain_dataitem");

  if(!bcflags.empty()){
    // Check whether bcflag exists. If so, retain only the panes with flag<=1.
    int bcflag = COM_get_dataitem_handle((bufwin+".bcflag").c_str());
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
	if ( flag==NULL )
	  COM_delete_pane( bufwin.c_str(),pane_ids[i]);
	bool delite = true;
	vector<int>::iterator bcfi = bcflags.begin();
	while(bcfi != bcflags.end() && delite)
	  if(*bcfi++ == *flag)
	    delite = false;
	if(delite)
	  COM_delete_pane( bufwin.c_str(), pane_ids[i]);
      }
      // remove buffers.
      COM_free_buffer( &pane_ids);
    }
  }
  if(apply_disp){
    // This is NOT correct for problems with regression.
    int disp_hndl = COM_get_dataitem_handle((bufwin+".uhat").c_str());
    if(disp_hndl > 0){
      std::cout << "Applying total displacements..." << std::endl;
      COM_LOAD_MODULE_STATIC_DYNAMIC( Simpal, "BLAS");
      int add = COM_get_function_handle( "BLAS.add");
      COM_call_function(IN_obtain,&disp_hndl,&disp_hndl);
      int nc_hndl = COM_get_dataitem_handle( bufwin + ".nc");
      COM_call_function( add, &disp_hndl, &nc_hndl, &nc_hndl);
      COM_UNLOAD_MODULE_STATIC_DYNAMIC(Simpal,"BLAS");
    }
  }
  // Read in the mesh.
  int buf_atts;
  if(all)
    buf_atts = COM_get_dataitem_handle((bufwin+".all").c_str());
  else
    buf_atts = COM_get_dataitem_handle((bufwin+".mesh").c_str());
  COM_call_function( IN_obtain, &buf_atts, &buf_atts);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC( SimIN, "TRAILIN");

  if(!all)
    // Remove all dataitems except for the mesh
    COM_delete_dataitem(  (bufwin+".data").c_str());

  
  std::cout << "Obtained window " << wname
	    << " from file " << fname << std::endl;

  // Change the memory layout to contiguous.
  if(all)
    COM_clone_dataitem( (wname+".all").c_str(), 
			 (bufwin+".all").c_str(), (with_ghost ? 1 : 0));
  else
    COM_clone_dataitem( (wname+".mesh").c_str(), 
			 (bufwin+".mesh").c_str(), (with_ghost ? 1 : 0));
  
  COM_delete_window( bufwin.c_str());
}


// A serial program, do not call this on mulitiple procs simultaneously
bool
TRAIL_AutoSurfer(GEM_Partition &gp,
		const string  &src,
		double        t,
		MPI_Comm      comm)
{
  // src will help specify the source window's file names:
  // src_in_<timestring>.txt
  std::string timestring(TRAIL_TimeString(t));
  std::string srcfile(src + "_in_" + timestring + ".txt");
  // the target window's files are specified by 
  // gp.surface_window_in_<timestring>.txt
  std::string trgfile(gp.surface_window+"_in_"+timestring+".txt");
  int rank = 0;
  MPI_Comm_rank(comm,&rank);
  std::string srcwin(src+"_coup");
  std::string trgwin(gp.surface_window+"_coup");
  vector<int> bcflags(2);
  bcflags[0] = 0;
  bcflags[1] = 1;
  const char *format = "HDF";
  

  
  // Serial step to create the common refinement
  if(!rank){
    read_file(trgfile,trgwin,bcflags,MPI_COMM_NULL,false,false,false);
    read_file(srcfile,srcwin,bcflags,MPI_COMM_NULL,false,false,false);

    COM_LOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");
    
    int RFC_overlay = COM_get_function_handle( "RFC.overlay");
    int RFC_write = COM_get_function_handle( "RFC.write_overlay");
    int src_mesh = COM_get_dataitem_handle( (srcwin+".mesh").c_str());
    int trg_mesh = COM_get_dataitem_handle( (trgwin+".mesh").c_str());
  
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_AutoSurfer: Creating overlay for coupled surfaces."
	       << std::endl;
    COM_call_function( RFC_overlay, &src_mesh, &trg_mesh);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_AutoSurfer: Writing overlay for coupled surfaces."
	       << std::endl;
    COM_call_function( RFC_write, &src_mesh, &trg_mesh, 
		       srcwin.c_str(), trgwin.c_str(), format);

    COM_UNLOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");
  }
  

  // Now we have a common refinement, use it to transfer data from the old
  // surface to the new surface
  COM_LOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");
  int RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
  int RFC_read = COM_get_function_handle( "RFC.read_overlay");
  
  string r_trgwin("new_surf_b");
  read_file(trgfile,trgwin,bcflags,comm,false,true,false);
  read_file(srcfile,srcwin,bcflags,comm,false,true,false);
  read_file(trgfile,r_trgwin,bcflags,comm,false,true,true);

  int srcmesh = COM_get_dataitem_handle( (srcwin+".mesh").c_str());
  int trgmesh = COM_get_dataitem_handle( (trgwin+".mesh").c_str());
  std::vector<int> pane_id;
    
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_AutoSurfer: Reading mesh overlay for coupled surfaces." 
	     << std::endl;
  COM_call_function( RFC_read, &srcmesh, &trgmesh, &comm,
		     srcwin.c_str(),trgwin.c_str(),format);
  
  int num_dataitems;
  string names;
  COM_get_dataitems( srcwin.c_str(),&num_dataitems,names);
  char loc;
  COM_Type comtype;
  int ncomp;
  std::string unit;
  istringstream Istr(names);
  for(int i = 0;i < num_dataitems;i++){
    string aname;
    Istr >> aname;
    COM_get_dataitem(srcwin+"."+aname,&loc,&comtype,&ncomp,&unit);
    if((loc == 'e' || loc == 'n') && comtype == COM_DOUBLE){
      if(!rank)
	cout << "Transferring dataitem: " << aname << " on " 
	     << (loc == 'e' ? "elements" : "nodes") << "." << endl; 
      COM_resize_array((srcwin+"."+aname).c_str());
      COM_new_dataitem((trgwin+"."+aname).c_str(),(char)loc,
			COM_DOUBLE,(int)ncomp,unit.c_str());
      COM_new_dataitem((r_trgwin+"."+aname).c_str(),(char)loc,
			COM_DOUBLE,(int)ncomp,unit.c_str());
      COM_resize_array((r_trgwin+"."+aname).c_str());
      COM_resize_array((trgwin+"."+aname).c_str());
      int src_ahdl  = COM_get_dataitem_handle((srcwin+"."+aname).c_str());
      int trg_ahdl  = COM_get_dataitem_handle((trgwin+"."+aname).c_str());
      COM_call_function( RFC_transfer, &src_ahdl, &trg_ahdl);
      int *srcpane_ids;
      int npanes;
      COM_get_panes( trgwin.c_str(), &npanes, &srcpane_ids);
      pane_id.resize(npanes);
      for(int i = 0;i < npanes;i++)
	pane_id[i] = srcpane_ids[i];
      // These are no longer necessary as we've duped the info into
      // a locally allocated array
      COM_free_buffer( &srcpane_ids);
      for(int p = 0;p < npanes;p++){
	void *src_ptr = NULL;
	int src_std = 0;
	int src_cap = 0;
	void *trg_ptr = NULL;
	int trg_std = 0;
	int trg_cap = 0;
	COM_get_array((trgwin+"."+aname).c_str(),pane_id[p],
		      &src_ptr,&src_std,&src_cap);
	COM_get_array((r_trgwin+"."+aname).c_str(),pane_id[p],
		      &trg_ptr,&trg_std,&trg_cap);
	if(src_ptr && trg_ptr && (trg_std*trg_cap == src_std*src_cap))
	  memcpy(trg_ptr,src_ptr,sizeof(double)*src_std*src_cap);
	else
	  if(gp._out)
	    *gp._out << "TRAIL_AutoSurfer: WARNING: non matching sizes for "
		     << aname << " on pane " << pane_id[p] << "." 
		     << std::endl; 
      }
    }
  } 
  
  // Write new surface and SimIN control files for it
  COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT,"Rocout");
  int OUT_set_option = COM_get_function_handle( "Rocout.set_option");
  std::string rankstr("0");
  COM_call_function( OUT_set_option, "rankwidth", rankstr.c_str());
  std::ostringstream Ostr;
  Ostr << r_trgwin << "_" << timestring << "_" << setw(5) 
       << setfill('0') << rank+1;
  int whand = COM_get_function_handle("Rocout.write_dataitem");
  int all = COM_get_dataitem_handle((r_trgwin+".all"));
  COM_call_function(whand,Ostr.str().c_str(),&all,
		    r_trgwin.c_str(),timestring.c_str());
  std::ofstream Ouf;
  string controlfilename(Ostr.str() + "_in.txt");
  Ouf.open(controlfilename.c_str());
  Ouf << "@Proc: " << rank << endl
      << "@Files: " << Ostr.str() << ".hdf" << std::endl;
  Ouf.clear();
  Ouf << "@Panes: ";
  std::vector<int>::iterator pii = pane_id.begin();
  while(pii != pane_id.end())
    Ouf << *pii++ << " ";
  Ouf << endl;
  Ouf.close();
  Ouf.clear();

  // Merge SimIN control files
  if(!rank){
    Ostr.clear();
    Ostr.str("");
    Ostr << r_trgwin << "_" << timestring << "_";
    ostringstream OutS;
    unsigned int id = 1;
    while(id <= gp._npart){
      ifstream Inf;
      ostringstream Ofn;
      Ofn << Ostr.str() << setw(5) << setfill('0')
	  << id++;
      string filename(Ofn.str() + "_in.txt");
      Inf.open(filename.c_str());
      if(Inf){
	OutS << Inf.rdbuf() << endl;
	Inf.close();
	unlink(filename.c_str());
      }
    }
    ofstream Ouf2;
    Ouf2.open((r_trgwin+"_in_"+timestring+".txt").c_str());
    Ouf2 << OutS.str();
    Ouf2.close();
    Ostr.clear();
    Ostr.str("");
  }
  
  COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimOUT,"Rocout");
  COM_UNLOAD_MODULE_STATIC_DYNAMIC(SurfX,"RFC");
  return(true);
}






