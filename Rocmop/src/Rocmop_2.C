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
// $Id: Rocmop_2.C,v 1.5 2008/12/06 08:45:24 mtcampbe Exp $

/****************************************************************
 * Rocmop_2.C (previously named Rocmop_1_1.C)
 * Originally written by Phil Alexandar
 *
 * Modified by Pornput Suriyamongkol - 03/22/07
 * --------------------------------------------
 * - Rename this file to Rocmop_2.C
 *
 * - Make it links to Mesquite 0.9.5
 *   - Disable smooth_boeing() 
 *
 * - Fixes in smooth()
 *   - Check total number of panes
 *   - Coordinate update algorithm
 *   - Displacement update algorithm
 *   - Constrain maximum displacement
 *
 ****************************************************************/

#ifdef MESQUITE
#define USE_STD_INCLUDES 1
#define USE_C_PREFIX_INCLUDES 1
#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "ShapeImprovementWrapper.hpp"

// algorithms
//#include "MeanRatioQualityMetric.hpp"
#include "MeanRatioFunctions.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
//#include "MsqMessage.hpp"
#include "MesqPane_95.h"

using namespace Mesquite;
#endif

#include "Rocblas.h"
#include "Rocmop_2.h"
#include "Rocmap.h"
#include "roccom.h"
#include "Pane_communicator.h"
#include "Pane_connectivity.h"
#include "Pane_ghost_connectivity.h"
#include "Geometric_Metrics_3.h"
#include "PN_patch.h"
#include "geometry.h"
#include "mapbasic.h"

// for debugging
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
// end for debugging

// To enable profiling, use ROCPROF=1 in your make command. 
#ifdef ROCPROF
#include "Rocprof.H"
#endif

MOP_BEGIN_NAMESPACE

using MAP::Pane_connectivity;
using MAP::Pane_communicator;
using MAP::Rocmap;
using MAP::Pane_ghost_connectivity;
using std::cout;
using std::endl;

Rocmop::~Rocmop() { 

  print_legible(0,"Entering Rocmop::~Rocmop");

  if (_buf_window){
    delete _buf_window; 
    _buf_window = NULL;
  }

  for(int i =0, ni=_dcs.size(); i<ni; ++i){
    delete _dcs[i];
    _dcs[i] = NULL;
  }
  print_legible(1,"Exiting Rocmop::~Rocmop");
}

// Utility function to fetch a non-empty, non-commented line from the
// file passed in "Inf". The noop char is the "comment" character, all
// lines beginning with 'noop', or empty lines will be ignored.  All 
// whitespace is also ignored.  Comments may go on the same line as 
// data lines.  Data lines must not start with whitespace.  Valid data
// should be separated by newlines.
std::ifstream &get_next_line(std::ifstream &Inf,
			     std::string &line,
			     const char noop)
{
  if(!Inf)
    return(Inf);
  std::getline(Inf,line);
  while((line.empty() || line[0] == ' ' || line[0] == noop) && Inf)
    std::getline(Inf,line);
  return(Inf);
}

// New Rocmop function to read the configuration file.  Currently, all 
// processors open and read this file (as opposed to process 0 bc'ing).
// If the file does not exist, or has formatting errors, the file is 
// ignored and Rocmop will continue to function with it's defaults.
void Rocmop::read_config_file(const std::string &cfname)
{
  std::ifstream Inf;
  Inf.open(cfname.c_str());

  // We can't use the usual print_legible function here because the
  // buffer winodw may not be created yet.
  int rank =0;
  if(COMMPI_Initialized()){
    int ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    assert( ierr == 0);
  }

  // If the file doesn't exist, fire off a warning and continue
  if(!Inf && rank==0){
    std::cerr << "Rocmop: Warning: Could not open configfile, " 
	      << cfname << "." << std::endl;
    return;
  }

  std::string line;
  std::istringstream Istr;

  // Read and set the verbosity
  int verbosity;
  get_next_line(Inf,line,'#');
  //     At this point, "line" should have data in it.  By using the 
  //     istringstream object this way, we eliminate all whitespace
  //     and extract only the valid part of the line.  If the file 
  //     has already closed, line will be empty, and verbosity will
  //     usually take on some nonsense value.  
  Istr.str(line);
  Istr >> verbosity;
  //     If the line had only the one chunk of data in it, Istr will 
  //     have a bit set that needs to be cleared with a call to "clear()"
  //     before populating it with another string.
  Istr.clear();
  if(verbosity < 0 || verbosity > 10){
    std::cerr << "Rocmop: Warning: Invalid verbosity from configfile. " 
	      << "Giving up on " << cfname << "." << std::endl;
    Inf.close();
    return;
  }
  _verb = verbosity;
  // Speak up if verbosity is turned on!
  if (rank==0 && _verb){
    std::cout << "Rocmop: Found configfile, " << cfname << "." << std::endl
	      << "Rocmop: Setting verbosity to " << _verb << "." << std::endl;
  }  

  // Load the next line from configfile and set smoothing method
  get_next_line(Inf,line,'#');
  Istr.str(line);
  Istr >> _method;
  Istr.clear();
  if(_method < 0 || _method > SMOOTH_NONE){
    std::cerr << "Rocmop: Warning: Invalid smoothing method from configfile. " 
	      << "Giving up on " << cfname << "." << std::endl;
    Inf.close();
    return;
  }
  // I love these ternary operators
  if (rank==0 && _verb){
    std::cout << "Rocmop: Setting method to "
	      << (_method == SMOOTH_VOL_MESQ_WG ? "SMOOTH_VOL_MESQ_WG" :
		  (_method == SMOOTH_VOL_MESQ_NG ? "SMOOTH_VOL_MESQ_NG" :
		   (_method == SMOOTH_SURF_MEDIAL ? "SMOOTH_SURF_MEDIAL" :
		    "SMOOTH_NONE"))) << "." << std::endl;
  }

  // Load the next line and set laziness
  int lazy;
  get_next_line(Inf,line,'#');
  Istr.str(line);
  Istr >> lazy;
  Istr.clear();
  if(lazy > 0) 
    _lazy = 1;
  if (rank==0 && _verb){
    std::cout << "Rocmop: Setting lazy to " << _lazy << "." << std::endl;
  }  

  // Load the next line and set tolerance for lazy threshold
  double tol;
  get_next_line(Inf,line,'#');
  Istr.str(line);
  Istr >> tol;
  Istr.clear();
  if(tol < 0. || tol > 180.){
    std::cerr << "Rocmop: Warning: Invalid dihedral angle tolerance"
	      << " from configfile. " 
	      << "Giving up on " << cfname << "." << std::endl;
    Inf.close();
    return;
  }
  _tol = tol;
  if (rank==0 && _verb){
    std::cout << "Rocmop: Setting tolerance to " << _tol << "." << std::endl;
  }  

  // Load the next line and set node displacement constraint for Mesquite smoothing 
  double max_disp;
  get_next_line(Inf,line,'#');
  Istr.str(line);
  Istr >> max_disp;
  Istr.clear();
  if(max_disp < 0. || max_disp > 10.0){
    std::cerr << "Rocmop: Warning: Invalid displacement constraint from configfile. " 
	      << "Giving up on " << cfname << "." << std::endl;
    Inf.close();
    return;
  }
  _maxdisp = max_disp;
  if (rank==0 && _verb){
    std::cout << "Rocmop: Setting displacement constraint to " 
	      << _maxdisp << "." << std::endl;
  }  

  // Load the next line and set N, which forces Rocmop to smooth every Nth step
  get_next_line(Inf,line,'#');
  Istr.str(line);
  Istr >> _smoothfreq;
  Istr.clear();
  if(_smoothfreq <= 0 || _method == SMOOTH_NONE)
    _smoothfreq = 0;
  if (rank==0 && _verb){
    if(_method == SMOOTH_NONE){
      std::cout << "Rocmop: No method selected, setting N to 0"
		<< ",disabling smoothing." << std::endl;
    }
    else{
      std::cout << "Rocmop: Setting N to " << _smoothfreq 
		<< (_smoothfreq==0 ? ", disabling smoothing." : ".") 
		<< std::endl;
     
    }
  }  
  
  // Displacement threshold.  When the max displacement among all processors 
  // exceeds this value, smooth.
  get_next_line(Inf,line,'#');
  Istr.str(line);
  Istr >> _disp_thresh;
  Istr.clear();
  if(_disp_thresh < 0.0)
    _disp_thresh = 0.0;
  else if ((_smoothfreq > 1) && (_disp_thresh > 0.0) ){
    _smoothfreq = 1;
    if(rank==0 && _verb)
      std::cout << "Rocmop: WARNING: N reset to 1 to enable displacement thresholding."
		<< std::endl;
  }
  if (rank==0 && _verb){
    std::cout << "Rocmop: Setting displacement threshold to " 
	      << _disp_thresh << "." << std::endl;
  }  
  
  // Close up the file, we're done.
  Inf.close();
  
}

void Rocmop::load( const std::string &mname) {

  Rocmop *mop = new Rocmop();

  // Hardcoded the control file location for now.  I didn't really see a clean
  // way to allow this to be dynamically set without adding additional init
  // calls to codes outside Rocmop. 
  mop->read_config_file("Rocmop/RocmopControl.txt");

  COM_new_window( mname.c_str());

  std::string glb=mname+".global";

  COM_new_attribute( glb.c_str(), 'w', COM_VOID, 1, "");
  COM_set_object( glb.c_str(), 0, mop);

  COM_Type types[4];
  types[0] = COM_RAWDATA; 
  types[1] = COM_METADATA;
  types[2] = COM_METADATA;
  types[3] = COM_DOUBLE;

  COM_set_member_function( (mname+".smooth").c_str(), 
			   (Member_func_ptr)(&Rocmop::smooth),
			   glb.c_str(), "bibB", types);

  types[1] = COM_STRING;
  types[2] = COM_VOID;
  COM_set_member_function( (mname+".set_value").c_str(), 
			   (Member_func_ptr)(&Rocmop::set_value),
			   glb.c_str(), "bii", types);    

  types[1] = COM_METADATA;
  types[2] = COM_INT;

  COM_set_member_function( (mname+".smooth_boeing").c_str(),
			   (Member_func_ptr)(&Rocmop::smooth_boeing),
			   glb.c_str(), "bbi", types);

  types[1] = COM_METADATA;
  types[2] = COM_METADATA;

  COM_set_member_function( (mname+".add_aspect_ratios").c_str(),
			   (Member_func_ptr)(&Rocmop::add_aspect_ratios),
			   glb.c_str(), "bbB", types);

  types[2] = COM_DOUBLE;
  types[3] = COM_DOUBLE;

  COM_set_member_function( (mname+".obtain_extremal_dihedrals").c_str(),
			   (Member_func_ptr)(&Rocmop::obtain_extremal_dihedrals),
			   glb.c_str(), "bioo", types);
  
  COM_window_init_done( mname.c_str());

}

void Rocmop::unload( const std::string &mname) {

  Rocmop *mop;
  std::string glb=mname+".global";

  COM_get_object( glb.c_str(), 0, &mop);
  delete mop;

  COM_delete_window( mname.c_str());
}

// Finds the maximum displacement over all processors
bool
Rocmop::check_displacements(COM::Attribute *w_disp)
{
#ifdef ROCPROF
  Profile_begin("Rocmop::check_disp");
#endif
  
  int disp_id = w_disp->id();
  double max_norm = 0.0;
  static double disp_tally = 0.0;
  bool retval = true;
  //  std::ofstream Ouf;
  //  Ouf.open("displacements.txt",ios::app);
  std::vector<Pane*> allpanes;
  const_cast<COM::Window*>(_usr_window)->panes(allpanes);
    
  // Find the maximum displacement on this processor
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    
    Vector_3<double> *ptr = NULL;
    COM::Attribute* ptr_att = allpanes[i]->attribute(disp_id);
    void* void_ptr = ptr_att->pointer();
    ptr = reinterpret_cast<Vector_3<double>*>(void_ptr);

    for(int j=0,nj = allpanes[i]->size_of_real_nodes(); j<nj; ++j){	
      double norm = ptr[j].norm();
      // Re-assign them to zero to maintain original Rocmop behavior
      // just in case.  I'm pretty sure that any values in this array
      // already get wiped out when the Mesquite displacements are 
      // copied in.
      ptr[j][0] = 0.0;
      ptr[j][1] = 0.0;
      ptr[j][2] = 0.0;
      //      Ouf << j << ": " << norm << std::endl;
      if(norm > max_norm)
	max_norm = norm;      
    }
  }
  // Maximize the displacement across all processors
  agree_double(max_norm,MPI_MAX);

  // Update the tally
  disp_tally += max_norm;

  // Build a blurb about displacement
  std::ostringstream Ostr;
  Ostr << "Mesh Displacement: " << disp_tally;

  // Determine whether we have exceeded the threshold
  if(disp_tally  > _disp_thresh){
    disp_tally = 0.0;
    retval = true;
    // If we intend to smooth - also print out the displacement
    print_legible(0,Ostr.str().c_str());
  }
  else {
    retval = false;
    // If not smoothing, print out displacement if verbosity is higher
    print_legible(1,Ostr.str().c_str());
  }
#ifdef ROCPROF
  Profile_end("Rocmop::check_disp");
#endif
  //  Ouf.close();
  return(retval);
}

// Zero displacements array
void
Rocmop::zero_displacements(COM::Attribute *w_disp)
{
#ifdef ROCPROF
  Profile_begin("Rocmop::zero_disp");
#endif
  
  int disp_id = w_disp->id();
  std::vector<Pane*> allpanes;
  const_cast<COM::Window*>(_usr_window)->panes(allpanes);
    
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    Vector_3<double> *ptr = NULL;
    COM::Attribute* ptr_att = allpanes[i]->attribute(disp_id);
    void* void_ptr = ptr_att->pointer();
    ptr = reinterpret_cast<Vector_3<double>*>(void_ptr);
    for(int j=0,nj = allpanes[i]->size_of_real_nodes(); j<nj; ++j){	
      // Zero to maintain original physics module behavior
      ptr[j][0] = 0.0;
      ptr[j][1] = 0.0;
      ptr[j][2] = 0.0;
    }
  }
#ifdef ROCPROF
  Profile_end("Rocmop::zero_disp");
#endif
}

// Add aspect ratio measures to the user's mesh
void Rocmop::add_aspect_ratios(COM::Attribute *usr_att,
			       COM::Attribute *buf_att){
  
  COM::Window* usr_window = usr_att->window();
  COM::Window* buf_window = NULL;

  bool del_buffer = false;

  if(!buf_att){      
    del_buffer = true;
    // Create a buffer window by cloning from the user window
    std::string buf_name(usr_window->name()+"-Rocmop_add_aspect");
    buf_window = new COM::Window(buf_name, 
				 usr_window->get_communicator());
    buf_window->inherit( usr_att, "", 
			 COM::Pane::INHERIT_CLONE, true, NULL,0);
    buf_window->init_done();
  }	      
  else
    buf_window = _buf_window;

  std::vector<Pane*> allpanes;

  buf_window->panes(allpanes);

  // Create and resize window level buffer attributes
#if 0
  COM::Attribute * w_buff_R =
    buf_window->new_attribute("Aspect Ratio:  R (circumradius)",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_buff_r =
    buf_window->new_attribute("Aspect Ratio:  r (inradius)",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_buff_l =
    buf_window->new_attribute("Aspect Ratio:  l (shortest edge)",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_buff_rR =
    buf_window->new_attribute("Aspect Ratio:  3r/R",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_buff_Rr =
    buf_window->new_attribute("Aspect Ratio:  R/r",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_buff_Rl =
    buf_window->new_attribute("Aspect Ratio:  R/l",'e',COM_DOUBLE,1,"");
#endif
  COM::Attribute * w_buff_min =
    buf_window->new_attribute("Aspect Ratio:  min dih",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_buff_max =
    buf_window->new_attribute("Aspect Ratio:  max dih",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_buff_type = 
    buf_window->new_attribute("Element type",'e',COM_INT,1,"");
#if 0  
  buf_window->resize_array(w_buff_R,0);
  buf_window->resize_array(w_buff_r,0);
  buf_window->resize_array(w_buff_l,0);
  buf_window->resize_array(w_buff_rR,0);
  buf_window->resize_array(w_buff_Rr,0);
  buf_window->resize_array(w_buff_Rl,0);
#endif
  buf_window->resize_array(w_buff_min,0);
  buf_window->resize_array(w_buff_max,0);
  buf_window->resize_array(w_buff_type,0);
  buf_window->init_done();

  // Obtain buffer attribute ids
#if 0
  int R_id   = w_buff_R->id();
  int r_id   = w_buff_r->id();
  int l_id   = w_buff_l->id();
  int rR_id  = w_buff_rR->id();
  int Rr_id  = w_buff_Rr->id();
  int Rl_id  = w_buff_Rl->id();
#endif
  int min_id = w_buff_min->id();
  int max_id = w_buff_max->id();
  int type_id = w_buff_type->id();

   // Create and resize window level user attributes
#if 0
  COM::Attribute * w_usr_R =
    usr_window->new_attribute("Aspect Ratio:  R",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_usr_r =
    usr_window->new_attribute("Aspect Ratio:  r",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_usr_l =
    usr_window->new_attribute("Aspect Ratio:  l",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_usr_rR =
    usr_window->new_attribute("Aspect Ratio:  3r/R",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_usr_Rr =
    usr_window->new_attribute("Aspect Ratio:  R/r",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_usr_Rl =
    usr_window->new_attribute("Aspect Ratio:  R/l",'e',COM_DOUBLE,1,"");
#endif
  COM::Attribute * w_usr_min =
    usr_window->new_attribute("Aspect Ratio:  min dih",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_usr_max =
    usr_window->new_attribute("Aspect Ratio:  max dih",'e',COM_DOUBLE,1,"");
  COM::Attribute * w_usr_type = 
    usr_window->new_attribute("Element type",'e',COM_INT,1,"");
  
#if 0
  usr_window->resize_array(w_usr_R,0);
  usr_window->resize_array(w_usr_r,0);
  usr_window->resize_array(w_usr_l,0);
  usr_window->resize_array(w_usr_rR,0);
  usr_window->resize_array(w_usr_Rr,0);
  usr_window->resize_array(w_usr_Rl,0);
#endif
  usr_window->resize_array(w_usr_min,0);
  usr_window->resize_array(w_usr_max,0);
  usr_window->resize_array(w_usr_type,0);
  usr_window->init_done();

#if 0    
  double R,r,l;
#endif
  double angles[] = {0.0,0.0};

  for(int i=0,ni=(int)allpanes.size();i<ni;++i){
#if 0
    // get Pane level buffer attributes
    double *R_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(R_id)->pointer());    
    double *r_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(r_id)->pointer());    
    double *l_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(l_id)->pointer());    
    double *rR_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(rR_id)->pointer());    
    double *Rr_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(Rr_id)->pointer());    
    double *Rl_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(Rl_id)->pointer());    
#endif
    double *min_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(min_id)->pointer());    
    double *max_ptr = reinterpret_cast<double *>
      (allpanes[i]->attribute(max_id)->pointer());    
    int * type_ptr = reinterpret_cast<int *>
      (allpanes[i]->attribute(type_id)->pointer());
    
    for(int j=0,nj=allpanes[i]->size_of_elements();j<nj;++j){
      
      Element_node_enumerator ene(allpanes[i],j+1);
      
#if 0
      Aspect_Metric_3 aspect_met;	
      aspect_met.initialize(ene);
      aspect_met.getAspects(R,r,l);
#endif
      Angle_Metric_3 angle_met;
      angle_met.initialize(ene);
      angle_met.compute(angles);
#if 0     
      R_ptr[j] = R;
      r_ptr[j] = r;
      l_ptr[j] = l;
      COM_assertion_msg(R != 0.0,
			"Cannot calculate aspect ratios with zero magnitude circumradius");
      rR_ptr[j] = (3.0*r)/R;
      COM_assertion_msg(r != 0.0,
			"Cannot calculate aspect ratios with zero magnitude circumradius");
      Rr_ptr[j] = R/r;
      COM_assertion_msg(l != 0.0,
			"Cannot calculate aspect ratios with zero magnitude circumradius");
      Rl_ptr[j] = R/l;
#endif
      min_ptr[j] = angles[0];
      max_ptr[j] = angles[1];
      type_ptr[j] = ene.type();
    }
  }
#if 0
  Rocmap::update_ghosts(w_buff_R);
  Rocmap::update_ghosts(w_buff_r);
  Rocmap::update_ghosts(w_buff_l);
  Rocmap::update_ghosts(w_buff_rR);
  Rocmap::update_ghosts(w_buff_Rr);
  Rocmap::update_ghosts(w_buff_Rl);
#endif
  Rocmap::update_ghosts(w_buff_min);
  Rocmap::update_ghosts(w_buff_max);
  Rocmap::update_ghosts(w_buff_type);
#if 0
  Rocblas::copy(w_buff_R,w_usr_R);
  Rocblas::copy(w_buff_r,w_usr_r);
  Rocblas::copy(w_buff_l,w_usr_l);
  Rocblas::copy(w_buff_rR,w_usr_rR);
  Rocblas::copy(w_buff_Rr,w_usr_Rr);
  Rocblas::copy(w_buff_Rl,w_usr_Rl);
#endif
  Rocblas::copy(w_buff_min,w_usr_min);
  Rocblas::copy(w_buff_max,w_usr_max);
  Rocblas::copy(w_buff_type,w_usr_type);

  if(del_buffer)
    delete buf_window;
}

#if 0
void Rocmop::get_user_coords(){
  std::vector<Pane*> allbufpanes;
  std::vector<Pane*> allusrpanes;
  const_cast<COM::Window*>(_usr_window)->panes(allusrpanes);
  _buf_window->panes(allbufpanes);
  COM_assertion_msg(allbufpanes.size() == allusrpanes.size(),
		    "Buffer window and user window have different number of panes.");
  
  double *usr_x = NULL;
  double *usr_y = NULL;
  double *usr_z = NULL;

  int x_space = 1;
  int y_space = 1;
  int z_space = 1;
  
  for(int i=0, ni = allusrpanes.size(); i<ni; ++i){

    int stride = allusrpanes[i]->attribute(COM_NC)->stride();
    // if stride == 0, then X,Y,Z coords stored in separate arrays
    if(stride>=3){
      x_space = y_space = z_space = allusrpanes[i]->attribute(COM_NC)->
	size_of_components();
      usr_x = allusrpanes[i]->attribute(COM_NC)->coordinates();
      usr_y = usr_x +1;
      usr_z = usr_y +1;
    }
    // else, if stride == 1, all X coords, then all Y, then all Z in 1 array
    else if(stride==1){
    }
  }
}

void Rocmop::get_real_disp(){
}
#endif

void Rocmop::smooth(const COM::Attribute *pmesh,
		    COM::Attribute *disp,
		    double* timestep){

  // The static var N will statically store the current value of 
  // N.  When N == 0, we will smooth and set N back to (_smoothfreq - 1).
  // What this amounts to is that we smooth only every _smoothfreqth call
  // to Rocmop.  _smoothfreq is specified as "N" in the config file.
  static int N = (_smoothfreq-1);
  
  if(N == 0){ // Go ahead and actually smooth if this is true
    N = (_smoothfreq-1); // Reset N 
    
    // We want to start profiling from here as it's smoothing we
    // want to time, not the NOOP calls.
#ifdef ROCPROF
    Profile_begin("Rocmop::smooth");
#endif
    
    COM_assertion_msg( validate_object()==0, "Invalid object");
    
    // Verify that we are working with a mesh or pmesh
    int pmesh_id = pmesh->id();
    COM_assertion_msg( pmesh &&
		       (pmesh_id==COM::COM_PMESH || pmesh_id==COM::COM_MESH),
		       "Input to Rocmop::smooth must be a mesh or pmesh attribute.");
    _is_pmesh = (pmesh_id == COM_PMESH) ? true : false;
    
    // If buffer window exists and was cloned from the current user
    // window, then just update the coordinates.
   
    if(_buf_window && (_usr_window == pmesh->window()))
    {
      // update nodel coordinates (real part) from _usr_window
      update_buf_real_nc();    
    }
    
    // Else, inherit(clone) a buffer window from the user's mesh.
    else
    {
      // Eliminate any old buffers
      if(_buf_window)
	delete _buf_window;
      
      _usr_window = pmesh->window();

      // First, check if pconn of _usr_window is complete
      int pconnErr = 0;
      pconnErr = check_input_pconn();
      agree_int(pconnErr, MPI_MAX);
      if (pconnErr == 1)
	COM_assertion_msg(pconnErr == 0, "*** Rocmop:ERROR: PCONN incomplete ***");

      // Create a buffer window by cloning from the user window
      std::string buf_name(_usr_window->name()+"-Rocmopbuf");
      _buf_window = new COM::Window(buf_name, 
				    _usr_window->get_communicator());
      _buf_window->inherit( const_cast<COM::Attribute*>(_usr_window->attribute(COM::COM_MESH)), "", 
			    COM::Pane::INHERIT_CLONE, false, NULL, 0);
      _buf_window->init_done();

      // Add the pconn
      if(_buf_window->size_of_panes_global() > 1){
	Pane_ghost_connectivity pgc(_buf_window);
	pgc.build_pconn();
      }
     
      // Perform initialization specific to the smoother.
      smoother_specific_init();
    }
   
    // Max possible worst element quality is 180.0, so this
    // default should result in mesh smoothing if _lazy==0
    double mesh_qual = 190.0;
    
    // Check initial quality if we are in lazy mode
    if(_lazy){
      
      // Fill a vector with pointers to the local panes
      std::vector<const Pane*> allpanes;
      _buf_window->panes(allpanes);
      mesh_qual = check_all_elem_quality(allpanes);
    }
    
    // Mechanism to trigger smoothing on tallied displacements
    bool exceeded = true;
    if(_disp_thresh > 0.0)
      exceeded = check_displacements(disp);
    
    if(exceeded || ((mesh_qual > _tol) && _lazy)){     
      print_legible(0,"Smoothing...");
      perform_smoothing();    
      print_legible(0,"Smoothing complete.");

      // At this point NC in _buf_window has been smoothed.
      // Get displacement.
      get_usr_disp(pmesh, disp, timestep);
    }
    else
    {
      _usr_window = pmesh->window();
      zero_displacements(disp);
    }
    
#ifdef ROCPROF
    Profile_end("Rocmop::smooth");
#endif
  } // if (N says it's smoothing time
  else{
    if(_smoothfreq > 0)
      N--;
    _usr_window = pmesh->window();
    zero_displacements(disp);
  }
}

/*************************************************************************
 *
 * Rocmop::update_buf_real_nc()
 *
 * Update update real NC of _buf_window from _usr_window
 *
 *************************************************************************/

void Rocmop::update_buf_real_nc()
{
  std::vector<const Pane*> allbufpanes;
  _buf_window->panes(allbufpanes);   
  std::vector<const Pane*> allusrpanes;
  _usr_window->panes(allusrpanes);

  COM_assertion_msg(allbufpanes.size() == allusrpanes.size(),
		    "Different number of panes on buffer and user windows.");

  for(int i = 0, ni = allbufpanes.size(); i < ni; ++i)
  {
    // Obtain pane level attributes
    const COM::Attribute *usr_nc = allusrpanes[i]->attribute( COM::COM_NC);
    const COM::Attribute *buf_nc = allbufpanes[i]->attribute( COM::COM_NC);
    
    COM_assertion_msg(usr_nc->stride() != 0,
		      "Rocmop can not operate on meshes with stride == 0");

    COM_assertion_msg(usr_nc->size_of_real_items() == buf_nc->size_of_real_items(),
		      "Number of real items differs between buffer and user windows");

    double* usr_nc_ptr = (double*)usr_nc->pointer();
    double* buf_nc_ptr = (double*)buf_nc->pointer();

    // get stride and distance_to_next_component
    int usr_nc_stride = usr_nc->stride();
    int buf_nc_stride = buf_nc->stride();
    int usr_nc_next_comp;
    int buf_nc_next_comp;
	  
    if (usr_nc_stride == 1)
      usr_nc_next_comp = usr_nc->size_of_items();
    else
      usr_nc_next_comp = 1;
	  
    if (buf_nc_stride == 1)
      buf_nc_next_comp = buf_nc->size_of_items();
    else
      buf_nc_next_comp = 1;
	
    // update real node coordinates
    for(int j=0, nj = usr_nc->size_of_real_items(); j<nj; ++j){
      *buf_nc_ptr                       = *usr_nc_ptr;
      *(buf_nc_ptr+ 1*buf_nc_next_comp) = *(usr_nc_ptr + 1*usr_nc_next_comp);
      *(buf_nc_ptr+ 2*buf_nc_next_comp) = *(usr_nc_ptr + 2*usr_nc_next_comp);
      
      buf_nc_ptr += buf_nc_stride;
      usr_nc_ptr += usr_nc_stride;
    }
  }
}


/*************************************************************************
 *
 * Rocmop::check_input_pconn()
 *
 * Check if all ghost nodes are listed in pconn block 3 of _usr_window
 *
 *************************************************************************/

int Rocmop::check_input_pconn()
{
  int rank =0;
  if(COMMPI_Initialized())
  {
    int ierr = MPI_Comm_rank(_usr_window->get_communicator(), &rank); 
    assert( ierr == 0);
  }

  std::vector<const Pane*> allusrpanes;
  _usr_window->panes(allusrpanes);

  for(int i = 0, ni = allusrpanes.size(); i < ni; ++i)
  {
    // Obtain pane level attributes
    const COM::Attribute *usr_nc = allusrpanes[i]->attribute(COM::COM_NC);
    const COM::Attribute *usr_pconn = allusrpanes[i]->attribute(COM::COM_PCONN);
    int minGhostNodeID = usr_nc->size_of_real_items() + 1;
    int maxGhostNodeID = usr_nc->size_of_items();
    int numGhostNodes  = maxGhostNodeID - minGhostNodeID + 1;
    int * pconnArray = (int *)usr_pconn->pointer();

    // Check pconn
    bool OORGhostNodes  = false;    // found out of range ghost node
    bool missGhostNodes = false;    // ghost nodes missed from pconn

    // Allocate marking array
    int * isMarked = new int [numGhostNodes];
    for (int j = 0; j < numGhostNodes; j++)
      isMarked[j] = 0;

    // Jump to pconn block 3
    int index = 0;
    int numNbrs = 0;
    int numItems = 0;

    for (int j = 0; j < 2; j++)
    {
      numNbrs = pconnArray[index];
      index++;

      for (int k = 0; k < numNbrs; k++)
      {
	index++;                         // step over proc id
	numItems = pconnArray[index];
	index++;
	index += numItems;
      }
    }

    // Begin marking
    numNbrs = pconnArray[index];
    index++;

    for (int j = 0; j < numNbrs; j++)
    {
      index++;                         // step over proc id
      numItems = pconnArray[index];
      index++;

      for (int k = 0; k < numItems; k++)
      {
	int lid = pconnArray[index];
	index++;

	if (lid < minGhostNodeID || lid > maxGhostNodeID)
	{
	  std::cerr << "Rocmop: Error in check_input_pconn(): Rank "
		    << rank << ", Ghost node " << lid 
		    << " is out of ghost node range.\n" << std::endl;
	  OORGhostNodes = true;
	}
	else if (isMarked[lid - minGhostNodeID] == 1)
	{
	  if (_verb > 0)
	    std::cerr << "Rocmop: Warning in check_input_pconn(): Rank "
		      << rank << ", Ghost node " << lid 
		      << " has > 1 owners.\n" << std::endl;
	}
	else
	  isMarked[lid - minGhostNodeID] = 1;
      }
    }

    // Check if any node is not marked
    for (int j = 0; j < numGhostNodes; j++)
      if (isMarked[j] == 0)
      {
	missGhostNodes = true;
	std::cerr << "Rocmop: Error in check_input_pconn(): Rank "
		  << rank << ", Ghost node " << (j + minGhostNodeID)
		  << " is not listed in pconn block 3.\n" << std::endl;
      }

    // Free up memory
    delete isMarked;

    // Return the result
    if (OORGhostNodes == true || missGhostNodes == true)
      return 1;
  }

  return 0;
}


/*************************************************************************
 *
 * Rocmop::get_usr_disp()
 *
 * Get displacement in _usr_window based on
 * nodal coordinates in _buf_window
 *
 *************************************************************************/

void Rocmop::get_usr_disp(const COM::Attribute *pmesh,
			  COM::Attribute *disp,
			  double* timestep)
{
  /*
   *    Disp array contains all nodal displacements. 
   *    For the real part, all displacements must not exceed 
   *    a given constraint. For the ghost part, they can since 
   *    they need to make up for incorrect initial positions 
   *    given by Rocflu (Rocflu never updates ghost nodes).
   *    
   *    The following steps are done.
   *
   * 1. Store displacements of REAL nodes in disp
   *
   * 2. Check if any displacements exceed the constraint
   *    2.1 If so, scale down all REAL displacements
   *
   * 3. Put updated NC in disp
   *
   * 4. Update ghost NC using Rocmap::update_ghosts(disp)
   *
   * 5. Disp now can be correctly calculated using 
   *    Rocblas::sub(disp, orig_nc, disp)
   *
   *
   *    Pornput Suriyamongkol 
   *         02/20/2007
   */
  
  std::vector<const Pane*> allusrpanes;
  std::vector<const Pane*> allbufpanes;
  _usr_window->panes(allusrpanes);
  _buf_window->panes(allbufpanes);
  COM_assertion_msg(allbufpanes.size() == allusrpanes.size(),
		    "Different number of panes on buffer and user windows.");
  
  // 1.  Store REAL displacment in disp 
  
  int disp_id = disp->id();
  for(int i = 0, ni = allusrpanes.size(); i < ni; ++i){
    
    const COM::Attribute *old_nc = allusrpanes[i]->attribute(COM::COM_NC);
    const COM::Attribute *new_nc = allbufpanes[i]->attribute(COM::COM_NC);
    COM::Attribute *p_disp       = const_cast<COM::Attribute*>
      (allusrpanes[i]->attribute(disp_id));
    
    COM_assertion_msg((p_disp->size_of_real_items() == new_nc->size_of_real_items()) &&
		      (p_disp->size_of_real_items() == old_nc->size_of_real_items()),
		      "Number of real items differs between buffer and user windows");
    
    double* old_nc_ptr = (double*)old_nc->pointer();
    double* new_nc_ptr = (double*)new_nc->pointer();
    double* p_disp_ptr = (double*)p_disp->pointer();
    
    int old_nc_stride = old_nc->stride();
    int new_nc_stride = new_nc->stride();
    int p_disp_stride = p_disp->stride();
    
    // Set distances to the next component for each data structure
    int old_nc_next_comp;
    int new_nc_next_comp;
    int p_disp_next_comp;
    
    if (old_nc_stride == 1)
      old_nc_next_comp = old_nc->size_of_items();
    else
      old_nc_next_comp = 1;
    
    if (new_nc_stride == 1)
      new_nc_next_comp = new_nc->size_of_items();
    else
      new_nc_next_comp = 1;
    
    if (p_disp_stride == 1)
      p_disp_next_comp = p_disp->size_of_items();
    else
      p_disp_next_comp = 1;
    
    // Put displacements for real nodes in p_disp
    for(int j = 0, nj = new_nc->size_of_real_items(); j < nj; ++j){
      
      *p_disp_ptr                      = *new_nc_ptr                      - *old_nc_ptr;
      *(p_disp_ptr+1*p_disp_next_comp) = *(new_nc_ptr+1*new_nc_next_comp) - *(old_nc_ptr+1*old_nc_next_comp);
      *(p_disp_ptr+2*p_disp_next_comp) = *(new_nc_ptr+2*new_nc_next_comp) - *(old_nc_ptr+2*old_nc_next_comp);
      
      // go to next item
      old_nc_ptr += old_nc_stride;
      new_nc_ptr += new_nc_stride;
      p_disp_ptr += p_disp_stride;
    }
  }
  
  // 2. Constrain displacements
  
  constrain_displacements(disp);
  
  // 3. Update real NC in disp
  
  for(int i = 0, ni = allusrpanes.size(); i < ni; ++i){
    
    const COM::Attribute *old_nc = allusrpanes[i]->attribute(COM::COM_NC);
    COM::Attribute *p_disp       = const_cast<COM::Attribute*>
      (allusrpanes[i]->attribute(disp_id));
    
    double* old_nc_ptr = (double*)old_nc->pointer();
    double* p_disp_ptr = (double*)p_disp->pointer();
    
    int old_nc_stride = old_nc->stride();
    int p_disp_stride = p_disp->stride();
    
    // Set distances to the next component for each data structure
    
    int old_nc_next_comp;
    int p_disp_next_comp;
    
    if (old_nc_stride == 1)
      old_nc_next_comp = old_nc->size_of_items();
    else
      old_nc_next_comp = 1;
    
    if (p_disp_stride == 1)
      p_disp_next_comp = p_disp->size_of_items();
    else
      p_disp_next_comp = 1;
    
    // Update real part of NC in disp
    for(int j = 0, nj = p_disp->size_of_real_items(); j < nj; ++j){
      
      *p_disp_ptr                      = *old_nc_ptr                      + *p_disp_ptr;
      *(p_disp_ptr+1*p_disp_next_comp) = *(old_nc_ptr+1*old_nc_next_comp) + *(p_disp_ptr+1*p_disp_next_comp);
      *(p_disp_ptr+2*p_disp_next_comp) = *(old_nc_ptr+2*old_nc_next_comp) + *(p_disp_ptr+2*p_disp_next_comp);
      
      // go to next item
      old_nc_ptr += old_nc_stride;
      p_disp_ptr += p_disp_stride;
    }
  }
  
  // 4. Update ghost part of disp 
  
  Rocmap::update_ghosts(disp);  
  
  // 5. Get displacment
  
  const COM::Attribute *orig_nc = pmesh->window()->attribute(COM::COM_NC);
  Rocblas::sub(disp, orig_nc, disp);
  
  if(timestep)
    Rocblas::div_scalar (disp,(const void*)timestep, disp);
}



void Rocmop::perform_smoothing(){
#ifdef ROCPROF
  Profile_begin("Rocmop::perform_smoothing");
#endif
  print_legible(1,"  Entering Rocmop::perform_smoothing");

  COM_assertion_msg(_buf_window, "Unexpected NULL pointer encountered.");

  // Select iterative or noniterative option depending on the smoother.
  if (0) // No noniterative smoothers currently
    perform_noniterative_smoothing();
  else if( _method < SMOOTH_NONE)
    perform_iterative_smoothing();
  else
    COM_assertion_msg(0, "No valid smoothing method selected");

  print_legible(1,"  Exiting Rocmop::perform_smoothing");
#ifdef ROCPROF
  Profile_end("Rocmop::perform_smoothing");
#endif
}

void Rocmop::perform_iterative_smoothing(){
  // All internal iterative behavior has been removed.
#ifdef ROCPROF
  Profile_begin("Rocmop::perform_itersmooth");
#endif

  print_legible(1,"    Entering Rocmop::perform_iterative_smoothing");

  COM_assertion_msg(_buf_window, "Unexpected NULL pointer encountered.");

  std::vector<const Pane*> allpanes;
  _buf_window->panes(allpanes);
  
#ifdef MESQUITE

  if(_method==SMOOTH_VOL_MESQ_WG)
    {
      smooth_vol_mesq_wg();
    }
  else if(_method==SMOOTH_VOL_MESQ_NG)
    {
      smooth_vol_mesq_ng();
    }
  
#else

  if((_method==SMOOTH_VOL_MESQ_WG) || (_method==SMOOTH_VOL_MESQ_NG))
    COM_assertion_msg(0,"Rocmop not compiled with MESQUITE");

#endif
  
  else COM_assertion_msg(0, "No valid iterative smoothing method selected");  
  
  print_legible(1,"    Exiting Rocmop::perform_iterative_smoothing");
#ifdef ROCPROF
  Profile_end("Rocmop::perform_itersmooth");
#endif
}

void Rocmop::perform_noniterative_smoothing(){
#ifdef ROCPROF
  Profile_begin("Rocmop::perform_nonitersmooth");
#endif

  print_legible(1,"    Entering Rocmop::perform_noniterative_smoothing");

  if(_niter != 1){
    std::cerr << "Although the maximum number of iterations > 1 is selected,\n"
	      << "the smoothing method is noniterative, and will run once.\n";
  }

  if (0)
    ;// Currently, no noniterative methods exist.
  else COM_assertion_msg(0, "No valid noniterative smoothing method selected");

  print_legible(1,"    Exiting Rocmop::perform_noniterative_smoothing");
#ifdef ROCPROF
  Profile_end("Rocmop::perform_nonitersmooth");
#endif
}

// Perform smoother specific initializing, for example initializing the
// Window_manifold_2 for a surface mesh, adding smoother specific attributes
// to the window, etc.
void Rocmop::smoother_specific_init(){
#ifdef ROCPROF
  Profile_begin("Rocmop::perform_smooth_init");
#endif

  print_legible(1,"    Entering Rocmop::smoother_specific_init");

  // get rid of any old data
  if(_wm){ delete _wm; _wm = NULL; }

  // perform initialization common to surface smoothers
  if(_method==SMOOTH_SURF_MEDIAL){
    // Initialize the Window_manifold
    if(_wm == NULL)
      Rocsurf::initialize(_buf_window->attribute(_is_pmesh ? COM_PMESH : COM_MESH));    
  }
  
  // perform smoother specific initialization 
  switch (_method){
    
#ifdef MESQUITE
  case SMOOTH_VOL_MESQ_WG: {
    // Obtain a list of elements containing shared nodes for each pane.
    determine_shared_border();
    if(_invert_tets){
      invert_elements(COM::Connectivity::TET4);
      _invert_tets = 0;
    }
    if(_invert_hexes){
      invert_elements(COM::Connectivity::HEX8);
      _invert_hexes = 0;
    }
    break;
  }
  case SMOOTH_VOL_MESQ_NG: {

    // Check to see if the physical surface boundary exists
    const std::string surf_attr("is_surface");
    const COM::Attribute *w_is_surface = _usr_window->attribute(surf_attr);
    if(w_is_surface){
      
      COM_assertion_msg( COM_compatible_types( w_is_surface->data_type(), COM_INT),
			 "Surface-list must have integer type");
      COM_assertion_msg( w_is_surface->is_nodal() == 1,
			 "Surface-list must be nodal");
      COM_assertion_msg( w_is_surface->size_of_components() == 1,
			 "Surface-list must have a single component");
      COM_assertion_msg( w_is_surface->initialized() == 1,
			 "Surface-list must be initialized");      

      // Clone the attribute
      COM::Attribute * new_attr = 
	_buf_window->inherit( const_cast<COM::Attribute *>(w_is_surface), 
			      surf_attr, COM::Pane::INHERIT_CLONE, true, NULL, 0);
      COM_assertion_msg(new_attr, "Attribute could not be allocated.");
    }
    // else, detect the physical boundary ourselves
    else{
      COM::Attribute* w_surf_attr =
	_buf_window->new_attribute( "is_surface", 'n', COM_INT, 1, "");
      _buf_window->resize_array( w_surf_attr, 0);
      
      determine_physical_border(w_surf_attr);
    }
    _buf_window->init_done();
    
    if(_invert_tets)
      invert_elements(COM::Connectivity::TET4);
    break;
  }
#else
  case SMOOTH_VOL_MESQ_WG:
  case SMOOTH_VOL_MESQ_NG:
    COM_assertion_msg(0, "Not compiled with MESQUITE");
    break;
#endif

  case SMOOTH_SURF_MEDIAL: {

    // Extend buffer window
    COM::Attribute* w_disps =
      _buf_window->new_attribute( "disps", 'n', COM_DOUBLE, 3, "");
    _buf_window->resize_array( w_disps, 0);

    COM::Attribute* w_facenormals = 
      _buf_window->new_attribute( "facenormals", 'e', COM_DOUBLE, 3, "");
    _buf_window->resize_array( w_facenormals, 0);
    
    COM::Attribute* w_facecenters = 
      _buf_window->new_attribute( "facecenters", 'e', COM_DOUBLE, 3, "");
    _buf_window->resize_array( w_facecenters, 0);

    COM::Attribute* w_eigvalues = 
    _buf_window->new_attribute( "lambda", 'n', COM_DOUBLE, 3, "");
    _buf_window->resize_array( w_eigvalues, 0);

    COM::Attribute* w_vnormals = 
    _buf_window->new_attribute( "vnormals", 'n', COM_DOUBLE, 3, "");
    _buf_window->resize_array( w_vnormals, 0);

    COM::Attribute* w_awnormals = 
    _buf_window->new_attribute( "awnormals", 'n', COM_DOUBLE, 3, "");
    _buf_window->resize_array( w_awnormals, 0);

    COM::Attribute* w_uwnormals = 
    _buf_window->new_attribute( "uwnormals", 'n', COM_DOUBLE, 3, "");
    _buf_window->resize_array( w_uwnormals, 0);

    COM::Attribute* w_eigvecs = 
    _buf_window->new_attribute( "eigvecs", 'n', COM_DOUBLE, 9, "");
    _buf_window->resize_array( w_eigvecs, 0);

    COM::Attribute* w_tangranks = 
    _buf_window->new_attribute( "tangranks", 'n', COM_INT, 1, "");
    _buf_window->resize_array( w_tangranks, 0);

    COM::Attribute* w_cntnranks = 
    _buf_window->new_attribute( "cntnranks", 'n', COM_INT, 1, "");
    _buf_window->resize_array( w_cntnranks, 0);
    
    COM::Attribute* w_cntnvecs =
    _buf_window->new_attribute( "cntnvecs", 'n', COM_DOUBLE, 6, "");
    _buf_window->resize_array( w_cntnvecs, 0);

    COM::Attribute* w_scales = 
    _buf_window->new_attribute( "scales", 'n', COM_DOUBLE, 1, "");
    _buf_window->resize_array( w_scales, 0);
    
    COM::Attribute* w_weights = 
    _buf_window->new_attribute( "weights", 'n', COM_DOUBLE, 1, "");
    _buf_window->resize_array( w_weights, 0);

    COM::Attribute* w_weights2 = 
    _buf_window->new_attribute( "weights2", 'n', COM_DOUBLE, 1, "");
    _buf_window->resize_array( w_weights2, 0);

    COM::Attribute* w_barycrds = 
    _buf_window->new_attribute( "barycrds", 'n', COM_DOUBLE, 2, "");
    _buf_window->resize_array( w_barycrds, 0);

    COM::Attribute* w_PNelemids = 
    _buf_window->new_attribute( "PNelemids", 'n', COM_INT, 1, "");
    _buf_window->resize_array( w_PNelemids, 0);

    // Extend the buffer window to hold local contributions to new placement
    // and the number of contributing faces for Laplacian smoothing.
    COM::Attribute * w_pnt_contrib = 
      _buf_window->new_attribute("pnt_contrib", 'n', COM_DOUBLE, 3, "");
    _buf_window->resize_array(w_pnt_contrib, 0);

    COM::Attribute * w_disp_count =
      _buf_window->new_attribute("disp_count", 'n', COM_DOUBLE, 1, "");
    _buf_window->resize_array(w_disp_count, 0);

    _buf_window->init_done();

    break;
  }
    //  case SMOOTH_LAPLACE : {
    //break;
    //}
  default:
    COM_assertion_msg(0, "Can't initialize for invalid smoother.");
    break;
  }

  COM_assertion_msg(_buf_window, "Unexpected NULL pointer encountered.");

  print_legible(1,"    Exiting Rocmop::smoother_specific_init");
#ifdef ROCPROF
  Profile_end("Rocmop::perform_smooth_init");
#endif
}

void Rocmop::set_value(const char* opt, const void* value)
{
  
  //print_legible(0,"Entering Rocmop::set_value");

  COM_assertion_msg( validate_object()==0, "Invalid object");

  COM_assertion_msg( opt && value,
		      "Rocmop::set_value does not except NULL parameters");
  std::string option;
  if(opt) option = opt;
  if ( option == "method") {
    COM_assertion_msg( *((int*)value) <= SMOOTH_NONE && *((int*)value)>=0 
                      ,"Illegal value for 'method' option");
    _method = *((int*)value);
  }
  else if ( option == "wrapper"){
    int wrapper = *((int*)value);
    COM_assertion_msg( wrapper < WRAPPER_MAX && wrapper >= 0
		       , "Illegal value for 'wrapper' option");
    _wrapper = wrapper;
  }
  else if ( option == "verbose"){ 
    _verb = *((int*)value); }
  else if ( option == "lazy"){ 
    _lazy = *((int*)value); }
  else if ( option == "tol"){ 
    COM_assertion_msg( *((float*)value) <= 180. && *((float*)value)>=0. 
                      ,"Illegal value for 'tol' option");
    _tol = *((float*)value); }
  else if ( option == "maxdisp"){ 
    COM_assertion_msg( *((float*)value)>0. 
                      ,"Illegal value for 'maxdisp' option");
    _maxdisp = *((float*)value); }
  else if ( option == "niter"){ 
    _niter = *((int*)value); }
  else if ( option == "ctol"){ 
    COM_assertion_msg( *((float*)value) <= 1. && *((float*)value)>=0. 
                      ,"Illegal value for '_ctol' option");
    _ctol = *((float*)value); }
  else if ( option == "ncycle"){
    _ncycle = *((int*)value); }
  else if ( option == "inverted"){
    _invert_tets = *((int*)value);
  }
  else if ( option == "invert_tets"){
    _invert_tets = *((int*)value); 
  }
  else if ( option == "invert_hexes"){
    _invert_hexes = *((int*)value);
  }
  else COM_assertion_msg( false, "Unknown option");

  //print_legible(1,"Exiting Rocmop::set_value");
}

void Rocmop::smooth_boeing(COM::Attribute* att, int* niter){
  /*
  std::vector<Pane*> allpanes;
  att->window()->panes(allpanes);

  COM_assertion_msg(allpanes.size() == 1,
		    "This function only supports winows with a single pane\n");

  MesqPane mp(allpanes[0],false);

  if(_wrapper == WRAPPER_SHAPE){
    
    Mesquite::MsqError err;
    ShapeImprovementWrapper mesh_quality_algorithm(err);
    MSQ_ERRRTN(err);
    for(int i=0; i< *niter; ++i){
      mesh_quality_algorithm.run_instructions(&mp,err);    
      MSQ_ERRRTN(err);
    }
  }
  else if(_wrapper == WRAPPER_BOEING){
    Mesquite::MsqError err;
    CGWrapper mesh_quality_algorithm(err);
    MSQ_ERRRTN(err);
    for(int i=0; i< *niter; ++i){
      mesh_quality_algorithm.run_instructions(&mp,err);    
      MSQ_ERRRTN(err);
    }
  }
  */
}

void Rocmop::smooth_mesquite(std::vector<COM::Pane*> &allpanes,
			     int ghost_level){

#ifdef ROCPROF
  Profile_begin("Rocmop::smooth_mesquite");
#endif
  Mesquite::MsqError err;
  ShapeImprovementWrapper mesh_quality_algorithm(err);

  int total_npanes = (int)allpanes.size();
  bool wg = (ghost_level == 0) ? false : true;
  
  // std::cout << "wg for Rocmop 114 is: " << wg << std::endl;

  for(int j=0; j < total_npanes; ++j){
    MesqPane mp(allpanes[j],wg);
    MeshSet mesh_set1;
    if(_verb > 4) 
      mp.set_verb(_verb - 4);

    mesh_set1.add_mesh(&mp, err);
    MSQ_ERRRTN(err);
#ifdef ROCPROF
    Profile_begin("Rocmop::Mesquite");
#endif
    mesh_quality_algorithm.run_instructions(mesh_set1,err);
#ifdef ROCPROF
    Profile_end("Rocmop::Mesquite");
#endif
    MSQ_ERRRTN(err);
  }
#ifdef ROCPROF
  Profile_end("Rocmop::smooth_mesquite");
#endif
}

void Rocmop::reduce_sum_on_shared_nodes(COM::Attribute *att){
  Pane_communicator pc(att->window(), att->window()->get_communicator());
  pc.init(att);
  pc.begin_update_shared_nodes();
  pc.reduce_on_shared_nodes(MPI_SUM);
  pc.end_update_shared_nodes();
}

void Rocmop::determine_pane_border(){
#ifdef ROCPROF
  Profile_begin("Rocmop::pane_border");
#endif

  print_legible(1,"Entering Rocmop::determine_pane_border");
  
  std::vector<const COM::Pane*> allpanes;
  _buf_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();

  _is_pane_bnd_node.resize(local_npanes);
  _is_pane_bnd_elem.resize(local_npanes);

  for(int i=0; i< local_npanes; ++i){
    int size_of_real_nodes = allpanes[i]->size_of_real_nodes();
    int size_of_real_elems = allpanes[i]->size_of_real_elements();
    _is_pane_bnd_node[i].resize(size_of_real_nodes,0);
    _is_pane_bnd_elem[i].resize(size_of_real_elems,0);

    std::vector<bool> is_isolated; // is a node isolated?
    MAP::Pane_boundary pb (allpanes[i]);
    pb.determine_border_nodes(_is_pane_bnd_node[i], is_isolated);
  }

  mark_elems_from_nodes(_is_pane_bnd_node,_is_pane_bnd_elem);
  
  print_legible(1,"Exiting Rocmop::determine_pane_border");
#ifdef ROCPROF
  Profile_end("Rocmop::pane_border");
#endif
}


void Rocmop::determine_shared_border(){
#ifdef ROCPROF
  Profile_begin("Rocmop::shared_border");
#endif
  
  print_legible(1,"Entering Rocmop::determine_shared_nodes");

  std::vector<const COM::Pane*> allpanes;
  _buf_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();
    
  _is_shared_node.resize(local_npanes);
  
  //First, get the list of shared nodes.
  for(int i=0; i < (int)(local_npanes); ++i){
    // Obtain the pane connectivity of the local pane.
    const COM::Attribute *pconn = allpanes[i]->attribute(COM::COM_PCONN);
    // Use the pconn offset
    const int *vs = (const int*)pconn->pointer()+Pane_connectivity::pconn_offset();
    int vs_size=pconn->size_of_real_items()-Pane_connectivity::pconn_offset();    
    _is_shared_node[i].resize(allpanes[i]->size_of_real_nodes(),0);
    
    // Determine the number of communicating panes for shared nodes.
    int count=0;
    for (int j=0, nj=vs_size; j<nj; j+=vs[j+1]+2) {
      if (_buf_window->owner_rank( vs[j]) >=0) ++count;
    }
    
    int index = 0;
    // Loop through communicating panes for shared nodes.
    for ( int j=0; j<count; ++j, index+=vs[index+1]+2) {
      // We skip the panes that are not in the current window 
      while ( _buf_window->owner_rank(vs[index])<0) {
	index+=vs[index+1]+2;
	COM_assertion_msg( index<=vs_size, "Invalid communication map");
      }	
      // Add shared nodes to the list
      for(int k=0; k<vs[index+1]; ++k){
	_is_shared_node[i][vs[index+2+k]-1] = 1;
      }
    }
  }

  mark_elems_from_nodes(_is_shared_node,_is_shared_elem);

  print_legible(1,"Exiting Rocmop::determine_shared_nodes");
#ifdef ROCPROF
  Profile_end("Rocmop::shared_border");
#endif
}

void Rocmop::determine_physical_border(){
#ifdef ROCPROF
  Profile_begin("Rocmop::phys_border");
#endif
  print_legible(1,"Entering Rocmop::determine_physical_border()");

  const std::string surf_attr("is_surface");
  COM::Attribute* w_is_surface = _buf_window->attribute(surf_attr);
  COM_assertion_msg( w_is_surface, "Unexpected NULL pointer");
  int is_surface_id = w_is_surface->id();

  std::vector<const COM::Pane*> allpanes;
  _buf_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();
    
  _is_phys_bnd_node.resize(local_npanes);

  for(int i=0; i < local_npanes; ++i){
    _is_phys_bnd_node[i].resize(allpanes[i]->size_of_real_nodes());

    // get pane level pointer to physical border property.
    const COM::Attribute *p_is_surface = allpanes[i]->attribute(is_surface_id);
    int *is_surface_ptr = (int*)p_is_surface->pointer();
    
    // loop through real nodes
    for(int j=0, nj =(int) allpanes[i]->size_of_real_nodes();j<nj; ++j){
      if (is_surface_ptr[j])
	_is_phys_bnd_node[i][j] = true;
    }
  }

  mark_elems_from_nodes(_is_phys_bnd_node,_is_phys_bnd_elem);


  print_legible(1,"Exiting Rocmop::determine_physical_border()");
#ifdef ROCPROF
  Profile_end("Rocmop::phys_border");
#endif
}

void Rocmop::mark_elems_from_nodes(std::vector<std::vector<bool> > &marked_nodes,
				   std::vector<std::vector<bool> > &marked_elems){
#ifdef ROCPROF
  Profile_begin("Rocmop::mark_elem_node");
#endif

  print_legible(1,"Entering Rocmop::mark_elems_from_nodes");

  std::vector<const COM::Pane*> allpanes;
  _buf_window->panes(allpanes);
  int local_npanes = (int)allpanes.size();

  marked_elems.clear();
  marked_elems.resize(local_npanes);

  //Loop through panes
  for(int i=0; i < (int)(local_npanes); ++i){

    marked_elems[i].clear();
    marked_elems[i].resize(allpanes[i]->size_of_real_elements(),false);

    // Loop through real elements.
    // Mark for quality check if they contain shared nodes.
    int s_real_elems = allpanes[i]->size_of_real_elements();
    std::vector<int> nodes;
    for(int j=1; j<= s_real_elems; ++j){
      Element_node_enumerator ene(allpanes[i],j);
      ene.get_nodes(nodes);
      for(int k=0, nk=nodes.size(); k<nk; ++k){
	if (marked_nodes[i][nodes[k]-1])
	  marked_elems[i][j-1] = true;
      }
    }
  }

  print_legible(1,"Exiting Rocmop::mark_elems_from_nodes");
#ifdef ROCPROF
  Profile_end("Rocmop::mark_elem_node");
#endif
}

void Rocmop::invert_elements(int conn_type){
#ifdef ROCPROF
  Profile_begin("Rocmop::invert_elements");
#endif
  print_legible(1,"Entering Rocmop::invert_elements");
  std::vector<Pane*> allpanes;
  _buf_window->panes(allpanes);
  for(int i=0, ni = allpanes.size(); i<ni; ++i){
    MesqPane* mp = new MesqPane(allpanes[i]);
    mp->invert();
    //mp->invert(conn_type);
    if(mp)
      delete mp;
    mp = NULL;
  }  
  print_legible(1,"Exiting Rocmop::invert_elements");
#ifdef ROCPROF
  Profile_end("Rocmop::invert_tets");
#endif
}

double Rocmop::check_marked_elem_quality(std::vector<std::vector<bool> > &marked_elems,
					 std::vector<COM::Pane*> &allpanes){
#ifdef ROCPROF
  Profile_begin("Rocmop::check_marked_quality");
#endif
  print_legible(1,"Entering Rocmop::check_marked_elems");

  double worst_angle = 0.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk=allpanes[i]->size_of_real_elements(); k<nk; ++k){
      if(marked_elems[i][k]){
	Element_node_enumerator ene(allpanes[i],k+1);
	Angle_Metric_3 am;	
	am.initialize(ene);
	am.compute(angles);
	if(angles[1]>worst_angle)
	  worst_angle = angles[1];
      }
    }
  }
  return worst_angle;

  print_legible(1,"Exiting Rocmop::check_marked_elems");
#ifdef ROCPROF
  Profile_end("Rocmop::check_marked_quality");
#endif
}

double Rocmop::check_all_elem_quality(std::vector<const COM::Pane*> &allpanes,
				      bool with_ghosts){

#ifdef ROCPROF
  Profile_begin("Rocmop::all_quality_const");
#endif
  print_legible(1,"Entering Rocmop::check_all_elem_quality");

  int rank =0;
  int ierr = MPI_Comm_rank( _buf_window->get_communicator(),
			    &rank); assert( ierr == 0);
  double worst_angle = 0.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    int nk=allpanes[i]->size_of_real_elements();
    if(with_ghosts)
      nk = allpanes[i]->size_of_elements();
    for(int k =0; k<nk; ++k){
      Element_node_enumerator ene(allpanes[i],k+1);
      Angle_Metric_3 am;
      am.initialize(ene);
      am.compute(angles);
      if(angles[1]>worst_angle)
	worst_angle = angles[1];
    }
  }

  agree_double(worst_angle, MPI_MAX);

#ifdef ROCPROF
  Profile_end("Rocmop::all_quality_const");
#endif
  return worst_angle;

  print_legible(1,"Exiting Rocmop::check_all_elem_quality");
}

double Rocmop::check_all_elem_quality(std::vector<COM::Pane*> &allpanes,
				      bool with_ghosts){
#ifdef ROCPROF
  Profile_begin("Rocmop::all_quality");
#endif
 
  print_legible(1,"Entering Rocmop::check_all_elem_quality");

  if(COMMPI_Initialized()){
    int rank =0;
    int ierr = MPI_Comm_rank( _buf_window->get_communicator(),
			      &rank); assert( ierr == 0);
  }

  double worst_angle = 0.0;
  double angles[] = {0.0, 0.0};
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    int nk=allpanes[i]->size_of_real_elements();
    if(with_ghosts){
      nk = allpanes[i]->size_of_elements();
    }
    for(int k =0; k<nk; ++k){
      Element_node_enumerator ene(allpanes[i],k+1);
      Angle_Metric_3 am;
      am.initialize(ene);
      am.compute(angles);

      if(angles[1]>worst_angle)
	worst_angle = angles[1];
    }
  }
#ifdef ROCPROF
  Profile_end("Rocmop::all_quality");
#endif

  return worst_angle;

  print_legible(1,"Exiting Rocmop::check_all_elem_quality");
}

void Rocmop::print_legible(int verb, const char *msg){
  if(_verb > verb){

  int rank =0;

  if(COMMPI_Initialized()){
    int ierr = MPI_Comm_rank( _buf_window->get_communicator(),
			      &rank); 
    assert( ierr == 0);
  }

  if (rank==0)
    std::cout << "Rocmop: " << msg << std::endl;
  }
}

void Rocmop::constrain_displacements(COM::Attribute * w_disp){
#ifdef ROCPROF
  Profile_begin("Rocmop::const_disp");
#endif
  
  if(_maxdisp > 0.0){

    int disp_id = w_disp->id();
    double max_norm = 0.0;

    std::vector<Pane*> allpanes;
    const_cast<COM::Window*>(_usr_window)->panes(allpanes);
    
    for(int i=0,ni = allpanes.size(); i<ni; ++i){
      
      Vector_3<double> *ptr = NULL;
      COM::Attribute* ptr_att = allpanes[i]->attribute(disp_id);
      void* void_ptr = ptr_att->pointer();
      ptr = reinterpret_cast<Vector_3<double>*>(void_ptr);

      for(int j=0,nj = allpanes[i]->size_of_real_nodes(); j<nj; ++j){	
	double norm = ptr[j].norm();
	if(norm > max_norm)
	  max_norm = norm;      
      }
    }

    agree_double(max_norm,MPI_MAX);
    
    if(max_norm > _maxdisp){
      double div = max_norm/_maxdisp;
      Rocblas::div_scalar(const_cast<const COM::Attribute*>(w_disp),&div,w_disp);
    }
  }
#ifdef ROCPROF
  Profile_end("Rocmop::const_disp");
#endif
}

void Rocmop::print_extremal_dihedrals(COM::Window * window){
#ifdef ROCPROF
  Profile_begin("Rocmop::ext_dihedrals");
#endif

  // Find the rank of this processor, making sure to only make
  // MPI calls on parallel runs
  int myrank =0;
  if(COMMPI_Initialized()){ // true if in parallel
    int ierr = MPI_Comm_rank( window->get_communicator(),
			      &myrank); 
    assert( ierr == 0);
  }

  // Calculate extremal angles on local panes, and
  // record their locations by pane id and element id

  double max_angle = 0.0;
  double min_angle = 180.0; 
  int max_pane = -1; // id of pane w/ the largest dihedral angle
  int min_pane = -1; // id of pane w/ the smallest dihedral angle
  int max_elem = -1; // id of element w/ the largest dihedral angle
  int min_elem = -1; // id of element w/ the smallest dihedral angle
  double angles[] = {0.0, 0.0};
    
  std::vector<Pane*> allpanes;
  window->panes(allpanes);
  
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk=allpanes[i]->size_of_real_elements(); k<nk; ++k){
      Element_node_enumerator ene(allpanes[i],k+1);
      Angle_Metric_3 am;
      am.initialize(ene);
      am.compute(angles);
      if(angles[1]>max_angle){
	max_angle = angles[1];
	max_pane = allpanes[i]->id();
	max_elem = k+1;
      }
      if(angles[0]<min_angle){
	min_angle = angles[0];
	min_pane = allpanes[i]->id();
	min_elem = k+1;
      }
    }
  }

  // Find minimum angle across all panes.

  // create MPI_DOUBLE_INT data types for send/recv buffs
  struct {
    double value;
    int rank;
  } send, recv;
  send.value = min_angle;
  send.rank = myrank;

  if(COMMPI_Initialized()){
    MPI_Allreduce(&send,&recv,1,MPI_DOUBLE_INT,MPI_MINLOC,
	       window->get_communicator());
  }
    if(recv.rank == myrank && _verb > 1){
      std::cout << "Rocmop:" << std::endl << "Rocmop: " 
		<< std::setw(10) << min_angle << " on element "
		<< min_elem << " of pane " << min_pane 
		<< "." << std::endl;
    }
    
  // Sync. output.
  if(COMMPI_Initialized()){
    MPI_Barrier(window->get_communicator());
  }

  // Find maximum angle across all panes.

  send.value = max_angle;

  if(COMMPI_Initialized()){
    MPI_Allreduce(&send,&recv,1,MPI_DOUBLE_INT,MPI_MAXLOC,
	       window->get_communicator());
  }
  if(recv.rank == myrank && _verb > 1){
    std::cout << "Rocmop:" << std::endl << "Rocmop: " 
	      << std::setw(10) << max_angle << " on element "
	      << max_elem << " of pane " << max_pane << std::endl;
  }
  // Sync. output.
  if(COMMPI_Initialized()){
    MPI_Barrier(window->get_communicator());
  }
#ifdef ROCPROF
  Profile_end("Rocmop::ext_dihedrals");
#endif

}

void Rocmop::obtain_extremal_dihedrals(const COM::Attribute *att, double *min_angle, double *max_angle){

  const COM::Window *win = att->window();

  bool cominit = COMMPI_Initialized();

  int rank =0;
  if(cominit){
    int ierr = MPI_Comm_rank( win->get_communicator(),
			      &rank); assert( ierr == 0);
  }

  std::vector<const Pane*> allpanes;
  win->panes(allpanes);

  max_angle[0] = 0.0;
  min_angle[0] = 180.0;
  double angles[] = {0.0, 0.0};

  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk=allpanes[i]->size_of_real_elements(); k<nk; ++k){
      Element_node_enumerator ene(allpanes[i],k+1);
      Angle_Metric_3 am;
      am.initialize(ene);
      am.compute(angles);

      if(angles[1]>max_angle[0])
	max_angle[0] = angles[1];
      if(angles[0]<min_angle[0])
	min_angle[0] = angles[0];
    }
  }

  if(cominit){
    double temp = max_angle[0];
    MPI_Allreduce(&max_angle[0], &temp, 1, MPI_DOUBLE, MPI_MAX,
		  win->get_communicator());
    max_angle[0] = temp;

    temp = min_angle[0];
    MPI_Allreduce(&min_angle[0], &temp, 1, MPI_DOUBLE, MPI_MIN,
		  win->get_communicator());
    min_angle[0] = temp;
  }
}

void Rocmop::print_quality(std::string &s,const std::string &s2){

  string outstr("angles.txt");


  ofstream file;
  ofstream file2;

  bool cominit = COMMPI_Initialized();

  int rank =0;
  if(cominit){
    int ierr = MPI_Comm_rank( _buf_window->get_communicator(),
			      &rank); assert( ierr == 0);
  }

  std::vector<Pane*> allpanes;
  _buf_window->panes(allpanes);

  double max_angle = 0.0;
  double min_angle = 180.0;
  double angles[] = {0.0, 0.0};

  file2.open(s2.c_str());
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk=allpanes[i]->size_of_real_elements(); k<nk; ++k){
      Element_node_enumerator ene(allpanes[i],k+1);
      Angle_Metric_3 am;
      am.initialize(ene);
      am.compute(angles);
      file2 << angles[0] << " " << angles[1] << std::endl;
      if(angles[1]>max_angle)
	max_angle = angles[1];
      if(angles[0]<min_angle)
	min_angle = angles[0];
    }
  }
  file2.close();

  agree_double(max_angle,MPI_MAX);
  agree_double(min_angle,MPI_MIN);

  if (rank==0){
    file.open(outstr.c_str(), std::ios::app);
    COM_assertion_msg( file, "File failed to open\n");
    file << std::left << std::setw(30) << s << std::setw(0)
	 << "(" << min_angle << " , " << max_angle << ")\n";
    file.close();
  }
}

void Rocmop::print_mquality(std::string &s,
			    std::vector<std::vector<bool> > &to_check){

  int ierr = 0, rank =0;

  if(COMMPI_Initialized()){
    ierr = MPI_Comm_rank( _buf_window->get_communicator()
			  , &rank); 
  }  
  assert( ierr == 0);
  

  std::vector<std::vector<bool> > elem_to_check;
  mark_elems_from_nodes(to_check, elem_to_check);

  std::vector<Pane*> allpanes;
  _buf_window->panes(allpanes);
  
  double max_angle = 0.0;
  double min_angle = 180.0;
  double angles[] = {0.0, 0.0};
  int id = -1;
  for(int i=0,ni = allpanes.size(); i<ni; ++i){
    for(int k =0,nk = elem_to_check[i].size(); k<nk; ++k){
      if(elem_to_check[i][k]){
	Element_node_enumerator ene(allpanes[i],k+1);
	Angle_Metric_3 am;
	am.initialize(ene);
	am.compute(angles);
	if(angles[1]>max_angle)
	  max_angle = angles[1];
	if(angles[0]<min_angle){
	  id = k;
	  min_angle = angles[0];
	}
      }
    }
  }

  double temp = min_angle;

  agree_double(max_angle,MPI_MAX);
  agree_double(min_angle,MPI_MIN);
  
  if (rank==0){
    string outstr("angles.txt");
    ofstream file (outstr.c_str(), std::ios::app);
    COM_assertion_msg( file, "File failed to open\n");

    file << std::left << std::setw(30) << s << std::setw(0)
	 << "(" << min_angle << " , " << max_angle << ")";
 
    if(_verb > 1)
      std::cout << std::left << std::setw(30) << "Rocmop: " << s 
                << std::setw(0) << "(" << min_angle << " , " 
	        << max_angle << ")" << std::endl;

    file.close();
  }

  if(COMMPI_Initialized())
    ierr = MPI_Barrier(_buf_window->get_communicator());
  assert( ierr == 0);

  if(min_angle == temp){
    
    string outstr("angles.txt");
    ofstream file (outstr.c_str(), std::ios::app);
    COM_assertion_msg( file, "File failed to open\n");
    
    file << "  worst = (" << rank << " , " << id << ")\n";
    file.close();
  }

  if(COMMPI_Initialized())
    ierr = MPI_Barrier(_buf_window->get_communicator());
  assert( ierr == 0);
}

void Rocmop::perturb_stationary(){
#ifdef ROCPROF
  Profile_begin("Rocmop::perturb_stat");
#endif
  // first calculate current quality of all nodes.
  std::vector<const Pane*> allpanes;
  std::vector<const Pane*> allpanes_usr;
  _buf_window->panes(allpanes);
  _usr_window->panes(allpanes_usr);

  COM::Attribute *w_fgpn_bnd = 
    _buf_window->attribute("is_fgpn_bnd");
  COM::Attribute *w_qual_b_perturb =
    _buf_window->attribute("qual_b_perturb");
  COM::Attribute *w_qual_a_perturb =
    _buf_window->attribute("qual_a_perturb");

  double angles[] = {0.0, 0.0};
  
  // Get worst quality on local panes
  for(int i=0,ni=allpanes.size(); i<ni; ++i){
  
    COM::Attribute *p_fgpn_bnd = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_fgpn_bnd->id()));
    COM::Attribute *p_qual_b_perturb = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_qual_b_perturb->id()));
 
    int* fgpn_bnd_ptr = reinterpret_cast<int*>(p_fgpn_bnd->pointer());
    double* qual_b_perturb_ptr = 
      reinterpret_cast<double*>(p_qual_b_perturb->pointer());

    std::vector<int> elist;
    // only worry about real nodes, ghost nodes should be fixed      
    for(int j=0, nj= allpanes[i]->size_of_real_nodes(); j<nj; ++j){
      if(fgpn_bnd_ptr[j]){
	_dcs[i]->incident_elements(j+1,elist);
	qual_b_perturb_ptr[j] = 0.0;
	for(uint k=0, nk=elist.size(); k<nk; ++k){
	  Element_node_enumerator ene(allpanes[i],k+1);
	  Angle_Metric_3 am;
	  am.initialize(ene);
	  am.compute(angles);
	  if(angles[1]>qual_b_perturb_ptr[j])
	    qual_b_perturb_ptr[j] = angles[1];	
	}
      }
    }
  }

  // Get the worst quality across all panes
  if(COMMPI_Initialized()){
    MAP::Pane_communicator pc(_buf_window, _buf_window->get_communicator());
    pc.init(w_qual_b_perturb);
    pc.begin_update_shared_nodes();
    pc.reduce_on_shared_nodes(MPI_MAX);
    pc.end_update_shared_nodes();
  }    

  // Now generate random perturbations
  for(int i=0,ni=allpanes.size(); i<ni; ++i){
  
    COM::Attribute *p_fgpn_bnd = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_fgpn_bnd->id()));
    COM::Attribute *p_nc = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(COM::COM_NC));
 
    int* fgpn_bnd_ptr = 
      reinterpret_cast<int*>(p_fgpn_bnd->pointer());
    Vector_3<double>* nc_ptr = 
      reinterpret_cast<Vector_3<double>*>(p_nc->pointer());

    std::vector<int> elist;
    // only worry about real nodes, ghost nodes should be fixed      
    for(int j=0, nj= (int)allpanes[i]->size_of_real_nodes(); j<nj; ++j){

      if(fgpn_bnd_ptr[j]){
	_dcs[i]->incident_elements(j+1,elist);

        if(_verb > 1)	
	std::cout << "Rocmop: Perturbing node " << j+1 << std::endl;
	
	//select a random element
	int rand_el = (std::rand()%elist.size());
	Element_node_enumerator ene(allpanes[i],elist[rand_el]);
	std::vector<int> enodes;
	ene.get_nodes(enodes);
	
	// get length of a randomly selected incident edge of the element
	int rand_ed = (std::rand()%3)+1;
	int nindex = 0;
	for(int nk= (int)enodes.size(); nindex<nk; ++nindex){
	  if(enodes[nindex]==j+1)	 
	    break;
	}      
	
	COM_assertion_msg((nindex>=0 && nindex<= (int)enodes.size()),
			  "Node not found in supposedly adjacent element\n");
	double length = (nc_ptr[j]-nc_ptr[enodes[(nindex+rand_ed)%4]-1]).norm();
	
	Vector_3<double> perturb(length,length,length);
	for(int k=0; k<3; ++k){
	  perturb[0] *= (std::rand()%2 == 0) ?
	    (double)((std::rand()%1000)+1)/1000.0 :
	    -1.0*((double)((std::rand()%1000)+1)/1000.0);	  
	}
	nc_ptr[j] += perturb;
      }
    }
  }

    // Get perturbed worst quality on local panes
  for(int i=0,ni=(int)allpanes.size(); i<ni; ++i){
  
    COM::Attribute *p_fgpn_bnd = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_fgpn_bnd->id()));
    COM::Attribute *p_qual_a_perturb = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_qual_a_perturb->id()));
 
    int* fgpn_bnd_ptr = 
      reinterpret_cast<int*>(p_fgpn_bnd->pointer());
    double* qual_a_perturb_ptr = 
      reinterpret_cast<double*>(p_qual_a_perturb->pointer());

    std::vector<int> elist;
    // only worry about real nodes, ghost nodes should be fixed      
    for(int j=0, nj= (int)allpanes[i]->size_of_real_nodes(); j<nj; ++j){
      if(fgpn_bnd_ptr[j]){
	_dcs[i]->incident_elements(j+1,elist);
	qual_a_perturb_ptr[j] = 0.0;
	for(int k=0, nk=(int)elist.size(); k<nk; ++k){
	  Element_node_enumerator ene(allpanes[i],k+1);
	  Angle_Metric_3 am;
	  am.initialize(ene);
	  am.compute(angles);
	  if(angles[1]>qual_a_perturb_ptr[j])
	    qual_a_perturb_ptr[j] = angles[1];	
	}
      }
    }
  }

  // Get the perturbed worst quality across all panes
  if(COMMPI_Initialized()){
    MAP::Pane_communicator pc(_buf_window, _buf_window->get_communicator());
    pc.init(w_qual_a_perturb);
    pc.begin_update_shared_nodes();
    pc.reduce_on_shared_nodes(MPI_MAX);
    pc.end_update_shared_nodes();
  }    

  // Reset positions of any nodes whose worst adj. element quality
  // has gotten worse.
  for(int i=0,ni=allpanes.size(); i<ni; ++i){
    
    COM::Attribute *p_fgpn_bnd = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_fgpn_bnd->id()));
    COM::Attribute *p_nc = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(COM::COM_NC));
    COM::Attribute *p_nc_usr = 
      const_cast<COM::Attribute*>(allpanes_usr[i]->attribute(COM::COM_NC));
      
    COM::Attribute *p_qual_b_perturb = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_qual_b_perturb->id()));
    COM::Attribute *p_qual_a_perturb = 
      const_cast<COM::Attribute*>(allpanes[i]->attribute(w_qual_a_perturb->id()));
 
    int* fgpn_bnd_ptr = 
      reinterpret_cast<int*>(p_fgpn_bnd->pointer());
    Vector_3<double>* nc_ptr = 
      reinterpret_cast<Vector_3<double>*>(p_nc->pointer());
    Vector_3<double>* nc_ptr_usr = 
      reinterpret_cast<Vector_3<double>*>(p_nc_usr->pointer());
    double* qual_b_perturb_ptr = 
      reinterpret_cast<double*>(p_qual_b_perturb->pointer());
    double* qual_a_perturb_ptr = 
      reinterpret_cast<double*>(p_qual_a_perturb->pointer());

    // only worry about real nodes, ghost nodes should be fixed      
for(int j=0, nj= (int)allpanes[i]->size_of_real_nodes(); j<nj; ++j){
      if(fgpn_bnd_ptr[j] &&
	 (qual_a_perturb_ptr[j] > qual_b_perturb_ptr[j]))
	nc_ptr[j] = nc_ptr_usr[j];
    }
  }  
#ifdef ROCPROF
  Profile_end("Rocmop::perturb_stat");
#endif
}


extern "C" void Rocmop_load_module( const char *mname) 
{ Rocmop::load( mname); }

extern "C" void Rocmop_unload_module( const char *mname) 
{ Rocmop::unload( mname); }

// Fortran bindings
extern "C" void COM_F_FUNC2(rocmop_load_module, ROCMOP_LOAD_MODULE)( const char *mname, long int length) 
{ Rocmop::load( std::string(mname, length)); }

extern "C" void COM_F_FUNC2(rocmop_unload_module, ROCMOP_UNLOAD_MODULE)( const char *mname, long int length) 
{ Rocmop::unload( std::string(mname, length)); }

MOP_END_NAMESPACE






