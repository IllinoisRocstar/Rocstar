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
/** \file Agent.h
 *  * Contains declaration of the base class for Agent implementations.
 *   * @see Agent.C
 *    */
/* Author: Gengbin Zheng */

#ifndef _AGENT_H_
#define _AGENT_H_

using namespace std;
#include <vector>
#include <cstring>

#include "com.h"
#include "commpi.h"

class Scheduler;
class Action;
class Agent;
class Coupling;
class InterpolateBase;

class AttributeBase {
public:
  Agent *agent;
  std::string  target_window;
  std::string  attr;
public:
  AttributeBase( Agent *ag, std::string target_window_, std::string attr_): 
		agent(ag), target_window(target_window_), attr(attr_) {}
  virtual ~AttributeBase() {}
  virtual void create(std::string bufname);
  virtual void assign(std::string bufname);
};

class NewAttribute : public AttributeBase {
protected:
  // new_dataitem
  char loc;
  int type;
  int ncomp;
  const char *unit;
public:
  NewAttribute(Agent *ag, std::string target_window_, std::string attr_, char loc_, int type_, int ncomp_, const char * unit_);
  void create(std::string bufname);
};

// clone_dataitem
class CloneAttribute : public AttributeBase {
protected:
  std::string  parent_window;
  std::string  parent_attr;  
  int          wg; 
  const char * ptnname;
  int          val;
  int          condition;
public:
  CloneAttribute(Agent *ag, int cond, std::string target_window_, std::string attr_, std::string parent_window_, std::string parent_attr_, int wg_=1, const char *ptnname_=0, int val_=0);
  void create(std::string bufname);
};

// use_dataitem
class UseAttribute : public AttributeBase {
protected:
  std::string  parent_window;
  std::string  parent_attr;  
  int          wg; 
  const char * ptnname;
  int          val;
public:
  UseAttribute(Agent *ag, std::string target_window_, std::string attr_, std::string parent_window_, std::string parent_attr_, int wg_=1, const char *ptnname_=0, int val_=0);
  void create(std::string bufname);
};

typedef vector< UserScheduler* > SchdulerList;

class PhysicsAction: public Action {
private:
  Agent *agent;
public:
  PhysicsAction(Agent *ag): Action(0, (const char **)NULL, NULL, NULL, (char *)"PhysicsAction"), agent(ag) {}
  void run(double t, double dt, double alpha);
  void declare(Scheduler &) {}
  virtual void print(FILE *f);
  virtual char *name();
};

class Agent: public COM_Object {
friend class PhysicsAction;
protected:
  const std::string     agent_name;
  MPI_Comm        communicator;
  Coupling	  *coupling;            // 
  PhysicsAction   action;
  SchdulerList    bcScheduler;
  UserScheduler   bcInitScheduler;
  UserScheduler   gmScheduler;
  UserScheduler   icScheduler;

  // Rocin and Rocout function handles
  static int read_by_control_handle;
  static int read_files_handle;
  static int obtain_attr_handle;
  static int write_attr_handle;
  static int write_ctrl_handle;

  // Handles to routines provided by physics modules
  int comm_rank;
  int init_handle, update_handle;
  int pre_hdf_handle, post_hdf_handle;
  int finalize_handle, compute_integrals_handle;  // fluid code handles

  // Hanldes to routines provided by agent to physics modules
  int  ic_handle, bc_handle, gm_handle;
  bool withgm;

  int dobackup;     // do backup on interpolation

  string rocmod_name, mod_instance;
  string inDir, outDir, suffix;

  double initial_time, timestamp, current_time, current_deltatime, old_dt;

  // surf_window can be solidWin of INITIALIZE_SOLID() in SolidAgent ,
  //  or fluidSurf of INITIALIZE_FLUID() in FluidAgent
  std::string surf_window, vol_window;
  std::string tmp_window;
  void *option_data;

  typedef vector< InterpolateBase* > InterpolateList;
  InterpolateList interpolateList;

  typedef vector< AttributeBase* > AttributeList;
  AttributeList dataitemList;
  
  int pc_hdls[2][3];   // Handles for dataitems to be 
                       // stored/restored for PC-iterations
  int pc_count;	       // burn only has one pc_hdls

  double integrals[MAN_INTEG_SIZE];       // Array for storing integrals

  bool skipInputIO;
public:
  Agent(Coupling *cp, std::string mod, std::string obj, 
	const char *agent_name, MPI_Comm com, bool wgm=false, bool skipio=false);
  virtual ~Agent();
  virtual void load_module() {}
  virtual void unload_module() {}
  virtual void init_module(double t, double dt);
  void init_subscheduler(double t);
  void callMethod(Scheduler_voidfn1_t fn, double t);
  void schedule();
  virtual void finalize();

  virtual void input(double t) = 0;
  virtual void output_restart_files(double t) = 0;
  virtual void output_visualization_files(double t) = 0;

  double max_timestep(double t, double dt);

  Action *get_main_action() { return &action; }

  // for creating dataitems by Actions
  void register_new_dataitem(std::string target_window_, std::string attr_, char loc_, int type_, int ncomp_, const char* unit_);
  void register_clone_dataitem(int cond, std::string target_window_, std::string attr_, std::string parent_window_, std::string parent_attr_, int wg_=1, const char *ptnname_=0, int val_=0);
  void register_use_dataitem(std::string target_window_, std::string attr_, std::string parent_window_, std::string parent_attr_, int wg_=1, const char *ptnname_=0, int val_=0);
  void create_registered_dataitems(std::string tmpBuf);
  void create_registered_window_dataitems(std::string target_window);
  virtual void create_buffer_all();
  void assign_dataitems();
  virtual void read_restart_data()  {}

  void add_data() { /* TODO */ }
  void add_icaction(Action *act);
  void add_bcaction(Action *act, int l=1);
  void add_bcinitaction(Action *act);
  void add_gmaction(Action *act);

  void init_callback( const char *surf_win, const char *vol_win, 
		      void *option=NULL);

  void init_bcactions(double t);
  void obtain_bc(double *a, int *l=NULL);

  void init_gmactions(double t);
  void obtain_gm(double *da);

  void init_bcinitaction(double t);
  void run_bcinitaction(double t, double dt);

  Coupling *get_coupling() { return coupling; }

  void register_interpolate(InterpolateBase *ip); 

  std::string get_surface_window() const { return surf_window; }
  std::string get_volume_window() const { return vol_window; }

  MPI_Comm  get_communicator() const { return communicator; }
  std::string get_rocmod_name() const { return rocmod_name; }
  std::string get_modinstance_name() const { return mod_instance; }
  std::string get_agent_name() const { return agent_name; }
  int  get_comm_rank() const { return comm_rank; }
  
  double get_old_dt() const { return old_dt; }

  virtual void init_convergence( int iPredCorr);
  virtual int check_convergence( double tolerMass, double tolerTract, double tolerVelo) { return 1; }

  virtual int compute_integrals() { return 0; }
  double * get_integrals() { return integrals; }

  void print(FILE *f);

  void get_time_string( double t, std::string &s);

protected:
  Scheduler * get_bcScheduler(unsigned int level) { return level-1<bcScheduler.size()?bcScheduler[level-1]:NULL; }

  void init_function_handles();
  void create_window(const char *window_name);

  // Split surface window into separate ones for output.
  void split_surface_window( const std::string surfAll,
			   const std::string surf_i, 
			   const std::string surf_nb, 
			   const std::string surf_b, 
			   const std::string surf_ni);

  // Read in HDF/CGNS files by control file. Return 0 if succeeds.
  // Return -1 if control file was not found.
  int  read_by_control_file( double t, const std::string base, const std::string window);

  // Write out the given dataitem into HDF/CGNS files with the given base
  // file name. If ref exist, then refer to it for mesh data.
  void write_data_files( double t, const std::string base, 
			 const std::string attr, const char *ref=NULL);

  // Write out Rocin file with the given base file name and window name.
  void write_control_file( double t, const std::string base, const std::string window);

  virtual void store_solutions( int converged);
  int check_convergence_help( int vcur, int vpre, double tol, std::string str);
};


#endif






