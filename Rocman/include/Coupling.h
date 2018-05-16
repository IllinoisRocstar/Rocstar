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
// $Id: Coupling.h,v 1.38 2009/11/04 15:15:14 mtcampbe Exp $

#ifndef _COUPLING_H_
#define _COUPLING_H_

using namespace std;
#include <vector>
#include <string>

#include "FluidAgent.h"
#include "SolidAgent.h"

class Scheduler;
class BurnAgent;
class RocmanControl_parameters;

typedef vector< Agent* > AgentList;

#define MAXLEN 128

class Control_parameters {
public:
  char coupling_scheme[MAXLEN];  //< Name of the coupling scheme
  char fluid_module[MAXLEN];     //< Names of fluid module
  char solid_module[MAXLEN];     //< Names of solid module
  char burn_module[MAXLEN];      //< Names of burn module
  char output_module[MAXLEN];     //< Names of IO module

  MPI_Comm    communicator;
  double  simue_time;           //< Maximum simulation time
  double  init_time;            //< initial time stamp
  double  current_time;         //< Current time stamp, init_time at beginning
  double  time_step;            //< current time step
  double  maxwalltime;
  int  maxNumDumps;
  int  current_dump;

  int  maxNumPredCorrCycles;
  int  maxNumTimeSteps;

  double  outputIntervalTime;
  double  zoomFactor;

  double  tolerTract;		// < Tolerances for convergence check
  double  tolerVelo;
  double  tolerMass;
  double  tolerDisp;

  char  timingDataDir[MAXLEN];
  std::string timingDataFile;

  int     cur_step;

  double  LastOutputTime;
  int     LastOutputStep;
  int     InterruptFlag;
  int     AutoRestart;

  int     controlVerb;
  bool    controlDebug;

  // internal
  int     iOutput;
  int     myRank;
  double  startTime;
public:
  Control_parameters();
  void read();
  void print();
  void update_communicator();
  void update_start_time(int step, double t);
};

// coupling
class Coupling : public COM_Object
{
protected:
  std::string coupling_name;
  vector<std::string>  modules;
  AgentList agents;
  UserScheduler scheduler;       //< NOTE: special simple user defined scheduler
  UserScheduler init_scheduler;   //<  for initialization
  int comm_rank;
  int init_started;                    // initial_start
  int restarting;		       // in restarting
  int init_remeshed;                   // initialization after remeshing
  // PC iteration
  int iPredCorr;		  // class
  int maxPredCorr;		    // modified

  Control_parameters *param;
  const RocmanControl_parameters *rocmanparam;

  std::string restartInfo;

    // compute integrals
  int overwrite_integ;
  std::string integFname;

  int overwrite_dist;
  std::string distFname;
private:
  void baseInit();
public:
  /// Constructor. Derived class will add actions for the coupling scheme
  Coupling(const char *coupl_name, const char *name, Control_parameters *p, const RocmanControl_parameters *mp);
  Coupling(const char *coupl_name, const char *fluidname, const char *solidname, Control_parameters *p, const RocmanControl_parameters *mp);
  Coupling(const char *coupl_name, const char *fluidname, const char *solidname, const char *burnname, Control_parameters *p, const RocmanControl_parameters *mp);

  /// Destructor
  virtual ~Coupling();

  const char *name() { return coupling_name.c_str(); }

  /// Add new agent.
  Agent *add_agent( Agent *);

  /// Schedule the top-level actions of the coupling scheme and 
  /// the actions of the agents.
  void schedule();

  /// Invoke initialization of the actions in the scheduler and the agents
  void init(double t, double dt, int reinit=0);

  void initialize(int reinit=0);

  /// Invoke finalization of the actions in the scheduler and the agents
  void finalize();

  /// Invoke the scheduler
  double run(double t, double dt, int iPredCorr, double zoom);
  void run_initactions( double t, double dt);

  /// Invoke input functions of the agents
  void input(double t);

  ///
  int get_ipc() const { return iPredCorr; }
  int get_max_ipc() const { return maxPredCorr; }

  int initial_start() const { return init_started; }
  int in_restart() const { return restarting; }

  /// true if in initialization step and remeshed is true
  int initial_remeshed() const { return init_remeshed; }

  int new_start(double t) const;

  void init_convergence( int iPredCorr);
  int check_convergence();

  virtual void update_integrals(double currentTime) {}
  virtual void update_distances(double currentTime) {}

  // Write out restart files (including visualization data)
  void output_restart_files(double t);

  // Write out visualization files
  void output_visualization_files(double t);

  const Control_parameters *get_control_param() { return param; }
  const RocmanControl_parameters *get_rocmancontrol_param() { return rocmanparam; }

  void read_restart_info();
  void write_restart_info(double CurrentTime, int iStep);

  void restart_at_time(double t, int step);
  virtual void reload_rocface(const RocmanControl_parameters *param) {}

  // print for visualization
  void print(const char *fname);

  // Experimental Interrupt
  void Interrupt(int *,const char *);
  int  ProcessInterrupt();

protected:
  void callMethod(Scheduler_voidfn1_t fn, double t);
  std::string normalize_modname(const char* mod);
};

// subclass fully coupled involving fluid and solid, and maybe burn
class FullyCoupling: public Coupling {
protected:
  FluidAgent *fluid_agent;
  SolidAgent *solid_agent;
  BurnAgent  *burn_agent;
public:
  FullyCoupling(const char *coupl_name, const char *fluidname, const char *solidname, Control_parameters *p, const RocmanControl_parameters *mp):
      Coupling(coupl_name, fluidname, solidname, p, mp), burn_agent(NULL) {}
  FullyCoupling(const char *coupl_name, const char *fluidname, const char *solidname, const char *burnname, Control_parameters *p, const RocmanControl_parameters *mp):
      Coupling(coupl_name, fluidname, solidname, burnname, p, mp) {}

  virtual void update_integrals(double currentTime);
  virtual void update_distances(double currentTime);
  virtual void reload_rocface(const RocmanControl_parameters *param);
};

// fully coupled
#define DECLARE_NEW_FULLY_COUPLING_SCHEME( New_scheme)			\
  class New_scheme : public FullyCoupling {					\
  public:								\
    New_scheme(const char *, const char *, MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp);		\
    New_scheme(const char *, const char *, const char*, MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp);	\
  };

#endif






