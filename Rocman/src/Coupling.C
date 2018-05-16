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
// $Id: Coupling.C,v 1.53 2008/12/06 08:45:22 mtcampbe Exp $

#include "rocman.h"

COM_EXTERN_MODULE(Simpal);
COM_EXTERN_MODULE(SurfX);
COM_EXTERN_MODULE(SurfUtil);
COM_EXTERN_MODULE(Rocprop);

class SolidAgent;
class FluidAgent;

void RocstarShutdown(int = 0);
#ifdef ARM
int Rocrem_remesh(const char* solver, const char* indir, const char* outdir,
		  double time, bool pt,bool remesh_surf, 
		  bool transfer_surf, double scaleFactor, MPI_Comm myComm,
		  int fieldWidth, int fileBase, int debug, int ngLayers);
#endif

int RocBlas::copy_scalar = 0;
int RocBlas::copy = 0;
int RocBlas::add = 0;
int RocBlas::sub = 0;
int RocBlas::div_scalar = 0;
int RocBlas::limit1 = 0;
int RocBlas::sub_scalar = 0;
int RocBlas::axpy_scalar = 0;
int RocBlas::mul = 0;
int RocBlas::div = 0;
int RocBlas::neg = 0;
int RocBlas::axpy = 0;
int RocBlas::nrm2 = 0;
int RocBlas::mul_scalar = 0;
int RocBlas::max_scalar_MPI = 0;
int RocBlas::min_scalar_MPI = 0;
int RocBlas::sum_scalar_MPI = 0;
int RocBlas::nrm2_scalar_MPI = 0;
int RocBlas::maxof_scalar = 0;


void RocBlas::initHandles()
{
  copy_scalar = COM_get_function_handle( "BLAS.copy_scalar");
  sub_scalar = COM_get_function_handle( "BLAS.sub_scalar");
  copy = COM_get_function_handle( "BLAS.copy");
  add = COM_get_function_handle( "BLAS.add");
  sub = COM_get_function_handle( "BLAS.sub");
  axpy = COM_get_function_handle( "BLAS.axpy");
  limit1 = COM_get_function_handle( "BLAS.limit1");
  mul = COM_get_function_handle( "BLAS.mul");
  div = COM_get_function_handle( "BLAS.div");
  neg = COM_get_function_handle( "BLAS.neg");
  nrm2 = COM_get_function_handle( "BLAS.nrm2");
  div_scalar = COM_get_function_handle( "BLAS.div_scalar");
  axpy_scalar = COM_get_function_handle( "BLAS.axpy_scalar");
  mul_scalar = COM_get_function_handle( "BLAS.mul_scalar");
  max_scalar_MPI = COM_get_function_handle( "BLAS.max_scalar_MPI");
  min_scalar_MPI = COM_get_function_handle( "BLAS.min_scalar_MPI");
  sum_scalar_MPI = COM_get_function_handle( "BLAS.sum_scalar_MPI");
  nrm2_scalar_MPI = COM_get_function_handle( "BLAS.nrm2_scalar_MPI");
  maxof_scalar = COM_get_function_handle( "BLAS.maxof_scalar");
}

void RocBlas::init()
{
  COM_LOAD_MODULE_STATIC_DYNAMIC( Simpal, "BLAS");
  initHandles();
}

void Coupling::baseInit()
{
  iPredCorr = 0;
  init_started = param->current_time==0.0;
  restarting = 0;
  init_remeshed = rocmanparam->remeshed;
  comm_rank = COMMPI_Comm_rank(param->communicator);
  restartInfo = "Restart.txt";
  integFname = "Rocman/Modout/GENX_integ.txt";
  distFname = "Rocman/Modout/GENX_dist.txt";
  overwrite_integ = overwrite_dist = 0;
  if ( param->current_time != 0.0) {
    FILE *fp = fopen(integFname.c_str(), "r");
    if (!fp) overwrite_integ = 1;
    fp = fopen(distFname.c_str(), "r");
    if (!fp) overwrite_dist = 1;
  }
  else {
    overwrite_integ = 1;
    overwrite_dist = 1;
  }
  // experimental Interrupt
  std::string manwinname = "Rocman";
  std::string manfuncname = manwinname+".interrupt";
  std::string coupobjname = manwinname+".coup";
  COM_Type types[3];
  types[0] = COM_RAWDATA;
  types[1] = COM_INT;
  types[2] = COM_STRING;
  COM_new_window(manwinname);
  COM_new_dataitem(coupobjname, 'w',COM_VOID, 1, "");
  COM_set_object(coupobjname.c_str(),0,this);
  COM_set_member_function(manfuncname.c_str(),
			  (Member_func_ptr)&Coupling::Interrupt,
			  coupobjname.c_str(),"bii",types);
}

Coupling::Coupling(const char *coupl_name, const char *name, 
                   Control_parameters *param_, 
                   const RocmanControl_parameters *mp): 
          coupling_name(coupl_name), param(param_), rocmanparam(mp)
{
  modules.push_back(name);
  baseInit();
}

Coupling::Coupling(const char *coupl_name, const char *fluidname, 
                   const char *solidname, 
                   Control_parameters *param_,
                   const RocmanControl_parameters *mp):
          coupling_name(coupl_name), param(param_), rocmanparam(mp)
{
  modules.push_back(fluidname);
  modules.push_back(solidname);
  baseInit();
}

Coupling::Coupling(const char *coupl_name, const char *fluidname, 
                   const char *solidname, const char *burnname,
                   Control_parameters *param_,
                   const RocmanControl_parameters *mp):
          coupling_name(coupl_name), param(param_), rocmanparam(mp)
{
  modules.push_back(fluidname);
  modules.push_back(solidname);
  modules.push_back(burnname);
  baseInit();
}

// Delete all agents
Coupling::~Coupling() {
  for ( int i=0, n=agents.size(); i<n; ++i) {
    delete agents[i];
  }
  agents.clear();
}

int Coupling::new_start(double t) const { 
  return t==0.0 || (rocmanparam->remeshed && t == param->init_time); 
}

Agent *Coupling::add_agent( Agent *agent)
{
  agents.push_back(agent);
  return agent;
}

// Schedule the actions for the scheduler and the agents
void Coupling::schedule() {
  init_scheduler.schedule();
  scheduler.schedule();

  for ( int i=0, n=agents.size(); i<n; ++i) {
    agents[i]->schedule();
  }
}

// Invoke initialization of the actions in the scheduler and the agents
void Coupling::init( double t, double dt, int reinit) {
  int i, n;

#if 0
  for ( i=0, n=agents.size(); i<n; ++i) 
    agents[i]->init_module( t);

  for ( i=0, n=agents.size(); i<n; ++i) 
    agents[i]->init_buffers( t);
#else
//  void _load_rocface(FluidAgent *fagent, SolidAgent *sagent, const RocmanControl_parameters *param);
  for ( i=0, n=agents.size(); i<n; ++i)  {
//    if (i==2) _load_rocface(agents[0], agents[1], rocmanparam);
    agents[i]->init_module( t, dt);
//    agents[i]->init_buffers( t);
  }
#endif

    // in warm restart, reset scheduler so that it can be re-inited
  if (reinit)
    callMethod(&Scheduler::restarting, t);

  callMethod(&Scheduler::init_actions, t);
}

// invoke method fn to all schedulers, and all schedulers of agents
void Coupling::callMethod(Scheduler_voidfn1_t fn, double t)
{
  int i, n;

  (init_scheduler.*fn)( t);

  for ( i=0, n=agents.size(); i<n; ++i) 
    agents[i]->callMethod(fn, t);

  (scheduler.*fn)( t);
}

// Invoke finalize of the actions in the scheduler and the agents
void Coupling::finalize() {
  init_scheduler.finalize_actions();
  scheduler.finalize_actions();

  for ( int i=0, n=agents.size(); i<n; ++i) 
    agents[i]->finalize();
}

double Coupling::run( double t, double dt, int pred, double zoom) 
{ 
  iPredCorr = pred;

  // Change dt to the maximum time setp allowed.
  for ( int i=0, n=agents.size(); i<n; ++i) {
    dt = std::min( dt, agents[i]->max_timestep(t, dt));
  }
  
  MAN_DEBUG(3, ("[%d] Rocstar: Coupling::run with t: %f dt: %f.\n", comm_rank, t, dt));

  scheduler.run_actions( t, dt); 

  // Return the new time as t + dt
  if (zoom > 0.0)
    return t + dt*zoom;
  else
    return t + dt;
}

void Coupling::run_initactions( double t, double dt)
{
  MAN_DEBUG(3, ("[%d] Rocstar: %s::run_initactions called with t:%e dt:%e.\n", comm_rank, name(), t, dt));
  init_scheduler.set_alpha(0.0);
  init_scheduler.run_actions( t, dt); 

  init_started = 0;
  init_remeshed = 0;
}

// Invoke input functions of the agents
void Coupling::input( double t) {
  for ( int i=0, n=agents.size(); i<n; ++i) {
    agents[i]->input( t);
  }
}

void Coupling::init_convergence( int iPredCorr)
{
  if ( maxPredCorr>1 && iPredCorr > 0) {    
    for ( int i=0, n=agents.size(); i<n; ++i) {
      agents[i]->init_convergence( iPredCorr);
    }
  }
}

int Coupling::check_convergence()
{
  int InterfaceConverged;
  if ( maxPredCorr>1) {
    InterfaceConverged = 0;
    for ( int i=0, n=agents.size(); i<n; ++i)
      if (!agents[i]->check_convergence( param->tolerMass, param->tolerTract, param->tolerVelo)) return InterfaceConverged;
    InterfaceConverged = 1;
  }
  else
    InterfaceConverged = 1;
  return InterfaceConverged;
}

// Write out restart files (including visualization data)
void Coupling::output_restart_files(double t) {
  for ( int i=0, n=agents.size(); i<n; ++i) {
    agents[i]->output_restart_files( t);
  }
}

// Write out visualization files
void Coupling::output_visualization_files(double t) {
  for ( int i=0, n=agents.size(); i<n; ++i) {
    agents[i]->output_visualization_files( t);
  }
}

// find restart time from Restart.txt
void Coupling::read_restart_info()
{
  if (param->current_time != 0.0) {
    FILE *fp = fopen(restartInfo.c_str(), "r");
    if (fp == NULL) {
      std::cerr << "Rocstar: Error: READ_RESTART: Failed to read file " << restartInfo << std::endl;
      MPI_Abort( MPI_COMM_WORLD, -1);
    }
    int curStep;
    double initialTime;
    while (!feof(fp)) {
      fscanf(fp, "%d %le", &curStep, &initialTime);
    }
    fclose(fp);
    if (comm_rank == 0)
      MAN_DEBUG(1,("Rocstar: This run is a restart continued at iteration %d with initial time %f.",curStep, initialTime));
    param->update_start_time(curStep - 1, initialTime);   // subtle - it starts from 0
  }
  else {
    if (comm_rank == 0)
      MAN_DEBUG(2,("Rocstar: This run is not a restart.\n"));
  }
}

void Coupling::write_restart_info(double CurrentTime, int iStep)
{
  FILE *fp;
  if ( comm_rank == 0 ) {
       if ( CurrentTime == 0.0)
         fp = fopen(restartInfo.c_str(), "w");
       else
         fp = fopen(restartInfo.c_str(), "a");
       if (fp == NULL) {
         std::cerr << "Rocstar: Error: Failed to open restart info file, " 
		   << restartInfo << "." << std::endl;
         MPI_Abort( MPI_COMM_WORLD, -1);
       }
       fprintf( fp, "%d %.20le \n", iStep, CurrentTime);
       fclose( fp);
  }
}

void Coupling::initialize(int reinit)
{
  // invoke input of coupling object (Input of initial data of physics modules)
  input( param->current_time);

  // invoke initialization of coupling object, load modules
  init( param->current_time, param->time_step, reinit);

  // reset rocface if restart
  if (reinit)
     reload_rocface(rocmanparam);   

  // initialize each modules by executing init schedulers
  run_initactions( param->current_time, param->time_step);
}

// experimental code
void Coupling::restart_at_time(double t, int step)
{
  int i, n;

  MAN_DEBUG(3, ("Rocstar: Coupling::restart_at_time with t: %f step:%d.\n", t, step));

#if 0
  // HACK
  system("rm RocburnAPN/Rocout/*");
  if (comm_rank == 0) {
    system("rm RocburnAPN/Rocout/*");
  }
  MPI_Barrier(get_control_param()->communicator);
#endif

  restarting = 1;

    // reset windows and buffers
  finalize();

    // unload modules, only call unload to zero out data, skip dlclose
  for ( i=0, n=agents.size(); i<n; ++i) {
    agents[i]->unload_module();
  }

  if (COM_get_window_handle("RFC") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");
  }
  if (COM_get_window_handle("SURF") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC( SurfUtil, "SURF");
  }
  if (COM_get_window_handle("PROP") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC( Rocprop, "PROP");
  }
  if (COM_get_window_handle("PROPCON") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocon,"PROPCON");
  }
  if (COM_get_window_handle("BLAS") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(Simpal,"BLAS");
  }
  if (COM_get_window_handle("MAP") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SurfMap,"MAP");
  }
    // update to given time
  param->update_start_time( step, t);

    // load module and call physics routine
  for ( i=0, n=agents.size(); i<n; ++i) {
    agents[i]->load_module();
  }

    // proceed to initialization phase
  initialize(1);    // reinit

/*
    // read restart data
  for ( i=0, n=agents.size(); i<n; ++i) {
    agents[i]->read_restart_data();
  }
*/
  restarting = 0;
}

// print in GDL
void Coupling::print(const char *fname)
{
  FILE *f = fopen(fname, "w");

  fprintf(f, "graph: { title: \"Coupling\" \n\
        display_edge_labels: yes \n\
        layoutalgorithm: tree   \n\
        scaling: maxspect   \n\
        node.color     : green   \n\
        node.textcolor : black   \n\
        node.bordercolor: black \n\
        node.borderwidth: 1    \n\
        edge.color     : blue   \n\
        edge.arrowsize : 7   \n\
        edge.thickness : 2   \n\
        edge.fontname:\"helvO08\"  \n\
        node.label: \"no type\" \n");

//  scheduler.print(f, name());
  for (unsigned int i=0; i<agents.size(); i++) {
    agents[i]->print(f);
    if (i>0) fprintf(f, "edge: { sourcename: \"%s\" targetname: \"%s\" label: \"\"}\n", agents[i-1]->get_agent_name().c_str(), agents[i]->get_agent_name().c_str());
  }

  fprintf(f, "} \n");
  fclose(f);
}

std::string Coupling::normalize_modname(const char* modname)
{
  std::string mod = modname;
  if (mod == "RocfluDummy")
    return "Rocflu";
  if (mod == "RocfloDummy")
    return "Rocflo";
  if (mod == "RocfracDummy")
    return "Rocfrac";
  if (mod == "RocsolidDummy")
    return "Rocsolid";
  return mod;
}


// experimental code
void Coupling::Interrupt(int *act,const char *message)
{
  int action = *act;
  if(!param->myRank){
    std::cout << "Rocstar::Interrupt invoked." << std::endl
	      << "Rocstar: ***************************************************" 
	      << std::endl
	      << "Rocstar: " << (message ? message : "") << std::endl
	      << "Rocstar: ***************************************************" 
	      << std::endl;
  }
  switch(action){
  case 0:
    // Just stop
    if(!param->myRank)
      std::cout << "Rocstar: Halting simulation." << std::endl;
    finalize();
    RocstarShutdown();
    break;
  case 1:
    // Output solutions and stop.
    if(!param->myRank)
      std::cout << "Rocstar: Writing restart files and halting simulation." 
		<< std::endl;
    output_restart_files( param->current_time);
    finalize();
    RocstarShutdown();
    break;
  case 2:
    //
    // Restart from last restart dump and increase the output
    // frequency for debugging.
    //
    // 0. Call warm restart function to junk current data, and 
    //    re-initialize with last successful solution dump
    // 1. Reset the max number of timesteps and output interval
    // 2. Code will stop and dump a soln when it reaches the step just prior 
    //    to this one
    //
    if(!param->myRank){
      std::cout << "Rocstar: Directed to warm restart from last dump "
		<< "and increase dump frequency." << std::endl;
    }
    param->InterruptFlag = 2;
    break;
  case 3:
    if(!param->myRank)
      std::cout << "Rocstar: Directed to remesh and warm restart from "
		<< "(time/step): ("
		<< param->LastOutputTime << "," << param->LastOutputStep 
		<< ")" << std::endl;
    param->InterruptFlag = 3;
    break;
  case 4:
    if(!param->myRank)
      std::cout << "Rocstar: Directed to restart from "
               << "(time/step): (" << param->LastOutputTime << ","
               << param->LastOutputStep << ")" << std::endl;
    if(message)
      std::cout << "MESSAGE: " << message << std::endl;
    param->InterruptFlag = 4;
    break;
  case 5:
    if(!param->myRank)
      std::cout << "Rocstar: Directed to dump at (time/step): ("
               << param->current_time << "," << param->cur_step
               << ")" << std::endl;
    output_restart_files( param->current_time);
    write_restart_info( param->current_time, param->cur_step+1);
    break;
  case 6:
    if(!param->myRank)
      std::cout << "Rocstar: Directed to dump and restart at (time/step): ("
               << param->current_time << "," << param->cur_step
               << ")" << std::endl;
    output_restart_files( param->current_time);
    write_restart_info( param->current_time, param->cur_step+1);
    if(message)
      std::cout << "MESSAGE: " << message << std::endl;
    param->InterruptFlag = 4;
    break;
  default:
    if(param->myRank == 0)
      std::cout << "Rocstar: Unknown interrupt action." << std::endl;
    MAN_DEBUG(0,("Rocstar: Coupling::Interrupt(unknown action): %s.\n",
		 message));
    return;
  }
}


int
Coupling::ProcessInterrupt()
{
  int save_step = 0;
  bool done = false;
  int fluid_agent_index = 0;
  int action = param->InterruptFlag;
  param->InterruptFlag = 0;
  if(action > 0){
    if(param->myRank == 0)
      std::cout << "Rocstar: Processing interrupt." << std::endl;
    switch(action){
    case 2:
      save_step = param->cur_step;
      restart_at_time(param->LastOutputTime,param->LastOutputStep);
      param->maxNumTimeSteps = save_step;
      param->outputIntervalTime /= 10.0;
      if(!param->myRank){
	std::cout << "Rocstar: Restarting from (time/step): ("
		  << param->current_time << "," << param->cur_step 
		  << ") with output"
		  << " interval, " << param->outputIntervalTime << "." 
		  << std::endl;
      }
      return(1);
      break;
    case 3:
      for(int i = 0,n=agents.size();i<n && !done;n++){
	std::string::size_type x = agents[i]->get_agent_name().find("lu");
	if( x != string::npos){
	  fluid_agent_index = i;
	  done = true;
	}
      }
      if(!done){
	std::cerr << "Rocstar: Could not find fluid agent.  Dying." 
		  << std::endl;
	RocstarShutdown(1);
      }
      // char* solver: Path to Rocout directory (i.e. blah/blah/Rocflu)
      // double time: timestep to read
      // bool use_parallel_transfer
      // bool remesh_surf: true=remesh surface false=leave surface alone
      // bool transfer_surf: true (only valid if remesh_surf = false)
      // double scaleFactor: scale mesh granularity
      // int fieldWidth: digits in filename partition id
      // int fileBase: 0 is partition ids start at 0, 1 otherwise
      // int debug: 0=off, 1 = some, more = more
      // int ngLayers: number of ghost layers; 0-9
      // returns 1 if success; 0 if fail
      // Rocrem_remesh(const char* solver, double time, bool remesh_surf, 
      //	       bool transfer_surf, double scaleFactor, MPI_Comm myComm,
      //	       int fieldWidth, int fileBase, int debug, int ngLayers) 
      //
      // Remeshing test case
      //
#ifdef ARM
      if(!Rocrem_remesh("Rocflu","Rocout","",param->LastOutputTime,true,true,false,1.0,
			agents[fluid_agent_index]->get_communicator(),
			4,0,2,2)){
	if(!param->myRank)
	  std::cerr << "Rocstar: Remeshing failed.  Stopping simulation." 
		    << std::endl;
	RocstarShutdown(1);
      }
#else
      if(!param->myRank){
	//	std::ofstream Ouf;
	//	Ouf.open("needs_remesh");
	//	Ouf.close();
	std::cerr << "Rocstar: Online remeshing not enabled. Shutting down." << std::endl;
      }
      RocstarShutdown(1);
#endif
      MPI_Barrier(MPI_COMM_WORLD);
      restart_at_time(param->LastOutputTime,param->LastOutputStep);
      MPI_Barrier(MPI_COMM_WORLD);
      if(!param->myRank)
	std::cout << "Rocstar: Restarting after remesh at (time/step): ("
		  << param->current_time << "," << param->cur_step << ")"
		  << std::endl;
      return(1);
      break;
    case 4:
      MPI_Barrier(MPI_COMM_WORLD);
      restart_at_time(param->LastOutputTime,param->LastOutputStep);
      MPI_Barrier(MPI_COMM_WORLD);
      if(!param->myRank)
       std::cout << "Rocstar: Warm restarting at (time/step): ("
                 << param->current_time << "," << param->cur_step << ")"
                 << std::endl;
      return(1);
      break;
    default:
      if(!param->myRank)
	std::cout << "Rocstar: Unknown interrupt action: " << action << "." 
		  << std::endl;
      return(0);
      break;
    }
  }
  return(0);
}

// Compute the integrals and dump into a file.
void FullyCoupling::update_integrals(double currentTime)
{
  int i, n, count=0;
  for ( i=0, n=agents.size(); i<n; ++i) {
    count += agents[i]->compute_integrals();
  }

  if ( comm_rank == 0 && count == 2) {
       FILE *fp;
       if ( overwrite_integ)
         fp = fopen(integFname.c_str(), "w");
       else
         fp = fopen(integFname.c_str(), "a");
       if (fp == NULL) {
         std::cerr << "Rocstar: Failed to open integral file, " 
		   << integFname << "." << std::endl;
         MPI_Abort( MPI_COMM_WORLD, -1);
       }
       MAN_DEBUG(3, ("Rocstar: FullyCoupling::update_integrals with t: %f.\n", currentTime));
       if ( overwrite_integ) {
         fprintf( fp, " time    f-volume        s-volume        f-mass          s-mass          f-burn area     s-burn area     f-non-burn area s-non-burn area s-volume-undef\n");
         overwrite_integ = 0;
       }
       double *f_integrals = fluid_agent->get_integrals();
       double *s_integrals = solid_agent->get_integrals();
       fprintf( fp, "%.18le %.18le %.18le %.18le %.18le %.18le %.18le %.18le %.18le %.18le\n", 
             currentTime,
             f_integrals[MAN_INTEG_VOL], s_integrals[MAN_INTEG_VOL],
             f_integrals[MAN_INTEG_MASS], s_integrals[MAN_INTEG_MASS],
             f_integrals[MAN_INTEG_IBAREA], s_integrals[MAN_INTEG_IBAREA],
             f_integrals[MAN_INTEG_INBAREA], s_integrals[MAN_INTEG_INBAREA],
             f_integrals[MAN_INTEG_VOL_UND]);
               
       fclose( fp);
  }
}

// Compute the distances between fluid and solid meshes and dump into a file.
void FullyCoupling::update_distances(double currentTime)
{
  double dist_max=-1, dist_min=-1, dist_nrm2=-1;

  MPI_Comm communicator = fluid_agent->get_communicator();

  int RFC_interpolate = COM_get_function_handle("RFC.interpolate");
  COM_call_function( RFC_interpolate, &solid_agent->y_hdl, &fluid_agent->nc_tmp_hdl);
  COM_call_function( RocBlas::sub, &fluid_agent->nc_tmp_hdl, &fluid_agent->nc_hdl, &fluid_agent->nc_tmp_hdl);
  COM_call_function( RocBlas::nrm2, &fluid_agent->nc_tmp_hdl, &fluid_agent->sq_dist_hdl);
  COM_call_function( RocBlas::max_scalar_MPI, &fluid_agent->sq_dist_hdl, &dist_max, &communicator);
  COM_call_function( RocBlas::min_scalar_MPI, &fluid_agent->sq_dist_hdl, &dist_min, &communicator);
  COM_call_function( RocBlas::sum_scalar_MPI, &fluid_agent->sq_dist_hdl, &dist_nrm2, &communicator);

  if ( comm_rank == 0) {
       FILE *fp;
       if ( overwrite_dist)
         fp = fopen(distFname.c_str(), "w");
       else
         fp = fopen(distFname.c_str(), "a");
       if (fp == NULL) {
         std::cerr << "Rocstar: Failed to open distance file, " << distFname 
		   << "." << std::endl;
         MPI_Abort( MPI_COMM_WORLD, -1);
       }
       MAN_DEBUG(3, ("Rocstar: FullyCoupling::update_distances with t: %f.\n", currentTime));
       if ( overwrite_dist) {
         fprintf( fp, "# time    distance-min    distance-max    distance-norm2  \n");
         overwrite_dist = 0;
       }
       fprintf( fp, "%.18le %.18le %.18le %.18le\n", 
             currentTime, sqrt(dist_min), sqrt(dist_max), sqrt(dist_nrm2));
               
       fclose( fp);
  }
}

void _load_rocface(FluidAgent *fagent, SolidAgent *sagent, const RocmanControl_parameters *param);

void FullyCoupling::reload_rocface(const RocmanControl_parameters *param)
{
  _load_rocface(fluid_agent, solid_agent, param);
}








