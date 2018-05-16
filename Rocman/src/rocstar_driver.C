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
// $Id: rocstar_driver.C,v 1.94 2009/08/12 04:17:05 mtcampbe Exp $

/** \file rocstar_driver.C
 *  * implements rocstar_driver() the main driver of the rocstar program
 *       */

/* Author: Gengbin Zheng */

#include <fstream>
#include <limits>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "rocman.h"
#if __CHARMC__
#include "charm++.h"
#endif
#include "builtin_couplings.h"
#include "derived_couplings.h"
#ifdef ROCPROF
#include "Rocprof.H"
#endif

COM_EXTERN_MODULE(SimIN);
COM_EXTERN_MODULE(SimOUT);
#ifdef WITH_PANDA
COM_EXTERN_MODULE(Rocpanda);
#endif
void init_profiling(Control_parameters &param, int comm_rank);

int man_verbose = 1;
void RocstarShutdown(int = 0);

void 
RocstarShutdown(int status)
{
  /* Finalize Roccom. */
  COM_finalize();  
#ifdef ROCPROF
  Rocprof_Finalize(true);
#endif
  /* Close down MPI. */
  MPI_Finalize();
  exit(status);
} 

double 
get_restart_time(const string &restart_file)
{
  int curStep = 0;
  double curTime = 0.0;
#if 0
  std::ifstream Inf;
  Inf.open(restart_file.c_str());
  if(!Inf)
    return(0.0);
  while(Inf)
    Inf >> curStep >> curTime;
  Inf.close();
#else
  FILE *fp = fopen(restart_file.c_str(), "r");
  if (fp == NULL) return 0.0;
  while (!feof(fp)) {
      fscanf(fp, "%d %le", &curStep, &curTime);
  }
  fclose(fp);
#endif
  return(curTime);
}

Coupling *create_coupling( Control_parameters &param, const RocmanControl_parameters &rocman_param) {
  const std::string &name = param.coupling_scheme;

  if ( name == "FluidAlone")
    return new FluidAlone( param.fluid_module, param.communicator, &param, &rocman_param);
  else if ( name == "FluidBurnAlone")
    return new FluidBurnAlone( param.fluid_module, param.burn_module, param.communicator, &param, &rocman_param);
  else if ( name == "SolidAlone") 
    return new SolidAlone( param.solid_module, param.communicator, &param, &rocman_param);
  else if ( name == "SolidBurn") 
    return new SolidBurn( param.solid_module, param.burn_module, param.communicator, &param, &rocman_param);
  else if ( name == "SolidFluidSPC") 
    return new SolidFluidSPC( param.fluid_module, param.solid_module, param.communicator, &param, &rocman_param);
  else if ( name == "FluidSolidISS") 
    return new FluidSolidISS( param.fluid_module, param.solid_module, param.communicator, &param, &rocman_param);
  else if ( name == "SolidFluidBurnSPC") 
    return new SolidFluidBurnSPC( param.fluid_module, param.solid_module, param.burn_module, param.communicator, &param, &rocman_param);
  else if ( name == "SolidFluidBurnEnergySPC") 
    return new SolidFluidBurnEnergySPC( param.fluid_module, param.solid_module, param.burn_module, param.communicator, &param, &rocman_param);
  else if ( name == "Test" || name.empty()){
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT,"OUTTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimIN,"INTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(SurfX,"FACETEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocon,"PROPCONTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocmop,"MOPTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocflu,"FLUTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocflo,"FLOTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocprop,"PROPTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocfrac,"FRACTEST");
  }
  else  {
    std::cerr << "Rocstar: ERROR: Unknown coupling scheme: " << name << std::endl;
    std::cerr << std::endl;
    std::cerr << "Rocstar: Rocstar Supported coupling schemes: " << std::endl;
    std::cerr << "Rocstar:    FluidAlone:              Fluid alone with no burn" << std::endl;
    std::cerr << "Rocstar:    FluidBurnAlone:          Fluid alone with burn" << std::endl;
    std::cerr << "Rocstar:    SolidAlone:              Solid alone with no burn" << std::endl;
    std::cerr << "Rocstar:    SolidFluidBurnSPC:       FullyCoupled simple staggered scheme with P-C" << std::endl;
    std::cerr << "Rocstar:    SolidFluidSPC:           FullyCoupled no burn" << std::endl;
    std::cerr << "Rocstar:    FluidSolidISS:           FullyCoupled no burn with improved staggered scheme" << std::endl;
    std::cerr << "Rocstar:    SolidFluidBurnEnergySPC: FullyCoupled with burn energy" << std::endl;
    COM_assertion_msg(0, "ERROR: Unknown coupling scheme!");
  }
  return NULL;
}

#define COM_DOUBLE_ATTRIBUTE(attrname, varname)  \
  attr = winname+"."+attrname; 		\
  COM_new_dataitem( attr.c_str(),'w',COM_DOUBLE,1,"");		\
  COM_set_size( attr.c_str(),0,1);				\
  COM_set_array( attr.c_str(),0, &varname);

#define COM_INT_ATTRIBUTE(attrname, varname)  \
  attr = winname+"."+attrname; 		\
  COM_new_dataitem( attr.c_str(),'w',COM_INT,1,"");	\
  COM_set_size( attr.c_str(),0,1);	\
  COM_set_array( attr.c_str(),0, &varname);

#define COM_STRING_ATTRIBUTE(attrname, varname)  \
  attr = winname+"."+attrname; 			\
  COM_new_dataitem( attr.c_str(),'w',COM_CHAR,1,"");	\
  COM_set_size( attr.c_str(),0,MAXLEN);			\
  COM_set_array( attr.c_str(),0,varname,MAXLEN);

#define COM_BOOL_ATTRIBUTE(attrname, varname)  \
  attr = winname+"."+attrname; 			\
  COM_new_dataitem( attr.c_str(),'w',COM_BOOL,1,"");	\
  COM_set_size( attr.c_str(),0,1);			\
  COM_set_array( attr.c_str(),0,&varname);


Control_parameters::Control_parameters()
{
  // default
  communicator = MPI_COMM_WORLD;
  zoomFactor = 0;
  cur_step = 0;
  maxwalltime = 3.15E7;
  simue_time = 1.0E12;
  strcpy(timingDataDir, "./");
  AutoRestart = 0;
  InterruptFlag = 0;
  maxNumDumps = 0;
  current_dump = 0;
  tolerTract = tolerMass = tolerVelo = tolerDisp = .001;
  maxNumPredCorrCycles = 1;
  maxNumTimeSteps = std::numeric_limits<int>::max();
  strcpy(output_module,"Rocout");
  strcpy(timingDataDir,"Rocman/Profiles");
  controlVerb = 5;
}

void Control_parameters::read()
{
  // Load I/O modules
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");

  string winname = "RocmanParam"; 
  COM_new_window(winname.c_str());

  string attr;
  COM_STRING_ATTRIBUTE("CouplingScheme", coupling_scheme);

  COM_STRING_ATTRIBUTE("FluidModule", fluid_module);

  COM_STRING_ATTRIBUTE("SolidModule", solid_module);

  COM_STRING_ATTRIBUTE("BurnModule", burn_module);

  COM_STRING_ATTRIBUTE("OutputModule", output_module);

  COM_DOUBLE_ATTRIBUTE("InitialTime", current_time);

  COM_DOUBLE_ATTRIBUTE("MaximumTime", simue_time);

  COM_DOUBLE_ATTRIBUTE("MaxWallTime", maxwalltime);

  COM_INT_ATTRIBUTE("MaxNumPredCorrCycles", maxNumPredCorrCycles);

  COM_INT_ATTRIBUTE("MaxNumTimeSteps", maxNumTimeSteps);

  COM_DOUBLE_ATTRIBUTE("CurrentTimeStep", time_step);

  COM_DOUBLE_ATTRIBUTE("OutputIntervalTime", outputIntervalTime);

  COM_DOUBLE_ATTRIBUTE("ZoomFactor", zoomFactor);

  // Tolerances for convergence check
  COM_DOUBLE_ATTRIBUTE("TolerTract", tolerTract);

  COM_DOUBLE_ATTRIBUTE("TolerVelo", tolerVelo);

  COM_DOUBLE_ATTRIBUTE("TolerMass", tolerMass);

  COM_DOUBLE_ATTRIBUTE("TolerDisp", tolerDisp);

  COM_STRING_ATTRIBUTE("GENXTimingDataDir", timingDataDir);
  COM_STRING_ATTRIBUTE("ProfileDir", timingDataDir);
  COM_BOOL_ATTRIBUTE("AutoRestart",AutoRestart);
  COM_INT_ATTRIBUTE("MaxNumDumps",maxNumDumps);
  // 
  COM_window_init_done(winname.c_str());

  //===== Read in parameter file using SimIN
  int IN_param = COM_get_function_handle( "IN.read_parameter_file");

  COM_call_function(IN_param, "RocstarControl.txt", (winname).c_str());

  // post process
  strcat(timingDataDir, "/");

  if(AutoRestart)
    current_time = get_restart_time("Restart.txt");

  init_time = current_time;       // init time, remain constant
  iOutput = (int)((1.000001*current_time)/outputIntervalTime)+1;

  COM_delete_window(winname.c_str());

    // load output module
  const std::string &outputModule = output_module;
  if (outputModule == "Rocout") {
    MAN_DEBUG(3, ("[%d] Rocstar: load_module SimOUT.\n", COMMPI_Comm_rank(communicator)));
    COM_LOAD_MODULE_STATIC_DYNAMIC( SimOUT, "OUT");
  }
#ifdef WITH_PANDA
  else if (outputModule == "Rocpanda") {
    MAN_DEBUG(3, ("[%d] Rocstar: load_module Rocpanda.\n", COMMPI_Comm_rank(communicator)));
    COM_LOAD_MODULE_STATIC_DYNAMIC( Rocpanda, "OUT");
      // Rocpanda server processes won't return
  }
#endif
  else {
    if(COMMPI_Comm_rank(communicator) == 0)
      printf("Rocstar: Error: Unknown output module %s!\n", output_module);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

    // get the MPI communicator set by Rocpanda
    // this is important, update communicator if it is Rocpanda
  update_communicator();

  myRank = COMMPI_Comm_rank(communicator);

    // initialize profiling file
  init_profiling(*this, myRank);

  if (myRank == 0) startTime = MPI_Wtime();
  MPI_Bcast( &startTime, 1, MPI_DOUBLE, 0, communicator);
}

// Rocpanda may require change to communicator
void Control_parameters::update_communicator()
{
  int get_comm = COM_get_function_handle( "OUT.get_comm");
  if (get_comm > 0) {
    COM_call_function( get_comm, &communicator);
    COM_set_default_communicator( communicator);
  }
}

void Control_parameters::print()
{
  printf("================  Rocstar Control file ================\n");
  printf("Rocstar: CouplingScheme: '%s'\n", coupling_scheme);
  printf("Rocstar: InitialTime: %f\n", current_time);
  printf("Rocstar: MaximumTime: %f\n", simue_time);
  printf("Rocstar: MaxWallTime: %f\n", maxwalltime);
  printf("Rocstar: AutoRestart: %s\n",(AutoRestart ? "Yes" : "No")); 
  printf("Rocstar: MaxNumPredCorrCycles: %d\n", maxNumPredCorrCycles);
  printf("Rocstar: MaxNumTimeSteps: %d\n", maxNumTimeSteps);
  if(maxNumDumps > 0)
    printf("Rocstar: MaxNumOutputDumps: %d\n", maxNumDumps);
  printf("Rocstar: CurrentTimeStep: %e\n", time_step);
  printf("Rocstar: OutputIntervalTime: %e\n", outputIntervalTime);
  printf("Rocstar: ZoomFactor: %f\n", zoomFactor);
  printf("Rocstar: TolerTract: %f\n", tolerTract);
  printf("Rocstar: TolerVelo: %f\n", tolerVelo);
  printf("Rocstar: TolerMass: %f\n", tolerMass);
  printf("Rocstar: TolerDisp: %f\n", tolerDisp);
  printf("Rocstar: ProfileDir: %s\n", timingDataDir);
  printf("Rocstar: ProfileFile: %s\n", timingDataFile.c_str());
  printf("=======================================================\n");
}

// 
void Control_parameters::update_start_time(int step, double t)
{
  cur_step = step;
  current_time = t;
  LastOutputTime = t;
  LastOutputStep = cur_step;
    // set output seq no
  iOutput = (int)((1.000001*current_time)/outputIntervalTime)+1;
}

RocmanControl_parameters::RocmanControl_parameters()
{
  // default
  verbose = 1;
  separate_out = 0;

  order = 0;
  traction_mode = 1;    // NO_SHEER;
  rhoc = 1703.0;
  pressure = 6.8e6;
  burn_rate = 0.01;
  P_ambient = 0.0;
  PROP_fom = 0;
  PROP_rediter = 2;
  PROP_fangle = 35.;
  PROPCON_enabled = 0;
  PROPCON_ndiv    = 100;
  async_in = 0;
  async_out = 0;
  rfc_verb = 1;
  rfc_order = 2;
  rfc_iter = 100;
  rfc_tol=1.e-6;

    // internal
  remeshed = 0;
}

void RocmanControl_parameters::read( MPI_Comm comm, int comm_rank) 
{
  const std::string filename = "Rocman/RocmanControl.txt";

  struct stat statBuf;

  string winname = "RocmanControlParam"; 
  COM_new_window(winname.c_str());

  string attr;
  COM_INT_ATTRIBUTE("Verbose", verbose);
  COM_INT_ATTRIBUTE("Separate_out", separate_out);

  COM_INT_ATTRIBUTE("InterpolationOrder", order);
  COM_INT_ATTRIBUTE("TractionMode", traction_mode);
  COM_DOUBLE_ATTRIBUTE("P_ambient", P_ambient);

  COM_DOUBLE_ATTRIBUTE("Rhoc", rhoc);
  COM_DOUBLE_ATTRIBUTE("Pressure", pressure);
  COM_DOUBLE_ATTRIBUTE("BurnRate", burn_rate);

  COM_INT_ATTRIBUTE("RFC_verb", rfc_verb);
  COM_INT_ATTRIBUTE("RFC_order", rfc_order);
  COM_INT_ATTRIBUTE("RFC_iteration", rfc_iter);
  COM_DOUBLE_ATTRIBUTE("RFC_tolerance", rfc_tol);

  COM_BOOL_ATTRIBUTE("Face-offsetting", PROP_fom);
  COM_BOOL_ATTRIBUTE("PROPCON_enabled", PROPCON_enabled);
  COM_INT_ATTRIBUTE("PROP_rediter", PROP_rediter);
  COM_INT_ATTRIBUTE("PROPCON_ndiv", PROPCON_ndiv);
  COM_DOUBLE_ATTRIBUTE("PROP_fangle", PROP_fangle);

  COM_BOOL_ATTRIBUTE("AsyncInput", async_in);
  COM_BOOL_ATTRIBUTE("AsyncOutput", async_out);

      // 
  COM_window_init_done(winname.c_str());

  if (stat(filename.c_str(), &statBuf) == 0) 
  {
    if ( comm_rank == 0) {			// only rank 0 reads
      //===== Read in parameter file using SimIN
      int IN_param = COM_get_function_handle( "IN.read_parameter_file");

      COM_call_function(IN_param, filename.c_str(), (winname).c_str());
    }

    MPI_Bcast(&verbose, 1, MPI_INT, 0, comm);
    MPI_Bcast(&separate_out, 1, MPI_INT, 0, comm);
    MPI_Bcast(&order, 1, MPI_INT, 0, comm);
    MPI_Bcast(&traction_mode, 1, MPI_INT, 0, comm);
    MPI_Bcast(&P_ambient, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&rhoc, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&pressure, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&burn_rate, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&rfc_verb, 1, MPI_INT, 0, comm);
    MPI_Bcast(&rfc_order, 1, MPI_INT, 0, comm);
    MPI_Bcast(&rfc_iter, 1, MPI_INT, 0, comm);
    MPI_Bcast(&rfc_tol, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&PROP_fom, 1, MPI_CHAR, 0, comm);
    MPI_Bcast(&PROP_rediter, 1, MPI_INT, 0, comm);
    MPI_Bcast(&PROPCON_enabled, 1, MPI_CHAR, 0, comm);
    MPI_Bcast(&PROPCON_ndiv, 1, MPI_INT, 0, comm);
    MPI_Bcast(&PROP_fangle, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&async_in, 1, MPI_CHAR, 0, comm);
    MPI_Bcast(&async_out, 1, MPI_CHAR, 0, comm);
  }

  COM_delete_window(winname.c_str());

  man_verbose = verbose;

   // set options for Rocout
  int OUT_set_option = COM_get_function_handle( "OUT.set_option");
  const char *rankWidth = "4";
  COM_call_function( OUT_set_option, "rankwidth", rankWidth);
  COM_call_function( OUT_set_option, "pnidwidth", "0");
  //const char *ioFormat = "HDF";
  const char *ioFormat = "CGNS";
  COM_call_function( OUT_set_option, "format", ioFormat);
  if ( async_out) 
       COM_call_function( OUT_set_option, "async", "on");
  else
       COM_call_function( OUT_set_option, "async", "off");
  if ( separate_out)
       COM_call_function( OUT_set_option, "rankdir", "on");
  else
       COM_call_function( OUT_set_option, "rankdir", "off");
}

void RocmanControl_parameters::print()
{
  printf("==========  Rocman Parameter file read ==========\n");
  printf("Rocstar: verbosity level is %d\n", verbose);
  printf("Rocstar: The order of interpolation is %d\n", order);
  printf("Rocstar: Traction mode is %d (1 for no sheer, 2 for with sheer)\n", traction_mode);
  printf("Rocstar: ambient pressure is %f\n", P_ambient);
  printf("Rocstar: Solid density (Rhoc) is %f kg/m^3\n", rhoc);
  printf("Rocstar: Pressure is %f Pa\n", pressure);
  printf("Rocstar: Burning rate is %f m/s\n", burn_rate);
  printf("Rocstar: RFC_verb: %d\n", rfc_verb);
  printf("Rocstar: Order of quadrature rule (RFC_order): %d\n", rfc_order);
  printf("Rocstar: Max iterations for iterative solver: %d\n", rfc_iter);
  printf("Rocstar: tolerance for iterative solver (RFC_tolerance): %f\n", rfc_tol);
  if (PROP_fom)
    printf("Rocstar: Using face-offsetting method for surface propagation.\n");
  else
    printf("Rocstar: Using marker-particle method for surface propagation.\n");
  if(PROPCON_enabled)
    printf("Rocstar: Using Rocon propagation constraints, (ndiv = %d).\n",PROPCON_ndiv);
  printf("Rocstar: Number of smoothing iterations in Rocprop: %d\n", PROP_rediter);
  printf("Rocstar: Feature-angle threshold in Rocprop: %f\n", PROP_fangle);
  printf("Rocstar: Async Input: %c\n", async_in?'T':'F');
  printf("Rocstar: Async Output: %c\n", async_out?'T':'F');
  printf("==================================================\n");
}

// MaximumTime, maxNumTimeSteps
bool reached_simulation_time( const Control_parameters &param) {
    // if exceed max num of steps
  if (param.cur_step>=param.maxNumTimeSteps) {
    MAN_DEBUG(2, ("\nRocstar: reached_simulation_time returns TRUE with cur_step: %d\n", param.cur_step));
    return true;
  }

  int verb = param.controlVerb;

    // if exceed max wall time
  double startTimeLoopTime;
  if (param.myRank == 0) startTimeLoopTime = MPI_Wtime();
  MPI_Bcast( &startTimeLoopTime, 1, MPI_DOUBLE, 0, param.communicator);
  const double elaped_time = startTimeLoopTime - param.startTime;
  if (elaped_time >= param.maxwalltime ) {
    if (param.myRank == 0 && verb > 1) 
      printf("Rocstar: Quitting; elapsed wall clock time = %f\n", elaped_time);
    return true;
  }
  else {
    if (param.myRank == 0) {
      MAN_DEBUG(2, ("\nRocstar: reached_simulation_time: elapsed wall time = %f \n", elaped_time));
    }
  }

    // if exceed MaximumTime the simulation time
  if (param.current_time >= param.simue_time) {   // MaximumTime
    if (param.myRank == 0 && verb > 0)
      printf("Rocstar: Quitting; CurrentTime = %e, MaximumTime = %e\n", param.current_time, param.simue_time);
    return true;
  }

  // if exceed the Maximum Number of dumps
  if( (param.maxNumDumps > 0 ) && (param.current_dump >= param.maxNumDumps)){
    if(param.myRank == 0 && verb > 0)
      printf("Rocstar: Quitting; Reached the maximum number of dumps for this run. (%d)\n",param.current_dump);
    return true;
  }

  return false;
}

bool reached_restartdump_time( Control_parameters &param) {
  bool result =  param.current_time + param.time_step*0.4 >= 
           param.iOutput*param.outputIntervalTime;
  if (result == true) {
    param.iOutput++;
    if(param.myRank == 0) 
      MAN_DEBUG(1, ("\nRocstar: reached_restartdump_time with current_time: %e time_step: %e iOutput: %d outputIntervalTime: %e\n", param.current_time, param.time_step, param.iOutput, param.outputIntervalTime));
  }
  return result;
}

bool reached_visdump_time( const Control_parameters &param) {
  // TODO: Implement this
  return false;
}
int check_for_interrupt(Coupling *coup,const Control_parameters &param)
{
  int interrupt_code = 0;
  std::string message("Interrupted by file");
  //  std::string therest;
  if(param.myRank == 0) {
    std::ifstream InterruptFile;
    InterruptFile.open("RocstarInterrupt.txt");
    if(InterruptFile){
      std::cout << "Rocstar: Processing interrupt from file." << std::endl;
      InterruptFile >> interrupt_code;
      std::getline(InterruptFile,message);
      InterruptFile.close();
      unlink("RocstarInterrupt.txt");
      std::ofstream InterruptOut;
      InterruptOut.open("RocstarInterrupt.txt.processed");
      InterruptOut << interrupt_code << " " << message << std::endl;
      InterruptOut.close();
    }
  }
  MPI_Bcast(&interrupt_code,1,MPI_INTEGER,0,param.communicator);
  if(interrupt_code > 0){
    coup->Interrupt(&interrupt_code,message.c_str());
    return(interrupt_code);
  }
  return(0);
}

void rocstar_driver( int verb, int remeshed, bool debug) {
  // Read in config files to initialize param
  Control_parameters  *param_ptr = new Control_parameters;
  Control_parameters  &param = *param_ptr;
  param.read();

  int comm_rank = param.myRank;


  param.controlVerb = verb;
  param.controlDebug = debug;

  // Read in rocman config files to initialize param
  RocmanControl_parameters  *rocman_param_ptr = new RocmanControl_parameters;
  RocmanControl_parameters  &rocman_param = *rocman_param_ptr;
  rocman_param.read( param.communicator, comm_rank);
  rocman_param.remeshed = remeshed;

  if ( comm_rank == 0 ) {
    if(debug)
      printf("Rocstar: Call CouplingInitialize\n");
    if(verb > 1){
     param.print();
     rocman_param.print();
    }
    COM_set_verbose( verb);  // shift by 2 to depress first level printing
    //    COM_set_debug( debug);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  //  if(!comm_rank)
  //    std::cout << "Rocstar: creating coupling." << std::endl;

  // Construct coupling object and schedule the actions
  Coupling *coup = create_coupling( param, rocman_param);
  if(!coup){
    if(comm_rank == 0)
      std::cerr << "Rocstar: Error: No coupling created. Exiting." << std::endl;
    return;
  }
  coup->schedule();

  //  std::cout << "Rocstar(" << comm_rank << "): coupling created." << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  //  if(!comm_rank)
  //    std::cout << "Rocstar: All couplings created." << std::endl;
  // visualization
  //  std::string fname = std::string(param.coupling_scheme) + ".gdl";
  //  coup->print(fname.c_str());

  // Load other service modules
  RocBlas::init();
  MPI_Barrier(MPI_COMM_WORLD);
  //  if(!comm_rank)
  //    std::cout << "Rocstar: Rocblas initd." << std::endl;

   // Have Roccom profile all the calls to COM_call_function.
  COM_set_profiling( 1);
  MPI_Barrier(MPI_COMM_WORLD);
  //  if(!comm_rank)
  //    std::cout << "Rocstar: profiling set." << std::endl;

  int with_pciter = (param.maxNumPredCorrCycles > 1);
  if ( comm_rank == 0 && verb > 1) {
    std::cout << "Rocstar: The maximum number of PC-iterations is " << param.maxNumPredCorrCycles << std::endl;
  }
  
  // in case of restarting, reset cur_step and current_time from Restart.txt
  coup->read_restart_info();
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(!comm_rank && debug)
    std::cout << "Rocstar: Read restart information, initializing" << std::endl;

  // call input, init and run_initactions
  coup->initialize();
  
  //  std::cout << "Rocstar(" << comm_rank << "): Coupling initialized." << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  //  if(!comm_rank)
  //    std::cout << "Rocstar: All couplings initd" << std::endl;
  
  if ( param.current_time == 0.0) {
    MPI_Barrier(param.communicator);
    if(comm_rank == 0 && verb > 0)
      std::cout << "Rocstar: Performing time 0 dump..." << std::endl;
    coup->update_integrals( param.current_time);
    coup->update_distances( param.current_time);
    coup->output_restart_files( param.current_time);
    MPI_Barrier(param.communicator);
    if(comm_rank == 0 && debug)
      std::cout << "Rocstar: done." << std::endl;
  }
  
  MPI_Barrier(param.communicator);
  param.LastOutputTime = param.current_time;
  param.LastOutputStep = param.cur_step;
  
  if (comm_rank == 0 && verb > 0) {
    std::cout << "Rocstar: Starting with step " << param.cur_step+1 << ", at time " << param.current_time << std::endl;
  }
  
  char  header[100];
  sprintf(header, "************** Initialization times at time step %d *****************\n", param.cur_step);
  COM_print_profile( param.timingDataFile.c_str(), header);
  COM_set_profiling( 1);     // Reset profiler;

  if (param.current_time == 0.0) 
    coup->write_restart_info( param.current_time, 1);

    // March through time until either the specified maximum time or the maximum
    // number of time steps have been exceeded (whichever comes first).
  int InterfaceConverged;
  MPI_Barrier(param.communicator);
  //  if(comm_rank == 0)
  //    std::cout << "Rocstar: Preparig to step..." << std::endl;
  while ( !reached_simulation_time(param) ) {
    param.cur_step ++;

    for (int iPredCorr = 1; iPredCorr <= param.maxNumPredCorrCycles; iPredCorr++) 
    {
      if ( comm_rank == 0  && verb > 0) {
        std::cout << "Rocstar:" << std::endl
		  << "Rocstar: ================================================================" << std::endl
		  << "Rocstar: System Time Step : " << param.cur_step << "   PC(" << iPredCorr << ")" << std::endl
		  << "Rocstar: ================================================================" << std::endl
		  << "Rocstar:" << std::endl
		  << "Rocstar: CurrentTime, CurrentTimeStep";
        if(verb > 1)
          std::cout << ", ZoomFactor";
        std::cout << ": " 
		  << param.current_time << " " 
		  << param.time_step << " ";
        if(verb > 1)     
	  std::cout << param.zoomFactor;
	std::cout << std::endl << "Rocstar:" << std::endl;
      }

      MPI_Barrier(param.communicator);
      //      if(comm_rank == 0)
      //	std::cout << "Rocstar: initializing convergence checking" << std::endl;
      if ( with_pciter) 
        coup->init_convergence( iPredCorr);
      MPI_Barrier(param.communicator);
      //      if(comm_rank == 0)
      //	std::cout << "Rocstar: convergence checking initialized" << std::endl;

      // invoke time integration of coupling scheme by passing in 
      // current time and desired time step and obtaining new time
      MPI_Barrier(param.communicator);
      //      if(comm_rank == 0)
      //	std::cout << "Rocstar: Running coupling.." << std::endl;
      double new_current_time = coup->run( param.current_time, param.time_step, 
                                           iPredCorr, param.zoomFactor);

      // Check for interrupt
      check_for_interrupt(coup,param);

      //      std::cout << comm_rank << ": solver done on this rank." << std::endl;
      MPI_Barrier(param.communicator);
      //      if(comm_rank == 0)
      //	std::cout << "Rocstar: coupling done." << std::endl;

      InterfaceConverged = 1;
      if ( with_pciter) 
        InterfaceConverged = coup->check_convergence();

      if ( comm_rank == 0 && debug )
        std::cout << "Rocstar:" << std::endl 
		  << "Rocstar: iPredCorr = " << iPredCorr << " is done" << std::endl;
      if ( InterfaceConverged ) {
          param.current_time = new_current_time;
          if ( comm_rank == 0 && verb > 1) 
              std::cout << "Rocstar: Success: predictor-corrector converged at time " << param.current_time << std::endl;
          break;
      }
      else { 
          if ( comm_rank == 0 && debug ) {
	    std::cout << "Rocstar:" << std::endl
		      << "Rocstar: Interface has -NOT- converged!" << std::endl
		      << "Rocstar: iPredCorr = " << iPredCorr << " is done" << std::endl;
	  }
      }
    }   // end of iPredCorr
    
    // New function to parse flags set by interrupt returns false if no action
    if(coup->ProcessInterrupt()) 
      continue;

    if ( !InterfaceConverged) {
      if (comm_rank == 0)
        std::cerr << "Rocstar: Disaster: predictor-corrector did not converge" << std::endl;
      MPI_Abort( MPI_COMM_WORLD, -1);
    }

#if 0
       // warm restart at step 3
    if (param.cur_step == 320) {
        static int didit = 0;     // do it only once
        if (!didit) {
          coup->restart_at_time(param.LastOutputTime, param.LastOutputStep);
          didit = 1;
          continue;		  // skip the rest !
        }
    }
#endif
    
    // if reached time for restart dump
    if ( reached_restartdump_time( param) ) {
      MPI_Barrier(param.communicator);
      // Compute integrals and write them into files for conservation-check
      coup->update_integrals( param.current_time);
      coup->update_distances( param.current_time);
      
      // invoke restart output (as well as visualization data)
      if(comm_rank == 0 && verb > 0)
	std::cout << "Rocstar: Dumping restart files... " << std::endl;
      coup->output_restart_files( param.current_time);
      coup->write_restart_info( param.current_time, param.cur_step+1);
      param.LastOutputTime = param.current_time;
      param.LastOutputStep = param.cur_step;
      param.current_dump++;
      MPI_Barrier(param.communicator);
      //if(comm_rank == 0)
	//std::cout << "done." << std::endl;
    }
    else if ( reached_visdump_time( param)) {
      if(comm_rank == 0 && verb > 0)
	std::cout << "Rocstar: Dumping coupling viz files." << std::endl;
      // output visualization files
      coup->output_visualization_files( param.current_time);
    }

      // Write out profile data
    if ( param.cur_step <= 100 || param.LastOutputTime == param.current_time) {
      sprintf(header, "************** Solver times up to time step %d since last output *********\n", param.cur_step);
      COM_print_profile( param.timingDataFile.c_str(), header);
      COM_set_profiling( 1);     // Reset profiler;
    }
    //if(comm_rank == 0)
      //std::cout << "Rocstar:" << std::endl;
  }   // end time step

    // Write final solution.
  if (param.current_time != param.LastOutputTime) {
    MPI_Barrier(param.communicator);
    if(comm_rank == 0 && verb > 0)
      std::cout << "Rocstar: Performing final dump..." << std::endl;
    coup->output_restart_files( param.current_time);
    coup->write_restart_info( param.current_time, std::min(param.cur_step, param.maxNumTimeSteps)+1);
    MPI_Barrier(param.communicator);
    if(comm_rank == 0 && debug)
      std::cout << "Rocstar: done." << std::endl;
  }

    // Finalize coupling
  coup->finalize();

#ifdef WITH_PANDA
    // unload Rocpanda to exit gracefully
  if ( strcasecmp(param.output_module, "Rocpanda") == 0) {
    MAN_DEBUG(3, ("[%d] Rocstar: unload_module Rocpanda.\n", param.myRank));
    COM_UNLOAD_MODULE_STATIC_DYNAMIC( Rocpanda, "OUT");
  }
#endif

    // Write out the final profile data
  sprintf(header, "************** Finalization times after time step %d *****************\n", param.cur_step);
  COM_print_profile( param.timingDataFile.c_str(), header);
  
  delete coup;
}

void init_profiling(Control_parameters &param, int comm_rank)
{
  std::string fname = "RocstarProfile";
  int data = comm_rank;
  int u = 10;   // <0, 99>
  int c = 2;
  int numProcs = COMMPI_Comm_size( param.communicator);
  while (numProcs >= u*10) {
    u *= 10;
    c ++;
  }
  for (int i=0; i<c; i++) {
    int d = data / u;
    data %= u;
    u /= 10;
    fname += ('0'+d);
  }
  fname += ".txt";
  param.timingDataFile = param.timingDataDir;
  param.timingDataFile += fname;

  int status = 1;  // new
  FILE *fd;
  if (param.current_time == 0.0) {
    fd = fopen( param.timingDataFile.c_str(), "w");
  }
  else {
    fd = fopen( param.timingDataFile.c_str(), "r");
    if (fd == NULL)
      fd = fopen( param.timingDataFile.c_str(), "w");
    else
      status = 0;    // old
  }
  if (fd == NULL) {
    param.timingDataFile = "./";
    param.timingDataFile += fname;
    if (param.current_time == 0) {
      fd = fopen( param.timingDataFile.c_str(), "w");
    }
    else {
      fd = fopen( param.timingDataFile.c_str(), "r");
      if (fd == NULL)
        fd = fopen( param.timingDataFile.c_str(), "w");
      else
        status = 0;
    }
    if (fd == NULL) {
      std::cerr << "Rocstar: Rank " << comm_rank << " could not open " 
                << param.timingDataFile.c_str() << "!" << std::endl;
      MPI_Abort( MPI_COMM_WORLD, -1);
    }
    else if(param.controlVerb > 1) 
      printf("Rocstar: Rank %d using %s instead for timing data!\n", comm_rank, param.timingDataFile.c_str());
  }

  int NumProcs;
  NumProcs = COMMPI_Comm_size( param.communicator);

  if (param.current_time == 0.0 || status == 1) {
    fprintf(fd, "\n");
    fprintf(fd, "Number of processors = %d\n", NumProcs);
    fprintf(fd, "\n");
    fprintf(fd, "MyId = %d\n", comm_rank);
    fprintf(fd, "\n");
  }
  fclose(fd);
}

// print double data of an dataitem on a pane
void debug_print(const std::string str, int pane, int pe, MPI_Comm comm, const char *memo)
{
  int comm_rank = COMMPI_Comm_rank(comm);
  if (comm_rank == pe) {
    double *vm;
    int strid, cap;
    printf("%s %s: before %p\n", str.c_str(), memo?memo:"", vm);
    COM_get_array(str.c_str(), pane, &vm, &strid, &cap);
    printf("%s %s: after %p \n", str.c_str(), memo?memo:"", vm);
    for (int i=0; i<strid*cap; i++) printf("%.17e ", vm[i]);
    printf("\n");
  }
}

void debug_int_print(const std::string str, int pane, int pe, MPI_Comm comm, const char *memo)
{
  int comm_rank = COMMPI_Comm_rank(comm);
  if (comm_rank == pe) {
    int *vm;
    int strid, cap;
    printf("%s %s: before %p\n", str.c_str(), memo?memo:"", vm);
    COM_get_array(str.c_str(), pane, &vm, &strid, &cap);
    printf("%s %s: after %p \n", str.c_str(), memo?memo:"", vm);
    for (int i=0; i<strid*cap; i++) printf("%d ", vm[i]);
    printf("\n");
  }
}






