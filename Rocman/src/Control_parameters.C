#include "Control_parameters.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <sys/stat.h>

#include "rocman.h"

/******************************************************************************
 *      Helper Macros
 ******************************************************************************/

#define COM_DOUBLE_ATTRIBUTE(attrname, varname)                                \
  do {                                                                         \
    attr = winname + "." + attrname;                                           \
    COM_new_dataitem(attr.c_str(), 'w', COM_DOUBLE, 1, "");                    \
    COM_set_size(attr.c_str(), 0, 1);                                          \
    COM_set_array(attr.c_str(), 0, &varname);                                  \
  } while (0)

#define COM_INT_ATTRIBUTE(attrname, varname)                                   \
  do {                                                                         \
    attr = winname + "." + attrname;                                           \
    COM_new_dataitem(attr.c_str(), 'w', COM_INT, 1, "");                       \
    COM_set_size(attr.c_str(), 0, 1);                                          \
    COM_set_array(attr.c_str(), 0, &varname);                                  \
  } while (0)

#define COM_STRING_ATTRIBUTE(attrname, varname)                                \
  do {                                                                         \
    attr = winname + "." + attrname;                                           \
    COM_new_dataitem(attr.c_str(), 'w', COM_CHAR, 1, "");                      \
    COM_set_size(attr.c_str(), 0, MAXLEN);                                     \
    COM_set_array(attr.c_str(), 0, varname, MAXLEN);                           \
  } while (0)

#define COM_BOOL_ATTRIBUTE(attrname, varname)                                  \
  do {                                                                         \
    attr = winname + "." + attrname;                                           \
    COM_new_dataitem(attr.c_str(), 'w', COM_BOOL, 1, "");                      \
    COM_set_size(attr.c_str(), 0, 1);                                          \
    COM_set_array(attr.c_str(), 0, &varname);                                  \
  } while (0)

/******************************************************************************
 *      Helper Functions
 ******************************************************************************/

double get_restart_time(const std::string &restart_file) {
  int curStep = 0;
  double curTime = 0.0;
  FILE *fp = fopen(restart_file.c_str(), "r");
  if (fp == NULL)
    return 0.0;
  while (!feof(fp))
    fscanf(fp, "%d %le", &curStep, &curTime);
  fclose(fp);
  return curTime;
}

void init_profiling(Control_parameters &param, int comm_rank) {
  std::string fname = "RocstarProfile";
  int data = comm_rank;
  int u = 10; // <0, 99>
  int c = 2;
  int numProcs = COMMPI_Comm_size(param.communicator);
  while (numProcs >= u * 10) {
    u *= 10;
    ++c;
  }
  for (int i = 0; i < c; ++i) {
    int d = data / u;
    data %= u;
    u /= 10;
    fname += '0' + d;
  }
  fname += ".txt";
  param.timingDataFile = param.timingDataDir;
  param.timingDataFile += fname;

  int status = 1; // new
  FILE *fd;
  if (param.current_time == 0.0) {
    fd = fopen(param.timingDataFile.c_str(), "w");
  } else {
    fd = fopen(param.timingDataFile.c_str(), "r");
    if (fd == NULL)
      fd = fopen(param.timingDataFile.c_str(), "w");
    else
      status = 0; // old
  }
  if (fd == NULL) {
    param.timingDataFile = "./";
    param.timingDataFile += fname;
    if (param.current_time == 0) {
      fd = fopen(param.timingDataFile.c_str(), "w");
    } else {
      fd = fopen(param.timingDataFile.c_str(), "r");
      if (fd == NULL)
        fd = fopen(param.timingDataFile.c_str(), "w");
      else
        status = 0;
    }
    if (fd == NULL) {
      std::cerr << "Rocstar: Rank " << comm_rank << " could not open "
                << param.timingDataFile.c_str() << "!" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    } else if (param.controlVerb > 1)
      printf("Rocstar: Rank %d using %s instead for timing data!\n", comm_rank,
             param.timingDataFile.c_str());
  }

  int NumProcs;
  NumProcs = COMMPI_Comm_size(param.communicator);

  if (param.current_time == 0.0 || status == 1) {
    fprintf(fd, "\n");
    fprintf(fd, "Number of processors = %d\n", NumProcs);
    fprintf(fd, "\n");
    fprintf(fd, "MyId = %d\n", comm_rank);
    fprintf(fd, "\n");
  }
  fclose(fd);
}

/******************************************************************************
 *      Control_parameters
 ******************************************************************************/

Control_parameters::Control_parameters() {
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
  current_time = 0.0;
  tolerTrac = tolerMass = tolerVelo = tolerDisp = .001;
  maxNumPredCorrCycles = 1;
  maxNumTimeSteps = std::numeric_limits<int>::max();
  strcpy(output_module, "Rocout");
  strcpy(timingDataDir, "Rocman/Profiles");
  controlVerb = 5;
}

void Control_parameters::read() {
  // Load I/O modules
  COM_LOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");

  std::string winname = "RocmanParam";
  COM_new_window(winname.c_str());

  std::string attr;
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
  COM_DOUBLE_ATTRIBUTE("TolerTract", tolerTrac);

  COM_DOUBLE_ATTRIBUTE("TolerVelo", tolerVelo);

  COM_DOUBLE_ATTRIBUTE("TolerMass", tolerMass);

  COM_DOUBLE_ATTRIBUTE("TolerDisp", tolerDisp);

  COM_STRING_ATTRIBUTE("GENXTimingDataDir", timingDataDir);
  COM_STRING_ATTRIBUTE("ProfileDir", timingDataDir);

  COM_BOOL_ATTRIBUTE("AutoRestart", AutoRestart);

  COM_INT_ATTRIBUTE("MaxNumDumps", maxNumDumps);
  //
  COM_window_init_done(winname.c_str());

  //===== Read in parameter file using SimIN
  int IN_param = COM_get_function_handle("IN.read_parameter_file");

  COM_call_function(IN_param, "RocstarControl.txt", winname.c_str());

  // post process
  strcat(timingDataDir, "/");

  if (AutoRestart)
    current_time = get_restart_time("Restart.txt");

  init_time = current_time; // init time, remain constant
  iOutput = (int)((1.000001 * current_time) / outputIntervalTime) + 1;

  COM_delete_window(winname.c_str());

  // load output module
  const std::string &outputModule = output_module;
  if (outputModule == "Rocout") {
    MAN_DEBUG(3, ("[%d] Rocstar: load_module SimOUT.\n",
                  COMMPI_Comm_rank(communicator)));
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
#ifdef WITH_PANDA
  } else if (outputModule == "Rocpanda") {
    MAN_DEBUG(3, ("[%d] Rocstar: load_module Rocpanda.\n",
                  COMMPI_Comm_rank(communicator)));
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocpanda, "OUT");
    // Rocpanda server processes won't return
#endif
  } else {
    if (COMMPI_Comm_rank(communicator) == 0)
      printf("Rocstar: Error: Unknown output module %s!\n", output_module);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // get the MPI communicator set by Rocpanda
  // this is important, update communicator if it is Rocpanda
  update_communicator();

  myRank = COMMPI_Comm_rank(communicator);

  // initialize profiling file
  init_profiling(*this, myRank);

  if (myRank == 0)
    startTime = MPI_Wtime();
  MPI_Bcast(&startTime, 1, MPI_DOUBLE, 0, communicator);
}

// Rocpanda may require change to communicator
void Control_parameters::update_communicator() {
  int get_comm = COM_get_function_handle("OUT.get_comm");
  if (get_comm > 0) {
    COM_call_function(get_comm, &communicator);
    COM_set_default_communicator(communicator);
  }
}

void Control_parameters::print() {
  printf("================  Rocstar Control file ================\n");
  printf("Rocstar: CouplingScheme: '%s'\n", coupling_scheme);
  printf("Rocstar: InitialTime: %f\n", current_time);
  printf("Rocstar: MaximumTime: %f\n", simue_time);
  printf("Rocstar: MaxWallTime: %f\n", maxwalltime);
  printf("Rocstar: AutoRestart: %s\n", (AutoRestart ? "Yes" : "No"));
  printf("Rocstar: MaxNumPredCorrCycles: %d\n", maxNumPredCorrCycles);
  printf("Rocstar: MaxNumTimeSteps: %d\n", maxNumTimeSteps);
  if (maxNumDumps > 0)
    printf("Rocstar: MaxNumOutputDumps: %d\n", maxNumDumps);
  printf("Rocstar: CurrentTimeStep: %e\n", time_step);
  printf("Rocstar: OutputIntervalTime: %e\n", outputIntervalTime);
  printf("Rocstar: ZoomFactor: %f\n", zoomFactor);
  printf("Rocstar: TolerTract: %f\n", tolerTrac);
  printf("Rocstar: TolerVelo: %f\n", tolerVelo);
  printf("Rocstar: TolerMass: %f\n", tolerMass);
  printf("Rocstar: TolerDisp: %f\n", tolerDisp);
  printf("Rocstar: ProfileDir: %s\n", timingDataDir);
  printf("Rocstar: ProfileFile: %s\n", timingDataFile.c_str());
  printf("=======================================================\n");
}

//
void Control_parameters::update_start_time(int step, double t) {
  cur_step = step;
  current_time = t;
  LastOutputTime = t;
  LastOutputStep = cur_step;
  // set output seq no
  iOutput = (int)((1.000001 * current_time) / outputIntervalTime) + 1;
}

/******************************************************************************
 *      RocmanControl_parameters
 ******************************************************************************/

RocmanControl_parameters::RocmanControl_parameters() {
  // default
  verbose = 1;
  separate_out = 0;

  order = 0;
  traction_mode = 1; // NO_SHEER;
  rhoc = 1703.0;
  pressure = 6.8e6;
  burn_rate = 0.01;
  P_ambient = 0.0;
  PROP_fom = 0;
  PROP_rediter = 2;
  PROP_fangle = 35.;
  PROPCON_enabled = 0;
  PROPCON_ndiv = 100;
  async_in = 0;
  async_out = 0;
  rfc_verb = 1;
  rfc_order = 2;
  rfc_iter = 100;
  rfc_tol = 1.e-6;

  // internal
  remeshed = 0;
}

void RocmanControl_parameters::read(MPI_Comm comm, int comm_rank) {
  const std::string filename = "Rocman/RocmanControl.txt";

  struct stat statBuf {};

  std::string winname = "RocmanControlParam";
  COM_new_window(winname.c_str());

  std::string attr;
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

  if (stat(filename.c_str(), &statBuf) == 0) {
    if (comm_rank == 0) { // only rank 0 reads
      //===== Read in parameter file using SimIN
      int IN_param = COM_get_function_handle("IN.read_parameter_file");

      COM_call_function(IN_param, filename.c_str(), winname.c_str());
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
  int OUT_set_option = COM_get_function_handle("OUT.set_option");
  const char *rankWidth = "4";
  COM_call_function(OUT_set_option, "rankwidth", rankWidth);
  COM_call_function(OUT_set_option, "pnidwidth", "0");
  // Let SimOUT determine the file format
  // COM_call_function(OUT_set_option, "format", ioFormat);
  if (async_out)
    COM_call_function(OUT_set_option, "async", "on");
  else
    COM_call_function(OUT_set_option, "async", "off");
  if (separate_out)
    COM_call_function(OUT_set_option, "rankdir", "on");
  else
    COM_call_function(OUT_set_option, "rankdir", "off");
}

void RocmanControl_parameters::print() {
  printf("==========  Rocman Parameter file read ==========\n");
  printf("Rocstar: verbosity level is %d\n", verbose);
  printf("Rocstar: The order of interpolation is %d\n", order);
  printf("Rocstar: Traction mode is %d (1 for no sheer, 2 for with sheer)\n",
         traction_mode);
  printf("Rocstar: ambient pressure is %f\n", P_ambient);
  printf("Rocstar: Solid density (Rhoc) is %f kg/m^3\n", rhoc);
  printf("Rocstar: Pressure is %f Pa\n", pressure);
  printf("Rocstar: Burning rate is %f m/s\n", burn_rate);
  printf("Rocstar: RFC_verb: %d\n", rfc_verb);
  printf("Rocstar: Order of quadrature rule (RFC_order): %d\n", rfc_order);
  printf("Rocstar: Max iterations for iterative solver: %d\n", rfc_iter);
  printf("Rocstar: tolerance for iterative solver (RFC_tolerance): %f\n",
         rfc_tol);
  if (PROP_fom)
    printf("Rocstar: Using face-offsetting method for surface propagation.\n");
  else
    printf("Rocstar: Using marker-particle method for surface propagation.\n");
  if (PROPCON_enabled)
    printf("Rocstar: Using Rocon propagation constraints, (ndiv = %d).\n",
           PROPCON_ndiv);
  printf("Rocstar: Number of smoothing iterations in Rocprop: %d\n",
         PROP_rediter);
  printf("Rocstar: Feature-angle threshold in Rocprop: %f\n", PROP_fangle);
  printf("Rocstar: Async Input: %c\n", async_in ? 'T' : 'F');
  printf("Rocstar: Async Output: %c\n", async_out ? 'T' : 'F');
  printf("==================================================\n");
}
