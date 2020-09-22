#ifndef _ROCSTAR_CONTROL_PARAMETERS_HPP_
#define _ROCSTAR_CONTROL_PARAMETERS_HPP_

#include "com.h"

const int MAXLEN = 128;

class Control_parameters {
public:
  char coupling_scheme[MAXLEN]; ///< Name of the coupling scheme
  char fluid_module[MAXLEN];    ///< Names of fluid module
  char solid_module[MAXLEN];    ///< Names of solid module
  char burn_module[MAXLEN];     ///< Names of burn module
  char output_module[MAXLEN];   ///< Names of IO module

  MPI_Comm communicator; ///< MPI communicator
  double simue_time;     ///< Maximum simulation time
  double init_time;      ///< Initial time stamp
  double current_time;   ///< Current time stamp, init_time at beginning
  double time_step;      ///< Current time step
  double maxwalltime;    ///< Maximum wall time
  int maxNumDumps;       ///< Maximum number of output dumps
  int current_dump;      ///< Current dump number

  int maxNumPredCorrCycles; ///< Maximum Number of predictor-corrector cycles
  int maxNumTimeSteps;      ///< Maximum number of time steps

  double outputIntervalTime; ///< Output interval time
  double zoomFactor;         ///< Time zooming factor

  double tolerTrac; ///< Traction tolerance for convergence check
  double tolerVelo; ///< Velocity tolerance for convergence check
  double tolerMass; ///< Mass tolerance for convergence check
  double tolerDisp; ///< Displacement tolerance for convergence check

  char timingDataDir[MAXLEN];
  std::string timingDataFile;

  int cur_step; ///< Current time step

  double LastOutputTime;
  int LastOutputStep;
  int InterruptFlag;
  int AutoRestart;

  int controlVerb;
  bool controlDebug;

  // internal
  int iOutput;
  int myRank;
  double startTime;

public:
  Control_parameters();
  void read();
  void print();
  void update_communicator();
  void update_start_time(int step, double t);
};

#endif //_ROCSTAR_CONTROL_PARAMETERS_HPP_
