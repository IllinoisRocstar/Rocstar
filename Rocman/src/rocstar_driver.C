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
#include <string>
#include <unistd.h>

#include "Control_parameters.h"
#include "RocBlas.h"
#include "rocman.h"

#include "builtin_couplings.h"
#include "derived_couplings.h"

#if __CHARMC__
#include "charm++.h"
#endif
#ifdef ROCPROF
#include "Rocprof.H"
#endif

COM_EXTERN_MODULE(SimIN)
COM_EXTERN_MODULE(SimOUT)
#ifdef WITH_PANDA
COM_EXTERN_MODULE(Rocpanda)
#endif

void init_profiling(Control_parameters &param, int comm_rank);

int man_verbose = 1;
void RocstarShutdown(int = 0);

void RocstarShutdown(int status) {
  /* Finalize Roccom. */
  COM_finalize();
#ifdef ROCPROF
  Rocprof_Finalize(true);
#endif
  /* Close down MPI. */
  MPI_Finalize();
  exit(status);
}

RocstarCoupling *create_coupling(Control_parameters &param,
                                 const RocmanControl_parameters &rocman_param) {
  const std::string &name = param.coupling_scheme;

  if (name == "FluidAlone")
    return new FluidAlone(param.communicator, &param, &rocman_param,
                          param.fluid_module);
  else if (name == "FluidBurnAlone")
    return new FluidBurnAlone(param.communicator, &param, &rocman_param,
                              param.fluid_module, param.burn_module);
  else if (name == "SolidAlone")
    return new SolidAlone(param.communicator, &param, &rocman_param,
                          param.solid_module);
  else if (name == "SolidBurnAlone" || name == "SolidBurn")
    return new SolidBurnAlone(param.communicator, &param, &rocman_param,
                              param.solid_module, param.burn_module);
  else if (name == "SolidFluidSPC")
    return new SolidFluidSPC(param.communicator, &param, &rocman_param,
                             param.fluid_module, param.solid_module);
  else if (name == "SolidFluidBurnSPC")
    return new SolidFluidBurnSPC(param.communicator, &param, &rocman_param,
                                 param.fluid_module, param.solid_module,
                                 param.burn_module);
  else if (name == "SolidFluidISS" || name == "FluidSolidISS")
    return new SolidFluidISS(param.communicator, &param, &rocman_param,
                             param.fluid_module, param.solid_module);
  else if (name == "SolidFluidBurnEnergySPC")
    return new SolidFluidBurnEnergySPC(param.communicator, &param,
                                       &rocman_param, param.fluid_module,
                                       param.solid_module, param.burn_module);
  else if (name == "Test" || name.empty()) {
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUTTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimIN, "INTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(SurfX, "FACETEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocon, "PROPCONTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocmop, "MOPTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocflu, "FLUTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocflo, "FLOTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocprop, "PROPTEST");
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocfrac, "FRACTEST");
  } else {
    std::cerr << "Rocstar: ERROR: Unknown coupling scheme: " << name << "\n"
              << R"(
Rocstar: Rocstar supported coupling schemes:
Rocstar: Built-in coupling schemes:
Rocstar:    FluidAlone:              Fluid alone with no burn
Rocstar:    FluidBurnAlone:          Fluid alone with burn
Rocstar:    SolidAlone:              Solid alone with no burn
Rocstar:    SolidBurnAlone:          Solid alone with burn
Rocstar:    SolidFluidSPC:           FullyCoupled no burn with
Rocstar:    SolidFluidBurnSPC:       FullyCoupled with burn with
                                     simple staggered scheme with P-C
Rocstar: Derived coupling schemes:
Rocstar:    SolidFluidISS:           FullyCoupled no burn with
                                     improved staggered scheme
Rocstar:    SolidFluidBurnEnergySPC: FullyCoupled with burn energy)"
              << std::endl;
    COM_abort_msg(EXIT_FAILURE, "ERROR: Unknown coupling scheme!");
  }
  return nullptr;
}

// MaximumTime, maxNumTimeSteps
bool reached_simulation_time(const Control_parameters &param) {
  int verb = param.controlVerb;

  // if exceed max num of steps
  if (param.cur_step >= param.maxNumTimeSteps) {
    MAN_DEBUG(
        2,
        ("\nRocstar: reached_simulation_time returns TRUE with cur_step: %d\n",
         param.cur_step));
    return true;
  }

  // if exceed max wall time
  double startTimeLoopTime;
  if (param.myRank == 0)
    startTimeLoopTime = MPI_Wtime();
  MPI_Bcast(&startTimeLoopTime, 1, MPI_DOUBLE, 0, param.communicator);
  const double elaped_time = startTimeLoopTime - param.startTime;
  if (elaped_time >= param.maxwalltime) {
    if (param.myRank == 0 && verb > 1)
      printf("Rocstar: Quitting; elapsed wall clock time = %f\n", elaped_time);
    return true;
  } else {
    if (param.myRank == 0) {
      MAN_DEBUG(
          2, ("\nRocstar: reached_simulation_time: elapsed wall time = %f \n",
              elaped_time));
    }
  }

  // if exceed MaximumTime the simulation time
  if (param.current_time >= param.simue_time) { // MaximumTime
    if (param.myRank == 0 && verb > 0)
      printf("Rocstar: Quitting; CurrentTime = %e, MaximumTime = %e\n",
             param.current_time, param.simue_time);
    return true;
  }

  // if exceed the Maximum Number of dumps
  if ((param.maxNumDumps > 0) && (param.current_dump >= param.maxNumDumps)) {
    if (param.myRank == 0 && verb > 0)
      printf("Rocstar: Quitting; Reached the maximum number of dumps for this "
             "run. (%d)\n",
             param.current_dump);
    return true;
  }

  return false;
}

bool reached_restartdump_time(Control_parameters &param) {
  bool result = param.current_time + param.time_step * 0.4 >=
                param.iOutput * param.outputIntervalTime;
  if (result) {
    param.iOutput++;
    if (param.myRank == 0)
      MAN_DEBUG(1, ("\nRocstar: reached_restartdump_time with current_time: %e "
                    "time_step: %e iOutput: %d outputIntervalTime: %e\n",
                    param.current_time, param.time_step, param.iOutput,
                    param.outputIntervalTime));
  }
  return result;
}

bool reached_visdump_time(const Control_parameters &param) {
  // TODO: Implement this
  return false;
}

/* AEG: disabling unused experimental code.
int check_for_interrupt(RocstarCoupling *coup,
                        const Control_parameters &param) {
  int interrupt_code = 0;
  std::string message("Interrupted by file");
  //  std::string therest;
  if (param.myRank == 0) {
    std::ifstream InterruptFile;
    InterruptFile.open("RocstarInterrupt.txt");
    if (InterruptFile) {
      std::cout << "Rocstar: Processing interrupt from file." << std::endl;
      InterruptFile >> interrupt_code;
      std::getline(InterruptFile, message);
      InterruptFile.close();
      unlink("RocstarInterrupt.txt");
      std::ofstream InterruptOut;
      InterruptOut.open("RocstarInterrupt.txt.processed");
      InterruptOut << interrupt_code << " " << message << std::endl;
      InterruptOut.close();
    }
  }
  MPI_Bcast(&interrupt_code, 1, MPI_INTEGER, 0, param.communicator);
  if (interrupt_code > 0) {
    coup->Interrupt(&interrupt_code, message.c_str());
    return interrupt_code;
  }
  return 0;
}
*/

void rocstar_driver(int verb, int remeshed, bool debug) {
  // Read in config files to initialize param
  auto *param_ptr = new Control_parameters;
  Control_parameters &param = *param_ptr;
  param.read();

  param.controlVerb = verb;
  param.controlDebug = debug;

  // Read in rocman config files to initialize param
  auto *rocman_param_ptr = new RocmanControl_parameters;
  RocmanControl_parameters &rocman_param = *rocman_param_ptr;
  rocman_param.read(param.communicator, param.myRank);
  rocman_param.remeshed = remeshed;

  if (param.myRank == 0) {
    if (debug)
      std::cout << "Rocstar: Call CouplingInitialize" << std::endl;
    if (verb > 1) {
      param.print();
      rocman_param.print();
    }
    COM_set_verbose(verb); // shift by 2 to depress first level printing
    // COM_set_debug(debug);
  }

  /*
  MPI_Barrier(MPI_COMM_WORLD);
  if (param.myRank == 0)
    std::cout << "Rocstar: creating coupling." << std::endl;
  */

  // Construct coupling object and schedule the actions
  RocstarCoupling *coup = create_coupling(param, rocman_param);
  if (!coup) {
    if (param.myRank == 0)
      std::cerr << "Rocstar: Error: No coupling created. Exiting." << std::endl;
    return;
  }

  /*
  std::cout << "Rocstar(" << param.myRank << "): coupling created."
            << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  if (param.myRank == 0)
    std::cout << "Rocstar: All couplings created." << std::endl;
  // visualization
  std::string fname = std::string(param.coupling_scheme) + ".gdl";
  coup->print(fname.c_str());
  */

  // Load other service modules
  RocBlas::init();
  /*
  MPI_Barrier(MPI_COMM_WORLD);
  if (param.myRank == 0)
    std::cout << "Rocstar: Rocblas initialized." << std::endl;
  */

  // Have Roccom profile all the calls to COM_call_function.
  COM_set_profiling(1);
  /*
  MPI_Barrier(MPI_COMM_WORLD);
  if (param.myRank == 0)
    std::cout << "Rocstar: profiling set." << std::endl;
  */

  bool with_pciter = param.maxNumPredCorrCycles > 1;
  if (param.myRank == 0 && verb > 1) {
    std::cout << "Rocstar: The maximum number of PC-iterations is "
              << param.maxNumPredCorrCycles << std::endl;
  }

  // in case of restarting, reset cur_step and current_time from Restart.txt
  coup->read_restart_info(param.current_time, param.cur_step);
  if (debug) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (param.myRank == 0)
      std::cout << "Rocstar: Read restart information, initializing"
                << std::endl;
  }

  // call input, init and run_initactions
  coup->init(param.current_time, param.time_step);
  /*
  std::cout << "Rocstar(" << param.myRank << "): Coupling initialized."
            << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  if (param.myRank == 0)
    std::cout << "Rocstar: All couplings initialized" << std::endl;
  */

  if (param.current_time == 0.0) {
    if (verb > 0) {
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: Performing time 0.0 dump..." << std::endl;
    }
    coup->update_integrals(param.current_time);
    coup->update_distances(param.current_time);
    coup->output_restart_files(param.current_time);
    if (debug) {
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: Performing time 0.0 dump...done." << std::endl;
    }
  }

  MPI_Barrier(param.communicator);
  param.LastOutputTime = param.current_time;
  param.LastOutputStep = param.cur_step;

  if (param.myRank == 0 && verb > 0) {
    std::cout << "Rocstar: Starting with step " << param.cur_step + 1
              << ", at time " << param.current_time << std::endl;
  }

  char header[100];
  sprintf(
      header,
      "************** Initialization times at time step %d *****************\n",
      param.cur_step);
  COM_print_profile(param.timingDataFile.c_str(), header);
  COM_set_profiling(1); // Reset profiler;

  if (param.current_time == 0.0)
    coup->write_restart_info(param.current_time, 1);

  // March through time until either the specified maximum time or the maximum
  // number of time steps have been exceeded (whichever comes first).
  bool InterfaceConverged = true;
  /*
  MPI_Barrier(param.communicator);
  if (param.myRank == 0)
    std::cout << "Rocstar: Preparing to step..." << std::endl;
  */
  while (!reached_simulation_time(param)) {
    param.cur_step++;

    for (int iPredCorr = 1; iPredCorr <= param.maxNumPredCorrCycles;
         iPredCorr++) {
      if (param.myRank == 0 && verb > 0) {
        std::cout << "Rocstar:\n"
                  << "Rocstar: " << std::string(80 - 9, '=') << "\n"
                  << "Rocstar: System Time Step : " << param.cur_step
                  << "   PC(" << iPredCorr << ")\n"
                  << "Rocstar: " << std::string(80 - 9, '=') << "\n"
                  << "Rocstar:\n"
                  << "Rocstar: CurrentTime, CurrentTimeStep";
        if (verb > 1)
          std::cout << ", ZoomFactor";
        std::cout << ": " << param.current_time << " " << param.time_step;
        if (verb > 1)
          std::cout << " " << param.zoomFactor;
        std::cout << "\n";
        std::cout << "Rocstar:" << std::endl;
      }

      /*
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: initializing convergence checking" << std::endl;
      */
      if (with_pciter)
        coup->init_convergence(iPredCorr);
      /*
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: convergence checking initialized" << std::endl;
      */

      // invoke time integration of coupling scheme by passing in
      // current time and desired time step and obtaining new time
      /*
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: Running coupling..." << std::endl;
      */
      double new_current_time = coup->run(param.current_time, param.time_step,
                                          iPredCorr, param.zoomFactor);

      /* AEG: disabling unused experimental code.
      // Check for interrupt
      check_for_interrupt(coup, param);
      */

      /*
      std::cout << "Rocstar(" << param.myRank << "): solver done on this rank."
                << std::endl;
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: coupling done." << std::endl;
      */

      InterfaceConverged = true;
      if (with_pciter)
        InterfaceConverged = coup->check_convergence();

      if (debug) {
        MPI_Barrier(param.communicator);
        if (param.myRank == 0)
          std::cout << "Rocstar:" << std::endl
                    << "Rocstar: iPredCorr = " << iPredCorr << " is done"
                    << std::endl;
      }

      if (InterfaceConverged) {
        param.current_time = new_current_time;
        if (verb > 1) {
          MPI_Barrier(param.communicator);
          if (param.myRank == 0)
            std::cout
                << "Rocstar: Success: predictor-corrector converged at time "
                << param.current_time << std::endl;
        }
        break;
      } else {
        if (debug) {
          MPI_Barrier(param.communicator);
          if (param.myRank == 0) {
            std::cout << "Rocstar:" << std::endl
                      << "Rocstar: Interface has -NOT- converged!" << std::endl
                      << "Rocstar: iPredCorr = " << iPredCorr << " is done"
                      << std::endl;
          }
        }
      }
    } // end of iPredCorr

    /* AEG: disabling unused experimental code.
    // New function to parse flags set by interrupt returns false if no action
    if (coup->ProcessInterrupt())
      continue;
    */

    if (!InterfaceConverged) {
      if (param.myRank == 0)
        COM_abort_msg(
            EXIT_FAILURE,
            "Rocstar: Disaster: predictor-corrector did not converge");
    }

    // if reached time for restart dump
    if (reached_restartdump_time(param)) {
      // MPI_Barrier(param.communicator);
      // Compute integrals and write them into files for conservation-check
      coup->update_integrals(param.current_time);
      coup->update_distances(param.current_time);

      // invoke restart output (as well as visualization data)
      if (verb > 0) {
        MPI_Barrier(param.communicator);
        if (param.myRank == 0)
          std::cout << "Rocstar: Dumping restart files... " << std::endl;
      }
      coup->output_restart_files(param.current_time);
      coup->write_restart_info(param.current_time, param.cur_step + 1);
      param.LastOutputTime = param.current_time;
      param.LastOutputStep = param.cur_step;
      param.current_dump++;
      /*
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: Dumping restart files...done." << std::endl;
      */
    } else if (reached_visdump_time(param)) {
      if (verb > 0) {
        MPI_Barrier(param.communicator);
        if (param.myRank == 0)
          std::cout << "Rocstar: Dumping coupling viz files." << std::endl;
      }
      // output visualization files
      coup->output_visualization_files(param.current_time);
    }

    // Write out profile data
    if (param.cur_step <= 100 || param.LastOutputTime == param.current_time) {
      sprintf(header,
              "************** Solver times up to time step %d since last "
              "output *********\n",
              param.cur_step);
      COM_print_profile(param.timingDataFile.c_str(), header);
      COM_set_profiling(1); // Reset profiler;
    }
    /*
    MPI_Barrier(param.communicator);
    if (param.myRank == 0)
      std::cout << "Rocstar:" << std::endl;
    */
  } // end time step

  // Write final solution.
  if (param.current_time != param.LastOutputTime) {
    if (verb > 0) {
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: Performing final dump..." << std::endl;
    }
    coup->output_restart_files(param.current_time);
    coup->write_restart_info(param.current_time,
                             std::min(param.cur_step, param.maxNumTimeSteps) +
                                 1);
    if (debug) {
      MPI_Barrier(param.communicator);
      if (param.myRank == 0)
        std::cout << "Rocstar: Performing final dump...done." << std::endl;
    }
  }

  // Finalize coupling
  coup->finalize();

#ifdef WITH_PANDA
  // unload Rocpanda to exit gracefully
  if (strcasecmp(param.output_module, "Rocpanda") == 0) {
    MAN_DEBUG(3, ("[%d] Rocstar: unload_module Rocpanda.\n", param.myRank));
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocpanda, "OUT");
  }
#endif

  // Write out the final profile data
  sprintf(header,
          "************** Finalization times after time step %d "
          "*****************\n",
          param.cur_step);
  COM_print_profile(param.timingDataFile.c_str(), header);

  delete coup;
}

// print double data of an dataitem on a pane
void debug_print(const std::string &str, int pane, int pe, MPI_Comm comm,
                 const char *memo) {
  if (COMMPI_Comm_rank(comm) == pe) {
    double *vm = nullptr;
    int strid, cap;
    printf("%s %s: before %p\n", str.c_str(), memo ? memo : "",
           static_cast<void *>(vm));
    COM_get_array(str.c_str(), pane, &vm, &strid, &cap);
    printf("%s %s: after  %p\n", str.c_str(), memo ? memo : "",
           static_cast<void *>(vm));
    for (int i = 0; i < strid * cap; i++)
      printf("%.17e ", vm[i]);
    printf("\n");
  }
}

void debug_int_print(const std::string &str, int pane, int pe, MPI_Comm comm,
                     const char *memo) {
  if (COMMPI_Comm_rank(comm) == pe) {
    int *vm = nullptr;
    int strid, cap;
    printf("%s %s: before %p\n", str.c_str(), memo ? memo : "",
           static_cast<void *>(vm));
    COM_get_array(str.c_str(), pane, &vm, &strid, &cap);
    printf("%s %s: after  %p\n", str.c_str(), memo ? memo : "",
           static_cast<void *>(vm));
    for (int i = 0; i < strid * cap; i++)
      printf("%d ", vm[i]);
    printf("\n");
  }
}
