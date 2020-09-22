#include "RocstarCoupling.h"

#include <cmath>
#include <iostream>

#include "Control_parameters.h"
#include "rocman.h"

/*
COM_EXTERN_MODULE(SurfMap)
COM_EXTERN_MODULE(SurfX)
COM_EXTERN_MODULE(SurfUtil)
COM_EXTERN_MODULE(Simpal)
COM_EXTERN_MODULE(Rocprop)
*/

/*
void RocstarShutdown(int = 0);
#ifdef ARM
int Rocrem_remesh(const char *solver, const char *indir, const char *outdir,
                  double time, bool pt, bool remesh_surf, bool transfer_surf,
                  double scaleFactor, MPI_Comm myComm, int fieldWidth,
                  int fileBase, int debug, int ngLayers);
#endif
*/

RocstarCoupling::RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                                 const Control_parameters *p,
                                 const RocmanControl_parameters *mp)
    : RocstarCoupling(coupl_name, com, p, mp, "", "", "") {}

RocstarCoupling::RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                                 const Control_parameters *p,
                                 const RocmanControl_parameters *mp,
                                 const std::string &name)
    : RocstarCoupling(coupl_name, com, p, mp, name, "", "") {}

RocstarCoupling::RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                                 const Control_parameters *p,
                                 const RocmanControl_parameters *mp,
                                 const std::string &fluidname,
                                 const std::string &solidname)
    : RocstarCoupling(coupl_name, com, p, mp, fluidname, solidname, "") {}

RocstarCoupling::RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                                 const Control_parameters *p,
                                 const RocmanControl_parameters *mp,
                                 const std::string &fluidname,
                                 const std::string &solidname,
                                 const std::string &burnname)
    : Coupling(coupl_name, com), param(p), rocmanparam(mp) {
  if (!fluidname.empty())
    modules.push_back(fluidname);
  if (!solidname.empty())
    modules.push_back(solidname);
  if (!burnname.empty())
    modules.push_back(burnname);

  init_started = param->current_time == 0.0;
  init_remeshed = rocmanparam->remeshed;
  restartInfo = "Restart.txt";

  /* AEG: disabling unused experimental code.
  // experimental Interrupt
  std::string manwinname = "Rocman";
  std::string manfuncname = manwinname + ".interrupt";
  std::string coupobjname = manwinname + ".coup";
  COM_new_window(manwinname);
  COM_new_dataitem(coupobjname, 'w', COM_VOID, 1, "");
  COM_set_object(coupobjname, 0, this);
  COM_set_member_function(
      manfuncname, (Member_func_ptr)&RocstarCoupling::Interrupt, coupobjname,
      "bii", {COM_RAWDATA, COM_INT, COM_STRING});
  */
}

int RocstarCoupling::new_start(double t) const {
  return Coupling::new_start(t) ||
         (rocmanparam->remeshed && t == param->init_time);
}

/* AEG: Functionality moved to parent class.
void RocstarCoupling::initialize(int reinit) {
  // reset rocface if restart
  if (reinit)
    reload_rocface(rocmanparam->rfc_verb);
}
*/

/* AEG: disabled unused method
std::string RocstarCoupling::normalize_modname(const std::string &modname) {
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
*/

/* AEG: disabling unused experimental code.
// experimental code
void RocstarCoupling::restart_at_time(double t, int step) {
  int i, n;

  restarting = true;

  // reset windows and buffers
  finalize();

  // unload modules, only call unload to zero out data, skip dlclose
  for (i = 0, n = agents.size(); i < n; ++i) {
    agents[i]->unload_module();
  }

  if (COM_get_window_handle("RFC") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SurfX, "RFC");
  }
  if (COM_get_window_handle("SURF") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SurfUtil, "SURF");
  }
  // if (COM_get_window_handle("PROP") > 0) {
  //   COM_UNLOAD_MODULE_STATIC_DYNAMIC( SurfProp, "PROP");
  // }
  // if (COM_get_window_handle("PROPCON") > 0) {
  //   COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocon,"PROPCON");
  // }
  if (COM_get_window_handle("BLAS") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(Simpal, "BLAS");
  }
  if (COM_get_window_handle("MAP") > 0) {
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SurfMap, "MAP");
  }
  // update to given time
  //  param->update_start_time( step, t);

  // load module and call physics routine
  for (i = 0, n = agents.size(); i < n; ++i) {
    agents[i]->load_module();
  }

  // proceed to initialization phase
  initialize(1);    // reinit


  // read restart data
  // for ( i=0, n=agents.size(); i<n; ++i) {
  //   agents[i]->read_restart_data();
  // }

  restarting = false;
}

// experimental code
void RocstarCoupling::Interrupt(int *act, const char *message) {
  int action = *act;
  if (!param->myRank) {
    std::cout << "Rocman::Interrupt invoked." << std::endl
              << "Rocman: ***************************************************"
              << std::endl
              << "Rocman: " << (message ? message : "") << std::endl
              << "Rocman: ***************************************************"
              << std::endl;
  }
  switch (action) {
    case 0:
      // Just stop
      if (!param->myRank)
        std::cout << "Rocman: Halting simulation." << std::endl;
      finalize();
      //    RocstarShutdown();
      break;
    case 1:
      // Output solutions and stop.
      if (!param->myRank)
        std::cout << "Rocman: Writing restart files and halting simulation."
                  << std::endl;
      output_restart_files(param->current_time);
      finalize();
      //    RocstarShutdown();
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
      if (!param->myRank) {
        std::cout << "Rocman: Directed to warm restart from last dump "
                  << "and increase dump frequency." << std::endl;
      }
      param->InterruptFlag = 2;
      break;
    case 3:
      if (!param->myRank)
        std::cout << "Rocman: Directed to remesh and warm restart from "
                  << "(time/step): ("
                  << param->LastOutputTime << "," << param->LastOutputStep
                  << ")" << std::endl;
      param->InterruptFlag = 3;
      break;
    case 4:
      if (!param->myRank)
        std::cout << "Rocman: Directed to restart from "
                  << "(time/step): (" << param->LastOutputTime << ","
                  << param->LastOutputStep << ")" << std::endl;
      if (message)
        std::cout << "MESSAGE: " << message << std::endl;
      param->InterruptFlag = 4;
      break;
    case 5:
      if (!param->myRank)
        std::cout << "Rocman: Directed to dump at (time/step): ("
                  << param->current_time << "," << param->cur_step
                  << ")" << std::endl;
      output_restart_files(param->current_time);
      write_restart_info(param->current_time, param->cur_step + 1);
      break;
    case 6:
      if (!param->myRank)
        std::cout << "Rocman: Directed to dump and restart at (time/step): ("
                  << param->current_time << "," << param->cur_step
                  << ")" << std::endl;
      output_restart_files(param->current_time);
      write_restart_info(param->current_time, param->cur_step + 1);
      if (message)
        std::cout << "MESSAGE: " << message << std::endl;
      param->InterruptFlag = 4;
      break;
    default:
      if (param->myRank == 0)
        std::cout << "Rocman: Unknown interrupt action." << std::endl;
      return;
  }
}

int RocstarCoupling::ProcessInterrupt() {
  int save_step = 0;
  bool done = false;
  int fluid_agent_index = 0;
  int action = param->InterruptFlag;
  param->InterruptFlag = 0;
//  if (action > 0) {
//    if (param->myRank == 0)
//      std::cout << "Rocman: Processing interrupt." << std::endl;
//    switch (action) {
//      case 2:save_step = param->cur_step;
//        restart_at_time(param->LastOutputTime, param->LastOutputStep);
//        param->maxNumTimeSteps = save_step;
//        param->outputIntervalTime /= 10.0;
//        if (!param->myRank) {
//          std::cout << "Rocman: Restarting from (time/step): ("
//                    << param->current_time << "," << param->cur_step
//                    << ") with output"
//                    << " interval, " << param->outputIntervalTime << "."
//                    << std::endl;
//        }
//        return (1);
//        break;
//      case 3:
//        for (int i = 0, n = agents.size(); i < n && !done; n++) {
//          std::string::size_type x = agents[i]->get_agent_name().find("lu");
//          if (x != string::npos) {
//            fluid_agent_index = i;
//            done = true;
//          }
//        }
//        if (!done) {
//          std::cerr << "Rocman: Could not find fluid agent.  Dying."
//                    << std::endl;
//        // RocstarShutdown(1);
//        }
//        // char* solver: Path to SimOUT directory (i.e. blah/blah/Rocflu)
//        // double time: timestep to read
//        // bool use_parallel_transfer
//        // bool remesh_surf: true=remesh surface false=leave surface alone
//        // bool transfer_surf: true (only valid if remesh_surf = false)
//        // double scaleFactor: scale mesh granularity
//        // int fieldWidth: digits in filename partition id
//        // int fileBase: 0 is partition ids start at 0, 1 otherwise
//        // int debug: 0=off, 1 = some, more = more
//        // int ngLayers: number of ghost layers; 0-9
//        // returns 1 if success; 0 if fail
//        // Rocrem_remesh(const char* solver, double time, bool remesh_surf,
//        //               bool transfer_surf, double scaleFactor, MPI_Comm myComm,
//        //               int fieldWidth, int fileBase, int debug, int ngLayers)
//        //
//        // Remeshing test case
//        //
//#ifdef ARM
//        if (!Rocrem_remesh("Rocflu", "SimOUT", "", param->LastOutputTime, true, true, false, 1.0,
//                           agents[fluid_agent_index]->get_communicator(),
//                           4, 0, 2, 2
//        )) {
//          if (!param->myRank)
//            std::cout << "Rocman: Remeshing failed.  Stopping simulation."
//                      << std::endl;
//          // RocstarShutdown(1);
//        }
//#else
//        if (!param->myRank) {
//          // std::ofstream Ouf;
//          // Ouf.open("needs_remesh");
//          // Ouf.close();
//          std::cout << "Rocman: Online remeshing not enabled. Shutting down." << std::endl;
//        }
//        // RocstarShutdown(1);
//#endif
//        MPI_Barrier(MPI_COMM_WORLD);
//        restart_at_time(param->LastOutputTime, param->LastOutputStep);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (!param->myRank)
//          std::cout << "Rocman: Restarting after remesh at (time/step): ("
//                    << param->current_time << "," << param->cur_step << ")"
//                    << std::endl;
//        return (1);
//        break;
//      case 4:MPI_Barrier(MPI_COMM_WORLD);
//        restart_at_time(param->LastOutputTime, param->LastOutputStep);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (!param->myRank)
//          std::cout << "Rocman: Warm restarting at (time/step): ("
//                    << param->current_time << "," << param->cur_step << ")"
//                    << std::endl;
//        return (1);
//        break;
//      default:
//        if (!param->myRank)
//          std::cout << "Rocman: Unknown interrupt action: " << action << "."
//                    << std::endl;
//        return (0);
//        break;
//    }
//  }
  return (0);
}
*/

// Compute the integrals and dump into a file.
void FullyCoupling::update_integrals(double currentTime) {
  int count = 0;
  for (auto &&agent : agents) {
    count += agent->compute_integrals();
  }

  if (comm_rank == 0 && count == 2) {
    FILE *fp = fopen(integFname.c_str(), overwrite_integ ? "w" : "a");
    if (fp == nullptr) {
      COM_abort_msg(EXIT_FAILURE, "Rocman: Failed to open integral file, " +
                                      integFname + ".");
    }
    MAN_DEBUG(2, ("Rocman: FullyCoupling::update_integrals with t: %f.\n",
                  currentTime));
    if (overwrite_integ) {
      fprintf(fp, "time                     "
                  "f-volume                 "
                  "s-volume                 "
                  "f-mass                   "
                  "s-mass                   "
                  "f-burn area              "
                  "s-burn area              "
                  "f-non-burn               "
                  "area                     "
                  "s-non-burn               "
                  "area                     "
                  "s-volume-undef"
                  "\n");
      overwrite_integ = false;
    }
    double *f_integrals = fluid_agent->get_integrals();
    double *s_integrals = solid_agent->get_integrals();
    fprintf(fp,
            "%.18le %.18le %.18le %.18le %.18le %.18le %.18le %.18le %.18le "
            "%.18le\n",
            currentTime, f_integrals[MAN_INTEG_VOL], s_integrals[MAN_INTEG_VOL],
            f_integrals[MAN_INTEG_MASS], s_integrals[MAN_INTEG_MASS],
            f_integrals[MAN_INTEG_IBAREA], s_integrals[MAN_INTEG_IBAREA],
            f_integrals[MAN_INTEG_INBAREA], s_integrals[MAN_INTEG_INBAREA],
            f_integrals[MAN_INTEG_VOL_UND]);

    fclose(fp);
  }
}

// Compute the distances between fluid and solid meshes and dump into a file.
void FullyCoupling::update_distances(double currentTime) {
  double dist_max = -1, dist_min = -1, dist_nrm2 = -1;

//  MPI_Comm communicator = fluid_agent->get_communicator();

//  int RFC_interpolate = COM_get_function_handle("RFC.interpolate");
//  COM_call_function(RFC_interpolate, &solid_agent->y_hdl, &fluid_agent->nc_tmp_hdl);
//  COM_call_function(RocBlas::sub, &fluid_agent->nc_tmp_hdl, &fluid_agent->nc_hdl, &fluid_agent->nc_tmp_hdl);
//  COM_call_function(RocBlas::nrm2, &fluid_agent->nc_tmp_hdl, &fluid_agent->sq_dist_hdl);
//  COM_call_function(RocBlas::max_scalar_MPI, &fluid_agent->sq_dist_hdl, &dist_max, &communicator);
//  COM_call_function(RocBlas::min_scalar_MPI, &fluid_agent->sq_dist_hdl, &dist_min, &communicator);
//  COM_call_function(RocBlas::sum_scalar_MPI, &fluid_agent->sq_dist_hdl, &dist_nrm2, &communicator);

  if (comm_rank == 0) {
    FILE *fp;
    fp = fopen(distFname.c_str(), overwrite_dist ? "w" : "a");
    if (fp == NULL) {
      std::cout << "Rocman: Failed to open distance file, " << distFname
                << "." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MAN_DEBUG(2, ("Rocman: FullyCoupling::update_distances with t: %f.\n", currentTime));
    if (overwrite_dist) {
      fprintf(fp, "# time    distance-min    distance-max    distance-norm2  \n");
      overwrite_dist = false;
    }
    fprintf(fp, "%.18le %.18le %.18le %.18le\n",
            currentTime, sqrt(dist_min), sqrt(dist_max), sqrt(dist_nrm2));

    fclose(fp);
  }
}

//void _load_rocface(FluidAgent *fagent, SolidAgent *sagent, int rfc_verb);

void FullyCoupling::reload_rocface(int rfc_verb) {
//    _load_rocface(fluid_agent, solid_agent, rfc_verb);
}

FullyCoupling::FullyCoupling(const std::string &coupl_name, MPI_Comm com,
                             const Control_parameters *p,
                             const RocmanControl_parameters *mp,
                             const std::string &fluidname,
                             const std::string &solidname)
    : FullyCoupling(coupl_name, com, p, mp, fluidname, solidname, "") {}

FullyCoupling::FullyCoupling(const std::string &coupl_name, MPI_Comm com,
                             const Control_parameters *p,
                             const RocmanControl_parameters *mp,
                             const std::string &fluidname,
                             const std::string &solidname,
                             const std::string &burnname)
    : RocstarCoupling(coupl_name, com, p, mp, fluidname, solidname, burnname),
      fluid_agent(nullptr), solid_agent(nullptr), burn_agent(nullptr) {
  integFname = "Rocman/Modout/ROCSTAR_integ.txt";
  distFname = "Rocman/Modout/ROCSTAR_dist.txt";
  if (!init_started) {
    FILE *fp = fopen(integFname.c_str(), "r");
    if (!fp)
      overwrite_integ = true;
    fp = fopen(distFname.c_str(), "r");
    if (!fp)
      overwrite_dist = true;
  } else {
    overwrite_integ = true;
    overwrite_dist = true;
  }
}
