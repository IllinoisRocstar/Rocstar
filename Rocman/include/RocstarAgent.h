//
// Created by agondolo on 8/31/18.
//

#ifndef _ROCSTAR_ROCSTARAGENT_H_
#define _ROCSTAR_ROCSTARAGENT_H_

#include "Agent.h"
#include "UserScheduler.h"

#define MAN_INTEG_SIZE 9

class RocstarCoupling;

class RocstarAgent : public Agent {
public:
  RocstarAgent(RocstarCoupling *cp, std::string mod, std::string obj,
               std::string agent_name, std::string vol_name,
               std::string surf_name, MPI_Comm com, bool wgm = false);

  void init_subscheduler(double t);

  virtual int compute_integrals() { return 0; }
  double *get_integrals() { return integrals; }

  RocstarCoupling *get_rocstar_coupling() { return coupling; };

  //  void set_dobackup(bool dobackup_) { dobackup = dobackup_; };

protected:
  RocstarCoupling *coupling;

  int compute_integrals_handle; ///< fluid code handles

  double integrals[MAN_INTEG_SIZE]; ///< Array for storing integrals
};

#endif //_ROCSTAR_ROCSTARAGENT_H_
