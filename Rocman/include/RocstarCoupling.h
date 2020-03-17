//
// Created by agondolo on 9/4/18.
//

#ifndef _ROCSTAR_ROCSTARCOUPLING_H_
#define _ROCSTAR_ROCSTARCOUPLING_H_

#include "Coupling.h"

#include "BurnAgent.h"
#include "FluidAgent.h"
#include "SolidAgent.h"
#include "UserScheduler.h"

class RocstarAgent;
class Control_parameters;
class RocmanControl_parameters;

class RocstarCoupling : public Coupling {
public:
  /// Constructor. Derived class will add actions for the coupling scheme
  RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                  const Control_parameters *p,
                  const RocmanControl_parameters *mp);
  RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                  const Control_parameters *p,
                  const RocmanControl_parameters *mp, const std::string &name);
  RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                  const Control_parameters *p,
                  const RocmanControl_parameters *mp,
                  const std::string &fluidname, const std::string &solidname);
  RocstarCoupling(const std::string &coupl_name, MPI_Comm com,
                  const Control_parameters *p,
                  const RocmanControl_parameters *mp,
                  const std::string &fluidname, const std::string &solidname,
                  const std::string &burnname);

  const Control_parameters *get_control_param() const { return param; }
  const RocmanControl_parameters *get_rocmancontrol_param() const {
    return rocmanparam;
  }

  /*
  /// Invoke initialization of the actions in the scheduler and the agents
  void initialize(int reinit = 0);

  static std::string normalize_modname(const std::string &mod);

  // Experimental Interrupt
  void Interrupt(int *, const char *);
  int ProcessInterrupt();
  void restart_at_time(double t, int step);
  */

  int new_start(double t) const override;

  // Virtual methods used in FullyCoupling
  virtual void reload_rocface(int rfc_verb){};
  virtual void update_integrals(double currentTime) {}
  virtual void update_distances(double currentTime) {}

protected:
  typedef std::vector<RocstarAgent *> RocstarAgentList;

  const Control_parameters *param;
  const RocmanControl_parameters *rocmanparam;

  RocstarAgentList agents;
};

// subclass fully coupled involving fluid and solid, and maybe burn
class FullyCoupling : public RocstarCoupling {
protected:
  FluidAgent *fluid_agent;
  SolidAgent *solid_agent;
  BurnAgent *burn_agent;

public:
  FullyCoupling(const std::string &coupl_name, MPI_Comm com,
                const Control_parameters *p, const RocmanControl_parameters *mp,
                const std::string &fluidname, const std::string &solidname);
  FullyCoupling(const std::string &coupl_name, MPI_Comm com,
                const Control_parameters *p, const RocmanControl_parameters *mp,
                const std::string &fluidname, const std::string &solidname,
                const std::string &burnname);

  void update_integrals(double currentTime) override;
  void update_distances(double currentTime) override;
  void reload_rocface(int rfc_verb) override;

private:
  // compute integrals
  bool overwrite_integ{false};
  std::string integFname;

  bool overwrite_dist{false};
  std::string distFname;
};

#define DECLARE_NEW_COUPLING_SCHEME_NO_BURN(New_scheme)                        \
  class New_scheme : public RocstarCoupling {                                  \
  public:                                                                      \
    New_scheme(MPI_Comm com, const Control_parameters *p,                      \
               const RocmanControl_parameters *mp, const std::string &);       \
  }

#define DECLARE_NEW_COUPLING_SCHEME_WITH_BURN(New_scheme)                      \
  class New_scheme : public RocstarCoupling {                                  \
  public:                                                                      \
    New_scheme(MPI_Comm com, const Control_parameters *p,                      \
               const RocmanControl_parameters *mp, const std::string &,        \
               const std::string &);                                           \
  }

#define DECLARE_NEW_FULLY_COUPLING_SCHEME_NO_BURN(New_scheme)                  \
  class New_scheme : public FullyCoupling {                                    \
  public:                                                                      \
    New_scheme(MPI_Comm com, const Control_parameters *p,                      \
               const RocmanControl_parameters *mp, const std::string &,        \
               const std::string &);                                           \
  }

#define DECLARE_NEW_FULLY_COUPLING_SCHEME_WITH_BURN(New_scheme)                \
  class New_scheme : public FullyCoupling {                                    \
  public:                                                                      \
    New_scheme(MPI_Comm com, const Control_parameters *p,                      \
               const RocmanControl_parameters *mp, const std::string &,        \
               const std::string &, const std::string &);                      \
  }

#endif //_ROCSTAR_ROCSTARCOUPLING_H_
