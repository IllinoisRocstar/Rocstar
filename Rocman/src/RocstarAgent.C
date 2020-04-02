#include <algorithm>
#include <utility>

#include "RocstarAgent.h"

#include "Control_parameters.h"
#include "RocstarCoupling.h"
#include "rocman.h"

RocstarAgent::RocstarAgent(RocstarCoupling *cp, std::string mod,
                           std::string obj, std::string agentname,
                           std::string volname, std::string surfname,
                           MPI_Comm com, bool wgm)
    : Agent(cp, std::move(agentname), com, std::move(mod), std::move(obj),
            std::move(surfname), std::move(volname)),
      coupling(cp) {
  MAN_DEBUG(1, ("Rocman: Agent::Agent create window %s.\n",
                get_agent_name().c_str()));

  timestamp = cp->get_control_param()->current_time;
  with_gm = wgm;
  dobackup = true;
}

/*
RocstarAgent::~RocstarAgent() {
}
*/

void RocstarAgent::init_subscheduler(double t) {
  MAN_DEBUG(1, ("Rocman: %sAgent::init_subscheduler called.\n",
                get_agent_name().c_str()));

  callMethod(&UserScheduler::init_actions, t);

  MAN_DEBUG(1, ("Rocman: %sAgent::init_subscheduler done.\n",
                get_agent_name().c_str()));
}
