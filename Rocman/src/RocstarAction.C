#include "RocstarAction.h"
#include "com_assertion.h"

RocstarAction::RocstarAction(int n, const char **at, int *id, void *p,
                             const char *name)
    : Action(name) {
  for (int i = 0; i < n; ++i) {
    action_data.emplace_back(at ? at[i] : "", id ? id[i] : 0,
                             ActionDataIntent::IN);
  }
}

RocstarAction::RocstarAction(int n, const std::string *at, int *id, void *p,
                             const char *name)
    : Action(name) {
  for (int i = 0; i < n; ++i) {
    action_data.emplace_back(at ? at[i] : "", id ? id[i] : 0,
                             ActionDataIntent::IN);
  }
}

void RocstarAction::set_attr(int n, const char **at, const int *id) {
  COM_assertion_msg(n == static_cast<int>(action_data.size()), name().c_str());
  for (int i = 0; i < n; ++i) {
    action_data[i].attr = at[i];
    action_data[i].idx = id ? id[i] : 0;
  }
}

void RocstarAction::set_attr(int n, const std::string *at, const int *id) {
  COM_assertion_msg(n == static_cast<int>(action_data.size()), name().c_str());
  for (int i = 0; i < n; ++i) {
    action_data[i].attr = at[i];
    action_data[i].idx = id ? id[i] : 0;
  }
}

void RocstarAction::set_io(int n, const int *io) {
  COM_assertion_msg(n == static_cast<int>(action_data.size()), name().c_str());
  for (int i = 0; i < n; ++i)
    action_data[i].inout = static_cast<ActionDataIntent>(io[i]);
}