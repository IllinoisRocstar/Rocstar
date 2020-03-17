//
// Created by agondolo on 2/7/20.
//

#ifndef ROCSTAR_ROCSTARACTION_H_
#define ROCSTAR_ROCSTARACTION_H_

#include "Action.h"

class RocstarAction : public Action {
public:
  RocstarAction(int n, const char **at, int *i = nullptr, void *p = nullptr,
                const char *name = nullptr);
  RocstarAction(int n, const std::string at[], int *i = NULL, void *p = 0,
                const char *name = NULL);

public:
  void finalize() override{};

protected:
  void set_attr(int n, const char *at[], const int *id = nullptr);
  void set_attr(int n, const std::string at[], const int *id = nullptr);

  // Set I/O properties of arguments
  void set_io(int n, const int *io);
  int get_io(int i) { return action_data[i].inout; }
};

#endif // ROCSTAR_ROCSTARACTION_H_
