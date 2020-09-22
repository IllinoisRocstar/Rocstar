#include "rocman.h"

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
