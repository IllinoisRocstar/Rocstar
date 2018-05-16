
#include "mpi.h"
#include "roccom.h"

extern "C" {
  MPI_Comm GetCommunicator()
  {
    return(COM_get_default_communicator());
  }
}
