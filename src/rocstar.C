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

#include <iostream>
#include <cstring>
#include <cstdlib>
//#include <mpi.h>
#include "com.h"
#include "splash.h"

//##include "makeopts.h"

using namespace std;

#ifdef ROCPROF
#include "Rocprof.H"
#endif

void rocstar_driver( int, int, bool);

struct Command_Options {
  Command_Options() : verbose(1),nrun(0),debug(false) {}

  int  verbose;
  int  remeshed;
  int  nrun;
  bool debug;
};

void
ReportUsage()
{
  std::cout << "rocstar [-h -d -v [0-2] -r [n] ]" << std::endl
	    << std::endl
	    << "Rocstar Help" << std::endl
	    << "------------------------" << std::endl
	    << " -h: help - prints this message" << std::endl
	    << " -v: verbosity - optional verbosity level (default = 1)" << std::endl
	    << " -r: runs - number of times to automatically restart (default = 0)" << std::endl
            << " -d: debug - turns on debugging messages for the run" << std::endl
	    << "------------------------" << std::endl;

}

void
parse_commandline( int argc, char *argv[], Command_Options &opts) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  for ( int i=1; i<argc; ++i) {
    if ( std::strcmp(argv[i], "-v") == 0) {
      if ( argc>i+1 && argv[i+1][0]>='0' && argv[i+1][0]<='9')
      { opts.verbose = std::atoi( argv[i+1]); ++i; }
    }
    else if ( std::strcmp(argv[i], "-remeshed") == 0) {
      opts.remeshed = 1;
    }
    else if ( std::strcmp(argv[i], "-r") == 0) {
      opts.nrun = 1;
      if ( argc>i+1 && argv[i+1][0]>='1' && argv[i+1][0]<='9')
      { opts.nrun = std::atoi( argv[i+1]); ++i; }
    }
    else if ( std::strcmp(argv[i], "-d") == 0){
      opts.debug = true;
    }
    else if ( std::strcmp(argv[i], "-debug") == 0){
      opts.debug = true;
    }
    else if (rank == 0) {
      if(std::strcmp(argv[i],"-h"))
	std::cerr << "Rocstar: Unknown option " << argv[i] << std::endl;
      ReportUsage();
    }
  }
  if(opts.verbose > 0 && rank == 0)
    std::cout << RocstarSplash << std::endl;
}

int
main( int argc, char *argv[]) {


  /* Start MPI. */
  MPI_Init( &argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef ROCPROF
  Rocprof_Init("Rocstar",rank);
#endif
  /* Initialize Roccom */
  COM_init( &argc, &argv);

  Command_Options opts;
  parse_commandline( argc, argv, opts);

  int nrun = 1;
  rocstar_driver( opts.verbose, opts.remeshed, opts.debug);
  while (nrun < opts.nrun){
    if(rank == 0) std::cout << " Rocstar automatic restart " << nrun++ << std::endl;
    rocstar_driver(opts.verbose,opts.remeshed, opts.debug);
  }
  /* Finalize Roccom. */
  COM_finalize();

#ifdef ROCPROF
  Rocprof_Finalize(true);
#endif
  /* Close down MPI. */
  MPI_Finalize();

  return 0;
}






