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
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <cassert>
#include <cstdlib>

using namespace std;

//#ifdef _TRAIL_MPI_
#include "mpi.h"
//#ifdef _ROCSTAR_X_
#include "com.h"
// COM_EXTERN_MODULE(Rocout);
//#endif
//#endif

#include "GEM.H"
#include "TRAIL_Flu.H"
#include "TRAIL_GeoPart.H"


int
main(int argc,char *argv[])
{
  bool debug = true;
  bool map_wrote = false;
  if(argc < 4){
    cerr << "part2flu <source_prefix> <target_prefix> <npartitions>" << endl;
    return(1);
  }
  int nproc = 1;
  int rank = 0;
  //#ifdef _TRAIL_MPI_
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //#ifdef _ROCSTAR_X_
  COM_init(&argc,&argv);
  //#endif
  //#endif
  std::ofstream DebugFile;
  std::ostringstream DBFileName;
  DBFileName << "gp2rstar_" << rank;
  DebugFile.open(DBFileName.str().c_str());
  bool verbose = true;
  string source_prefix(argv[1]);
  string targ_prefix(argv[2]);
  unsigned int npart;
  string snpart(argv[3]);
  istringstream Istr(snpart);
  Istr >> npart;
  unsigned int N;
  int R;
  unsigned int first;

  if(npart < (unsigned int)nproc){
    N = 1;
    first = rank;
  }
  else{
    N = npart/nproc;
    R = npart - (nproc*N);
    if(rank < R){
      first = rank * N + rank;
      N++;
    }
    else 
      first = rank * N + R;
  }   
  if(debug)
    DebugFile << "Processing " << npart << " partitions." << std::endl
	      << "Local partitions: (" << N << "," << first 
	      << ")" << std::endl;
  int part_index = first;
  unsigned int n = 0;
  // Now that we know which partitions we will process on this processor, the
  // task at hand is to populate the GEM_ structures from Bill's output. That's
  // all the following code does.
  while(n < N){
    if(!rank && verbose)
      std::cout << "Processing " << part_index << "." << std::endl;
    if(debug)
      DebugFile << "Processing " << part_index << "." << std::endl;
    std::ifstream Inf;
    std::ostringstream FileName;
    FileName << source_prefix << "." << part_index << ".asc";
    Inf.open(FileName.str().c_str());
    if(!Inf){
      std::cerr << "Cannot open " << FileName.str() << ", aborting." << std::endl;
      return(false);
    }
    GEM_Partition my_partition(part_index+1);
    my_partition._cell_ordering[0] = 1;
    my_partition._cell_ordering[1] = 4;
    my_partition._cell_ordering[2] = 3;
    my_partition._cell_ordering[3] = 2;
    my_partition._npart = npart;
    if(!TRAIL_FluInitSolver(my_partition,"Rocflu")){
      std::cerr << rank << ": TRAIL_FluInitSolver() failed, aborting." << std::endl;
      return(1);
    }
    if(debug){
      my_partition.debug();
      my_partition._out = &DebugFile;
    }
    //#ifdef _TRAIL_MPI_
    my_partition._comm = MPI_COMM_WORLD;
    //#endif
    if(!TRAIL_GeoPartReadASCII(my_partition,Inf)){
      cerr << rank << ": Cannot process " << FileName.str() << ", aborting." << endl;
      return(1);
    }
    if(verbose && !rank)
      cout << "Done reading Bill's partition, populating Rocflu data." << endl;
    if(debug)
      DebugFile << "Done reading Bill's partition, populating Rocflu data." << std::endl;

    // This function populates the Rocflu boundary patches.  You'll need a
    // *.cgi file for this (needed when translating a grid generator's BC id
    // to Rocflu patch id's)  
    if(!TRAIL_FluPopulatePatches(my_partition))
      std::cerr << rank << ": TRAIL_FluPopulatePatches failed, aborting." << std::endl;
    if(debug)
      DebugFile << "Patches done." << std::endl;


    // When populating from scratch (like we do with Bill's meshes), we need
    // to initialize actual solutions on the mesh.  We don't need to do this
    // after remeshing, or if data already exists.
    TRAIL_FluInitVolSoln(my_partition);
    if(debug)
      DebugFile << "Volume solution done." << std::endl;
    TRAIL_FluInitSurfSoln(my_partition);
    if(debug){
      DebugFile << "Surface solution done." << std::endl;
      my_partition.report();
      DebugFile << "Writing Rocflu native data." << std::endl;
    }

    // Writes the Rocflu native text files that are required by both Rocflu
    // standalone and Rocstar
    if(!TRAIL_FluWriteNative(my_partition,"Rocflu/Modin")){
      std::cerr << rank << ": TRAIL_FluWriteNative failed, aborting." << std::endl;
      return(1);
    }
    // Write this file only once
    if(!rank && !map_wrote){
      TRAIL_FluWriteMAP(my_partition,npart,npart,"Rocflu/Modin");
      map_wrote = true;
    }
    // FIXME! - Still need to populate remote border indices!

    // FIXME
    // Writes the Grid and Soln files needed by standalone Rocflu
    //    fluregion.WriteFluGridASCII(targ_prefix,0.0);
    //    fluregion.WriteFluSolnASCII(targ_prefix,0,0.0,true);
    // FIXME

    if(debug){
      DebugFile << "Done writing native files." << std::endl
		<< "Writing Rocstar data." << std::endl;
    }
    //#ifdef _ROCSTAR_X_
    if(debug)
      COM_set_verbose(6);
    // The next few calls set up, and populate a Rocflu window
    my_partition.InitRoccomWindows("flu");
    if(debug && !rank)
      cout << "  window created." << endl;
    TRAIL_FluRegisterSurfMesh(my_partition);
    TRAIL_FluRegisterVolSoln(my_partition);
    TRAIL_FluRegisterSurfSoln(my_partition);
    my_partition.WindowInitDone();
    if(debug && !rank)
      cout << "  window data registered." << endl;
    my_partition.WriteRocstar("Rocflu/Rocin",0.0);
    // These calls identify the target directory and write all the files
    // Rocstar needs.
//     string::size_type y = targ_prefix.find_last_of("/");
//     string targ_dir = targ_prefix.substr(0,y);
//     fluregion.WriteRocstar(targ_dir);
    if(debug && !rank)
      cout << "  hdf files written." << endl;
    
    // Clean up after ourselves in Roccom
    my_partition.DestroyWindows();
    if(debug && !rank)
      cout << "  window data destroyed." << endl;
    //#endif
    if(debug && !rank)
      cout << "Local partition " << n << " processed." << endl;
    n++;
    part_index++;
  }
  if(verbose && !rank)
    std::cout 
      << "All local partitions processed. Waiting for other processors...";
  MPI_Barrier(MPI_COMM_WORLD);
  if(verbose && !rank)
    std::cout << "done." << std::endl;
  
  //#ifdef _TRAIL_MPI_
  MPI_Finalize();
  //#endif
  
    //#ifdef _ROCSTAR_X_
  COM_finalize();
  //#endif
  return(0);
}








