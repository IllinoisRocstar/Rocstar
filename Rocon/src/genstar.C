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

#include "GEM.H"
//#include "Partition.H"
#include "FluRegion.H"

//#ifdef _TRAIL_MPI_
#include "mpi.h"
//#ifdef _ROCSTAR_X_
#include "roccom.h"
COM_EXTERN_MODULE(Rocout);
//#endif
//#endif

int
main(int argc,char *argv[])
{
  if(argc < 4){
    cerr << "genstar <source_prefix> <target_prefix> <npartitions>" << endl;
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
  if(!rank)
    cout << rank << ": Processing " << npart << " partitions." << endl;
  if(verbose)
    cout << rank << ": Local partitions: (" << N << "," << first 
	 << ")" << endl;
  int part_index = first;
  unsigned int n = 0;
  while(n < N){
    if(!rank && verbose)
      cout << "Processing " << part_index++ << "." << endl;
    FluRegion fluregion;
    fluregion.debug();
    bool unsteady = true;
    double t = 0.0;
    unsigned int niter = 0;
    fluregion.ReadRegionASCII(source_prefix,part_index,niter,t,unsteady);
    fluregion.WriteFluNative(targ_prefix);
    fluregion.WriteFluGridASCII(targ_prefix,t,unsteady);
    fluregion.WriteFluSolnASCII(targ_prefix,niter,t,unsteady);
    fluregion.CreateRegionMapFile(targ_prefix,npart,npart);
    //#ifdef _ROCSTAR_X_
    fluregion.InitRoccomWindows("flu",fluregion._borders,fluregion._patches);
    fluregion.RegisterFluSurfaceMesh();
    fluregion.RegisterVolumeSoln(true);
    fluregion.RegisterSurfaceSoln(false);
    fluregion.InitDone();
    string::size_type y = targ_prefix.find_last_of("/");
    string targ_dir = targ_prefix.substr(0,y);
    fluregion.WriteRocstar(targ_dir);
    fluregion.DestroyWindows();
    //#endif

    n++;
  }
  //#ifdef _TRAIL_MPI_
  MPI_Finalize();
  //#endif

    //#ifdef _ROCSTAR_X_
  COM_finalize();
  //#endif
  return(0);
}







