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
#include <fstream>
#include <string>
#include <list>
#include <map>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cassert>

#include "TRAIL_UnixUtils.H"

#include "mpi.h"

#include "GEM.H"
#include "TRAIL.H"
#include "TRAIL_Flu.H"
#include "TRAIL_Remesh.H"

#ifdef _ROCSTAR_X_
#include "roccom.h"
#endif

COM_EXTERN_MODULE( Rocface);
COM_EXTERN_MODULE( Rocout);
COM_EXTERN_MODULE( Rocin);
COM_EXTERN_MODULE( Rocblas);
COM_EXTERN_MODULE( Rocsurf);

using namespace std;

void
TRAIL_RemeshRunDirSetup(const std::string &path,double t,MPI_Comm comm, 
		       bool shuffle)
{
  int rank = 0;
  MPI_Comm_rank(comm,&rank);

  std::string rocstarrundir(TRAIL_CWD());
  std::string rocremdir("Rocrem");

  // For now, only these two are needed
  std::string solver(path);
  std::string rocman("Rocman");
  std::string ddrocm("../"+rocman);
  std::string rocburnapn("RocburnAPN");
  std::string ddsolver("../"+solver);
  std::string ddrocb("../"+rocburnapn);
  struct stat fstat;
  if(!rank){
    if(stat(solver.c_str(),&fstat)){
      std::cerr << "TRAIL_RemeshRunDirSetup: ERROR: Solver directory, " 
		<< solver << ", does not exist. Exiting." << std::endl;
      exit(1);
    }
    if(stat(rocremdir.c_str(),&fstat))
      mkdir(rocremdir.c_str(),S_IRWXU | S_IRWXG);
  }
  chdir(rocremdir.c_str());
  if(!rank){
    unlink(solver.c_str());
    symlink(ddsolver.c_str(),solver.c_str());
    if(!stat(ddrocm.c_str(),&fstat)){
      unlink(rocman.c_str());
      symlink(ddrocm.c_str(),rocman.c_str());
    }
    if(!stat(ddrocb.c_str(),&fstat)){
      unlink(rocburnapn.c_str());
      symlink(ddrocb.c_str(),rocburnapn.c_str());
    }
    unlink("rocstardir");
    symlink(rocstarrundir.c_str(),"rocstardir");
    //    if(shuffle)
    //      TRAIL_RemeshShuffle(solver,t);
  }
}

void
TRAIL_RemeshShutdown(MPI_Comm comm)
{
  chdir("rocstardir");
}

bool 
TRAIL_RemeshInitFluSurfData(GEM_Partition &gp)
{
  unsigned int ndb = gp._db.size();
  unsigned int n = 0;
  if(!TRAIL_FluInitSurfSoln(gp))
    return false;
  while(n < ndb)
    if(!TRAIL_RemeshInitFluSurfData(gp._db[n++],0))
      return false;
  return true;
}

bool
TRAIL_RemeshInitFluSurfData(GEM_DomainBoundary &db,int src_index)
{
  if(db._debug && db._out)
    *db._out << "TRAIL_RemeshInitFluSurfData:Enter" << std::endl;
  unsigned int ntri  = db._triconn.size()/3;
  unsigned int nquad = db._quadconn.size()/4;
  unsigned int nnodes = db.NNodes();
  unsigned int n_elem = ntri + nquad;
  if(db._debug && db._out)
    *db._out << "TRAIL_RemeshInitFluSurfData: Ntri = " << ntri 
	     << " Nquad = " << nquad << endl
	     << "TRAIL_RemeshInitFluSurfData: Solution sizes(" 
	     << nnodes << "," << n_elem << ")" << endl;
  unsigned int ind = 0;
  unsigned int cdind = 0;
  unsigned int vind = 0;
  // Deallocate unused datum
  db._data._field_data[0].resize(0);  // gsp's are 0 after remeshing
  db._data._stride_field[0] = 0;
  db._data._field_data[1].resize(0);  // rhofvf_alp not used
  db._data._stride_field[1] = 0;
  db._data._field_data[2].resize(0);  // nf_alp not used
  db._data._stride_field[2] = 0;
  db._data._field_data[3].resize(0);  // rhof_alp not used
  db._data._stride_field[3] = 0;
  db._data._field_data[4].resize(0);  // pf
  db._data._stride_field[4] = 0;
  db._data._field_data[5].resize(0);  // qc
  db._data._stride_field[5] = 0;
  db._data._field_data[6].resize(0);  // qr
  db._data._stride_field[6] = 0;
  db._data._field_data[7].resize(0);  // tf
  db._data._stride_field[7] = 0;
  db._data._field_data[8].resize(0);  // Tb_alp
  db._data._stride_field[8] = 0;
  db._data._field_data[11].resize(0); // Tf
  db._data._stride_field[11] = 0;
  db._data._field_data[12].resize(0); // vm
  db._data._stride_field[12] = 0;
  db._data._field_data[13].resize(0); // vs
  db._data._stride_field[13] = 0;
  db._data._field_data[14].resize(0); // vs_old
  db._data._stride_field[14] = 0;
  db._data._field_data[16].resize(0); // mdot
  db._data._stride_field[16] = 0;
  db._data._field_data[17].resize(0); // mdot_old
  db._data._stride_field[17] = 0;
  db._data._field_data[18].resize(0); // du_alp
  db._data._stride_field[18] = 0;
  db._data._field_data[20].resize(0); // sq_dist
  db._data._stride_field[20] = 0;
  if (db._solver_data._field_data.size() == 0){
    db._data._field_data[9].resize(0);
    db._data._field_data[10].resize(0);
    db._data._field_data[15].resize(0);
    db._data._int_data[0].resize(0);
    db._data._field_data[19].resize(0);
  } 
  else if((db._solver_data._stride_field[0] == 2 || 
	   db._solver_data._stride_field[0] == 3)) {
    assert(db._solver_data._stride_field[1] == 0 ||
	   db._solver_data._stride_field[1] == 3);
    db._data._int_data[0].resize(n_elem,0);
    while(ind < n_elem){
      // Keep mdot, Tflm, bflag, and ts
      db._data._field_data[9][ind]      = 
	db._solver_data._field_data[src_index][cdind++]; // mdot_alp
      db._data._field_data[10][ind]     = 
	db._solver_data._field_data[src_index][cdind++]; // Tflm_alp
      if(db._data._field_data[9][0] > 0) // if faces on fire, set bflag
	db._data._int_data[0][ind] = 1;
      if(db._solver_data._stride_field[0] == 3) 
	db._data._field_data[15][ind]            = 
	  db._solver_data._field_data[src_index][cdind++]; // ts
      //    if(db_debug && db_out)
      //      *db_out << "Element[" << ind+1 << "] rhof_alp(" 
      //	    << _soln._rhof_alp[ind] << ")  _nf_alp<"
      //	    << _soln._nf_alp[vind] << "," 
      //	    << _soln._nf_alp[vind+1] << "," 
      //	    << _soln._nf_alp[vind+2] << ">" << endl;
      ind++;
      vind+=3;
    }
    db._solver_data._field_data[src_index].resize(0);
    db._solver_data._stride_field[src_index] = 0;
    if(db._data._field_data[9][0] == -9999){
      db._data._field_data[9].resize(0);
      db._data._stride_field[9] = 0;
      db._data._field_data[10].resize(0);
      db._data._stride_field[10] = 0;
      db._data._int_data[0].resize(0);
      db._data._stride_int[0] = 0;
    }
    if(db._data._field_data[15][0] == -9999){
      db._data._field_data[15].resize(0);
      db._data._stride_field[15] = 0;
    }
    src_index++;
    ind = 0;
    cdind = 0;
    if(db._solver_data._stride_field[1] == 3){
      while(ind < nnodes){
	db._data._field_data[19][ind*3]    = 
	  db._solver_data._field_data[src_index][cdind++]; // nc_t0-x
	db._data._field_data[19][ind*3+1]  = 
	  db._solver_data._field_data[src_index][cdind++]; // nc_t0-y
	db._data._field_data[19][ind*3+2]  = 
	  db._solver_data._field_data[src_index][cdind++]; // nc_t0-z
	ind++;
      }
    }
    db._solver_data._field_data[src_index].resize(0);
    db._solver_data._stride_field[src_index] = 0;
    if(db._data._field_data[19][0] == -9999){
      db._data._field_data[19].resize(0);
      db._data._stride_field[19] = 0;
    }
  }
  else
    if(db._debug && db._out)
      *db._out << "TRAIL_RemeshInitFluSurfData: No data to process." << std::endl;
  if(db._debug && db._out)
    *db._out << "TRAIL_RemeshInitFluSurfData:Exit" << std::endl;
  return(true);
}

bool
TRAIL_RemeshInitFluVolData(GEM_Partition &gp,int src_index)
{
  //  assert((gp._solver_data._stride_field[0] == 8) &&
  //	 (gp._solver_data._stride_field[1] == 3 || 
  //	  gp._solver_data._stride_field[1] == 0));
  unsigned int ntet   = gp._tetconn.size()/4;
  unsigned int nhex   = gp._hexconn.size()/8;
  unsigned int npris  = gp._prisconn.size()/6;
  unsigned int npyr   = gp._pyrconn.size()/5;
  unsigned int nnodes = gp._nc.size()/3;
  unsigned int ncells = ntet + nhex + npris + npyr;
  //  assert(gp._solver_data._field_data[src_index].size() == 8*ncells);
  //  assert(gp._solver_data._field_data[src_index+1].size() == 3*nnodes);
  if(!TRAIL_FluInitVolSoln(gp))
    return false;
  unsigned int ind = 0;
  unsigned int cdind = 0;
  unsigned int vind = 0;
  while(ind < ncells){
    // rhof
    gp._data._field_data[2][ind] = 
      gp._solver_data._field_data[src_index][cdind++]; 
    // rhovf_x
    gp._data._field_data[3][vind++] = 
      gp._solver_data._field_data[src_index][cdind++]; 
    // rhovf_y
    gp._data._field_data[3][vind++] = 
      gp._solver_data._field_data[src_index][cdind++]; 
    // rhovf_z
    gp._data._field_data[3][vind++] = 
      gp._solver_data._field_data[src_index][cdind++];
    // rhoEf 
    gp._data._field_data[4][ind] = 
      gp._solver_data._field_data[src_index][cdind++]; 
    // pf
    gp._data._field_data[5][ind] = 
      gp._solver_data._field_data[src_index][cdind++]; 
    // Tf
    gp._data._field_data[6][ind] = 
      gp._solver_data._field_data[src_index][cdind++]; 
    // af
    gp._data._field_data[7][ind] = 
      gp._solver_data._field_data[src_index][cdind++]; 
    ind++;
  }
  ind = 0;
  gp._solver_data._field_data[src_index].resize(0);
  gp._solver_data._stride_field[src_index] = 0;
  src_index++;
  // nodal displacements (disp)
  if(gp._solver_data._field_data[src_index].size() == 3*nnodes){
    while(ind < nnodes){
      gp._data._field_data[8][ind*3]   = 
	gp._solver_data._field_data[src_index][ind*3];
      gp._data._field_data[8][ind*3+1] = 
	gp._solver_data._field_data[src_index][ind*3+1];
      gp._data._field_data[8][ind*3+2] = 
	gp._solver_data._field_data[src_index][ind*3+2];
      ind++;
    }
  }
  else{
    while(ind < nnodes){
      gp._data._field_data[8][ind*3]   = 0.0;
      gp._data._field_data[8][ind*3+1] = 0.0;
      gp._data._field_data[8][ind*3+2] = 0.0;
      ind++;
    }
  }
  // Grid speeds are 0 after remeshing
  gp._data._field_data[9].resize(0);
  gp._data._stride_field[9] = 0;
  gp._solver_data._field_data[src_index].resize(0);
  gp._solver_data._stride_field[src_index] = 0;
  return(true);
}


void
TRAIL_RemeshFixRocstarFiles(GEM_Partition &gp,const std::string &path,double t)
{
  // Read in all the files and output to a single Rocin control file
  // Volume first:
  // set up base filename
  TRAIL_MergeRocinFiles(gp.volume_window,gp.volume_window,path,t,gp._npart);
  TRAIL_MergeRocinFiles(gp.surface_window,gp.surface_window,path,t,gp._npart);
}

bool
TRAIL_RemeshAutoSurfer(GEM_Partition &gp,
		      const std::string  &src,
		      const std::string  &srcpath,
		      const std::string  &trgpath,
		      const std::string  &destpath,
		      double        t,
		      MPI_Comm      comm)
{
  std::string timestring(TRAIL_TimeString(t));
  std::string srcfile(src + "_in_" + timestring + ".txt");
  std::string trgfile(gp.surface_window+"_in_"+timestring+".txt");
  std::string srcwin(src);
  std::string trailwin(src+"_trail");
  std::string trgwin(gp.surface_window);
  std::string homedir(TRAIL_CWD());
  std::string format("HDF");
  std::string crpath(trgpath+"/AutoSurf");
  struct stat fstat;
  int rank = 0;
  MPI_Comm_rank(comm,&rank);

  if(!rank){
    if(!stat(crpath.c_str(),&fstat)){
      rename(crpath.c_str(),(crpath+"_save").c_str());
      std::cerr << "TRAIL_RemeshAutoSurfer: WARNING: " << crpath << " already "
		<< "existed.  Renaming to " << crpath << "_save." << std::endl;
    }
    TRAIL_CreateDirectory(crpath);
  }
  MPI_Barrier(comm);
  // Serial step to create the common refinement
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RemeshAutoSurfer: Creating overlay for coupled surfaces."
	     << std::endl;

  if(!rank)
    TRAIL_AutoSurfer(srcwin,trgwin,srcpath,trgpath,crpath,t,comm,gp._out);
  MPI_Barrier(comm);
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RemeshAutoSurfer: Overlay for coupled surfaces complete."
	     << std::endl
	     << "TRAIL_RemeshAutoSurfer: Transferring data....";
  MPI_Barrier(comm);
  if(!rank && gp._debug)
    std::cout << "Roctrail> Transferring data in parallel..." << std::endl;
  
  
  TRAIL_TransferSurfDataFILE(srcwin,trgwin,"new_surf",srcpath,trgpath,
			    destpath,crpath,t,gp._id,comm,gp._out);
  
  
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_AutoSurfer: Transfer for coupled surfaces complete." 
	     << std::endl
	     << "TRAIL_AutoSurfer: Writing populated surfaces." 
	     << std::endl;
  if(!rank && gp._debug)
    std::cout << "Roctrail> Transfer done." << std::endl;

  
  MPI_Barrier(comm);
  return(true);
}

// SERIAL- DO NOT INVOKE WITH MULTIPROCS
void
TRAIL_RemeshShuffle(const std::string &solver,double t,bool debug)
{
  if(t==0.0)
    return;
  std::string timestring(TRAIL_TimeString(t));
  std::string time0string(TRAIL_TimeString(0.0));
  std::string ssd("Modin");
  std::string rsd("Rocout");
  std::string strg("Modin.remesh_"+timestring);
  std::string strg0("Modin.remesh_"+time0string);
  std::string rtrg("Rocout.remesh_"+timestring);
  std::string rtrg0("Rocout.remesh_"+time0string);
  std::string solversrc(solver+"/"+ssd);
  std::string solvertrg0(solver+"/"+strg0);
  std::string solvertrg(solver+"/"+strg);
  std::string rocstarsrc(solver+"/"+rsd);
  std::string rocstartrg(solver+"/"+rtrg);
  std::string rocstartrg0(solver+"/"+rtrg0);
  std::string rocstarlast(solver+"/Rocout.remesh_last");
  std::string solverlast(solver+"/Modin.remesh_last");
  struct stat fstat0;
  struct stat fstat1;
  if(lstat(rocstarsrc.c_str(),&fstat0)){
    std::cerr << "TRAIL_Remesh::DirectoryShuffling: ERROR: " << rocstarsrc 
	      << " does not exist.  Exiting." << std::endl;
    exit(1);
  }
  // If the <solver>/Rocout is a directory, then this must be the first remesh
  // The plan is to mv <solver>/Rocout to <solver>/Rocout.remesh_00.000000
  if(S_ISDIR(fstat0.st_mode)){
    if(debug)
      std::cout << "TRAIL_RemeshShuffle: " << rocstarsrc << " directory exists."
		<< std::endl;
    // bummer, <solver>/Rocout.remesh_00.000000 already existed.  this should 
    // probably never happen, but if it does, warn and back it up as  
    // <solver>/Rocout.remesh_save_00.000000
    if(!stat(rocstartrg0.c_str(),&fstat1)){
      std::string savepath(rocstarsrc+".remesh_save_"+time0string);
      rename(rocstartrg0.c_str(),savepath.c_str());
      std::cerr << "TRAIL_Remesh::DirectoryShuffling: WARNING: " << rocstartrg0
		<< " already existed.  Saving it as " << savepath << "." 
		<< std::endl;
    }
    rename(rocstarsrc.c_str(),rocstartrg0.c_str());
    symlink(rtrg0.c_str(),rocstarsrc.c_str());
  }
  // <solver>/Rocout is not a directory, it is a link, that's okay unless it's
  // already pointing at <solver>/Rocout.remesh_<timestamp>.  In this case, the
  // most likely scenario is the one where you are remeshing the same timestamp
  // multiple times without running rocstar.
  else{
    std::string pointsto(ResolveLink(rocstarsrc.c_str()));
    if(debug)
      std::cout << "TRAIL_RemeshShuffle: " << rocstarsrc 
		<< " is a link to " << pointsto << "." << std::endl;
    if(pointsto == rtrg){
      if(debug)
	std::cout << "TRAIL_RemeshShuffle: " << rocstarsrc 
		  << " directory did not exist." << std::endl;
      if(stat(rocstarlast.c_str(),&fstat1) ||
	 stat(solverlast.c_str(),&fstat1)){
	std::cerr << "TRAIL_Remesh::DirectoryShuffling: ERROR: Detected "
		  << "multiple remesh of same dump, but last source cannot "
		  << "be found. Exiting." << std::endl;
	exit(1);
      }
      unlink(rocstarsrc.c_str());
      rename(rocstarlast.c_str(),rocstarsrc.c_str());
      unlink(solversrc.c_str());
      rename(solverlast.c_str(),solversrc.c_str());
    }
    else{ // it points to the last remeshing
      unlink(rocstarlast.c_str());
      unlink(solverlast.c_str());
    }
  }
  // If rocstartrg already exists, warn and back it up 
  if(!stat(rocstartrg.c_str(),&fstat0)){
    std::string savepath(rocstarsrc+".remesh_save_"+timestring);
    rename(rocstartrg.c_str(),savepath.c_str());
    std::cerr << "TRAIL_Remesh::DirectoryShuffling: WARNING: " << rocstartrg
	      << " already existed.  Saving it as " << savepath << "." 
	      << std::endl;
  }
  mkdir(rocstartrg.c_str(),S_IRWXG | S_IRWXU);
  symlink(ResolveLink(rocstarsrc.c_str()).c_str(),rocstarlast.c_str());
  if(lstat(solversrc.c_str(),&fstat0)){
    std::cerr << "TRAIL_Remesh::DirectoryShuffling: ERROR: " << solversrc 
	      << " does not exist.  Exiting." << std::endl;
    exit(1);
  }
  // If the <solver>/Modin is a directory, then this must be the first remesh
  // The plan is to mv <solver>/Modin to <solver>/Modin.remesh_00.000000
  if(S_ISDIR(fstat0.st_mode)){
    // bummer, <solver>/Modin.remesh_00.000000 already existed.  this should 
    // probably never happen, but if it does, warn and back it up as  
    // <solver>/Modin.remesh_save_00.000000
    if(!stat(solvertrg0.c_str(),&fstat1)){
      std::string savepath(solversrc+".remesh_save_"+time0string);
      rename(solvertrg0.c_str(),savepath.c_str());
      std::cerr << "TRAIL_Remesh::DirectoryShuffling: WARNING: " << solvertrg0
		<< " already existed.  Saving it as " << savepath << "." 
		<< std::endl;
    }
    rename(solversrc.c_str(),solvertrg0.c_str());
    symlink(strg0.c_str(),solversrc.c_str());
  }
  // If solvertrg already exists, warn and back it up 
  if(!stat(solvertrg.c_str(),&fstat0)){
    std::string savepath(solversrc+".remesh_save_"+timestring);
    rename(solvertrg.c_str(),savepath.c_str());
    std::cerr << "TRAIL_Remesh::DirectoryShuffling: WARNING: " << solvertrg
	      << " already existed.  Saving it as " << savepath << "." 
	      << std::endl;
  }
  mkdir(solvertrg.c_str(),S_IRWXG | S_IRWXU);
  symlink(ResolveLink(solversrc.c_str()).c_str(),solverlast.c_str());
  // Deal with RocburnAPN
  std::string rocburnapn("RocburnAPN");
  std::string rocburnapnsrc(rocburnapn+"/"+rsd);
  std::string rocburnapntrg0(rocburnapn+"/"+rtrg0);
  std::string rocburnapntrg(rocburnapn+"/"+rtrg);
  if(!lstat(rocburnapnsrc.c_str(),&fstat0)){
    if(S_ISDIR(fstat0.st_mode)){
      if(!stat(rocburnapntrg0.c_str(),&fstat1)){
	// RocburnAPN/Rocout is a directory, but RocburnAPN/Rocout.remesh_0 
	// already exists.  Warn and back up the previous one.
	std::string savepath(rocburnapnsrc+".remesh_save_"+time0string);
	rename(rocburnapntrg0.c_str(),savepath.c_str());
	std::cerr << "TRAIL_Remesh::DirectoryShuffling: WARNING: " 
		  << rocburnapntrg0 << " already existed, backup it up as "
		  << savepath << "." << std::endl;
      }
      rename(rocburnapnsrc.c_str(),rocburnapntrg0.c_str());
      symlink(rtrg0.c_str(),rocburnapnsrc.c_str());
    }
    if(!stat(rocburnapntrg.c_str(),&fstat1)){
      // RocburnAPN/Rocout.remesh_<timestamp> already exists, warn and back it
      // up.
      std::string savepath(rocburnapnsrc+".remesh_save_"+timestring);
      rename(rocburnapntrg.c_str(),savepath.c_str());
      std::cerr << "TRAIL_Remesh::DirectoryShuffling: WARNING: " 
		<< rocburnapntrg << " already existed, backup it up as "
		<< savepath << "." << std::endl;
    }
    mkdir(rocburnapntrg.c_str(), S_IRWXG | S_IRWXU);
    unlink(rocburnapnsrc.c_str());
    symlink(rtrg.c_str(),rocburnapnsrc.c_str());
  }
}

// bool
// TRAIL_RemeshWrite(GEM_Partition &gp,const string &path,double t,MPI_Comm comm,
// 		 bool transfer_surface_data)
// {
//   transfer_surface_data = false;
//   std::string homedir(CWD());
//   std::string timestring(TRAIL_TimeString(t));
//   std::string solver(path);
//   std::string solversrc(solver+"/Modin");
//   std::string strg("Modin.remesh_"+timestring);
//   std::string solvertrg(solver+"/"+strg);
//   std::string rocstarsrc(solver+"/Rocout");
//   std::string rtrg("Rocout.remesh_"+timestring);
//   std::string rocstartrg(solver+"/"+rtrg);
//   std::string volcntl("remesh_vol_in_"+timestring+".txt");
//   std::string surfcntl("remesh_surf_in_"+timestring+".txt");
//   std::string fluvolcntl(rocstarsrc+"/fluid_in_"+timestring+".txt");
//   std::string flusurfcntl(rocstarsrc+"/ifluid_in_"+timestring+".txt");
//   COM_set_profiling(0);
//   int rank = 0;
//   MPI_Comm_rank(comm,&rank);
//   int nproc = 0;
//   MPI_Comm_size(comm,&nproc);
//   gp._npart = nproc;
//   //  if(!rank)
//   //    DirectoryShuffling(path,t);
//   //  transfer_surface_data = true;
//   if(comm != MPI_COMM_WORLD && gp._out)
//     *gp._out << "TRAIL_RemeshWrite: WARNING: communicator is "
// 	     << (comm == MPI_COMM_NULL ? "null." : "not worldly.") 
// 	     << std::endl;
//   if(gp._debug && gp._out)
//     *gp._out << "TRAIL_RemeshWrite: Preparing to write remeshed dataset for " 
// 	     << path << "." << std::endl;
//   MPI_Barrier(comm);
//   if(!rank)
//     std::cout << "Initializing solver data for " << path << "...";
//   if(gp._debug){
//     if(gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Window test..." << std::endl;
//     MPI_Barrier(comm);
//     string test("test");
//     if(gp._out)
//       *gp._out << "   creation..." << std::endl;
//     MPI_Barrier(comm);
//     COM_new_window(test);
//     MPI_Barrier(comm);
//     if(gp._out)
//       *gp._out << "   destruction..." << std::endl;
//     MPI_Barrier(comm);
//     COM_delete_window(test);
//     MPI_Barrier(comm);
//     if(gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Window test passed, graduating..." 
// 	       << std::endl;
//   }
//   std::string pre(path+"/Rocout");
//   // Write Rocflu remesh
//   if(path == "Rocflu"){
    
//     // COM_delete_window("ifluid");
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Initializing solver data..." << std::endl;
//     if(!TRAIL_FluInitSolver(gp,path))
//       return(false);
//     // Now that we have Rocflu's paths defined and target directories created,
//     // we can copy in the necessary Rocflu casefiles from the solver source to
//     // the solver target directories.
//     TRAIL_FluCopyCaseFiles(gp,solvertrg);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Solver data init done." << std::endl
// 	       << "TRAIL_RemeshWrite: Writing native files (pass 1)..." 
// 	       << std::endl;
//     if(!TRAIL_FluWriteNative(gp,solvertrg))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Solver native files written (pass 1)." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Populating remote border indices..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!TRAIL_FluPopRemBordIndFILE(gp,0.0,true,solvertrg))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Remote border indices populated." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Writing native files (pass 2)..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!TRAIL_FluWriteNative(gp,solvertrg))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Solver native files written (pass 2)." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Initializing volume data..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!TRAIL_RemeshInitFluVolData(gp,0))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Volume data initialized." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Initializing surface data..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!TRAIL_RemeshInitFluSurfData(gp))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Surface data initialized." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Initializing Roccom Windows..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!rank)
//       std::cout << "done." << std::endl
// 		<< "Creating Roccom Windows...";
//     if(!gp.InitRoccomWindows("remesh"))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Roccom windows initialized." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Registering surface meshes..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!TRAIL_FluRegisterSurfMesh(gp))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Surface meshes registered." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Registering volume solution..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!TRAIL_FluRegisterVolSoln(gp))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Volume soln registered." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Registering surface solution..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!TRAIL_FluRegisterSurfSoln(gp))
//       return(false);
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Surface soln registered." 
// 	       << std::endl
// 	       << "TRAIL_RemeshWrite: Preparing to write..."
// 	       << std::endl;
//     MPI_Barrier(comm);
//     if(!gp.WindowInitDone())
//       return(false);
//     if(!rank)
//       std::cout << "done." << std::endl
// 		<< "Writing windows to Rocstar format..." 
// 		<< std::endl;
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << " writing..." << std::endl;
//     MPI_Barrier(comm);
//     if(!gp.WriteRocstar(rocstartrg,t))
//       return(false);
//     MPI_Barrier(comm);
//     gp.DestroyWindows();
//     // Take care of the directory thing
//     // Record the name of the real source directory for surface xfer
//     // rm <solver>/Rocout link and relink to rocstartrg
//     // rm <solver>/Modin link and relink to solvertrg
//     if(!rank){
//       //      rename(rocstarsrc.c_str(),(solver+"/Rocout_old").c_str());
//       unlink(rocstarsrc.c_str());
//       symlink(rtrg.c_str(),rocstarsrc.c_str());
//       unlink(solversrc.c_str());
//       symlink(strg.c_str(),solversrc.c_str());
//     };
//     if(gp._debug && gp._out)
//       *gp._out << " done writing" << std::endl;
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Merging Rocin files..." << std::endl;
//     if(!rank){
//       TRAIL_MergeRocinFiles(gp.volume_window,gp.volume_window,rocstartrg,t,gp._npart);
//       TRAIL_MergeRocinFiles(gp.surface_window,gp.surface_window,rocstartrg,t,gp._npart);
//     }
//       TRAIL_RemeshFixRocstarFiles(gp,rocstartrg,t);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: Rocin files merged." << std::endl;
//     MPI_Barrier(comm);
//     if(!rank){
//       symlink(volcntl.c_str(),fluvolcntl.c_str());
//       symlink(surfcntl.c_str(),flusurfcntl.c_str());
//     }
//     if(transfer_surface_data){
//       if(!rank)
// 	std::cout << "Proceeding with surface data transfer..." << std::endl;
//       if(gp._debug && gp._out)
// 	*gp._out << "TRAIL_RemeshWrite: Testing surface transfer..." 
// 		 << std::endl;
//       TRAIL_RemeshAutoSurfer(gp,"ifluid",(solver+"/Rocout.remesh_last"),
// 			    rocstartrg,rocstartrg,t,comm);
//       MPI_Barrier(comm);
//       if(!rank){
// 	unlink(flusurfcntl.c_str());
// 	surfcntl.assign("new_surf_in_"+timestring+".txt");
// 	symlink(surfcntl.c_str(),flusurfcntl.c_str());
//       }
//       if(gp._debug && gp._out)
// 	*gp._out << "TRAIL_RemeshWrite: Surface transfer complete." 
// 		 << std::endl;
//       MPI_Barrier(comm);
      
//       if(!rank)
// 	std::cout << "Surface data transfer succeeded." << std::endl;
//     }
//     //    if(!rank)
//     //      unlink((solver+"/Rocout_old").c_str());
//     MPI_Barrier(comm);
//     if(gp._debug && gp._out)
//       *gp._out << "TRAIL_RemeshWrite: All done." << std::endl
// 	       << "TRAIL_RemeshWrite: Exit" << std::endl;    
    
//     MPI_Barrier(comm);
//     if(!rank)
//       std::cout << "Rocstar restart data done." << std::endl;
//     MPI_Barrier(comm);
//     return(true);
//   }
//   if(gp._out)				      
//     *gp._out << "TRAIL_RemeshWrite::Error: Unknown solver, " << path << "." 
// 	     << std::endl
// 	     << (gp._debug ? "TRAIL_RemeshWrite: Exit" : "") << std::endl;
//   return(false);
// }
bool
TRAIL_RemeshWrite(GEM_Partition &gp,
		 const std::string &solver,
		 const std::string &srcd,
		 const std::string &trgd,
		 double t,MPI_Comm comm,
		 bool surftrans,bool rocstarshuffle)
{
  //  transfer_surface_data = false;
  std::string source_directory(srcd);
  std::string target_directory(trgd);
  std::string rocstarrunsrc("Rocout");
  std::string defaultrocstartrg("Rocout.remesh_"+TRAIL_TimeString(t));
  if(source_directory.empty())
    source_directory = rocstarrunsrc;
  if(target_directory.empty())
    target_directory = defaultrocstartrg;
  std::string homedir(TRAIL_CWD());
  std::string timestring(TRAIL_TimeString(t));
  std::string solversrc(solver+"/Modin");
  std::string strg("Modin.remesh_"+timestring);
  std::string solvertrg(solver+"/"+strg);
  std::string rocstarsrc(solver+"/"+source_directory);
  std::string rtrg(target_directory);
  std::string rocstartrg(solver+"/"+rtrg);
  std::string volcntl("remesh_vol_in_"+timestring+".txt");
  std::string surfcntl("remesh_surf_in_"+timestring+".txt");
  std::string fluvolcntl(rocstartrg+"/fluid_in_"+timestring+".txt");
  std::string flusurfcntl(rocstartrg+"/ifluid_in_"+timestring+".txt");
  
  COM_set_profiling(0);
  int rank = 0;
  MPI_Comm_rank(comm,&rank);
  int nproc = 0;
  MPI_Comm_size(comm,&nproc);
  gp._npart = nproc;
  std::cout << flush;
  MPI_Barrier(comm);
  std::string rocstardir("..");
  if(!rank){
    TRAIL_CD(rocstardir,gp._out);
    rocstardir.assign(TRAIL_CWD());
    TRAIL_CD(homedir,gp._out);
  } 
  std::string solverdir(rocstardir+"/"+solver);
  if(!TRAIL_FILEEXISTS(rocstarsrc)){
    if(!rank)
      std::cerr << "TRAIL_RemeshWrite: Fatal Error: Source directory, " 
		<< rocstarsrc << " does not exist.  Exiting." << std::endl;
    if(gp._out)
      *gp._out << "TRAIL_RemeshWrite: Fatal Error: Souce directory, " 
	       << rocstarsrc << " does not exist.  Exiting." << std::endl;
    exit(1);
  } 
  if(!TRAIL_FILEEXISTS(solversrc)){
    if(!rank)
      std::cerr << "TRAIL_RemeshWrite: Fatal Error: Source directory, " 
		<< solversrc << " does not exist.  Exiting." << std::endl;
    if(gp._out)
      *gp._out << "TRAIL_RemeshWrite: Fatal Error: Souce directory, " 
	       << solversrc << " does not exist.  Exiting." << std::endl;
    exit(1);
  }   
  if(!TRAIL_FILEEXISTS(rocstartrg))
    TRAIL_CreateDirectory(rocstartrg);
  if(!TRAIL_FILEEXISTS(solvertrg))
    TRAIL_CreateDirectory(solvertrg);
  if(comm != MPI_COMM_WORLD && gp._out)
    *gp._out << "TRAIL_RemeshWrite: WARNING: communicator is "
	     << (comm == MPI_COMM_NULL ? "null." : "not worldly.") 
	     << std::endl;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RemeshWrite: Preparing to write remeshed dataset for " 
	     << solver << "." << std::endl;
  MPI_Barrier(comm);
  if(!rank)
    std::cout << std::endl
	      << "Roctrail> Processing solver data for " << solver 
	      << std::endl;
  MPI_Barrier(comm);
  if(gp._debug){
    if(gp._out)
      *gp._out << "TRAIL_RemeshWrite: Window test..." << std::endl;
    MPI_Barrier(comm);
    string test("test");
    if(gp._out)
      *gp._out << "   creation..." << std::endl;
    MPI_Barrier(comm);
    COM_new_window(test);
    MPI_Barrier(comm);
    if(gp._out)
      *gp._out << "   destruction..." << std::endl;
    MPI_Barrier(comm);
    COM_delete_window(test);
    MPI_Barrier(comm);
    if(gp._out)
      *gp._out << "TRAIL_RemeshWrite: Window test passed, graduating..." 
	       << std::endl;
  }
  // Write Rocflu remesh
  if(solver == "Rocflu"){
    
    // COM_delete_window("ifluid");
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Initializing solver data..." << std::endl;
    if(!TRAIL_FluInitSolver(gp,solver))
      return(false);
    // Now that we have Rocflu's paths defined and target directories created,
    // we can copy in the necessary Rocflu casefiles from the solver source to
    // the solver target directories.
    TRAIL_FluCopyCaseFiles(gp,solvertrg);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Solver data init done." << std::endl
	       << "TRAIL_RemeshWrite: Writing native files (pass 1)..." 
	       << std::endl;
    if(!TRAIL_FluWriteNative(gp,solvertrg))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Solver native files written (pass 1)." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Populating remote border indices..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!TRAIL_FluPopRemBordIndFILE(gp,0.0,true,solvertrg))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Remote border indices populated." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Writing native files (pass 2)..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!TRAIL_FluWriteNative(gp,solvertrg))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Solver native files written (pass 2)." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Initializing volume data..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!TRAIL_RemeshInitFluVolData(gp,0))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Volume data initialized." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Initializing surface data..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!TRAIL_RemeshInitFluSurfData(gp))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Surface data initialized." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Initializing Roccom Windows..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!rank)
      std::cout << "Roctrail> Creating Roccom Windows" << std::endl;
    if(!gp.InitRoccomWindows("remesh"))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Roccom windows initialized." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Registering surface meshes..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!TRAIL_FluRegisterSurfMesh(gp))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Surface meshes registered." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Registering volume solution..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!TRAIL_FluRegisterVolSoln(gp))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Volume soln registered." 
	       << std::endl
	       << "TRAIL_RemeshWrite: Registering surface solution..."
	       << std::endl;
    MPI_Barrier(comm);
    if(!TRAIL_FluRegisterSurfSoln(gp))
      return(false);
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: Surface soln registered." 
	       << std::endl;
    MPI_Barrier(comm);
    if(!gp.WindowInitDone())
      return(false);
    if(!rank)
      std::cout << "Roctrail> Writing windows to Rocstar format"
		<< std::endl;
    MPI_Barrier(comm);
    if(gp._debug && gp._out)
      *gp._out << " writing..." << std::endl;
    MPI_Barrier(comm);
    if(!gp.WriteRocstar(rocstartrg,t))
      return(false);
    MPI_Barrier(comm);
    TRAIL_CD(homedir,gp._out);
    gp.DestroyWindows();
    if(gp._debug && gp._out)
      *gp._out << " done writing" << std::endl;
    MPI_Barrier(comm);

    // Set up Rocin Control Files for Rocstar
    if(rocstarshuffle){
      if(!rank){
	std::cout << "Roctrail> Setting up Rocin control files for Rocstar" 
		  << std::endl;
	if(!TRAIL_FILEEXISTS(fluvolcntl+".orig") && TRAIL_FILEEXISTS(fluvolcntl))
	  rename(fluvolcntl.c_str(),(fluvolcntl+".orig").c_str());
	if(!TRAIL_FILEEXISTS(flusurfcntl+".orig") && 
	   TRAIL_FILEEXISTS(flusurfcntl))
	  rename(flusurfcntl.c_str(),(flusurfcntl+".orig").c_str());
	symlink(volcntl.c_str(),fluvolcntl.c_str());
	symlink(surfcntl.c_str(),flusurfcntl.c_str());
      }
    }
    if(surftrans){
      if(!rank)
	std::cout << "Roctrail> Proceeding with surface data transfer..." 
		  << std::endl;
      if(gp._debug && gp._out)
	*gp._out << "TRAIL_RemeshWrite: Performing surface data transfer..." 
		 << std::endl;
      TRAIL_RemeshAutoSurfer(gp,"ifluid",rocstarsrc,
			    rocstartrg,rocstartrg,t,comm);
      MPI_Barrier(comm);
      // Reset Rocin Control Files for populated surface
      if(rocstarshuffle){
	if(!rank){
	  std::cout << "Roctrail> Resetting Rocin control files"
		    << std::endl;
	  unlink(flusurfcntl.c_str());
	  surfcntl.assign("new_surf_in_"+timestring+".txt");
	  symlink(surfcntl.c_str(),flusurfcntl.c_str());
	}
      }
      if(!rank)
	std::cout << "Roctrail> Surface data transfer done." << std::endl;
      if(gp._debug && gp._out)
	*gp._out << "TRAIL_RemeshWrite: Surface transfer complete." 
		 << std::endl;
      MPI_Barrier(comm);
    }
    MPI_Barrier(comm);
    TRAIL_CD(homedir,gp._out);
    // Prepare Rocstar restart files
    if(rocstarshuffle){
      if(!rank){
	std::cout << "Roctrail> Preparing Rocstar restart directories." 
		  << std::endl;
	std::string rocstarbackup0(rocstardir+"/"+solver+
				   "/Rocout.remesh_00.000000");
	if(rocstartrg != rocstarrunsrc){
	  // <solver>/Rocout is an actual directory, back it up as 
	  // Rocout.remesh_00.000000
	  if(TRAIL_ISDIR(solverdir+"/"+rocstarrunsrc) && 
	     !TRAIL_ISLINK(solverdir+"/"+rocstarrunsrc)){
	    TRAIL_SafeRemove(rocstarbackup0,".save");
	    rename((solverdir+"/"+rocstarrunsrc).c_str(),
		   rocstarbackup0.c_str());
	  }
	  else if(TRAIL_ISLINK(solverdir+"/"+rocstarrunsrc))
	    unlink((solverdir+"/"+rocstarrunsrc).c_str());
	  if(gp._debug)
	    std::cout << "Roctrail> Attempting to link " << rocstarrunsrc 
		      << " to " << rtrg << "." << std::endl;
	  symlink(rtrg.c_str(),(solverdir+"/"+rocstarrunsrc).c_str());
	  // Modin is a directory, act accordingly
	  std::string solverbackup0(rocstardir+"/"+solver+
				    "/Modin.remesh_00.000000");
	  if(TRAIL_ISDIR(rocstardir+"/"+solversrc) && 
	     !TRAIL_ISLINK(rocstardir+"/"+solversrc)){
	    TRAIL_SafeRemove(solverbackup0,".save");
	    rename((rocstardir+"/"+solversrc).c_str(),solverbackup0.c_str());
	  }
	  else if(TRAIL_ISLINK(solversrc))
	    unlink((rocstardir+"/"+solversrc).c_str());
	  if(gp._debug)
	    std::cout << "Roctrail> Attempting to link " << solversrc << " to "
		      << strg << "." << std::endl;
	  symlink(strg.c_str(),(rocstardir+"/"+solversrc).c_str());
	}
      }
      if(!rank)
	std::cout << "Roctrail> Rocstar restart data done." << std::endl;
    }
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RemeshWrite: All done." << std::endl
	       << "TRAIL_RemeshWrite: Exit" << std::endl;    
    if(!rank)
      std::cout << "Roctrail> Goodbye" << std::endl;
    MPI_Barrier(comm);
    return(true);
  }
  if(gp._out)				      
    *gp._out << "TRAIL_RemeshWrite::Error: Unknown solver, " << solver << "." 
	     << std::endl
	     << (gp._debug ? "TRAIL_RemeshWrite: Exit" : "") << std::endl;
  return(false);
}







