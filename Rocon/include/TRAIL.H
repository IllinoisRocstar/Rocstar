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
#ifndef _TRAIL_H_
#define _TRAIL_H_
#include "PMesh.H"
#include "BSMesh.H"
#include "GEM.H"
#include <algorithm>
#include "com.h"

int TRAIL_FD2FE_WinCreate(const std::string &wname,const std::string &outwname,std::ostream *ouf=NULL);
int TRAIL_FD2FE_WinCreate2(const std::string &wname,const std::string &outwname,std::ostream *ouf=NULL);
int TRAIL_FD2FE_Transfer(const std::string &,const std::string &,const std::string &,std::ostream *ouf=NULL);
int TRAIL_FE2FD_Transfer(const std::string &fewin,const std::string &fdwin,
			const std::string &attlist,MPI_Comm communicator,std::ostream *ouf=NULL);
int TRAIL_SurfaceMesh2Window(const std::string &wname,int pane_id,Mesh::NodalCoordinates &,Mesh::Connectivity &);

void TRAIL_AutoSurfer(const std::string &src, const std::string &trg,
		     const std::string &srcpath=".",
		     const std::string &trgpath=".",
		     const std::string &destpath=".", double t=0,
		     MPI_Comm comm=MPI_COMM_NULL,
		     std::ostream * = NULL);
void
TRAIL_CreateRobustFC(const std::string &wname,const std::string &path);

void TRAIL_MergeRocinFiles(const std::string &srcname,
			  const std::string &trgname,
			  const std::string &path="./",
			  double t=0,unsigned int np=1,
			  std::ostream * = NULL);


void TRAIL_HDF2Window( const std::string &fname, const std::string &wname,int verb=0);
void TRAIL_GetRocstarDumpStrings(const std::string &filename,std::string &wname,
				std::string &timestring,std::string &rankstring);


void
TRAIL_File2Window( const std::string &fname,
		  const std::string &wname,
		  std::vector<int>  &bcflags,
		  MPI_Comm     comm=MPI_COMM_NULL,
		  bool         apply_disp=false,
		  bool         all=false,
		  bool         with_ghost=false);
bool
TRAIL_WriteWindow(const std::string &wname,const std::string &path,
		 const std::string &twin,const std::string &tpath,
		 double t,unsigned int id,MPI_Comm comm,std::ostream * = NULL);

void
TRAIL_WriteRocinControl(std::vector<int> &pane_id,const std::string &pre,
		       int rank);
std::string TRAIL_TimeString( double t);
double TRAIL_TimeString( const std::string &);

//std::string TRAIL_CWD();

bool TRAIL_TransferSurfDataFILE(const std::string &src,
			       const std::string &trg,
			       const std::string &dest,
			       const std::string &srcpath,
			       const std::string &trgpath,
			       const std::string &destpath,
			       const std::string &crpath,
			       double t,unsigned int id,
			       MPI_Comm comm,
			       std::ostream * = NULL);

//std::string TRAIL_ResolveLink(const std::string &path);
void TRAIL_Debug(GEM_Partition &gp);
// class TRAIL_Directory : public std::vector<std::string>
// {
// protected:
//   std::string _path;
//   bool   _good;
//   DIR    *_dir;
// public:
//   TRAIL_Directory(const std::string &s = "");
//   ~TRAIL_Directory();
//   int open(const std::string &s = "");
//   void close();
//   operator void* ();
//   bool operator ! ();
// };

void TRAIL_RocmopSmooth(GEM_Partition &gp, unsigned int niter);
void
TRAIL_RocpropSmoothSurf(double *nc,unsigned int nnodes,
		       unsigned int *ec,unsigned int nel,
		       unsigned int *cnstr_type,
		       unsigned int niter);

class SolnMetaData {
public:
  std::string name;
  char loc;
  std::string unit;
  unsigned int ncomp;
};

int TRAIL_UnstructuredMesh2Pane( const std::string &wname,int pane_id,
				Mesh::UnstructuredMesh &mesh,
				SolnMetaData &smdv,
				std::vector<std::vector<double> > &soln_data,
				int verblevel	);

void
TRAIL_GetWindowSolnMetaData(const std::string &wname,std::vector<SolnMetaData> &smdv,int verblevel=0);
void TRAIL_GetWindowSolnData(const std::string &wname,std::vector<std::vector<std::vector<double> > > &soln_data,
			    std::vector<SolnMetaData> &smdv,int verblevel = 0);
void TRAIL_Window2UnstructuredMesh(const std::string &wname,std::vector<Mesh::UnstructuredMesh> &meshes,
				  std::vector<SolnMetaData> &smdv,
				  std::vector<std::vector<std::vector<double> > > &soln_data,
				  int verb=0,bool no_ghost=false);
#endif


