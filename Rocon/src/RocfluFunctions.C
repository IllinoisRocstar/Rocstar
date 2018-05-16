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
#include <cmath>
#include <cassert>

//#ifdef _TRAIL_MPI_
#include "mpi.h"
//#endif

#include "GEM.H"
#include "TRAIL.H"
#include "TRAIL_Flu.H"

//#ifdef _ROCSTAR_X_
#include "com.h"
//#endif

using namespace std;

istream &
SkipLines(istream &ioo,unsigned int n)
{
  string line;
  unsigned int l = 0;
  while(l++ < n)
    getline(ioo,line);
  return(ioo);
}

bool
TRAIL_FluInitSolver(GEM_Partition &gp,const string &prefix)
{
  if(!TRAIL_FluReadControlFile(gp,prefix))
    return(false);
  gp._solver_data._int_data.resize(2);
  gp._solver_data._stride_int.resize(2);
  gp._solver_data._int_data[1].resize(1);
  gp._solver_data._stride_int[1] = 1;
  gp._solver_data._int_data[1][0] = TRAIL_FluNumPatches(gp);
  if(gp._solver_data._int_data[1][0] == 0)
    return(false);
  return true;
}

bool 
TRAIL_FluReadControlFile(GEM_Partition &gp,const string &prefix)
{
  ifstream Inf;
  string fname;
  if(!prefix.empty())
    fname = prefix + "/RocfluControl.txt";
  Inf.open(fname.c_str());
  if(!Inf)
    return false;
  gp._solver_data._string_data.resize(3);
  gp._solver_data._string_data[2] = prefix; // solver name
  // casename and casepath
  Inf >> gp._solver_data._string_data[0] >> gp._solver_data._string_data[1];
  Inf.close();
  return(true);
}

bool
TRAIL_FluWriteCOM(GEM_Partition &gp,const std::string &path)
{
  if(gp._solver_data._string_data.size() < 3){
    if(gp._out)
      *gp._out << "TRAIL_FluWriteCOM::Error: No casename or path found."
	       << std::endl;
    return(false);
  }
  string pre("./" + path + "/" + gp._solver_data._string_data[0]);
  if(gp._solver_data._int_data.empty()){
    if(gp._out)
      *gp._out 
	<< "TRAIL_FluWriteCOM: Error, called with empty rbid's."
	<< endl;
    return(false);
  }
  assert(gp._solver_data._int_data[0].size() == gp._pb.size());
  ofstream Ouf;
  ostringstream Ostr;
  Ostr << pre << ".com_";
  Ostr << setw(5) << setfill('0');
  Ostr << gp._id;
  Ouf.open(Ostr.str().c_str());
  if(!Ouf){
    if(gp._out)
      *gp._out 
	<< "TRAIL_FluWriteCOM: Error: Could not open " << Ostr.str() << "."
	<< std::endl;
    return(false);
  }
  Ouf << "# ROCFLU communication lists file" << endl
      << "# Dimensions"  << endl
      << setw(8) << gp._pb.size() << endl
      << "# Information" << endl;
  int pbindex = 0;
  vector<GEM_PartitionBoundary>::iterator fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    Ouf << setw(8) << fbi->_rpart << setw(8) 
	<< gp._solver_data._int_data[0][pbindex++] 
	<< endl;
    fbi++;
  }
  Ouf << "# Cells" << endl;
  fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    Ouf << setw(8) << fbi->_sendcells.size() << setw(8) 
	<< fbi->_recvcells.size() << endl;
    int line = 0;
    vector<unsigned int>::iterator sci = fbi->_sendcells.begin();
    while(sci != fbi->_sendcells.end()){
      Ouf << setw(8) << *sci++;
      line++;
      if(line == 10){ 
	Ouf << endl;
	line = 0;
      }
    }
    if(line){
      Ouf << endl;
      line = 0;
    }
    sci = fbi->_recvcells.begin();
    while(sci != fbi->_recvcells.end()){
      Ouf << setw(8) << *sci++;
      line++;
      if(line == 10){ 
	Ouf << endl;
	line = 0;
      }
    }
    if(line){
      Ouf << endl;
      line = 0;
    }
    fbi++;
  }    
  Ouf << "# Vertices" << endl;
  fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    Ouf << setw(8) << fbi->_sendnodes.size() << setw(8) 
	<< fbi->_recvnodes.size() << setw(8) << fbi->_sharenodes.size() 
	<< endl;
    vector<unsigned int>::iterator ni = fbi->_sendnodes.begin();
    int line = 0;
    while(ni != fbi->_sendnodes.end()){
      Ouf << setw(8) << *ni++;
      line++;
      if(line == 10){
	Ouf << endl;
	line = 0;
      }
    } 
    if(line || fbi->_sendnodes.empty()){
      Ouf << endl;
      line = 0;
    }
    ni = fbi->_recvnodes.begin();
    line = 0;
    while(ni != fbi->_recvnodes.end()){
      Ouf << setw(8) << *ni++;
      line++;
      if(line == 10){
	Ouf << endl;
	line = 0;
      }
    } 
    if(line || fbi->_recvnodes.empty()){
      Ouf << endl;
      line = 0;
    }
    ni = fbi->_sharenodes.begin();
    line = 0;
    while(ni != fbi->_sharenodes.end()){
      Ouf << setw(8) << *ni++;
      line++;
      if(line == 10){
	Ouf << endl;
	line = 0;
      }
    } 
    if(line || fbi->_sharenodes.empty()){
      Ouf << endl;
      line = 0;
    }
    fbi++;
  }
  Ouf << "# End" << endl;
  return(true);
}

bool
TRAIL_FluReadCOM(GEM_Partition &gp)
{
  if(gp._solver_data._string_data.size() < 3){
    if(gp._out)
      *gp._out << "TRAIL_FluReadCOM::Error: No casename or path found."
	       << std::endl;
    return(false);
  }
  
  string pre(gp._solver_data._string_data[1] + "/" + 
	     gp._solver_data._string_data[0]);
  if(gp._solver_data._int_data.empty()){
    if(gp._out)
      *gp._out 
	<< "TRAIL_FluReadCOM: Error,called with empty _int_data."
	<< endl;

    return(false);
  }
  ifstream Inf;
  ostringstream Ostr;
  Ostr << pre << ".com_";
  Ostr << setw(5) << setfill('0');
  Ostr << gp._id;
  Inf.open(Ostr.str().c_str());
  if(!Inf){
    if(gp._out)
      *gp._out 
	<< "TRAIL_FluReadCOM: Error: Cannot open " << Ostr.str() << "."
	<< std::endl;
    return(false);
  }    
  string line;
  SkipLines(Inf,2);
  unsigned int nborders = 0;
  Inf >> nborders;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluReadCOM(" << gp._id << "): nborders: " 
	     << nborders << "\n";
  assert(nborders == gp._pb.size() && nborders == 
	 gp._solver_data._int_data[0].size());
  SkipLines(Inf,2);
  unsigned int pbindex = 0;
  vector<GEM_PartitionBoundary>::iterator fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    unsigned int rpart;
    int rbid;
    Inf >> rpart >> rbid;
    assert(fbi->_rpart == rpart && 
	   gp._solver_data._int_data[0][pbindex++] == rbid);
    fbi++;
  }
  SkipLines(Inf,2);
  fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    unsigned int rcells,scells;
    Inf >> scells >> rcells;
    assert(scells == fbi->_sendcells.size() &&
	   rcells == fbi->_recvcells.size());
    int cell = 0;
    int ncells = fbi->_sendcells.size();
    while(cell < ncells)
      Inf >> fbi->_sendcells[cell++];
    cell = 0;
    ncells = fbi->_recvcells.size();
    while(cell < ncells)
      Inf >> fbi->_recvcells[cell++];
    fbi++;
  }    
  SkipLines(Inf,2);
  fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    unsigned int nshare,nsend,nrecv;
    Inf >> nsend >> nrecv >> nshare;
    assert(nsend  == fbi->_sendnodes.size() &&
	   nrecv  == fbi->_recvnodes.size() &&
	   nshare == fbi->_sharenodes.size());
    unsigned int node = 0;
    unsigned int nnodes = fbi->_sendnodes.size();
    while(node < nnodes)
      Inf >> fbi->_sendnodes[node++];
    node = 0;
    nnodes = fbi->_recvnodes.size();
    while(node < nnodes)
      Inf >> fbi->_recvnodes[node++];
    nnodes = fbi->_sharenodes.size();
    node = 0;
    while(node < nnodes)
      Inf >> fbi->_sharenodes[node++];
    fbi++;
  }
  Inf.close();
  return(true);
}

bool
TRAIL_FluWriteMAP(GEM_Partition &gp,unsigned int nproc,
		 unsigned int nregions,const std::string &path)
{
  if(gp._solver_data._string_data.size() < 3){
    if(gp._out)
      *gp._out << "TRAIL_FluWriteMAP::Error: No casename or path found."
	       << std::endl;
    return(false);
  }
  std::string pre("./" + path + "/" + gp._solver_data._string_data[0]);
  ofstream Ouf;
  ostringstream Ostr;
  Ostr << pre << ".map";
  Ouf.open(Ostr.str().c_str());
  if(!Ouf)
    return(false);
  Ouf << "# ROCFLU region mapping file" << endl
      << "# Number of regions" << endl
      << setw(8) << nregions << endl
      << "# Number of processes" << endl
      << setw(8) << nproc << endl;
  unsigned int regpproc = nregions/nproc;
  if(regpproc == 0)
    return(false);
  unsigned int left     = nregions%nproc;
  unsigned int proc = 0;
  unsigned int reg = 0;
  while(proc < nproc){
    unsigned int nreg = regpproc;
    if(left > 0){
      nreg++;
      left--;
    }
    Ouf << "# Process " << setw(6) << setfill('0') << proc+1 << endl
	<< setfill(' ') << setw(8) << nreg << endl;
    unsigned int uplimit = reg + nreg;
    while(reg < uplimit)
      Ouf << setw(8) << (1+reg++) << endl;
    proc++;
  }
  Ouf << "# End" << endl;
  Ouf.close();
  return(true);
}

bool
TRAIL_FluWriteDIM(GEM_Partition &gp,double t, bool unsteady,
		 const std::string &path)
{
  if(gp._solver_data._string_data.size() < 3){
    if(gp._out)
      *gp._out << "TRAIL_FluWriteDIM::Error: No casename or path found."
	       << std::endl;
    return(false);
  }
  std::string pre("./" + path + "/" + gp._solver_data._string_data[0]);
  std::ofstream Ouf;
  std::ostringstream Ostr;
  Ostr << pre << ".dim_";
  Ostr << setw(5) << setfill('0');
  Ostr << gp._id;
  std::string filebase(Ostr.str());
  if(unsteady){
    Ostr.clear();
    Ostr.str("");
    Ostr << "_" << scientific << setprecision(5) << t;
    std::string timestring(Ostr.str());
    timestring.replace(8,1,"E");
    filebase+=timestring;
  }
  Ouf.open(filebase.c_str());
  if(!Ouf){
    if(gp._out)
      *gp._out << "TRAIL_FluWriteDIM::Error: Could not open " << filebase << "." 
	       << std::endl;
    return(false);
  }
  unsigned int ntet = gp._tetconn.size()/4;
  unsigned int nnodes = gp._nc.size()/3;
  unsigned int nhex = gp._hexconn.size()/8;
  unsigned int npris = gp._prisconn.size()/6;
  unsigned int npyr = gp._pyrconn.size()/5;
  unsigned int nelem = ntet + nhex + npris + npyr;
  Ouf << "# ROCFLU dimensions file" << endl
      << "# Vertices" << endl
      << setw(8) << nnodes - gp._ngnodes << setw(8) << nnodes
      << setw(8) << nnodes+(nnodes/5) << endl
      << "# Cells" << endl
      << setw(8) << (nelem-gp._ngtet-gp._ngpyr-gp._ngpris-gp._nghex)
      << setw(8) << nelem 
      << setw(8) << (nelem+(nelem/5)) << endl
      << "# Tetrahedra" << endl
      << setw(8) << ntet - gp._ngtet << setw(8) << ntet << setw(8) 
      << (ntet+(ntet/5))
      << endl
      << "# Hexahedra" << endl
      << setw(8) << nhex - gp._nghex << setw(8) << nhex << setw(8) 
      << (nhex+(nhex/5))
      << endl
      << "# Prisms" << endl
      << setw(8) << npris - gp._ngpris << setw(8) << npris << setw(8) 
      << npris+(npris/5)
      << endl
      << "# Pyramids" << endl
      << setw(8) << npyr - gp._ngpyr << setw(8) << npyr 
      << setw(8) << npyr+(npyr/5)
      << endl
      << "# Patches" << endl
      << setw(8) << gp._db.size() << setw(8) 
      << gp._solver_data._int_data[1][0] << endl; // total_num_patches
  vector<GEM_DomainBoundary>::iterator fpi = gp._db.begin();
  while(fpi != gp._db.end()){
    Ouf << setw(8) << fpi->_id << setw(8) 
	<< fpi->_triconn.size()/3 - fpi->_ngtri
	<< setw(8) << fpi->_triconn.size()/3 << setw(8) 
	<< fpi->_quadconn.size()/4 - fpi->_ngquad
	<< setw(8) << fpi->_quadconn.size()/4 << endl;
    fpi++;
  }
  Ouf << "# Borders" << endl
      << setw(8) << gp._pb.size() << endl;
  vector<GEM_PartitionBoundary>::iterator fbi = gp._pb.begin();
  unsigned int pbindex = 0;
  assert(gp._solver_data._int_data[0].size() == gp._pb.size());
  while(fbi != gp._pb.end()){
    Ouf << setw(8) << fbi->_rpart << setw(8) 
	<< gp._solver_data._int_data[0][pbindex++] 
	<< setw(8) 
	<< fbi->_sendcells.size() << setw(8) << fbi->_recvcells.size() 
	<< setw(8) 
	<< fbi->_sendnodes.size() << setw(8) << fbi->_recvnodes.size() 
	<< setw(8) 
	<< fbi->_sharenodes.size() << endl;
    fbi++;
  }
  Ouf << "# End" << endl;
  return(true);
}

bool
TRAIL_FluReadDIM(GEM_Partition &gp,double t,bool unsteady)
{
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluReadDIM: Entry" << "\n";
  if(gp._solver_data._string_data.size() < 3){
    if(gp._out)
      *gp._out << "TRAIL_FluReadDIM::Error: No casename or path found."
	       << std::endl;
    return(false);
  }
  string pre(gp._solver_data._string_data[1] + "/" + 
	     gp._solver_data._string_data[0]);
  ifstream Inf;
  ostringstream Ostr;
  Ostr << pre << ".dim_";
  Ostr << setw(5) << setfill('0');
  Ostr << gp._id;
  string filebase(Ostr.str());
  if(unsteady){
    Ostr.clear();
    Ostr.str("");
    Ostr << "_" << scientific << setprecision(5) << t;
    string timestring(Ostr.str());
    timestring.replace(8,1,"E");
    filebase+=timestring;
  }
  Inf.open(filebase.c_str());
  if(!Inf){
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_FluReadDIM: Unable to open " << filebase 
	       << " for reading." << endl;
    return(false);
  }
  string line;
  unsigned int nreal_cells = 0;
  unsigned int ncells = 0;
  unsigned int npatches = 0;
  unsigned int nreal_tets = 0;
  unsigned int nreal_hex = 0;
  unsigned int nreal_pris = 0;
  unsigned int nreal_pyr = 0;
  unsigned int nreal_nodes = 0;
  unsigned int ntet = 0;
  unsigned int nhex = 0;
  unsigned int npris = 0;
  unsigned int npyr = 0;
  unsigned int nnodes = 0;
  SkipLines(Inf,2);
  Inf >> nreal_nodes >> nnodes;
  SkipLines(Inf,2);
  Inf >> nreal_cells >> ncells;
  SkipLines(Inf,2);
  Inf >> nreal_tets >> ntet;
  SkipLines(Inf,2);
  Inf >> nreal_hex >> nhex;
  SkipLines(Inf,2);
  Inf >> nreal_pris >> npris;
  SkipLines(Inf,2);
  Inf >> nreal_pyr >> npyr;
  SkipLines(Inf,2);
  unsigned int nelem = ntet + nhex + npyr + npris;
  Inf >> npatches >> gp._solver_data._int_data[1][0]; // npatches_total
  gp._ngnodes = nnodes - nreal_nodes;
  gp._ngtet =  ntet - nreal_tets;
  gp._nghex =  nhex - nreal_hex;
  gp._ngpyr =  npyr - nreal_pyr;
  gp._ngpris = npris - nreal_pris;
  assert(nelem == ncells);
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluReadDIM: Elems: (" << nelem << "," 
	     << gp._ngtet+gp._nghex+gp._ngpyr+gp._ngpris << ")" << endl
	     << "TRAIL_FluReadDIM: Nodes: (" << nnodes << "," 
	     << gp._ngnodes << ")" << endl
	     << "TRAIL_FluReadDIM: Patches: (" << npatches << "," 
	     << gp._solver_data._int_data[1][0] << ")" << endl;  
  gp._db.resize(npatches);
  gp._tetconn.resize(4*ntet);
  gp._hexconn.resize(8*nhex);
  gp._prisconn.resize(6*npris);
  gp._pyrconn.resize(5*npyr);
  gp._nc.resize(3*nnodes);
  vector<GEM_DomainBoundary>::iterator fpi = gp._db.begin();
  while(fpi != gp._db.end()){
    unsigned int nreal_tri = 0;
    unsigned int nreal_quad = 0;
    unsigned int ntri = 0;
    unsigned int nquad = 0;
    Inf >> fpi->_id >> nreal_tri >> ntri >> nreal_quad >> nquad;
    fpi->_ngtri  = ntri - nreal_tri;
    fpi->_ngquad = nquad - nreal_quad;
    fpi->_triconn.resize(ntri*3);
    fpi->_quadconn.resize(nquad*4);
    fpi++;
  }
  unsigned int nborders = 0;
  SkipLines(Inf,2);
  Inf >> nborders;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluReadDIM: Number of borders: " << nborders 
	     << endl;
  gp._pb.resize(nborders);
  gp._solver_data._int_data[0].resize(nborders);
  unsigned int pbindex = 0;
  vector<GEM_PartitionBoundary>::iterator fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    unsigned int nsend,nrecv,nshared,csend,crecv;
    Inf >> fbi->_rpart 
	>> gp._solver_data._int_data[0][pbindex++] >> csend >> crecv
	>> nsend >> nrecv >> nshared;
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_FluReadDIM: Border: (" 
	       << fbi->_rpart << "," 
	       << gp._solver_data._int_data[0][pbindex - 1] << "," 
	       << csend << "," << crecv << "," << nsend << "," << nrecv 
	       << "," << nshared << ")" << endl;
    fbi->_sharenodes.resize(nshared);
    fbi->_sendnodes.resize(nsend);
    fbi->_recvnodes.resize(nrecv);
    fbi->_sendcells.resize(csend);
    fbi->_recvcells.resize(crecv);
    fbi++;
  }
  Inf.close();
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluReadDIM: Exit" << endl;
  return(true);
}

bool
TRAIL_FluPopRemBordIndFILE(GEM_Partition &gp,double t,bool unsteady,
			  const std::string &path)
{
  if(gp._solver_data._string_data.size() < 3){
    if(gp._out)
      *gp._out << "TRAIL_FluPopRemBordIndFILE:Error: No casename or path found."
	       << std::endl;
    return(false);
  }
  std::string pre("./" + path + "/" + gp._solver_data._string_data[0]);
  unsigned int pbindex = 0;
  std::vector<GEM_PartitionBoundary>::iterator fbi = gp._pb.begin();
  while(fbi != gp._pb.end()){
    std::ifstream Inf;
    std::ostringstream Ostr;
    Ostr << pre << ".dim_";
    Ostr << setw(5) << setfill('0');
    Ostr << fbi->_rpart;
    std::string filebase(Ostr.str());
    if(unsteady){
      Ostr.clear();
      Ostr.str("");
      Ostr << "_" << scientific << setprecision(5) << t;
      std::string timestring(Ostr.str());
      timestring.replace(8,1,"E");
      filebase+=timestring;
    }
    Inf.open(filebase.c_str());
    if(!Inf){
      if(gp._debug && gp._out)
	*gp._out 
	  << "TRAIL_FluPopRemBordIndFILE: Unable to open "
	  << filebase << " for reading.\n";
      return(false);
    }
    string line;
    unsigned int fakeout = 0;
    unsigned int nreal_cells = 0;
    unsigned int nnodes = 0;
    unsigned int nhex = 0;
    unsigned int ntet = 0;
    unsigned int npris = 0;
    unsigned int npyr = 0;
    unsigned int ncells = 0;
    unsigned int npatches = 0;
    unsigned int nreal_tets = 0;
    unsigned int nreal_hex = 0;
    unsigned int nreal_pris = 0;
    unsigned int nreal_pyr = 0;
    unsigned int nreal_nodes = 0;
    unsigned int npatches_total = 0;
    unsigned int ngnodes = 0;
    unsigned int ngtet = 0;
    unsigned int nghex = 0;
    unsigned int ngpris = 0;
    unsigned int ngpyr = 0;
    SkipLines(Inf,2);
    Inf >> nreal_nodes >> nnodes;
    SkipLines(Inf,2);
    Inf >> nreal_cells >> ncells;
    SkipLines(Inf,2);
    Inf >> nreal_tets >> ntet;
    SkipLines(Inf,2);
    Inf >> nreal_hex >> nhex;
    SkipLines(Inf,2);
    Inf >> nreal_pris >> npris;
    SkipLines(Inf,2);
    Inf >> nreal_pyr >> npyr;
    SkipLines(Inf,2);
    Inf >> npatches >> npatches_total;
    ngnodes = nnodes - nreal_nodes;
    ngtet = ntet - nreal_tets;
    nghex = nhex - nreal_hex;
    ngpyr = npyr - nreal_pyr;
    ngpris = npris - nreal_pris;
    if(gp._debug && gp._out)
      *gp._out  << "GEM_Partition(" << gp._id 
		<< ")::PopRemBordIndFILE::Elems: (" 
		<< ntet+npris+nhex+npyr
		<< "," << ngtet+nghex+ngpyr+ngpris << ")\n"
		<< "GEM_Partition(" << gp._id << ")::PopRemBordIndFILE" 
		<< "::Nodes: (" 
		<< nnodes << "," 
		<< ngnodes << ")\n"
		<< "GEM_Partition(" << gp._id 
		<< ")::PopRemBordIndFILE:: Patches: (" 
		<< npatches << "," 
		<< npatches_total << ")\n";
    for(unsigned int pind = 0;pind < npatches;pind++){
      unsigned int nreal_tri = 0;
      unsigned int nreal_quad = 0;
      Inf >> fakeout >> nreal_tri >> fakeout >> nreal_quad >> fakeout;
    }
    unsigned int nborders = 0;
    SkipLines(Inf,2);
    Inf >> nborders;
    if(gp._debug && gp._out)
      *gp._out  << "GEM_Partition(" << gp._id << ")::PopRemBordIndFILE:: " 
		<< "Number of borders: " << nborders 
		<< "\n";
    bool done = false;
    for(unsigned int bind = 0;(bind < nborders && !done);bind++){
      unsigned int rpart;
      Inf >> rpart;
      if(rpart == gp._id){
	gp._solver_data._int_data[0][pbindex++] = bind+1;
	done = true;
      }
      else
	Inf >> fakeout >> fakeout >> fakeout >> fakeout >> fakeout >> fakeout;
    }
    if(!done)
      return(false);
    Inf.close();
    //    this->WriteFluDIM(prefix,t,unsteady);
    fbi++;
  }
  return(true);
}

bool 
TRAIL_FluWriteCMP(GEM_Partition &gp,const std::string &path)
{
  if(gp._solver_data._string_data.size() < 3){
    if(gp._out)
      *gp._out << "TRAIL_FluWriteCMP:Error: No casename or path found."
	       << std::endl;
    return(false);
  }
  std::string pre("./" + path + "/" + gp._solver_data._string_data[0]);
  std::ofstream Ouf;
  std::ostringstream Ostr;
  Ostr << pre << ".cmp_";
  Ostr << setw(5) << setfill('0');
  Ostr << gp._id;
  Ouf.open(Ostr.str().c_str());
  if(!Ouf)
    return(false);
  unsigned int ntet  = gp._tetconn.size()/4;
  unsigned int nhex  = gp._hexconn.size()/8;
  unsigned int npris = gp._prisconn.size()/6;
  unsigned int npyr  = gp._pyrconn.size()/5;
  //  BuildCellMapping();
  Ouf << "# ROCFLU cell mapping file" << endl
      << "# Dimensions" << endl
      << setw(8) << ntet << setw(8) << nhex << setw(8) << npris 
      << setw(8) << npyr << endl;
  //  vector<pair<unsigned int,unsigned int> >::iterator ci = _cellmap.begin();
  //  int line = 0;
  unsigned int el = 0;
  unsigned int npl = 10;
  if(ntet > 0){
    el = 0;
    Ouf << "# Tetrahedra" << endl;
    while(el < ntet){
      Ouf << setw(8) << gp.Elem2Cell(make_pair(1,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  if(nhex > 0){
    el = 0;
    Ouf << "# Hexahedra" << endl;
    while(el < nhex){
      Ouf << setw(8) << gp.Elem2Cell(make_pair(4,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  if(npris > 0){
    el = 0;
    Ouf << "# Prisms" << endl;
    while(el < npris){
      Ouf << setw(8) << gp.Elem2Cell(make_pair(3,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  if(npyr > 0){
    el = 0;
    Ouf << "# Pyramids" << endl;
    while(el<npyr){
      Ouf << setw(8) << gp.Elem2Cell(make_pair(2,++el));
      if(!(el%npl))
	Ouf << endl;
    }
    if(el%npl) Ouf << endl;
  }
  Ouf << "# End" << endl;
  Ouf.close();
  return(true);
}

void
TRAIL_FluResizeVolSoln(GEM_Partition &gp)
{
  unsigned int ncells = gp._tetconn.size()/4 + gp._prisconn.size()/6 +
    gp._pyrconn.size()/5 + gp._hexconn.size()/8;
  unsigned int nnodes = gp._nc.size()/3;
  unsigned int nvfaces = gp._nvface;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluResizeVolSoln: Allocating for " << ncells << " cells"
	     << ", " << nnodes << " nodes and " << nvfaces << " volume faces."
	     << std::endl; 
  gp._data._int_data.resize(1);
  gp._data._stride_int.resize(1);
  gp._data._int_data[0].resize(1); // unsteady flag
  gp._data._int_data[0][0] = 1;  // unsteady yes 
  gp._data._field_data.resize(10);
  gp._data._stride_field.resize(10);
  gp._data._field_data[0].resize(1); // _current_time
  gp._data._stride_field[0] = 1;
  gp._data._field_data[0][0] = 0.0;
  gp._data._field_data[1].resize(1); // _residual
  gp._data._field_data[1][0] = 0.0;
  gp._data._stride_field[1] = 1;
  gp._data._field_data[2].resize(ncells,0.0); // rhof
  gp._data._stride_field[2] = 1;
  gp._data._field_data[3].resize(3*ncells,0.0); // rhovf
  gp._data._stride_field[3] = 3;
  gp._data._field_data[4].resize(ncells,0.0); //rhoEf
  gp._data._stride_field[4] = 1;
  gp._data._field_data[5].resize(ncells,0.0); // pf
  gp._data._stride_field[5] = 1;
  gp._data._field_data[6].resize(ncells,0.0); // Tf
  gp._data._stride_field[6] = 1;
  gp._data._field_data[7].resize(ncells,0.0); // af
  gp._data._stride_field[7] = 1;
  gp._data._field_data[8].resize(nnodes*3,0.0); // disp
  gp._data._stride_field[8] = 3;
  gp._data._field_data[9].resize(nvfaces,0.0); // gridspeeds
  gp._data._stride_field[9] = 1;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluResizeVolSoln: Volume data allocated." << std::endl;
}

bool
TRAIL_FluInitVolSoln(GEM_Partition &gp)
{
  if(gp._solver_data._string_data.size() < 3 ||
     gp._solver_data._string_data[1].empty()){
    if(gp._out)
      *gp._out << "  TRAIL_FluInitVolSoln::Error: Solver data contains no path."
	       << std::endl;
    
    return(false);
  }
  string pre(gp._solver_data._string_data[1]+"/"+
	     gp._solver_data._string_data[0]);
  string fname(pre+".inp");
  ifstream Inf;
  Inf.open(fname.c_str());
  if(!Inf){
    if(gp._out)
      *gp._out << " TRAIL_FluInitVolSoln::Error: Could not open " << fname 
	       << "." << std::endl;
    
    return false;
  }
  unsigned int ntet  = gp._tetconn.size()/4;
  unsigned int npyr  = gp._pyrconn.size()/5;
  unsigned int nhex  = gp._hexconn.size()/8;
  unsigned int npris = gp._prisconn.size()/6;
  unsigned int nelem = ntet + npyr + npris + nhex;
  //  unsigned int nvert = gp._nc.size()/3;
  TRAIL_FluResizeVolSoln(gp);
  string line;
  double cp    = -1.0;
  double gamma = -1.0;
  double velx  =  0.0;
  double vely  =  0.0;
  double velz  =  0.0;
  double press = -1.0;
  double rho   = -1.0;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluInitVolSoln: Initializing volume data." << std::endl;
  while(getline(Inf,line)){
    if(line[0] != '#'){
      string::size_type xtime  = line.find("STARTTIME");
      string::size_type ftype  = line.find("FLOWTYPE");
      string::size_type xgamma = line.find("GAMMA");
      string::size_type xpress = line.find("PRESS");
      string::size_type x_dens = line.find("DENS");
      string::size_type x_velx = line.find("VELX");
      string::size_type x_vely = line.find("VELY");
      string::size_type x_velz = line.find("VELZ");
      string::size_type x_cp   = line.find("CP");
      string tstr;
      istringstream Istr(line);
      if(xtime  != string::npos)
	Istr >> tstr >> gp._data._field_data[0][0]; // current time
      if(ftype  != string::npos)
	Istr >> tstr >> gp._data._int_data[0][0]; // flowtype (steady/unsteady)
      if(xgamma != string::npos)
	Istr >> tstr >> gamma;
      if(xpress != string::npos)
	Istr >> tstr >> press;
      if(x_dens != string::npos)
	Istr >> tstr >> rho;
      if(x_velx != string::npos)
	Istr >> tstr >> velx;
      if(x_vely != string::npos)
	Istr >> tstr >> vely;
      if(x_velz != string::npos)
	Istr >> tstr >> velz;
      if(x_cp   != string::npos)
	Istr >> tstr >> cp;
    }
  }
  if(gamma < 0 || press < 0 || rho < 0 || cp < 0)
    return(false);
  gp._data._field_data[2].resize(nelem,rho); // rhof
  gp._data._stride_field[2] = 1;
  double E = press/(rho*(gamma - 1.0)) + 
    (velx*velx + vely*vely + velz*velz)/2.0;
  gp._data._field_data[3].resize(3*nelem,0.0); // rhovf
  gp._data._stride_field[3] = 3;
  if(velx != 0.0 || vely != 0.0 || velz != 0.0){
    unsigned int index = 0;
    for(unsigned int el = 0;el < nelem;el++){
      gp._data._field_data[3][index++] = rho*velx;
      gp._data._field_data[3][index++] = rho*vely;
      gp._data._field_data[3][index++] = rho*velz;
    }
  }
  double R = (cp*(gamma-1.0))/gamma;
  double Ti = press/(rho*R);
  double C  = sqrt(gamma*press/rho);
  gp._data._field_data[4].resize(nelem,rho*E); // rhoEf
  gp._data._stride_field[4] = 1;
  gp._data._field_data[5].resize(nelem,press); // pf
  gp._data._stride_field[5] = 1;
  gp._data._field_data[6].resize(nelem,Ti); // Tf
  gp._data._stride_field[6] = 1;
  gp._data._field_data[7].resize(nelem,C); // af
  gp._data._stride_field[7] = 1;
  //  _data._field_data[8].resize(gp._nc.size(),0.0); // disp
  //  _data._stride_field[8] = 3;
  //  gp._data._field_data[9].resize(gp._nvface,0.0); // gridspeeds (gsp)
  //  gp._data._stride_field[9] = 1;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluInitVolSoln: Volume data initialized." << std::endl;
  return(true);
}


void
TRAIL_FluResizeSurfSoln(GEM_DomainBoundary &db)
{
  unsigned int nfaces = db._triconn.size()/3 + db._quadconn.size()/4;
  unsigned int nvert = db.NNodes();
  if(db._debug && db._out)
    *db._out << "TRAIL_FluResizeSurfSoln: Resizing with " << nfaces << " faces"
	     << " and " << nvert << " vertices." << std::endl;
  db._data._field_data.resize(21);
  db._data._stride_field.resize(21);
  db._data._int_data.resize(3);
  db._data._stride_int.resize(3);
  db._data._int_data[1].resize(1); // bcflag
  db._data._stride_int[1] = 1;
  db._data._int_data[2].resize(1); // cnstr_type
  db._data._stride_int[2] = 1;
  db._data._field_data[0].resize(nfaces,0.0);    // gsp
  db._data._stride_field[0] = 1;
  db._data._field_data[1].resize(3*nfaces,0.0);  // rhofvf_alp
  db._data._stride_field[1] = 3;
  db._data._field_data[2].resize(3*nfaces,0.0);  // nf_alp
  db._data._stride_field[2] = 3;
  db._data._field_data[3].resize(nfaces,0.0);    // rhof_alp
  db._data._stride_field[3] = 1;
  db._data._field_data[4].resize(nfaces,0.0);    // pf
  db._data._stride_field[4] = 1;
  db._data._field_data[5].resize(nfaces,0.0);    // qc
  db._data._stride_field[5] = 1;
  db._data._field_data[6].resize(nfaces,0.0);    // qr
  db._data._stride_field[6] = 1;
  db._data._field_data[7].resize(3*nfaces,0.0);  // tf
  db._data._stride_field[7] = 3;
  db._data._field_data[8].resize(nfaces,0.0);    // Tb_alp
  db._data._stride_field[8] = 1;
  db._data._field_data[9].resize(nfaces,0.0);    // mdot_alp
  db._data._stride_field[9] = 1;
  db._data._field_data[10].resize(nfaces,0.0);   // Tflm_alp
  db._data._stride_field[10] = 1;
  db._data._int_data[0].resize(nfaces,0);     // bflag
  db._data._stride_int[0] = 1;
  db._data._field_data[11].resize(nfaces,0.0);   // Tf
  db._data._stride_field[11] = 1;
  db._data._field_data[12].resize(3*nvert,0.0); // vm
  db._data._stride_field[12] = 3;
  db._data._field_data[13].resize(3*nfaces,0.0); // vs
  db._data._stride_field[13] = 3;
  db._data._field_data[14].resize(3*nfaces,0.0); // vs_old
  db._data._stride_field[14] = 3;
  db._data._field_data[15].resize(nfaces,0.0);   // ts
  db._data._stride_field[15] = 1;
  db._data._field_data[16].resize(nfaces,0.0);   // mdot
  db._data._stride_field[16] = 1;
  db._data._field_data[17].resize(nfaces,0.0);   // mdot_old
  db._data._stride_field[17] = 1;
  db._data._field_data[18].resize(nvert*3,0.0);  // du_alp
  db._data._stride_field[18] = 3;
  db._data._field_data[19].resize(3*nvert,0.0);  // nc_t0
  db._data._stride_field[19] = 3;
  db._data._field_data[20].resize(nvert,0.0);    // sq_dist
  db._data._stride_field[20] = 1;
  if(db._debug && db._out)
    *db._out << "TRAIL_FluResizeSurfSoln: Resizing done." << std::endl;
}

bool 
TRAIL_FluInitSurfSoln(GEM_DomainBoundary &db,const std::string &prefix,bool all)
{
  string pre(prefix);
  if(pre.empty())
    return(false);
  unsigned int ntri   = db._triconn.size()/3;
  unsigned int nquad  = db._quadconn.size()/4;
  unsigned int nfaces = ntri + nquad;
  //  unsigned int nvert  = db.NNodes();
  //  if(all){
  // FIXME FIXME
  //#ifdef _ROCSTAR_X_
  TRAIL_FluResizeSurfSoln(db);
  //#else
  // Rocflu only needs only gsp for field data
  //    nfaces = (ntri+nquad)- (_ngtri + _ngquad);
  //    TRAIL_FluResizeSurfSoln(nfaces,nvert,db);
  //    for(int i = 1;i < 22;i++){
  //      db._data._field_data[i].resize(0);
  //      db._data._stride_field[i] = 0;
  //    }
  //#endif
  //  }
  //#ifdef _ROCSTAR_X_
  int constraint_map[] = { 2 , 120, 121, 122, -122, -121, -120, 0};
  // Populate the Rocstar data living on the surface
  string fname(pre+".bc");
  ifstream Inf;
  Inf.open(fname.c_str());
  if(!Inf){
    if(db._out)
      *db._out << "TRAIL_FluInitSurfSoln: Error: could not locate bc file, "
	       << fname << "." << std::endl;
    return false;
  }
  string line;
  bool found_patch = false;
  while(getline(Inf,line) && !found_patch){
    if(line[0] != '#'){
      string::size_type patchx = line.find("PATCH");
      if( patchx != string::npos){
	string tstr;
	unsigned int ubound,lbound;
	istringstream Istr(line);
	Istr >> tstr >> lbound >> ubound;
	if(db._id >= lbound && db._id <= ubound)
	  found_patch = true;
      }
    }
  }
  if(!found_patch){
    if(db._out)
      *db._out << "TRAIL_FluInitSurfSoln: Error: "
	       << "Did not find patch " << db._id << " in bc file." << endl;
    return(false);
  }
  // Now we want to search forward for the relevant data until we encounter #
  bool end_section = false;
  db._data._int_data[1][0] = 2; //_bcflag      = 2;
  db._data._int_data[2][0] = 2; // _constr_type = 2;
  db._data._stride_int[1] = 1;
  db._data._stride_int[2] = 1;
  while(getline(Inf,line) && !end_section){
    if(line[0] == '#') end_section = true;
    else{
      string::size_type movx  = line.find("MOVEDIR");
      string::size_type coupx = line.find("COUPLED");
      string::size_type burn  = line.find("BFLAG");
      if(movx != string::npos){
	string tstr;
	int mvdir;
	istringstream Istr(line);
	Istr >> tstr >> mvdir;
	db._data._int_data[2][0] = constraint_map[mvdir];
      }
      if(coupx != string::npos){
	string tstr;
	istringstream Istr(line);
	Istr >> tstr >> db._data._int_data[1][0];
      }
      if(burn != string::npos){
	string tstr;
	unsigned int burnf = 0;
	istringstream Istr(line);
	Istr >> tstr >> burnf;
	db._data._int_data[0].resize(nfaces,burnf);
	db._data._stride_int[0] = 1;
      }
    }
  }
  //#endif
  Inf.close();
  if(db._debug && db._out)
    *db._out << "TRAIL_FluInitSurfSoln: Exit" << endl;
  return true;
}

int
TRAIL_FluNumPatches(GEM_Partition &gp)
{
  if(gp._solver_data._string_data.size() < 3 ||
     gp._solver_data._string_data[0].empty() ||
     gp._solver_data._string_data[1].empty()){
    if(gp._out)
      *gp._out << "TRAIL_FluNumPatches:Error: solver data has no casename"
	       << " or path." << std::endl;
    return 0;
  } 
  string pre(gp._solver_data._string_data[1]+ "/" +
	     gp._solver_data._string_data[0]);
  string filename(pre+".bc");
  ifstream Inf;
  Inf.open(filename.c_str());
  string line;
  int npatches = 0;
  if(!Inf) return 0;
  while(!Inf.eof()){
    bool found = false;
    while(!found){
      getline(Inf,line);
      if(line[0] == '#')
	found = true;
    }
    found = false;
    while(getline(Inf,line) && !found){
      if(line[0] != '#'){
	string::size_type patchx = line.find("PATCH");
	if( patchx != string::npos){
	  string tstr;
	  unsigned int ubound,lbound;
	  istringstream Istr(line);
	  Istr >> tstr >> lbound >> ubound;
	  npatches += ((ubound - lbound) + 1);
	  found = true;
	}
      }
    }
  }
  Inf.close();
  return npatches;  
}

bool
TRAIL_FluInitSurfSoln(GEM_Partition &gp)
{
  if(gp._solver_data._string_data.size() < 3 || 
     gp._solver_data._string_data[1].empty()){
    if(gp._out)
      *gp._out << "TRAIL_FluInitSurfSoln::Error: Solver data contains no"
	       << " path." << std::endl;
    return(false);
  }
  unsigned int ndb = gp._db.size();
  unsigned int n = 0;
  std::string pre(gp._solver_data._string_data[1]+"/"+
		  gp._solver_data._string_data[0]);
  while(n < ndb){
    if(!TRAIL_FluInitSurfSoln(gp._db[n++],pre,true))
      return(false);
  }
  return(true);
}

bool
TRAIL_FluRegisterVolSoln(GEM_Partition &gp,bool all)
{
  gp.Create_com_volsoln("rhof",gp._data._field_data[2],1,"kg/m^3");
  gp.Create_com_volsoln("rhovf",gp._data._field_data[3],3,"kg/(m^2s)");
  gp.Create_com_volsoln("rhoEf",gp._data._field_data[4],1,"J/m^3");
  COM_new_dataitem((gp.volume_window+".gs"),'p',COM_DOUBLE,1,"m/s");
  COM_set_size((gp.volume_window+".gs"),gp.pane_id,gp._nvface);
  if(gp._data._field_data[9].size() != 0)
    COM_set_array((gp.volume_window+".gs"),gp.pane_id,
		  &(gp._data._field_data[9][0]));
  gp.Create_com_volsoln("pf",gp._data._field_data[5],1,"Pa");
  gp.Create_com_volsoln("Tf",gp._data._field_data[6],1,"K");
  gp.Create_com_volsoln("af",gp._data._field_data[7],1,"m/s");
  return(true);
}

bool
TRAIL_FluRegisterSurfSoln(GEM_DomainBoundary &db,const string &wname,bool all)
{
  if(db._debug && db._out)
    *db._out << "TRAIL_FluRegisterSoln: "
	     << "Registering pane " << db.pane_id << "." << endl;
  if(db._data._field_data[0].size() != 0){
    COM_set_array((wname+".gs"),db.pane_id,&(db._data._field_data[0][0]));
  }
  COM_set_size((wname+".patchNo"),db.pane_id,1);
  COM_set_array((wname+".patchNo"),db.pane_id,&db._id);
  COM_set_size((wname+".bcflag"),db.pane_id,1);
  COM_set_array((wname+".bcflag"),db.pane_id,&db._data._int_data[1][0]);
  COM_set_size((wname+".cnstr_type"),db.pane_id,1);
  COM_set_array((wname+".cnstr_type"),db.pane_id,&db._data._int_data[2][0]);
  if(db._data._int_data[0].size() > 0)
    COM_set_array((wname+".bflag"),db.pane_id,&(db._data._int_data[0]));
  if(all){
    db.Create_com_surfsoln(wname,"rhofvf_alp",db._data._field_data[1],
			   3,"kg/(m^2s)");
    db.Create_com_surfsoln(wname,"nf_alp",db._data._field_data[2],3,"");
    db.Create_com_surfsoln(wname,"rhof_alp",db._data._field_data[3],
			   1,"kg/m^3");
    db.Create_com_surfsoln(wname,"pf",db._data._field_data[4],1,"Pa");
    db.Create_com_surfsoln(wname,"qc",db._data._field_data[5],1,"W/m^2");
    db.Create_com_surfsoln(wname,"qr",db._data._field_data[6],1,"W/m^2");
    db.Create_com_surfsoln(wname,"tf",db._data._field_data[7],3,"Pa");
    db.Create_com_surfsoln(wname,"Tb_alp",db._data._field_data[8],1,"K");
    db.Create_com_surfsoln(wname,"mdot_alp",db._data._field_data[9],
			   1,"kg/(m^2s)");
    db.Create_com_surfsoln(wname,"Tflm_alp",db._data._field_data[10],1,"K");
    db.Create_com_surfsoln(wname,"Tf",db._data._field_data[11],1,"K");
    db.Create_com_surfsoln(wname,"vs",db._data._field_data[13],3,"m/s");
    db.Create_com_surfsoln(wname,"vs_old",db._data._field_data[14],3,"m/s");
    db.Create_com_surfsoln(wname,"ts",db._data._field_data[15],1,"Pa");
    db.Create_com_surfsoln(wname,"mdot",db._data._field_data[16],
			   1,"kg/(m^2s)");
    db.Create_com_surfsoln(wname,"mdot_old",db._data._field_data[17],
			   1,"kg/(m^2s)");
    // Nodal Data
    if(db._data._field_data[19].size() > 0)
      COM_set_array((wname+".nc_t0"),db.pane_id,&db._data._field_data[19][0]);
    if(db._data._field_data[18].size() > 0)
      COM_set_array((wname+".du_alp"),db.pane_id,&db._data._field_data[18][0]);
    if(db._data._field_data[12].size() > 0)
      COM_set_array((wname+".vm"),db.pane_id,&db._data._field_data[12][0]);
    if(db._data._field_data[20].size() > 0)
      COM_set_array((wname+".sq_dist"),db.pane_id,
		    &db._data._field_data[20][0]);

  }
  return(true);
}

bool
TRAIL_FluRegisterSurfSoln(GEM_Partition &gp,bool all)
{
  
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluRegisterSurfSoln Enter" 
	     << endl;
  COM_new_dataitem((gp.surface_window+".gs"),'e',COM_DOUBLE,1,"m/s");
  COM_new_dataitem((gp.surface_window+".bcflag"),'p',COM_INT,1,"");
  COM_new_dataitem((gp.surface_window+".cnstr_type"),'p',COM_INT,1,"");
  COM_new_dataitem((gp.surface_window+".patchNo"),'p',COM_INT,1,"");
  COM_new_dataitem((gp.surface_window+".bflag"),'e',COM_INT,1,"");
  if(all){
    COM_new_dataitem(gp.surface_window+".rhofvf_alp",'e',
		      COM_DOUBLE,3,"kg/(m^2s)");
    COM_new_dataitem(gp.surface_window+".nf_alp",'e',COM_DOUBLE,3,"");
    COM_new_dataitem(gp.surface_window+".rhof_alp",'e',COM_DOUBLE,1,"kg/m^3");
    COM_new_dataitem(gp.surface_window+".pf",'e',COM_DOUBLE,1,"Pa");
    COM_new_dataitem(gp.surface_window+".qc",'e',COM_DOUBLE,1,"W/m^2");
    COM_new_dataitem(gp.surface_window+".qr",'e',COM_DOUBLE,1,"W/m^2");
    COM_new_dataitem(gp.surface_window+".tf",'e',COM_DOUBLE,3,"Pa");
    COM_new_dataitem(gp.surface_window+".Tb_alp",'e',COM_DOUBLE,1,"K");
    COM_new_dataitem(gp.surface_window+".mdot_alp",'e',COM_DOUBLE,1,
		      "kg/(m^2s)");
    COM_new_dataitem(gp.surface_window+".Tflm_alp",'e',COM_DOUBLE,1,"K");
    COM_new_dataitem(gp.surface_window+".Tf",'e',COM_DOUBLE,1,"K");
    COM_new_dataitem(gp.surface_window+".ts",'e',COM_DOUBLE,1,"Pa");
    COM_new_dataitem(gp.surface_window+".mdot",'e',COM_DOUBLE,1,"kg/(m^2s)");
    COM_new_dataitem(gp.surface_window+".mdot_old",'e',COM_DOUBLE,1,
		      "kg/(m^2s)");
    COM_new_dataitem(gp.surface_window+".vs",'e',COM_DOUBLE,3,"m/s");
    COM_new_dataitem(gp.surface_window+".vs_old",'e',COM_DOUBLE,3,"m/s");
    
    // Nodal data
    COM_new_dataitem((gp.surface_window+".nc_t0"),'n',COM_DOUBLE,3,"m");
    COM_new_dataitem((gp.surface_window+".vm"),'n',COM_DOUBLE,3,"m/s");
    COM_new_dataitem((gp.surface_window+".sq_dist"),'n',COM_DOUBLE,1,"m");
    COM_new_dataitem((gp.surface_window+".du_alp"),'n',COM_DOUBLE,3,"m");
  }
  vector<GEM_DomainBoundary>::iterator fpi = gp._db.begin();
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluRegisterSurfSoln: "
	     << "Registering domain boundaries..." << endl;
  while(fpi != gp._db.end()){
    GEM_DomainBoundary &fp = *fpi++;
    TRAIL_FluRegisterSurfSoln(fp,gp.surface_window,all);
  }
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_FluRegisterSurfSoln: "
	     << "Registering domain boundaries done." << endl;
  return(true);
}

bool
TRAIL_FluRegisterSurfMesh(GEM_Partition &gp)
{
  COM_new_dataitem((gp.surface_window+".t3g:real"),'p',COM_INT,3,"");
  COM_new_dataitem((gp.surface_window+".t3g:virtual"),'p',COM_INT,3,"");      
  COM_new_dataitem((gp.surface_window+".q4g:real"),'p',COM_INT,4,"");
  COM_new_dataitem((gp.surface_window+".q4g:virtual"),'p',COM_INT,4,"");      

  vector<GEM_DomainBoundary>::iterator fpi = gp._db.begin();
  while(fpi != gp._db.end()){
    GEM_DomainBoundary &fp = *fpi++;
    unsigned int nreal = fp._triconn.size()/3 - fp._ngtri;
    if(nreal > 0){
      COM_set_size((gp.surface_window+".t3g:real"),fp.pane_id,nreal);    
      COM_set_array((gp.surface_window+".t3g:real"),fp.pane_id,
		    &(fp._triconn[0]),3);
    }
    nreal = fp._quadconn.size()/4 - fp._ngquad;
    if(nreal > 0){
      COM_set_size((gp.surface_window+".q4g:real"),fp.pane_id,nreal);
      COM_set_array((gp.surface_window+".q4g:real"),fp.pane_id,
		    &(fp._quadconn[0]),4);
    }
    if(fp._ngtri > 0){
      nreal = fp._triconn.size()/3 - fp._ngtri;
      COM_set_size((gp.surface_window+".t3g:virtual"),fp.pane_id,
		   fp._ngtri,fp._ngtri);
      COM_set_array((gp.surface_window+".t3g:virtual"),fp.pane_id,
		    &(fp._triconn[nreal*3]),3);
    }
    if(fp._ngquad > 0){
      nreal = fp._quadconn.size()/4 - fp._ngquad;
      COM_set_size((gp.surface_window+".q4g:virtual"),fp.pane_id,
		   fp._ngquad,fp._ngquad);
      COM_set_array((gp.surface_window+".q4g:virtual"),fp.pane_id,
		    &(fp._quadconn[nreal*4]),4);
    }
  }
  return(true);
}

bool 
TRAIL_FluBuildPatchMapping(map<unsigned int,unsigned int> &pmap,
			  const string &prefix)
{
  // Build the patch mapping from the hand produced patch mapping file
  ostringstream Ostr;
  string pre(prefix);
  //  if(pre.empty()){
  //     if(_casename.empty())
  //       pre = "./";
  //     else
  //       pre = _casepath + "/" + _casename;
  //}
  Ostr << pre << ".cgi";
  ifstream bcfile;
  bcfile.open(Ostr.str().c_str());
  if(!bcfile)
    return(false);
  //  map< unsigned int,unsigned int > pmap;
  unsigned int nggpmap_lines;
  unsigned int npatches_total;
  // _solver_data._int_data[1][0] is _npatches_total
  bcfile >> npatches_total >> nggpmap_lines;
  unsigned int patch_line = 0;
  while(patch_line < nggpmap_lines){
    unsigned int low, high, flupatch;
    bcfile >> low >> high >> flupatch;
    unsigned int p = low;
    while(p <= high)
      pmap[p++] = flupatch;
    patch_line++;
  }
  //  if(_debug && _out){
  //    *_out << "BuildFluPatchMapping: Patch Mapping:\n";
  //    map<unsigned int,unsigned int>::iterator mi = pmap.begin();
  //    while(mi != pmap.end()){
  //      *_out << "BuildFluPatchMapping: (" << mi->first << "," 
  //	   << mi->second << ")\n";
  //      mi++;
  //    }
  //  }
  bcfile.close();
  return(true);
}  


bool 
TRAIL_FluPopulatePatches(GEM_Partition &gp)
{
  if(gp._solver_data._string_data.size() < 3 ||
     gp._solver_data._string_data[0].empty() ||
     gp._solver_data._string_data[1].empty()){
    if(gp._out)
      *gp._out << "TRAIL_FluPopulatePatches:Error: solver data has no casename"
	       << " or path." << std::endl;
    return(false);
  } 
  string pre(gp._solver_data._string_data[1]+ "/" +
	     gp._solver_data._string_data[0]);
  map<unsigned int,unsigned int> patch_mapping;
  TRAIL_FluBuildPatchMapping(patch_mapping,pre);
  gp.MapDomainBoundaries(patch_mapping);
  return(true);
}

bool 
TRAIL_FluWriteNative(GEM_Partition &gp,const std::string &path)
{
  bool error = false;
  if(gp._solver_data._string_data.size() < 3)
    error = true;
  if(gp._solver_data._string_data[0].empty() || 
     gp._solver_data._string_data[1].empty())
    error = true;
  if(gp._solver_data._int_data.empty())
    error = true;
  if(error){
    if(gp._out)
      *(gp._out) << "TRAIL_FluWriteNative: Error.  Called with empty" 
		 << " solver data." << endl;
    return false;
  }
  if(!TRAIL_FluWriteCMP(gp,path))
    error = true;
  if(gp._solver_data._int_data[0].empty()){
    gp._solver_data._int_data[0].resize(gp._pb.size(),0);
    gp._solver_data._stride_int[0] = gp._pb.size();
  }
  if(!TRAIL_FluWriteDIM(gp,0.0,true,path))
    error = true;
  if(!TRAIL_FluWriteCOM(gp,path))
    error = true;
  return(!error);
} 



// Populates our connectivites from the ASCII grid files
// pre should contain the path + casename
// Ex: /my/directory/StarSlice
bool
TRAIL_FluReadGridASCII(GEM_Partition &gp,double t,bool unsteady)
{
  if(gp._solver_data._string_data.size() < 3 ||
     gp._solver_data._string_data[0].empty() ||
     gp._solver_data._string_data[1].empty()){
    if(gp._out)
      *gp._out << "TRAIL_FluReadGridASCII:Error: solver data has no casename"
	       << " or path." << std::endl;
    return 0;
  } 
  string pre(gp._solver_data._string_data[1]+ "/" +
	     gp._solver_data._string_data[0]);
  ostringstream Ostr;
  Ostr << pre << ".grda_" << setw(5) << setfill('0') << gp._id;
  string filebase(Ostr.str());
  Ostr.clear();
  Ostr.str("");
  Ostr << "_" << scientific << setprecision(5) << t;
  string timestring(Ostr.str());
  timestring.replace(8,1,"E"); 
  if(unsteady) filebase+=timestring;
  ifstream Inf;
  Inf.open(filebase.c_str());
  if(!Inf){
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_FluReadGridASCII: Unable to open " << filebase 
	       << " for reading.\n";
    return(false);
  }
  string line;
  SkipLines(Inf,6);
  unsigned int nnode,ntet,nhex,npris,npyr;
  Inf >> nnode >> ntet >> nhex >> npris >> npyr;
  assert(nnode == gp._nc.size()/3       && 
	 ntet  == gp._tetconn.size()/4  && 
  	 nhex  == gp._hexconn.size()/8  && 
	 npris == gp._prisconn.size()/6 &&
  	 npyr  == gp._pyrconn.size()/5);
  SkipLines(Inf,2);
  for(int i = 0; i < 3;i++){
    unsigned int node = 0;
    while(node < nnode)
      Inf >> gp._nc[(node++)*3+i];
  }
  if(ntet > 0){
    gp._tetconn.resize(ntet*4);
    SkipLines(Inf,2);
    for(int i = 0;i < 4;i++){
      unsigned int el = 0;
      while(el < ntet)
	Inf >> gp._tetconn[(el++)*4+i];
    }
  }
  if(nhex > 0){
    gp._hexconn.resize(nhex*8);
    SkipLines(Inf,2);
    for(int i = 0;i < 8;i++){
      unsigned int el = 0;
      while(el < nhex)
	Inf >> gp._hexconn[(el++)*8+i];
    }
  }
  if(npris > 0){
    gp._prisconn.resize(npris*6);
    SkipLines(Inf,2);
    for(int i = 0;i < 6;i++){
      unsigned int el = 0;
      while(el < npris)
	Inf >> gp._prisconn[(el++)*6+i];
    }
  }
  if(npyr > 0){
    gp._pyrconn.resize(npyr*5);
    SkipLines(Inf,2);
    for(int i = 0;i < 5;i++){
      unsigned int el = 0;
      while(el < npyr)
	Inf >> gp._pyrconn[(el++)*5+i];
    }
  }
  SkipLines(Inf,2);
  unsigned int npatches = 0;
  unsigned int patch    = 0;
  Inf >> npatches;
  assert(npatches == gp._db.size());
  while(patch < npatches){
    GEM_DomainBoundary &fp = gp._db[patch++];
    unsigned int ntri,nquad;
    Inf >> ntri >> nquad;
    //    assert(ntri == fp._ntri && nquad == fp._nquad);
    if(ntri > 0){
      fp._triconn.resize(3*ntri);
      for(int i = 0;i < 3;i++){
	unsigned int tri = 0;
	while(tri < ntri)
	  Inf >> fp._triconn[(tri++)*3+i];
      }
    }
    if(nquad > 0){
      fp._quadconn.resize(4*nquad);
      for(int i = 0;i < 4;i++){
	unsigned int quad = 0;
	while(quad < nquad)
	  Inf >> fp._quadconn[(quad++)*4+i];
      }
    }
  }
  Inf.close();
  return true;
}

void
TRAIL_FluCopyCaseFiles(GEM_Partition &gp,const std::string &path)
{
  std::string casename(gp._solver_data._string_data[0]);
  std::string casepath(gp._solver_data._string_data[1]);
  std::ifstream Inf;
  std::ofstream Ouf;
  std::string filename;
  filename.assign(casename+".bc");
  Inf.open((casepath+"/"+filename).c_str());
  Ouf.open((path+"/"+filename).c_str());
  Ouf << Inf.rdbuf();
  Ouf.close();
  Inf.close();
  Ouf.clear();
  Inf.clear();
  filename.assign(casename+".cgi");
  Inf.open((casepath+"/"+filename).c_str());
  Ouf.open((path+"/"+filename).c_str());
  Ouf << Inf.rdbuf();
  Ouf.close();
  Inf.close();
  Ouf.clear();
  Inf.clear();
  filename.assign(casename+".inp");
  Inf.open((casepath+"/"+filename).c_str());
  Ouf.open((path+"/"+filename).c_str());
  Ouf << Inf.rdbuf();
  Ouf.close();
  Inf.close();
  Ouf.clear();
  Inf.clear();
  filename.assign(casename+".map");
  Inf.open((casepath+"/"+filename).c_str());
  Ouf.open((path+"/"+filename).c_str());
  Ouf << Inf.rdbuf();
  Ouf.close();
  Inf.close();
}

// bool
// ReadFluSolnASCII(const string &prefix,unsigned int n,
// 		 double t,bool unsteady,GEM_Partition &gp)
// {
//   string pre(prefix);
//   if(pre.empty()){
// //     if(_casename.empty())
// //       pre = "./";
// //     else
// //       pre = _casepath + "/" + _casename;
//     pre = "rocflu";
//   }
//   ostringstream Ostr;
//   Ostr << pre << ".mixt.cva_" << setw(5) << setfill('0') << gp._id << "_";
//   string filebase(Ostr.str());
//   Ostr.str("");
//   Ostr << setw(6) << setfill('0') << n;
//   string iterstring(Ostr.str());
//   Ostr.clear();
//   Ostr.str("");
//   Ostr << scientific << setprecision(5) << t;
//   string timestring(Ostr.str());
//   timestring.replace(7,1,"E");
//   ifstream Inf;
//   if(unsteady) filebase+=timestring;
//   else filebase+=iterstring;
//   Inf.open(filebase.c_str());
//   if(!Inf){
//     if(gp._debug && gp._out)
//       *gp._out << "FluRegion::ReadFluSolnASCII: Unable to open " << filebase 
// 	       << " for reading.\n";
//     return(false);
//   }
//   string line;
//   SkipLines(Inf,4);
//   Inf >> _soln._resid;
//   SkipLines(Inf,2);
//   Inf >> _soln._current_time;
//   SkipLines(Inf,2);
//   unsigned int ncells = 0;
//   Inf >> ncells;
//   assert(ncells == nelem());
//   SkipLines(Inf,2);
//   unsigned int cell = 0;
//   _soln._rhof.resize(ncells,0.0);
//   _soln._rhovf.resize(3*ncells,0.0);
//   _soln._rhoEf.resize(ncells,0.0);
//   while(cell < ncells)
//     Inf >> _soln._rhof[cell++];
//   SkipLines(Inf,2);  
//   for(int i = 0;i < 3;i++){
//     cell = 0;
//     while(cell < ncells)
//       Inf >> _soln._rhovf[(cell++)*3+i];
//     SkipLines(Inf,2);
//   }
//   cell = 0;
//   while(cell < ncells)
//     Inf >> _soln._rhoEf[cell++];
//   if(!ReadGSPASCII(pre,n,t,unsteady)){
//     if(_debug && _out)
//       *_out 
// 	<< "FluRegion::ReadFluSolnASCII: Unable to process grid speed file."
// 	<< "\n";
//     return(false);
//   }
//   Inf.close();
//   return(true);
// }
// bool
// FluRegion::WriteFluSolnASCII(const string &prefix,unsigned int n,
// 			     double t,bool unsteady)
// {
//   if(_debug && _out)
//     *_out << "FluRegion::WriteFluSolnASCII: Enter\n";
//   string pre(prefix);
//   if(pre.empty()){
//     if(_casename.empty())
//       pre = "./";
//     else
//       pre = _casepath + "/" + _casename + ".";
//   }
//   else
//     pre += ".";
//   ostringstream Ostr;
//   Ostr << pre << "mixt.cva_" << setw(5) << setfill('0') << _id << "_";
//   string filebase(Ostr.str());
//   Ostr.str("");
//   Ostr << setw(6) << setfill('0') << n;
//   string iterstring(Ostr.str());
//   Ostr.clear();
//   Ostr.str("");
//   Ostr << scientific << setprecision(5) << t;
//   string timestring(Ostr.str());
//   timestring.replace(7,1,"E");
//   ofstream Ouf;
//   if(!unsteady)
//     filebase+=iterstring;
//   else
//     filebase+=timestring;
//   Ouf.open(filebase.c_str());
//   if(!Ouf)
//     return false;
//   unsigned int ncells = nelem();
//   unsigned int n_cvar = 5;
//   Ouf << "# ROCFLU flow file" << endl
//       << "# Precision and range" << endl
//       << setw(8) << 15 << setw(8) << 307 << endl
//       << "# Initial residual" << endl
//       << setw(23) << scientific << setprecision(16)
//       << _soln._resid << endl
//       << "# Physical time" << endl
//       << setw(23) << scientific << setprecision(16)
//       << _soln._current_time << endl
//       << "# Dimensions" << endl
//       << setw(8) << ncells << setw(8) << n_cvar << endl
//       << "# Mixture density" << endl;
//   unsigned int cell = 0;
//   unsigned int ll = 5;
//   while(cell < ncells){
//     Ouf <<  setw(23) << scientific << setprecision(16)
// 	<< _soln._rhof[cell++];
//     if(!(cell%ll))
//       Ouf << endl;
//   }
//   if(cell%ll)
//     Ouf << endl;
//   for(int i = 0;i < 3;i++){
//     Ouf << "# Mixture " 
// 	<< (i == 0 ? "x-momentum" :
// 	    i == 1 ? "y-momentum" : "z-momentum") 
// 	<< endl; 
//     cell = 0;
//     while(cell < ncells){
//       Ouf << setw(23) << scientific << setprecision(16) 
// 	  << _soln._rhovf[(cell++)*3+i];
//       if(!(cell%ll))
// 	Ouf << endl;
//     }
//     if(cell%ll)
//       Ouf << endl;
//   }
//   Ouf << "# Mixture total internal energy" << endl;
//   cell = 0;
//   while(cell < ncells){
//     Ouf << setw(23) << scientific << setprecision(16) 
// 	<< _soln._rhoEf[cell++];
//     if(!(cell%ll))
//       Ouf << endl;
//   }
//   if(cell%ll)
//     Ouf << endl;
//   Ouf << "# End" << endl;
//   Ouf.close();
//   WriteGSPASCII(pre,n,t);
//   if(_debug && _out)
//     *_out << "FluRegion::WriteFluSolnASCII: Exit\n";
//   return(true);
// }

// bool
// FluRegion::WriteGSPASCII(const string &prefix,unsigned int n,
// 			 double t,bool unsteady)
// {
//   if(_debug && _out)
//     *_out << "FluRegion::WriteGSPASCII: Enter\n";
//   string pre(prefix);
//   if(pre.empty()){
//     if(_casename.empty())
//       pre = "./";
//     else
//       pre = _casepath + "/" + _casename + ".";
//   }
//   else
//     pre += ".";
//   ostringstream Ostr;
//   Ostr << pre << "gspa_" << setw(5) << setfill('0') << _id << "_";
//   string filebase(Ostr.str());
//   Ostr.str("");
//   Ostr << setw(6) << setfill('0') << n;
//   string iterstring(Ostr.str());
//   Ostr.clear();
//   Ostr.str("");
//   Ostr << scientific << setprecision(5) << t;
//   string timestring(Ostr.str());
//   timestring.replace(7,1,"E");
//   ofstream Ouf;
//   if(!unsteady) filebase+=iterstring;
//   else filebase+=timestring;
//   Ouf.open(filebase.c_str());
//   if(!Ouf)
//     return false;
//   Ouf << "# ROCFLU grid speeds file" << endl
//       << "# Precision and range" << endl
//       << setw(8) << 15 << setw(8) << 307 << endl
//       << "# Physical time" << endl
//       << setw(23) << scientific << setprecision(16) 
//       << t << endl
//       << "# Dimensions" << endl
//       << setw(8) << _nvface << endl
//       << "# Grid speeds" << endl;
//   unsigned int face = 0;
//   unsigned int ll = 5;
//   bool gsp = (_soln._gsp.size() != 0);
//   if(_debug && _out){
//     *_out << "FluRegion::WriteGSPASCII: nvface = " << _nvface << "\n";
//     if(!gsp) 
//       *_out 
// 	<< "FluRegion::WriteGSPASCII:no gridspeeds for the volume faces.\n";
//   }
//   while(face < _nvface){
//     Ouf << setw(23) << scientific << setprecision(16)
// 	<< (gsp ? _soln._gsp[face] : 0);
//     face++;
//     if(!(face%ll))
//       Ouf << endl;
//   }
//   if(face%ll)
//     Ouf << endl;
//   if(_debug && _out)
//     *_out << "FluRegion::WriteGSPASCII: npatches = " << _patches.size() 
// 	  << "\n";
//   vector<FluPatch>::iterator fpi = _patches.begin();
//   while(fpi != _patches.end()){
//     FluPatch &fp = *fpi++;
//     gsp = (fp._soln._gsp.size() != 0);
//     unsigned int nfaces = fp._ntri+fp._nquad-(fp._ngtri+fp._ngquad);
//     if(_debug && _out){
//     *_out << "FluRegion::WriteGSPASCII: npatch_faces = " << nfaces << "\n";
//       if(!gsp) 
// 	*_out << "FluRegion::WriteGSPASCII: no gridspeeds for patch faces." 
// 	     << "\n";
//     }
//     face = 0;
//     while(face < nfaces){
//       Ouf << setw(23) << scientific << setprecision(16)
// 	  << (gsp ? fp._soln._gsp[face] : 0);
//       face++;
//       if(!(face%ll))
// 	Ouf << endl;
//     }
//     if(face%ll)
//       Ouf << endl;
//   }
//   Ouf << "# End" << endl;
//   Ouf.close();  
//   if(_debug && _out)
//     *_out << "FluRegion::WriteGSPASCII: Exit\n";
//   return true;
// }

// bool
// FluRegion::ReadGSPASCII(const string &prefix,unsigned int n,
// 			double t,bool unsteady)
// {
//   if(_debug && _out)
//     *_out << "FluRegion::ReadGSPASCII: Enter\n";
//   string pre(prefix);
//   if(pre.empty()){
//     if(_casename.empty())
//       pre = "./";
//     else
//       pre = _casepath + "/" + _casename;
//   }
//   ostringstream Ostr;
//   Ostr << pre << ".gspa_" << setw(5) << setfill('0') << _id << "_";
//   string filebase(Ostr.str());
//   Ostr.str("");
//   Ostr << setw(6) << setfill('0') << n;
//   string iterstring(Ostr.str());
//   Ostr.clear();
//   Ostr.str("");
//   Ostr << scientific << setprecision(5) << t;
//   string timestring(Ostr.str());
//   timestring.replace(7,1,"E");
//   ifstream Inf;
//   if(!unsteady)
//     Inf.open((filebase+iterstring).c_str());
//   else
//     Inf.open((filebase+timestring).c_str());
//   if(!Inf)
//     return false;
//   string line;
//   SkipLines(Inf,6);
//   Inf >> _nvface;
//   _soln._gsp.resize(_nvface,0.0);
//   SkipLines(Inf,2);
//   unsigned int face = 0;
//   while(face < _nvface)
//     Inf >> _soln._gsp[face++];
//   vector<FluPatch>::iterator fpi = _patches.begin();
//   while(fpi != _patches.end()){
//     FluPatch &fp = *fpi++;
//     unsigned int nfaces = fp._ntri+fp._nquad-(fp._ngtri+fp._ngquad);
//     fp._soln._gsp.resize(nfaces,0.0);
//     face = 0;
//     while(face < nfaces)
//       Inf >> fp._soln._gsp[face++];
//     fp.InitSurfaceSoln(pre,false);
//   }
//   Inf.close();
//   return true;
// }
// bool
// FluRegion::ReadRegionASCII(const string &prefix,unsigned int region_id,
// 			   unsigned int n,double t,bool unsteady)
// {
//   if(_debug && _out)
//     *_out << "FluRegion::ReadRegionASCII:: Reading Rocflu Region "
// 	 << "from Rocflu native files...\n";
//   string pre(prefix);
//   if(pre.empty()){
//     if(_casename.empty())
//       pre = "./";
//     else
//       pre = _casepath + "/" + _casename;
//   }
//   _id = region_id;
//   if(!ReadFluDIM(pre,t,unsteady))
//     return(false);
//   if(!ReadFluCOM(pre)){
//     if(_debug && _out)
//       *_out << "FluRegion::ReadRegionASCII: ReadFluCOM failed.\n";
//     return(false);
//   }
//   if(!ReadFluGridASCII(pre,t,unsteady)){
//     if(_debug && _out)
//       *_out << "FluRegion::ReadRegionASCII: ReadFluGridASCII failed.\n";
//     return(false);
//   }
//   if(!ReadFluSolnASCII(pre,n,t,unsteady)){
//     if(_debug && _out)
//       *_out << "FluRegion::ReadRegionASCII: ReadFluSolnASCII failed.\n";
//     return(false);
//   }
//   if(_debug && _out)
//     *_out << "FluRegion::ReadRegionASCII:: Done reading Rocflu Region "
// 	 << "from Rocflu native files.\n";
//   return(true);
// }










