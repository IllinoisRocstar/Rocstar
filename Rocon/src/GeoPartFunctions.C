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
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <map>
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

bool
TRAIL_GeoPartReadPatchASCII(GEM_DomainBoundary &db,std::ifstream &Inf)
{
  unsigned int nfaces  = 0;
  unsigned int ngfaces = 0;
  db._ngtri = 0;
  db._ngquad = 0;
  Inf >> db._id >> nfaces >> ngfaces;
  unsigned int nfaces_accum = 0;
  nfaces_accum += nfaces;
  db._ngtri = ngfaces;
  if(db._id < 0)
    db._id *= -1;
  db._triconn.resize(3*nfaces);
  unsigned int face = 0;
  while(face < nfaces){
    unsigned int node1,node2,node3;
    Inf >> node1 >> node2 >> node3;
    db._triconn[face*3]   = node1;
    db._triconn[face*3+1] = node2;
    db._triconn[face*3+2] = node3;
    face++;
  }
  Inf >> nfaces >> ngfaces;
  nfaces_accum += nfaces;
  db._ngquad = ngfaces;
  db._quadconn.resize(4*nfaces);
  face = 0;
  while(face < nfaces){
    unsigned int node1,node2,node3,node4;
    Inf >> node1 >> node2 >> node3 >> node4;
    db._quadconn[face*4]   = (node1);
    db._quadconn[face*4+1] = (node2);
    db._quadconn[face*4+2] = (node3);
    db._quadconn[face*4+3] = (node4);
    face++;
  }
  return(true);
}

bool
TRAIL_GeoPartReadPatchesASCII(GEM_Partition &gp,std::ifstream &Inf)
{
  //  map<unsigned int, unsigned int> patch_mapping;
  //  if(!build_patch_mapping(pre,patch_mapping))
  //    return false;
  int nlocal_patches;
  Inf >> nlocal_patches;
  gp._db.resize(nlocal_patches);
  int bound = 0;
  while(bound < nlocal_patches)
    if(!TRAIL_GeoPartReadPatchASCII(gp._db[bound++],Inf))
      return(false);
  return(true);
}

bool 
TRAIL_GeoPartReadPartitionBoundaryASCII(GEM_Partition &gp,
				       unsigned int pb_index,
				       std::ifstream &Inf)
{

  GEM_PartitionBoundary &pb = gp._pb[pb_index];
  // Just read the nodes and elements from the list
  int nnodes  = 0;
  int ngnodes = 0;
  int nshared = 0;
  int nsend   = 0;
  int nreal_nodes = 0;
  int rbid = 0;
  int dummy = 0;
  Inf >> pb._rpart >> gp._solver_data._int_data[0][pb_index]
      >> nnodes >> nshared >> ngnodes;
  nreal_nodes = nnodes - ngnodes;
  nsend = nreal_nodes - nshared;
  pb._rpart += 1; // This converts to the required 1-based id's
  pb._sendnodes.resize(nsend);
  pb._sharenodes.resize(nshared);
  pb._recvnodes.resize(ngnodes);

  int node = 0;
  while(node < nsend)
    Inf >> pb._sendnodes[node++] >> dummy;
  node = 0;
  while(node < nshared)
    Inf >> pb._sharenodes[node++] >> dummy;
  node = 0;
  while(node < ngnodes)
    Inf >> pb._recvnodes[node++] >> dummy;


  unsigned int ncells_total = 0;
  unsigned int nreal_cells = 0;
  unsigned int nghost_cells = 0;

  unsigned int ntets = 0;
  unsigned int ngtets = 0;
  unsigned int nreal_tets = 0;
  Inf >> ntets >> ngtets;
  std::vector<unsigned int> btets;
  btets.resize(ntets + ngtets);
  unsigned int elem = 0;
  while(elem < ntets)
    Inf >> btets[elem++] >> dummy;
  ncells_total += ntets;
  nghost_cells += ngtets;
  nreal_tets = ntets - ngtets;

  unsigned int npyrs = 0;
  unsigned int ngpyrs = 0;
  unsigned int nreal_pyrs = 0;
  std::vector< unsigned int > bpyrs;
  Inf >> npyrs >> ngpyrs;
  bpyrs.resize(npyrs);
  elem = 0;
  while(elem < npyrs)
    Inf >> bpyrs[elem++] >> dummy;
  ncells_total += npyrs;
  nghost_cells += ngpyrs;
  nreal_pyrs = npyrs - ngpyrs;

  unsigned int npris = 0;
  unsigned int ngpris = 0;
  unsigned int nreal_pris = 0;
  std::vector< unsigned int > bpris;
  Inf >> npris >> ngpris;
  bpris.resize(npris);
  elem = 0;
  while(elem < npris)
    Inf >> bpris[elem++] >> dummy;
  ncells_total += npris;
  nghost_cells += ngpris;
  nreal_pris = npris - ngpris;

  unsigned int nhexs = 0;
  unsigned int nghexs = 0;
  unsigned int nreal_hexs = 0;
  std::vector< unsigned int > bhexs;
  Inf >> nhexs >> nghexs;
  bhexs.resize(nhexs);
  elem = 0;
  while(elem < nhexs)
    Inf >> bhexs[elem++] >> dummy;
  ncells_total += nhexs;
  nghost_cells += nghexs;
  nreal_hexs = nhexs - nghexs;

  nreal_cells = ncells_total - nghost_cells;
  pb._recvcells.resize(nghost_cells);
  pb._sendcells.resize(nreal_cells);

  elem = 0;
  unsigned int send_cell = 0;
  unsigned int recv_cell = 0;
  while(elem < ntets)
    if(elem < nreal_tets)
      pb._sendcells[send_cell++] = gp.Elem2Cell(std::make_pair(1,btets[elem++]));
    else
      pb._recvcells[recv_cell++] = gp.Elem2Cell(std::make_pair(1,btets[elem++]));
  elem = 0;
  while(elem < npyrs)
    if(elem < nreal_pyrs)
      pb._sendcells[send_cell++] = gp.Elem2Cell(std::make_pair(2,bpyrs[elem++]));
    else
      pb._recvcells[recv_cell++] = gp.Elem2Cell(std::make_pair(2,bpyrs[elem++]));
  elem = 0;
  while(elem < npris)
    if(elem < nreal_pris)
      pb._sendcells[send_cell++] = gp.Elem2Cell(std::make_pair(3,bpris[elem++]));
    else
      pb._recvcells[recv_cell++] = gp.Elem2Cell(std::make_pair(3,bpris[elem++]));
  elem = 0;
  while(elem < nhexs)
    if(elem < nreal_hexs)
      pb._sendcells[send_cell++] = gp.Elem2Cell(std::make_pair(4,bhexs[elem++]));
    else
      pb._recvcells[recv_cell++] = gp.Elem2Cell(std::make_pair(4,bhexs[elem++]));

  assert(send_cell == nreal_cells && recv_cell == nghost_cells);
  return(true);
}

bool 
TRAIL_GeoPartReadPartitionBoundariesASCII(GEM_Partition &gp,std::ifstream &Inf)
{
  int nbound = gp._pb.size();
  int bound = 0;
  while(bound < nbound)
    if(!TRAIL_GeoPartReadPartitionBoundaryASCII(gp,bound++,Inf))
      return(false);
  return(true);
}

bool 
TRAIL_GeoPartReadASCII(GEM_Partition &gp,std::ifstream &Inf)
{

  unsigned int nvfaces,nbfaces,nnodes;
  Inf >> nvfaces >> nbfaces >> nnodes >> gp._ngnodes;
  gp._nvface =  nvfaces + nbfaces;
  unsigned int node = 0;
  gp._nc.resize(3*nnodes);
  while(node < nnodes){
    Inf >> gp._nc[3*node]   >> gp._nc[3*node+1] 
	>> gp._nc[3*node+2];
    node++;
  }
  unsigned int ntet = 0;
  Inf >> ntet >> gp._ngtet;
  unsigned int tet = 0;
  gp._tetconn.resize(4*ntet);
  while(tet < ntet){
    Inf >> gp._tetconn[4*tet]   >> gp._tetconn[4*tet+1] 
	>> gp._tetconn[4*tet+2] >> gp._tetconn[4*tet+3];
    tet++;
  }
  unsigned int npyr = 0;
  Inf >> npyr >> gp._ngpyr;
  unsigned int pyr = 0;
  gp._pyrconn.resize(5*npyr);
  while(pyr < npyr){
    Inf >> gp._pyrconn[5*pyr]   >> gp._pyrconn[5*pyr+1] 
	>> gp._pyrconn[5*pyr+2] >> gp._pyrconn[5*pyr+3] 
	>> gp._pyrconn[5*pyr+4];
    pyr++;
  }

  unsigned int npris = 0;
  Inf >> npris >> gp._ngpris;
  unsigned int pris = 0;
  gp._prisconn.resize(6*npris);
  while(pris < npris){
    Inf >> gp._prisconn[6*pris]   >> gp._prisconn[6*pris+1] 
	>> gp._prisconn[6*pris+2] >> gp._prisconn[6*pris+3] 
	>> gp._prisconn[6*pris+4] >> gp._prisconn[6*pris+5];
    pris++;
  }

  unsigned int nhex = 0;
  Inf >> nhex >> gp._nghex;
  unsigned int hex = 0;
  gp._hexconn.resize(8*nhex);
  while(hex < nhex){
    Inf >> gp._hexconn[8*hex]   >> gp._hexconn[8*hex+1] 
	>> gp._hexconn[8*hex+2] >> gp._hexconn[8*hex+3] 
	>> gp._hexconn[8*hex+4] >> gp._hexconn[8*hex+5] 
	>> gp._hexconn[8*hex+6] >> gp._hexconn[8*hex+7];
    hex++;
  }

  //  _nelem = _nhex + _ntet + _npyr + _npris;
  unsigned int nbound = 0;
  Inf >> nbound;
  gp._pb.resize(nbound);
  // Note: the gp._cell_ordering[] array should already be 
  // set up properly before these routines are called.  To
  // stay general, we require that this is done outside of 
  // this function
  if(!TRAIL_GeoPartReadPartitionBoundariesASCII(gp,Inf))
    return(false);
  if(!TRAIL_GeoPartReadPatchesASCII(gp,Inf))
    return(false);
  Inf.close();
  return(true);
}







