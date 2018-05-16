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

using namespace std;

#include "GEM.H"
#include "Partition.H"


bool
PartitionPatch::ReadPartitionPatch(ifstream &Inf)
{
  int buf[3];
  Inf.read(reinterpret_cast<char *>(buf),12);
  _id      = -buf[0];
  _nfaces  = buf[1];
  _ngfaces = buf[2];
  //  _faceconn.resize(4*_nfaces);
  //  Inf.read(reinterpret_cast<char *>(&_faceconn[0]),4 * nfaces);
  // skip the face connectivity until Bill tell us whats up with it
  Inf.seekg(16 * _nfaces,ios::cur);
  return(true);
}
  
bool
PartitionPatch::ReadPartitionPatchASCII(ifstream &Inf)
{
  _ngtri = 0;
  _ngquad = 0;
  Inf >> _id >> _nfaces >> _ngfaces;
  unsigned int nfaces_accum = 0;
  nfaces_accum += _nfaces;
  _ngtri = _ngfaces;
  _id *= -1;
  _triconn.resize(3*_nfaces);
  unsigned int face = 0;
  while(face < _nfaces){
    unsigned int node1,node2,node3;
    Inf >> node1 >> node2 >> node3;
    _triconn[face*3]   = node1;
    _triconn[face*3+1] = node2;
    _triconn[face*3+2] = node3;
    face++;
  }
  Inf >> _nfaces >> _ngfaces;
  nfaces_accum += _nfaces;
  _ngquad = _ngfaces;
  _quadconn.resize(4*_nfaces);
  face = 0;
  while(face < _nfaces){
    unsigned int node1,node2,node3,node4;
    Inf >> node1 >> node2 >> node3 >> node4;
    _quadconn[face*4]   = (node1);
    _quadconn[face*4+1] = (node2);
    _quadconn[face*4+2] = (node3);
    _quadconn[face*4+3] = (node4);
    face++;
  }
  _nfaces = nfaces_accum;
  _ngfaces = _ngquad + _ngtri;
  return(true);
}
  
bool
PartitionPatch::WritePartitionPatchASCII(ofstream &Ouf)
{
  unsigned int nfaces = _triconn.size()/3;
  unsigned int face = 0;
  Ouf << setw(12) << (_id < 0 ? _id : -1*_id) << endl 
      << setw(12) << nfaces << endl
      << setw(12) << _ngtri << endl;
  while(face < nfaces){
    Ouf << setw(12) << _triconn[face*3]
	<< setw(12) << _triconn[face*3+1]
	<< setw(12) << _triconn[face*3+2] << endl;
    face++;
  }
  nfaces = _quadconn.size()/4;
  face = 0;
  Ouf << setw(12) << nfaces << endl
      << setw(12) << _ngquad << endl;
  while(face < nfaces){
    Ouf << setw(12) << _quadconn[face*4]
	<< setw(12) << _quadconn[face*4+1]
	<< setw(12) << _quadconn[face*4+2]
	<< setw(12) << _quadconn[face*4+3] << endl;
    face++;
  }
  return(true);
}
  

bool
Partition::ReadPartitionPatches(ifstream &Inf)
{
  //  map<unsigned int, unsigned int> patch_mapping;
  //  if(!build_patch_mapping(pre,patch_mapping))
  //    return false;
  int nlocal_patches;
  Inf.read(reinterpret_cast<char *>(&nlocal_patches),4);
  _patches.resize(nlocal_patches);
  int bound = 0;
  while(bound < nlocal_patches){
    PartitionPatch *pp = &_patches[bound++];
    if(!pp->ReadPartitionPatch(Inf))
      return(false);
  }
  return(true);
}

bool
Partition::ReadPartitionPatchesASCII(ifstream &Inf)
{
  int nlocal_patches;
  Inf >> nlocal_patches;
  _patches.resize(nlocal_patches);
  int bound = 0;
  while(bound < nlocal_patches){
    PartitionPatch *pp = &_patches[bound++];
    if(!pp->ReadPartitionPatchASCII(Inf))
      return(false);
  }
  return(true);
}


bool
Partition::WritePartitionPatchesASCII(ofstream &Ouf)
{
  Ouf << setw(12) << _patches.size() << endl;
  unsigned int bound = 0;
  while(bound < _patches.size()){
    PartitionPatch *pp = &_patches[bound++];
    if(!pp->WritePartitionPatchASCII(Ouf))
      return(false);
  }
  return(true);
}
  
bool
PartitionBoundary::ReadPartitionBoundary(ifstream &Inf)
{
    int buf[4];
    Inf.read(reinterpret_cast<char *>(buf),16);
    _rpart      = buf[0];
    _rbid       = buf[1];
    int nnodes  = buf[2];
    _ngnodes    = buf[3];

    _bnodes.resize(nnodes);
    int *nbuf = new int [2 * nnodes];
    Inf.read(reinterpret_cast<char *>(nbuf),8*nnodes);

    int node = 0;
    while(node < nnodes){
      _bnodes[node] = make_pair(nbuf[2*node],nbuf[2*node+1]);
      node++;
    }
    delete [] nbuf;

    Inf.read(reinterpret_cast<char *>(buf),8);
    int nelem = buf[0];
    _ngtets   = buf[1];
    _btets.resize(nelem);
    nbuf = new int [2*nelem];
    Inf.read(reinterpret_cast<char *>(nbuf),nelem*8);
    int elem = 0;
    while(elem < nelem){
      _btets[elem] = make_pair(nbuf[2*elem],nbuf[2*elem+1]);
      elem++;
    }
    delete [] nbuf;

    Inf.read(reinterpret_cast<char *>(buf),8);
    nelem = buf[0];
    _ngpyr   = buf[1];
    _bpyr.resize(nelem);
    nbuf = new int [2*nelem];
    Inf.read(reinterpret_cast<char *>(nbuf),nelem*8);
    elem = 0;
    while(elem < nelem){
      _bpyr[elem] = make_pair(nbuf[2*elem],nbuf[2*elem+1]);
      elem++;
    }
    delete [] nbuf;

    Inf.read(reinterpret_cast<char *>(buf),8);
    nelem = buf[0];
    _ngpris   = buf[1];
    _bpris.resize(nelem);
    nbuf = new int [2*nelem];
    Inf.read(reinterpret_cast<char *>(nbuf),nelem*8);
    elem = 0;
    while(elem < nelem){
      _bpris[elem] = make_pair(nbuf[2*elem],nbuf[2*elem+1]);
      elem++;
    }
    delete [] nbuf;

    Inf.read(reinterpret_cast<char *>(buf),8);
    nelem = buf[0];
    _nghex   = buf[1];
    _bhex.resize(nelem);
    nbuf = new int [2*nelem];
    Inf.read(reinterpret_cast<char *>(nbuf),nelem*8);
    elem = 0;
    while(elem < nelem){
      _bhex[elem] = make_pair(nbuf[2*elem],nbuf[2*elem+1]);
      elem++;
    }
    delete [] nbuf;
    
    return(true);
}

bool 
PartitionBoundary::ReadPartitionBoundaryASCII(ifstream &Inf)
{
  int nnodes  = 0;
  Inf >> _rpart >> _rbid >> nnodes >> _ngnodes;
  _rpart += 1; // This converts to the required 1-based id's
  _bnodes.resize(nnodes);
  int node = 0;
  int l,r;
  while(node < nnodes){
    Inf >> l >> r;
    _bnodes[node++] = make_pair(l,r);
  }

  int nelem = 0;
  Inf >> nelem >> _ngtets;
  _btets.resize(nelem);
  int elem = 0;
  while(elem < nelem){
    Inf >> l >> r;
    _btets[elem++] = make_pair(l,r);
  }
  
  nelem = 0;
  Inf >> nelem >> _ngpyr;
  _bpyr.resize(nelem);
  elem = 0;
  while(elem < nelem){
    Inf >> l >> r;
    _bpyr[elem++] = make_pair(l,r);
  }
  
  nelem = 0;
  Inf >> nelem >> _ngpris;
  _bpris.resize(nelem);
  elem = 0;
  while(elem < nelem){
    Inf >> l >> r;
    _bpris[elem++] = make_pair(l,r);
  }
  
  nelem = 0;
  Inf >> nelem >> _nghex;
  _bhex.resize(nelem);
  elem = 0;
  while(elem < nelem){
    Inf >> l >> r;
    _bhex[elem++] = make_pair(l,r);
  }
  
  return(true);
}

bool 
PartitionBoundary::WritePartitionBoundaryASCII(ofstream &Ouf)
{
  unsigned int nelem = _bnodes.size();
  Ouf << setw(12) << _rpart-1 << endl
      << setw(12) << _rbid << endl
      << setw(12) << nelem << endl 
      << setw(12) << _ngnodes << endl;
  unsigned int node = 0;
  while(node < nelem){
    Ouf << setw(12) << _bnodes[node].first 
	<< setw(12) << _bnodes[node].second
	<< endl;
    node++;
  }
  nelem = _btets.size();
  Ouf << setw(12) << nelem << endl
      << setw(12) << _ngtets << endl;
  unsigned int elem = 0;
  while(elem < nelem){
    Ouf << setw(12) << _btets[elem].first
	<< setw(12) << _btets[elem].second
	<< endl;
    elem++;
  }
  
  nelem = _bpyr.size();
  Ouf <<  setw(12) << nelem << endl
      << setw(12) << _ngpyr << endl;
  elem = 0;
  while(elem < nelem){
    Ouf << setw(12) << _bpyr[elem].first
	<< setw(12) << _bpyr[elem].second
	<< endl;
    elem++;
  }
  nelem = _bpris.size();
  Ouf << setw(12) << nelem << endl
      << setw(12) << _ngpris << endl;
  elem = 0;
  while(elem < nelem){
    Ouf << setw(12) << _bpris[elem].first
	<< setw(12) << _bpris[elem].second
	<< endl;
    elem++;
  }
  nelem = _bhex.size();
  Ouf << setw(12) << nelem << endl
      << setw(12) << _nghex << endl;
  elem = 0;
  while(elem < nelem){
    Ouf << setw(12) << _bhex[elem].first
	<< setw(12) << _bhex[elem].second
	<< endl;
    elem++;
  }
  return(true);
}

bool 
Partition::ReadPartitionBoundariesASCII(ifstream &Inf)
{
  int nbound = _boundaries.size();
  int bound = 0;
  while(bound < nbound){
    PartitionBoundary *pb = &_boundaries[bound++];
    if(!pb->ReadPartitionBoundaryASCII(Inf))
      return(false);
    //    pb->MyRank(_rank);
  }
  return(true);
}

bool 
Partition::WritePartitionBoundariesASCII(ofstream &Ouf)
{
  int nbound = _boundaries.size();
  int bound = 0;
  while(bound < nbound){
    PartitionBoundary *pb = &_boundaries[bound++];
    if(!pb->WritePartitionBoundaryASCII(Ouf))
      return(false);
  }
  return(true);
}

bool 
Partition::ReadPartitionBoundaries(ifstream &Inf)
{
  int nbound = _boundaries.size();
  int bound = 0;
  while(bound < nbound){
    PartitionBoundary *pb = &_boundaries[bound++];
    if(!pb->ReadPartitionBoundary(Inf))
      return(false);
    pb->MyRank(_id);
  }
  return(true);
}


bool 
Partition::ReadPartition(const string &pre,unsigned int rank)
{
  ifstream Inf;
  ostringstream Ostr;

  Ostr << pre << "." << rank;
  Inf.open(Ostr.str().c_str());
  if(!Inf) return(false);
  
  long buf[2];
  this->_id = rank;
  Inf.read(reinterpret_cast<char *>(&buf[0]),8);
  this->_nnodes = buf[0];
  cout << "buf[0] = " << buf[0] << endl;
  this->_ngnodes = buf[1];
  // skip over coordinates 3 of 4 bytes for each node
  Inf.seekg(this->_nnodes * 3 * 4,ios::cur);
  Inf.read(reinterpret_cast<char *>(buf),8);
  this->_ntet = buf[0];
  this->_ngtet = buf[1];
  // skip over connectivity 4 of 4 bytes for each tet
  Inf.seekg(this->_ntet * 4 * 4,ios::cur);
  Inf.read(reinterpret_cast<char *>(buf),8);
  this->_npyr = buf[0];
  this->_ngpyr = buf[1];
  // skip over connectivity 5 of 4 bytes for each pyr
  Inf.seekg(this->_npyr * 5 * 4,ios::cur);
  Inf.read(reinterpret_cast<char *>(buf),8);
  this->_npris = buf[0];
  this->_ngpris = buf[1];
  // skip over connectivity 6 of 4 bytes for each pris
  Inf.seekg(this->_npris * 6 * 4,ios::cur);
  Inf.read(reinterpret_cast<char *>(buf),8);
  this->_nhex = buf[0];
  this->_nghex = buf[1];
  // skip over connectivity 8 of 4 bytes for each hex
  Inf.seekg(this->_nhex * 8 * 4,ios::cur);
  this->_nelem = this->_nhex + this->_ntet + this->_npyr + this->_npris;
  int nbound = 0;
  Inf.read(reinterpret_cast<char *>(&nbound),4);
  _boundaries.resize(nbound);
  if(!ReadPartitionBoundaries(Inf))
    return(false);
  if(!ReadPartitionPatches(Inf))
    return(false);
  Inf.close();
  return(true);
}
  
bool 
Partition::ReadPartitionASCII(const string &pre,unsigned int rank)
{
  ifstream Inf;
  ostringstream Ostr;

  Ostr << pre << "." << rank << ".asc";
  Inf.open(Ostr.str().c_str());
  if(!Inf){
    cerr << "Cannot open " << Ostr.str() << "." << endl;
    return(false);
  }
 
  _id = rank + 1;
  unsigned int nvfaces,nbfaces;
  Inf >> nvfaces >> nbfaces >> _nnodes >> _ngnodes;
  _nvface =  nvfaces + nbfaces;
  unsigned int node = 0;
  _nc.resize(3*_nnodes);
  while(node < _nnodes){
    Inf >> _nc[3*node] >> _nc[3*node+1] >> _nc[3*node+2];
    node++;
  }
  Inf >> _ntet >> _ngtet;
  unsigned int tet = 0;
  _tetconn.resize(4*_ntet);
  while(tet < _ntet){
    Inf >> _tetconn[4*tet] >> _tetconn[4*tet+1] >> _tetconn[4*tet+2]
	>> _tetconn[4*tet+3];
    tet++;
  }
  Inf >> _npyr >> _ngpyr;
  unsigned int pyr = 0;
  _pyrconn.resize(5*_npyr);
  while(pyr < _npyr){
    Inf >> _pyrconn[5*pyr] >> _pyrconn[5*pyr+1] >> _pyrconn[5*pyr+2]
	>> _pyrconn[5*pyr+3] >> _pyrconn[5*pyr+4];
    pyr++;
  }

  Inf >> _npris >> _ngpris;
  unsigned int pris = 0;
  _prisconn.resize(6*_npris);
  while(pris < _npris){
    Inf >> _prisconn[6*pris] >> _prisconn[6*pris+1] >> _prisconn[6*pris+2]
	>> _prisconn[6*pris+3] >> _prisconn[6*pris+4] 
	>> _prisconn[6*pris+5];
    pris++;
  }

  Inf >> _nhex >> _nghex;
  unsigned int hex = 0;
  _hexconn.resize(8*_nhex);
  while(hex < _nhex){
    Inf >> _hexconn[8*hex] >> _hexconn[8*hex+1] >> _hexconn[8*hex+2]
	>> _hexconn[8*hex+3] >> _hexconn[8*hex+4] 
	>> _hexconn[8*hex+5] >> _hexconn[8*hex+6] >> _hexconn[8*hex+7];
    hex++;
  }
  _nelem = _nhex + _ntet + _npyr + _npris;
  unsigned int nbound = 0;
  Inf >> nbound;
  _boundaries.resize(nbound);
  if(!ReadPartitionBoundariesASCII(Inf))
    return(false);
  if(!ReadPartitionPatchesASCII(Inf))
    return(false);
  Inf.close();
  return(true);
}

bool
Partition::BuildPartitionBoundaries(const vector<int> &Pconn)
{
  unsigned int index = 0;
  unsigned int npb = Pconn[index++];
  _boundaries.clear();
  _boundaries.resize(npb);
  
  // Shared Nodes
  unsigned int pb_index = 0;
  while(pb_index < npb){
    PartitionBoundary &pb = _boundaries[pb_index++];
    unsigned int pane_id = Pconn[index++];
    pb._rpart  = pane_id/100;
    //    pb._mypart = _part;
    //    pb._rbid   = pane_id%100;
    unsigned int nshared = Pconn[index++];
    pb._sharenodes.resize(nshared);
    unsigned int node = 0;
    while(node < nshared)
      pb._sharenodes[node++] = Pconn[index++];
  }
  index++;

  // Send Nodes
  pb_index = 0;
  while(pb_index < npb){
    PartitionBoundary &pb = _boundaries[pb_index++];
    unsigned int pane_id = Pconn[index++];
    assert(pb._rpart == pane_id/100);
    unsigned int nsend = Pconn[index++];
    pb._sendnodes.resize(nsend);
    unsigned int node = 0;
    while(node < nsend)
      pb._sendnodes[node++] = Pconn[index++];
  }
  index++;

  // Recv Nodes
  pb_index = 0;
  while(pb_index < npb){
    PartitionBoundary &pb = _boundaries[pb_index++];
    unsigned int pane_id = Pconn[index++];
    assert(pb._rpart == pane_id/100);
    unsigned int nrecv = Pconn[index++];
    pb._recvnodes.resize(nrecv);
    unsigned int node = 0;
    while(node < nrecv)
      pb._recvnodes[node++] = Pconn[index++];
  }
  index++;

  // Send Elems - first, we have to populate a temprorary array
  // then, we have to sort through the array and determine how
  // many of each element type there are.
  pb_index = 0;
  while(pb_index < npb){
    PartitionBoundary &pb = _boundaries[pb_index++];
    unsigned int pane_id = Pconn[index++];
    assert(pb._rpart == pane_id/100);
    unsigned int nsend = Pconn[index++];
    pb._sendcells.resize(nsend);
    unsigned int node = 0;
    while(node < nsend)
      pb._sendcells[node++] = Pconn[index++];
    // Now temp_id contains the mapped element id
    // We need to know whether each element is a 
    // tet, pyr, pris, or hex and keep a count.
    node = 0;
    unsigned int ntet = 0;
    unsigned int npyr = 0;
    unsigned int npris = 0;
    unsigned int nhex = 0;
    while(node < nsend){
      if(pb._sendcells[node] <= _ntet)
	ntet++;
      else if (pb._sendcells[node] <= (_ntet+_npyr))
	npyr++;
      else if (pb._sendcells[node] <= (_ntet+_npyr+_npris))
	npris++;
      else if (pb._sendcells[node] <= (_ntet+_npyr+_npris+_nhex))
	nhex++;
      else
	assert(pb._sendcells[node] <= (_ntet+_npyr+_npris+_nhex));
      node++;
    }
    pb._belem_send[0].resize(ntet);
    pb._belem_send[1].resize(npyr);
    pb._belem_send[2].resize(npris);
    pb._belem_send[3].resize(nhex);
    node = 0;
    ntet = 0;
    npyr = 0;
    npris = 0;
    nhex = 0;
    while(node < nsend){
      if(pb._sendcells[node] <= _ntet)
	pb._belem_send[0][ntet++] = pb._sendcells[node];
      else if (pb._sendcells[node] <= (_ntet+_npyr))
	pb._belem_send[1][npyr++] = pb._sendcells[node] - _ntet;
      else if (pb._sendcells[node] <= (_ntet+_npyr+_npris))
	pb._belem_send[2][npris++] = pb._sendcells[node] - (_ntet+_npyr);
      else 
	pb._belem_send[3][nhex++]  = pb._sendcells[node] - (_ntet+_npyr+_npris);
      node++;
    }
  }
  index++;

  // Recv Elems - first, we have to populate a temprorary array
  // then, we have to sort through the array and determine how
  // many of each element type there are.  Then we populate the 
  // array
  pb_index = 0;
  while(pb_index < npb){
    PartitionBoundary &pb = _boundaries[pb_index++];
    unsigned int pane_id = Pconn[index++];
    assert(pb._rpart == pane_id/100);
    unsigned int nrecv = Pconn[index++];
    pb._recvcells.resize(nrecv);
    unsigned int node = 0;
    while(node < nrecv)
      pb._recvcells[node++] = Pconn[index++];
    // Now pb._recvcells contains the mapped element id
    // We need to know whether each element is a 
    // tet, pyr, pris, or hex and keep a count.
    node = 0;
    unsigned int ntet = 0;
    unsigned int npyr = 0;
    unsigned int npris = 0;
    unsigned int nhex = 0;
    while(node < nrecv){
      if(pb._recvcells[node] <= _ntet)
	ntet++;
      else if (pb._recvcells[node] <= (_ntet+_npyr))
	npyr++;
      else if (pb._recvcells[node] <= (_ntet+_npyr+_npris))
	npris++;
      else if (pb._recvcells[node] <= (_ntet+_npyr+_npris+_nhex))
	nhex++;
      else
	assert(pb._recvcells[node] <= (_ntet+_npyr+_npris+_nhex));
      node++;
    }
    pb._belem_recv[0].resize(ntet);
    pb._belem_recv[1].resize(npyr);
    pb._belem_recv[2].resize(npris);
    pb._belem_recv[3].resize(nhex);
    node = 0;
    ntet = 0;
    npyr = 0;
    npris = 0;
    nhex = 0;
    while(node < nrecv){
      if(pb._recvcells[node] <= _ntet)
	pb._belem_recv[0][ntet++] = pb._recvcells[node];
      else if (pb._recvcells[node] <= (_ntet+_npyr))
	pb._belem_recv[1][npyr++] = pb._recvcells[node] - _ntet;
      else if (pb._recvcells[node] <= (_ntet+_npyr+_npris))
	pb._belem_recv[2][npris++] = pb._recvcells[node] - (_ntet+_npyr);
      else 
	pb._belem_recv[3][nhex++]  = pb._recvcells[node] - (_ntet+_npyr+_npris);
      node++;
    }
  }
  return(true);
}

bool 
Partition::WritePartitionASCII(const string &pre)
{
  ofstream Ouf;
  ostringstream Ostr;

  Ostr << pre << "." << _id-1 << ".asc";
  Ouf.open(Ostr.str().c_str());
  if(!Ouf){
    cerr << "Cannot open " << Ostr.str() << "." << endl;
    return(false);
  }
  

  Ouf << setw(12) << _nvface
      << setw(12) << 0
      << setw(12) << _nnodes 
      << setw(12) << _ngnodes
      << endl;
  unsigned int node = 0;
  while(node < _nnodes){
    Ouf << setw(20) <<  _nc[3*node] 
	<< setw(20) <<  _nc[3*node+1] 
	<< setw(20) <<  _nc[3*node+2]
	<< endl;
    node++;
  }
  Ouf << setw(12) << _ntet 
      << setw(12) << _ngtet
      << endl;
  unsigned int tet = 0;
  while(tet < _ntet){
    Ouf << setw(12) << _tetconn[4*tet] 
	<< setw(12) << _tetconn[4*tet+1]
	<< setw(12) << _tetconn[4*tet+2]
	<< setw(12) <<  _tetconn[4*tet+3]
	<< endl;
    tet++;
  }
  Ouf << setw(12) << _npyr
      << setw(12) << _ngpyr 
      << endl;
  unsigned int pyr = 0;
  while(pyr < _npyr){
    Ouf << setw(12) << _pyrconn[5*pyr] 
	<< setw(12) << _pyrconn[5*pyr+1] 
	<< setw(12) << _pyrconn[5*pyr+2]
	<< setw(12) << _pyrconn[5*pyr+3] 
	<< setw(12) << _pyrconn[5*pyr+4]
	<< endl;
    pyr++;
  }

  Ouf << setw(12) << _npris 
      << setw(12) << _ngpris
      << endl;
  unsigned int pris = 0;
  while(pris < _npris){
    Ouf << setw(12) << _prisconn[6*pris] 
	<< setw(12) << _prisconn[6*pris+1] 
	<< setw(12) << _prisconn[6*pris+2]
	<< setw(12) << _prisconn[6*pris+3] 
	<< setw(12) << _prisconn[6*pris+4] 
	<< setw(12) << _prisconn[6*pris+5]
	<< endl;
    pris++;
  }

  Ouf << setw(12) <<  _nhex 
      << setw(12) <<  _nghex
      << endl;
  unsigned int hex = 0;
  while(hex < _nhex){
    Ouf << setw(12) << _hexconn[8*hex] 
	<< setw(12) << _hexconn[8*hex+1]
	<< setw(12) << _hexconn[8*hex+2]
	<< setw(12) << _hexconn[8*hex+3]
	<< setw(12) << _hexconn[8*hex+4] 
	<< setw(12) << _hexconn[8*hex+5]
	<< setw(12) << _hexconn[8*hex+6] 
	<< setw(12) << _hexconn[8*hex+7]
	<< endl;
    hex++;
  }
  Ouf << _boundaries.size() << endl;
  if(!WritePartitionBoundariesASCII(Ouf))
    return(false);
  if(!WritePartitionPatchesASCII(Ouf))
    return(false);
  Ouf.close();
  return(true);
}

void
report_border_array(const vector<pair<unsigned int,unsigned int> > &v)
{
  vector<pair<unsigned int,unsigned int> >::const_iterator vi = v.begin();
  while(vi != v.end()){
    cerr << vi->first << " " << vi->second << endl;
    vi++;
  }
}

void
PartitionBoundary::populate_local_arrays(const PartitionBoundary &rpb,
					 GEM_Partition &gp)
{
  bool debug = gp._debug;
  if(debug)
    cerr << "PartitionBoundary::populate_local_arrays entry" << endl;
  // Total number of nodes on this boundary (shared+send+recv)
  unsigned int nnodes = _bnodes.size();
  // Number of local reals is local total minus local ghosts
  unsigned int nreal = _bnodes.size() - _ngnodes;
  // Number of remote recv'd nodes
  unsigned int nremote_gnodes = rpb._ngnodes;
  // Number of local shared is local real - remote ghosts
  unsigned int nshared = nreal - nremote_gnodes;
  // Num remote real is remote total minus remote ghosts
  unsigned int nremote_rnodes = rpb._bnodes.size() - nremote_gnodes;
  // Num local send nodes must equal num remote ghosts
  unsigned int nsend = nreal - nshared;
  assert(nsend == nremote_gnodes);
  _recvnodes.resize(_ngnodes);
  _sharenodes.resize(nshared);
  _sendnodes.resize(nsend);
  //
  // We always let the sender determine the node ordering,
  // so the recv arrays need to be re-ordered so that they
  // match the remote border's node ordering.  Also, we
  // let the lesser of the two partition_id's determine the
  // shared node order.  
  //
  // The format for the data in _bnodes[] at this point is:
  //
  // vector<pair<LocalNodeIndex,RemoteBoundaryNodeIndex> >
  // 
  // And we want to populate the send, recv, and shared 
  // node arrays in the correct order.  For those that need
  // the remote ordering, we swap the pairs above and sort
  // the array based on the RemoteBoundaryNodeIndex (which
  // indexes into the remote partition's boundary node 
  // array).  Then we take a second pass to populate the 
  // local indices in that ordering.
  //
  // These hold the inverted local<-->remote node indices
  vector<pair<unsigned int,unsigned int> > recv_iindex;
  vector<pair<unsigned int,unsigned int> > shared_iindex;
  recv_iindex.resize(_ngnodes);
  shared_iindex.resize(nshared);
  unsigned int node = 0;
  unsigned int shnode = 0;
  unsigned int snode = 0;
  unsigned int rnode = 0;
  while(node < nnodes){
    // All ghost nodes are processed here
    if(node >= nreal)
      recv_iindex[rnode++] = make_pair(_bnodes[node].second,
				       _bnodes[node].first);
    //      recv_iindex[rnode++] = _bnodes[node];
    else {
      if(_bnodes[node].second > nremote_rnodes)
	_sendnodes[snode++] = _bnodes[node].first;
      else {
	if(_mypart < _rpart) 
	  _sharenodes[shnode++] = _bnodes[node].first;
	else
	  shared_iindex[shnode++] = make_pair(_bnodes[node].second,
					      _bnodes[node].first);
      }
    }
    node++;
  }
  sort(recv_iindex.begin(),recv_iindex.end());
  node = 0;
  while(node < rnode){
    _recvnodes[node] = recv_iindex[node].second;
    node++;
  }
  // Check partition_id (encoded as _myrank and _rrank) and order shared
  // nodes if needed.
  if(_mypart > _rpart){
    node = 0; 
    sort(shared_iindex.begin(),shared_iindex.end());
    while(node < shnode){
      _sharenodes[node] = shared_iindex[node].second;
      node++;
    }
  }
  if(debug)
    cerr << "PartitionBoundary::populate_local_arrays: Beginning elements."
	 << endl;
  // Nodes are done.  Now we need send and recv element arrays. These will be
  // done independently for each element type.  
  unsigned int ntets = _btets.size();
  unsigned int nhex  = _bhex.size();
  unsigned int npyr  = _bpyr.size();
  unsigned int npris = _bpris.size();
  unsigned int nreal_tets = ntets - _ngtets;
  unsigned int nreal_hex = nhex - _nghex;
  unsigned int nreal_pyr = npyr - _ngpyr;
  unsigned int nreal_pris = npris - _ngpris;
  _belem_send[0].resize(nreal_tets);
  _belem_send[1].resize(nreal_pyr);
  _belem_send[2].resize(nreal_pris);
  _belem_send[3].resize(nreal_hex);
  _belem_recv[0].resize(_ngtets);
  _belem_recv[1].resize(_ngpyr);
  _belem_recv[2].resize(_ngpris);
  _belem_recv[3].resize(_nghex);
  _sendcells.resize(nreal_tets+nreal_pyr+nreal_pris+nreal_hex);
  _recvcells.resize(_ngtets+_ngpyr+_ngpris+_nghex);
  unsigned int elem = 0;
  unsigned int selem = 0;
  unsigned int relem = 0;
  unsigned int scell = 0;
  unsigned int rcell = 0;
  recv_iindex.resize(_ngtets);
  while(elem < ntets){
    if(elem < nreal_tets){
      _belem_send[0][selem++] = _btets[elem].first;
      _sendcells[scell++] = gp.Elem2Cell(make_pair((unsigned int)1,
						   _btets[elem].first));
    }
    else
      recv_iindex[relem++] = make_pair(_btets[elem].second,
				       _btets[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[0][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = gp.Elem2Cell(make_pair((unsigned int)1,
						 recv_iindex[elem].second));
    elem++;
  }
  assert((relem == _ngtets) && (selem == nreal_tets)); 
  elem = 0;
  selem = 0;
  relem = 0;
  recv_iindex.resize(_ngpyr);
  while(elem < npyr){
    if(elem < nreal_pyr){
      _belem_send[1][selem++] = _bpyr[elem].first;
      _sendcells[scell++] = gp.Elem2Cell(make_pair((unsigned int)2,
						   _bpyr[elem].first));
    }
    else
      recv_iindex[relem++] = make_pair(_bpyr[elem].second,
				       _bpyr[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[1][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = gp.Elem2Cell(make_pair((unsigned int)2,
						 recv_iindex[elem].second));
    elem++;
  }
  assert((relem == _ngpyr) && (selem == nreal_pyr)); 
  elem = 0;
  selem = 0;
  relem = 0;
  recv_iindex.resize(_ngpris);
  while(elem < npris){
    if(elem < nreal_pris){
      _belem_send[2][selem++] = _bpris[elem].first;
      _sendcells[scell++] = gp.Elem2Cell(make_pair((unsigned int)3,
						   _bpris[elem].first));
    }
    else
      recv_iindex[relem++] = make_pair(_bpris[elem].second,
				       _bpris[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[2][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = gp.Elem2Cell(make_pair((unsigned int)3,
						 recv_iindex[elem].second));
    elem++;
  }
  assert((relem == _ngpris) && (selem == nreal_pris)); 
  elem = 0;
  selem = 0;
  relem = 0;
  recv_iindex.resize(_nghex);
  while(elem < nhex){
    if(elem < nreal_hex){
      _belem_send[3][selem++] = _bhex[elem].first;
      _sendcells[scell++] = gp.Elem2Cell(make_pair((unsigned int)4,
						   _bhex[elem].first));
    }
    else
      recv_iindex[relem++] = make_pair(_bhex[elem].second,
				       _bhex[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[3][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = gp.Elem2Cell(make_pair((unsigned int)4,
						 recv_iindex[elem].second));
    elem++;
  }
  assert((relem == _nghex) && (selem == nreal_hex)); 
}


void
PartitionBoundary::populate_local_arrays(const PartitionBoundary &rpb)
{
  // Total number of nodes on this boundary (shared+send+recv)
  unsigned int nnodes = _bnodes.size();
  // Number of local reals is local total minus local ghosts
  unsigned int nreal = _bnodes.size() - _ngnodes;
  // Number of remote recv'd nodes
  unsigned int nremote_gnodes = rpb._ngnodes;
  // Number of local shared is local real - remote ghosts
  unsigned int nshared = nreal - nremote_gnodes;
  // Num remote real is remote total minus remote ghosts
  unsigned int nremote_rnodes = rpb._bnodes.size() - nremote_gnodes;
  // Num local send nodes must equal num remote ghosts
  unsigned int nsend = nreal - nshared;
  assert(nsend == nremote_gnodes);
  _recvnodes.resize(_ngnodes);
  _sharenodes.resize(nshared);
  _sendnodes.resize(nsend);
  //
  // We always let the sender determine the node ordering,
  // so the recv arrays need to be re-ordered so that they
  // match the remote border's node ordering.  Also, we
  // let the lesser of the two partition_id's determine the
  // shared node order.  
  //
  // The format for the data in _bnodes[] at this point is:
  //
  // vector<pair<LocalNodeIndex,RemoteBoundaryNodeIndex> >
  // 
  // And we want to populate the send, recv, and shared 
  // node arrays in the correct order.  For those that need
  // the remote ordering, we swap the pairs above and sort
  // the array based on the RemoteBoundaryNodeIndex (which
  // indexes into the remote partition's boundary node 
  // array).  Then we take a second pass to populate the 
  // local indices in that ordering.
  //
  // These hold the inverted local<-->remote node indices
  vector<pair<unsigned int,unsigned int> > recv_iindex;
  vector<pair<unsigned int,unsigned int> > shared_iindex;
  recv_iindex.resize(_ngnodes);
  shared_iindex.resize(nshared);
  unsigned int node = 0;
  unsigned int shnode = 0;
  unsigned int snode = 0;
  unsigned int rnode = 0;
  while(node < nnodes){
    // All ghost nodes are processed here
    if(node >= nreal)
      recv_iindex[rnode++] = make_pair(_bnodes[node].second,
				       _bnodes[node].first);
    //      recv_iindex[rnode++] = _bnodes[node];
    else {
      if(_bnodes[node].second > nremote_rnodes)
	_sendnodes[snode++] = _bnodes[node].first;
      else {
	if(_mypart < _rpart) 
	  _sharenodes[shnode++] = _bnodes[node].first;
	else
	  shared_iindex[shnode++] = make_pair(_bnodes[node].second,
					      _bnodes[node].first);
      }
    }
    node++;
  }
  sort(recv_iindex.begin(),recv_iindex.end());
  node = 0;
  while(node < rnode){
    _recvnodes[node] = recv_iindex[node].second;
    node++;
  }
  // Check partition_id (encoded as _myrank and _rrank) and order shared
  // nodes if needed.
  if(_mypart > _rpart){
    node = 0; 
    sort(shared_iindex.begin(),shared_iindex.end());
    while(node < shnode){
      _sharenodes[node] = shared_iindex[node].second;
      node++;
    }
  }

  // Nodes are done.  Now we need send and recv element arrays. These will be
  // done independently for each element type.  
  unsigned int ntets = _btets.size();
  unsigned int nhex  = _bhex.size();
  unsigned int npyr  = _bpyr.size();
  unsigned int npris = _bpris.size();
  unsigned int nreal_tets = ntets - _ngtets;
  unsigned int nreal_hex = nhex - _nghex;
  unsigned int nreal_pyr = npyr - _ngpyr;
  unsigned int nreal_pris = npris - _ngpris;
  _belem_send[0].resize(nreal_tets);
  _belem_send[1].resize(nreal_pyr);
  _belem_send[2].resize(nreal_pris);
  _belem_send[3].resize(nreal_hex);
  _belem_recv[0].resize(_ngtets);
  _belem_recv[1].resize(_ngpyr);
  _belem_recv[2].resize(_ngpris);
  _belem_recv[3].resize(_nghex);
  _sendcells.resize(nreal_tets+nreal_pyr+nreal_pris+nreal_hex);
  _recvcells.resize(_ngtets+_ngpyr+_ngpris+_nghex);
  unsigned int elem = 0;
  unsigned int selem = 0;
  unsigned int relem = 0;
  unsigned int scell = 0;
  unsigned int rcell = 0;
  recv_iindex.resize(_ngtets);
  while(elem < ntets){
    if(elem < nreal_tets){
      _belem_send[0][selem++] = _btets[elem].first;
      _sendcells[scell++] = _btets[elem].first;
    }
    else
      recv_iindex[relem++] = make_pair(_btets[elem].second,
				       _btets[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[0][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = recv_iindex[elem].second;
    elem++;
  }
  assert((relem == _ngtets) && (selem == nreal_tets)); 
  elem = 0;
  selem = 0;
  relem = 0;
  recv_iindex.resize(_ngpyr);
  while(elem < npyr){
    if(elem < nreal_pyr){
      _belem_send[1][selem++] = _bpyr[elem].first;
      _sendcells[scell++] = _bpyr[elem].first + ntets;
    }
    else
      recv_iindex[relem++] = make_pair(_bpyr[elem].second,
				       _bpyr[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[1][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = recv_iindex[elem].second + ntets;
    elem++;
  }
  assert((relem == _ngpyr) && (selem == nreal_pyr)); 
  elem = 0;
  selem = 0;
  relem = 0;
  recv_iindex.resize(_ngpris);
  while(elem < npris){
    if(elem < nreal_pris){
      _belem_send[2][selem++] = _bpris[elem].first;
      _sendcells[scell++] = _bpris[elem].first + ntets + npyr;
    }
    else
      recv_iindex[relem++] = make_pair(_bpris[elem].second,
				       _bpris[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[2][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = recv_iindex[elem].second + ntets + npyr;
    elem++;
  }
  assert((relem == _ngpris) && (selem == nreal_pris)); 
  elem = 0;
  selem = 0;
  relem = 0;
  recv_iindex.resize(_nghex);
  while(elem < nhex){
    if(elem < nreal_hex){
      _belem_send[3][selem++] = _bhex[elem].first;
      _sendcells[scell++] = _bhex[elem].first+ntets+npyr+npris;
    }
    else
      recv_iindex[relem++] = make_pair(_bhex[elem].second,
				       _bhex[elem].first);
    elem++;
  }
  elem = 0;
  sort(recv_iindex.begin(),recv_iindex.end());
  while(elem < relem){
    _belem_recv[3][elem] = recv_iindex[elem].second;
    _recvcells[rcell++] = recv_iindex[elem].second+ntets+npyr+npris;
    elem++;
  }
  assert((relem == _nghex) && (selem == nreal_hex)); 
}


#if 0
Partition::Partition(const GEM_Partition &reg)
  : GEM_Partition(reg)
{
  //  _part = _id;
  //  reg._id = part._rank+1;
  //  reg._ngnodes = part._ngnodes;
  //  reg._ngtet = part._ngtet;
  //  reg._nghex = part._nghex;
  //  reg._ngpyr = part._ngpyr;
  //  reg._ngpris = part._ngpris;
  _patches.resize(_com
  reg._domain_b.resize(part._patches.size());
  reg._part_b.resize(part._boundaries.size());
  reg._nc = part._nc;
  reg._tetconn = part._tetconn;
  reg._hexconn = part._hexconn;
  reg._pyrconn = part._pyrconn;
  reg._prisconn = part._prisconn;

  // Populate all the partition boundary structures.  Lots of index juggling here to 
  // to match the mapping CELLS[0:NELEM-1]. NELEM = NTETS+NPYR+NPRIS+NHEX, in that 
  // exact order.  The Partitioner output has an array for each element type.
  unsigned int nbord = reg._part_b.size();
  unsigned int bord = 0;
  while(bord < nbord){

    // Nodes are straightforward copy of the data
    reg._part_b[bord]._rpbindex   = part._boundaries[bord]._rbid;
    reg._part_b[bord]._rpart      = part._boundaries[bord]._rrank + 1;
    reg._part_b[bord]._sendnodes  = part._boundaries[bord]._sendnodes;
    reg._part_b[bord]._recvnodes  = part._boundaries[bord]._recvnodes;
    reg._part_b[bord]._sharenodes = part._boundaries[bord]._sharenodes;

    // Cells must be re-arranged to the above described ordering
    unsigned int nsend = 0;
    unsigned int nrecv = 0;

    unsigned int ntets       = part._boundaries[bord]._btets.size();
    unsigned int nreal_tets  = ntets - part._boundaries[bord]._ngtets;
    unsigned int npyrs       = part._boundaries[bord]._bpyr.size();
    unsigned int nreal_pyrs  = npyrs - part._boundaries[bord]._ngpyr;
    unsigned int npriss      = part._boundaries[bord]._bpris.size();
    unsigned int nreal_priss = npriss - part._boundaries[bord]._ngpris;
    unsigned int nhexs       = part._boundaries[bord]._bhex.size();
    unsigned int nreal_hexs  = nhexs - part._boundaries[bord]._nghex;

    nsend = nreal_tets + nreal_pyrs + nreal_priss + nreal_hexs;
    nrecv = ntets - nreal_tets + npyrs - nreal_pyrs + npriss - nreal_priss 
      + nhexs - nreal_hexs;
    reg._part_b[bord]._sendcells.resize(nsend);
    reg._part_b[bord]._recvcells.resize(nrecv);

    vector<unsigned int> &sendv   = part._boundaries[bord]._belem_send[0];
    vector<unsigned int> &recvv   = part._boundaries[bord]._belem_recv[0];
    vector<unsigned int> &regsend = reg._part_b[bord]._sendcells;
    vector<unsigned int> &regrecv = reg._part_b[bord]._recvcells;

    unsigned int sendi = 0;
    unsigned int nerecv = ntets - nreal_tets;
    unsigned int nesend = nreal_tets;
    unsigned int curr_recv_pos = sendi;
    unsigned int curr_send_pos = sendi;

    while(sendi < nesend)
      regsend[curr_send_pos++] = sendv[sendi++];
    sendi = 0;
    while(sendi < nerecv)
      regrecv[curr_recv_pos++] = recvv[sendi++];

    sendv = part._boundaries[bord]._belem_send[1];
    recvv = part._boundaries[bord]._belem_recv[1];
    sendi = 0;
    nerecv = npyrs - nreal_pyrs;
    nesend = nreal_pyrs;
    while(sendi < nesend)
      regsend[curr_send_pos++] = sendv[sendi++];
    sendi = 0;
    while(sendi < nerecv)
      regrecv[curr_recv_pos++] = recvv[sendi++];

    sendv = part._boundaries[bord]._belem_send[2];
    recvv = part._boundaries[bord]._belem_recv[2];
    sendi = 0;
    nerecv = npriss - nreal_priss;
    nesend = nreal_priss;
    while(sendi < nesend)
      regsend[curr_send_pos++] = sendv[sendi++];
    sendi = 0;
    while(sendi < nerecv)
      regrecv[curr_recv_pos++] = recvv[sendi++];

    sendv = part._boundaries[bord]._belem_send[3];
    recvv = part._boundaries[bord]._belem_recv[3];
    sendi = 0;
    nerecv = nhexs - nreal_hexs;
    nesend = nreal_hexs;
    while(sendi < nesend)
      regsend[curr_send_pos++] = sendv[sendi++];
    sendi = 0;
    while(sendi < nerecv)
      regrecv[curr_recv_pos++] = recvv[sendi++];

    assert(curr_recv_pos == regrecv.size() && curr_send_pos == regsend.size());
    bord++;
  }      

  // Populate the domain boundaries - straightforward data copy
  unsigned int local_patch = 0;
  unsigned int nlocal_patches = reg._domain_b.size();
  while(local_patch < nlocal_patches){
    PartitionPatch *pp = &part._patches[local_patch];
    reg._domain_b[local_patch]._id       = (pp->_id < 0 ? (-1*pp->_id) : pp->_id);
    reg._domain_b[local_patch]._ngtri    = pp->_ngtri;
    reg._domain_b[local_patch]._ngquad   = pp->_ngquad;
    reg._domain_b[local_patch]._triconn  = pp->_triconn;
    reg._domain_b[local_patch]._quadconn = pp->_quadconn;
    local_patch++;
  }
  return;
}
#endif


  






