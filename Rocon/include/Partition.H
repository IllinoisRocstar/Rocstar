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
#ifndef _GEM_PARTITION_H_
#define _GEM_PARTITION_H_


// Partition boundary refers to an interface between two partitions
// where communication is needed
class PartitionBoundary : public GEM_PartitionBoundary {
public:
  unsigned int _mypart;
  unsigned int _rbid;
  vector< pair<unsigned int,unsigned int> >  _bnodes;
  unsigned int _ngnodes;
  vector< pair<unsigned int,unsigned int> >  _btets;
  unsigned int _ngtets;
  vector< pair<unsigned int, unsigned int> > _bhex;
  unsigned int _nghex;
  vector< pair<unsigned int, unsigned int> > _bpyr;
  unsigned int _ngpyr;
  vector< pair<unsigned int, unsigned int> > _bpris;
  unsigned int _ngpris;
  vector<unsigned int> _belem_send[4];
  vector<unsigned int> _belem_recv[4];
  PartitionBoundary()
    : GEM_PartitionBoundary()
  {};
  PartitionBoundary(const GEM_PartitionBoundary &pb)
    : GEM_PartitionBoundary(pb)
  {
    _mypart  = 0;
    _rbid    = 0;
    _ngnodes = pb._recvnodes.size();
  };
  bool ReadPartitionBoundary(ifstream &Inf);
  bool ReadPartitionBoundaryASCII(ifstream &Inf);
  bool WritePartitionBoundaryASCII(ofstream &Ouf);
  unsigned int MyRank(void) const { return _mypart; };
  void MyRank(unsigned int n) {_mypart = n;};
  void populate_local_arrays(const PartitionBoundary &);
  void populate_local_arrays(const PartitionBoundary &rpb,
			     GEM_Partition &gp);
  //  void populate_element_arrays(const GEM_Partition &);
};

// Partition patch refers to the solver's domain boundary - or it's 
// interface to the rest of the physical universe.
class PartitionPatch : public GEM_DomainBoundary {
public:
  unsigned int _nfaces;
  unsigned int _ngfaces;
  vector<unsigned int> _faceconn;
  bool ReadPartitionPatchASCII(ifstream &Inf);
  bool WritePartitionPatchASCII(ofstream &Ouf);
  bool ReadPartitionPatch(ifstream &Inf);
  PartitionPatch() 
    : GEM_DomainBoundary()
  {
    _nfaces = 0;
    _ngfaces = 0;
  };
  PartitionPatch(const GEM_DomainBoundary &db)
    : GEM_DomainBoundary(db)
  {
    _nfaces  = db._triconn.size()/3 + db._quadconn.size()/4;
    _ngfaces = db._ngtri + db._ngquad;
  }; 
};

// A partition is whole or part of a mesh representing a solver's domain
class Partition : public GEM_Partition {
public:
  unsigned int _nnodes;
  unsigned int _nelem;
  unsigned int _ntet;
  unsigned int _nhex;
  unsigned int _npyr;
  unsigned int _npris;
  unsigned int _npatch;
  vector<PartitionBoundary> _boundaries;
  vector<PartitionPatch>    _patches;
  bool ReadPartition(const string &pre,unsigned int rank);
  bool ReadPartitionBoundaries(ifstream &Inf);
  bool ReadPartitionPatches(ifstream &Inf);
  bool ReadPartitionASCII(const string &pre,unsigned int rank);
  bool ReadPartitionBoundariesASCII(ifstream &Inf);
  bool ReadPartitionPatchesASCII(ifstream &Inf);
  bool WritePartitionASCII(const string &pre);
  bool WritePartitionBoundariesASCII(ofstream &Ouf);
  bool WritePartitionPatchesASCII(ofstream &Ouf);
  bool BuildPartitionBoundaries(const vector<int> &Pconn);
  Partition(const GEM_Partition &gp);
  Partition()
  {};
};
  
#endif






