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
#ifndef _FLUREGION_H_
#define _FLUREGION_H_

class FluVolumeSoln {
public:
  double _current_time;
  double _resid;
  bool   _unsteady;
  vector<double> _rhof;  // Density
  vector<double> _rhovf; // Momentum density
  vector<double> _rhoEf; // Energy density
  vector<double> _pf;    // Pressure
  vector<double> _Tf;    // Temperature
  vector<double> _af;    // Speed of sound
  vector<double> _disp;  // nodal displacements
  vector<double> _gsp;   // Speed of each face
  FluVolumeSoln(unsigned int ncells=0,unsigned int nfaces=0,
		unsigned int nnodes=0)
  {
    Resize(ncells,nfaces,nnodes);
  };
  void
  Resize(unsigned int ncells=0,unsigned int nfaces=0,unsigned int nnodes=0)
  {
    _rhof.resize(ncells,0.0);
    _rhovf.resize(3*ncells,0.0);
    _rhoEf.resize(ncells,0.0);
    _pf.resize(ncells,0.0);
    _Tf.resize(ncells,0.0);
    _af.resize(ncells,0.0);
    _disp.resize(nnodes*3,0.0);
    _gsp.resize(nfaces,0.0);
    _unsteady = false;
    _resid = 0.0;
    _current_time = 0.0;
  };
};

class FluSurfSoln {
public:
  vector<double> _du_alp;     // Vertex displacement
  vector<double> _rhofvf_alp; // Momentum flux
  vector<double> _nf_alp;     // face normal
  vector<double> _rhof_alp;   // density
  vector<double> _pf;         // pressure
  vector<double> _qc;         // convective heat flux
  vector<double> _qr;         // radiative heat flux
  vector<double> _tf;         // traction
  vector<double> _Tb_alp;     // Boundary temperature
  vector<double> _mdot_alp;   // Injection mass flux
  vector<double> _Tflm_alp;   // Flame temperature
  vector<double> _gsp;        // Grid Speeds
  vector<int>    _bflag;      // Burning flag
  
  FluSurfSoln(unsigned int nfaces=0,unsigned int nvert=0)
  {
    Resize(nfaces,nvert);
  };
  void Resize(unsigned int nfaces=0,unsigned int nvert=0)
  {
    _gsp.resize(nfaces,0);
    _du_alp.resize(nvert*3,0.0);
    _rhofvf_alp.resize(3*nfaces,0.0);
    _nf_alp.resize(3*nfaces,0.0);
    _rhof_alp.resize(nfaces,0.0);
    _pf.resize(nfaces,0.0);
    _qc.resize(nfaces,0.0);
    _qr.resize(nfaces,0.0);
    _tf.resize(3*nfaces,0.0);
    _Tb_alp.resize(nfaces,0.0);
    _mdot_alp.resize(nfaces,0.0);
    _Tflm_alp.resize(nfaces,0.0);
    _bflag.resize(nfaces,0);
  };  
};

class FluPatch : public GEM_DomainBoundary {
public:
  //  unsigned int _patchNo; = _id
  unsigned int _ntri;
  unsigned int _nquad;
  FluSurfSoln _soln;
  bool InitSurfaceSoln(const string &prefix="",bool=false);
  //#ifdef _ROCSTAR_X_
  unsigned int _bcflag; // 0:nonburning, 1:burning, 2:noninteracting
  int _constr_type;
  bool RegisterSoln(const string &,bool);
  //#endif
};

class FluBorder : public GEM_PartitionBoundary {
public:
  unsigned int _rbid;
};
  

class FluRegion : public GEM_Partition {
public:
  string               _casename;
  string               _casepath;
  //#ifdef _ROCSTAR_X_
  //  string               _rocstar_outpath;
  //#endif
  unsigned int         _nnodes;
  unsigned int         _ntet;
  unsigned int         _nhex;
  unsigned int         _npyr;
  unsigned int         _npris;
  unsigned int         _npatches_total;
  vector<FluPatch>     _patches;
  vector<FluBorder>    _borders;
  FluVolumeSoln        _soln;
  vector<pair<unsigned int,unsigned int> > _cellmap;
  void report();
  unsigned int nelem();
  bool BuildCellMapping();
  bool Cell2Elem(unsigned int, unsigned  
		 int&, unsigned int&);
  bool BuildPatchMapping(map<unsigned int,unsigned int> &,const string &p="");
  bool WriteFluNative(const string &p="");
  bool WriteFluCMP(const string &p="");
  bool ReadFluCOM(const string &p="");
  bool WriteFluCOM(const string &p="");
  bool ReadFluDIM(const string &p="",double=0.0,bool=true);
  bool WriteFluDIM(const string &p="",double=0.0,bool=true);
  bool ReadFluGridASCII(const std::string&p="", double=0.0,bool=true);
  bool WriteFluGridASCII(const string &pre="", double=0.0,bool=true);
  bool ReadFluSolnASCII(const std::string&p="",unsigned int=0, 
			double=0.0,bool=true);
  bool WriteFluSolnASCII(const std::string&p="",unsigned int=0, 
			 double=0.0,bool=true);
  bool WriteGSPASCII(const std::string&p="",unsigned int=0,
		     double=0.0,bool=true);
  bool ReadGSPASCII(const std::string&p="",unsigned int=0,
		    double=0.0,bool=true);
  bool ReadRegionASCII(const string &p="",unsigned int=0,unsigned int=0,
		       double=0.0,bool=true);
  bool InitSurfaceSoln(const string &p="");
  bool InitVolumeSoln(const string &p="");
  bool CreateRegionMapFile(const string &p="",unsigned int=0,unsigned int=0);
  void PopRemBordIndMPI();
  bool PopRemBordIndFILE(const string &p="",double=0.0,bool=true);
  bool ComRemeshInitData(const string &wname,double *cell_data,int nval_cells,
			 double *node_data,int nval_nodes);
  bool ReadControlFile();
  //  bool validate_comm_list(int ncsend,int ncrecv,int *csend,int *crecv);
  //  bool PopulateFluVolumeSoln(double *data);
  FluRegion()
    : GEM_Partition()
  {
    _casename = "";
    _nnodes = 0;
    _ntet   = 0;
    _npyr   = 0;
    _npris  = 0;
    _nhex   = 0;
    _npatches_total = 0;
    _patches.resize(0);
    _borders.resize(0);
    _cell_ordering[0] = 1;
    _cell_ordering[1] = 4;
    _cell_ordering[2] = 3;
    _cell_ordering[3] = 2;
  };
  FluRegion(const GEM_Partition &gp)
    : GEM_Partition(gp)
  {
    _casename = "";
    _nnodes  = gp._nc.size()/3;
    _ntet    = gp._tetconn.size()/4;
    _npyr    = gp._pyrconn.size()/5;
    _npris   = gp._prisconn.size()/6;
    _nhex    = gp._hexconn.size()/8;
    _cell_ordering[0] = 1;
    _cell_ordering[1] = 4;
    _cell_ordering[2] = 3;
    _cell_ordering[3] = 2;
  };
  template<typename DB>
  bool PopulateFluPatches(const string &pre,vector<DB> &db)
  {
    map<unsigned int,unsigned int> patch_mapping;
    BuildPatchMapping(patch_mapping,pre);
    MapDomainBoundaries(db,patch_mapping,_patches);
    vector<FluPatch>::iterator fpi = _patches.begin();
    while(fpi != _patches.end()){
      fpi->_ntri  = fpi->_triconn.size()/3;
      fpi->_nquad = fpi->_quadconn.size()/4;
      fpi++;
    }
    return(true);
  };

  template<typename PB>
  bool PopulateFluBorders(vector<PB> &pb)
  {
    unsigned int nborders = pb.size();
    _borders.resize(nborders);
    unsigned int border = 0;
    while(border < nborders){
      _borders[border]._rpart      = pb[border]._rpart;
      _borders[border]._sendcells  = pb[border]._sendcells;
      _borders[border]._recvcells  = pb[border]._recvcells;
      _borders[border]._sharenodes = pb[border]._sharenodes;
      _borders[border]._sendnodes  = pb[border]._sendnodes;
      _borders[border]._recvnodes  = pb[border]._recvnodes;
      _borders[border]._rbid       = pb[border]._rbid;
      _borders[border]._out        = pb[border]._out;
      border++;
    }
    return(true);
  };
  bool
  validate_comm_list(int ncsend,int ncrecv,int *csend,int *crecv)
  {
    int index = 0;
    int nreal_cell = _tetconn.size()/4 + _prisconn.size()/6 +
      _pyrconn.size()/5 + _hexconn.size()/8 - (_ngtet + _ngpris +
					       _ngpyr + _nghex);
    bool rval = true;
    while(index < ncsend){
      int ind = index++;
      if(!(csend[ind] <= nreal_cell)){
	if(_out)
	  *_out << "SEND CELL " << index << " is a ghost cell!!" << endl;
	rval = false;
      }
      if(!(csend[ind] > 0)){
	if(_out)
	  *_out << "SEND CELL " << index << " is zero or negative!" << endl;
	rval = false;
      }
    }
    index = 0;
    list<int> recvcell_list;
    while(index < ncrecv) {
      int ind = index++;
      if(!(crecv[ind] > nreal_cell)){
	if(_out)
	  *_out << "RECV CELL " << index << " is a real cell!!" << endl;
	rval = false;
      }
      if(!(crecv[ind] > 0)){
	if(_out)
	  *_out << "RECV CELL " << index << " is zero or negative!" << endl;
	rval = false;
      }
      bool duped = false;
      list<int>::iterator rci = recvcell_list.begin();
      while(rci != recvcell_list.end() && !duped){
	if(crecv[ind] == *rci++){
	  if(_out)
	    *_out << "RECV CELL " << index 
		  << " is duplicated in the receive list!" 
		  << endl;
	  duped = true;
	}
      }
      if(!duped)
	recvcell_list.push_back(crecv[ind]);
    }
    return(rval);
  };
  void AddFluBorder(int rpid,int rpin,int nnshare, int nnsend,
		    int nnrecv,int ncsend,int ncrecv,
		    int *nshared,int *nsend,int *nrecv,
		    int *csend,int *crecv){
    assert(rpid > 0 && rpin >= 0);
    if(_debug && _out)
      *_out << "FluRegion(" << _id << ")::AddFluBorder: Adding Border with"
	    << " region " << rpid << "." << endl;
    FluBorder new_border;
    new_border._rbid = rpin;
    new_border._out = _out;
    if(!validate_comm_list(ncsend,ncrecv,csend,crecv)){
      if(_out)
	*_out << "FluRegion(" << _id << ")::AddFluBorder: Validation of "
	      << "communication arrays failed, aborting." << endl;
      exit(-1);
    }
    new_border.populate(rpid,nnshare,nnsend,nnrecv,ncsend,ncrecv,
			nshared,nsend,nrecv,csend,crecv);
    _borders.push_back(new_border);
  };
  void AddFluPatch(int patch_id,int ntri, int ngtri, int *tris,
		   int nquad,int ngquad, int *quads)
  {
    assert(ntri >= ngtri && nquad >= ngquad);
    if(_debug && _out)
      *_out << "FluRegion(" << _id << ")::AddFluPatch: Adding Patch with"
	    << " id " << patch_id << "." << endl;
    FluPatch new_patch;
    new_patch._ntri = ntri;
    new_patch._nquad = nquad;
    new_patch._id = patch_id;
    new_patch._ngtri = ngtri;
    new_patch._ngquad = ngquad;
    int indy = 0;
    new_patch._triconn.resize(3*ntri);
    new_patch._quadconn.resize(4*nquad);
    new_patch._out = _out;
    while(indy < 3*ntri){
      assert(tris[indy] != 0);
      new_patch._triconn[indy] = tris[indy];
      assert(new_patch._triconn[indy] != 0);
      indy++;
    }
    indy = 0;
    while(indy < 4*nquad){
      assert(quads[indy] != 0);
      new_patch._quadconn[indy] = quads[indy];
      indy++;
    }
    _patches.push_back(new_patch);
  };
  //#ifdef _ROCSTAR_X_
  bool RegisterFluSurfaceMesh();
  bool RegisterSurfaceSoln(bool);
  bool RegisterVolumeSoln(bool);
  //#endif
};

istream &SkipLines(istream &,unsigned int = 1);

#endif






