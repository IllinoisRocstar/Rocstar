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
#ifndef _FLOGRID_H_
#define _FLOGRID_H_


class FloGridBlock {
protected:
  unsigned int              _ncells[3];
  std::vector<double>       _incoords;
  unsigned int              _ng;
  std::vector<double>       _soln;
  unsigned int              _grid_no;
private:
  bool                      flolayout;
public:
  FloGridBlock() 
  { 
    _ncells[0] = _ncells[1] = _ncells[2] = _ng = 0;
    _incoords.resize(0); 
    _soln.resize(0);
    _grid_no = 0;
    flolayout = false;
  };
  FloGridBlock(unsigned int i,unsigned int j,unsigned int k)
  {
    _ncells[0] = i;
    _ncells[1] = j;
    _ncells[2] = k;
    _ng = 0;
    _grid_no = 0;
    _incoords.resize(3*(i+1)*(j+1)*(k+1));
    flolayout = false;
  };
  FloGridBlock(const FloGridBlock &inblock)
  {
    _ncells[0] = inblock._ncells[0];
    _ncells[1] = inblock._ncells[1];
    _ncells[2] = inblock._ncells[2];
    _grid_no = inblock._grid_no;
    _incoords.resize(3*(inblock._ncells[0]+1)*(inblock._ncells[1]+1)*
		     (inblock._ncells[2]+1));
    std::vector<double>::const_iterator ii = inblock._incoords.begin();
    std::vector<double>::iterator ni = _incoords.begin();
    while(ni != _incoords.end())
      *ni++ = *ii++;
    flolayout = inblock.flolayout;
    _soln = inblock._soln;
    _ng = inblock._ng;
  };
  void DestroyGrid()
  {
    _incoords.resize(0);
    // release the memory for real...
    std::vector<double>(_incoords).swap(_incoords);
  };
  void SetGhostLayers(unsigned int n)
  {
    _ng = n;
  };
  unsigned int NGhostLayers() const {return (_ng); };
  int CreateSoln()
  {
    _soln.resize(5*NTotalCells());
  };
  unsigned int NNodes() const { return ((_ncells[0]+1)*(_ncells[1]+1)*(_ncells[2]+1)); };
  unsigned int NumCells() const { return (_ncells[0]*_ncells[1]*_ncells[2]); };
  unsigned int NTotalCells() const { return ((_ncells[0]+(2*_ng))*(_ncells[1]+(2*_ng))*(_ncells[2]+(2*_ng))); };
  std::vector<double> &Solution() { return (_soln); };
  const std::vector<double> &Solution() const { return (_soln); };
  void FloLayout()
  {
    if(!flolayout){
      flolayout = true;
      unsigned int nnodes = NNodes();
      std::vector<double> coordcopy(3*nnodes);
      std::vector<double>::iterator ci = _incoords.begin();
      std::vector<double>::iterator cci = coordcopy.begin();
      while(ci != _incoords.end())
	*cci++ = *ci++;
      for(int i = 0;i < nnodes;i++){
	_incoords[i]   = coordcopy[3*i];
	_incoords[nnodes+i] = coordcopy[(3*i) + 1];
	_incoords[2*nnodes+i] = coordcopy[(3*i) + 2];
      }
      
    }
  };
  void MeshLayout()
  {
    if(flolayout){
      flolayout = false;
      unsigned int nnodes = NNodes();
      std::vector<double> coordcopy(3*nnodes);
      std::vector<double>::iterator ci = _incoords.begin();
      std::vector<double>::iterator cci = coordcopy.begin();
      while(ci != _incoords.end())
	*cci++ = *ci++;
      for(int i = 0;i < nnodes;i++){
	_incoords[3*i]   = coordcopy[i];
	_incoords[3*i+1] = coordcopy[nnodes + i];
	_incoords[3*i+2] = coordcopy[(2*nnodes) + i];
      }
      
    }
  };

  FloGridBlock &operator=(const FloGridBlock &inblock)
  {
    _ncells[0] = inblock._ncells[0];
    _ncells[1] = inblock._ncells[1];
    _ncells[2] = inblock._ncells[2];
    
    _incoords.resize(3*(inblock._ncells[0]+1)*(inblock._ncells[1]+1)*
		     (inblock._ncells[2]+1));
    std::vector<double>::const_iterator ii = inblock._incoords.begin();
    std::vector<double>::iterator ni = _incoords.begin();
    while(ni != _incoords.end())
      *ni++ = *ii++;
    flolayout = inblock.flolayout;

  };
  unsigned int isize() const { return _ncells[0]; };
  unsigned int jsize() const { return _ncells[1]; };
  unsigned int ksize() const { return _ncells[2]; };
  unsigned int *NCells() { return (_ncells); };
  const unsigned int *NCells() const { return(_ncells); };
  void Resize(unsigned int i,unsigned int j,unsigned int k)
  {
    _ncells[0] = i; _ncells[1] = j; _ncells[2] = k;
    _incoords.resize(3*(i+1)*(j+1)*(k+1));
    flolayout = false;
  };
  void ResizeSoln(unsigned int i,unsigned int j,unsigned int k,unsigned int ng)
  {
    _ncells[0] = i; _ncells[1] = j; _ncells[2] = k;
    _ng = ng;
    _soln.resize(5*((i+(2*ng))*(j+(2*ng))*(k+(2*ng))));
    //    flolayout = false;
  };
  std::vector<double> &Coords() { return _incoords; };
  const std::vector<double> &Coords() const { return _incoords; };
  int ReadBlock(std::istream &BlockIn)
  {
    std::string line;
    std::getline(BlockIn,line);
    if(line.empty())
      return(1);
    std::istringstream Istr(line);
    unsigned int ni      = 0;
    unsigned int nj      = 0;
    unsigned int nk      = 0;
    Istr >> _grid_no >> ni >> nj >> nk;
    if((_grid_no <= 0) || (ni <= 0) || (nj <= 0) || (nk <= 0)){
      std::cout << "FloGridBlock::ReadBlock: Error reading block" << std::endl
		<< "FloGridBlock::ReadBlock: line:" << line << std::endl
		<< "FloGridBlock::ReadBlock: grid no: " << _grid_no << "ni = " << ni 
		<< " nj = " << nj << " nk = " << nk << std::endl;
      return(1);
    }
    unsigned int number_of_nodes = (ni+1)*(nj+1)*(nk+1);
    //    std::cout << "Partition " << grid_no << ": (" << ni
    //	      << "," << nj << "," << nk << ")" << std::endl
    //	      << "npts = " << number_of_nodes << std::endl;
    this->Resize(ni,nj,nk);
    double VERYSMALL = std::numeric_limits<double>::min();
    double VERYLARGE = std::numeric_limits<double>::max();
    std::vector<double>::iterator ici = _incoords.begin();
    while(ici != _incoords.end())
      BlockIn >> *ici++;
    std::getline(BlockIn,line);
    GeoPrim::CPoint minp(VERYLARGE,VERYLARGE,VERYLARGE);
    GeoPrim::CPoint maxp(-1.0*VERYLARGE,-1.0*VERYLARGE,-1.0*VERYLARGE);
    ici = _incoords.begin();
    unsigned int ncount = 0;
    while(ici != _incoords.end()){
      double *minval = NULL;
      double *maxval = NULL;
      double value = *ici++;
      if(ncount < number_of_nodes){
	minval = &minp.x();
	maxval = &maxp.x();
      }
      else if(ncount < 2*number_of_nodes){
	minval = &minp.y();
	maxval = &maxp.y();
      }
      else {
	minval = &minp.z();
	maxval = &maxp.z();
      }
      ncount++;
      if (value < *minval)
	*minval = value;
      if (value > *maxval)
	*maxval = value;
    }
    //    std::cout << "Mesh limits: " << std::endl
    //	      << "Minimum point: " << minp << std::endl
    //      << "Maximum point: " << maxp << std::endl;
    flolayout = true;
    return(0);
  };
  int ReadBlockSoln(std::istream &BlockSolnIn)
  {
    std::string line;
    std::getline(BlockSolnIn,line);
    if(line.empty())
      return(1);
    std::istringstream Istr(line);
    unsigned int grid_no = 0;
    unsigned int ni      = 0;
    unsigned int nj      = 0;
    unsigned int nk      = 0;
    unsigned int ng      = 0;
    Istr >> grid_no >> ni >> nj >> nk >> ng;
    if(_grid_no != grid_no) std::cout << "grid " << _grid_no << " tried to read soln for block " << grid_no << std::endl;
    assert(_grid_no == grid_no);
    if((grid_no <= 0) || (ni <= 0) || (nj <= 0) || (nk <= 0)){
      std::cout << "FloGridBlock::ReadBlockSoln: Error reading block" << std::endl
		<< "FloGridBlock::ReadBlockSoln: line:" << line << std::endl
		<< "FloGridBlock::ReadBlockSoln: grid no: " << grid_no << "ni = " << ni 
		<< " nj = " << nj << " nk = " << nk << std::endl;
      return(1);
    }
    unsigned int number_of_total_cells = (ni+(2*ng))*(nj+(2*ng))*(nk+(2*ng));
    unsigned int number_of_real_cells  = (ni*nj*nk);
    assert(number_of_real_cells == NumCells());
    //    std::cout << "Partition " << grid_no << ": (" << ni
    //	      << "," << nj << "," << nk << ")" << std::endl
    //	      << "npts = " << number_of_nodes << std::endl;
    this->ResizeSoln(ni,nj,nk,ng);
    std::vector<double>::iterator ici = _soln.begin();
    while(ici != _soln.end())
      BlockSolnIn >> *ici++;
    std::getline(BlockSolnIn,line);
    return(0);
  };
  int WriteBlockSoln(std::ostream &Ouf)
  {
    std::ostringstream Ostr;
    Ostr << _grid_no << " " << _ncells[0] << " " << _ncells[1]
	<< " " << _ncells[2] << " " << _ng << std::endl;
    std::vector<double>::iterator si = _soln.begin();
    while(si != _soln.end()){
      Ostr << std::scientific << std::setprecision(16) << *si++ << " ";
    }
    Ostr << std::endl;
    Ouf << Ostr.str();
    return(Ostr.str().size());
  };
  
};


class FloGrid {
protected:
  std::string               _casename;
  std::vector<FloGridBlock> _blocks;
  std::ofstream _solnfile;
  double _time;
  double _unknown_number;
public:
  FloGrid() : _casename("")
  { _blocks.resize(0); _time=0.0;};
  FloGrid(const std::string &casename) 
    : _casename(casename)
  {_blocks.resize(0); _time=0.0;};
  std::string Casename() const { return (_casename); };
  double Time() const { return (_time); };
  std::vector<FloGridBlock> &Blocks() { return (_blocks); };
  const std::vector<FloGridBlock> &Blocks() const { return (_blocks); };
  double UnknownNumber() const { return (_unknown_number); };
  int ReadAllBlocks()
  {
    std::string sourcegridfile(_casename+".grda_0.00000E-00");
    //    std::cout << "Reading grid file: " << sourcegridfile << std::endl;
    
    std::ifstream Inf;
    Inf.open(sourcegridfile.c_str());
    if(!Inf){
      std::cout << "FloGrid::ReadAllBlocks: Could not open source grid file, " 
		<< sourcegridfile << std::endl;
      return(1);
    }
    
    std::string ftimestamp;
    std::getline(Inf,ftimestamp);
    while(Inf){
      FloGridBlock block;
      if(block.ReadBlock(Inf) == 0){
	_blocks.push_back(block);
      }
      
    }
    Inf.close();
    return(0);
  };
  int ReadNBlocks(int first,int nblocks)
  {
    std::string sourcegridfile(_casename+".grda_0.00000E-00");
    //    std::cout << "Reading grid file: " << sourcegridfile << std::endl;
    
    std::ifstream Inf;
    Inf.open(sourcegridfile.c_str());
    if(!Inf){
      std::cout << "FloGrid::ReadAllBlocks: Could not open source grid file, " 
		<< sourcegridfile << std::endl;
      return(1);
    }
    
    std::string ftimestamp;
    std::getline(Inf,ftimestamp);
    int n = 1;
    while(n < first){
      FloGridBlock block;
      block.ReadBlock(Inf);
      n++;
    }
    while(n < (first + nblocks)){
      FloGridBlock block;
      if(block.ReadBlock(Inf) == 0){
	_blocks.push_back(block);
      }
      n++;
    }
    Inf.close();
    return(0);
  };
  int BlockCount()
  {
    std::string sourcegridfile(_casename+".grda_0.00000E-00");
    std::cout << "Reading grid file: " << sourcegridfile << std::endl;
    
    std::ifstream Inf;
    Inf.open(sourcegridfile.c_str());
    if(!Inf){
      std::cout << "FloGrid::ReadAllBlocks: Could not open source grid file, " 
		<< sourcegridfile << std::endl;
      return(1);
    }
    
    std::string ftimestamp;
    std::getline(Inf,ftimestamp);
    int n = 0;
    while(Inf){
      FloGridBlock block;
      block.ReadBlock(Inf);
      n++;
    }
    Inf.close();
    return(n-1);
  };
  void DestroyGrids()
  {
    std::vector<FloGridBlock>::iterator fgbi = _blocks.begin();
    while(fgbi != _blocks.end())
      fgbi++->DestroyGrid();
  };
  void SetGhostLayers(unsigned int n)
  {
    std::vector<FloGridBlock>::iterator fgbi = _blocks.begin();
    while(fgbi != _blocks.end())
      fgbi++->SetGhostLayers(n);
  };
  void CreateSolutions()
  {
    std::vector<FloGridBlock>::iterator fgbi = _blocks.begin();
    while(fgbi != _blocks.end())
      fgbi++->CreateSoln();    
  };
  int OpenSolutionFile(const std::string &timestamp,double time,double unknown_number)
  {
    std::string solnfile(_casename+".sola_"+timestamp);
    
    //    std::ofstream Ouf;
    _solnfile.open(solnfile.c_str());
    if(!_solnfile){
      std::cout << "FloGrid::OpenSolutionFile Could not open source soln file, " 
		<< solnfile << std::endl;
      return(1);
    }
    _solnfile << std::scientific << std::setprecision(16);
    _solnfile << time << "\t" << unknown_number << std::endl;
    return(0);
  };
  void CloseSolutionFile() { _solnfile.close(); };
  std::ofstream &SolnFile() { return (_solnfile); };
  int WriteSolutionFile(const std::string &timestamp,double time,double unknown_number)
  {
    if(OpenSolutionFile(timestamp,time,unknown_number))
      return(1);
    std::vector<FloGridBlock>::iterator fgbi = _blocks.begin();
    while(fgbi != _blocks.end()){
      fgbi++->WriteBlockSoln(_solnfile);
      if(!_solnfile && (fgbi != _blocks.end())) 
	return(1);
      
    }
    CloseSolutionFile();
    return(0);
    
  };
  int WriteBlocks(std::ostream &Ouf)
  {
    //    std::string sourcesolnfile(_casename+".sola_"+timestamp);
    //    std::cout << "Writing soln file: " << sourcesolnfile << std::endl;
    
    //    std::ofstream Ouf;
    //    Ouf.open(sourcesolnfile.c_str());
    //    if(!Ouf){
    //      std::cout << "FloGrid::WriteSolutionFile Could not open source soln file, " 
    //		<< sourcesolnfile << std::endl;
    //      return(1);
    //    }
    //    Ouf << std::scientific << std::setprecision(16);
    //    Ouf << time << "\t" << unknown_number << std::endl;
    int datasize = 0;
    std::vector<FloGridBlock>::iterator fgbi = _blocks.begin();
    while(fgbi != _blocks.end())
      datasize += (fgbi++->WriteBlockSoln(Ouf));
    return(datasize);    
  };
  int ReadAllSolutions(const std::string &timestamp)
  {
    std::string sourcesolnfile(_casename+".sola_"+timestamp);
    //    std::cout << "Reading soln file: " << sourcesolnfile << std::endl;
    
    std::ifstream Inf;
    Inf.open(sourcesolnfile.c_str());
    if(!Inf){
      std::cout << "FloGrid::ReadAllSoln: Could not open source soln file, " 
		<< sourcesolnfile << std::endl;
      return(1);
    }
    
    std::string ftimestamp;
    std::getline(Inf,ftimestamp);
    std::istringstream Istr(ftimestamp);
    Istr >> _time >> _unknown_number;
    std::vector<FloGridBlock>::iterator fgbi = _blocks.begin();
    while(fgbi != _blocks.end()){
      fgbi++->ReadBlockSoln(Inf);
      if(!Inf && (fgbi != _blocks.end())) 
	return(1);
      
    }
    Inf.close();
    return(0);
  };
};

class FloVolumeSoln {
public:
  double _current_time;
  double _resid;
  bool   _unsteady;
  std::vector<double> _rhof;  // Density
  std::vector<double> _rhouf;
  std::vector<double> _rhovf; // Momentum density
  std::vector<double> _rhowf;
  std::vector<double> _rhoEf; // Energy density
  std::vector<double> _pf;    // Pressure
  std::vector<double> _Tf;    // Temperature
  std::vector<double> _af;    // Speed of sound
  std::vector<double> _disp;  // nodal displacements
  //  std::vector<double> _gsp;   // Speed of each face
  FloVolumeSoln(unsigned int ncell_i=0, unsigned int ncell_j=0, unsigned int ncell_k=0)
    
  {
    Resize(ncell_i,ncell_j,ncell_k);
  };
  void
  Resize(unsigned int ncell_i,unsigned int ncell_j,unsigned int ncell_k)
  {
    unsigned int ncells = ncell_i * ncell_j * ncell_k;
    unsigned int number_of_nodes = (ncell_i+1)*(ncell_j+1)*(ncell_k+1);
    _rhof.resize(ncells,0.0);
    _rhouf.resize(ncells,0.0);
    _rhovf.resize(ncells,0.0);
    _rhowf.resize(ncells,0.0);
    _rhoEf.resize(ncells,0.0);
    _pf.resize(ncells,0.0);
    _Tf.resize(ncells,0.0);
    _af.resize(ncells,0.0);
    _disp.resize(number_of_nodes*3,0.0);
    //    _gsp.resize(nfaces,0.0);
    _unsteady = false;
    _resid = 0.0;
    _current_time = 0.0;
  };
};

//class FloBorder : public GEM_PartitionBoundary 
//{
//public:
//  unsigned int _rbid;
//};


// class FloRegion : public GEM_Partition {
// public:
//   string               _casename;
//   string               _casepath;
//   unsigned int         _nnodes;
//   unsigned int         _ncells;
//   unsigned int         _ni;
//   unsigned int         _nj;
//   unsigned int         _nk;
//   //  unsigned int         _npatches_total;
//   //  std::vector<FluPatch>     _patches;
//   std::vector<FloBorder>    _borders;
//   FloVolumeSoln        _soln;
//   //  std::vector<pair<unsigned int,unsigned int> > _cellmap;
//   //   void report();
//   //   unsigned int nelem();
//   //   bool BuildCellMapping();
//   //   bool Cell2Elem(unsigned int, unsigned  
//   // 		 int&, unsigned int&);
//   //   bool BuildPatchMapping(map<unsigned int,unsigned int> &,const string &p="");
//   //   bool WriteFluNative(const string &p="");
//   //   bool WriteFluCMP(const string &p="");
//   //   bool ReadFluCOM(const string &p="");
//   //   bool WriteFluCOM(const string &p="");
//   //   bool ReadFluDIM(const string &p="",double=0.0,bool=true);
//   //   bool WriteFluDIM(const string &p="",double=0.0,bool=true);
//   //   bool ReadFluGridASCII(const std::string&p="", double=0.0,bool=true);
//   //   bool WriteFluGridASCII(const string &pre="", double=0.0,bool=true);
//   //   bool ReadFluSolnASCII(const std::string&p="",unsigned int=0, 
//   // 			double=0.0,bool=true);
//   //   bool WriteFluSolnASCII(const std::string&p="",unsigned int=0, 
//   // 			 double=0.0,bool=true);
//   //   bool WriteGSPASCII(const std::string&p="",unsigned int=0,
//   // 		     double=0.0,bool=true);
//   //   bool ReadGSPASCII(const std::string&p="",unsigned int=0,
//   // 		    double=0.0,bool=true);
//   //   bool ReadRegionASCII(const string &p="",unsigned int=0,unsigned int=0,
//   // 		       double=0.0,bool=true);
//   //   bool InitSurfaceSoln(const string &p="");
//   //   bool InitVolumeSoln(const string &p="");
//   //   bool CreateRegionMapFile(const string &p="",unsigned int=0,unsigned int=0);
//   //   void PopRemBordIndMPI();
//   //   bool PopRemBordIndFILE(const string &p="",double=0.0,bool=true);
//   //   bool ComRemeshInitData(const string &wname,double *cell_data,int nval_cells,
// // 			 double *node_data,int nval_nodes);
// //   bool ReadControlFile();
//   //  bool validate_comm_list(int ncsend,int ncrecv,int *csend,int *crecv);
//   //  bool PopulateFluVolumeSoln(double *data);
//   FloRegion(): GEM_Partition()
//   {
//     _casename = "";
//     _nnodes = 0;
//     _ni   = 0;
//     _nj   = 0;
//     _nk   = 0;
//     _ncells = 0;
//     _borders.resize(0);
//   };
  
//   FloRegion(const GEM_Partition &gp) : GEM_Partition(gp)
//   {
//     _casename = "";
//     _nnodes  = gp._nc.size()/3;
//     _ncells    = gp._hexconn.size()/8;
//   };
  
  
  
//   template<typename PB>
//   bool PopulateFloBorders(std::vector<PB> &pb)
//   {
//     unsigned int nborders = pb.size();
//     _borders.resize(nborders);
//     unsigned int border = 0;
//     while(border < nborders){
//       _borders[border]._rpart      = pb[border]._rpart;
//       _borders[border]._sendcells  = pb[border]._sendcells;
//       _borders[border]._recvcells  = pb[border]._recvcells;
//       _borders[border]._sharenodes = pb[border]._sharenodes;
//       _borders[border]._sendnodes  = pb[border]._sendnodes;
//       _borders[border]._recvnodes  = pb[border]._recvnodes;
//       _borders[border]._rbid       = pb[border]._rbid;
//       _borders[border]._out        = pb[border]._out;
//       border++;
//     }
//     return(true);
//   };

// //   bool
// //   validate_comm_list(int ncsend,int ncrecv,int *csend,int *crecv)
// //   {
// //     int index = 0;
// //     int nreal_cell = _tetconn.size()/4 + _prisconn.size()/6 +
// //       _pyrconn.size()/5 + _hexconn.size()/8 - (_ngtet + _ngpris +
// // 					       _ngpyr + _nghex);
// //     bool rval = true;
// //     while(index < ncsend){
// //       int ind = index++;
// //       if(!(csend[ind] <= nreal_cell)){
// // 	if(_out)
// // 	  *_out << "SEND CELL " << index << " is a ghost cell!!" << endl;
// // 	rval = false;
// //       }
// //       if(!(csend[ind] > 0)){
// // 	if(_out)
// // 	  *_out << "SEND CELL " << index << " is zero or negative!" << endl;
// // 	rval = false;
// //       }
// //     }
// //     index = 0;
// //     list<int> recvcell_list;
// //     while(index < ncrecv) {
// //       int ind = index++;
// //       if(!(crecv[ind] > nreal_cell)){
// // 	if(_out)
// // 	  *_out << "RECV CELL " << index << " is a real cell!!" << endl;
// // 	rval = false;
// //       }
// //       if(!(crecv[ind] > 0)){
// // 	if(_out)
// // 	  *_out << "RECV CELL " << index << " is zero or negative!" << endl;
// // 	rval = false;
// //       }
// //       bool duped = false;
// //       list<int>::iterator rci = recvcell_list.begin();
// //       while(rci != recvcell_list.end() && !duped){
// // 	if(crecv[ind] == *rci++){
// // 	  if(_out)
// // 	    *_out << "RECV CELL " << index 
// // 		  << " is duplicated in the receive list!" 
// // 		  << endl;
// // 	  duped = true;
// // 	}
// //       }
// //       if(!duped)
// // 	recvcell_list.push_back(crecv[ind]);
// //     }
// //     return(rval);
// //   };
//   void AddFluBorder(int rpid,int rpin,int nnshare, int nnsend,
// // 		    int nnrecv,int ncsend,int ncrecv,
// // 		    int *nshared,int *nsend,int *nrecv,
// // 		    int *csend,int *crecv){
// //     assert(rpid > 0 && rpin >= 0);
// //     if(_debug && _out)
// //       *_out << "FluRegion(" << _id << ")::AddFluBorder: Adding Border with"
// // 	    << " region " << rpid << "." << endl;
// //     FluBorder new_border;
// //     new_border._rbid = rpin;
// //     new_border._out = _out;
// //     if(!validate_comm_list(ncsend,ncrecv,csend,crecv)){
// //       if(_out)
// // 	*_out << "FluRegion(" << _id << ")::AddFluBorder: Validation of "
// // 	      << "communication arrays failed, aborting." << endl;
// //       exit(-1);
// //     }
// //     new_border.populate(rpid,nnshare,nnsend,nnrecv,ncsend,ncrecv,
// // 			nshared,nsend,nrecv,csend,crecv);
// //     _borders.push_back(new_border);
// //   };
// //   void AddFluPatch(int patch_id,int ntri, int ngtri, int *tris,
// // 		   int nquad,int ngquad, int *quads)
// //   {
// //     assert(ntri >= ngtri && nquad >= ngquad);
// //     if(_debug && _out)
// //       *_out << "FluRegion(" << _id << ")::AddFluPatch: Adding Patch with"
// // 	    << " id " << patch_id << "." << endl;
// //     FluPatch new_patch;
// //     new_patch._ntri = ntri;
// //     new_patch._nquad = nquad;
// //     new_patch._id = patch_id;
// //     new_patch._ngtri = ngtri;
// //     new_patch._ngquad = ngquad;
// //     int indy = 0;
// //     new_patch._triconn.resize(3*ntri);
// //     new_patch._quadconn.resize(4*nquad);
// //     new_patch._out = _out;
// //     while(indy < 3*ntri){
// //       assert(tris[indy] != 0);
// //       new_patch._triconn[indy] = tris[indy];
// //       assert(new_patch._triconn[indy] != 0);
// //       indy++;
// //     }
// //     indy = 0;
// //     while(indy < 4*nquad){
// //       assert(quads[indy] != 0);
// //       new_patch._quadconn[indy] = quads[indy];
// //       indy++;
// //     }
// //     _patches.push_back(new_patch);
// //   };
// // #ifdef _ROCSTAR_X_
// //   bool RegisterFluSurfaceMesh();
// //   bool RegisterSurfaceSoln(bool);
// //   bool RegisterVolumeSoln(bool);
// //#endif
// };

//istream &SkipLines(istream &,unsigned int = 1);
#endif






