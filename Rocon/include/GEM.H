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
#ifndef _GEM_H_
#define _GEM_H_
//#ifdef _ROCSTAR_X_
//#include "roccom.h"
//#endif
#include "mpi.h"
// ============================================================================
// Generic Interface Mesh Header
//
// The aim is to provide a simple mesh structure which should be easily
// populated from a number of working mesh representations. We will write
// translators to/from this simple/generic mesh representation to Rocstar 
// and other solver/code native mesh representations.
//=============================================================================
//
// Not yet used, but will be needed for quadratic elements
class GEM_Connectivity : public std::vector<unsigned int>
{
public:
  unsigned int _esize;
  GEM_Connectivity(const std::vector<unsigned int> &v)
    : std::vector<unsigned int>(v)
  {
    _esize = 0;
  };
};

template<class T>
class GEM_DataItem : public std::vector<T>
{
protected:
  std::string _name;
  int         _stride;
public:
  const std::string &Name() const { return(_name); };
  std::string &Name() {return(_name);};
  std::string &Name(const std::string &name) { return(_name.assign(name));};
  int Stride() const { return(_stride); };
  int Stride(int stride) { return(_stride = stride); };
  int &Stride() { return(_stride); };
  GEM_DataItem<T>(const std::string &name="",unsigned int size=0)
    : std::vector<T>(size)
  {
    _name.assign(name);
    _stride = 1;
  };
  GEM_DataItem<T>(const std::string &name,unsigned int size,const T &val)
    : std::vector<T>(size,val)
  {
    _name.assign(name);
    _stride = 1;
  };
  GEM_DataItem<T>(const std::vector<T> &vec,const std::string &name="",
		  int stride = 1)
    : std::vector<T>(vec)
  {
    _name.assign(name);
    _stride = stride;
  };
};
  
// To store field and integer data
class GEM_UserData
{
public:
  // Name of the dataset.  This data is used to identify
  // the window name in the case of Rocstar.
  std::string _name;
  // Example _field_data[0] = vector<double> of all elemental data
  //         _field_data[1] = vector<double> of all nodal data
  //         _field_data[2] = vector<double> of all facial data
  std::vector< GEM_DataItem<double> > _double_items;
  std::vector< std::vector<double> >  _field_data; 

  // Some codes have integer values that need to be associated with
  // cells and nodes.
  std::vector< GEM_DataItem<int> >    _int_items;
  std::vector< std::vector<int> >     _int_data;

  // Specifies the strides of the field and integer data
  std::vector<int>                    _stride_field;
  std::vector<int>                    _stride_int;
  std::vector<std::string>            _string_data;
  std::vector< GEM_DataItem<std::string> > _string_items;
  GEM_UserData()
  {
    _double_items.resize(0);
    _field_data.resize(0);
    _int_items.resize(0);
    _int_data.resize(0);
    _stride_field.resize(0);
    _stride_int.resize(0);
    _string_data.resize(0);
    _string_items.resize(0);
  };
};

//
// DomainBoundary = the physical boundary of the computational domain
// This structure describes the piece of the DomainBoundary that might exist
// on a given partition.
//
class GEM_DomainBoundary {
public:
  unsigned int _id;     // User data - but usually used to indicate unique BC
  unsigned int _ngtri;  // Number of ghost triangles
  unsigned int _ngquad; // Number of ghost quads 
  unsigned int _nnodes;
  std::vector<unsigned int> _triconn;  // Triangle connectivity
  std::vector<unsigned int> _quadconn; // Quad connectivity
  std::vector<double> surface_coordinates;
  unsigned int   surface_ngnodes;
  std::vector<unsigned int> surface_tri;
  std::vector<unsigned int> surface_quad;
  GEM_UserData         _data;
  GEM_UserData         _solver_data;
  bool                 _debug;
  std::ostream *_out;
  unsigned int NNodes(void);
  bool SetSolverDataBlock(const std::string &wname,
			  double *cell_data,
			  int nval_cells, 
			  double *node_data,int nval_nodes);
  // Rocflu Functions
  //  bool InitFluSurfaceSoln(const string &prefix="",bool all = true);
  //  bool RegisterFluSoln(const string &wname,bool all = true);
  //#ifdef _ROCSTAR_X_
  int pane_id;
  //#endif  
  GEM_DomainBoundary(const GEM_DomainBoundary &db)
  {
    _id       = db._id;
    _out      = db._out;
    _ngtri    = db._ngtri;
    _ngquad   = db._ngquad;
    _triconn  = db._triconn;
    _quadconn = db._quadconn;
    _debug    = db._debug;
    _data     = db._data;
    _nnodes   = db._nnodes;
  };
  GEM_DomainBoundary()
  {
    _id = 0;
    _out = NULL;
    _ngtri = 0;
    _ngquad = 0;
    _triconn.resize(0);
    _quadconn.resize(0);
    _nnodes = 0;
    //    _data.resize(0);
    _debug = false;
    surface_coordinates.resize(0);
  };
  GEM_DomainBoundary &operator=(const GEM_DomainBoundary &db)
  {
    _id       = db._id;
    _out      = db._out;
    _ngtri    = db._ngtri;
    _ngquad   = db._ngquad;
    _triconn  = db._triconn;
    _quadconn = db._quadconn;
    _data     = db._data;
    _debug    = db._debug;
    _nnodes   = db._nnodes;
    surface_coordinates = db.surface_coordinates;
    surface_tri = db.surface_tri;
    surface_quad = db.surface_quad;
    return(*this);
  };
  bool debug(bool s = true)
  {
    _debug = s;
    return(_debug);
  };
  void report();
  // Creates a nodal coordinate array [x y z x y z...] with ghosts
  // at the end and connectivity arrays that refer to the new nc.
  // This is just a convenience function to build those surface local arrays
  // that are required by Roccom.
  void PopulateSurfaceArrays(const std::vector<double> &,unsigned int);
  //#ifdef _ROCSTAR_X_
  bool Register_com_surfmesh(const std::string &wname);
  void Create_com_surfsoln(const std::string &wname,const std::string &fname,
			   std::vector<double> &fvec, unsigned int ncomp,
			   const std::string &unit);
  //#endif
};

//
// PartitionBoundary = non physical boundary introduced by domain partitioning.
// This structure describes a boundary for which information must be exchanged
// with a neighboring partition.
//
// Notes:
// The _sendcells and _recvcells refer to a partition cell mapping following 
// the convention of GEM_Partition::_cell_ordering[].  
//  
// It is up to the user to make sure that the ordering appears here matches
// that of the sending and receiving partitions.
//
class GEM_PartitionBoundary {
public:
  unsigned int _rpart;  // Remote partition id (1 based)
  std::vector<unsigned int> _sendcells;  // Real Cells associated with this PB
  std::vector<unsigned int> _recvcells;  // Ghost Cells associated with this PB
  std::vector<unsigned int> _sendnodes;  // Real nodes associated with this PB
  std::vector<unsigned int> _recvnodes;  // Ghost nodes associated with this PB
  std::vector<unsigned int> _sharenodes; // Nodes *on* the PB
  GEM_UserData         _data;
  std::ostream        *_out;
  bool                 _debug;
  void report();
  void populate(int rpid, int nnshared,int nnsend,int nnrecv,
		int ncsend,int ncrecv, int *sharedn, int *sendn,
		int *recvn, int *sendc,int *recvc);
  GEM_PartitionBoundary(const GEM_PartitionBoundary &pb)
  {
    _rpart      = pb._rpart;
    _out        = pb._out;
    _sendcells  = pb._sendcells;
    _recvcells  = pb._recvcells;
    _sendnodes  = pb._sendnodes;
    _recvnodes  = pb._recvnodes;
    _sharenodes = pb._sharenodes;
    _data       = pb._data;
  };
  GEM_PartitionBoundary()
  {
    _rpart = 0;
    _out = NULL;
    _sendcells.resize(0);
    _sendnodes.resize(0);
    _recvnodes.resize(0);
    _recvcells.resize(0);
    _sharenodes.resize(0);
    //    _data.resize(0);
  };
  GEM_PartitionBoundary &operator=(const GEM_PartitionBoundary &pb)
  {
    _rpart      = pb._rpart;
    _out        = pb._out;
    _sendcells  = pb._sendcells;
    _recvcells  = pb._recvcells;
    _sendnodes  = pb._sendnodes;
    _recvnodes  = pb._recvnodes;
    _sharenodes = pb._sharenodes;
    _data       = pb._data;
    return(*this);
  };
};


//
// Partition = A part of a mesh representation of some computational domain
// This structure completely describes a partition and provides the translator
// functionality.  The design choice was to make the translator functions
// member functions of this class (for now).
//
class GEM_Partition {
public:
  unsigned int                  _id;      // Unique id (should be 1 based)
  unsigned int                  _npart;   // Total number of partitions
  unsigned int                  _ngnodes; // Number of ghost nodes
  unsigned int                  _ngtet;   // Number of ghost tets
  unsigned int                  _nghex;   // Number of ghost hex
  unsigned int                  _ngpyr;   // Number of ghost pyr
  unsigned int                  _ngpris;  // Number of ghost pris
  unsigned int                  _nvface;  // Number volume faces (r-r,r-g)
  // Nodal Coordinates [x y z x y z...]
  std::vector<double>           _nc;
  // Tet connectivity [1 2 3 4 1 2 ...]
  std::vector<unsigned int>     _tetconn; 
  // Pyr conn [1 2 3 4 5 1 2 3 ...]
  std::vector<unsigned int>     _pyrconn; 
  // Pris conn[1 2 3 4 5 6 1 2 3 ...] 
  std::vector<unsigned int>     _prisconn;
  // Hex conn[1 2 3 4 5 6 7 8 1 2 3 ..]
  std::vector<unsigned int>     _hexconn;
  // 1234=tet,pyr,pris,hex
  unsigned int                  _cell_ordering[4];
  std::vector<GEM_PartitionBoundary> _pb;
  std::vector<GEM_DomainBoundary>    _db;
  GEM_UserData                  _data; // any volume solution data
  GEM_UserData                  _solver_data; // solver specific data
  std::ostream *_out;
  //  vector<std::pair<unsigned int,unsigned int> > _cellmap;
  bool                          _debug;
  //#ifdef _TRAIL_MPI_
  MPI_Comm _comm;
  //#endif
  //#ifdef _ROCSTAR_X_
  int                                pane_id;
  std::vector<int>                   pconn;
  unsigned int                       pconn_nghost;
  std::string                        volume_window;
  std::string                        surface_window;
  //#endif
  // Functions for Rocflu
  //  bool ReadFluControlFile();
  //  bool WriteFluNative();
  //  bool WriteFluCOM(const string &pre = "");
  //  bool WriteFluDIM(const string &pre = "",double t=0,bool unsteady=true);
  //  bool WriteFluCMP(const string &pre = "");
  //  bool ReadFluCOM(const string &pre = "");
  //  bool WriteFluMAP(const string &pre = "",unsigned int nproc = 1,
  //		   unsigned int nreg = 1);
  //  bool ReadFluDIM(const string &prefix,double t,bool unsteady);
  //  bool PopRemBordIndFILE(const string &p="",double=0.0,bool=true);
  //  bool InitFluVolumeSoln(const string &prefix = "");
  //  bool InitFluSurfaceSoln(const string &prefix = "");
  //  bool BuildFluPatchMapping(map<unsigned int,unsigned int> &,
  //			    const string &p="");
  //  bool RegisterFluVolumeSoln(bool all = true);
  //  bool RegisterFluSurfaceSoln(bool all = true);
  //  bool RegisterFluSurfaceMesh();
  //  bool PopulateFluPatches(const string &pre,vector<GEM_DomainBoundary> &db);
  bool validate_comm_list(int ncsend,int ncrecv,int *csend,int *crecv);
  bool ValidateMesh();
  void AddParitionBoundary(int rpid,int nnshare, int nnsend,
			   int nnrecv,int ncsend,int ncrecv,
			   int *nshared,int *nsend,int *nrecv,
			   int *csend,int *crecv);

  void AddDomainBoundary(int db_id,int ntri, int ngtri, int *tris,
			 int nquad,int ngquad, int *quads);
  bool SetSolverDataBlock(const std::string &wname,double *cell_data,
			  int nval_cells,double *node_data,
			  int nval_nodes);
  bool AddSolverDataBlock(const std::string &wname,double *cell_data,
			  int nval_cells,double *node_data,
			  int nval_nodes);
  bool debug(bool s = true);
  void ResolveCellMapping(GEM_Partition &sp);
  void SetNodalCoordinates(double *data,int nn,int ng);
  void SetVolumeElements(int *data,int ncells,int ng,int npe);
  bool PopulatePartitionBoundaries(std::vector<GEM_PartitionBoundary> &pb);
  void MapDomainBoundaries(std::map<unsigned int,unsigned int> &bcmap);
  // Dumps a description of the partition to stdout
  void report(); 
  void report_partition_boundaries();
  void report_domain_boundaries();
  // Returns total number of elements
  unsigned int nelem() 
  {
    return(_tetconn.size()/4+_pyrconn.size()/5+
	   _prisconn.size()/6+_hexconn.size()/8);
  };
  // The following functions deal with codes that map their cells in some 
  // other order besides the default. The default cell mapping goes like:
  // real tets  : cells 1 - N
  // real pyrs  : cells N+1 - N'
  // real pris  : cells N'+1 - N''
  // real hex   : cells N''+1 - nreal_elem
  // ghost tets : cells nreal_elem + 1 - M
  // ghost pyrs : cells M+1 - M'
  // ghost pris : cells M'+1 - M''
  // ghost hex  : cells M''+1 - (nelem_total)
  // This mapping is realized by using the following two functions:
  // Elem2Cell - takes a key-value element pair (explained below), and returns
  // the cell id. 
  // Example: Elem2Cell(<1,2>) = 2  (cell 2 is tet 2)
  unsigned int Elem2Cell(std::pair<unsigned int,unsigned int>);
  // Cell2Elem - will return the key-value pair of the element type, and 
  // id.  Element type is always defined like:
  //      1 = tet, 2 = pyr, 3 = pris, 4 = hex
  // Example: <1,4> = the 4th tet element   
  std::pair<unsigned int,unsigned int> Cell2Elem(unsigned int);
  
  //#ifdef _ROCSTAR_X_
  void Register_com_volconn(const std::string &wname,int paneid,unsigned int nel,
			    unsigned int ngel,std::vector<unsigned int> &conn,
			    unsigned int esize,bool ghost_part = false);
  
  bool WindowInitDone();
  bool DestroyWindows();
  // Write a Roccom window to whatever format is supported by Rocstar 
  // (CGNS, HDF4, etc) to the specified directory.
  bool WriteRocstar(const std::string &,double t = 0.0);
  bool ReadRocstar(const std::string &,double t = 0.0);
  bool InitRoccomWindows(const std::string &wname);
  void Create_com_pconn(std::vector<unsigned int> rpids,
			std::vector< std::vector< 
			std::vector<unsigned int> > > &index_vectors,
			unsigned int &nreal,unsigned int &ng,
			std::vector<int> &pc);
  bool CreatePconn(const std::string &wname);
  bool PopulateVolumeWindow(const std::string &wname);
  bool PopulateSurfaceWindow(const std::string &wname);
  //  bool ComRemeshInitData(const std::string &wname,double *cell_data,
  //			 int nval_cells, double *node_data,int nval_nodes);
  //  bool ComRemeshInitFluData(const std::string &wname,double *cell_data,
  //			    int nval_cells, double *node_data, int nval_nodes);
  //  bool ComRemeshInitFluData();
  void Create_com_volsoln(const std::string &fname,std::vector<double> &fvec,
			  unsigned int ncomp,const std::string &unit);
  //#endif
  GEM_Partition(const GEM_Partition &gp)
  {
    _out      = gp._out;
    _id       = gp._id;
    _ngnodes  = gp._ngnodes;
    _ngtet    = gp._ngtet;
    _nghex    = gp._nghex;
    _ngpyr    = gp._ngpyr;
    _ngpris   = gp._ngpris;
    _nc       = gp._nc;
    _tetconn  = gp._tetconn;
    _pyrconn  = gp._pyrconn;
    _prisconn = gp._prisconn;
    _hexconn  = gp._hexconn;
    _nvface   = gp._nvface;
    _debug    = gp._debug;
    _cell_ordering[0] = gp._cell_ordering[0];
    _cell_ordering[1] = gp._cell_ordering[1];
    _cell_ordering[2] = gp._cell_ordering[2];
    _cell_ordering[3] = gp._cell_ordering[3];
    //#ifdef _TRAIL_MPI_
    _comm = gp._comm;
    //#endif
  };
  GEM_Partition()
  {
    _id       = 0;
    _out      = NULL;
    _ngnodes  = 0;
    _ngtet    = 0;
    _nghex    = 0;
    _ngpyr    = 0;
    _ngpris   = 0;
    _nvface   = 0;
    _nc.resize(0);
    _debug    = false;
    _tetconn.resize(0);
    _pyrconn.resize(0);
    _prisconn.resize(0);
    _hexconn.resize(0);
    _cell_ordering[0] = 1;
    _cell_ordering[1] = 2;
    _cell_ordering[2] = 3;
    _cell_ordering[3] = 4;
    //#ifdef _TRAIL_MPI_
    _comm = MPI_COMM_NULL;
    //#endif
  };
  //#ifdef _TRAIL_MPI_
  GEM_Partition(unsigned int id,MPI_Comm comm = MPI_COMM_NULL)
    : _id(id), _comm(comm)
  {
    _out      = NULL;
    _ngnodes  = 0;
    _ngtet    = 0;
    _nghex    = 0;
    _ngpyr    = 0;
    _ngpris   = 0;
    _nvface   = 0;
    _nc.resize(0);
    _debug    = false;
    _tetconn.resize(0);
    _pyrconn.resize(0);
    _prisconn.resize(0);
    _hexconn.resize(0);
    _cell_ordering[0] = 1;
    _cell_ordering[1] = 2;
    _cell_ordering[2] = 3;
    _cell_ordering[3] = 4;
  };
  //#endif
  GEM_Partition &operator=(const GEM_Partition &gp)
  {
    _id       = gp._id;
    _out      = gp._out;
    _ngnodes  = gp._ngnodes;
    _ngtet    = gp._ngtet;
    _nghex    = gp._nghex;
    _ngpyr    = gp._ngpyr;
    _ngpris   = gp._ngpris;
    _nc       = gp._nc;
    _nvface   = gp._nvface;
    _tetconn  = gp._tetconn;
    _pyrconn  = gp._pyrconn;
    _prisconn = gp._prisconn;
    _hexconn  = gp._hexconn;
    _debug    = gp._debug;
    //#ifdef _TRAIL_MPI_
    _comm     = gp._comm;
    //#endif
    _cell_ordering[0] = gp._cell_ordering[0];
    _cell_ordering[1] = gp._cell_ordering[1];
    _cell_ordering[2] = gp._cell_ordering[2];
    _cell_ordering[3] = gp._cell_ordering[3];
    return(*this);
  };
};

// Utilities
bool
flip_elements(std::vector<unsigned int> &,unsigned int);

#endif 






