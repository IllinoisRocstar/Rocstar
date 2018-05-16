/**
 \file 
 \ingroup support
 \brief Parallel Mesh Stuff

 All mesh data structures (MDS) will follow the CGNS conventions:\n
 \verbatim
 Triangle Edges: 12 23 31    (counterclockwise)
 Square Edges:   12 23 34 41 (ccw)
 Tet Faces:      132 241 342 143
 Pyr Faces:      1432 251 352 453 154
 Prism Faces:    2541 3652 1463 132 456
 Hex Faces:      1432 2651 3762 4873 1584 5678 
 \endverbatim 
 \image latex elements.eps "Element Types"
 \image html elements.jpg "Element Types" width=10cm
*/
/*
                  Supported Cell Types

        4                                              6   
       /|\             5       8---------7            /|\ 
      / | \         /  |       | \       |\          / | \
     /  |  \       /  /        |  \      | \        /  |  \
    /   |   \     /   |        |   5--------\6     /   3   \
   1----|----3   4----|-3      4---|-----3   |    4---------5
    \   |   /     \  /   \      \  |      \  |    |  /   \  |
     \  |  /       \ |    \      \ |       \ |    | /     \ |
      \ | /         \|     \      \|        \|    |/       \|
        2            1------2      1---------2    1---------2
 
*/
#ifndef _PMESH_H_
#define _PMESH_H_
#include "Mesh.H"
#include "COMM.H"

namespace Mesh {  
  /** 
      \brief Border data in CSR format
      \verbatim
   This is the border data in CSR format.  Review of CSR:
   Ap: is the array of sizes (in number of nonzero entries)
       in the Sparse Row.  Example: 
       [0 5 14 27 84 85 100] (5 in first row (row 0) in Ai[0:4]) 
                             (9 in row[1], Ai[5:13]), etc        
  
   Ai: is the array containing the column numbers of the nz 
       entries in the sparse row.  Example: 
       [0 3 4 7 8 1 5 ......] (row 0): X00XX00XX000.....
                              (row 1): 0X000X000000.....
  
   Ax: would store the actual values of the sparse matrix 
  
   Each row corresponds to a border *node* - so the actual send and
   receive arrays are puffed up to contain a number of matrix values
   required for nremote_dof * nlocal_dof.
     \endverbatim
  */
  struct BorderData {
    /// \brief sizes in ndof sent
    std::vector<Mesh::IndexType> SendAp; 
    /// \brief dof id's of the sent dofs
    std::vector<Mesh::IndexType> SendAi;
    /// \brief sizes in ndof recv
    std::vector<Mesh::IndexType> RecvAp;
    /// \brief dof id's of recvd dofs
    std::vector<Mesh::IndexType> RecvAi; 
    std::vector<Mesh::IndexType> border_dof_to_global;
    std::vector<double> SendBuffer;
    std::vector<double> RecvBuffer;
  };
  struct Border {
    Mesh::IndexType rpart; /// remote partition id
    Mesh::IndexVec nsend;  /// local id's of remote nodes
    Mesh::IndexVec nrecv;  /// local id's of local nodes
    Mesh::IndexType recvsize; /// total num owned dofs
    Mesh::IndexType sendsize; /// hrm
    Mesh::IndexType remote_dofcount; /// total num remote dofs
    BorderData data; /// csr support/comm buffers
  };

  struct PartInfo {
    Mesh::IndexType npart;   /// total number of partitions
    Mesh::IndexType part;    /// partition id
    Mesh::IndexType nelem;   /// total number of elems
    Mesh::IndexType nnodes;  /// total number of nodes
    Mesh::IndexType nborder; /// number of borders
    Mesh::IndexType nshared; /// number of shared nodes
    Mesh::IndexType nown;    /// number of shared nodes owned
    Mesh::IndexType nlocal;  /// number of local nodes (convenience)
    Mesh::IndexType doffset; /// global dof offset
  };
  
  class Partition {
  private:
    IRAD::Comm::CommunicatorObject *_communicator;
  public:
    Partition() : _communicator(NULL), _mesh(),_info(),_borders() {};
    MeshUtilityObject _mesh;
    //    NodalCoordinates    _nc;
    //    Connectivity        _ec;
    PartInfo            _info;
    std::vector<Border> _borders;
    int Read(const std::string &meshname,int id);
    int Read(const std::string &MeshName,IRAD::Comm::CommunicatorObject &comm,
	     bool allow_n2m,std::ostream &ErrOut);
    //    IRAD::Comm::CommunicatorObject *CommunicatorPtr(){ return(_communicator); };
    //    std::ostream &Report(std::ostream &Ostr);
    MeshUtilityObject &Mesh() { return (_mesh); };
    const MeshUtilityObject &Mesh() const { return (_mesh); };
    PartInfo &Info() { return (_info); };
    const PartInfo &Info() const { return(_info); };
    std::vector<Border> &Borders() { return(_borders); };
    const std::vector<Border> &Borders() const { return(_borders); };
    IRAD::Comm::CommunicatorObject *CommunicatorPtr(){ return(_communicator); };
     //    std::ostream &Report(std::ostream &Ostr);
    int GetBorderElements(std::vector<Mesh::IndexType> &be) const;
  };
  
}
#endif
