///
/// \file
/// \ingroup support
/// \brief Parallel Mesh implementation
///
#include <iomanip>
#include <sstream>

#include "PMesh.H"

namespace Mesh {

  int Partition::GetBorderElements(std::vector<Mesh::IndexType> &be) const
  {
    std::list<Mesh::IndexType> belist;
    Mesh::IndexType number_of_nodes = _mesh.NumberOfNodes();
    Mesh::IndexType number_of_elements = _mesh.NumberOfElements();
    std::vector<bool> element_processed(number_of_elements,false);
    Mesh::IndexType nvolume_nodes = number_of_nodes - _info.nshared;
    const Mesh::Connectivity &dc = _mesh.GetCon(MeshUtilityObject::ENDUAL);
    for(Mesh::IndexType nn = nvolume_nodes;nn < number_of_nodes;nn++){
      std::vector<Mesh::IndexType>::const_iterator ei = 
	dc[nn].begin();
      while(ei != dc[nn].end()){
	if(!element_processed[*ei-1]){
	  element_processed[*ei-1] = true; // cheap uniquity
	  belist.push_back(*ei);
	}
	ei++;
      }
    }
    belist.sort();
    std::list<Mesh::IndexType>::iterator beli = belist.begin();
    be.resize(belist.size());
    std::vector<Mesh::IndexType>::iterator bei = be.begin();
    while(bei != be.end())
      *bei++ = *beli++;
    return(0);
  };

  int Partition::Read(const std::string &MeshName,IRAD::Comm::CommunicatorObject &comm,bool allow_n2m,std::ostream &ErrOut)
  {
    unsigned int nproc = comm.Size();
    unsigned int rank = comm.Rank();
    unsigned int id = rank + 1;
    // Read the partitioning info
    if(nproc > 1){
      std::ifstream InfoInf;
      std::ostringstream FNOstr;
      FNOstr << MeshName << "." << id << ".info";
      InfoInf.open(FNOstr.str().c_str());
      if(!InfoInf){
	comm.SetErr(1);
      }
      if(comm.Check()){
	ErrOut << "Partition::Read Could not find partition info for " << MeshName;
	return(1);
      }
      InfoInf >> _info.npart >> _info.part >> _info.nelem >> _info.nnodes >> _info.nborder 
	      >> _info.nshared >> _info.nown;
      if(_info.npart != nproc && !allow_n2m){
	comm.SetErr(1);
      }
      if(comm.Check()){
	ErrOut << "Partition::Read " << MeshName << " has partition/processor mismatch: ("
	       << _info.npart << "/" << nproc << ").\n";
	return(1);
      }
      _info.nlocal = _info.nnodes  - _info.nshared + _info.nown;
      InfoInf.close();
    }
    else{
      _info.nborder = 0;
      _info.nshared = 0;
    }
    // done reading partition info
    std::ifstream Inf;
    std::ostringstream Ostr;
    Ostr << MeshName;
    if(nproc > 1)
      Ostr << "." << id << ".pmesh";
    else
      Ostr << ".mesh";
    Inf.open(Ostr.str().c_str());
    // Don't blink or you'll miss the mesh actually being loaded here.
    Inf >> _mesh.NC() >> _mesh.ECon();
    _mesh.ECon().ShrinkWrap();
    //    Inf >> _nc >> _ec;
    //    _ec.ShrinkWrap();
    // done loading actual mesh  
    //    Mesh::IndexType number_of_nodes    = _nc.Size();
    //    Mesh::IndexType number_of_elements = _ec.Nelem();
    Mesh::IndexType number_of_nodes    = _mesh.NumberOfNodes();
    Mesh::IndexType number_of_elements = _mesh.NumberOfElements();
    if(nproc == 1){
      _info.nnodes = number_of_nodes;
      _info.nown   = 0;
      _info.nlocal = number_of_nodes;
      _info.nelem  = number_of_elements;
    }
    else{
      Inf >> _info.nborder;
      _borders.resize(_info.nborder);
      for(Mesh::IndexType nn = 0;nn < _info.nborder;nn++){
	Mesh::IndexType nrecv = 0;
	Mesh::IndexType nsend = 0;
	Inf >> _borders[nn].rpart >> nrecv >> nsend;
	_borders[nn].nrecv.resize(nrecv);
	for(Mesh::IndexType ii = 0;ii < nrecv;ii++){
	  Mesh::IndexType nodeid = 0;
	  Inf >> _borders[nn].nrecv[ii] >> nodeid;
	}
	_borders[nn].nsend.resize(nsend);
	for(Mesh::IndexType ii = 0;ii < nsend;ii++){
	  Mesh::IndexType nodeid = 0;
	  Inf >> _borders[nn].nsend[ii] >> nodeid;
	}
      }
    }
    Inf.close();
    _communicator = &comm;
    return(0);
  }


  int
  Partition::Read(const std::string &MeshName,int id)
  {
    // rank = id - 1
    std::ifstream Inf;
    std::ostringstream FNOstr;
    if(id > 0){
      FNOstr << MeshName << "." << id << ".info";
      Inf.open(FNOstr.str().c_str());
      if(!Inf)
	return(2);
      Inf >> _info.npart >> _info.part >> _info.nelem >> _info.nnodes >> _info.nborder >> _info.nshared >> _info.nown;
      Inf.close();
      _info.nlocal = _info.nnodes - _info.nshared + _info.nown;
      FNOstr.str("");
    }
    _info.nborder = 0;
    _info.nshared = 0;
    _info.nown = 0;
    if(id > 0)
      FNOstr << MeshName << "." << id << ".pmesh";
    else 
      FNOstr << MeshName;
    Inf.open(FNOstr.str().c_str());
    if(!Inf)
      return(3);
    //    Inf >> _nc;
    //    Inf >> _ec;
    //    _ec.ShrinkWrap();
    //    _ec.Sync();
    Inf >> _mesh.NC();
    Inf >> _mesh.ECon();
    _mesh.ECon().ShrinkWrap();
    _mesh.ECon().Sync();
    if(id > 0){
      Inf >> _info.nborder;
      _borders.resize(_info.nborder);
      for(Mesh::IndexType nn = 0;nn < _info.nborder;nn++){
	Mesh::IndexType nrecv = 0;
	Mesh::IndexType nsend = 0;
	Inf >> _borders[nn].rpart >> nrecv >> nsend;
	for(Mesh::IndexType ii = 0;ii < nrecv;ii++){
	  Mesh::IndexType nodeid = 0;
	  Inf >> nodeid;
	  _borders[nn].nrecv.push_back(nodeid);
	  Inf >> nodeid;
	}
	for(Mesh::IndexType ii = 0;ii < nsend;ii++){
	  Mesh::IndexType nodeid = 0;
	  Inf >> nodeid;
	  _borders[nn].nsend.push_back(nodeid);
	  Inf >> nodeid;
	}
      }
    }
    Inf.close();
    return(0);
  };
  
}
