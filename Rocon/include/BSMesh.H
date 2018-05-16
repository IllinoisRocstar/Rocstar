///
/// \file
/// \ingroup support
/// \brief Block Structured Mesh stuff
///
#ifndef _BSMESH_H_
#define _BSMESH_H_
#include "Mesh.H"
namespace Mesh {
  ///
  /// \brief Simple Block Structured Mesh object
  ///
  /// Block structured mesh object holds irange, jrange, and krange.
  /// Elements are ordered in dictionary ordering of lower left node
  template<typename T>
  class BSExtent : public std::vector<std::vector<T> > 
  {
  private:
    size_t _nd;
    std::vector<size_t> _Np;
    std::vector<size_t> _N;
  public:
    BSExtent()
    {
      this->resize(0);
      _nd = 0;
      _Np.resize(0);
      _N.resize(0);
    }
    BSExtent(const std::vector<std::vector<T> > &inextent)
    {
      _nd = inextent.size();
      this->resize(_nd); 
      for(size_t i = 0; i < _nd; i++){
	typename std::vector<T>::const_iterator iei = inextent[i].begin();
	while(iei != inextent[i].end())
	  (*this)[i].push_back(*iei++);
      }     
      Sync();
    };
    BSExtent(const std::vector<T> &inflatextent){
      _nd = inflatextent.size()/2;
      this->resize(_nd);
      size_t nd = 0;
      typename BSExtent<T>::iterator bsi = this->begin();
      while(bsi != this->end()){
	bsi->push_back(inflatextent[nd++]);
	bsi->push_back(inflatextent[nd++]);
	bsi++;
      }
      Sync();
    };
    BSExtent(const T *src,int nd=3){
      _nd = nd;
      this->resize(_nd);
      unsigned int nindex = 0;
      typename BSExtent<T>::iterator bsi = this->begin();
      while(bsi != this->end()){
	bsi->push_back(src[nindex++]);
	bsi->push_back(src[nindex++]);
	bsi++;
      }
      Sync();
    };
    void destroy()
    {
      typename BSExtent<T>::iterator bsi = this->begin();
      while(bsi != this->end())
	bsi++->resize(0);
      this->resize(0);
      _nd = 0;
      _Np.resize(0);
      _N.resize(0);
    };
    void Init(const std::vector<T> &inflatextent){
      destroy();
      _nd = inflatextent.size()/2;
      this->resize(_nd);
      size_t nd = 0;
      typename BSExtent<T>::iterator bsi = this->begin();
      while(bsi != this->end()){
	bsi->push_back(inflatextent[nd++]);
	bsi->push_back(inflatextent[nd++]);
	bsi++;
      }
      Sync();
    };
    void Sync()
    {
      _nd = this->size();
      _N.resize(_nd);
      _Np.resize(_nd);
      for(size_t i = 0; i < _nd; i++){
	_N[i] = (*this)[i][1] - (*this)[i][0] + 1;
      }     
      _Np[0] = 1;
      for(size_t i = 1; i < _nd;i++){
	_Np[i] = _N[i-1] * _Np[i-1];
      }
    };
    T NNodes(){
      T nnodes = 1;
      typename BSExtent<T>::iterator bsi = this->begin();
      while(bsi != this->end()){
	nnodes *= (((*bsi)[1] - (*bsi)[0])+1);
	bsi++;
      }
      return(nnodes);
    };
    int ND(){return(this->size());};
    void Flatten(std::vector<T> &output)
    {
      output.resize(0);
      typename BSExtent<T>::iterator bsi = this->begin();
      while(bsi != this->end()){
	output.push_back((*bsi)[0]);
	output.push_back((*bsi)[1]);
	bsi++;
      }
    };
    void dir_loop(int nd,T indoff,T plane,std::vector<T> &N,
		  std::vector<T> &NP,const BSExtent<T> &inextent,
		  std::vector<T> &indices)
    {
      int dindex = nd - 1;
      plane += NP[dindex];
      T start = inextent[dindex][0];
      T end   = inextent[dindex][1];
      for(T dd = start;dd <= end;dd++){
	if(dindex == 0)
	  indices.push_back(dd - (*this)[dindex][0]+1+indoff);
	else{
	  T ndoff = indoff + (dd - (*this)[dindex][0])*NP[dindex];
	  dir_loop(dindex,ndoff,plane,N,NP,inextent,indices);
	}
      }
	
    };
    void GetFlatIndices(const BSExtent<T> &extent,std::vector<T> &indices)
    {
      int nd = ND();
      indices.resize(0);
      std::vector<T> N(nd);
      std::vector<T> NP(nd);
      for(int nc = 0;nc < nd;nc++){
	N[nc] = (*this)[nc][1] - (*this)[nc][0] + 1;
	NP[nc] = (nc == 0 ? 1 : N[nc-1] * NP[nc-1]);
      }
      NP[0] = 0;
      dir_loop(nd,0,0,N,NP,extent,indices);
    }; 
    void Overlap(const BSExtent<T> &inextent,BSExtent<T> &outextent)
    {
      size_t nd = ND();
      outextent.resize(0);
      bool match = true;
      for(size_t j = 0;j < nd && match;j++){
	if(!(((*this)[j][0] >= inextent[j][0] &&
	      (*this)[j][0] <= inextent[j][1]) ||
	     ((*this)[j][1] >= inextent[j][0] &&
	      (*this)[j][1] <= inextent[j][1])))
	  match = false;
      }
      if(match){ // then this extent has some overlap with the searched
	outextent.resize(nd);
	for(size_t i = 0;i < nd;i++)
	  outextent[i].resize(2);
	for(size_t i = 0;i < nd;i++){
	  outextent[i][0] = std::max((*this)[i][0],inextent[i][0]);
	  outextent[i][1] = std::min((*this)[i][1],inextent[i][1]);
	}
      }
    };
    void FindSharedNodes(std::vector<BSExtent<T> > &extent_pool,std::vector<BSExtent<T> > &shared_extents,
			 std::vector<T> &neighbors)
    {
      shared_extents.resize(0);
      T nindex = 0;
      typename std::vector<BSExtent<T> >::iterator epi = extent_pool.begin();
      while(epi != extent_pool.end()){
	if(*epi != *this){
	  BSExtent<T> overlap;
	  this->Overlap(*epi,overlap);
	  if(!overlap.empty()){
	    neighbors.push_back(nindex);
	    shared_extents.push_back(overlap);
	  }
	}
	nindex++;
	epi++;
      }
    };
    inline T NodeNum(std::vector<T> &index)
    {
      T l0 = 1;
      for(size_t i = (_nd-1);i >= 0;i--)
	l0 += (index[i] - (*this)[i][0])*_Np[i]; 
    };
    void CreateUnstructuredMesh(Connectivity &conn)
    {
      // get coordinates of lower left hand coordinates:
      // istart to (iend - 1), or istart
      // jstart to (jend - 1), or jstart
      // kstart to (kend - 1), or kstart
      assert(_nd == 3);
      BSExtent<T> lowerleft(*this);
      int mesh_nd = 3;
      for(int i = 0;i < 3;i++){
	if(_N[i] > 1)
	  lowerleft[i][1]--;
	else
	  mesh_nd--;
      }
      assert(mesh_nd > 1);
      std::vector<T> indices;
      GetFlatIndices(lowerleft,indices);
      IndexType offset = (_N[0] > 1 ? _N[0] : _N[1]);
      IndexType planeoffset = (mesh_nd == 2 ? 0 : offset*_N[1]);
      IndexType number_of_elements = indices.size();
      for(IndexType i = 0;i < number_of_elements;i++){
	std::vector<IndexType> element;
	element.push_back(indices[i]);
	element.push_back(indices[i]+1);
	element.push_back(indices[i]+1+offset);
	element.push_back(indices[i]+offset);
	if(planeoffset > 0){
	  element.push_back(indices[i]+planeoffset);
	  element.push_back(indices[i]+1+planeoffset);
	  element.push_back(indices[i]+planeoffset+offset+1);
	  element.push_back(indices[i]+planeoffset+offset);
	}
	conn.AddElement(element);
      }
      conn.Sync();
    };
  };

  int GenerateCartesianGrid(Mesh::NodalCoordinates &nc,
			    Mesh::BSExtent<Mesh::IndexType> &gridextent,
			    std::vector<Mesh::IndexType> &gridsizes,
			    GeoPrim::CBox &box);

}
#endif
