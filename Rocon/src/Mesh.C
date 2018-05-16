///
/// \file
/// \ingroup support
/// \brief Mesh stuff implementation
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include "Mesh.H"

namespace Mesh {

  Mesh::NodalCoordinates::NodalCoordinates()
  {
    ncdata = NULL;
    nnodes = 0;
    mydata = false;
  };
  Mesh::NodalCoordinates::NodalCoordinates(IndexType n)
  {
    if(n > 0){
      ncdata = new double [3 * n];
      nnodes = n;
      mydata = true;
    }
    else {
      ncdata = NULL;
      nnodes = 0;
      mydata = false;
    }
  };
  Mesh::NodalCoordinates::NodalCoordinates(IndexType n,double *data)
  {
    if((n > 0) && (data != NULL)){
      mydata = false;
      nnodes = n;
      ncdata = data;
    }
    else {
      ncdata = NULL;
      nnodes = 0;
      mydata = false;
    }
  };
  Mesh::NodalCoordinates::~NodalCoordinates(){destroy();};
  bool Mesh::NodalCoordinates::good() const { return(ncdata ? true : false);};
  Mesh::IndexType Mesh::NodalCoordinates::size() const { return(nnodes);};
  Mesh::IndexType Mesh::NodalCoordinates::Size() const { return(nnodes);};
  void Mesh::NodalCoordinates::destroy(){
    if(ncdata && mydata){
      delete [] ncdata;
      nnodes = 0;
      ncdata = NULL;
      mydata = false;
    }
  };
  void Mesh::NodalCoordinates::init(){destroy();};
  void Mesh::NodalCoordinates::init(IndexType n)
  {
    destroy();
    if(n > 0){
      ncdata = new double [3 * n];
      nnodes = n;
      mydata = true;
    }
  };
  void Mesh::NodalCoordinates::init(IndexType n, double *data)
  {
    destroy();
    if(n > 0){
      ncdata = data;
      nnodes = n;
      mydata = false;
    }
  };

  void Mesh::NodalCoordinates::init_node(IndexType n,const GeoPrim::CPoint &point)
  {
    this->x(n) = point.x();
    this->y(n) = point.y();
    this->z(n) = point.z();
  };

  const GeoPrim::CPoint NodalCoordinates::closest_point(const GeoPrim::CPoint &p) const
  {
    double dist = 1000000;
    GeoPrim::CPoint retval;
    for(unsigned int i = 0;i < nnodes;i++){
      GeoPrim::CPoint testp(&ncdata[3*i]);
      double testd = (p-testp).norm();
      if(testd < dist){
	dist = testd;
	retval = testp;
      }
    }
    return(retval);
  };

  Mesh::IndexType NodalCoordinates::closest_node(const GeoPrim::CPoint &p,double *dist_ptr) const
  {
    double dist = 1000000;
    Mesh::IndexType reti = 0;
    for(unsigned int i = 0;i < nnodes;i++){
      GeoPrim::CPoint testp(&ncdata[3*i]);
      double testd = (p-testp).norm();
      if(testd < dist){
	dist = testd;
	reti = i + 1;
      }
    }
    if(dist_ptr)
      *dist_ptr = dist;
    return(reti);
  };

  void NodalCoordinates::init_copy(IndexType n, double *data)
  {
    destroy();
    if(n > 0 && data){
      ncdata = new double [3*n];
      nnodes = n;
      mydata = true;
      std::memcpy(ncdata,data,n*sizeof(double)*3);
    }
  };
 
  int CollideCellsWithBox(Mesh::NodalCoordinates &nc,
			  Mesh::Connectivity     &conn,
			  GeoPrim::CBox          &box,
			  std::vector<Mesh::IndexType> &candidates,
			  std::vector<Mesh::IndexType> &cells)
 {
    std::vector<IndexType>::iterator ci = candidates.begin();
    while(ci != candidates.end()){
      unsigned int cell_id = *ci++;
      unsigned int esize = conn.Esize(cell_id);
      std::vector<GeoPrim::CPoint> element_points(esize);
      // Make a vector of the element points
      for(IndexType node = 1;node <= esize;node++){
	GeoPrim::CPoint ep(nc[conn.Node(cell_id,node)]);
	element_points[node-1] = ep;
      }
      // Make the bounding box for the element
      GeoPrim::CBox ebox(element_points);
      // collide the ebox with the box
      GeoPrim::CBox cbox(box.intersect(ebox));
      if(!cbox.empty())
	cells.push_back(cell_id);
    }
    
  };

  int CollideMeshWithBox(Mesh::NodalCoordinates &nc,
			 Mesh::Connectivity     &conn,
			 GeoPrim::CBox          &box,
			 std::vector<Mesh::IndexType> &cells)
  {
    unsigned int nelem = conn.size();
    IndexType n = 1;
    while(n <= nelem){
      unsigned int esize = conn.Esize(n);
      std::vector<GeoPrim::CPoint> element_points(esize);
      // Make a vector of the element points
      for(IndexType node = 1;node <= esize;node++){
	GeoPrim::CPoint ep(nc[conn.Node(n,node)]);
	element_points[node-1] = ep;
      }
      // Make the bounding box for the element
      GeoPrim::CBox ebox(element_points);
      // collide the ebox with the box
      GeoPrim::CBox cbox = box.intersect(ebox);
      if(!cbox.empty())
	cells.push_back(n);
      n++;
    }
    
  };

  int GetMeshCentroids(Mesh::NodalCoordinates &nc,
		       Mesh::Connectivity     &conn,
		       std::vector<double>    &centroids)
  {
    unsigned int nelem = conn.size();
    centroids.resize(3*nelem,0.0);
    for(int i = 0;i < nelem;i++){
      GeoPrim::C3Point p(&(centroids[3*i]));
      Mesh::GenericElement e;
      e.Centroid(conn[i],nc,p);
      //      std::cout << "a" << i+1 << ": " << p << std::endl;
    }
  };
  
  void GetCoordinateBounds(NodalCoordinates &nc,std::vector<double> &bnds)
  {
    bnds.resize(6);
    bnds[0] = bnds[2] = bnds[4] = std::numeric_limits<double>::max();
    bnds[1] = bnds[3] = bnds[5] = std::numeric_limits<double>::min();
    Mesh::IndexType number_of_nodes = nc.Size();
    for(Mesh::IndexType nn = 1;nn <= number_of_nodes;nn++){
      if(nc.x(nn) < bnds[0])
	bnds[0] = nc.x(nn);
      if(nc.x(nn) > bnds[1])
	bnds[1] = nc.x(nn);
      if(nc.y(nn) < bnds[2])
	bnds[2] = nc.y(nn);
      if(nc.y(nn) > bnds[3])
	bnds[3] = nc.y(nn);
      if(nc.z(nn) < bnds[4])
	bnds[4] = nc.z(nn);
      if(nc.z(nn) > bnds[5])
	bnds[5] = nc.z(nn);
    }
  };


  bool TestFace::operator<(const TestFace &f)
  {
    return(IRAD::Util::LessThan(*this,f));
  };

  Mesh::Connectivity::Connectivity() : std::vector<std::vector<IndexType> >()
  {
  	this ->_nelem = 0;
    this->resize(0);
    std::vector<std::vector<Mesh::IndexType> >(*this).swap(*this);
  };

  Mesh::Connectivity::Connectivity(unsigned int N) : std::vector<std::vector<IndexType> >(N)
  {
    std::vector<std::vector<Mesh::IndexType> >(*this).swap(*this);
  };
  void Mesh::Connectivity::Resize(unsigned int N){this->resize(N);};
  Mesh::Connectivity::~Connectivity() { destroy(); };
  void Mesh::Connectivity::AddElement(const std::vector<IndexType> &elem)
  {
    _nelem++;
    this->push_back(elem);
  };
  void Mesh::Connectivity::AddElement()
  {
    _nelem++;
    std::vector<IndexType> a;
    this->push_back(a);
  };
  void Mesh::Connectivity::Sync(){ _nelem = this->size();};
  void Mesh::Connectivity::ShrinkWrap()
  {
    _nelem = this->size();
    Connectivity::iterator ti = this->begin();
    while(ti != this->end()){
      std::vector<Mesh::IndexType>(*ti).swap(*ti);
      ti++;
    }
    std::vector<std::vector<Mesh::IndexType> >(*this).swap(*this);
  };
  void Mesh::Connectivity::Truncate()
  {
    Connectivity::iterator ti = this->begin();
    while(ti != this->end()){
      if(ti->empty()){
	ti = this->erase(ti);
	_nelem--;
      }
      else
	ti++;
    }
  };
  void Mesh::Connectivity::destroy(){
    DestroySizes();
    Connectivity::iterator ti = this->begin();
    while(ti != this->end()){
      ti->resize(0);
      std::vector<IndexType> tmp(*ti);
      ti->swap(tmp);
      ti++;
    }
    this->resize(0);
    std::vector<std::vector<IndexType> >(*this).swap(*this);
  };
  void Mesh::Connectivity::SyncSizes(){
    _sizes.resize(_nelem);
    for(IndexType i = 0;i < _nelem;i++){
      _sizes[i] = (*this)[i].size();
      //	std::vector<IndexType> tmp((*this)[i]);
      //	((*this)[i].swap(tmp));
    }
  };
  void Mesh::Connectivity::DestroySizes(){ _sizes.resize(0);};

  // warning, if called with unsync'd Connectivity table, result
  // is silently garbage.
  void Mesh::Connectivity::InverseDegenerate(Connectivity &rc,IndexType nnodes) const
  {
    if(nnodes <= 0)
      nnodes = MaxNodeId<Connectivity,std::vector<IndexType> >(*this);
    rc.Resize(nnodes);
    rc.Sync();
    for(IndexType i = 0;i < _nelem;i++){
      if((*this)[i].size() > 1){
	std::vector<IndexType>::const_iterator ii = (*this)[i].begin();
	while(ii != (*this)[i].end()){
	  rc[*ii-1].push_back(i+1);
	  ii++;
	}
      }
    }
  };

  // warning, if called with unsync'd Connectivity table, result
  // is silently garbage.
  void Connectivity::Inverse(Connectivity &rc,IndexType nnodes) const
  {
    if(nnodes <= 0)
      nnodes = MaxNodeId<Connectivity,std::vector<IndexType> >(*this);
    //      Connectivity::const_iterator ci = this->begin();
    //      while(ci != this->end()){
    //	std::vector<IndexType>::const_iterator ii = ci->begin();
    //	while(ii != ci->end()){
    //	  nnodes = (nnodes > *ii ? nnodes : *ii);
    //	  ii++;
    //	}
    //	ci++;
    //      }
    //    }
    rc.Resize(nnodes);
    rc.Sync();
    for(IndexType i = 0;i < _nelem;i++){
      std::vector<IndexType>::const_iterator ii = (*this)[i].begin();
      while(ii != (*this)[i].end()){
	rc[*ii-1].push_back(i+1);
	ii++;
      }
    }
  };

  void Connectivity::BuildFaceConnectivity(Connectivity &fcon,Connectivity &ef,
					   std::vector<SymbolicFace> &sf,Connectivity &dc) const
  {
    //    this->Sync();
    //    this->SyncSizes();
    ef.Resize(this->Nelem());
    Mesh::IndexType number_of_elements = this->Nelem();
    Mesh::IndexType nface_estimate = static_cast<Mesh::IndexType>(2.2 * number_of_elements);
    fcon.Resize(0);
    fcon.reserve(nface_estimate);
    sf.reserve(nface_estimate);
    Mesh::IndexType number_of_faces = 0;
    std::vector<Mesh::Connectivity> all_face_conn;
    all_face_conn.resize(this->Nelem());
    //    std::vector<Mesh::IndexType> element_nfaces;
    //    element_nfaces.resize(this->Nelem());
    // This loop sizes the C[F] array (i.e. for each cell, which faces)
    // They're init'd to 0 so that we can tell which faces have been processed
    for(Mesh::IndexType element_being_processed = 1;
	element_being_processed <= number_of_elements;
	element_being_processed++)
      {
	Mesh::IndexType index = element_being_processed - 1;
	Mesh::IndexType size_of_element = this->Esize(element_being_processed);
	Mesh::GenericElement ge(size_of_element);
	//	element_nfaces[index] = nfaces;
	ef[index].resize(ge.nfaces(),0);
	//	std::vector<Mesh::IndexType>(ef[index]).swap(ef[index]);
	ge.get_face_connectivities(all_face_conn[index],(*this)[index]);
      }
    ef.Sync();
    ef.SyncSizes();
    // This loop populates the F[N] (i.e. for each face, which nodes), and
    // the C[F] arrays.
    for(Mesh::IndexType element_being_processed = 1;
	element_being_processed <= number_of_elements;
	element_being_processed++)
      {
	Mesh::IndexType index = element_being_processed - 1;
	Mesh::IndexType size_of_element = this->Esize(element_being_processed);
	Mesh::GenericElement ge(size_of_element);
	Mesh::IndexType nfaces = ef.Esize(element_being_processed);
	//	Mesh::Connectivity efc;
	// This call populates the F[N] for the faces local
	// to the element_being_processed.
	//	ge.get_face_connectivities(efc,(*this)[index]);
	// This loop goes through each face of the element_being_processed
	// and updates the F[N] and C[F] arrays.
	for(Mesh::IndexType face_being_processed = 1;
	    face_being_processed <= nfaces;
	    face_being_processed++){
	  Mesh::IndexType findex = face_being_processed - 1;
	  if(!ef[index][findex]){ // 0 if face hasn't yet been processed
	    fcon.push_back(all_face_conn[index][findex]); // add the new face to the F[N] array
	    ef[index][findex] = number_of_faces + 1; // add the face id to the C[F] array
	    SubEntityId seid1;
	    SubEntityId seid2;
	    sf.push_back(std::make_pair(seid1,seid2));
	    sf[number_of_faces].first.first = element_being_processed;
	    sf[number_of_faces].first.second = face_being_processed;
	    // Now look at each cell containing each node of the current face and determine
	    // which one (if any) have the same face.   This will be the face neighbor of the
	    // element_being_processed for the face identified by (number_of_faces + 1).
	    //	    std::list<Mesh::IndexType> element_nbrlist;
	    bool found = false;
	    std::vector<Mesh::IndexType>::iterator fni = fcon[number_of_faces].begin();
	    while(fni != fcon[number_of_faces].end() && !found){ // for every node in this face
	      Mesh::IndexType face_node_index = *fni++ - 1;
	      std::vector<Mesh::IndexType>::iterator enbri = dc[face_node_index].begin();
	      while(enbri != dc[face_node_index].end() && !found){ // look thru all the node's cells
		Mesh::IndexType enbr_index = *enbri++ - 1;
		if(enbr_index+1 > element_being_processed){
		  Mesh::IndexType nbr_nfaces = ef.Esize(enbr_index+1);
		  for(Mesh::IndexType nbrface = 1;nbrface <= nbr_nfaces && !found;nbrface++){
		    Mesh::IndexType nbr_face_index = nbrface -1;
		    if(!ef[enbr_index][nbr_face_index]){
		      if(IRAD::Util::HaveOppositeOrientation(all_face_conn[index][findex],
						 all_face_conn[enbr_index][nbr_face_index])){
			found = true;
			ef[enbr_index][nbr_face_index] = number_of_faces + 1;
			sf[number_of_faces].second.first = enbr_index+1;
			sf[number_of_faces].second.second = nbr_face_index+1;
		      } // if faces match
		    }
		  }
		}
	      }
	    }
	    number_of_faces++;
	  }
	}
      }
    fcon.Sync();
  };

  // warning, if called with unsync'd Connectivity table, result
  // is silently garbage.
  void Connectivity::graph_mode(IndexType offset)
  {
    Connectivity::iterator ci = this->begin();
    IndexType n = 1;
    while(ci != this->end()){
      ci->erase(std::lower_bound(ci->begin(),ci->end(),offset+n++));
      ci++;
    }
  };

  void Connectivity::matrix_mode(IndexType offset)
  {
    Connectivity::iterator ci = this->begin();
    IndexType n = 1;
    while(ci != this->end()){
      ci->insert(std::lower_bound(ci->begin(),ci->end(),offset+n),1,offset+n);
      n++;
      ci++;
    }
  };

  // input is dual connectivity (i.e. for every node, which elements)
  void Connectivity::GetNeighborhood(NeighborHood &rl,
				     const Connectivity &dc,
				     IndexType size)
  {
    rl.resize(_nelem);
    //    rl._nelem = _nelem;
    for(IndexType i = 0; i < _nelem;i++){
      std::vector<IndexType>::const_iterator ni = (*this)[i].begin();
      while(ni != (*this)[i].end()){
	IndexType index = *ni - 1;
	std::vector<IndexType>::const_iterator dci = dc[index].begin();
	while(dci != dc[index].end())
	  rl[i].insert(*dci++);
	ni++;
      }
      rl[i].erase(i+1);
    }
  };

  // input is dual connectivity (i.e. for every node, which elements)
  void Connectivity::GetNeighborhood(Connectivity &rl,
				     Connectivity &dc,
				     bool exclude_self,
				     bool sortit)
  {
    rl.Resize(_nelem);
    rl._nelem = _nelem;
    std::vector<bool> added(_nelem,false);
    for(IndexType i = 0; i < _nelem;i++){
      IndexType current_element = i + 1;
      std::vector<IndexType>::iterator ni = (*this)[i].begin();
      std::list<IndexType> nbrlist;
      while(ni != (*this)[i].end()){
	IndexType index = *ni - 1;
	std::vector<IndexType>::iterator dci = dc[index].begin();
	while(dci != dc[index].end()){
	  IndexType ai = *dci - 1;
	  if(!added[ai]){
	    nbrlist.push_back(*dci);
	    added[ai] = true;
	  }
	  dci++;
	}
	ni++;
      }
      if(sortit)
	nbrlist.sort();
      //      nbrlist.unique();
      if(exclude_self){
	nbrlist.remove(current_element);
        added[i] = false;
      }
      rl[i].resize(nbrlist.size());
      IndexType j = 0;
      std::list<IndexType>::iterator si = nbrlist.begin();
      while(si != nbrlist.end()){
	rl[i][j++] = *si;
	added[*si-1] = false;
	si++;
      }
    }
  };

  // input is dual connectivity (i.e. for every node, which elements)
  void Connectivity::GetAdjacent(Connectivity &rl,
				 Connectivity &dc,
				 IndexType n,
				 bool sortit)
  {
    rl.Resize(_nelem);
    rl._nelem = _nelem;
    IndexType nadj = n;
    if(n == 0){
      std::cout << "GETADJACENT: Determining MaxNodeId" << std::endl;
      nadj = MaxNodeId<Connectivity,std::vector<IndexType> >(dc);
    }
    std::vector<bool> added(nadj,false);
    for(IndexType i = 0; i < _nelem;i++){
      //      IndexType current_element = i + 1;
      std::vector<IndexType>::iterator ni = (*this)[i].begin();
      std::list<IndexType> nbrlist;
      while(ni != (*this)[i].end()){
	IndexType index = *ni - 1;
	std::vector<IndexType>::iterator dci = dc[index].begin();
	while(dci != dc[index].end()){
	  IndexType ai = *dci - 1;
	  if(!added[ai]){
	    nbrlist.push_back(*dci);
	    added[ai] = true;
	  }
	  dci++;
	}
	ni++;
      }
      if(sortit)
	nbrlist.sort();
      //      nbrlist.unique();
      rl[i].resize(nbrlist.size());
      IndexType j = 0;
      std::list<IndexType>::iterator si = nbrlist.begin();
      while(si != nbrlist.end()){
	rl[i][j++] = *si;
	added[*si-1] = false;
	si++;
      }
    }
  };


  // Doesn't have to be nodes - this is a graph type operation
  // Takes a list of nodes in nodes and determines a list of 
  // unique elements they touch as indicated by the dual 
  // connectivity - can be removed from Connectivity class
  void Connectivity::ElementsOn(std::vector<Mesh::IndexType> &nodes,
				Connectivity &dc,
				std::vector<Mesh::IndexType> &subset)
  {
    std::list<Mesh::IndexType> alist;
    std::vector<Mesh::IndexType>::iterator ni = nodes.begin();
    while(ni != nodes.end()){
      Mesh::IndexType node_index = *ni++ - 1;
      std::vector<Mesh::IndexType>::iterator dci = dc[node_index].begin();
      while(dci != dc[node_index].end())
	alist.push_back(*dci++);
    }
    alist.sort();
    alist.unique();
    subset.resize(alist.size());
    std::vector<Mesh::IndexType>::iterator ssi = subset.begin();
    std::list<Mesh::IndexType>::iterator li = alist.begin();
    while(li != alist.end())
      *ssi++ = *li++;
  };

				
  // Does breadth first renumbering and produces the remap: 
  // remap[old_id] = new_id
  void Connectivity::BreadthFirstRenumber(std::vector<Mesh::IndexType> &remap) 
  {
    remap.resize(_nelem,0);
    IndexType renumber = 1;
    for(IndexType i = 0;i < _nelem && renumber <= _nelem;i++){
      //    while(renumber <= _nelem){
      if(remap[i] == 0){
	remap[i] = renumber++; 
	std::list<IndexType> processing_queue;
	std::vector<IndexType>::iterator ni = (*this)[i].begin();
	while(ni != (*this)[i].end()){
	  IndexType index = *ni++ - 1;
	  if(remap[index] == 0){
	    processing_queue.push_back(index);
	    remap[index] = renumber++;
	  }
	}
	std::list<IndexType>::iterator pqi = processing_queue.begin();
	while(pqi != processing_queue.end()){
	  IndexType index = *pqi;
	  std::vector<IndexType>::iterator ini = (*this)[index].begin();
	  while(ini != (*this)[index].end()){
	    IndexType iindex = *ini++ - 1;
	    if(remap[iindex] == 0){
	      processing_queue.push_back(iindex);
	      remap[iindex] = renumber++;
	    }
	  }
	  pqi++;
	}
      }
    }
    assert((renumber == (_nelem+1)));
  };


  GeoPrim::C3Point
  GenericCell_2::Centroid(std::vector<IndexType> &ec,NodalCoordinates &nc) const
  {
    GeoPrim::C3Point centroid(0,0,0);
    IndexType i = 0;
    IndexType esize = ec.size();
    while(i < esize)
      centroid += GeoPrim::C3Point(nc[ec[i++]]);
    double scal = 1.0/(static_cast<double>(esize));
    return(centroid *= scal);
  };

  // Populate SF with the shape function of the element at the nat.c.
  void
  GenericCell_2::shape_func(const GeoPrim::CVector &natc,double SF[]) const
  {
    switch (_size) {
    case 3: {
      SF[0] = 1. - natc[0] - natc[1];
      SF[1] = natc[0];
      SF[2] = natc[1];
      return;
    }
    case 4: {
      const double xi = natc[0];
      const double xi_minus = 1. - xi;
      const double eta = natc[1];
      const double eta_minus = 1. - eta;

      SF[0] = xi_minus * eta_minus;
      SF[1] = xi * eta_minus;
      SF[2] = xi * eta;
      SF[3] = xi_minus * eta;
      return;
    }
    default:
      std::cerr << "GenericCell_2::shape_func:Error: unkown element type."
		<< std::endl;
      exit(1);
    }
  };


  // This function evaluates the shape function of an element. It
  // computes the barycentric coordinates from the natural coordinates.
  void
  GenericCell_2::dshape_func( const GeoPrim::CVector &natc, double Np[][2]) const {
    switch ( _size) {
    case 3: {
      Np[0][0] = -1;    Np[0][1] = -1;
      Np[1][0] =  1;    Np[1][1] =  0;
      Np[2][0] =  0;    Np[2][1] =  1;
      return;
    }
    case 4: {
      const double xi = natc[0];
      const double xi_minus = 1. - xi;
      const double eta = natc[1];
      const double eta_minus = 1. - eta;

      Np[0][0] = -eta_minus;   Np[0][1] = -xi_minus;
      Np[1][0] = eta_minus;    Np[1][1] = -xi;
      Np[2][0] = eta;          Np[2][1] = xi;
      Np[3][0] = -eta;         Np[3][1] = xi_minus;

      return;
    }
    default:
      std::cerr << "GenericCell_2::dshape_func:error: Unknown element type."
 		<< std::endl;
      exit(1);
    }
  };
  // This function evaluates the Jacobian at a given point.
  // J and nc are both size 2 (nc = natural coordinates)
  void GenericCell_2::jacobian( const GeoPrim::CPoint p[],
 				const GeoPrim::CVector &nc,
 				GeoPrim::CVector J[]) const {
    switch (_size) {
    case 3:
      J[0] = p[1]-p[0]; J[1] = p[2]-p[0];
      return;
    case 4: {
      const double xi        = nc[0];
      const double xi_minus  = 1. - xi;
      const double eta       = nc[1];
      const double eta_minus = 1. - eta;
      J[0] = ((p[1]-p[0]) *= eta_minus) += (( p[2]-p[3]) *= eta );
      J[1] = ((p[3]-p[0]) *= xi_minus) += (( p[2]-p[1]) *= xi);
      return;
    }
    case 6: {
      const double xi   = nc[0];
      const double eta  = nc[1];
      const double zeta = 1.-xi-eta;

      J[0] = ((p[1]-p[0]) *= 4.*xi -1.)
	+= ((p[3]-p[0]) *= 4.*zeta-4*xi)
	+= ((p[4]-p[5]) *= 4.*eta);
      J[1] = ((p[2]-p[0]) *= 4.*eta -1.)
	+= ((p[5]-p[0]) *= 4*zeta-4*eta)
	+= ((p[4]-p[3]) *= 4.*xi);
      return;
    }
    case 8: {
      const double xi        = nc[0];
      const double xi_minus  = 1. - xi;
      const double eta       = nc[1];
      const double eta_minus = 1. - eta;

      J[0] = ((p[1]-p[0]) *= eta_minus) += ((p[2]-p[3]) *= eta)
	-= ((((p[0]-p[4])+=(p[1]-p[4])) *= 2.*eta_minus*(xi_minus-xi))
 	    += (((p[2]-p[6])+=(p[3]-p[6])) *= 2.*eta*(xi_minus-xi))
 	    += (((p[1]-p[5])+=(p[2]-p[5])-=((p[0]-p[7])+=(p[3]-p[7])))
 		*= 2.*eta*eta_minus));

      ((J[1] = p[3]-p[0]) *= xi_minus) += ((p[2]-p[1]) *= xi)
 	-= ((((p[1]-p[5])+=(p[2]-p[5])) *= 2.*xi*(eta_minus-eta))
 	    += (((p[0]-p[7])+=(p[3]-p[7])) *= 2.*xi_minus*(eta_minus-eta))
	    += (((p[2]-p[6])+=(p[1]-p[4])-=((p[3]-p[6])+=(p[0]-p[4])))
 		*= 2.*xi*xi_minus));
      return;
    }
    default:
      abort(); // Should never reach here
    }
  };

  void GenericCell_2::GetNormalSet(GeoPrim::CVector ns[],
				   IndexType elnum,
				   const Connectivity &ec,
				   const NodalCoordinates &nc)
     {
       GeoPrim::CPoint P[3];
       for(IndexType i = 0;i < _size;i++){
	 IndexType i_minus = (i==0 ? _size -1 : i - 1);
	 P[i_minus] = nc[ec.Node(elnum,i_minus+1)];
	 P[i] = nc[ec.Node(elnum,i+1)];
	 IndexType i_plus = (i == (_size-1) ? 0 : i+1);
	 P[i_plus] = nc[ec.Node(elnum,i_plus+1)];
	 ns[i] = GeoPrim::CVector(P[i] - P[i_minus])%GeoPrim::CVector(P[i_plus] - P[i]);
       }
     };

  void GenericCell_2::GetPointSet(GeoPrim::CPoint P[],
				  IndexType elnum,
				  const Connectivity &ec,
				  const NodalCoordinates &nc)
     {
       for(IndexType i = 0;i < _size;i++)
	 P[i].init(nc[ec.Node(elnum,i+1)]);
     };

  // Evaluates the shape function and jacobian at a given point
  void
  GenericCell_2::shapef_jacobian_at(const GeoPrim::CPoint &p,
 				    GeoPrim::CVector &natc,
 				    IndexType elnum,
 				    const Connectivity &ec,
 				    const NodalCoordinates &nc,
 				    GeoPrim::CVector &fvec,
 				    GeoPrim::CVector jac[]) const // fjac is 3x2
  {
    GeoPrim::CPoint P[(const IndexType)_size];
    double SF[(const IndexType)_size];
    this->shape_func(natc,SF);
    fvec.init(-1.0*p);
    for(IndexType i = 0;i < _size;i++){
      P[i].init(nc[ec.Node(elnum,i+1)]);
      fvec += SF[i]*P[i];
    }
    GeoPrim::CVector v1(P[0],P[1]);
    GeoPrim::CVector v2(P[0],P[2]);
    this->jacobian(P,natc,jac);
    jac[2] = v1%v2;
    jac[2].normalize();
    Transpose(jac);
  }

  // Interpolate a field value from the nodes of the element to the point
  // having natural coordinates nc.
  void GenericCell_2::interpolate(const GeoPrim::CVector f[],
				  const GeoPrim::CVector &nc,
				  GeoPrim::CVector &v) const
  {
    const double xi = nc.x();
    const double eta = nc.y();
    v = f[0];
    switch(_size) { 
    case 3: {
      v += ((f[1]-f[0])*=xi) += ((f[2]-f[0])*=eta);
      return;
    }
    case 4: {
      v += ((f[1]-f[0]) *= (xi * (1.-eta)))
	+= ((f[3]-f[0]) *= eta)
	+= ((f[2]-f[3]) *= xi*eta);
      return;
    }
   case 6: {
      const double zeta=1.-xi-eta;

      v += ((f[1]-f[0]) *= xi * (2.*xi-1.))
	+= ((f[3]-f[0]) *= 4.* xi * zeta)
	+= ((f[2]-f[0]) *= eta * (2.*eta-1.))
	+= ((f[5]-f[0]) *= 4. * eta * zeta)
	+= ((f[4]-f[0]) *= 4.*xi*eta);
      return;
    }
    case 8: {
      const double xi_minus = 1. - xi;
      const double eta_minus = 1. - eta;
      v += ((f[1]-f[0]) *= eta * (2.*eta-1.))
	+= ((f[3]-f[0]) *= eta)
	+= ((f[2]-f[3]) *= xi * eta)
	-= ((((f[0]-f[4])+=(f[1]-f[4])) *= 2.*xi*xi_minus*eta_minus)
	    += (((f[2]-f[6])+=(f[3]-f[6])) *= 2.*xi*xi_minus*eta)
	    += (((f[1]-f[5])+=(f[2]-f[5])) *= 2.*xi*eta*eta_minus)
	    += (((f[0]-f[7])+=(f[3]-f[7])) *= 2.*xi_minus*eta*eta_minus));
      return;
    }
    default:
      assert(false);
    }
    return;
  }

  void GenericCell_2::ReOrient(std::vector<IndexType> &ec)
  {
    switch(_size) { 
    case 3: {
      Mesh::IndexType temp = ec[0];
      ec[0] = ec[1];
      ec[1] = temp;
      return;
    }
    case 4: {
      Mesh::IndexType temp = ec[0];
      ec[0] = ec[1];
      ec[1] = temp;
      temp = ec[2];
      ec[2] = ec[3];
      ec[3] = temp;
      return;
    }
    default:
      assert(false);
    }
    return;
    
  }

  // Evaluates the shape function and jacobian at a given point
  void
  GenericElement::shapef_jacobian_at(const GeoPrim::CPoint &p,
				     GeoPrim::CVector &natc,
				     IndexType elnum,
				     const Connectivity &ec,
				     const NodalCoordinates &nc,
				     GeoPrim::CVector &fvec,
				     GeoPrim::CVector fjac[]) const
  {
    GeoPrim::CPoint P[(const IndexType)_size];
    double SF[(const IndexType)_size];
    this->shape_func(natc,SF);
    fvec.init(-1.0*p);
    for(IndexType i = 0;i < _size;i++){
      //      std::cout << "debugging: " << GeoPrim::CPoint(nc[ec.Node(elnum,i+1)]) << std::endl;
      P[i].init(nc[ec.Node(elnum,i+1)]);
      fvec += SF[i]*P[i];
    }
    this->jacobian(P,natc,fjac);
    Transpose(fjac);
  }

  // Populate SF with the shape function of the element at the nat.c.
  void
  GenericElement::shape_func(const GeoPrim::CVector &nc,double SF[]) const
  {
    switch (_size) {
      // Tets are 4 nodes or 10
    case 4: {
      SF[0] = 1. - nc.x() - nc.y() - nc.z();
      SF[1] = nc.x();
      SF[2] = nc.y();
      SF[3] = nc.z();
      break;
    }
    case 10: {
      const double xi = nc.x();
      const double eta = nc.y();
      const double zeta = nc.z();
      const double alpha = (1. - xi - eta - zeta);
      SF[0] = alpha*(1. - 2.*(xi + eta + zeta));
      SF[1] = xi *(2.* xi - 1.);
      SF[2] = eta *(2. * eta - 1.);
      SF[3] = zeta *(2. * zeta - 1.);
      SF[4] = 4.* xi * alpha;
      SF[5] = 4.* eta * alpha;
      SF[6] = 4.* zeta * alpha;
      SF[7] = 4. * xi * eta;
      SF[8] = 4. * eta * zeta;
      SF[9] = 4. * xi * zeta;
      break;
    }
      // Hex's are 8 nodes or 20
    case 8: {
      const double xi = nc.x();
      const double xi_minus = 1. - xi;
      const double eta = nc.y();
      const double eta_minus = 1. - eta;
      const double zeta = nc.z();
      const double zeta_minus = 1. - zeta;
      SF[0] = xi_minus * eta_minus * zeta_minus;
      SF[1] = xi * eta_minus * zeta_minus;
      SF[2] = xi * eta * zeta_minus;
      SF[3] = xi_minus * eta * zeta_minus;
      SF[4] = xi_minus * eta_minus * zeta;
      SF[5] = xi * eta_minus * zeta;
      SF[6] = xi * eta * zeta;
      SF[7] = xi_minus * eta * zeta;
      break;
    }
    default:
      std::cerr << "GenericElement::shape_func:Error: unkown element type."
		<< std::endl;
      exit(1);
    }
  }

  // Populate dSF with the derivative of the shape function for the element
  // at nat.c.
  void
  GenericElement::dshape_func(const GeoPrim::CVector &nc,double dSF[][3]) const
  {
    switch (_size) {
    case 4: {
      dSF[0][0] = -1;  dSF[0][1] = -1;  dSF[0][2] = -1;
      dSF[1][0] =  1;  dSF[1][1] =  0;  dSF[1][2] =  0;
      dSF[2][0] =  0;  dSF[2][1] =  1;  dSF[2][2] =  0;
      dSF[3][0] =  0;  dSF[3][1] =  0;  dSF[3][2] =  1;
      break;
    }
    case 10:{
      const double xi = nc.x();
      const double eta = nc.y();
      const double zeta = nc.z();
      const double alpha = (1. - xi - eta - zeta);
      dSF[0][0] = (4.*(xi+eta+zeta)-3.);      dSF[0][1] = dSF[0][0];                 dSF[0][2] = dSF[0][0];
      dSF[1][0] = 4.*xi - 1.;                 dSF[1][1] = 0;                         dSF[1][2] = 0;
      dSF[2][0] = 0;                          dSF[2][1] = 4.*eta - 1.;               dSF[2][2] = 0;
      dSF[3][0] = 0;                          dSF[3][1] = 0;                         dSF[3][2] = 4.*zeta - 1.;
      dSF[4][0] = 4.*(alpha - xi);            dSF[4][1] = -4.*xi;                    dSF[4][2] = -4.*xi;
      dSF[5][0] = -4.*eta;                    dSF[5][1] = 4.*(alpha - eta);          dSF[5][2] = -4.*eta;
      dSF[6][0] = -4.*zeta;                   dSF[6][1] = -4.*zeta;                  dSF[6][2] = 4.*(alpha - zeta);
      dSF[7][0] = 4.*eta;                     dSF[7][1] = 4.*xi;                     dSF[7][2] = 0;
      dSF[8][0] = 0;                          dSF[8][1] = 4.*zeta;                   dSF[8][2] = 4.*eta;
      dSF[9][0] = 4.*zeta;                    dSF[9][1] = 0;                         dSF[9][2] = 4.*xi;
      break;
    }
    case 8: {
      const double xi = nc.x();
      const double xi_minus = 1. - xi;
      const double eta = nc.y();
      const double eta_minus = 1. - eta;
      const double zeta = nc.z();
      const double zeta_minus = 1. - zeta;
      dSF[0][0] = -1.*eta_minus*zeta_minus;  dSF[0][1] = -1.*xi_minus*zeta_minus;  dSF[0][2] = -1.*xi_minus*eta_minus;
      dSF[1][0] = eta_minus*zeta_minus;      dSF[1][1] = -1.*xi*zeta_minus;        dSF[1][2] = -1.*xi*eta_minus;
      dSF[2][0] = eta*zeta_minus;            dSF[2][1] = xi*zeta_minus;            dSF[2][2] = -1.*xi*eta;
      dSF[3][0] = -1.*eta*zeta_minus;        dSF[3][1] = xi_minus*zeta_minus;      dSF[3][2] = -1.*xi_minus*eta;
      dSF[4][0] = -1.*eta_minus*zeta;        dSF[4][1] = -1.*xi_minus*zeta;        dSF[4][2] = xi_minus*eta_minus;
      dSF[5][0] = eta_minus*zeta;            dSF[5][1] = -1.*xi*zeta;              dSF[5][2] = xi*eta_minus;
      dSF[6][0] = eta*zeta;                  dSF[6][1] = xi*zeta;                  dSF[6][2] = xi*eta;
      dSF[7][0] = -1.*eta*zeta;              dSF[7][1] = xi_minus * zeta;          dSF[7][2] = xi_minus*eta;
      break;
    }
    default:
      std::cerr << "GenericElement::dshape_func:error: Unknown element type."
		<< std::endl;
      exit(1);
    }
  }

  // Populate J with the jacobian for the point having natural
  // coordinates nc.  p are the nodal coordinates of the vertices
  void
  GenericElement::jacobian(const GeoPrim::CPoint p[],
			   const GeoPrim::CVector &nc,
			   GeoPrim::CVector J[]) const
  {
    switch(_size){
    case 4: {
      J[0] = p[1] - p[0];
      J[1] = p[2] - p[0];
      J[2] = p[3] - p[0];
      break;
    }
    case 10: {
      const double xi = nc.x();
      const double eta = nc.y();
      const double zeta = nc.z();
      const double alpha = (1. - xi - eta - zeta);
      GeoPrim::CPoint P(p[0]*(4.*(xi+eta+zeta)-3.));
      J[0] = ((p[9]-p[6])*=4.*zeta)+=((p[7]-p[5])*=4.*eta)+=
	(p[4]*(4.*(alpha-xi))+p[1]*(4.*xi-1.)+P);
      J[1] = ((p[8]-p[6])*=4.*zeta)+=((p[7]-p[4])*=4.*xi)+=
	(p[5]*(4.*(alpha-eta))+p[2]*(4.*eta-1.)+P);
      J[2] = ((p[9]-p[4])*=4.*xi)+=((p[8]-p[5])*=4.*eta)+=
	(p[6]*(4.*(alpha-zeta))+p[3]*(4.*zeta-1.)+P);
      break;
    }
    case 8: {
      const double xi = nc.x();
      const double xi_minus = 1. - xi;
      const double eta = nc.y();
      const double eta_minus = 1. - eta;
      const double zeta = nc.z();
      const double zeta_minus = 1. - zeta;
      J[0] = ((p[6]-p[7])*=eta*zeta)+=((p[5]-p[4])*=eta_minus*zeta)+=
	((p[2]-p[3])*=eta*zeta_minus)+=((p[1]-p[0])*=eta_minus*zeta_minus);
      J[1] = ((p[7]-p[4])*=xi_minus*zeta)+=((p[6]-p[5])*=xi*zeta)+=
	((p[3]-p[0])*=xi_minus*zeta_minus)+=((p[2]-p[1])*=xi*zeta_minus);
      J[2] = ((p[7]-p[3])*=xi_minus*eta)+=((p[6]-p[2])*=xi*eta)+=
	((p[5]-p[1])*=xi*eta_minus)+=((p[4]-p[0])*=xi_minus*eta_minus);
      break;
    }
    default:
      std::cerr << "GenericElement::jacobian:Error: Cannot handle this"
		<< " element size (yet)." << std::endl;
      exit(1);
    }
  }

  // Interpolate a field value from the nodes of the element to the point
  // having natural coordinates nc.
  void
  GenericElement::interpolate(const GeoPrim::CVector f[],
			      const GeoPrim::CVector &nc,
			      GeoPrim::CVector &v) const
  {
    const double xi = nc.x();
    const double eta = nc.y();
    const double zeta = nc.z();
    switch(_size) {
    case 4: {
      v = f[0];
      v += (((f[1]-f[0])*=xi) += ((f[2]-f[0])*=eta) += ((f[3] - f[0])*=zeta));
      //    v += (f[1]-f[0])*=xi+=(f[2]-f[0])*=eta+=(f[3]-f[0])*=zeta;
      break;
    }
    case 10: {
      const double alpha = (1.-xi-eta-zeta);
      v = (alpha*(1.-2.*(xi+eta+zeta))*f[0] +
	   xi*(2.*xi-1.)*f[1] +
	   eta*(2.*eta-1.)*f[2] +
	   zeta*(2.*zeta-1.)*f[3] +
	   4.*xi*alpha*f[4] +
	   4.*eta*alpha*f[5] +
	   4.*zeta*alpha*f[6] +
	   4.*xi*eta*f[7] +
	   4.*eta*zeta*f[8] +
	   4.*zeta*xi*f[9]);
      break;
    }
    case 8: {
      const double xi = nc.x();
      const double xi_minus = 1. - xi;
      const double eta = nc.y();
      const double eta_minus = 1. - eta;
      const double zeta = nc.z();
      const double zeta_minus = 1. - zeta;
      v = (xi_minus*eta_minus*zeta_minus*f[0] +
	   xi*eta_minus*zeta_minus*f[1] +
	   xi*eta*zeta_minus*f[2] +
	   xi_minus*eta*zeta_minus*f[3] +
	   xi_minus*eta_minus*zeta*f[4] +
	   xi*eta_minus*zeta*f[5] +
	   xi*eta*zeta*f[6] +
	   xi_minus*eta*zeta*f[7]);
      break;
    }
    default:
     std::cerr << "interpolate::error Cannot handle this element "
	   << "type (yet)." << std::endl;
      exit(1);
    }
  }

///
/// \brief Bounding boxes for a mesh
///
/// This function will determine some bounding boxes for the mesh.
/// Namely, it gets a box for the entire mesh, one for the smallest
/// element and one for the largest.
  void
  GetMeshBoxes(const NodalCoordinates &nc,
	       const Connectivity &ec,
	       GeoPrim::CBox &mesh_box,
	       GeoPrim::CBox &small_box,
	       GeoPrim::CBox &large_box)
  {
    mesh_box.init(nc[1],nc.size());
    small_box = mesh_box;
    IndexType nelem = ec.Nelem();
    IndexType n = 1;
    while(n <= nelem){
      unsigned int esize = ec.Esize(n);
      std::vector<GeoPrim::CPoint> element_points(esize);
      // Make a vector of the element points
      for(IndexType node = 1;node <= esize;node++){
	GeoPrim::CPoint ep(nc[ec.Node(n,node)]);
	element_points[node-1] = ep;
      }
      // Make the bounding box for the element
      GeoPrim::CBox box(element_points);
      if(!box.empty()){
	if(box < small_box) small_box = box;
	if(box > large_box) large_box = box;
      }
      n++;
    }
  }


  /// \brief Get elements in box
  ///
  /// Given a box in a cartesian space, this function will populate a
  /// list<IndexType> with all the elements that lie within.
  /// \n
  /// Inputs: CBox object (see GeoPrimitives.H)
  ///         NodalCoordinates object
  ///         ElementConnectivity object
  ///         DualConnectivity object
  ///         list<IndexType> stores results
  ///
  ///
  void
  FindElementsInBox(const GeoPrim::CBox &box,
		    const NodalCoordinates &nc,
		    const Connectivity &dc, // dual connectivity
		    std::list<IndexType> &elements)
  {
    int nnodes = nc.Size();
    for(int i = 1;i <= nnodes;i++){
      GeoPrim::CPoint p(nc[i]);
      if(box.contains(p)){
	std::vector<IndexType>::const_iterator dci = dc[i-1].begin();
	while(dci != dc[i-1].end()){
	  elements.push_back(*dci++);
	}
      }
    }
    if(!elements.empty()){
      elements.sort();
      elements.unique();
    }
  }


  ///
  /// \brief Locate element containing given physical point
  ///
  /// This function will locate which element in a mesh contains
  /// a specified point by constructing a box around the point and
  /// solving (by Newton-Raphson) the system for each element that
  /// touches the box until the element is found.
  /// \n
  ///
  /// Inputs: A point in cartesian, ie (X,Y,Z), coordinates
  ///         A Nodal Coordinates object for the mesh
  ///         An Element Connectivity object for the mesh
  ///         Dual Connectivity object for the mesh
  ///         CBox specifying the parameters of the box use
  ///
  /// Returns: An integer indicating which element in the
  ///          connectivity object contains the point. A
  ///          0 will indicate that the point could not
  ///          be located.
  ///
  ///          The natural coordinates of the point in
  ///          the containing element is returned in natc
  ///
  Mesh::IndexType
  FindPointInCells(const GeoPrim::CPoint &p, // Target point
		   const NodalCoordinates &nc, // Source
		   const Connectivity &ec, // Source connectivity
		   const std::vector<Mesh::IndexType> &elements,    // candidate cells
		   GeoPrim::CVector &natc)  // Returns Targ nat
  {

    //    GeoPrim::CBox bounds(box.around(p));
    //    std::list<IndexType> elements;
    //    FindElementsInBox(bounds, nc, dc, elements);
    std::vector<IndexType>::const_iterator ei = elements.begin();
    while(ei != elements.end()){
      GeoPrim::CVector guess;
      unsigned int esize = ec.Esize(*ei);
      if(esize == 4 || esize == 10)
	guess.init(.25,.25,.25);
      else if (esize == 8 || esize == 20)
	guess.init(.5,.5,.5);
      else
	assert(false);
      natc.init(guess);
      // Solve the non-linear system using newton-raphson with an
      // initial guess as the center of the theoretical element.
      if(!NewtonRaphson(natc,*ei,GenericElement(esize),ec,nc,p)){
	//	double dist = 0.0;
	//	Mesh::IndexType closest_node = nc.closest_node(p,&dist);
	//	std::cout << "Mesh::FindPointInMesh: Error: Closest approach: " << dist
	//		  << " at node " << closest_node << "." << std::endl
	//		  << "Mesh::FindPointInMesh: Element(" << *ei << ") = (";
	//	IRAD::Util::DumpContents(std::cout,ec[*ei-1],",");
	//	std::cout << ")" << std::endl;
	//	std::cout << "Point: " << p << std::endl
	//	          << "Nodes: " << std::endl;
	//	std::vector<Mesh::IndexType>::const_iterator ni = ec[*ei-1].begin();
	//	while(ni != ec[*ei-1].end()){
	//	  std::cout << ni - ec[*ei-1].begin() << ": " << GeoPrim::C3Point(nc[*ni]) << std::endl;
	//	  ni++;
	//	}
	//	assert(false);
	return(0);
      }
      if(natc[0] >= LTOL && natc[0] <= HTOL &&
	 natc[1] >= LTOL && natc[1] <= HTOL &&
	 natc[2] >= LTOL && natc[2] <= HTOL){
	if(esize == 4 || esize == 10){
	  if((natc[0]+natc[1]+natc[2]) <= HTOL)
	    return (*ei);
	}
	else if(esize == 8 || esize == 20)
	  return(*ei);
	else{
	  std::cerr << "Mesh::FindPointInMesh: Error: Cannot handle"
		    << " this element type. (yet)" << std::endl;
	  exit(1);
	}
      }
      ei++;
    }
    return (0);
  }

  ///
  /// \brief Locate element containing given physical point
  ///
  /// This function will locate which element in a mesh contains
  /// a specified point by constructing a box around the point and
  /// solving (by Newton-Raphson) the system for each element that
  /// touches the box until the element is found.
  /// \n
  ///
  /// Inputs: A point in cartesian, ie (X,Y,Z), coordinates
  ///         A Nodal Coordinates object for the mesh
  ///         An Element Connectivity object for the mesh
  ///         Dual Connectivity object for the mesh
  ///         CBox specifying the parameters of the box use
  ///
  /// Returns: An integer indicating which element in the
  ///          connectivity object contains the point. A
  ///          0 will indicate that the point could not
  ///          be located.
  ///
  ///          The natural coordinates of the point in
  ///          the containing element is returned in natc
  ///
  Mesh::IndexType
  FindPointInMesh(const GeoPrim::CPoint &p, // Target Mesh point
		  const NodalCoordinates &nc, // Source
		  const Connectivity &ec, // Source connectivity
		  const Connectivity &dc,    // Source dual connectivity
		  const GeoPrim::CBox &box,    // neigborhood (typically defined by some source character)
		  GeoPrim::CVector &natc)  // Returns Targ nat
  {

    GeoPrim::CBox bounds(box.around(p));
    std::list<IndexType> elements;
    FindElementsInBox(bounds, nc, dc, elements);
    std::list<IndexType>::iterator ei = elements.begin();
    while(ei != elements.end()){
      GeoPrim::CVector guess;
      unsigned int esize = ec.Esize(*ei);
      if(esize == 4 || esize == 10)
	guess.init(.25,.25,.25);
      else if (esize == 8 || esize == 20)
	guess.init(.5,.5,.5);
      natc.init(guess);
      // Solve the non-linear system using newton-raphson with an
      // initial guess as the center of the theoretical element.
      if(!NewtonRaphson(natc,*ei,GenericElement(esize),ec,nc,p)){
	double dist = 0.0;
	Mesh::IndexType closest_node = nc.closest_node(p,&dist);
	std::cout << "Mesh::FindPointInMesh: Error: Closest approach: " << dist
		  << " at node " << closest_node << "." << std::endl
		  << "Mesh::FindPointInMesh: Element(" << *ei << ") = (";
	IRAD::Util::DumpContents(std::cout,ec[*ei-1],",");
	std::cout << ")" << std::endl;
	std::vector<Mesh::IndexType>::const_iterator ei2 = dc[closest_node-1].begin();
	GeoPrim::CBox new_bounds;
	while(ei2 != dc[closest_node-1].end()){
	  // For every element touching this node
	  std::vector<Mesh::IndexType>::const_iterator ni = ec[*ei2-1].begin();
	  while(ni != ec[*ei2-1].end())
	    new_bounds.AddPoint(nc[*ni++]);
	  ei2++;
	}
	std::cout << "Mesh::FindPointInMesh: Bounding box of failed search: "
		  << new_bounds << std::endl;
	return(0);
      }
      if(natc[0] >= LTOL && natc[0] <= HTOL &&
	 natc[1] >= LTOL && natc[1] <= HTOL &&
	 natc[2] >= LTOL && natc[2] <= HTOL){
	if(esize == 4 || esize == 10){
	  if((natc[0]+natc[1]+natc[2]) <= HTOL)
	    return (*ei);
	}
	else if(esize == 8 || esize == 20)
	  return(*ei);
	else{
	  std::cerr << "Mesh::FindPointInMesh: Error: Cannot handle"
		    << " this element type. (yet)" << std::endl;
	  exit(1);
	}
      }
      ei++;
    }
    return (0);
  }


  ///
  /// \brief Locate element containing given physical point
  ///
  /// This function will locate which element in a mesh contains
  /// a specified point by constructing a box around the point and
  /// solving (by Newton-Raphson) the system for each element that
  /// touches the box until the element is found.
  /// \n
  ///
  /// Inputs: A point in cartesian, ie (X,Y,Z), coordinates
  ///         A Nodal Coordinates object for the mesh
  ///         An Element Connectivity object for the mesh
  ///         Dual Connectivity object for the mesh
  ///         CBox specifying the parameters of the box use
  ///
  /// Returns: An integer indicating which element in the
  ///          connectivity object contains the point. A
  ///          0 will indicate that the point could not
  ///          be located.
  ///
  ///          The natural coordinates of the point in
  ///          the containing element is returned in natc
  ///
  Mesh::IndexType
  FindPointInMesh_2(const GeoPrim::CPoint &p, // Target Mesh point
		    const NodalCoordinates &nc, // Source
		    const Connectivity &ec, // Source connectivity
		    const Connectivity &dc,    // Source dual connectivity
		    const GeoPrim::CBox &box,    // Source
		    GeoPrim::CVector &natc)  // Returns Targ nat
  {

    GeoPrim::CBox bounds(box.around(p));
    std::list<IndexType> elements;
    FindElementsInBox(bounds, nc, dc, elements);
    std::list<IndexType>::iterator ei = elements.begin();
    while(ei != elements.end()){
      GeoPrim::CVector guess;
      unsigned int esize = ec.Esize(*ei);
      if(esize == 3 || esize == 6)
	guess.init(.25,.25,0);
      else if (esize == 4 || esize == 8)
	guess.init(.5,.5,0);
      natc.init(guess);
      // Solve the non-linear system using newton-raphson with an
      // initial guess as the center of the theoretical element.
      if(!NewtonRaphson_2(natc,*ei,GenericCell_2(esize),ec,nc,p)){
	double dist = 0.0;
	Mesh::IndexType closest_node = nc.closest_node(p,&dist);
	std::cout << "Mesh::FindPointInMesh: Error: Closest approach: " << dist
		  << " at node " << closest_node << "." << std::endl
		  << "Mesh::FindPointInMesh: Element(" << *ei << ") = (";
	IRAD::Util::DumpContents(std::cout,ec[*ei-1],",");
	std::cout << ")" << std::endl;
	std::vector<Mesh::IndexType>::const_iterator ei2 = dc[closest_node-1].begin();
	GeoPrim::CBox new_bounds;
	while(ei2 != dc[closest_node-1].end()){
	  // For every element touching this node
	  std::vector<Mesh::IndexType>::const_iterator ni = ec[*ei2-1].begin();
	  while(ni != ec[*ei2-1].end())
	    new_bounds.AddPoint(nc[*ni++]);
	  ei2++;
	}
	std::cout << "Mesh::FindPointInMesh: Bounding box of failed search: "
		  << new_bounds << std::endl;
	return(0);
      }
      if(natc[0] >= LTOL && natc[0] <= HTOL &&
	 natc[1] >= LTOL && natc[1] <= HTOL)
	{
	  if(esize == 3 || esize == 6){
	    if(((natc[0]+natc[1]) <= HTOL))
	      return (*ei);
	  }
	  else if(esize == 4 || esize == 8)
	    return(*ei);
	  else{
	    std::cerr << "Mesh::FindPointInMesh: Error: Cannot handle"
		      << " this element type. (yet)" << std::endl;
	    exit(1);
	  }
	}
      ei++;
    }
    return (0);
  }


  GeoPrim::C3Point
  GenericElement::Centroid(std::vector<IndexType> &ec,NodalCoordinates &nc) const
  {
    GeoPrim::C3Point centroid(0,0,0);
    IndexType i = 0;
    IndexType esize = ec.size();
    while(i < esize)
      centroid += GeoPrim::C3Point(nc[ec[i++]]);
    double scal = 1.0/(static_cast<double>(esize));
    return(centroid *= scal);
  };

  void
  GenericElement::Centroid(std::vector<IndexType> &ec,NodalCoordinates &nc,GeoPrim::C3Point &centroid) const
  {
    centroid.Init();
    IndexType i = 0;
    IndexType esize = ec.size();    
    while(i < esize)
      centroid += GeoPrim::C3Point(nc[ec[i++]]);
    double scal = 1.0/(static_cast<double>(esize));
    centroid *= scal;
  };

  bool
  GenericElement::Inverted(std::vector<IndexType> &ec,NodalCoordinates &nc) const
  {
    // Form the vector V1 from element centroid to first face centroid
    // Dot V1 with the first face normal
    // return(result < 0)
    if(_size == 4 || _size == 10){ // Tet Faces:      132 241 342 143
      GeoPrim::C3Point p1(nc[ec[0]]);
      GeoPrim::C3Point p2(nc[ec[2]]);
      GeoPrim::C3Point p3(nc[ec[1]]);
      GeoPrim::C3Facet aface(&p1,&p2,&p3);
      GeoPrim::C3Vector avec(Centroid(ec,nc),aface.Centroid());
      return((avec*aface.Normal()) < 0);
    }
    else if(_size == 5 || _size == 13){            // Pyr Faces:      1432 251 352 453 154
      GeoPrim::C3Point p1(nc[ec[0]]);
      GeoPrim::C3Point p2(nc[ec[3]]);
      GeoPrim::C3Point p3(nc[ec[2]]);
      GeoPrim::C3Point p4(nc[ec[1]]);
      GeoPrim::C3Facet bface(&p1,&p2,&p3,&p4);
      GeoPrim::C3Vector bvec(Centroid(ec,nc),bface.Centroid());
      return((bvec*bface.Normal()) < 0);
    }
    else if(_size == 6 || _size == 15){            // Prism Faces:    2541 3652 1463 132 456
      GeoPrim::C3Point p1(nc[ec[1]]);
      GeoPrim::C3Point p2(nc[ec[4]]);
      GeoPrim::C3Point p3(nc[ec[3]]);
      GeoPrim::C3Point p4(nc[ec[0]]);
      GeoPrim::C3Facet cface(&p1,&p2,&p3,&p4);
      GeoPrim::C3Vector cvec(Centroid(ec,nc),cface.Centroid());
      return((cvec*cface.Normal()) < 0);
    }
    else if(_size == 8 || _size == 20){            // Hex Faces:      1432 2651 3762 4873 1584 5678
      GeoPrim::C3Point p1(nc[ec[0]]);
      GeoPrim::C3Point p2(nc[ec[3]]);
      GeoPrim::C3Point p3(nc[ec[2]]);
      GeoPrim::C3Point p4(nc[ec[1]]);
      GeoPrim::C3Facet dface(&p1,&p2,&p3,&p4);
      GeoPrim::C3Vector dvec(Centroid(ec,nc),dface.Centroid());
      return((dvec*dface.Normal()) < 0);
    }
    else
      return(false);
    return(false);
  };

  bool
  GenericElement::ShapeOK(std::vector<IndexType> &ec,NodalCoordinates &nc) const
  {
    // ensure element is not poorly shaped
    //  make sure point 3 is not co-linear with points 1 and 2
    //  make sure 4th point not co-planar with base
    if(_size == 4){ // Tet Faces:      132 241 342 143
      GeoPrim::C3Point p1(nc[ec[0]]);
      GeoPrim::C3Point p2(nc[ec[2]]);
      GeoPrim::C3Point p3(nc[ec[1]]);
      GeoPrim::C3Point p4(nc[ec[3]]);
      GeoPrim::C3Plane plane(p1,p2,p3);
      return(!plane.contains_point(p4));
    }
    // ensure element in not poorly shaped
    //   make sure quad base is planar
    //   make sure area of quad base is positive
    //   make sure 5th point not co-planar with quad base
    else if(_size == 5){            // Pyr Faces:      1432 251 352 453 154
      GeoPrim::C3Point p1(nc[ec[0]]);
      GeoPrim::C3Point p2(nc[ec[3]]);
      GeoPrim::C3Point p3(nc[ec[2]]);
      GeoPrim::C3Point p4(nc[ec[1]]);
      GeoPrim::C3Point p5(nc[ec[4]]);
      GeoPrim::C3Facet bface(&p1,&p2,&p3,&p4);
      GeoPrim::C3Vector bvec(Centroid(ec,nc),bface.Centroid());
      return((bvec*bface.Normal()) < 0);
    }
    // ensure element is not poorly shaped
    //   make sure both triangular faces are triangles (i.e. not co-linear)
    //   make sure no points of opposing face is co-planar with base
    //   make sure quad faces are planar
    //   make sure quad faces have positive area
    else if(_size == 6){            // Prism Faces:    2541 3652 1463 132 456
      GeoPrim::C3Point p1(nc[ec[1]]);
      GeoPrim::C3Point p2(nc[ec[4]]);
      GeoPrim::C3Point p3(nc[ec[3]]);
      GeoPrim::C3Point p4(nc[ec[0]]);
      GeoPrim::C3Facet cface(&p1,&p2,&p3,&p4);
      GeoPrim::C3Vector cvec(Centroid(ec,nc),cface.Centroid());
      return((cvec*cface.Normal()) < 0);
    }
    // ensure element is not poorly shaped
    //   make sure each face is planar
    //   make sure each face has positive area
    else if(_size == 8){            // Hex Faces:      1432 2651 3762 4873 1584 5678
      GeoPrim::C3Point p1(nc[ec[0]]);
      GeoPrim::C3Point p2(nc[ec[1]]);
      GeoPrim::C3Point p3(nc[ec[2]]);
      GeoPrim::C3Point p4(nc[ec[3]]);

      GeoPrim::C3Point p5(nc[ec[4]]);
      GeoPrim::C3Point p6(nc[ec[5]]);
      GeoPrim::C3Point p7(nc[ec[6]]);
      GeoPrim::C3Point p8(nc[ec[7]]);
      
      GeoPrim::C3Plane plane1(p1,p4,p3);
      if(!plane1.contains_point(p2))
	return(false);
      GeoPrim::C3Plane plane2(p2,p6,p5);
      if(!plane2.contains_point(p1))
	return(false);
      GeoPrim::C3Plane plane3(p3,p7,p6);
      if(!plane3.contains_point(p2))
	return(false);
      GeoPrim::C3Plane plane4(p4,p8,p7);
      if(!plane4.contains_point(p3))
	return(false);
      GeoPrim::C3Plane plane5(p1,p5,p8);
      if(!plane5.contains_point(p4))
	return(false);
      GeoPrim::C3Plane plane6(p5,p6,p7);
      if(!plane6.contains_point(p8))
	return(false);
    }
    else
      return(true);
    return(true);
  };

  void
  GenericElement::ReOrient(std::vector<IndexType> &ec)
  {
    IndexType temp;
    switch(_size){
    case 4:   // Tet Faces:      132 241 342 143
    case 10:
      temp = ec[1];
      ec[1] = ec[0];
      ec[0] = temp;
      break;
    case 5:   // Pyr Faces:      1432 251 352 453 154
    case 13:
      temp = ec[2];
      ec[2] = ec[0];
      ec[0] = ec[2];
      break;
    case 6:   // Prism Faces:    2541 3652 1463 132 456
    case 15:
      temp = ec[1];
      ec[1] = ec[0];
      ec[0] = ec[1];
      temp = ec[4];
      ec[4] = ec[3];
      ec[3] = temp;
      break;
    case 8:   /// Hex Faces:      1432 2651 3762 4873 1584 5678
    case 20:
      temp = ec[2];
      ec[2] = ec[0];
      ec[0] = ec[2];
      temp = ec[6];
      ec[6] = ec[4];
      ec[4] = temp;
      break;
    }
  };



  // Searches _all_ elements one by one for a given point (last resort)
  Mesh::IndexType
  GlobalFindPointInMesh(const GeoPrim::CPoint &p,      // Target Mesh point
			const NodalCoordinates &nc,    // Source
			const Connectivity &ec, // Source
			const Connectivity &dc,    // Source
			const GeoPrim::CBox &box,      // Source
			GeoPrim::CVector &natc)  // Returns Targ nat
  {
    GeoPrim::CVector guess;
    int ein = 1;
    int ein_size = ec.Nelem();
    while(ein <= ein_size){
      unsigned int esize = ec.Esize(ein);
      if(esize == 4 || esize == 10)
	guess.init(.25,.25,.25);
      else if (esize == 8 || esize == 20)
	guess.init(.5,.5,.5);
      natc.init(guess);
      // Solve the non-linear system using newton-raphson with an
      // initial guess as the center of the theoretical element.
      if(!NewtonRaphson(natc,ein,GenericElement(esize),ec,nc,p)){
	std::cerr << "GlobalFindPointInMesh: error NewtonRaphson failed."
		  << std::endl;
	return(0);
      }
      if(natc[0] >= LTOL && natc[0] <= HTOL &&
	 natc[1] >= LTOL && natc[1] <= HTOL &&
	 natc[2] >= LTOL && natc[2] <= HTOL){
	if(esize == 4 || esize == 10){
	  if((natc[0]+natc[1]+natc[2]) <= HTOL){
	    return (ein);
	  }
	}
	else if(esize == 8 || esize == 20)
	  return(ein);
	else{
	  std::cerr << "GlobalFindPointInMesh: Error: Cannot handle"
		    << " this element type. (yet)" << std::endl;
	  exit(1);
	}
      }
      ein++;
    }
    return(0);
  }


  // Does not work for parallel meshes - does not work with BC's.  Need
  // to support T3D mesh formats as well:
  // T3D:
  //     tri,quad,mixed   topological entity counts
  //     [3,4,7]         [0,1,2,2,3,3,3]
  //
  int ReadMesh(const std::string &path,Mesh::UnstructuredMesh &mesh)
  {
    std::ifstream Inf;
    Inf.open(path.c_str());
    if(!Inf)
      return(1);
    Inf >> mesh.nc >> mesh.con;
    if(!Inf)
      return(1);
    Inf.close();
    return(0);
  }

  /// Tet Faces:      132 241 342 143
  /// Pyr Faces:      1432 251 352 453 154
  /// Prism Faces:    2541 3652 1463 132 456
  /// Hex Faces:      1432 2651 3762 4873 1584 5678
  // Input: Single Element Connectivity
  // Output: vector of face connectivities for this element
  void
  GenericElement::get_face_connectivities(Connectivity &rv,
					  const std::vector<IndexType> &e) const
  {
    //    std::vector<std::vector<IndexType> > rv;
    switch(_size){
    case 4:   // Tet Faces:      132 241 342 143
    case 10:
      rv.Resize(4);

      rv[0].resize(3);
      //      std::vector<Mesh::IndexType>(rv[0]).swap(rv[0]);
      rv[0][0] = e[0]; // 1
      rv[0][1] = e[2]; // 3
      rv[0][2] = e[1]; // 2

      rv[1].resize(3);
      //      std::vector<Mesh::IndexType>(rv[1]).swap(rv[1]);
      rv[1][0] = e[1]; // 2
      rv[1][1] = e[3]; // 4
      rv[1][2] = e[0]; // 1

      rv[2].resize(3);
      //      std::vector<Mesh::IndexType>(rv[2]).swap(rv[2]);
      rv[2][0] = e[2]; // 3
      rv[2][1] = e[3]; // 4
      rv[2][2] = e[1]; // 2

      rv[3].resize(3);
      //      std::vector<Mesh::IndexType>(rv[3]).swap(rv[3]);
      rv[3][0] = e[0]; // 1
      rv[3][1] = e[3]; // 4
      rv[3][2] = e[2]; // 3
      break;
    case 5:   // Pyr Faces:      1432 251 352 453 154
    case 13:
      rv.Resize(5);

      rv[0].resize(4);
      //      std::vector<Mesh::IndexType>(rv[0]).swap(rv[0]);
      rv[0][0] = e[0]; // 1
      rv[0][1] = e[3]; // 4
      rv[0][2] = e[2]; // 3
      rv[0][3] = e[1]; // 2

      rv[1].resize(3);
      //      std::vector<Mesh::IndexType>(rv[1]).swap(rv[1]);
      rv[1][0] = e[1]; // 2
      rv[1][1] = e[4]; // 5
      rv[1][2] = e[0]; // 1

      rv[2].resize(3);
      //      std::vector<Mesh::IndexType>(rv[2]).swap(rv[2]);
      rv[2][0] = e[2]; // 3
      rv[2][1] = e[4]; // 5
      rv[2][2] = e[1]; // 2

      rv[3].resize(3);
      //      std::vector<Mesh::IndexType>(rv[3]).swap(rv[3]);
      rv[3][0] = e[3]; // 4
      rv[3][1] = e[4]; // 5
      rv[3][2] = e[2]; // 3

      rv[4].resize(3);
      //      std::vector<Mesh::IndexType>(rv[4]).swap(rv[4]);
      rv[4][0] = e[0]; // 1
      rv[4][1] = e[4]; // 5
      rv[4][2] = e[3]; // 4
      break;
    case 6:   // Prism Faces:    2541 3652 1463 132 456
    case 15:
      rv.Resize(5);

      rv[0].resize(4);
      //      std::vector<Mesh::IndexType>(rv[0]).swap(rv[0]);
      rv[0][0] = e[1]; // 2
      rv[0][1] = e[4]; // 5
      rv[0][2] = e[3]; // 4
      rv[0][3] = e[0]; // 1

      rv[1].resize(4);
      //      std::vector<Mesh::IndexType>(rv[1]).swap(rv[1]);
      rv[1][0] = e[2]; // 3
      rv[1][1] = e[5]; // 6
      rv[1][2] = e[4]; // 5
      rv[1][3] = e[1]; // 2

      rv[2].resize(4);
      //      std::vector<Mesh::IndexType>(rv[2]).swap(rv[2]);
      rv[2][0] = e[0]; // 1
      rv[2][1] = e[3]; // 4
      rv[2][2] = e[5]; // 6
      rv[2][3] = e[2]; // 3

      rv[3].resize(3);
      //      std::vector<Mesh::IndexType>(rv[3]).swap(rv[3]);
      rv[3][0] = e[0]; // 1
      rv[3][1] = e[2]; // 3
      rv[3][2] = e[1]; // 2

      rv[4].resize(3);
      //      std::vector<Mesh::IndexType>(rv[4]).swap(rv[4]);
      rv[4][0] = e[4]; // 5
      rv[4][1] = e[5]; // 6
      rv[4][2] = e[3]; // 4
      break;
    case 8:   /// Hex Faces:      1432 2651 3762 4873 1584 5678
    case 20:
      rv.Resize(6);

      rv[0].resize(4);
      //      std::vector<Mesh::IndexType>(rv[0]).swap(rv[0]);
      rv[0][0] = e[0]; // 1
      rv[0][1] = e[3]; // 4
      rv[0][2] = e[2]; // 3
      rv[0][3] = e[1]; // 2

      rv[1].resize(4);
      //      std::vector<Mesh::IndexType>(rv[1]).swap(rv[1]);
      rv[1][0] = e[1]; // 2
      rv[1][1] = e[5]; // 6
      rv[1][2] = e[4]; // 5
      rv[1][3] = e[0]; // 1

      rv[2].resize(4);
      //      std::vector<Mesh::IndexType>(rv[2]).swap(rv[2]);
      rv[2][0] = e[2]; // 3
      rv[2][1] = e[6]; // 7
      rv[2][2] = e[5]; // 6
      rv[2][3] = e[1]; // 2

      rv[3].resize(4);
      //      std::vector<Mesh::IndexType>(rv[3]).swap(rv[3]);
      rv[3][0] = e[3]; // 4
      rv[3][1] = e[7]; // 8
      rv[3][2] = e[6]; // 7
      rv[3][3] = e[2]; // 3

      rv[4].resize(4);
      //      std::vector<Mesh::IndexType>(rv[4]).swap(rv[4]);
      rv[4][0] = e[0]; // 1
      rv[4][1] = e[4]; // 5
      rv[4][2] = e[7]; // 8
      rv[4][3] = e[3]; // 4

      rv[5].resize(4);
      //      std::vector<Mesh::IndexType>(rv[5]).swap(rv[5]);
      rv[5][0] = e[4]; // 5
      rv[5][1] = e[5]; // 6
      rv[5][2] = e[6]; // 7
      rv[5][3] = e[7]; // 8
      break;
    default:
      rv.Resize(0);
    }
    //    return(rv);
  };

  ///
  /// Newton-Raphson method, customized for using the GeoPrimitives
  /// and computational mesh constructs
  ///
  /// Parameters:
  /// double natc[] = initial guess at the natural coordinates
  /// IndexType elnum = index of the mesh element to use
  /// const ElementConnectivity &ec = The element connectivity for the mesh
  /// const NodalCoordinates &nc = The nodal coordinates for the mesh
  /// const CPoint &point = the point at which we wish a solution
  ///
  bool
  NewtonRaphson(GeoPrim::CVector &natc,
		IndexType elnum,
		const GenericElement &el,
		const Connectivity &ec,
		const NodalCoordinates &nc,
		const GeoPrim::CPoint &point)
  {
    int k,i;
    int ntrial = 100; // Number of iterations before giving up
    double errx,errf;
    int indx[3];
    double errtol = 1e-12;
    GeoPrim::CVector p;
    GeoPrim::CVector fvec;
    GeoPrim::CVector fjac[3];
    for (k=0;k<ntrial;k++) {
      //      std::cout << k << std::endl;
      el.shapef_jacobian_at(point,natc,elnum,ec,nc,fvec,fjac);
      errf=0.0;
      for (i=0;i<3;i++){
	//	if(fvec[i] < errtol)
	//	  fvec[i] = 0.0;
	errf += fabs(fvec[i]);
      }
      if (errf <= errtol)
	return (true);
      p = -1.0 * fvec;
      //      std::cout << "Jac: "  << fjac[0] << std::endl
      //		<< "Jac: "  << fjac[1] << std::endl
      //		<< "Jac: "  << fjac[2] << std::endl
      //		<< "fv : "  << fvec    << std::endl
      //		<< "nc : "  << natc    << std::endl
      //	<< "err: "  << errf    << std::endl;
      if(!LUDcmp(fjac,indx)){
	//	std::cerr << "NewtonRaphson::error: LUDcmp failed." << std::endl;
	//	std::cout << "Jac: " << fjac[0] << std::endl
	//		  << "Jac: " << fjac[1] << std::endl
	//		  << "Jac: " << fjac[2] << std::endl;
	return(false);
      }
      LUBksb(fjac,indx,p);
      errx=0.0;
      for (i=0;i<3;i++)
	errx += fabs(p[i]);
      natc += p;
      if (errx <= errtol)
	return (true);
    }
    //    std::cerr << "NewtonRaphson::warning: reached maximum iterations" << std::endl;
    return (false);
  }
  ///
  /// Newton-Raphson method, customized for using the GeoPrimitives
  /// and computational mesh constructs
  ///
  /// Parameters:
  /// double natc[] = initial guess at the natural coordinates
  /// IndexType elnum = index of the mesh element to use
  /// const ElementConnectivity &ec = The element connectivity for the mesh
  /// const NodalCoordinates &nc = The nodal coordinates for the mesh
  /// const CPoint &point = the point at which we wish a solution
  ///
  bool
  NewtonRaphson_2(GeoPrim::CVector &natc,
		  IndexType elnum,
		  const GenericCell_2 &el,
		  const Connectivity &ec,
		  const NodalCoordinates &nc,
		  const GeoPrim::CPoint &point)
  {
    int k,i;
    int ntrial = 10; // Number of iterations before giving up
    double errx,errf;
    int indx[3];
    GeoPrim::CVector p;
    GeoPrim::CVector fvec;
    GeoPrim::CVector fjac[3];
    for (k=0;k<ntrial;k++) {
      el.shapef_jacobian_at(point,natc,elnum,ec,nc,fvec,fjac);
      errf=0.0;
      for (i=0;i<3;i++)
	errf += fabs(fvec[i]);
      if (errf <= TOL)
	return (true);
      p = -1.0 * fvec;
      if(!LUDcmp(fjac,indx)){
	//	std::cerr << "NewtonRaphson::error: LUDcmp failed." << std::endl;
	return(false);
      }
      LUBksb(fjac,indx,p);
      errx=0.0;
      for (i=0;i<3;i++)
	errx += fabs(p[i]);
      natc += p;
      if (errx <= TOL)
	return (true);
    }
    //    std::cerr << "NewtonRaphson::warning: reached maximum iterations" << std::endl;
    return (false);
  }


  // LU Decomp
#define LUDTINY 1.0e-24
  bool
  LUDcmp(GeoPrim::CVector a[], int indx[])
  {
    int i,j,k;
    int imax = 0;
    double big,dum,sum,temp;
    GeoPrim::CVector vv;



    for (i=0;i<3;i++) {
      big=0.0;
      for (j=0;j<3;j++)
	if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big == 0.0){
	//	std::cerr << "LUDcmp::error: Singular matrix" << std::endl;
	return(false);
      }
      vv[i]=1.0/big;
    }
    for (j=0;j<3;j++) {
      for (i=0;i<j;i++) {
	sum=a[i][j];
	for (k=0;k<i;k++)
	  sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
      }
      big=0.0;
      for (i=j;i<3;i++) {
	sum=a[i][j];
	for (k=0;k<j;k++)
	  sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
	if ( (dum=vv[i]*fabs(sum)) >= big) {
	  big=dum;
	  imax=i;
	}
      }
      if (j != imax) {
	for (k=0;k<3;k++) {
	  dum=a[imax][k];
	  a[imax][k]=a[j][k];
	  a[j][k]=dum;
	}
	vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=LUDTINY;
      if (j != 3) {
	dum=1.0/(a[j][j]);
	for (i=j+1;i<3;i++)
	  a[i][j] *= dum;
      }
    }
    return (true);
  }
#undef LUDTINY

  // LU Back substitution
  void
  LUBksb(GeoPrim::CVector a[],
	 int indx[],
	 GeoPrim::CVector &b)
  {
    int i,ii=0,ip,j;
    double sum;

    for (i=0;i<3;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii){
	for (j=ii;j<=i-1;j++)
	  sum -= a[i][j]*b[j];
      }
      else if (sum)
	ii=i;
      b[i]=sum;
    }
    for (i=2;i>=0;i--) {
      sum=b[i];
      for (j=i+1;j<3;j++)
	sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
    }
  }

  std::istream &operator>>(std::istream &iSt, Mesh::GeometricEntity &ge)
  {
    std::string line;
    std::getline(iSt,line);
    std::istringstream Istr(line);
    Istr >> ge.first;
    iSt >> ge._collections;
    ge.second = ge._collections.size();
    return(iSt);
  }

  std::istream &operator>>(std::istream &iSt, Mesh::Connectivity &ec)
  {
    std::string line;
    std::getline(iSt,line);
    std::istringstream Istr(line);
    IndexType nelem = 0;
    Istr >> nelem;
    if(nelem > 0){
      ec.Resize(nelem);
      ec._nelem = nelem;
      for(IndexType i = 0;i < nelem;i++){
	std::getline(iSt,line);
	std::istringstream Istr2(line);
	IndexType node;
	ec[i].reserve(8);
	while(Istr2 >> node)
	  ec[i].push_back(node);
	std::vector<Mesh::IndexType>(ec[i]).swap(ec[i]);
      }
    }
    return(iSt);
  }

  std::istream &operator>>(std::istream &iSt, Mesh::NodalCoordinates &nc)
  {
    //  if(nc.size() == 0) return (iSt);
    std::string line;
    IndexType nnodes;
    std::getline(iSt,line);
    std::istringstream Istr(line);
    Istr >> nnodes;
    if(nnodes > 0){
      nc.init(nnodes);
      for(IndexType n = 1;n <= nc.size();n++){
	double *inPtr = nc[n];
	std::getline(iSt,line);
	std::istringstream Istr2(line);
	Istr2 >> inPtr[0] >> inPtr[1] >> inPtr[2];
      }
    }
    return(iSt);
  }

  std::ostream &operator<<(std::ostream &oSt, const Mesh::NodalCoordinates &nc)
  {
    //  if(nc.size() == 0) return(oSt);
    oSt.setf(std::ios::scientific);
    oSt.setf(std::ios::showpoint);
    oSt << std::setprecision(10) << std::setiosflags(std::ios::left);
    oSt << nc.size();
    if(nc.size() > 0) oSt << "\n";
    for(IndexType n = 1;n <= nc.size() ;n++){
      oSt << std::setw(16) << nc[n][0] << "\t"
	  << std::setw(16) << nc[n][1] << "\t"
	  << std::setw(16) << nc[n][2];
      if(n != nc.size())
	oSt << std::endl;
    }
    return (oSt);
  }

  std::ostream &operator<<(std::ostream &oSt, const Mesh::Connectivity &ec)
  {
    oSt << std::setiosflags(std::ios::left) << ec.size();
    if(ec.size() == 0) return(oSt);
    oSt << "\n";
    Mesh::Connectivity::const_iterator ci = ec.begin();
    while(ci != ec.end()){
      std::vector<IndexType>::const_iterator ni = ci->begin();
      while(ni != ci->end()){
	oSt << *ni++;
	if(ni != ci->end())
	  oSt << "\t";
      }
      ci++;
      if(ci != ec.end())
	oSt << "\n";
    }
    return(oSt);
  }

  std::ostream &operator<<(std::ostream &oSt, const Mesh::GeometricEntity &ge)
  {
    if(ge.first.empty()) return(oSt);
    oSt << std::setiosflags(std::ios::left) << ge.first << "\n";
    oSt << ge._collections;
    return(oSt);
  }


} // namespace Mesh

