#include <string>
#include <cmath>
#include <iostream>
#include <sstream>

#include "Rocon.H"
#include "TRAIL.H"
#include "Parameters.H"

COM_EXTERN_MODULE(SimIN);

/// Init Rocon from a given input mesh
void Rocon::init_from_file( const char *inp, const int *ndiv )
{
  if(!inp){
    std::cout << "Rocon::Error: Called without an input file." << std::endl;
    return;
  }
  std::string configfile(inp);
  ifstream ConfigFile;
  ConfigFile.open(configfile.c_str());
  if(!ConfigFile){
    std::cout << "Rocon::Error: Could not open configuration file, " 
	      << configfile << std::endl;
    return;
  }
  IRAD::Util::Parameters configparams;
  configparams.ReadFromStream(ConfigFile);
  ConfigFile.close();
  std::string inpath(configparams.Param("constraint_surface"));
  int verbosity = configparams.GetValue<int>("verbosity");
  _transform_x = 1.0;
  _transform_y = 1.0;
  _transform_z = 1.0;
  _transform = false;
  if(!configparams.Param("transform_x").empty()){
    _transform_x = configparams.GetValue<double>("transform_x");
    _transform = true;
  }
  if(!configparams.Param("transform_y").empty()){
    _transform_y = configparams.GetValue<double>("transform_y");
    _transform = true;
  }
  if(!configparams.Param("transform_z").empty()){
    _transform_z = configparams.GetValue<double>("transform_z");
    _transform = true;
  }
  if(!configparams.Param("transform").empty()){
    _transform_x = _transform_y = _transform_z = configparams.GetValue<double>("transform");
    _transform = true;
  }
  if(verbosity > 1)
    std::cout << "Rocon: Reading " << inpath << std::endl;
  int numdiv = 100;
  if(ndiv)
    numdiv = *ndiv;
  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocin, "PROPCONIN");
  int IN_read;
  IN_read = COM_get_function_handle( "PROPCONIN.read_window");
  MPI_Comm comm_null = MPI_COMM_NULL;  
  std::string bufwin("bufwin");
  COM_call_function( IN_read, inpath.c_str(), bufwin.c_str(), &comm_null);
  int IN_obtain   = COM_get_function_handle( "PROPCONIN.obtain_dataitem");
  int mesh_handle = COM_get_dataitem_handle((bufwin+".mesh").c_str());
  COM_call_function( IN_obtain, &mesh_handle, &mesh_handle);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC( Rocin, "PROPCONIN");
  COM_window_init_done(bufwin.c_str());
  int init_handle = COM_get_function_handle((_wname+".initialize").c_str());
  COM_get_dataitem_handle((bufwin+".mesh").c_str());
  COM_call_function(init_handle,&mesh_handle,&numdiv);
  COM_delete_window(bufwin.c_str());
}

void CheckCoordinates(double *coords,int number_of_points,double tol)
{
  for(int i = 0; i < number_of_points*3;i++)
    {
      if(std::fabs(coords[i]) < tol && (coords[i] != 0.0)){
	//	std::cout << "Rocon> WARNING: node " << (i/3)+1 << "/" 
	//		  << number_of_points << " "   
	//	  << (((i%3)==0 ? "x " : ((i%3) < 2) ? "y " : "z "))
	//		  << "coordinate = " << coords[i] << std::endl;
	coords[i] = 0.0;
      }
    }
}
void TransformCoordinates(double *coords,int number_of_points,double xfac,double yfac,double zfac)
{
  for(int i = 0; i < number_of_points;i++)
    {
      coords[i*3] *= xfac;
      coords[i*3 + 1] *= yfac;
      coords[i*3 + 2] *= zfac;
    }
}

/// Initialize Rocon with given mesh.
void Rocon::initialize( const COM::DataItem *pmesh,const int *ndiv)
{
  COM_assertion_msg( validate_object()==0, "Invalid object");
  COM_assertion_msg( pmesh, "Mesh must be present");
  _checked = false;

  COM::Window *win = const_cast<COM::Window*>(pmesh->window());
  this->_comm = win->get_communicator();
  if ( COMMPI_Initialized())  this->_rank = COMMPI_Comm_rank(this->_comm);
  else this->_rank = 0;

  if(false){
    double test_coords[21];
    test_coords[0] = test_coords[1] = test_coords[2] = 0.0;
    test_coords[3] = test_coords[7] = 1.0;
    test_coords[4] = test_coords[5] = test_coords[8] = 0.0;
    test_coords[6] = 0.5;
    test_coords[9] = test_coords[11] = test_coords[14] = 0.0;
    test_coords[10] = test_coords[12] = test_coords[13] = 1.0;
    test_coords[15] = test_coords[18] = test_coords[19] = 0.0;
    test_coords[16] = 1.0;
    test_coords[17] = test_coords[20] = -1.0;
    double test_point[3];
    test_point[0] = test_point[1] = 0.5;
    test_point[2] = 0.0;
    double test_disp[3];
    test_disp[0] = test_disp[1] = 0.0;
    test_disp[2] = -1.0;
    int test_conn[15];
    test_conn[0] = 1;
    test_conn[1] = 0;
    test_conn[2] = 2;
    test_conn[3] = 3;
    test_conn[4] = 0;
    test_conn[5] = 2;
    test_conn[6] = 2;
    test_conn[7] = 1;
    test_conn[8] = 4;
    test_conn[9] = 5;
    test_conn[10] = 0;
    test_conn[11] = 3;
    test_conn[12] = 5;
    test_conn[13] = 6;
    test_conn[14] = 0;
    double test_position[3];
    MeshBndSurf testmbs;
    testmbs.initialize(test_coords,test_conn,7,5,1);
    if(testmbs.intersection(test_point,test_disp,test_position) && !_rank)
      std::cout << "Rocon> Passed the test, found <" << test_position[0] << ","
		<< test_position[1] << "," << test_position[2] << " >" << std::endl;
    else if(!_rank)
      std::cout << "Rocon> FAILED THE TEST" << std::endl;
  }


  COM_new_window("tempwin");
  // ensures contiguous layout
  COM_clone_dataitem("tempwin.mesh", (win->name()+".mesh").c_str(), 0);
  COM_window_init_done("tempwin");
  
  int    *connectivity      = NULL;
  double *nodal_coordinates = NULL;
  std::vector<int> pane_ids;
  COM_get_panes("tempwin",pane_ids);
  COM_assertion_msg((pane_ids.size() == 1),"Too many panes in input mesh.");
  COM_get_array("tempwin.nc",pane_ids[0],&nodal_coordinates);
  COM_assertion_msg((nodal_coordinates != NULL),
		    "Failed to extract coordinates from input.");
  COM_get_array("tempwin.:t3:real",pane_ids[0],&connectivity);
  COM_assertion_msg((connectivity != NULL),
		    "Failed to extract surface triangles from input.");
  int number_of_divisions = *ndiv;
  int number_of_nodes = 0;
  int number_of_elements = 0;
  COM_get_size("tempwin.nc", pane_ids[0],&number_of_nodes);
  COM_get_size("tempwin.:t3:real", pane_ids[0],&number_of_elements);
  if(_transform){
    if(!_rank)
      std::cout << "Rocon> Transforming coordinates" << std::endl;
    TransformCoordinates(nodal_coordinates,number_of_nodes,
			 _transform_x,_transform_y,_transform_z);
  }
  if(!_rank){
    std::cout << "Rocon> Checking input coordinates. " << std::endl;
  }
  
  CheckCoordinates(nodal_coordinates,number_of_nodes,1.0e-9);
  // Connectivity from Window is 1-based, change it to 0-based.
  for(int i = 0;i < number_of_elements*3;i++)
    connectivity[i] -= 1;
  this->_mbs = new MeshBndSurf;
  this->_mbs->initialize(nodal_coordinates,connectivity,number_of_nodes,
			 number_of_elements,number_of_divisions);
  if(!_rank)
    std::cout << "Rocon> Initialized constraint surface with " 
	      << number_of_nodes << " nodes and " 
	      << number_of_elements << " triangles. (ndiv = " 
	      << *ndiv << ")" << std::endl;
  if(false){
    double test_point[6];
    test_point[0] = test_point[1] = 0.0;
    test_point[2] = 0.0599999999999;
    test_point[3] = test_point[4] = 0.0;
    test_point[5] = 0.06;
    double test_disp[6];
    test_disp[0] = test_disp[1] = 0.0;
    test_disp[2] = 2e-9;
    test_disp[3] = 2e-3;
    test_disp[4] = 0.0;
    test_disp[5] = 2e-3;
    double test_position[3];
    if(this->_mbs->intersection(test_point,test_disp,test_position) && !_rank)
      std::cout << "Rocon> Passed the test2, found <" << test_position[0] << ","
		<< test_position[1] << "," << test_position[2] << " >" << std::endl;
    else if(!_rank)
      std::cout << "Rocon> FAILED THE TEST 2" << std::endl;

    if(this->_mbs->intersection(&(test_point[3]),&(test_disp[3]),test_position) && !_rank)
      std::cout << "Rocon> Passed the test3, found <" << test_position[0] << ","
		<< test_position[1] << "," << test_position[2] << " >" << std::endl;
    else if(!_rank)
      std::cout << "Rocon> FAILED THE TEST 3" << std::endl;

  }    
  COM_delete_window("tempwin");
}

/// set bflag to 0 for every element whos nodes all
/// have cflag = 1.
void Rocon::burnout( const COM::DataItem *pmesh, 
		     const COM::DataItem *cflag,
		     COM::DataItem *bflag)
{
  COM_assertion_msg( validate_object()==0, "Invalid object");
  COM_assertion_msg( pmesh, "Mesh must be present");
  COM_assertion_msg( cflag,  "Cflag must be present");
  COM_assertion_msg( bflag,   "Bflag must be present");

  COM::Window *meshwin = const_cast<COM::Window*>(pmesh->window());
  COM::Window *cflagwin = const_cast<COM::Window*>(cflag->window());
  COM::Window *bflagwin  = const_cast<COM::Window*>(bflag->window());

  std::string meshwinname(meshwin->name());
  std::string cflagwinname(cflagwin->name());
  std::string bflagwinname(bflagwin->name());
  
  std::string cflagname(cflag->name());
  std::string bflagname(bflag->name());

  unsigned int number_of_panes = 0;
  std::vector<int> pane_ids;
  COM_get_panes(meshwinname.c_str(),pane_ids);
  number_of_panes = pane_ids.size();
  int total_burned_out = 0;
  int local_burned_out = 0;

  // First, make sure that nodes on nonburning panes are
  // seen as "burned out" by the burnout processing.  
  // We need to be sure and remove this 
  std::vector<int>::iterator pii = pane_ids.begin();
  while(pii != pane_ids.end()){
    int pane_id = *pii++;
    int *bcflag = NULL;
    int nreal = 0;
    int nghost = 0;
    COM_get_array((cflagwinname+".bcflag").c_str(),pane_id,&bcflag);
    if(bcflag != NULL){
      if(*bcflag == 0 || *bcflag == 2){
	int *constrained = NULL;
	COM_get_array((cflagwinname+"."+cflagname).c_str(),pane_id,&constrained);
	if(constrained){
	  COM_get_size((cflagwinname+"."+cflagname).c_str(),pane_id,&nreal,&nghost);
	  if(nreal > 0){
	    for(int i = 0;i < nreal;i++)
	      constrained[i] = 1;
	  }
	}
      }	
    }
  }
  pii = pane_ids.begin();
  //  std::vector<int>::iterator pii = pane_ids.begin();
  while(pii != pane_ids.end()){
    int    *connectivity          = NULL;
    int    *constrained           = NULL;
    int    *burning               = NULL;
    int nreal                     = 0;
    int nghost                    = 0;
    int pane_id = *pii++;
    COM_get_array((cflagwinname+"."+cflagname).c_str(),pane_id,&constrained);
    if(constrained == NULL){
      if(!_rank && _verbose)
	std::cout << "Rocon> Skipping pane " << pane_id 
		  << " which has no cflag array."
		  << std::endl;
      continue;
    }
    int bflag_size   = 0;
    int bflag_ghosts = 0;
    int bflag_stride = 0;
    int bflag_cap    = 0;
    int *bcflag      = NULL;
    COM_get_array((bflagwinname+"."+bflagname).c_str(),pane_id,
		  &burning,&bflag_stride,&bflag_cap);
    //    COM_get_array((bflagwinname+"."+bflagname).c_str(),pane_id,
    //		  &bcflag);
    bool burning_pane = (burning != NULL);
    //    if(bcflag && burning)
    //      if(*bcflag == 1)
    //	burning_pane = true;
    //    if(!burning_pane){
    //      if(!_rank && _verbose)
    //	std::cout << "Rocon> Skipping pane " << pane_id 
    //		  << " which has no matching bflag array." << std::endl;
    //      continue;
    //    }
    if(burning_pane)
      COM_get_size((bflagwinname+"."+bflagname).c_str(),pane_id,
		   &bflag_size,&bflag_ghosts);
    // has a cflag array, should be processed
    Mesh::Connectivity pane_conn;
    COM_get_size((cflagwinname+".nc").c_str(),pane_id,&nreal,&nghost);
    if(_rank==0 && _verbose){
      std::cout << "Rocon> Burning out window " << cflagwinname 
		<< " with " << nghost << " ghost nodes." << std::endl
		<< "Rocon> Bflag stats: size = " << bflag_size 
		<< ", nghosts = " << bflag_ghosts
		<< ", stride = " << bflag_stride << ", cap = " 
		<< bflag_cap << std::endl;
    }
    int nnodes = nreal;
    int nreal_nodes  = nreal-nghost;
    int nghost_nodes = nghost;
    nreal = nghost = 0;
    
    COM_get_size((cflagwinname+".conn").c_str(),pane_id,&nreal,&nghost);
    if(nghost > 0 && !_rank){
      std::cout << "Rocon> Burning out window " << cflagwinname << " with "
		<< nghost << " ghost elements." << std::endl;
    }
    int nelem = nreal;
    int nreal_elem  = nreal - nghost;
    int nghost_elem = nghost;
    nreal = nghost = 0;
    if((nelem == 0 || nnodes == 0) && !_rank &&  _verbose){
      std::cout << "Rocon> Empty pane " << pane_id << " on window " 
		<< cflagwinname << "." << std::endl;
      continue;
    }
    // Obtain the connectivity tables
    int nmeshtypes;       // Number of connectivity tables
    std::string connNames; // Names of connectivity tables separated by space
    COM_get_connectivities(cflagwinname.c_str(), pane_id, 
			   &nmeshtypes, connNames);
    int   ndims;
    const int* dims = NULL;
    char elemType = '\0';
    int eTotal = 0;
    //      std::vector<ConnInfo> connInfo(nConn);
    if (nmeshtypes == 1 && (connNames.find(":st") != std::string::npos)) 
      { 
	// Structured mesh
	if(!_rank && _verbose){
	  std::cout << "Rocon> Processing block structured mesh." 
		    << std::endl;
	}
	COM_get_size((cflagwinname+"."+connNames).c_str(), pane_id, 
		     &nreal, &nghost);
	COM_get_array_const((cflagwinname+"."+connNames).c_str(), 
			    pane_id, &dims);
	int ndims = nreal;
	int nghost_layers = nghost;
	// nghost is wrt element layers, dims are wrt nodes
	if(!_rank && _verbose){
	  std::cout << "Rocon> Conn size (dimension,ghost) = (" 
		    << ndims << "," << nghost_layers << std::endl
		    << "Rocon> dims: ";
	  for(int i = 0;i < ndims;i++){
	    std::cout << dims[i]; 
	    if(i != (ndims-1))
	      std::cout << ",";
	  }
	  std::cout << std::endl;
	}
	std::vector<Mesh::IndexType> flat_extent(6,0);
	for(int i = 0;i < ndims;i++){
	  flat_extent[i*2]   = 1;
	  flat_extent[i*2+1] = dims[i];
	}
	for(int i = 0;i < 6;i++)
	  if(flat_extent[i] == 0)
	    flat_extent[i] = 1;
	Mesh::BSExtent<Mesh::IndexType> bsextent(flat_extent);
	bsextent.CreateUnstructuredMesh(pane_conn);
	if(!_rank && _verbose){
	  std::cout << "Rocon> Pane Connectivity has " << pane_conn.Nelem() 
		    << " elements of size " << pane_conn.Esize(1) << std::endl
		    << "Rocon> First elem: (" << pane_conn.Node(1,1) << ","
		    << pane_conn.Node(1,2) << "," << pane_conn.Node(1,3) << ","
		    << pane_conn.Node(1,4) << ")" << std::endl;
	}
      }
    else 
      {
	if(!_rank && _verbose)
	  std::cout << "Rocon> Processing unstructured mesh" << std::endl;
	Mesh::Connectivity ghost_conn;
	std::istringstream Istr(connNames);
	std::string ucname;
	int nnodes_per_element;
	while(Istr >> ucname){
	  int *conn  = NULL;
	  int stride = 0;
	  int cap    = 0;
	  int nelem  = 0;
	  nghost = 0;
	  std::string fullname(cflagwinname + "." + ucname);
	  COM_get_size(fullname.c_str(),pane_id,&nelem,&nghost);
	  COM_get_array(fullname.c_str(),pane_id,&conn,&stride,&cap);
	  if(!_rank && _verbose)
	    std::cout << "Rocon> Found connectivity " << ucname 
		      << " for pane " << pane_id 
		      << ", with size: " << nelem << ", and " 
		      << nghost << " ghosts. Stride = " 
	      	      << stride << ", and cap = " << cap << std::endl;
	  int npe = 9;
	  if(ucname.find("t3") != std::string::npos)
	    npe = 3;
	  else if(ucname.find("t6") != std::string::npos)
	    npe = 6;
	  else if(ucname.find("q4") != std::string::npos)
	    npe = 4;
	  else if(ucname.find("q8") != std::string::npos)
	    npe = 8;
	  else if(ucname.find("q9") != std::string::npos)
	    npe = 9;
	  else{
	    std::string msg("Rocon> Unknown connectivity type: "+ucname);
	    COM_assertion_msg(false,msg.c_str());
	  }
	  if(stride == 0)
	    stride = 1;
	  for(int j = 0;j < nelem;j++){
	    std::vector<Mesh::IndexType> new_element(npe,0);
	    if(stride == 1)
	      for(int k = 0;k < npe;k++)
		new_element.push_back(conn[j*npe+k]);
	    else
	      for(int k = 0;k < npe;k++)
		new_element.push_back(conn[j+k*stride]);
	    pane_conn.AddElement(new_element);
	  }
	}
      }
    pane_conn.Sync();
    // Now Burnout on this pane if it is burning
    if(burning_pane){
      for(int j = 0;j < pane_conn.Nelem();j++){
	std::vector<Mesh::IndexType>::iterator ei = pane_conn[j].begin();
	bool all_burned_out = true;
	while(ei != pane_conn[j].end() && all_burned_out)
	  if(constrained[*ei++-1] == 0)
	    all_burned_out = false;
	if(all_burned_out){
	  burning[j] = 0;
	  local_burned_out++;
	}
      }
    }
    // experimental, set stupid cflag to 1 for nonburning panes
    else {
      for(int j = 0;j < pane_conn.Nelem();j++){
	std::vector<Mesh::IndexType>::iterator ei = pane_conn[j].begin();
	while(ei != pane_conn[j].end())
	  constrained[*ei++-1] = 1;
      }
    }
  }    
  MPI_Reduce(&local_burned_out,&total_burned_out,
	     1,MPI_INTEGER,MPI_SUM,0,_comm);
  if(!_rank && _verbose)
    std::cout << "Rocon> Burned out " << total_burned_out << " elements."
	      << std::endl;
  //      if(nghost > 0)
  //	number_of_nodes -= nghost;
  //      if(!_checked){
  //	std::cout << "Rocon> Checking coords for pane " << pane_id << std::endl;
  //	CheckCoordinates(original_coordinates,number_of_nodes,1e-9);
  //      }
}

/// set bflag to 0 for every element whos nodes all
/// have cflag = 1.
void Rocon::burnout_filter( const COM::DataItem *bflag,
			    COM::DataItem *target)
{
  COM_assertion_msg( validate_object()==0, "Invalid object");
  COM_assertion_msg( bflag,   "Bflag must be present");
  COM_assertion_msg( target,  "Target attribute must be present");

  COM::Window *targwin = const_cast<COM::Window*>(target->window());
  COM::Window *bflagwin  = const_cast<COM::Window*>(bflag->window());

  std::string targwinname(targwin->name());
  std::string bflagwinname(bflagwin->name());
  
  std::string targname(target->name());
  std::string bflagname(bflag->name());

  unsigned int number_of_panes = 0;
  std::vector<int> pane_ids;
  COM_get_panes(targwinname.c_str(),pane_ids);
  number_of_panes = pane_ids.size();
  int local_nfiltered  = 0;
  int global_nfiltered = 0;
  std::vector<int>::iterator pii = pane_ids.begin();
  while(pii != pane_ids.end()){
    int    *burning               = NULL;
    double *targ                  = NULL;
    int targ_size                 = 0;
    int bflag_size                = 0;
    int bflag_ghosts              = 0;
    int targ_ghosts               = 0;
    int pane_id = *pii++;
    COM_get_array((bflagwinname+"."+bflagname).c_str(),pane_id,&burning);
    COM_get_array((targwinname+"."+targname).c_str(),pane_id,&targ);
    if(burning && targ){
      COM_get_size((bflagwinname+"."+bflagname).c_str(),pane_id,
		   &bflag_size,&bflag_ghosts);
      COM_get_size((targwinname+"."+targname).c_str(),pane_id,
		   &targ_size,&targ_ghosts);
      assert(targ_size == bflag_size);
      if(_rank==0 && _verbose){
	std::cout << "Rocon> Burning out " << targname << " on " 
		  << targwinname << std::endl
		  << "Rocon> Target stats: size = " << targ_size 
		  << ", nghosts = " << targ_ghosts << std::endl;
      }
      for(int i = 0;i < targ_size;i++)
	if(burning[i] == 0){
	  local_nfiltered++;
	  targ[i] = 0.0;
	}
    }
  }
  MPI_Reduce(&local_nfiltered,&global_nfiltered,
	     1,MPI_INTEGER,MPI_SUM,0,_comm);
  if(!_rank && _verbose)
    std::cout << "Rocon> Burnout filtered " << global_nfiltered << " elements."
	      << std::endl;
}

/// Displace the points in pmesh by disp and calcuate the intersections
/// (if any) with the constraint mesh set in "initialize".  The int IO 
/// attribute "constr" indicates (by a 1 or 0) which nodes are already 
/// on the constraint surface on input.  On output, "constr" will be
/// updated to include the newly constrained nodes as well. On output,
/// pos contains the intersection position of any nodes which collide
/// with the constraint surface.
void Rocon::find_intersections( const COM::DataItem *pmesh, 
				const COM::DataItem *disp,
				COM::DataItem *pos,
				COM::DataItem *constr)
{
  COM_assertion_msg( validate_object()==0, "Invalid object");
  COM_assertion_msg( pmesh, "Mesh must be present");
  COM_assertion_msg( disp,  "Displacements must be present");
  COM_assertion_msg( pos,   "Positions att must be present");

  COM::Window *meshwin = const_cast<COM::Window*>(pmesh->window());
  COM::Window *dispwin = const_cast<COM::Window*>(disp->window());
  COM::Window *poswin  = const_cast<COM::Window*>(pos->window());
  COM::Window *conwin  = const_cast<COM::Window*>(constr->window());

  std::string meshwinname(meshwin->name());
  std::string dispwinname(dispwin->name());
  std::string poswinname(poswin->name());
  std::string conwinname(conwin->name());
  
  std::string dispname(disp->name());
  std::string posname(pos->name());
  std::string conname(constr->name());

  unsigned int number_of_panes = 0;
  std::vector<int> pane_ids;
  COM_get_panes(meshwinname.c_str(),pane_ids);
  number_of_panes = pane_ids.size();

  std::vector<int> paneids;
  COM_get_panes(dispwinname.c_str(),paneids);
  COM_assertion_msg(paneids.size() == number_of_panes,
		    "Number of panes must match.");
  COM_get_panes(poswinname.c_str(),paneids);
  COM_assertion_msg(paneids.size() == number_of_panes,
		    "Number of panes must match.");
  COM_get_panes(conwinname.c_str(),paneids);
  COM_assertion_msg(paneids.size() == number_of_panes,
		    "Number of panes must match.");
  
  std::vector<int>::iterator pii = pane_ids.begin();
  int local_number_constrained        = 0;
  int local_number_checked            = 0;
  int global_number_constrained       = 0;
  int global_number_checked           = 0;
  while(pii != pane_ids.end()){
    double *original_coordinates  = NULL;
    double *displacements         = NULL;
    double *positions             = NULL;
    int    *constrained           = NULL;
    int number_of_nodes           = 0;
    int nghost                    = 0;
    int pane_id = *pii++;
    COM_get_array((meshwinname+".nc").c_str(),pane_id,&original_coordinates);
    COM_get_size((meshwinname+".nc"),pane_id,&number_of_nodes,&nghost);
    if(nghost > 0)
      number_of_nodes -= nghost;
    if(!_checked){
      //      std::cout << "Rocon> Checking coords for pane " << pane_id << std::endl;
      CheckCoordinates(original_coordinates,number_of_nodes,1e-9);
    }
    COM_get_array((dispwinname+"."+dispname).c_str(),pane_id,&displacements);
    COM_get_array((poswinname+"."+posname).c_str(),pane_id,&positions);
    COM_get_array((conwinname+"."+conname).c_str(),pane_id,&constrained);
    if(constrained == NULL && !_rank && _verbose) {
      std::cout << "Rocon> Skipping pane " << pane_id << " for no cflag." << std::endl;
    }
    if((constrained != NULL) && (number_of_nodes > 0)){
      // In an attempt to fix some problems, a check is added here to see if this
      // is a burning pane.   If not, then it'll set all nodes to "constrained" on
      // that pane.   Hopefully this will help sticky elements finally burn out
      // and keep nonburning panes from getting jerked around anomalously.
      int *burning_flag = NULL;
      int *bcflag       = NULL;
      COM_get_array((conwinname+".bflag").c_str(),pane_id,&burning_flag);
      //      COM_get_array((conwinname+".bcflag").c_str(),pane_id,&bcflag);
      //      bool pane_is_burning = true;
      //      if((bcflag != NULL) && (burning_flag != NULL))
      //	if(*bcflag == 1)
      //	  pane_is_burning = true;
      bool pane_is_burning = (burning_flag != NULL);
      if(pane_is_burning)
	local_number_checked += number_of_nodes;
      for(int i = 0;i < number_of_nodes;i++){
	// could check here if they are already constrained going in, 
	// then just set disp to 0 else perform the check below. disp 
	// to 0 and/or positions to current pos.
	//	if(constrained[i] == 0){
	if(pane_is_burning){
	  double point[3];
	  double disp[3];
	  point[0] = original_coordinates[i*3]   - displacements[i*3];
	  point[1] = original_coordinates[i*3+1] - displacements[i*3+1];
	  point[2] = original_coordinates[i*3+2] - displacements[i*3+2];
	  disp[0] = 2.0*displacements[i*3];
	  disp[1] = 2.0*displacements[i*3+1];
	  disp[2] = 2.0*displacements[i*3+2];
		  //point[0] = original_coordinates[i*3]   - displacements[i*3];
		  // point[1] = original_coordinates[i*3+1] - displacements[i*3+1];
		  //point[2] = original_coordinates[i*3+2] - displacements[i*3+2];
		  //disp[0] = 2.0*displacements[i*3];
		  //disp[1] = 2.0*displacements[i*3+1];
		  //disp[2] = 2.0*displacements[i*3+2];
	  //	if(this->_mbs->intersection(&original_coordinates[i*3],
	  //				    &displacements[i*3],&positions[i*3])){
	  if(this->_mbs->intersection(point,disp,&positions[i*3])){
	    local_number_constrained++;
	    constrained[i] = 1;
	  }
	  else{
	    positions[i*3]   = original_coordinates[i*3]   + displacements[i*3];
	    positions[i*3+1] = original_coordinates[i*3+1] + displacements[i*3+1];
	    positions[i*3+2] = original_coordinates[i*3+2] + displacements[i*3+2];
	  }
	    
	}
	else{
	  constrained[i] = 1;
	  // Experimental - 0 out displacements
	  //	  displacements[i*3]   = 0.0;
	  //	  displacements[i*3+1] = 0.0;
	  //  displacements[i*3+2] = 0.0;
	}
	//	}
	//	else{
	//	  displacements[i*3]   = 0.0;
	//	  displacements[i*3+1] = 0.0;
	//	  displacements[i*3+2] = 0.0;
	//	}
      }
    }
  }
  this->_comm = meshwin->get_communicator();
  if ( COMMPI_Initialized())  this->_rank = COMMPI_Comm_rank(this->_comm);
  else this->_rank = 0;
  MPI_Reduce(&local_number_constrained,&global_number_constrained,
	     1,MPI_INTEGER,MPI_SUM,0,_comm);
  MPI_Reduce(&local_number_checked,&global_number_checked,
	     1,MPI_INTEGER,MPI_SUM,0,_comm);
  if(!_rank){
    std::cout << "Rocon> Checked " << global_number_checked 
	      << " nodes." << std::endl;
    std::cout << "Rocon> Found " << global_number_constrained 
	      << " intersections." << std::endl;
  }
  //  _checked = true;
}

/// Constrain the displacements, disp, to the appropriate positions, pos.
void Rocon::constrain_displacements( const COM::DataItem *pmesh, 
				     COM::DataItem *disp,
				     const COM::DataItem *pos, 
				     const COM::DataItem *constr)
{
  COM_assertion_msg( validate_object()==0, "Invalid object");
  COM_assertion_msg( pmesh, "Mesh must be present");
  COM_assertion_msg( disp,  "Displacements must be present");
  COM_assertion_msg( pos,   "Positions must be present");
  
  COM::Window *meshwin = const_cast<COM::Window*>(pmesh->window());
  COM::Window *dispwin = const_cast<COM::Window*>(disp->window());
  COM::Window *poswin  = const_cast<COM::Window*>(pos->window());
  COM::Window *conwin  = const_cast<COM::Window*>(constr->window());

  std::string meshwinname(meshwin->name());
  std::string dispwinname(dispwin->name());
  std::string poswinname(poswin->name());
  std::string conwinname(conwin->name());
  
  std::string dispname(disp->name());
  std::string posname(pos->name());
  std::string conname(constr->name());

  unsigned int number_of_panes = 0;
  std::vector<int> pane_ids;
  COM_get_panes(meshwinname.c_str(),pane_ids);
  number_of_panes = pane_ids.size();

  std::vector<int> paneids;
  COM_get_panes(dispwinname.c_str(),paneids);
  COM_assertion_msg(paneids.size() == number_of_panes,
		    "Number of panes must match.");
  COM_get_panes(poswinname.c_str(),paneids);
  COM_assertion_msg(paneids.size() == number_of_panes,
		    "Number of panes must match.");
  COM_get_panes(conwinname.c_str(),paneids);
  COM_assertion_msg(paneids.size() == number_of_panes,
		    "Number of panes must match.");
  
  std::vector<int>::iterator pii = pane_ids.begin();
  int local_number_checked = 0;
  int global_number_checked = 0;
  int local_number_constrained = 0;
  int global_number_constrained = 0;
  while(pii != pane_ids.end()){
    double *original_coordinates  = NULL;
    double *displacements         = NULL;
    double *positions             = NULL;
    int    *constrained           = NULL;
    int number_of_nodes           = 0;
    int nghost                    = 0;
    int pane_id = *pii++;
    COM_get_array((meshwinname+".nc").c_str(),pane_id,&original_coordinates);
    COM_get_size((meshwinname+".nc"),pane_id,&number_of_nodes,&nghost);
    if(nghost > 0)
      number_of_nodes -= nghost;
    COM_get_array((dispwinname+"."+dispname).c_str(),pane_id,&displacements);
    COM_get_array((poswinname +"."+posname).c_str(),pane_id,&positions);
    COM_get_array((conwinname +"."+conname).c_str(),pane_id,&constrained);
    if((constrained != NULL) && (number_of_nodes > 0)){
      int *burning_flag = NULL;
      int *bcflag  = NULL;
      COM_get_array((conwinname+".bflag").c_str(),pane_id,&burning_flag);
      //      COM_get_array((conwinname+".bcflag").c_str(),pane_id,&bcflag);
      bool pane_is_burning = (burning_flag != NULL);
      //      if(bcflag && burning_flag)
      //	if(*bcflag == 1)
      //	  pane_is_burning = true;
      if(pane_is_burning){
	local_number_checked += number_of_nodes;
	for(int i = 0;i < number_of_nodes;i++){
	  if(constrained[i] == 1){
	    local_number_constrained++;
	    displacements[i*3]   = positions[i*3]   - original_coordinates[i*3];
	    displacements[i*3+1] = positions[i*3+1] - original_coordinates[i*3+1];
	    displacements[i*3+2] = positions[i*3+2] - original_coordinates[i*3+2];
	    // experimental, just stop them
	    //	    displacements[i*3]   = 0.0;
	    //	    displacements[i*3+1] = 0.0;
	    //	    displacements[i*3+2] = 0.0;
	    if(std::fabs(displacements[i*3])   <= this->_TOL) 
	      displacements[i*3]   = 0;
	    if(std::fabs(displacements[i*3+1]) <= this->_TOL) 
	      displacements[i*3+1] = 0;
	    if(std::fabs(displacements[i*3+2]) <= this->_TOL)
	      displacements[i*3+2] = 0;
	  }
	}
      }
    }
  }
  this->_comm = meshwin->get_communicator();
  if ( COMMPI_Initialized())  this->_rank = COMMPI_Comm_rank(this->_comm);
  else this->_rank = 0;
  MPI_Reduce(&local_number_constrained,&global_number_constrained,
	     1,MPI_INTEGER,MPI_SUM,0,_comm);
  if(!_rank){
    std::cout << "Rocon> Constrained " << global_number_constrained 
	      << " nodes." << std::endl;
  }


}



void Rocon::load( const std::string &mname) {
  Rocon *rp   = new Rocon();
  COM_new_window( mname.c_str());
  rp->WindowName() = mname;
  rp->SetVerbosity(1);
  std::string glb=mname+".global";
  COM_new_dataitem( glb.c_str(), 'w', COM_VOID, 1, "");
  COM_set_object( glb.c_str(), 0, rp);
  
  COM_Type types[5];

  types[0] = COM_RAWDATA;
  types[2] = COM_INT;
  types[1] = types[3] = types[4] = COM_METADATA;
  COM_set_member_function( (mname+".initialize").c_str(), 
  			   (Member_func_ptr)(&Rocon::initialize), 
  			   glb.c_str(), "bii", types);

  types[1] = COM_STRING;
  COM_set_member_function( (mname+".init_from_file").c_str(),
			   (Member_func_ptr)(&Rocon::init_from_file),
			   glb.c_str(),"bii", types);
  types[1] = COM_METADATA;
  types[2] = COM_METADATA;
  COM_set_member_function( (mname+".find_intersections").c_str(), 
  			   (Member_func_ptr)(&Rocon::find_intersections), 
  			   glb.c_str(), "biiob", types);
  COM_set_member_function( (mname+".constrain_displacements").c_str(), 
  			   (Member_func_ptr)(&Rocon::constrain_displacements), 
  			   glb.c_str(), "bibii", types);

  COM_set_member_function( (mname+".burnout").c_str(), 
  			   (Member_func_ptr)(&Rocon::burnout), 
  			   glb.c_str(), "biib", types);
  
  COM_set_member_function( (mname+".burnout_filter").c_str(), 
  			   (Member_func_ptr)(&Rocon::burnout_filter), 
  			   glb.c_str(), "bib", types);
  
    
  COM_window_init_done( mname.c_str());
}

void Rocon::unload( const std::string &mname) {
  Rocon *rp;
  std::string glb=mname+".global";

  COM_get_object( glb.c_str(), 0, &rp);
  
  COM_assertion_msg( rp->validate_object()==0, "Invalid object");

  delete rp;
  COM_delete_window( mname.c_str());
}

// C/C++ binding
extern "C" void Rocon_load_module( const char *name)
{ Rocon::load( std::string(name)); }
extern "C" void Rocon_unload_module( const char *name) 
{ Rocon::unload( std::string(name)); }

// Fortran binding
extern "C" void rocon_load_module( const char *name, long int length)
{ Rocon::load( std::string(name, length)); }
extern "C" void rocon_unload_module( const char *name, long int length) 
{ Rocon::unload( std::string(name, length)); }

extern "C" void ROCON_LOAD_MODULE( const char *name, long int length)
{ Rocon::load( std::string(name, length)); }
extern "C" void ROCON_UNLOAD_MODULE( const char *name, long int length) 
{ Rocon::unload( std::string(name, length)); }

extern "C" void rocon_load_module_( const char *name, long int length)
{ Rocon::load( std::string(name, length)); }
extern "C" void rocon_unload_module_( const char *name, long int length) 
{ Rocon::unload( std::string(name, length)); }

extern "C" void ROCON_LOAD_MODULE_( const char *name, long int length)
{ Rocon::load( std::string(name, length)); }
extern "C" void ROCON_UNLOAD_MODULE_( const char *name, long int length) 
{ Rocon::unload( std::string(name, length)); }
