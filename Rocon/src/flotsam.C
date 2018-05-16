#include <iostream>
#include <iomanip>
#include <limits>

#include "Global.H"
#include "PMesh.H"
#include "BSMesh.H"
#include "Profiler.H"

#include "FloGrid.H"
#include "flotsam_cl.H"


// Mesh::FindPointInMesh: Error: Closest approach: 0.0254483 at node 146.
// Mesh::FindPointInMesh: Element(113) = (129,130,138,137,361,362,370,369)
//  Point: 0.595836  0.111593  0.0254
//  Nodes: 
//  0: 0.58992  0.107061  0
//  1: 0.5946  0.103897  0
//  2: 0.5946  0.107263  0
//  3: 0.58992  0.110554  0
//
//  4: 0.58992  0.107061  0.0508
//  5: 0.5946  0.103897  0.0508
//  6: 0.5946  0.107263  0.0508
//  7: 0.58992  0.110554  0.0508
// Mesh::FindPointInMesh: Error: Closest approach: 0.0255076 at node 139.
// Mesh::FindPointInMesh: Element(107) = (122,123,131,130,354,355,363,362)
//  Point: 0.60162  0.103916  0.0254
//  Nodes: 
//   0: 0.5946  0.100531  0
//   1: 0.59928  0.0974934  0
//   2: 0.59928  0.100733  0
//   3: 0.5946  0.103897  0
//
//   4: 0.5946  0.100531  0.0508
//   5: 0.59928  0.0974934  0.0508
//   6: 0.59928  0.100733  0.0508
//   7: 0.5946  0.103897  0.0508
// LUDcmp::error: Singular matrix
// Mesh::FindPointInMesh: Error: Closest approach: 0.0256659 at node 238.
// Mesh::FindPointInMesh: Element(13) = (14,15,23,22,246,247,255,254)
//   Point: 0.696963  -0.0204658  0.0762
//   Nodes:
//   0: 0.699524  -0.0237386  0
//   1: 0.694296  -0.0245241  0
//   2: 0.694257  -0.0253172  0
//   3: 0.699457  -0.0244292  0
//   4: 0.699524  -0.0237386  0.0508
//   5: 0.694296  -0.0245241  0.0508
//   6: 0.694257  -0.0253172  0.0508
//   7: 0.699457  -0.0244292  0.0508
//   Mesh::FindPointInMesh: Error: Closest approach: 0.0255667 at node 235.
//   Mesh::FindPointInMesh: Element(2) = (2,3,11,10,234,235,243,242)
//     Point: 0.712735  -0.0197399  0.0762
//   Nodes:
//   0: 0.720621  -0.020316  0
//   1: 0.715364  -0.020999  0
//   2: 0.715261  -0.0214581  0
//   3: 0.72052  -0.0207251  0
//   4: 0.720621  -0.020316  0.0508
//   5: 0.715364  -0.020999  0.0508
//   6: 0.715261  -0.0214581  0.0508
//   7: 0.72052  -0.0207251  0.0508
//   LUDcmp::error: Singular matrix
//   Mesh::FindPointInMesh: Error: Closest approach: 0.0255693 at node 236.
//   Mesh::FindPointInMesh: Element(19) = (21,22,30,29,253,254,262,261)
//     Point: 0.707478  -0.0203717  0.0762
//   Nodes:
//   0: 0.704675  -0.0235657  0
//   1: 0.699457  -0.0244292  0
//   2: 0.69939  -0.0251198  0
//   3: 0.704589  -0.024166  0
//   4: 0.704675  -0.0235657  0.0508
//   5: 0.699457  -0.0244292  0.0508
//   6: 0.69939  -0.0251198  0.0508
//   7: 0.704589  -0.024166  0.0508                                                                                                               
void test2(void)
{
  Mesh::BSExtent<Mesh::IndexType> ghost_extent;
  Mesh::BSExtent<Mesh::IndexType> real_extent;
  std::vector<Mesh::IndexType> flat_ghost(6,0);
  std::vector<Mesh::IndexType> flat_real(6,0);
  flat_ghost[0] = flat_ghost[2] = flat_ghost[4] = 1;
  flat_ghost[1] = flat_ghost[3] = flat_ghost[5] = 4;
  flat_real[0] = flat_real[2] = flat_real[4] = 2;
  flat_real[1] = flat_real[3] = flat_real[5] = 3;
  ghost_extent.Init(flat_ghost);
  real_extent.Init(flat_real);
  std::vector<Mesh::IndexType> flat_real_indices;
  ghost_extent.GetFlatIndices(real_extent,flat_real_indices);
  std::vector<Mesh::IndexType>::iterator fri = flat_real_indices.begin();
  while(fri != flat_real_indices.end()){
    std::cout << *fri++ << " ";
  }
  std::cout << std::endl;
}

void test(void)
{
  Mesh::Connectivity conn;
  Mesh::NodalCoordinates nc(24);
  GeoPrim::CPoint p;
  nc.init_node(1,p.init(.58992,.107061,0));
  nc.init_node(2,p.init(.5946,0.103897,0));
  nc.init_node(3,p.init(.5946,0.107263,0));
  nc.init_node(4,p.init(.58992,0.110554,0));
  nc.init_node(5,p.init(.58992,0.107061,.0508));
  nc.init_node(6,p.init(.5946,.103897,0.0508));
  nc.init_node(7,p.init(.5946,.107263,0.0508));
  nc.init_node(8,p.init(.58992,0.110554,0.0508));
  nc.init_node(9,p.init(0.5946,0.100531,0));
  nc.init_node(10,p.init(0.59928,0.0974934,0));
  nc.init_node(11,p.init(0.59928,0.100733 ,0));
  nc.init_node(12,p.init(0.5946,0.103897,0));
  nc.init_node(13,p.init(0.5946,0.100531,0.0508));
  nc.init_node(14,p.init(0.59928,0.0974934,0.0508));
  nc.init_node(15,p.init(0.59928,0.100733,0.0508));
  nc.init_node(16,p.init(0.5946,0.103897,0.0508));
  nc.init_node(17,p.init(0.699524,-0.0237386,0));
  nc.init_node(18,p.init(0.694296,-0.0245241,0));
  nc.init_node(19,p.init(0.694257,-0.0253172,0));
  nc.init_node(20,p.init(0.699457,-0.0244292,0));
  nc.init_node(21,p.init(0.699524,-0.0237386,0.0508));
  nc.init_node(22,p.init(0.694296,-0.0245241,0.0508));
  nc.init_node(23,p.init(0.694257,-0.0253172,0.0508));
  nc.init_node(24,p.init(0.699457,-0.0244292,0.0508));
  conn.resize(3);
  conn[0].push_back(1);
  conn[0].push_back(2);
  conn[0].push_back(3);
  conn[0].push_back(4);
  conn[0].push_back(5);
  conn[0].push_back(6);
  conn[0].push_back(7);
  conn[0].push_back(8);
  conn[1].push_back(9);
  conn[1].push_back(10);
  conn[1].push_back(11);
  conn[1].push_back(12);
  conn[1].push_back(13);
  conn[1].push_back(14);
  conn[1].push_back(15);
  conn[1].push_back(16);
  conn[2].push_back(17);
  conn[2].push_back(18);
  conn[2].push_back(19);
  conn[2].push_back(20);
  conn[2].push_back(21);
  conn[2].push_back(22);
  conn[2].push_back(23);
  conn[2].push_back(24);
  conn.Sync();
//  Point: 0.595836  0.111593  0.0254
//   Point: 0.696963  -0.0204658  0.0762
  std::vector<Mesh::IndexType> els(3);
  els[0] = 1;
  els[1] = 2;
  els[2] = 3;
  GeoPrim::CVector natc;
  unsigned int cell_id = 0;
  std::vector<GeoPrim::CPoint> test_points(3);
  test_points[0].init(.595836,.111593,0.0254);
  test_points[1].init( 0.60162,0.103916,0.0254 );
  test_points[2].init(0.696963 ,-0.0204658 ,0.0762);
  std::vector<GeoPrim::CPoint>::iterator tpi = test_points.begin();
  while(tpi != test_points.end())
    cell_id = Mesh::FindPointInCells(*tpi++,nc,conn,els,natc);
    
}

void writeVtkData( Mesh::NodalCoordinates &nc,Mesh::Connectivity &con,
		   const std::string &filename,std::vector<double> &soln,
		   std::vector<Mesh::IndexType> indices,
		   unsigned int stride)
{
  
  //  std::cout << "Writing VTK data...";
  //  std::cout.flush( );
  
  std::ofstream ofs;
  ofs.open( filename.c_str());
  
  if( !ofs.is_open( ) )
    {
      std::cerr << "Cannot write VTK file: " << filename << std::endl;
      std::cerr << "File: " << __FILE__ 		 << std::endl;
      std::cerr << "Line: " << __LINE__ 		 << std::endl;
      assert( false );
    }
  
  unsigned int nnodes = nc.NNodes();
  /* STEP 1: Write the header */
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << "Global Mesh Output" << std::endl;
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";
  ofs << "POINTS " << nnodes << " double\n";
  
  /* STEP 2: Write the nodes */
  for( int i=1; i <= nnodes; ++i )
    {
      ofs << GeoPrim::C3Point(nc[i]) << std::endl;
    }
  
  /* STEP 3: Write the max element size */
  int maxnodes = 0; /* the maximum number of nodes per element */
  unsigned int nelem = con.size();
  con.SyncSizes();
  for( int i=0; i < nelem; ++i )
    {
      unsigned int esize = con[i].size();
      if(maxnodes < esize)
	maxnodes = esize;
    }
  
  /* STEP 4: Write the element connectivity */
  ofs << "CELLS " << nelem << " " << nelem*( maxnodes+1) << "\n";
  for( int i=0; i < nelem; ++i  )
    {
      ofs <<  con.Esize(i+1) << " ";
      for( int j=0; j < con.Esize(i+1); ++j )
	ofs <<  con[ i ][ j ]-1 << " ";
      ofs << std::endl;
    }
  
  /* STEP 5: Write the cell types */
  ofs << "CELL_TYPES " << nelem << std::endl;
  for( int i=0; i < nelem; ++i )
    {
      
      if( con.Esize(i+1) == 8 )
	ofs << "12\n";
      else if( con.Esize(i+1) == 4 )
	ofs << "10\n";
      else {
	std::cerr << "Undefined element type!\n";
	std::cerr << "File: " << __FILE__ << std::endl;
	std::cerr << "Line: " << __LINE__ << std::endl;
      }
      
    }
  
  /* STEP 6: Write the solution data attached to the points */
  bool wroteHeader = false;
  std::vector<std::string> Name(5);
  Name[0] = "rho";
  Name[1] = "rho-u";
  Name[2] = "rho-v";
  Name[3] = "rho-w";
  Name[4] = "rho-E";
  if(!soln.empty()){
    for( int i=0; i < 5; ++i )
      {   
	if( !wroteHeader )
	  {
	    ofs << "CELL_DATA " << nelem << std::endl;
	    wroteHeader = true;
	  }
	ofs << "SCALARS " << Name[i] << " double" << std::endl;
	ofs << "LOOKUP_TABLE default\n";
	for( int j=0; j < nelem; ++j )
	  {
	    
	    ofs << soln[ i*stride + indices[j] - 1] << " ";
	  }
	ofs << std::endl;
      }
  }
  ofs.close( );
  //  std::cout << "[DONE]\n";
}

int main(int argc,char *argv[])
{
  //  test2();
  //  return(0);
  std::ostream* StdOut = NULL; // program output
  std::ostream* ErrOut = NULL; // program errors
  std::ostream* LogOut = NULL; // log for each proc (if debugging on)
  std::ofstream LogFile;
  IRAD::Global::GlobalObj<std::string,Mesh::IndexType,IRAD::Profiler::ProfilerObj> global;
  IRAD::Comm::CommunicatorObject comm(&argc,&argv);
  unsigned int rank  = comm.Rank();
  unsigned int nproc = comm.Size();
  //  std::cout << "rank/nproc: " << rank << "/" << nproc << std::endl;
  if(rank == 0){
    StdOut = &std::cout;
    ErrOut = &std::cerr;
  }
  global.Init("flotsam",rank);
  global.SetDebugLevel(0);
  global.FunctionEntry("main");
  
  FlotsamComLine comline((const char **)argv);
  comline.Initialize();
  int clerr = comline.ProcessOptions();
  std::string sverb    =  comline.GetOption("verb");
  bool debug           = !comline.GetOption("debug").empty();
  bool array_checking  = !comline.GetOption("checking").empty();

  IRAD::Comm::DataTypes IndexDT = (sizeof(Mesh::IndexType) == sizeof(size_t) ?
			     IRAD::Comm::DTSIZET : IRAD::Comm::DTUINT);
  if(!comline.GetOption("help").empty()){
    if(StdOut) 
      *StdOut << comline.LongUsage() << std::endl;
    comm.SetExit(1);
  }
  if(comm.Check())
    return(0);
  
  if(clerr){
    if(ErrOut) 
      *ErrOut << comline.ErrorReport() << std::endl
	      << std::endl << comline.ShortUsage() << std::endl;
    comm.SetExit(1);
  }
  if(comm.Check())
    return(1);
  int verblevel = 0;
  if(!sverb.empty()){
    verblevel = 1;
    if(sverb != ".true."){
      std::istringstream Istr(sverb);
      Istr >> verblevel;
    }
  }
  if(nproc > 1){
    if(rank == 0 && verblevel > 0)
      std::cout << "Flotsam running on " << nproc << " processors." << std::endl; 
    std::ostringstream debugfilename;
    debugfilename << comline.ProgramName() << ".output." << rank;
    LogFile.open(debugfilename.str().c_str());
    global.SetDebugStream(LogFile);
    LogOut = &LogFile;
    ErrOut = &LogFile;
    if(verblevel > 1 && rank > 0)
      StdOut = &LogFile;
  }
  else{
    global.SetDebugStream(std::cerr);
    LogOut=&std::cout;
    ErrOut=&std::cout;
  }
  if(verblevel > 1 || debug){
    global.SetDebugLevel(1);
  }
  if(debug && rank > 0)
    StdOut = LogOut;
  std::vector<std::string> infiles = comline.GetArgs();
  if(infiles.size()==0) {
    if(ErrOut) 
      *ErrOut << "flotsam::Error: No input specified." << std::endl
	      << std::endl << comline.ShortUsage() << std::endl;
    comm.SetExit(1);
  }
  if(comm.Check())
    return(1);

  // XXXXXXXXXXXXXXXXXXXXXXXXXX
  // Do not modify above here.  
  // XXXXXXXXXXXXXXXXXXXXXXXXXX

  std::string casename(infiles[0]);
  std::string timestamp(infiles[1]);
  std::string targname(infiles[2]);

  FloGrid sourcegrid(casename);
  FloGrid targetgrid(targname);

  unsigned int total_source_cells = 0;
  unsigned int total_source_nodes = 0;
  unsigned int total_target_cells = 0;
  unsigned int total_target_nodes = 0;
  unsigned int NSourceRegions = 0;
  unsigned int NTargetRegions = 0;
  unsigned int totalcells = 0;

  // This function reads the entire Rocflo grid from the file
  // casename.grda_0.00000E-00
  // The grid blocks are kept separated as they are in the input
  std::vector<std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> > > target_source_cells;
  std::vector<std::vector<Mesh::IndexType> > source_cell_ids;
  std::vector<std::vector<Mesh::IndexType> > target_cell_ids;

  std::vector<Mesh::Connectivity> source_mesh;
  std::vector<Mesh::Connectivity> target_mesh;

  bool do_search = true;
  if(do_search) {
    global.FunctionEntry("ReadFloGrid");
    sourcegrid.ReadAllBlocks();
    global.FunctionExit("ReadFloGrid");
    
    
    // The number of regions is the number of blocks in the input grid
    NSourceRegions = sourcegrid.Blocks().size();
    
    std::vector<std::vector<double> > source_cell_centers(NSourceRegions);

    // The extents are just an abstraction to represent indexing into a logically
    // rectangular dataset ijk.
    std::vector<Mesh::BSExtent<Mesh::IndexType> > source_extents(NSourceRegions);
    
    // The connectivities are the unstructured-like representation of the input
    // block structured mesh.
    source_mesh.resize(NSourceRegions);
    
    // The boxes are: boxes[region_index][0] = Box around entire region
    //                boxes[region_index][1] =  "    "    smallest element in region
    //                boxes[region_index][2] =  "    "    largest  element in region
    GeoPrim::CBox sourcebox;
    std::vector<std::vector<GeoPrim::CBox> > source_boxes(NSourceRegions);
    
    

    // This section loops through and populates the unstructured representations and 
    // block bounding boxes.
    std::vector<FloGridBlock>::iterator source_block_iterator = sourcegrid.Blocks().begin();
    std::vector<Mesh::BSExtent<Mesh::IndexType> >::iterator source_extent_iterator = source_extents.begin();
    std::vector<Mesh::Connectivity>::iterator source_mesh_iterator = source_mesh.begin();
    std::vector<std::vector<GeoPrim::CBox> >::iterator source_boxes_iterator = source_boxes.begin();
    while( source_block_iterator != sourcegrid.Blocks().end())
      {
	// Convert the block coordinates from
	// block format to sequential format, i.e.
	// [--x--,--y--,--z--] to [xyzxyzxyz]
	source_block_iterator->MeshLayout();
	
	// Block extents go from 1 to (size+1). Size indicates
	// number of *elements* in a given direction while
	// the number of nodes in that direction is size+1.
	std::vector<Mesh::IndexType> source_block_extent(6,0);
	unsigned int isize = source_block_iterator->isize();
	unsigned int jsize = source_block_iterator->jsize();
	unsigned int ksize = source_block_iterator->ksize();
	source_block_extent[0] = source_block_extent[2] = source_block_extent[4] = 1;
	source_block_extent[1] = isize+1;
	source_block_extent[3] = jsize+1;
	source_block_extent[5] = ksize+1;
	source_extent_iterator->Init(source_block_extent);
	
	// Use the extent object's canned routine for converting
	// to an unstructured connectivity array
	global.FunctionEntry("CreateUnsMesh");
	source_extent_iterator->CreateUnstructuredMesh(*source_mesh_iterator);
	source_mesh_iterator->Sync();
	source_mesh_iterator->SyncSizes();
	global.FunctionExit("CreateUnsMesh");
	
	
	// Create a lightweight NodalCoordinates object around 
	// the existing block coordinates in place.  This will
	// be used in the canned unstructured routines.
	Mesh::NodalCoordinates nc;
	unsigned int nnodes = source_block_iterator->NNodes();
	nc.init(nnodes,&((source_block_iterator->Coords())[0]));
	unsigned int source_region_index = source_block_iterator - sourcegrid.Blocks().begin();
	Mesh::GetMeshCentroids(nc,*source_mesh_iterator,source_cell_centers[source_region_index]);
	global.FunctionEntry("Orient elements");
	Mesh::IndexType invertedcount = 0;
	unsigned int number_of_elements = isize*jsize*ksize;
	for(Mesh::IndexType element_being_processed = 1;
	    element_being_processed <= number_of_elements;
	    element_being_processed++)
	  {
	    Mesh::IndexType index = element_being_processed - 1;
	    Mesh::IndexType size_of_element = (*source_mesh_iterator)[index].size();
	    Mesh::GenericElement ge(size_of_element);
	    if(ge.Inverted((*source_mesh_iterator)[index],nc)){
	      invertedcount++;
	      ge.ReOrient((*source_mesh_iterator)[index]);
	    }
	  }
	global.FunctionExit("Orient elements");
	//      if(StdOut && verblevel) 
	//	*StdOut << "Number of elements reoriented: " << invertedcount 
	//		<< std::endl;
	
	global.FunctionEntry("BoxProcessing");
	
	// Get the characteristic boxes on the new unstructured mesh
	source_boxes_iterator->resize(3);
	Mesh::GetMeshBoxes(nc,*source_mesh_iterator,(*source_boxes_iterator)[0],(*source_boxes_iterator)[1],
			   (*source_boxes_iterator)[2]);
	sourcebox.merge((*source_boxes_iterator)[0]); // for debugging - ensure sourcebox = region box
	global.FunctionExit("BoxProcessing");
	
	// track the total number of nodes and cells
	total_source_nodes += nnodes;
	total_source_cells += isize*jsize*ksize;
	
	source_boxes_iterator++;
	source_mesh_iterator++;
	source_extent_iterator++;
	source_block_iterator++;
      }
    
    if(StdOut && verblevel > 0)
      *StdOut << "Source grid has " << NSourceRegions << " regions, " 
	      << total_source_nodes << " nodes, and " << total_source_cells
	      << " cells." << std::endl
	      << "With bounding box: " << sourcebox << std::endl;
    
    // 
    // Repeat everything above on the target mesh
    //
    //    FloGrid targetgrid(targname);
    GeoPrim::CBox targetbox;
    if(!rank)
      NTargetRegions = targetgrid.BlockCount();
    comm.BroadCast(NTargetRegions,0);
    if(debug)
      std::cout << "ntarget regions = " << NTargetRegions << std::endl;
    int ntargregpproc = NTargetRegions/nproc;
    int nleftovers = NTargetRegions%nproc;
    int first = rank * ntargregpproc + 1;
    int ntargreg = ntargregpproc;
    if(rank < nleftovers){
      ntargreg++;
      first+=rank;
    }
    else
      first += nleftovers;
    if(debug){
      std::cout << "First/N: " << first << "/" << ntargreg << std::endl;
    }
    global.FunctionEntry("ReadFloGrid");
    //    targetgrid.ReadAllBlocks();
    targetgrid.ReadNBlocks(first,ntargreg);
    global.FunctionExit("ReadFloGrid");    
    unsigned int TargetGridSize = NTargetRegions;
    NTargetRegions = ntargreg;
    //    NTargetRegions = targetgrid.Blocks().size();
    std::vector<Mesh::BSExtent<Mesh::IndexType> > target_extents(NTargetRegions);
    target_mesh.resize(NTargetRegions);
    
    std::vector<std::vector<GeoPrim::CBox> > target_boxes(NTargetRegions);
    
    std::vector<FloGridBlock>::iterator target_block_iterator = targetgrid.Blocks().begin();
    std::vector<Mesh::BSExtent<Mesh::IndexType> >::iterator target_extent_iterator = target_extents.begin();
    std::vector<Mesh::Connectivity>::iterator target_mesh_iterator = target_mesh.begin();
    std::vector<std::vector<GeoPrim::CBox> >::iterator target_boxes_iterator = target_boxes.begin();
    while( target_block_iterator != targetgrid.Blocks().end())
      {
	target_block_iterator->MeshLayout();
	std::vector<Mesh::IndexType> target_block_extent(6,0);
	unsigned int isize = target_block_iterator->isize();
	unsigned int jsize = target_block_iterator->jsize();
	unsigned int ksize = target_block_iterator->ksize();
	target_block_extent[0] = target_block_extent[2] = target_block_extent[4] = 1;
	target_block_extent[1] = isize+1;
	target_block_extent[3] = jsize+1;
	target_block_extent[5] = ksize+1;
	target_extent_iterator->Init(target_block_extent);
	global.FunctionEntry("CreateUnsMesh");
	target_extent_iterator->CreateUnstructuredMesh(*target_mesh_iterator);
	target_mesh_iterator->Sync();
	target_mesh_iterator->SyncSizes();
	global.FunctionExit("CreateUnsMesh");
	target_boxes_iterator->resize(3);
	Mesh::NodalCoordinates nc;
	unsigned int nnodes = target_block_iterator->NNodes();
	total_target_nodes += nnodes;
	total_target_cells += (isize*jsize*ksize);
	nc.init(nnodes,&((target_block_iterator->Coords())[0]));
	global.FunctionEntry("Orient elements");
	Mesh::IndexType invertedcount = 0;
	unsigned int number_of_elements = isize*jsize*ksize;
	for(Mesh::IndexType element_being_processed = 1;
	    element_being_processed <= number_of_elements;
	    element_being_processed++)
	  {
	    Mesh::IndexType index = element_being_processed - 1;
	    Mesh::IndexType size_of_element = (*target_mesh_iterator)[index].size();
	    Mesh::GenericElement ge(size_of_element);
	    if(ge.Inverted((*target_mesh_iterator)[index],nc)){
	      invertedcount++;
	      ge.ReOrient((*target_mesh_iterator)[index]);
	    }
	  }
	global.FunctionExit("Orient elements");
	//      if(StdOut && verblevel) 
	//	*StdOut << "Number of elements reoriented: " << invertedcount 
	//		<< std::endl;
	global.FunctionEntry("BoxProcessing");
	Mesh::GetMeshBoxes(nc,*target_mesh_iterator,(*target_boxes_iterator)[0],(*target_boxes_iterator)[1],
			   (*target_boxes_iterator)[2]);
	targetbox.merge((*target_boxes_iterator)[0]);
	global.FunctionExit("BoxProcessing");
	target_boxes_iterator++;
	target_mesh_iterator++;
	target_extent_iterator++;
	target_block_iterator++;
      }
    
    unsigned int totalnodes = 0;
    unsigned int nregions = 0;
    double minx, maxx, miny, maxy, minz, maxz;
    minx = miny = minz = std::numeric_limits<double>::max();
    maxx = maxy = maxz = -1 * std::numeric_limits<double>::max();
    comm.Reduce(NTargetRegions,nregions,IRAD::Comm::DTUINT,IRAD::Comm::SUMOP,0);
    comm.Reduce(total_target_cells,totalcells,IRAD::Comm::DTUINT,IRAD::Comm::SUMOP,0);
    comm.Reduce(total_target_nodes,totalnodes,IRAD::Comm::DTUINT,IRAD::Comm::SUMOP,0);
    comm.Reduce(targetbox.P1().x(),minx,IRAD::Comm::DTDOUBLE,IRAD::Comm::MINOP,0);
    comm.Reduce(targetbox.P1().y(),miny,IRAD::Comm::DTDOUBLE,IRAD::Comm::MINOP,0);
    comm.Reduce(targetbox.P1().z(),minz,IRAD::Comm::DTDOUBLE,IRAD::Comm::MINOP,0);
    comm.Reduce(targetbox.P2().x(),maxx,IRAD::Comm::DTDOUBLE,IRAD::Comm::MAXOP,0);
    comm.Reduce(targetbox.P2().y(),maxy,IRAD::Comm::DTDOUBLE,IRAD::Comm::MAXOP,0);
    comm.Reduce(targetbox.P2().z(),maxz,IRAD::Comm::DTDOUBLE,IRAD::Comm::MAXOP,0);
    GeoPrim::CBox domainbox(GeoPrim::CPoint(minx,miny,minz),GeoPrim::CPoint(maxx,maxy,maxz));
    
    if(StdOut && verblevel > 0) {
      if(rank)
	*StdOut << "Target grid (" << rank << ") has " << NTargetRegions << " regions, " 
		<< total_target_nodes << " nodes, and " << total_target_cells
		<< " cells." << std::endl
		<< "With bounding box: " << targetbox << std::endl;
      else
	*StdOut << "Target grid has " << nregions << " regions, " 
		<< totalnodes << " nodes, and " << totalcells
		<< " cells." << std::endl
		<< "With bounding box: " << domainbox << std::endl;
    }
    
    
    // for every target region, which source regions
    // target_connectivity[target_region_index][1:nintersections] = source_region_id
    Mesh::Connectivity target_connectivity(NTargetRegions);
    
    // holds the intesections of the colliding target/source regions
    // target_intersections[target_region_index][1:nintersecting_regions] = intersecting_box
    std::vector< std::vector<GeoPrim::CBox> > target_intersections(NTargetRegions);
    
    
    Mesh::Connectivity::iterator target_connectivity_iterator = target_connectivity.begin();
    std::vector< std::vector<GeoPrim::CBox> >::iterator target_intersections_iterator = target_intersections.begin();
    target_boxes_iterator = target_boxes.begin();
    unsigned int target_id = 0;
    
    // This loop collides the target and source regions, thereby populating the 
    // target connectivity and target_intersections
    global.FunctionEntry("RegionIntersection");
    while(target_boxes_iterator != target_boxes.end()){
      target_id++;
      source_boxes_iterator = source_boxes.begin();
      unsigned int source_id = 0;
      //    std::cout << "Target region " << target_id << " collides with source regions: ";
      while(source_boxes_iterator != source_boxes.end()){
	source_id++;
	GeoPrim::CBox intersection = (*target_boxes_iterator)[0].intersect((*source_boxes_iterator)[0]);
	if(!intersection.empty()){
	  //	std::cout << source_id << " ";
	  target_connectivity_iterator->push_back(source_id);
	  target_intersections_iterator->push_back(intersection);
	}
	source_boxes_iterator++;
      }
      //    std::cout << std::endl;
      target_connectivity_iterator++;
      target_boxes_iterator++;
      target_intersections_iterator++;
    }
    global.FunctionExit("RegionIntersection");
    
    
    global.FunctionEntry("FindIntersectionCells");
    // Holds the nodes that participate in each intersecting region
    // target_intersection_nodes[1:ntargetregions][1:nintersection][1:nnodes_intersecting_region]
    std::vector<Mesh::Connectivity> target_intersection_cells(NTargetRegions);
    std::vector<std::vector<double> > target_region_cell_centers(NTargetRegions);
    
    // For every target region
    // For every target region's cell
    // a vector of pairs<source region id,cell id>
    // where the target cell is found
    std::vector<std::vector<std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> > > > target_source_matches(NTargetRegions);
    target_source_cells.resize(NTargetRegions);
    std::vector<Mesh::Connectivity>::iterator tint_cells_iterator = target_intersection_cells.begin();
    std::vector<std::vector<std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> > > >::iterator  tsmit = 
      target_source_matches.begin();
    target_connectivity_iterator = target_connectivity.begin();
    target_boxes_iterator = target_boxes.begin();
    target_block_iterator = targetgrid.Blocks().begin();
    target_mesh_iterator = target_mesh.begin();
    target_intersections_iterator = target_intersections.begin();
    while(tint_cells_iterator != target_intersection_cells.end())
      {
	//      std::ofstream Ouf;
	unsigned int target_region_index = tint_cells_iterator - target_intersection_cells.begin();
	unsigned int nnodes = target_block_iterator->NNodes();
	unsigned int ncells = target_block_iterator->NumCells();
	//      std::ostringstream Ostr;
	//      Ostr << "debug" << target_region_index+1 ;
	//      Ouf.open(Ostr.str().c_str());
	
	tsmit->resize(ncells);
	target_source_cells[target_region_index].resize(ncells);
	// std::vector<unsigned int> nfound(ncells,0);
	std::vector<bool> cell_found(ncells,false);
	Mesh::NodalCoordinates nc(nnodes,&(target_block_iterator->Coords()[0]));
	GeoPrim::CBox intersection_box;
	global.FunctionEntry("Centroids");
	Mesh::GetMeshCentroids(nc,*target_mesh_iterator,target_region_cell_centers[target_region_index]);
	global.FunctionExit("Centroids");
	tint_cells_iterator->Resize(target_intersections_iterator->size());
	Mesh::Connectivity::iterator intersection_cells_iterator = tint_cells_iterator->begin();
	std::vector<GeoPrim::CBox>::iterator intersections_iterator = target_intersections_iterator->begin();
	std::vector<Mesh::IndexType>::iterator source_region_ids = target_connectivity_iterator->begin();
	while(intersections_iterator != target_intersections_iterator->end())
	  {
	    // Initialize source mesh data structures
	    Mesh::IndexType source_region_index = *source_region_ids++ - 1;
	    Mesh::IndexType nsource_nodes = sourcegrid.Blocks()[source_region_index].NNodes();
	    Mesh::NodalCoordinates source_coords(nsource_nodes,&(sourcegrid.Blocks()[source_region_index].Coords()[0]));
	    Mesh::Connectivity &sourcemesh = source_mesh[source_region_index];
	    double source_edgelen = (source_boxes[source_region_index][1].P2() - source_boxes[source_region_index][1].P1()).norm();
	    // Obtain the source cells which live in this intersection
	    std::vector<Mesh::IndexType> source_cells;
	    global.FunctionEntry("CollideMeshBox");
	    Mesh::CollideMeshWithBox(source_coords,sourcemesh,*intersections_iterator,source_cells);
	    global.FunctionExit("CollideMeshBox");
	    std::vector<Mesh::IndexType> target_cells;
	    global.FunctionEntry("CollideMeshBox");
	    Mesh::CollideMeshWithBox(nc,*target_mesh_iterator,*intersections_iterator,target_cells);
	    global.FunctionExit("CollideMeshBox");
	    
	    if(source_cells.empty() || target_cells.empty()){
	      //	      std::cout << "found " << source_cells.size() << " source candidates." << std::endl
	      //			<< "found " << target_cells.size() << " target candidates." << std::endl
	      //			<< "source id = " << source_region_index+1 << ", box: " << *intersections_iterator 
	      //			<< std::endl << "source box = " << source_boxes[source_region_index][0] << std::endl;
	    }
	    else {
	      if(!rank)
		std::cout << ".";
	      std::vector<Mesh::IndexType>::iterator ttccii = target_cells.begin();
	      //	      for(int i = 0;i < ncells;i++){
	      while(ttccii != target_cells.end()){
		int i = *ttccii++ - 1;
		intersection_box.merge(*intersections_iterator);
		GeoPrim::CPoint p(&(target_region_cell_centers[target_region_index][3*i]));
		//		if(intersections_iterator->contains(p)){
		  // Ok, the target cell is in the intersection of the two regions
		  // Find source cell and get solution values
		  // 1. Extract the source cells close to the given point
		std::vector<Mesh::IndexType> candidates;
		GeoPrim::CBox little_box(source_boxes[source_region_index][1].around(p));
		// weed out the cells by checking distance from centroid
		global.FunctionEntry("WeedOut");
		std::vector<Mesh::IndexType>::iterator srccellit = source_cells.begin();
		while(srccellit != source_cells.end()){
		  GeoPrim::CPoint p2(&(source_cell_centers[source_region_index][3*(*srccellit-1)]));
		  if((p2-p).norm() <= source_edgelen)
		    candidates.push_back(*srccellit);
		  srccellit++;
		}
		global.FunctionExit("WeedOut");
		if(candidates.empty()){
		  global.FunctionEntry("CollideCellsBox");
		  Mesh::CollideCellsWithBox(source_coords,sourcemesh,little_box,source_cells,candidates);
		  global.FunctionExit("CollideCellsBox");
		}
		else{
		  std::vector<Mesh::IndexType> candcopy(candidates);
		  candidates.resize(0);
		  global.FunctionEntry("CollideCellsBoxMini");
		  Mesh::CollideCellsWithBox(source_coords,sourcemesh,little_box,candcopy,candidates);
		  global.FunctionExit("CollideCellsBoxMini");
		}
		GeoPrim::CVector natc;
		unsigned int source_cell_id = 0;
		if(candidates.empty()){
		  //		std::cout << "did not find candidates for target mesh point: (" << target_region_index+1
		  //			  << "," << i+1 << "," << p << ")  with box: " << little_box << " in source mesh: " 
		  //			  << source_region_index+1 << " with bounding box " << source_boxes[source_region_index][0]
		  //			  << std::endl;
		  global.FunctionEntry("FindPointInAllCells");
		  source_cell_id = Mesh::FindPointInCells(p,source_coords,sourcemesh,source_cells,natc);
		  global.FunctionExit("FindPointInAllCells");
		}
		else{
		  global.FunctionEntry("FindPointInCells");
		  source_cell_id = Mesh::FindPointInCells(p,source_coords,sourcemesh,candidates,natc);
		  global.FunctionExit("FindPointInCells");
		}
		if(source_cell_id > 0){
		  (*tsmit)[i].push_back(std::make_pair(source_region_index+1,source_cell_id));
		  intersection_cells_iterator->push_back(i+1);
		  cell_found[i] = true;
		}
		else{ // resort to something else terrible
		  //		Ouf << "Target Point: " << p << " with box [" << little_box << "] " << std::endl
		  //		    << "Source Region: " << source_region_index+1 << " with box [" 
		  //		    << source_boxes[source_region_index][0] << "]" << std::endl
		  //		    << "Source Elements: " << std::endl;
		  GeoPrim::CVector d(1000,1000,1000);
		  unsigned int closest_cell = 0;
		  std::vector<Mesh::IndexType> *thecells;
		  if(candidates.empty())
		    thecells = &source_cells;
		  else
		    thecells = &candidates;
		  //		  std::vector<Mesh::IndexType>::iterator ci = source_cells.begin();
		  //		  while(ci != source_cells.end())
		  std::vector<Mesh::IndexType>::iterator ci = thecells->begin();
		  while(ci != thecells->end())
		    {
		      source_cell_id = *ci++;
		      GeoPrim::CVector dist = p - GeoPrim::CPoint(&source_cell_centers[source_region_index][3*(source_cell_id-1)]);
		      if(dist.mag() < d.mag()){
			d = dist;
			closest_cell = source_cell_id;
		      }
		    }
		  //		Ouf << "Closest cell (" << source_region_index+1 << "," << closest_cell << ") : ("
		  //		    << d << ")" << std::endl;
		  (*tsmit)[i].push_back(std::make_pair(source_region_index+1,closest_cell));
		  intersection_cells_iterator->push_back(i+1);
		  // 		if(candidates.empty()){
		  // 		}
		  // 		else {
		  // 		  std::vector<Mesh::IndexType>::iterator ci = candidates.begin();
		  // 		  while(ci != candidates.end())
		  // 		    {
		  // 		      unsigned int source_cell_id = *ci++;
		  // 		      Ouf << "Element " << source_cell_id << ": " << std::endl;
		  // 		      std::vector<Mesh::IndexType>::iterator ni = sourcemesh[source_cell_id-1].begin();
		  // 		      while(ni != sourcemesh[source_cell_id-1].end())
		  // 			{
		  // 			  unsigned int node_id = *ni++;
		  // 			  Ouf << node_id << ": " 
		  // 			      << GeoPrim::C3Point(source_coords[node_id])
		  // 			      << std::endl;
		  // 			}
		  
		  // 		    }
		  
		  // 		}		
		}
	      }
	    }
	    intersection_cells_iterator++;
	    intersections_iterator++;
	  }
	    
	    
	//      std::vector<bool>::iterator cfi = cell_found.begin();
	//      unsigned int number_not_found = 0;
	//      while(cfi != cell_found.end())
	//	if(!*cfi++)
	//	  number_not_found++;
	
	if(!((*target_boxes_iterator)[0] == intersection_box) && StdOut)
	  *StdOut << "Mesh Box: " << (*target_boxes_iterator)[0] 
		  << " Intersecting region: " << intersection_box << std::endl;
	tsmit++;
	target_mesh_iterator++;
	target_block_iterator++;
	target_intersections_iterator++;
	tint_cells_iterator++;
	target_boxes_iterator++;
	target_connectivity_iterator++;
      }
    
    global.FunctionExit("FindIntersectionCells");
    //  }
    comm.Barrier();
    if(!rank && StdOut)
      std::cout << std::endl << "Preparing to write solution." << std::endl;
    global.FunctionEntry("Reporting");
    unsigned int number_not_found = 0;
    bool reporting = false;
    tsmit = target_source_matches.begin();
    while(tsmit != target_source_matches.end())
      {
	unsigned int target_region_index = tsmit - target_source_matches.begin();
	std::vector<std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> > >::iterator tci = tsmit->begin();
	while(tci != tsmit->end()){
	  unsigned int target_cell_index = tci - tsmit->begin();
	  if(tci->empty()){
	    number_not_found++;
	    if(debug && StdOut)
	      *StdOut << "(" << target_region_index+1 << "," << target_cell_index+1 << ") not found." << std::endl;
	  }
	  else {
	    if(debug && StdOut){
	      *StdOut << "(" << target_region_index+1 << "," << target_cell_index+1 << ") has " 
		      << tci->size() << " matches." << std::endl
		      << "target point: " 
		      << GeoPrim::C3Point(&target_region_cell_centers[target_region_index][3*target_cell_index]) << std::endl
		      << "source(s): " << std::endl;
	    }
	    if(tci->size() == 1)
	      target_source_cells[target_region_index][target_cell_index] = (*tci)[0];
	    else {
	      GeoPrim::CPoint p(&target_region_cell_centers[target_region_index][3*target_cell_index]);
	      std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> >::iterator tcit = tci->begin();
	      GeoPrim::CVector dist(1000,1000,1000);
	      std::pair<Mesh::IndexType,Mesh::IndexType> match;
	      while(tcit != tci->end()){
		Mesh::IndexType source_region_index = tcit->first - 1;
		Mesh::IndexType source_cell_index   = tcit->second - 1;
		GeoPrim::CVector d = p - GeoPrim::CPoint(&source_cell_centers[source_region_index][3*source_cell_index]);
		if(d.mag() < dist.mag()){
		  dist = d;
		  match = *tcit;
		}
		tcit++;
	      }
	      target_source_cells[target_region_index][target_cell_index] = match;
	    }
	    std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> >::iterator tcit = tci->begin();
	    while(tcit != tci->end())
	      {
		if(debug && StdOut){
		  *StdOut << "(" << tcit->first << "," << tcit->second << ") : (" 
			  << GeoPrim::C3Point(&source_cell_centers[tcit->first-1][3*(tcit->second-1)])
			  << ") : (" 
			  << GeoPrim::C3Point(&source_cell_centers[tcit->first-1][3*(tcit->second-1)]) - 
		    GeoPrim::C3Point(&target_region_cell_centers[target_region_index][3*target_cell_index])
			  << ")" << std::endl;
		}
		tcit++;
	      }
	    
	  }
	  tci++;
	}
	tsmit++;
      }
    if(number_not_found > 0 && StdOut)
      *StdOut << number_not_found << "/" << total_target_cells << " cells not found or processed." << std::endl;
    global.FunctionExit("Reporting");
  }
  // Searching complete... Now do the solution transfer part 

  // get rid of the rudding meshes
  //  sourcegrid.DestroyGrids();
  //  targetgrid.DestroyGrids();

  global.FunctionEntry("SolutionRead");
  // read the source solution
  if(sourcegrid.ReadAllSolutions(timestamp)){
    std::cout << "Could not read source solutions." << std::endl;
    exit(1);
  }
  global.FunctionExit("SolutionRead");

  // allocate space for the target solution
  unsigned int ng = sourcegrid.Blocks()[0].NGhostLayers();
  double unknown_number = sourcegrid.UnknownNumber();
  double time           = sourcegrid.Time();
  targetgrid.SetGhostLayers(ng);
  targetgrid.CreateSolutions();
  std::vector<FloGridBlock> &targetgrids = targetgrid.Blocks();

  // Populate an array of flat source cell ids for performance
  source_cell_ids.resize(NSourceRegions);
  std::vector<FloGridBlock> &sourcegrids = sourcegrid.Blocks();
  std::vector<FloGridBlock>::iterator sbi = sourcegrids.begin();
  while(sbi != sourcegrids.end()){
    unsigned int source_grid_index = sbi - sourcegrids.begin();
    Mesh::BSExtent<Mesh::IndexType> real_extent;
    Mesh::BSExtent<Mesh::IndexType> ghost_extent;
    std::vector<Mesh::IndexType> source_block_extent(6,0);
    std::vector<Mesh::IndexType> source_ghost_extent(6,0);
    unsigned int isize = sbi->isize();
    unsigned int jsize = sbi->jsize();
    unsigned int ksize = sbi->ksize();
    source_block_extent[0] = source_block_extent[2] = source_block_extent[4] = 1 + ng;
    source_ghost_extent[0] = source_ghost_extent[2] = source_ghost_extent[4] = 1;
    source_block_extent[1] = isize+ng;
    source_block_extent[3] = jsize+ng;
    source_block_extent[5] = ksize+ng;
    source_ghost_extent[1] = isize+(2*ng);
    source_ghost_extent[3] = jsize+(2*ng);
    source_ghost_extent[5] = ksize+(2*ng);
    real_extent.Init(source_block_extent);
    ghost_extent.Init(source_ghost_extent);
    ghost_extent.GetFlatIndices(real_extent,source_cell_ids[source_grid_index]);
    sbi++;
  }

  // Transfer the data
  global.FunctionEntry("DataTransfer");
  target_cell_ids.resize(NTargetRegions);
  std::vector<std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> > >::iterator tsci = target_source_cells.begin();
  while(tsci != target_source_cells.end()){
    unsigned int target_region_index = tsci - target_source_cells.begin();
    //    std::vector<Mesh::IndexType> &target_cell_ids = all_target_cell_ids[target_region_index];
    std::vector<Mesh::IndexType> &target_cell_ids_local = target_cell_ids[target_region_index];
    Mesh::BSExtent<Mesh::IndexType> real_extent;
    Mesh::BSExtent<Mesh::IndexType> ghost_extent;
    std::vector<Mesh::IndexType> target_block_extent(6,0);
    std::vector<Mesh::IndexType> target_ghost_extent(6,0);
    unsigned int isize = targetgrids[target_region_index].isize();
    unsigned int jsize = targetgrids[target_region_index].jsize();
    unsigned int ksize = targetgrids[target_region_index].ksize();
    target_block_extent[0] = target_block_extent[2] = target_block_extent[4] = 1 + ng;
    target_ghost_extent[0] = target_ghost_extent[2] = target_ghost_extent[4] = 1;
    target_block_extent[1] = isize+ng;
    target_block_extent[3] = jsize+ng;
    target_block_extent[5] = ksize+ng;
    target_ghost_extent[1] = isize+(2*ng);
    target_ghost_extent[3] = jsize+(2*ng);
    target_ghost_extent[5] = ksize+(2*ng);
    real_extent.Init(target_block_extent);
    ghost_extent.Init(target_ghost_extent);
    //    ghost_extent.GetFlatIndices(target_block_extent,target_cell_ids_local);
    ghost_extent.GetFlatIndices(real_extent,target_cell_ids_local);
    unsigned int target_block = (isize+(2*ng))*(jsize+(2*ng))*(ksize+(2*ng));
    //    ghost_extent.GetFlatIndices(real_extent,all_target_cell_ids[target_region_index]);
    std::vector<std::pair<Mesh::IndexType,Mesh::IndexType> >::iterator tci = tsci->begin();
    while(tci != tsci->end()){
      unsigned int target_cell_index = target_cell_ids_local[tci - tsci->begin()] - 1;
      unsigned int source_region_index = tci->first - 1;
      unsigned int source_cell_index = source_cell_ids[source_region_index][tci->second -1] - 1;
      unsigned int source_block = (sourcegrids[source_region_index].isize()+(2*ng))*
	(sourcegrids[source_region_index].jsize()+(2*ng))*(sourcegrids[source_region_index].ksize()+(2*ng));

      targetgrids[target_region_index].Solution()[target_cell_index] =
	sourcegrids[source_region_index].Solution()[source_cell_index];

      targetgrids[target_region_index].Solution()[target_block+target_cell_index] =
	sourcegrids[source_region_index].Solution()[source_block+source_cell_index];

      targetgrids[target_region_index].Solution()[(2*target_block) + target_cell_index] =
	sourcegrids[source_region_index].Solution()[(2*source_block) + source_cell_index];

      targetgrids[target_region_index].Solution()[(3*target_block) + target_cell_index] =
	sourcegrids[source_region_index].Solution()[(3*source_block) + source_cell_index];

      targetgrids[target_region_index].Solution()[(4*target_block) + target_cell_index] =
	sourcegrids[source_region_index].Solution()[(4*source_block) + source_cell_index];

      tci++;
    }
    tsci++;
  }
  global.FunctionExit("DataTransfer");
  global.FunctionEntry("GhostZoneFlood");
  std::vector<FloGridBlock>::iterator tgi = targetgrids.begin();
  while(tgi != targetgrids.end()){
    unsigned int target_grid_index = tgi - targetgrids.begin();
    Mesh::BSExtent<Mesh::IndexType> full_extent;
    std::vector<Mesh::IndexType> flat_full_extent(6,0);
    unsigned int isize = tgi->isize();
    unsigned int jsize = tgi->jsize();
    unsigned int ksize = tgi->ksize();
    flat_full_extent[0] = flat_full_extent[2] = flat_full_extent[4] = 1;
    flat_full_extent[1] = isize+(2*ng);
    flat_full_extent[3] = jsize+(2*ng);
    flat_full_extent[5] = ksize+(2*ng);
    full_extent.Init(flat_full_extent);
    unsigned int target_block = (isize+(2*ng))*(jsize+(2*ng))*(ksize+(2*ng));
    bool do_iplane = true;
    if(do_iplane){
      Mesh::BSExtent<Mesh::IndexType> plane_extent;
      std::vector<Mesh::IndexType> flat_plane_extent(6,0);
      flat_plane_extent[0] = flat_plane_extent[1] = ng+1;
      flat_plane_extent[2] = flat_plane_extent[4] = ng+1;
      flat_plane_extent[3] = jsize + ng;
      flat_plane_extent[5] = ksize + ng;
      plane_extent.Init(flat_plane_extent);
      std::vector<Mesh::IndexType> real_cell_ids;
      full_extent.GetFlatIndices(plane_extent,real_cell_ids);
      for(int ii = 1;ii <= ng;ii++){
	std::vector<Mesh::IndexType> flat_ghost_layer(6,0);
	flat_ghost_layer[0] = flat_ghost_layer[1] = ii;
	flat_ghost_layer[2] = flat_plane_extent[2];
	flat_ghost_layer[3] = flat_plane_extent[3];
	flat_ghost_layer[4] = flat_plane_extent[4];
	flat_ghost_layer[5] = flat_plane_extent[5];
	Mesh::BSExtent<Mesh::IndexType> ghost_plane_extent(flat_ghost_layer);
	std::vector<Mesh::IndexType> ghost_cell_ids;
	full_extent.GetFlatIndices(ghost_plane_extent,ghost_cell_ids);
	assert(ghost_cell_ids.size() == real_cell_ids.size());
	std::vector<Mesh::IndexType>::iterator rci = real_cell_ids.begin();
	std::vector<Mesh::IndexType>::iterator gci = ghost_cell_ids.begin();
	while(rci != real_cell_ids.end()){
	  unsigned int real_cell_index  = *rci++ - 1;
	  unsigned int ghost_cell_index = *gci++ - 1; 
	  targetgrids[target_grid_index].Solution()[ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[real_cell_index];
	  targetgrids[target_grid_index].Solution()[target_block+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[target_block+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(2*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(2*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(3*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(3*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(4*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(4*target_block)+real_cell_index];
	}
      }
    } 
    if(do_iplane){
      Mesh::BSExtent<Mesh::IndexType> plane_extent;
      std::vector<Mesh::IndexType> flat_plane_extent(6,0);
      flat_plane_extent[0] = flat_plane_extent[1] = isize+ng;
      flat_plane_extent[2] = flat_plane_extent[4] = ng+1;
      flat_plane_extent[3] = jsize + ng;
      flat_plane_extent[5] = ksize + ng;
      plane_extent.Init(flat_plane_extent);
      std::vector<Mesh::IndexType> real_cell_ids;
      full_extent.GetFlatIndices(plane_extent,real_cell_ids);
      for(int ii = 1;ii <= ng;ii++){
 	std::vector<Mesh::IndexType> flat_ghost_layer(6,0);
 	flat_ghost_layer[0] = flat_ghost_layer[1] = isize+ng+ii;
 	flat_ghost_layer[2] = flat_plane_extent[2];
 	flat_ghost_layer[3] = flat_plane_extent[3];
 	flat_ghost_layer[4] = flat_plane_extent[4];
 	flat_ghost_layer[5] = flat_plane_extent[5];
 	Mesh::BSExtent<Mesh::IndexType> ghost_plane_extent(flat_ghost_layer);
 	std::vector<Mesh::IndexType> ghost_cell_ids;
 	full_extent.GetFlatIndices(ghost_plane_extent,ghost_cell_ids);
 	std::vector<Mesh::IndexType>::iterator rci = real_cell_ids.begin();
 	std::vector<Mesh::IndexType>::iterator gci = ghost_cell_ids.begin();
	assert(ghost_cell_ids.size() == real_cell_ids.size());
	while(rci != real_cell_ids.end()){
	  unsigned int real_cell_index  = *rci++ - 1;
	  unsigned int ghost_cell_index = *gci++ - 1; 
	  targetgrids[target_grid_index].Solution()[ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[real_cell_index];
	  targetgrids[target_grid_index].Solution()[target_block+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[target_block+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(2*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(2*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(3*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(3*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(4*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(4*target_block)+real_cell_index];
	}
	// 	while(rci != real_cell_ids.end())
	// 	  *gci++ = *rci++;
      }
    } 
    bool do_jplane = true;
    if(do_jplane){
      Mesh::BSExtent<Mesh::IndexType> plane_extent;
      std::vector<Mesh::IndexType> flat_plane_extent(6,0);
      flat_plane_extent[0] = 1;
      flat_plane_extent[1] = isize+(2*ng);
      flat_plane_extent[2] = flat_plane_extent[3] = ng+1; 
      flat_plane_extent[4] = ng+1;
      flat_plane_extent[5] = ksize + ng;
      plane_extent.Init(flat_plane_extent);
      std::vector<Mesh::IndexType> real_cell_ids;
      full_extent.GetFlatIndices(plane_extent,real_cell_ids);
      for(int ii = 1;ii <= ng;ii++){
 	std::vector<Mesh::IndexType> flat_ghost_layer(6,0);
 	flat_ghost_layer[0] = flat_plane_extent[0];
 	flat_ghost_layer[1] = flat_plane_extent[1];
 	flat_ghost_layer[2] = flat_ghost_layer[3] = ii;
 	flat_ghost_layer[4] = flat_plane_extent[4];
 	flat_ghost_layer[5] = flat_plane_extent[5];
 	Mesh::BSExtent<Mesh::IndexType> ghost_plane_extent(flat_ghost_layer);
 	std::vector<Mesh::IndexType> ghost_cell_ids;
 	full_extent.GetFlatIndices(ghost_plane_extent,ghost_cell_ids);
 	std::vector<Mesh::IndexType>::iterator rci = real_cell_ids.begin();
 	std::vector<Mesh::IndexType>::iterator gci = ghost_cell_ids.begin();
	assert(ghost_cell_ids.size() == real_cell_ids.size());
	while(rci != real_cell_ids.end()){
	  unsigned int real_cell_index  = *rci++ - 1;
	  unsigned int ghost_cell_index = *gci++ - 1; 
	  targetgrids[target_grid_index].Solution()[ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[real_cell_index];
	  targetgrids[target_grid_index].Solution()[target_block+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[target_block+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(2*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(2*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(3*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(3*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(4*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(4*target_block)+real_cell_index];
	}
	// 	while(rci != real_cell_ids.end())
	// 	  *gci++ = *rci++;
      }
    } 
    if(do_jplane){
      Mesh::BSExtent<Mesh::IndexType> plane_extent;
      std::vector<Mesh::IndexType> flat_plane_extent(6,0);
      flat_plane_extent[0] = 1;
      flat_plane_extent[1] = isize+(2*ng);
      flat_plane_extent[2] = flat_plane_extent[3] = jsize+ng; 
      flat_plane_extent[4] = ng+1;
      flat_plane_extent[5] = ksize + ng;
      plane_extent.Init(flat_plane_extent);
      std::vector<Mesh::IndexType> real_cell_ids;
      full_extent.GetFlatIndices(plane_extent,real_cell_ids);
      for(int ii = 1;ii <= ng;ii++){
 	std::vector<Mesh::IndexType> flat_ghost_layer(6,0);
 	flat_ghost_layer[0] = flat_plane_extent[0];
 	flat_ghost_layer[1] = flat_plane_extent[1];
 	flat_ghost_layer[2] = flat_ghost_layer[3] = ii+jsize+ng;
 	flat_ghost_layer[4] = flat_plane_extent[4];
 	flat_ghost_layer[5] = flat_plane_extent[5];
 	Mesh::BSExtent<Mesh::IndexType> ghost_plane_extent(flat_ghost_layer);
 	std::vector<Mesh::IndexType> ghost_cell_ids;
 	full_extent.GetFlatIndices(ghost_plane_extent,ghost_cell_ids);
 	std::vector<Mesh::IndexType>::iterator rci = real_cell_ids.begin();
 	std::vector<Mesh::IndexType>::iterator gci = ghost_cell_ids.begin();
	assert(ghost_cell_ids.size() == real_cell_ids.size());
	while(rci != real_cell_ids.end()){
	  unsigned int real_cell_index  = *rci++ - 1;
	  unsigned int ghost_cell_index = *gci++ - 1; 
	  targetgrids[target_grid_index].Solution()[ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[real_cell_index];
	  targetgrids[target_grid_index].Solution()[target_block+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[target_block+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(2*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(2*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(3*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(3*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(4*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(4*target_block)+real_cell_index];
	}
	// 	while(rci != real_cell_ids.end())
	// 	  *gci++ = *rci++;
      }
    } 
    bool do_kplane = true;
    if(do_kplane){
      Mesh::BSExtent<Mesh::IndexType> plane_extent;
      std::vector<Mesh::IndexType> flat_plane_extent(6,0);
      flat_plane_extent[0] = 1;
      flat_plane_extent[1] = isize+(2*ng);
      flat_plane_extent[2] = 1;
      flat_plane_extent[3] = jsize+(2*ng);
      flat_plane_extent[4] = flat_plane_extent[5] = ng+1;
      plane_extent.Init(flat_plane_extent);
      std::vector<Mesh::IndexType> real_cell_ids;
      full_extent.GetFlatIndices(plane_extent,real_cell_ids);
      for(int ii = 1;ii <= ng;ii++){
 	std::vector<Mesh::IndexType> flat_ghost_layer(6,0);
 	flat_ghost_layer[0] = flat_plane_extent[0];
 	flat_ghost_layer[1] = flat_plane_extent[1];
 	flat_ghost_layer[2] = flat_plane_extent[2];
 	flat_ghost_layer[3] = flat_plane_extent[3];
 	flat_ghost_layer[4] = flat_ghost_layer[5] = ii;
 	Mesh::BSExtent<Mesh::IndexType> ghost_plane_extent(flat_ghost_layer);
 	std::vector<Mesh::IndexType> ghost_cell_ids;
 	full_extent.GetFlatIndices(ghost_plane_extent,ghost_cell_ids);
 	std::vector<Mesh::IndexType>::iterator rci = real_cell_ids.begin();
 	std::vector<Mesh::IndexType>::iterator gci = ghost_cell_ids.begin();
	assert(ghost_cell_ids.size() == real_cell_ids.size());
	while(rci != real_cell_ids.end()){
	  unsigned int real_cell_index  = *rci++ - 1;
	  unsigned int ghost_cell_index = *gci++ - 1; 
	  targetgrids[target_grid_index].Solution()[ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[real_cell_index];
	  targetgrids[target_grid_index].Solution()[target_block+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[target_block+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(2*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(2*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(3*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(3*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(4*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(4*target_block)+real_cell_index];
	}
	// 	while(rci != real_cell_ids.end())
	// 	  *gci++ = *rci++;
      }
    } 
    if(do_kplane){
      Mesh::BSExtent<Mesh::IndexType> plane_extent;
      std::vector<Mesh::IndexType> flat_plane_extent(6,0);
      flat_plane_extent[0] = 1;
      flat_plane_extent[1] = isize+(2*ng);
      flat_plane_extent[2] = 1;
      flat_plane_extent[3] = jsize+(2*ng);
      flat_plane_extent[4] = flat_plane_extent[5] = ng+ksize;
      plane_extent.Init(flat_plane_extent);
      std::vector<Mesh::IndexType> real_cell_ids;
      full_extent.GetFlatIndices(plane_extent,real_cell_ids);
      for(int ii = 1;ii <= ng;ii++){
 	std::vector<Mesh::IndexType> flat_ghost_layer(6,0);
 	flat_ghost_layer[0] = flat_plane_extent[0];
 	flat_ghost_layer[1] = flat_plane_extent[1];
 	flat_ghost_layer[2] = flat_plane_extent[2];
 	flat_ghost_layer[3] = flat_plane_extent[3];
 	flat_ghost_layer[4] = flat_ghost_layer[5] = ii+ng+ksize;
 	Mesh::BSExtent<Mesh::IndexType> ghost_plane_extent(flat_ghost_layer);
 	std::vector<Mesh::IndexType> ghost_cell_ids;
 	full_extent.GetFlatIndices(ghost_plane_extent,ghost_cell_ids);
 	std::vector<Mesh::IndexType>::iterator rci = real_cell_ids.begin();
 	std::vector<Mesh::IndexType>::iterator gci = ghost_cell_ids.begin();
	assert(ghost_cell_ids.size() == real_cell_ids.size());
	while(rci != real_cell_ids.end()){
	  unsigned int real_cell_index  = *rci++ - 1;
	  unsigned int ghost_cell_index = *gci++ - 1; 
	  targetgrids[target_grid_index].Solution()[ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[real_cell_index];
	  targetgrids[target_grid_index].Solution()[target_block+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[target_block+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(2*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(2*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(3*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(3*target_block)+real_cell_index];
	  targetgrids[target_grid_index].Solution()[(4*target_block)+ghost_cell_index] =
	    targetgrids[target_grid_index].Solution()[(4*target_block)+real_cell_index];
	}
	// 	while(rci != real_cell_ids.end())
	// 	  *gci++ = *rci++;
      }
    } 
    tgi++;
  }
  
  global.FunctionExit("GhostZoneFlood");

  unsigned int ncells = 0;
  if(rank == 0) {
    std::vector<double> sourcemin(5,std::numeric_limits<double>::max());
    std::vector<double> sourcemax(5,-1.0*std::numeric_limits<double>::max());
    std::vector<double> sourcemean(5,0);
    std::vector<std::vector<Mesh::IndexType> >::iterator sri = source_cell_ids.begin();
    while(sri != source_cell_ids.end()){
      unsigned int source_region_index = sri - source_cell_ids.begin();
      std::vector<Mesh::IndexType>::iterator sci = sri->begin();
      unsigned int blocksize = (sourcegrids[source_region_index].isize()+(2*ng)) *
	(sourcegrids[source_region_index].jsize()+(2*ng)) * (sourcegrids[source_region_index].ksize()+(2*ng));
      while(sci != sri->end()){
	ncells++;
	unsigned int source_cell_index = *sci++-1;
      
	double val = sourcegrids[source_region_index].Solution()[source_cell_index];
	sourcemean[0] += val;
	if(val > sourcemax[0])
	  sourcemax[0] = val;
	if(val < sourcemin[0])
	  sourcemin[0] = val;
	
	val = sourcegrids[source_region_index].Solution()[blocksize+source_cell_index];
	sourcemean[1] += val;
	if(val > sourcemax[1])
	  sourcemax[1] = val;
	if(val < sourcemin[1])
	  sourcemin[1] = val;
	
	val = sourcegrids[source_region_index].Solution()[(2*blocksize)+source_cell_index];
	sourcemean[2] += val;
	if(val > sourcemax[2])
	  sourcemax[2] = val;
	if(val < sourcemin[2])
	  sourcemin[2] = val;
      
	val = sourcegrids[source_region_index].Solution()[(3*blocksize)+source_cell_index];
	sourcemean[3] += val;
	if(val > sourcemax[3])
	  sourcemax[3] = val;
	if(val < sourcemin[3])
	  sourcemin[3] = val;
	
	val = sourcegrids[source_region_index].Solution()[(4*blocksize)+source_cell_index];
	sourcemean[4] += val;
	if(val > sourcemax[4])
	  sourcemax[4] = val;
	if(val < sourcemin[4])
	  sourcemin[4] = val;
	
      }
      sri++;
    }
    if(ncells != total_source_cells)
      std::cout << "Warning: Source soln extents indicate only " << ncells 
		<< "/" << total_source_cells << " real cells." << std::endl;
    double scale = 1.0/static_cast<double>(ncells);
    sourcemean[0] *= scale;
    sourcemean[1] *= scale;
    sourcemean[2] *= scale;
    sourcemean[3] *= scale;
    sourcemean[4] *= scale;
    std::cout << "Source Data (min,max,mean):" << std::endl
	      << "rho  : (" << sourcemin[0] << "," << sourcemax[0] 
	      << "," << sourcemean[0] << ")" << std::endl
	      << "rho-u: (" << sourcemin[1] << "," << sourcemax[1] 
	      << "," << sourcemean[1] << ")" << std::endl
	      << "rho-v: (" << sourcemin[2] << "," << sourcemax[2] 
	      << "," << sourcemean[2] << ")" << std::endl
	      << "rho-w: (" << sourcemin[3] << "," << sourcemax[3] 
	      << "," << sourcemean[3] << ")" << std::endl
	      << "rho-E: (" << sourcemin[4] << "," << sourcemax[4] 
	      << "," << sourcemean[4] << ")" << std::endl;
  }
  std::vector<double> targetmin(5,std::numeric_limits<double>::max());
  std::vector<double> ttargetmin(5,std::numeric_limits<double>::max());
  std::vector<double> targetmax(5,-1*std::numeric_limits<double>::max());
  std::vector<double> ttargetmax(5,-1*std::numeric_limits<double>::max());
  std::vector<double> targetmean(5,0);
  std::vector<double> ttargetmean(5,0);
  ncells = 0;
  std::vector<std::vector<Mesh::IndexType> >::iterator sri = target_cell_ids.begin();
  while(sri != target_cell_ids.end()){
    unsigned int target_region_index = sri - target_cell_ids.begin();
    std::vector<Mesh::IndexType>::iterator sci = sri->begin();
    unsigned int blocksize = (targetgrids[target_region_index].isize()+(2*ng)) *
      (targetgrids[target_region_index].jsize()+(2*ng))*(targetgrids[target_region_index].ksize()+(2*ng));
    while(sci != sri->end()){
      ncells++;
      unsigned int target_cell_index = *sci++-1;
      double val = targetgrids[target_region_index].Solution()[target_cell_index];
      targetmean[0] += val;
      if(val > targetmax[0])
	targetmax[0] = val;
      if(val < targetmin[0])
	targetmin[0] = val;
      val = targetgrids[target_region_index].Solution()[blocksize+ target_cell_index];
      targetmean[1] += val;
      if(val > targetmax[1])
	targetmax[1] = val;
      if(val < targetmin[1])
	targetmin[1] = val;
      val = targetgrids[target_region_index].Solution()[2*blocksize+target_cell_index];
      targetmean[2] += val;
      if(val > targetmax[2])
	targetmax[2] = val;
      if(val < targetmin[2])
	targetmin[2] = val;
      val = targetgrids[target_region_index].Solution()[3*blocksize+target_cell_index];
      targetmean[3] += val;
      if(val > targetmax[3])
	targetmax[3] = val;
      if(val < targetmin[3])
	targetmin[3] = val;
      val = targetgrids[target_region_index].Solution()[4*blocksize +target_cell_index];
      targetmean[4] += val;
      if(val > targetmax[4])
	targetmax[4] = val;
      if(val < targetmin[4])
	targetmin[4] = val;
      
    }
    sri++;
  }
  if(ncells != total_target_cells)
    std::cout << "Warning: Target soln extents indicate only " << ncells 
	      << "/" << total_target_cells << " real cells." << std::endl;

  comm.Reduce(targetmean,ttargetmean,IRAD::Comm::DTDOUBLE,IRAD::Comm::SUMOP,0);
  comm.Reduce(targetmin,ttargetmin,IRAD::Comm::DTDOUBLE,IRAD::Comm::MINOP,0);
  comm.Reduce(targetmax,ttargetmax,IRAD::Comm::DTDOUBLE,IRAD::Comm::MAXOP,0);
  double scale = 1.0/static_cast<double>(ncells);
  targetmean[0] *= scale;
  targetmean[1] *= scale;
  targetmean[2] *= scale;
  targetmean[3] *= scale;
  targetmean[4] *= scale;
  scale = 1.0/static_cast<double>(totalcells);
  ttargetmean[0] *= scale;
  ttargetmean[1] *= scale;
  ttargetmean[2] *= scale;
  ttargetmean[3] *= scale;
  ttargetmean[4] *= scale;
  if(verblevel > 0 && StdOut){
    if( (nproc==1) || ((rank > 0) && debug)){
      *StdOut << "Target Data (min,max,mean):" << std::endl
	      << "rho  : (" << targetmin[0] << "," << targetmax[0] 
	      << "," << targetmean[0] << ")" << std::endl
	      << "rho-u: (" << targetmin[1] << "," << targetmax[1] 
	      << "," << targetmean[1] << ")" << std::endl
	      << "rho-v: (" << targetmin[2] << "," << targetmax[2] 
	      << "," << targetmean[2] << ")" << std::endl
	      << "rho-w: (" << targetmin[3] << "," << targetmax[3] 
	      << "," << targetmean[3] << ")" << std::endl
	      << "rho-E: (" << targetmin[4] << "," << targetmax[4] 
	      << "," << targetmean[4] << ")" << std::endl;
    }
    else {
      if(debug){
	*StdOut << "Target Data (min,max,mean):" << std::endl
		<< "rho  : (" << targetmin[0] << "," << targetmax[0] 
		<< "," << targetmean[0] << ")" << std::endl
		<< "rho-u: (" << targetmin[1] << "," << targetmax[1] 
		<< "," << targetmean[1] << ")" << std::endl
		<< "rho-v: (" << targetmin[2] << "," << targetmax[2] 
		<< "," << targetmean[2] << ")" << std::endl
		<< "rho-w: (" << targetmin[3] << "," << targetmax[3] 
		<< "," << targetmean[3] << ")" << std::endl
		<< "rho-E: (" << targetmin[4] << "," << targetmax[4] 
		<< "," << targetmean[4] << ")" << std::endl;
      }
	
      *StdOut << "Target Data (min,max,mean):" << std::endl
	      << "rho  : (" << ttargetmin[0] << "," << ttargetmax[0] 
	      << "," << ttargetmean[0] << ")" << std::endl
	      << "rho-u: (" << ttargetmin[1] << "," << ttargetmax[1] 
	      << "," << ttargetmean[1] << ")" << std::endl
	      << "rho-v: (" << ttargetmin[2] << "," << ttargetmax[2] 
	      << "," << ttargetmean[2] << ")" << std::endl
	      << "rho-w: (" << ttargetmin[3] << "," << ttargetmax[3] 
	      << "," << ttargetmean[3] << ")" << std::endl
	      << "rho-E: (" << ttargetmin[4] << "," << ttargetmax[4] 
	      << "," << ttargetmean[4] << ")" << std::endl;
    }
    
  }
  
  global.FunctionEntry("SolutionWrite");
  if(rank == 0)
    targetgrid.OpenSolutionFile(timestamp,time,unknown_number);
  std::ostringstream Ostr;
  Ostr << std::scientific << std::setprecision(16);
  int datasize = targetgrid.WriteBlocks(Ostr);
  std::vector<int> datasizes(nproc,0);
  comm.Gather(datasize,datasizes);
  if(rank > 0)
    comm._Send(const_cast<char *>(Ostr.str().c_str()),datasize,0,0);
  if(!rank){
    std::ofstream &Ouf = targetgrid.SolnFile();
    targetgrid.WriteBlocks(Ouf);
    for(unsigned int npi = 1;npi < nproc;npi++)
      {
	std::vector<char> recvbuf(datasizes[npi]);
	comm._Recv(&recvbuf[0],datasizes[npi],npi,0);
	std::vector<char>::iterator rbi = recvbuf.begin();
	while(rbi != recvbuf.end())
	  Ouf << *rbi++;
      }
    targetgrid.CloseSolutionFile();
  }
  comm.Barrier();
  //  targetgrid.WriteSolutionFile(timestamp,time,unknown_number);
  global.FunctionExit("SolutionWrite");
  
//    tgi = targetgrids.begin();
//    while(tgi != targetgrids.end()){
//      Mesh::NodalCoordinates nc;
//      unsigned int nnodes = tgi->NNodes();
//      nc.init(nnodes,&((tgi->Coords())[0]));
//      unsigned int block_index = tgi - targetgrids.begin();
//      unsigned int blocksize = (targetgrids[block_index].isize()+(2*ng)) *
//        (targetgrids[block_index].jsize()+(2*ng))*(targetgrids[block_index].ksize()+(2*ng));
//      std::ostringstream Ostr;
//      Ostr << "target_" << block_index+1 << ".vtk";
//      std::vector<double> emptyv;
//      emptyv.resize(0);
//      writeVtkData(nc,target_mesh[block_index],Ostr.str(),tgi->Solution(),target_cell_ids[block_index],blocksize);
//      tgi++;
//    }
//    tgi = sourcegrids.begin();
//    while(tgi != sourcegrids.end()){
//      Mesh::NodalCoordinates nc;
//      unsigned int nnodes = tgi->NNodes();
//      nc.init(nnodes,&((tgi->Coords())[0]));
//      unsigned int block_index = tgi - sourcegrids.begin();
//      unsigned int blocksize = (sourcegrids[block_index].isize()+(2*ng)) *
//        (sourcegrids[block_index].jsize()+(2*ng))*(sourcegrids[block_index].ksize()+(2*ng));
//      std::ostringstream Ostr;
//      Ostr << "source_" << block_index+1 << ".vtk";
//      std::vector<double> emptyv;
//      emptyv.resize(0);
//      writeVtkData(nc,source_mesh[block_index],Ostr.str(),tgi->Solution(),source_cell_ids[block_index],blocksize);
//      tgi++;
//    }

  //    test();
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // This stuff goes at the very end of the program, don't modify
  // below here.
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  global.FunctionExit("main");
  comm.Barrier();
  if(LogOut && debug)
    *LogOut << "All processors made it to the end." << std::endl;
  global.Finalize();
  //  global.Profiler.Finalize();
  if(StdOut)
    //    global.Report(Profiler.SummarizeSerialExecution(*StdOut);
    global.Report(*StdOut);
  
  return(0);
}
