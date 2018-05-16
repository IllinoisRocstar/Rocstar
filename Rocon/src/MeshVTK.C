
#include "MeshVTK.H"

namespace Mesh {
  /*
   * @brief Writes VTK files for the corresponding meshes.
   * @param meshes List of meshes in the current window.
   */
  void writeVtkFiles( std::vector< Mesh::UnstructuredMesh > &meshes )
  {
    
    // TODO: Generalize this function for any type of mesh.
    for( int m=0; m < meshes.size( ); ++m )
      {
	std::ostringstream oss; oss.clear( );
	oss << Program.outputfile << "." << Program.paneIds[ m ] << ".vtk";
	
	std::ofstream ofs;
	ofs.open( oss.str( ).c_str( ) );
	
	if( !ofs.is_open( ) )
	  {
	    std::cerr << "Cannot write VTK file: " << oss.str( ) << std::endl;
	    std::cerr << "File: " << __FILE__ << std::endl;
	    std::cerr << "Line: " << __LINE__ << std::endl;
	    assert( false );
	  }
	
	ofs << "# vtk DataFile Version 3.0\n";
	ofs << oss.str( ) << std::endl;
	ofs << "ASCII\n";
	ofs << "DATASET UNSTRUCTURED_GRID\n";
	ofs << "POINTS " << meshes[ m ].nc.Size( ) << " double\n";
	
	for( int i=1; i <= meshes[ m ].nc.Size( ); ++i )
	  {
	    ofs << meshes[ m ].nc.x( i ) << " ";
	    ofs << meshes[ m ].nc.y( i ) << " ";
	    ofs << meshes[ m ].nc.z( i ) << "\n";
	  }
	
	ofs << "CELLS " << meshes[ m ].con.Nelem( ) << " " << meshes[ m ].con.Nelem( )*9 << "\n";
	for( int e=1; e <= meshes[ m ].con.Nelem( ); ++e )
	  {
	    ofs <<  "8 ";
	    for( int j=1; j <= 8; ++j )
	      ofs << meshes[ m ].con.Node( e, j )-1<< " ";
	    ofs << std::endl;
	  }
	
	ofs << "CELL_TYPES " << meshes[ m ].con.Nelem( ) << std::endl;
	for( int e=1; e <= meshes[ m ].con.Nelem( ); ++ e )
	  ofs << "12\n";
	
	ofs << "POINT_DATA " << Program.mask[ m ].size( ) << std::endl;
	ofs << "SCALARS SharedNodes double\nLOOKUP_TABLE default\n";
	for( int i=0; i < Program.mask[ m ].size( ); ++i )
	  {
	    if( Program.mask[ m ][ i ] )
	      ofs << "1\n";
	    else
	      ofs << "0\n";
	  }
	ofs.close( );
      }    
  } 
  /**
   * @brief This method prints the global grid in VTK file format.
   * @note Mostly used for debugging.
   * @param vlist the global vertex list.
   * @param elist the element list.
   */
  void printVtk( std::vector< double > &vlist, std::vector< std::vector< int > > &elist )
  {
    std::ofstream ofs;
    ofs.open( std::string( "global_mesh.vtk" ).c_str( ) );
    
    if( !ofs.is_open( ) )
      {
	std::cerr << "Cannot write global VTK file! " << std::endl;
	std::cerr << "File: " << __FILE__ 						<< std::endl;
	std::cerr << "Line: " << __LINE__ 						<< std::endl;
	assert( false );
      }
    
    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Global Mesh Output" << std::endl;
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << vlist.size( )/3 << " double\n";
    
    for( int i=0; i < vlist.size( )/3; ++i )
      {
	ofs << vlist[ i*3   ] << " ";
	ofs << vlist[ i*3+1 ] << " ";
	ofs << vlist[ i*3+2 ] << "\n";
      }
    
    int maxnodes = 0; /* the maximum number of nodes per element */
    for( int i=0; i < elist.size( ); ++i )
      {
	if( maxnodes < elist[ i ].size( ) )
	  maxnodes = elist[ i ].size( );
      }
    
    ofs << "CELLS " << elist.size( ) << " " << elist.size( )*( maxnodes+1) << "\n";
    for( int i=0; i < elist.size( ); ++i  )
      {
	ofs <<  elist[ i ].size( ) << " ";
	for( int j=0; j < elist[ i ].size( ); ++j )
	  ofs <<  elist[ i ][ j ]<< " ";
	ofs << std::endl;
      }
    
    ofs << "CELL_TYPES " << elist.size( ) << std::endl;
    for( int i=0; i < elist.size( ); ++i )
      {
	
	if( elist[ i ].size( ) == 8 )
	  ofs << "12\n";
	else if( elist[ i ].size( ) == 4 )
	  ofs << "10\n";
	else {
	  std::cerr << "Undefined element type!\n";
	  std::cerr << "File: " << __FILE__ << std::endl;
	  std::cerr << "Line: " << __LINE__ << std::endl;
	}
	
      }
    ofs.close( );
  }  
/**
 * @brief Writes the merged grid and its data in VTK UNSTRUCTURED_GRID file format.
 * @pre Program.globalNodeList.size( ) 	  > 0.
 * @pre Program.globalElementList.size( ) > 0.
 * @pre Program.smdv.size( ) > 0.
 * @pre Program.globalSolution.size( ) == Program.smdv.size( ).
 * @pre Program.solutionignore.size( ) == Program.smdv.size( ).
 * @pre Program.outputfile != ""
 * @post A VTK formatted file is written in Program.outputfile.vtk
 * @note The extension in the outputfile will be appended.
 */
  void writeVtkData( NodalCoordinates &nc,Connectivity &con,const std::string &filename,std::vector<double> &soln)
  {
    
    std::cout << "Writing VTK data...";
    std::cout.flush( );
    
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
    for( int i=1; i < nnodes; ++i )
      {
	ofs << GeoPrim::C3Point(&(nc[i])) << std:endl;
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
	ofs <<  con.Esize(i) << " ";
	for( int j=0; j < con.Esize(i); ++j )
	  ofs <<  con[ i ][ j ]-1 << " ";
	ofs << std::endl;
      }
    
    /* STEP 5: Write the cell types */
    ofs << "CELL_TYPES " << nelem << std::endl;
    for( int i=0; i < nelem; ++i )
      {
	
	if( con.Esize(i) == 8 )
	  ofs << "12\n";
	else if( con.Esize(i) == 4 )
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
	       ofs << soln[ i*nelem + j ] << " ";
	     }
	   ofs << std::endl;
	 }
     }
     ofs.close( );
     std::cout << "[DONE]\n";
  }
}
