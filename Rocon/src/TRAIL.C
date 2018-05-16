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
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <map>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "TRAIL_UnixUtils.H"

//#ifdef _TRAIL_MPI_
#include "mpi.h"
//#else
//typedef int MPI_Comm;
//#endif

#include "GEM.H"
#include "TRAIL.H"

//#ifdef _ROCSTAR_X_
#include "com.h"
//#endif

COM_EXTERN_MODULE( SurfX);
COM_EXTERN_MODULE( SimOUT);
COM_EXTERN_MODULE( SimIN);
COM_EXTERN_MODULE( Simpal);
COM_EXTERN_MODULE( SurfUtil);
COM_EXTERN_MODULE( Rocmop);
COM_EXTERN_MODULE( Rocprop);

void
TRAIL_Debug(GEM_Partition &gp)
{
  int rank = 0;
  //#ifdef _TRAIL_MPI_
  MPI_Comm_rank(gp._comm,&rank);
  //#endif
  if(!rank)
    TRAIL_CreateDirectory("Roctrail");
  //#ifdef _TRAIL_MPI_
  MPI_Barrier(gp._comm);
  //#endif
  gp.debug(true);
  std::ostringstream Ostr;
  Ostr << "Roctrail/Roctrail_debug_" << gp._id;
  std::ofstream *Ouf;
  Ouf = new std::ofstream;
  Ouf->open(Ostr.str().c_str());
  gp._out = Ouf;
}

//  Get a string that encodes the given time. If the string is
//    xx.yyyyyy, it means 0.yyyyyy*10^xx nanoseconds. This format
//    ensures that sorting the strings alphabetically gives the right
//    ordering of time. The string is null-terminated.
std::string
TRAIL_TimeString( double t)
{
  std::ostringstream Ostr;
  Ostr << std::scientific << std::setprecision(5) << t*1e10;
  std::string buf(Ostr.str());
  std::string timestring(buf.substr(9)+"."+buf.substr(0,1)+buf.substr(2,5));
  return(timestring);
}

double
TRAIL_TimeString(const std::string &ts)
{
  std::string xp_s(ts.substr(0,2));
  std::string num_s(ts.substr(3));
  std::istringstream X_inst(ts);
  int xp = 0;
  X_inst >> xp;
  xp = xp - 15;
  double x = pow(10,(double)xp);
  std::istringstream Num_inst(num_s);
  double num = 0;
  Num_inst >> num;
  return(num*=x);

}

void
TRAIL_GetRocstarDumpStrings(const std::string &filename,std::string &wname,std::string &timestring,
			   std::string &rankstring)
{
  std::string fname(filename.substr(0,filename.find_last_of(".")));
  std::string::size_type x = fname.find_last_of("_");
  std::string::size_type len = fname.length() - x;
  rankstring.assign(fname.substr(x+1,len-1));
  fname = fname.substr(0,x);
  x = fname.find_last_of("_");
  len = fname.length() - x;
  timestring.assign(fname.substr(x+1,len-1));
  wname = fname.substr(0,x);
}


//#ifdef _MESH_X_
int
TRAIL_ExtractPanes(const std::string &window_name,
		  const std::string &attribute_name,
		  int attribute_value,
		  std::vector<int> &pane_id)
{
  std::vector<int> pane_ids;
  pane_id.resize(0);
  COM_get_panes(window_name.c_str(),pane_ids);
  std::vector<int>::iterator pi = pane_ids.begin();
  while(pi != pane_ids.end()){
    int *attptr = NULL;
    COM_get_array((window_name+"."+attribute_name).c_str(),*pi,&attptr);
    if(*attptr == attribute_value)
      pane_id.push_back(*pi);
    pi++;
  }
  return 0;
}

int
TRAIL_GetPanelAttribute(const std::string &window_name,
		       const std::string &attribute_name,
		       const std::string &qual_name,
		       int qualval,
		       std::vector<int> &attvec)
{
  std::vector<int> pane_ids;
  attvec.resize(0);
  COM_get_panes(window_name.c_str(),pane_ids);
  std::vector<int>::iterator pi = pane_ids.begin();
  while(pi != pane_ids.end()){
    int *attptr = NULL;
    COM_get_array((window_name+"."+qual_name).c_str(),*pi,&attptr);
    if(*attptr == qualval){
      COM_get_array((window_name+"."+attribute_name).c_str(),*pi,&attptr);
      if(attptr)
	attvec.push_back(*attptr);
    }
    pi++;
  }
  return 0;
}

int
TRAIL_GetPanelAttribute(const std::string &window_name,
		       const std::string &attribute_name,
		       std::vector<int> &attvec)
{
  std::vector<int> pane_ids;
  attvec.resize(0);
  COM_get_panes(window_name.c_str(),pane_ids);
  std::vector<int>::iterator pi = pane_ids.begin();
  while(pi != pane_ids.end()){
    int *attptr = NULL;
    COM_get_array((window_name+"."+attribute_name).c_str(),*pi,&attptr);
    if(attptr)
      attvec.push_back(*attptr);
    pi++;
  }
  return 0;
}

int TRAIL_UniqueAcrossProcs(std::vector<int> &input_data,std::vector<int> &output_data,
			   MPI_Comm communicator)
{
  int nprocs = 1;
  MPI_Comm_size(communicator,&nprocs);
  int nlocal_items = input_data.size();
  std::vector<int> nitems_per_processor(nprocs);
  MPI_Allgather(&nlocal_items,1,MPI_INT,&nitems_per_processor[0],1,MPI_INT,communicator);
 std::vector<int> displacements(nprocs);
 displacements[0] = 0;
  int total_count = nitems_per_processor[0];
  for(int i = 1;i < nprocs;i++){
    total_count += nitems_per_processor[i];
    displacements[i] = nitems_per_processor[i-1]+displacements[i-1];
  }
  std::vector<int> all_items(total_count);
  MPI_Allgatherv(&input_data[0],nlocal_items,MPI_INT,
  		 &all_items[0],&nitems_per_processor[0],&displacements[0],
  		 MPI_INT,communicator);
  std::sort(all_items.begin(),all_items.end());
  std::vector<int>::iterator ui = std::unique(all_items.begin(),all_items.end());
  output_data.resize(ui - all_items.begin());
  std::copy(all_items.begin(),ui,output_data.begin());
  return 0;
}

/// Given an array of ranges:
/// search_extent[nd][2]
/// and one array of array of ranges:
/// partition_extent[number_of_partitions][nd][2]
/// Provide an array of range arrays:
/// neighbor_extent[number_of_neighbors][nd][2]
/// and a processor list:
/// neighbors[number_of_neighbors]
///
/// This is the basic search to find what elements of the
/// extent pool collide with the given search extent.
int
TRAIL_Search_Block_Structured_Pool(std::vector<std::vector<int> > &search_extent,
				  std::vector<std::vector<std::vector<int> > > &extent_pool,
				  std::vector<std::vector<std::vector<int> > > &neighbor_extent,
				  std::vector<int> &neighbors)
{
  size_t pool_depth = extent_pool.size();
  size_t nd = search_extent.size();
  neighbor_extent.resize(0);
  neighbors.resize(0);
  std::vector<std::vector<int> > extent(nd);
  for(size_t i = 0;i < nd;i++)
    extent[i].resize(2);
  for(size_t i = 0;i < pool_depth;i++){
    bool match = true;
    for(size_t j = 0;j < nd;j++){
      if(!((search_extent[j][0] >= extent_pool[i][j][0] &&
	    search_extent[j][0] <= extent_pool[i][j][1]) ||
	   (search_extent[j][1] >= extent_pool[i][j][0] &&
	    search_extent[j][1] <= extent_pool[i][j][1])))
	match = false;
    }
    if(match){ // then partition has some searched for mesh points
      neighbors.push_back(i);
      for(size_t j = 0;j < nd;j++){
	extent[j][0] = std::max(search_extent[j][0],extent_pool[i][j][0]);
	extent[j][1] = std::min(search_extent[j][1],extent_pool[i][j][1]);
      }
      neighbor_extent.push_back(extent);
    }
  }
  return(0);
}


/// Given two arrays of ranges:
/// global_extent[nd][2]
/// local_extent[nd][2]
/// an extent pool:
/// extent_pool[npool][nd][2]
/// provide the neighboring extent pool:
/// neighbor_extent[nneighbors][nd][2]
/// and a neighbor id list:
/// neighbors[nneighbors]
///
/// This will take the local extent and
/// create the full list of neighbors and
/// the range of remote indices.
int
TRAIL_Get_Block_Structured_Neighbors(std::vector<std::vector<int> > &local_extent,
				    std::vector<std::vector<int> > &global_extent,
				    std::vector<std::vector<std::vector<int> > > &extent_pool,
				    std::vector<std::vector<int> > &ghost_extent,
				    std::vector<std::vector<std::vector<int> > > &neighbor_extent,
				    std::vector<int> &neighbors)
{
  size_t pool_depth = extent_pool.size();
  size_t nd = local_extent.size();
  std::vector<std::vector<std::vector<int> > > pool_extents(pool_depth);
  neighbors.resize(0);
  neighbor_extent.resize(0);
  for(size_t i = 0;i < pool_depth;i++){
    pool_extents[i].resize(nd);
    for(size_t j = 0;j < nd;j++)
      pool_extents[i][j].resize(2,0);
  }
  ghost_extent = local_extent;
  for(size_t i = 0;i < nd;i++){
    // Left neighbors
    if(local_extent[i][0] > global_extent[i][0])
      ghost_extent[i][0]  = local_extent[i][0] - 1;
    if(local_extent[i][1] < global_extent[i][1])
      ghost_extent[i][1]  = local_extent[i][1] + 1;
  }
  for(size_t i = 0;i < nd;i++){
    // Left neighbors
    if(local_extent[i][0] > global_extent[i][0]){
      std::vector<std::vector<int> > search_extent(ghost_extent);
      search_extent[i][1] = search_extent[i][0];
      //      search_extent[i][0] = search_extent[i][1] = ghost_extent[i][0];
      std::vector<std::vector<std::vector<int> > > directional_extent;
      std::vector<int> directional_neighbors;
      TRAIL_Search_Block_Structured_Pool(search_extent,extent_pool,
					directional_extent,directional_neighbors);
      std::vector<int>::iterator dni = directional_neighbors.begin();
      size_t ncount = 0;
      while(dni != directional_neighbors.end()){
	int neighbor_index = *dni++;
	for(size_t j = 0;j < nd;j++){
	  if(pool_extents[neighbor_index][j][0] > 0) // This is done to avoid clash with previous setting
	    pool_extents[neighbor_index][j][0] = std::min(directional_extent[ncount][j][0],
							  pool_extents[neighbor_index][j][0]);
	  else
	    pool_extents[neighbor_index][j][0] = directional_extent[ncount][j][0];
	  if(pool_extents[neighbor_index][j][1] > 0)
	    pool_extents[neighbor_index][j][1] = std::max(directional_extent[ncount][j][1],
							  pool_extents[neighbor_index][j][1]);
	  else
	    pool_extents[neighbor_index][j][1] = directional_extent[ncount][j][1];

	}
	ncount++;
	//	neighbors.push_back(neighbor_index);
      }
    }
    // Right neighbors
    if(local_extent[i][1] < global_extent[i][1]){
      std::vector<std::vector<int> > search_extent(ghost_extent);
      search_extent[i][0] = search_extent[i][1] = ghost_extent[i][1];
      std::vector<std::vector<std::vector<int> > > directional_extent;
      std::vector<int> directional_neighbors;
      TRAIL_Search_Block_Structured_Pool(search_extent,extent_pool,
					directional_extent,directional_neighbors);
      std::vector<int>::iterator dni = directional_neighbors.begin();
      size_t ncount = 0;
      while(dni != directional_neighbors.end()){
	int neighbor_index = *dni++;
	for(size_t j = 0;j < nd;j++){
	  if(pool_extents[neighbor_index][j][0] > 0)
	    pool_extents[neighbor_index][j][0] = std::min(directional_extent[ncount][j][0],
							  pool_extents[neighbor_index][j][0]);
	  else
	    pool_extents[neighbor_index][j][0] = directional_extent[ncount][j][0];
	  if(pool_extents[neighbor_index][j][1] > 0)
	    pool_extents[neighbor_index][j][1] = std::max(directional_extent[ncount][j][1],
							  pool_extents[neighbor_index][j][1]);
	  else
	    pool_extents[neighbor_index][j][1] = directional_extent[ncount][j][1];

	}
	ncount++;
	//	neighbors.push_back(neighbor_index);
      }
    }
  }
  for(size_t i = 0;i < pool_depth;i++){
    if(pool_extents[i][0][0] > 0){
      neighbors.push_back(i);
      neighbor_extent.push_back(pool_extents[i]);
    }
  }
  return(0);
}

template<typename T>
void TRAIL_Copy2Attribute(const std::string &aname,const std::vector<T> &container,int pane_id,int asize = 1)
{

  COM_set_size(aname.c_str(),pane_id,asize);
  COM_resize_array(aname.c_str(),pane_id);
  T *leptr;
  COM_get_array(aname.c_str(),pane_id,&leptr);
  typename std::vector<T>::const_iterator ci = container.begin();
  while(ci != container.end())
    *leptr++ = *ci++;
}

template<typename T>
void TRAIL_SetAttribute(const std::string &aname,int pane_id,T &value)
{

  COM_set_size(aname.c_str(),pane_id,1);
  COM_resize_array(aname.c_str(),pane_id);
  T *leptr;
  COM_get_array(aname.c_str(),pane_id,&leptr);
  *leptr = value;
}

template<typename T>
void TRAIL_Copy2Attribute(const std::string &aname,const std::vector<std::vector<T> > &container,int pane_id)
{
  unsigned int asize = container.size();
  COM_set_size(aname.c_str(),pane_id,asize);
  COM_resize_array(aname.c_str(),pane_id);
  T *leptr;
  COM_get_array(aname.c_str(),pane_id,&leptr);
  typename std::vector<std::vector<T> >::const_iterator ci = container.begin();
  while(ci != container.end()){
    typename std::vector<T>::const_iterator vi = ci->begin();
    while(vi != ci->end())
      *leptr++ = *vi++;
    ci++;
  }
}

/// Creates a window from a Mesh object. (copies data)
int TRAIL_SurfaceMesh2Window(const std::string &wname,int pane_id,Mesh::NodalCoordinates &coords,
			    Mesh::Connectivity &conn)
{
  COM_new_window(wname);
  unsigned int number_of_nodes = coords.Size();
  COM_set_size((wname+".nc").c_str(),pane_id,number_of_nodes);
  COM_resize_array((wname+".nc").c_str(),pane_id);
  double *nc = NULL;
  COM_get_array((wname+".nc").c_str(),pane_id,&nc);
  std::memcpy(nc,coords[1],number_of_nodes*3*sizeof(double));
  Mesh::Connectivity tris;
  Mesh::Connectivity quads;
  Mesh::Connectivity::iterator ci = conn.begin();
  bool known_element_type = false;
  while(ci != conn.end()){
    if(ci->size() == 3)
      tris.AddElement(*ci++);
    else if (ci->size() == 4)
      quads.AddElement(*ci++);
    else
      assert(known_element_type);
  }
  unsigned int number_of_tris  = tris.Nelem();
  unsigned int number_of_quads = quads.Nelem();
  int *con_array = NULL;
  if(number_of_tris > 0){
    COM_set_size((wname+".:t3:real").c_str(),pane_id,number_of_tris);
    COM_resize_array((wname+".:t3:real").c_str(),pane_id);
    COM_get_array((wname+".:t3:real").c_str(),pane_id,&con_array);
    ci = tris.begin();
    while(ci != tris.end()){
      unsigned int element_id = ci - tris.begin() + 1;
      unsigned int element_index = element_id - 1;
      con_array[3*element_index]   = (*ci)[0];
      con_array[3*element_index+1] = (*ci)[1];
      con_array[3*element_index+2] = (*ci)[2];
      ci++;
    }
  }
  if(number_of_quads > 0){
    COM_set_size((wname+".:q4:").c_str(),pane_id,number_of_quads);
    COM_resize_array((wname+".:q4:").c_str(),pane_id);
    COM_get_array((wname+".:q4:").c_str(),pane_id,&con_array);
    ci = quads.begin();
    while(ci != quads.end()){
      unsigned int element_id = ci - quads.begin() + 1;
      unsigned int element_index = element_id - 1;
      con_array[4*element_index]   = (*ci)[0];
      con_array[4*element_index+1] = (*ci)[1];
      con_array[4*element_index+2] = (*ci)[2];
      con_array[4*element_index+3] = (*ci)[3];
      ci++;
    }
  }
  COM_window_init_done(wname.c_str());
  return(0);
}

/// Creates a window from a Mesh object. (copies data)
int TRAIL_UnstructuredMesh2Pane(const std::string &wname,int pane_id,
				 Mesh::UnstructuredMesh &mesh,
				 SolnMetaData &smdv,
				 std::vector<std::vector<double> > &soln_data,
				 int verblevel)
{


/*
 * Note: This throws an error that the window is duplicated. Looking
 * at Roccom_base::get_status it looks like
 * Now, the window should be created prior to calling this function.
 * George Zagaris (gzagaris@illinois.edu)
 */
//  if(COM_get_status(wname.c_str(),0)<0)
//    COM_new_window(wname);


  unsigned int number_of_nodes = mesh.nc.Size();
//  std::cout << "Number of nodes: " << mesh.nc.Size() << std::endl;
  COM_set_size((wname+".nc").c_str(),pane_id,number_of_nodes);
  COM_resize_array((wname+".nc").c_str(),pane_id);
  double *nc = NULL;
  COM_get_array((wname+".nc").c_str(),pane_id,&nc);
  std::memcpy(nc,mesh.nc[1],number_of_nodes*3*sizeof(double));
  Mesh::Connectivity tets;
  Mesh::Connectivity bricks;
  Mesh::Connectivity::iterator ci = mesh.con.begin();
//  std::cout << "Number of elements: " << mesh.con.Nelem( ) << std::endl;
  bool known_element_type = false;
  while(ci != mesh.con.end()){
    if(ci->size() == 4)
      tets.AddElement(*ci++);
    else if (ci->size() == 8)
      bricks.AddElement(*ci++);
    else
      assert(known_element_type);
  }
  unsigned int number_of_tets  = tets.Nelem();
  unsigned int number_of_bricks = bricks.Nelem();
//  std::cout << "Number of tets: " << number_of_tets << std::endl;
//  std::cout << "Number of bricks: " << number_of_bricks << std::endl;
  int *con_array = NULL;
  if(number_of_tets > 0){
    COM_set_size((wname+".:T4:real").c_str(),pane_id,number_of_tets);
    COM_resize_array((wname+".:T4:real").c_str(),pane_id);
    COM_get_array((wname+".:T4:real").c_str(),pane_id,&con_array);
    ci = tets.begin();
    while(ci != tets.end()){
      unsigned int element_id = ci - tets.begin() + 1;
      unsigned int element_index = element_id - 1;
      con_array[4*element_index]   = (*ci)[0];
      con_array[4*element_index+1] = (*ci)[1];
      con_array[4*element_index+2] = (*ci)[2];
      con_array[4*element_index+3] = (*ci)[3];
      ci++;
    }
  }
  if(number_of_bricks > 0){
    COM_set_size((wname+".:B8:real").c_str(),pane_id,number_of_bricks);
    COM_resize_array((wname+".:B8:real").c_str(),pane_id);
    COM_get_array((wname+".:B8:real").c_str(),pane_id,&con_array);
    ci = bricks.begin();
    while(ci != bricks.end()){
      unsigned int element_id = ci - bricks.begin() + 1;
      unsigned int element_index = element_id - 1;
      con_array[8*element_index]   = (*ci)[0];
      con_array[8*element_index+1] = (*ci)[1];
      con_array[8*element_index+2] = (*ci)[2];
      con_array[8*element_index+3] = (*ci)[3];
      con_array[8*element_index+4] = (*ci)[4];
      con_array[8*element_index+5] = (*ci)[5];
      con_array[8*element_index+6] = (*ci)[6];
      con_array[8*element_index+7] = (*ci)[7];
      ci++;
    }
  }

 // Window is created and initialized outside
 // COM_window_init_done(wname.c_str());

  return(0);
}

/// Adds ghost zones for block structured meshes to close gaps
/// in the interface surface mesh.
/// Needs the following data for each patch/pane:
/// Block ID: block_id
/// Patch ID: patch_id
/// Local Patch Extent: local_patch_extent [istart,iend,jstart,jend,kstart,kend] (global indices for local patch extent)
/// Global Patch Extent: global_patch_extent [istart,iend,jstart,jend,kstart,kend]
int
TRAIL_FD2FE_WinCreate(const std::string &wname,const std::string &outwname,std::ostream *ouf)
{
  // Get the window communicator
  MPI_Comm communicator = MPI_COMM_NULL;
  COM_get_communicator(wname.c_str(),&communicator);
  IRAD::Comm::CommunicatorObject BaseComm(communicator);
  int nprocs = BaseComm.Size();
  int rank   = BaseComm.Rank();
  if(ouf && !rank)
    *ouf << "TRAIL_AddBlockStructuredGhostZone::Nprocs = " << nprocs << std::endl;
  std::vector<int> pane_ids;
  COM_get_panes(wname.c_str(),pane_ids);
  if(ouf && !rank)
    *ouf << "Number of panes: " << pane_ids.size() << std::endl;
  // Form a list of unique block id's across all processors
  std::vector<int> local_block_ids;
  TRAIL_GetPanelAttribute(wname,"block_id",local_block_ids);
  if(ouf && !rank){
    *ouf << "Local Block_ids: ";
    IRAD::Util::DumpContents(*ouf,local_block_ids," ");
    *ouf << std::endl;
  }
  std::vector<int> global_blocks;
  TRAIL_UniqueAcrossProcs(local_block_ids,global_blocks,BaseComm.World());
  if(ouf && !rank){
    *ouf << "Global Block_ids: ";
    IRAD::Util::DumpContents(*ouf,local_block_ids," ");
    *ouf << std::endl;
  }
  // Now all_block_ids is an array of all the unique block ids across all procs
  // For each block
  std::vector<int>::iterator bii = global_blocks.begin();
  std::string wname2 = outwname; // wname+"_uns";
  bool window_created = false;
  if(!window_created){
    COM_new_window(wname2.c_str());
    COM_new_dataitem((wname2+".local_extent").c_str(), 'p',COM_INTEGER,6,"");
    COM_new_dataitem((wname2+".global_extent").c_str(),'p',COM_INTEGER,6,"");
    COM_new_dataitem((wname2+".send_extent").c_str(),  'p',COM_INTEGER,6,"");
    COM_new_dataitem((wname2+".send_panes").c_str(),   'p',COM_INTEGER,1,"");
    COM_new_dataitem((wname2+".recv_panes").c_str(),   'p',COM_INTEGER,1,"");
    COM_new_dataitem((wname2+".recv_extent").c_str(),  'p',COM_INTEGER,6,"");
    //    COM_new_dataitem((wname2+".sister_pane").c_str(),  'p',COM_INTEGER,1,"");
    COM_new_dataitem((wname2+".block_id").c_str(),     'p',COM_INTEGER,1,"");
    COM_new_dataitem((wname2+".patch_id").c_str(),     'p',COM_INTEGER,1,"");
    window_created = true;
    COM_new_dataitem((wname+".send_extent").c_str(),   'p',COM_INTEGER,6,"");
    COM_new_dataitem((wname+".send_panes").c_str(),    'p',COM_INTEGER,1,"");
    COM_new_dataitem((wname+".recv_extent").c_str(),   'p',COM_INTEGER,6,"");
    COM_new_dataitem((wname+".recv_panes").c_str(),    'p',COM_INTEGER,1,"");
    //   COM_new_dataitem((wname+".sister_pane").c_str(),   'p',COM_INTEGER,1,"");
  }
  while(bii != global_blocks.end()){
    int block_id = *bii++;
    std::vector<int>::iterator fi = std::find(local_block_ids.begin(),local_block_ids.end(),block_id);
    int block_color = 0;
    if(fi != local_block_ids.end()) // Then this processor has data for the given block
      block_color = 1;
    // Split the communicator into haves and have nots for this block
    if(ouf && !rank)
      *ouf << "Splitting communicator for blocks." << std::endl;
    IRAD::Comm::CommunicatorObject BlockComm;
    BaseComm.Split(block_color,rank,BlockComm);
    if(block_color == 1){ // all the guys with data on this block
      int block_nproc = BlockComm.Size();
      int block_rank  = BlockComm.Rank();
      // Form a list of unique patch id's across all processors
      std::vector<int> local_patch_ids;
      std::vector<int> panes;
      TRAIL_ExtractPanes(wname,"block_id",block_id,panes); // get pane id's of panes with block id = block_id
      std::vector<int>::iterator pi = panes.begin();
      while(pi != panes.end()){
	int *patch_id;
	COM_get_array((wname+".patch_id").c_str(),*pi,&patch_id);
	if(patch_id)
	  local_patch_ids.push_back(*patch_id);
	pi++;
      }
      std::vector<int> all_patches(block_nproc);
      TRAIL_UniqueAcrossProcs(local_patch_ids,all_patches,BlockComm.World());
      // Now all_patches is an array of all the unique patch id across all procs having the current block
      std::vector<int>::iterator pii = all_patches.begin();
      while(pii != all_patches.end()){ // For each patch
	int patch_id = *pii++;
	// Determine if local processor owns part of the patch
	// Split the communicator into haves and have nots
	int patch_color = 0;
	IRAD::Comm::CommunicatorObject PatchComm;
	std::vector<int>::iterator fp = std::find(local_patch_ids.begin(),local_patch_ids.end(),patch_id);
	if(fp != local_patch_ids.end())
	  patch_color = 1;
	if(ouf && !rank)
	  *ouf << "Splitting communicator for patches." << std::endl;
	BlockComm.Split(patch_color,block_rank,PatchComm);
	if(patch_color == 1) { // all of us that have data on the given block/patch
	  int patch_nproc = PatchComm.Size();
	  int patch_rank  = PatchComm.Rank();
	  std::vector<int> patch_pane;
	  TRAIL_ExtractPanes(wname,"patch_id",patch_id,patch_pane); // get the pane for this patch (hopefully only 1)
	  int *global_patch_extent = NULL;
	  COM_get_array((wname+".global_extent").c_str(),patch_pane[0],&global_patch_extent);
	  if(!global_patch_extent){
	    std::cerr << "ERROR: Window " << wname << " has no global_extent attribute." << std::endl;
	    return 1;
	  }
	  int *local_patch_extent_ptr = NULL;
	  COM_get_array((wname+".local_extent").c_str(),patch_pane[0],&local_patch_extent_ptr);
	  if(!local_patch_extent_ptr){
	    std::cerr << "ERROR: Window " << wname << " has no local_extent attribute." << std::endl;
	    return 1;
	  }
	  std::vector<int> local_patch_extent;
	  for(unsigned int aa = 0;aa < 6;aa++)
	    local_patch_extent.push_back(local_patch_extent_ptr[aa]);
	  if(ouf && !rank)
	    *ouf << "Getting Local Coordinates." << std::endl;
	  double *LocalCoordinates = NULL;
	  COM_get_array((wname+".nc").c_str(),patch_pane[0],&LocalCoordinates);
	  // Determine the total number of ghost nodes we need and augment the extended local patch extent
	  // and allocate an array for holding ALL the nodal coordinates
	  // Communicate the local patch extents
	  std::vector<int> all_patch_extents(6*patch_nproc);
	  PatchComm.AllGather(local_patch_extent,all_patch_extents,6,6);
	  std::vector<int> all_pane_ids(patch_nproc);
	  PatchComm.AllGather(patch_pane[0],all_pane_ids);
	  if(ouf && !rank)
	    *ouf << "All patch extents communicated." << std::endl;
	  std::vector<std::vector<std::vector<int> > > extent_pool(patch_nproc);
	  std::vector<std::vector<int> > ghost_extent;
	  std::vector<std::vector<std::vector<int> > > neighbor_extent;
	  std::vector<int> neighbors;
	  std::vector<std::vector<int> > local_extent(3);
	  std::vector<std::vector<int> > global_extent(3);
	  for(int dd = 0;dd < 3;dd++){
	    local_extent[dd].push_back(local_patch_extent[dd*2]);
	    local_extent[dd].push_back(local_patch_extent[dd*2+1]);
	    global_extent[dd].push_back(global_patch_extent[dd*2]);
	    global_extent[dd].push_back(global_patch_extent[dd*2+1]);
	  }
	  for(int pp = 0;pp < patch_nproc;pp++){
	    extent_pool[pp].resize(3);
	    for(int dd = 0;dd < 3;dd++){
	      extent_pool[pp][dd].push_back(all_patch_extents[pp*6+dd*2]);
	      extent_pool[pp][dd].push_back(all_patch_extents[pp*6+dd*2+1]);
	    }
	  }
	  if(ouf && !rank)
	    *ouf << "Calling TRAIL_Get_Block_Structured_Neighbors." << std::endl;
	  TRAIL_Get_Block_Structured_Neighbors(local_extent,
					      global_extent,
					      extent_pool,
					      ghost_extent,
					      neighbor_extent,
					      neighbors);

	  if(ouf && !rank)
	    *ouf << "Calling TRAIL_Get_Block_Structured_Neighbors done." << std::endl;
	  // Loop thru neighbors and post receives for the Nodal Coordinates and
	  // the local extents requested from each neighbor.
	  std::vector<int>::iterator ni = neighbors.begin();
	  int number_of_neighbors = neighbors.size();
	  std::vector<std::vector<int> > RemoteNeighborExtent(number_of_neighbors);
	  std::vector<std::vector<double> > RemoteNodalCoordinates(number_of_neighbors);
	  std::vector<std::vector<double> > SendCoordinates(number_of_neighbors);
	  std::vector<std::vector<int> > FlattenedNeighborExtents(number_of_neighbors);
	  std::vector<int> RecvMsgID(number_of_neighbors);
	  PatchComm.Barrier();
	  if(ouf && !rank)
	    *ouf << "All processors Doing ComLoop 1" << std::endl;
	  for(int nc = 0;nc < number_of_neighbors;nc++){
	    int nid = *ni++;
	    if(nid != patch_rank){
	      RemoteNeighborExtent[nc].resize(6);
	      RecvMsgID[nc] = PatchComm.ARecv(RemoteNeighborExtent[nc],nid,1);
	      Mesh::BSExtent<int> bs_extent(neighbor_extent[nc]);
	      RemoteNodalCoordinates[nc].resize(3*bs_extent.NNodes());
	      bs_extent.Flatten(FlattenedNeighborExtents[nc]);
	      PatchComm.ASend(FlattenedNeighborExtents[nc],nid,1);
	      if(ouf){
		*ouf << "Sending the extents I request from " << nid << ": ";
		IRAD::Util::DumpContents(*ouf,FlattenedNeighborExtents[nc]," ");
		*ouf << std::endl;
	      }
	      PatchComm.ARecv(RemoteNodalCoordinates[nc],nid,2);
	    }
	  }
	  // Loop thru the neighbors and send each one the indices we need and
	  // the nodal coordinates they need.
	  ni = neighbors.begin();
	  PatchComm.Barrier();
	  if(ouf && !rank)
	    *ouf << "Doing ComLoop 2." << std::endl;
	  for(int nc = 0;nc < number_of_neighbors;nc++){
	    int nid = *ni++;
	    if(nid != patch_rank){
	      PatchComm.WaitRecv(RecvMsgID[nc]);
	      Mesh::BSExtent<int> bs_extent(RemoteNeighborExtent[nc]);
	      int nnodes = bs_extent.NNodes();
	      if(ouf && !rank)
		*ouf << "Sending " << nnodes << " nodal coordinates to remote proc." << std::endl;
	      SendCoordinates[nc].resize(0);
	      std::vector<int> flat_indices;
	      Mesh::BSExtent<int> LocalExtent(local_extent);
	      LocalExtent.GetFlatIndices(bs_extent,flat_indices);
	      std::vector<int>::iterator gii = flat_indices.begin();
	      while(gii != flat_indices.end()){
		int nodeind = (*gii++ - 1)*3;
		SendCoordinates[nc].push_back(LocalCoordinates[nodeind]);
		SendCoordinates[nc].push_back(LocalCoordinates[nodeind + 1]);
		SendCoordinates[nc].push_back(LocalCoordinates[nodeind + 2]);
	      }
	      if(ouf && !rank)
		*ouf << "SendCoordinates.size == " << SendCoordinates[nc].size() << std::endl;
	      PatchComm.ASend(SendCoordinates[nc],nid,2);
	    }
	  }
	  PatchComm.Barrier();
	  if(ouf && !rank)
	    *ouf << "Waiting for all messages." << std::endl;
	  PatchComm.WaitAll();
	  PatchComm.ClearRequests();
	  if(ouf && !rank)
	    *ouf << "All messages received, requests cleared." << std::endl;
	  // Now all the coordinates we need to build the ghost zone are in RemoteNodalCoordinates array
	  // Here we register the new surface with Roccom.   This will be the surface collecting points
	  // greedily on the right and not at all on the left.
	  // Create and size the unstructured extent array
	  Mesh::BSExtent<int> uns_extent(local_extent);
	  Mesh::BSExtent<int>::iterator uei = uns_extent.begin();
	  Mesh::BSExtent<int>::iterator gei = ghost_extent.begin();
	  while(uei != uns_extent.end())
	    (*uei++)[1] = (*gei++)[1];
	  uns_extent.Sync();
// 	  if(ouf && !rank){
// 	    std::vector<int> tempout;
// 	    uns_extent.Flatten(tempout);
// 	    *ouf << "Unstructured extent: ";
// 	    IRAD::Util::DumpContents(*ouf,tempout," ");
// 	    *ouf << std::endl;
// 	  }
	  // Register it as the local_extent
	  std::vector<int> flatunsextent;
	  uns_extent.Flatten(flatunsextent);
	  TRAIL_Copy2Attribute((wname2+".local_extent"),flatunsextent,patch_pane[0]);
	  // Do the global
	  Mesh::BSExtent<int> GlobalExtent(global_extent);
	  std::vector<int> flatglobextent;
	  GlobalExtent.Flatten(flatglobextent);
	  TRAIL_Copy2Attribute((wname2+".global_extent"),flatunsextent,patch_pane[0]);
	  TRAIL_SetAttribute((wname2+".block_id"),patch_pane[0],block_id);
	  TRAIL_SetAttribute((wname2+".patch_id"),patch_pane[0],patch_id);
	  if(ouf && !rank)
	    *ouf << "Patch ID = " << patch_id << std::endl;
	  //	  TRAIL_SetAttribute((wname2+".sister_pane")),patch_pane[0],pane_id);
	  Mesh::Connectivity conn;
	  uns_extent.CreateUnstructuredMesh(conn);
// 	  if(ouf && !rank)
// 	    *ouf << "Root proc created " << conn.Nelem() << " elements." << std::endl
// 		 << "Connectivity:" << std::endl
// 		 << conn << std::endl;
	  std::vector<int> flatcon;
	  unsigned int nnodes_new = uns_extent.NNodes();
	  std::vector<double> NewCoords(3*nnodes_new);
	  conn.Flatten(flatcon);
	  COM_set_size((wname2+".nc").c_str(),patch_pane[0],nnodes_new);
	  COM_resize_array((wname2+".nc").c_str(),patch_pane[0]);
	  COM_set_size((wname2+".:q4:").c_str(),patch_pane[0],conn.Nelem());
	  COM_resize_array((wname2+".:q4:").c_str(),patch_pane[0]);
	  double *coords;
	  COM_get_array((wname2+".nc").c_str(),patch_pane[0],&coords);
	  int *conndat;
	  COM_get_array((wname2+".:q4:").c_str(),patch_pane[0],&conndat);
	  memcpy(conndat,&flatcon[0],sizeof(int)*4*conn.Nelem());
	  std::vector<int> localind;
	  Mesh::BSExtent<int> LocalExtent(local_extent);
	  unsigned int nlocalnodes = LocalExtent.NNodes();
	  uns_extent.GetFlatIndices(LocalExtent,localind);
	  for(unsigned int a = 0;a < nlocalnodes;a++){
	    coords[(localind[a]-1)*3]   = LocalCoordinates[a*3];
	    coords[(localind[a]-1)*3+1] = LocalCoordinates[a*3+1];
	    coords[(localind[a]-1)*3+2] = LocalCoordinates[a*3+2];
	  }
	  std::vector<int> send_panes;
	  std::vector<int> recv_panes;
	  std::vector< std::vector<int> > flat_recv_extents;
	  std::vector< std::vector<int> > flat_send_extents;
	  for(int nc = 0;nc < number_of_neighbors;nc++){
	    int nid = neighbors[nc];
	    if(nid != patch_rank){
	      if(FlattenedNeighborExtents[nc][1] == LocalExtent[0][1]+1 ||
		 FlattenedNeighborExtents[nc][3] == LocalExtent[1][1]+1 ||
		 FlattenedNeighborExtents[nc][5] == LocalExtent[2][1]+1){
		// Then it has coordinates on the right that I need
		// First, trim the flattened neighbor extent so that it doesn't
		// include any indices out of the local range
		Mesh::BSExtent<int> NeighborExtent(FlattenedNeighborExtents[nc]);
		Mesh::BSExtent<int> TrimmedNeighborExtent(NeighborExtent);
		TrimmedNeighborExtent[0][0] = std::max(TrimmedNeighborExtent[0][0],uns_extent[0][0]);
		TrimmedNeighborExtent[1][0] = std::max(TrimmedNeighborExtent[1][0],uns_extent[1][0]);
		TrimmedNeighborExtent[2][0] = std::max(TrimmedNeighborExtent[2][0],uns_extent[2][0]);
		std::vector<int> flatrecvextent;
		TrimmedNeighborExtent.Flatten(flatrecvextent);
		flat_recv_extents.push_back(flatrecvextent);
		recv_panes.push_back(all_pane_ids[nid]);
		// Now it's trimmed, and we want to copy the relavent coordinates from
		// the original untrimmed received coordinates array to the appropriate
		// spot in the new coordinates array.  To do so, we need the flat indices
		// of the TrimmedNeighborExtent wrt both the source and destination.
		std::vector<int> source_indices;
		std::vector<int> dest_indices;
		NeighborExtent.GetFlatIndices(TrimmedNeighborExtent,source_indices);
		uns_extent.GetFlatIndices(TrimmedNeighborExtent,dest_indices);
		unsigned int ncopynodes = source_indices.size();
		for(unsigned int a = 0;a < ncopynodes;a++){
		  unsigned int sind = (source_indices[a]-1)*3;
		  unsigned int dind = (dest_indices[a]-1)*3;
		  coords[dind]   = RemoteNodalCoordinates[nc][sind];
		  coords[dind+1] = RemoteNodalCoordinates[nc][sind+1];
		  coords[dind+2] = RemoteNodalCoordinates[nc][sind+2];
		}
	      }
	      if(FlattenedNeighborExtents[nc][1] == LocalExtent[0][0]-1 ||
		 FlattenedNeighborExtents[nc][3] == LocalExtent[1][0]-1 ||
		 FlattenedNeighborExtents[nc][5] == LocalExtent[2][0]-1){
		if(ouf)
		  *ouf << "Yeah, I send left to pane " << all_pane_ids[nid] << std::endl;
		// Then this processor is on the left
		// First, trim the flattened neighbor extent so that it doesn't
		// include any indices out of the local range
		send_panes.push_back(all_pane_ids[nid]);
		Mesh::BSExtent<int> NeighborExtent(RemoteNeighborExtent[nc]); // extent requested by remote proc
		//		Mesh::BSExtent<int> TrimmedNeighborExtent(NeighborExtent);
		//		TrimmedNeighborExtent[0][0] = std::min(TrimmedNeighborExtent[0][1],uns_extent[0][1]);
		//		TrimmedNeighborExtent[1][0] = std::min(TrimmedNeighborExtent[1][1],uns_extent[1][1]);
		//		TrimmedNeighborExtent[2][0] = std::min(TrimmedNeighborExtent[2][1],uns_extent[2][1]);
		std::vector<int> flatsendextent;
		NeighborExtent.Flatten(flatsendextent);
		if(ouf){
		  *ouf << "Send extents: ";
		  IRAD::Util::DumpContents(*ouf,flatsendextent," ");
		  *ouf << std::endl;
		}
		flat_send_extents.push_back(flatsendextent);
	      }
	    }
	  }
	  int nrecv = recv_panes.size();
	  int nsend = send_panes.size();
	  TRAIL_Copy2Attribute((wname+".send_panes"),send_panes,patch_pane[0],nsend);
	  TRAIL_Copy2Attribute((wname2+".send_panes"),send_panes,patch_pane[0],nsend);
	  TRAIL_Copy2Attribute((wname+".recv_panes"),recv_panes,patch_pane[0],nrecv);
	  TRAIL_Copy2Attribute((wname2+".recv_panes"),recv_panes,patch_pane[0],nrecv);
	  TRAIL_Copy2Attribute((wname+".send_extent"),flat_send_extents,patch_pane[0]);
	  TRAIL_Copy2Attribute((wname2+".recv_extent"),flat_recv_extents,patch_pane[0]);
	  //	  COM_window_init_done(wname.c_str());
	  //	  COM_window_init_done(wname2.c_str());
	  if(ouf && !rank){
	    *ouf << "Number of panes I send to: " << nsend << std::endl
		 << "Number of panes I recv from: " << nrecv << std::endl
		 << "Send panes: ";
	    IRAD::Util::DumpContents(*ouf,send_panes," ");
	    *ouf << std::endl
		 << "Recv panes: ";
	    IRAD::Util::DumpContents(*ouf,recv_panes," ");
	    *ouf << std::endl
		 << "Send Extents: " << std::endl;
	    for(int o = 0;o < nsend;o++){
	      IRAD::Util::DumpContents(*ouf,flat_send_extents[o]," ");
	      *ouf << std::endl;
	    }
	    *ouf << "Recv Extents:" << std::endl;
	    for(int o = 0;o < nrecv;o++){
	      IRAD::Util::DumpContents(*ouf,flat_recv_extents[o]," ");
	      *ouf << std::endl;
	    }
	  }
	}
      }
      // MPI_Free_Comm (if needed, happens automatically when CommunicatorObject is destroyed)
    }
    // MPI_Free_Comm (if needed, happens automatically when CommunicatorObject is destroyed)
  }
  BaseComm.Barrier();
  return(0);
}

template<typename T>
int TRAIL_Att2Vec(const std::string &att,int pane_id,std::vector<T> &dest)
{
  T *ptr = NULL;
  dest.resize(0);
  int ncomp = 0;
  std::string unit;
  int insize = 0;
  COM_Type type;
  char loc;
  COM_get_dataitem(att.c_str(),&loc,&type,&ncomp,&unit);
  COM_get_size(att.c_str(),pane_id,&insize);
  COM_get_array(att.c_str(),pane_id,&ptr);
  if(!ptr)
    return 1;
  int count = 0;
  int total = ncomp*insize;
  while(count < total)
    dest.push_back(ptr[count++]);
  return 0;
}

int
TRAIL_FE2FD_Transfer(const std::string &fewin,const std::string &fdwin,
		    const std::string &attlist,MPI_Comm communicator,
		    std::ostream *ouf)
{
  std::vector<std::string> atts;
  std::istringstream Istr(attlist);
  std::string attstr;
  if(ouf)
    *ouf << "FE2FD_Transfer: Enter" << std::endl;
  MPI_Comm oldcomm = COM_get_default_communicator();
  COM_set_default_communicator(communicator);
  while(Istr >> attstr)
    atts.push_back(attstr);
  // Get the window communicator
  //  MPI_Comm communicator = MPI_COMM_NULL;
  //  COM_get_communicator(fdwin.c_str(),&communicator);
  IRAD::Comm::CommunicatorObject BaseComm(communicator);
  int rank   = BaseComm.Rank();
  std::vector<std::string>::iterator atti = atts.begin();
  if(ouf)
    *ouf << "FE2FD: Reducing attributes over shared nodes." << std::endl;
  while(atti != atts.end()){
    std::string attstr = fewin+"."+*atti++;
    if(ouf)
      *ouf << "FE2FD: Reducing " << attstr << std::endl;
    int att_hndl = COM_get_dataitem_handle(attstr);
    int att_color = 0;
    if(att_hndl >= 0)
      att_color = 1;
    IRAD::Comm::CommunicatorObject AttComm;
    BaseComm.Split(att_color,rank,AttComm);
    if(att_color == 1){
      MPI_Comm ocomm = COM_get_default_communicator();
      COM_set_default_communicator(AttComm.GetCommunicator());
      COM_LOAD_MODULE_STATIC_DYNAMIC(SurfMap,"MAP");
      int mean_over_shared_nodes = COM_get_function_handle("MAP.reduce_average_on_shared_nodes");
      COM_call_function(mean_over_shared_nodes,&att_hndl);
      COM_UNLOAD_MODULE_STATIC_DYNAMIC(SurfMap,"MAP");
      COM_set_default_communicator(ocomm);
      if(ouf)
	*ouf << "FE2FD: Reduced " << attstr << std::endl;
    }
    else if (ouf)
      *ouf << "FE2FD: " << attstr << " did not exist." << std::endl;
  }
  //  COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocmap,"MAP");
  if(ouf)
    *ouf << "FE2FD: Waiting on all procs to finish reducing." << std::endl;
  BaseComm.Barrier();
  if(ouf)
    *ouf << "FE2FD: All procs done reducing attributes." << std::endl;
  //  int nprocs = BaseComm.Size();
  std::vector<int> pane_ids;
  COM_get_panes(fdwin.c_str(),pane_ids);
  if(ouf)
    *ouf << "FE2FD: Number of panes: " << pane_ids.size() << std::endl;
  // Form a list of unique block id's across all processors
  std::vector<int> local_block_ids;
  TRAIL_GetPanelAttribute(fdwin,"block_id",local_block_ids);
  if(ouf){
    *ouf << "FE2FD: Local Block_ids: ";
    IRAD::Util::DumpContents(*ouf,local_block_ids," ");
    *ouf << std::endl;
  }
  std::vector<int> global_blocks;
  TRAIL_UniqueAcrossProcs(local_block_ids,global_blocks,BaseComm.World());
  if(ouf){
    *ouf << "FE2FD: Global Block_ids: ";
    IRAD::Util::DumpContents(*ouf,local_block_ids," ");
    *ouf << std::endl;
  }
  // Now all_block_ids is an array of all the unique block ids across all procs
  // For each block
  std::vector<int>::iterator bii = global_blocks.begin();
  while(bii != global_blocks.end()){
    int block_id = *bii++;
    std::vector<int>::iterator fi = std::find(local_block_ids.begin(),local_block_ids.end(),block_id);
    int block_color = 0;
    if(fi != local_block_ids.end()) // Then this processor has data for the given block
      block_color = 1;
    // Split the communicator into haves and have nots for this block
    if(ouf)
      *ouf << "FE2FD: Splitting communicator for blocks." << std::endl;
    IRAD::Comm::CommunicatorObject BlockComm;
    BaseComm.Split(block_color,rank,BlockComm);
    if(block_color == 1){ // all the guys with data on this block
      int block_nproc = BlockComm.Size();
      int block_rank  = BlockComm.Rank();
      // Form a list of unique patch id's across all processors
      std::vector<int> local_patch_ids;
      std::vector<int> panes;
      TRAIL_ExtractPanes(fdwin,"block_id",block_id,panes); // get pane id's of panes with block id = block_id
      std::vector<int>::iterator pi = panes.begin();
      while(pi != panes.end()){
	int *patch_id;
	COM_get_array((fdwin+".patch_id").c_str(),*pi,&patch_id);
	if(patch_id)
	  local_patch_ids.push_back(*patch_id);
	pi++;
      }
      std::vector<int> all_patches(block_nproc);
      TRAIL_UniqueAcrossProcs(local_patch_ids,all_patches,BlockComm.World());
      // Now all_patches is an array of all the unique patch id across all procs having the current block
      std::vector<int>::iterator pii = all_patches.begin();
      while(pii != all_patches.end()){ // For each patch
	int patch_id = *pii++;
	// Determine if local processor owns part of the patch
	// Split the communicator into haves and have nots
	int patch_color = 0;
	IRAD::Comm::CommunicatorObject PatchComm;
	std::vector<int>::iterator fp = std::find(local_patch_ids.begin(),local_patch_ids.end(),patch_id);
	if(fp != local_patch_ids.end())
	  patch_color = 1;
	if(ouf)
	  *ouf << "FE2FD: Splitting communicator for patches." << std::endl;
	BlockComm.Split(patch_color,block_rank,PatchComm);
	if(patch_color == 1) { // all of us that have data on the given block/patch
	  //	  int patch_nproc = PatchComm.Size();
	  //	  int patch_rank  = PatchComm.Rank();
	  std::vector<int> patch_pane;
	  TRAIL_ExtractPanes(fdwin,"patch_id",patch_id,patch_pane); // get the pane for this patch (hopefully only 1)
	  //  std::vector<int>::iterator pii = pane_ids.begin();
	  //  while(pii != pane_ids.end()){
	  int pane_id = patch_pane[0];
	  int err = 0;
	  std::vector<int> flat_uns_extent;
	  err = TRAIL_Att2Vec(fewin+".local_extent",pane_id,flat_uns_extent);
	  std::vector<int> flat_fd_extent;
	  err = TRAIL_Att2Vec(fdwin+".local_extent",pane_id,flat_fd_extent);
	  assert(err == 0);
	  Mesh::BSExtent<int> UnsExtent(flat_uns_extent);
	  Mesh::BSExtent<int> FDExtent(flat_fd_extent);
	  std::vector<int> src_indices;
	  UnsExtent.GetFlatIndices(FDExtent,src_indices);
	  std::vector<int> asizes(atts.size(),0);
	  int attindex = 0;
	  std::vector<std::string>::iterator atti = atts.begin();
	  while(atti != atts.end()){
	    std::string att = *atti++;
	    int ncomp = 0;
	    char loc;
	    COM_Type type;
	    std::string unit;
	    COM_get_dataitem((fewin+"."+att).c_str(),&loc,&type,&ncomp,&unit);
	    asizes[attindex++] = COM_get_sizeof(type,ncomp);
	  }
	  if(ouf)
	    *ouf << "FE2FD: All attributes sized up." << std::endl;

	  atti = atts.begin();
	  attindex = 0;
	  if(ouf)
	    *ouf << "FE2FD: Copying buffers." << std::endl;
	  while(atti != atts.end()){
	    std::string attstr = fewin+"."+*atti;
	    std::string trgatt = fdwin+"."+*atti++;
	    char *srcptr = NULL;
	    char *trgptr = NULL;
	    COM_get_array(attstr.c_str(),pane_id,&srcptr);
	    COM_get_array(trgatt.c_str(),pane_id,&trgptr);
	    std::vector<int>::iterator sii = src_indices.begin();
	    int count = 0;
	    while(sii != src_indices.end()){
	      int src_index = *sii++ - 1;
	      memcpy(&trgptr[count++*asizes[attindex]],&srcptr[src_index*asizes[attindex]],asizes[attindex]);
	    }
	    attindex++;
	  }
	}
      }
    }
  }
  if(ouf)
    *ouf << "FD2FE: Waiting at final barrier" << std::endl;
  BaseComm.Barrier();
  if(ouf)
    *ouf << "FD2FE: Exit" << std::endl;
  COM_set_default_communicator(oldcomm);
  return(0);
}

int
TRAIL_Add_Attributes(const std::string &srcwin,
		    const std::string &destwin,
		    const std::vector<std::string> &atts)
{
  std::vector<std::string>::const_iterator ai = atts.begin();
  while(ai != atts.end()){
    std::string att = *ai++;
    std::string tatt = destwin+"."+att;
    std::string satt = srcwin+"."+att;
    // First, make sure the destination attribute exists, if not then create and size it
    {
      if(COM_get_dataitem_handle(tatt) <= 0){
	int ncomp = 0;
	char loc;
	COM_Type type;
	std::string unit;
	COM_get_dataitem(satt,&loc,&type,&ncomp,&unit);
	COM_new_dataitem(tatt,loc,type,ncomp,unit.c_str());
	//	COM_window_init_done(destwin,0);
      }
    }
  }
  COM_window_init_done(destwin,0);
  return(0);
}

int
TRAIL_FD2FE_Transfer(const std::string &srcwin,const std::string &destwin,
		    const std::string &attlist,std::ostream *ouf)
{
  std::vector<std::string> atts;
  std::istringstream Istr(attlist);
  std::string attstr;
  while(Istr >> attstr)
    atts.push_back(attstr);
  // Get the window communicator
  MPI_Comm communicator = MPI_COMM_NULL;
  COM_get_communicator(srcwin.c_str(),&communicator);
  IRAD::Comm::CommunicatorObject BaseComm(communicator);
  IRAD::Comm::CommunicatorObject BlockComm;
  IRAD::Comm::CommunicatorObject PatchComm;
  // Make sure the attribute exists
  TRAIL_Add_Attributes(srcwin,destwin,atts);
  //  int nprocs = BaseComm.Size();
  int rank   = BaseComm.Rank();
  std::vector<int> pane_ids;
  COM_get_panes(srcwin.c_str(),pane_ids);
  if(ouf && !rank)
    *ouf << "TRAIL: Number of panes: " << pane_ids.size() << std::endl;
  // Form a list of unique block id's across all processors
  std::vector<int> local_block_ids;
  TRAIL_GetPanelAttribute(srcwin,"block_id",local_block_ids);
  if(ouf && !rank){
    *ouf << "TRAIL: Local Block_ids: ";
    IRAD::Util::DumpContents(*ouf,local_block_ids," ");
    *ouf << std::endl;
  }
  std::vector<int> global_blocks;
  TRAIL_UniqueAcrossProcs(local_block_ids,global_blocks,BaseComm.World());
  if(ouf && !rank){
    *ouf << "TRAIL: Global Block_ids: ";
    IRAD::Util::DumpContents(*ouf,local_block_ids," ");
    *ouf << std::endl;
  }
  // Now all_block_ids is an array of all the unique block ids across all procs
  // For each block
  if(true){
    std::vector<int>::iterator bii = global_blocks.begin();
    while(bii != global_blocks.end()){
      int block_id = *bii++;
      std::vector<int>::iterator fi = std::find(local_block_ids.begin(),local_block_ids.end(),block_id);
      int block_color = 0;
      if(fi != local_block_ids.end()) // Then this processor has data for the given block
	block_color = 1;
      // Split the communicator into haves and have nots for this block
      if(ouf && !rank)
	*ouf << "TRAIL: Splitting communicator for blocks." << std::endl;
      //    IRAD::Comm::CommunicatorObject BlockComm;
      BaseComm.Split(block_color,rank,BlockComm);
      if(true){
	if(block_color == 1){ // all the guys with data on this block
	  int block_nproc = BlockComm.Size();
	  int block_rank  = BlockComm.Rank();
	  // Form a list of unique patch id's across all processors
	  std::vector<int> local_patch_ids;
	  std::vector<int> panes;
	  TRAIL_ExtractPanes(srcwin,"block_id",block_id,panes); // get pane id's of panes with block id = block_id
	  std::vector<int>::iterator pi = panes.begin();
	  while(pi != panes.end()){
	    int *patch_id;
	    COM_get_array((srcwin+".patch_id").c_str(),*pi,&patch_id);
	    if(patch_id)
	      local_patch_ids.push_back(*patch_id);
	    pi++;
	  }
	  std::vector<int> all_patches(block_nproc);
	  TRAIL_UniqueAcrossProcs(local_patch_ids,all_patches,BlockComm.World());
	  // Now all_patches is an array of all the unique patch id across all procs having the current block
	  std::vector<int>::iterator pii = all_patches.begin();
	  while(pii != all_patches.end()){ // For each patch
	    int patch_id = *pii++;
	    // Determine if local processor owns part of the patch
	    // Split the communicator into haves and have nots
	    int patch_color = 0;
	    //	IRAD::Comm::CommunicatorObject PatchComm;
	    std::vector<int>::iterator fp = std::find(local_patch_ids.begin(),local_patch_ids.end(),patch_id);
	    if(fp != local_patch_ids.end())
	      patch_color = 1;
	    if(ouf && !rank)
	      *ouf << "TRAIL: Splitting communicator for patches." << std::endl;
	    BlockComm.Split(patch_color,block_rank,PatchComm);
	    if(true){
	      if(patch_color == 1) { // all of us that have data on the given block/patch
		int patch_nproc = PatchComm.Size();
		//	  int patch_rank  = PatchComm.Rank();
		std::vector<int> patch_pane;
		TRAIL_ExtractPanes(srcwin,"patch_id",patch_id,patch_pane); // get the pane for this patch (hopefully only 1)
		//  std::vector<int>::iterator pii = pane_ids.begin();
		//  while(pii != pane_ids.end()){
		int pane_id = patch_pane[0];
		std::vector<int> all_pane_ids(patch_nproc);
		PatchComm.AllGather(pane_id,all_pane_ids);
		int err = 0;
		std::vector<int> flat_uns_extent;
		err = TRAIL_Att2Vec(destwin+".local_extent",pane_id,flat_uns_extent);
		std::vector<int> flat_global_extent;
		err = TRAIL_Att2Vec(destwin+".global_extent",pane_id,flat_global_extent);
		std::vector<int> shared_panes;
		err = TRAIL_Att2Vec(destwin+".shared_panes",pane_id,shared_panes);
		std::vector<int> flat_shared_extents;
		err = TRAIL_Att2Vec(destwin+".shared_extents",pane_id,flat_shared_extents);
		std::vector<int> flat_fd_extent;
		err = TRAIL_Att2Vec(srcwin+".local_extent",pane_id,flat_fd_extent);
		assert(err == 0);
		
		if(ouf){
		  *ouf << "fd_extent: (["
		       << flat_fd_extent[0] << "," << flat_fd_extent[1] << "],[" 
		       << flat_fd_extent[2] << "," << flat_fd_extent[3] << "],["
		       << flat_fd_extent[4] << "," << flat_fd_extent[5] << "])"
		       << std::endl;
		  *ouf << "unstructured_extent: (["
		       << flat_uns_extent[0] << "," << flat_uns_extent[1] << "],[" 
		       << flat_uns_extent[2] << "," << flat_uns_extent[3] << "],["
		       << flat_uns_extent[4] << "," << flat_uns_extent[5] << "])"
		       << std::endl;
		  *ouf << "global_extent: (["
		       << flat_global_extent[0] << "," << flat_global_extent[1] << "],[" 
		       << flat_global_extent[2] << "," << flat_global_extent[3] << "],["
		       << flat_global_extent[4] << "," << flat_global_extent[5] << "])"
		       << std::endl;
		  for(int iii = 0;iii < flat_shared_extents.size()/6;iii++){
		    unsigned int index = iii*6;
		    *ouf << "shared_extents (" << shared_panes[iii] << "): (["
			 << flat_shared_extents[index+0] << "," << flat_shared_extents[index+1] << "],[" 
			 << flat_shared_extents[index+2] << "," << flat_shared_extents[index+3] << "],["
			 << flat_shared_extents[index+4] << "," << flat_shared_extents[index+5] << "])"
			 << std::endl;
		  }
		}
		
		std::vector<int> all_fd_extents(6*patch_nproc);
		//	  PatchComm.AllGather<std::vector<int>,int>(flat_fd_extent,all_fd_extents,6,6);
		PatchComm.AllGather(flat_fd_extent,all_fd_extents,6,6);
		if(ouf && !rank)
		  *ouf << "TRAIL: All FD Extents communicated." << std::endl;
		// Now need to determine the interpane communication lists so that
		// we can gather the FD representation data needed to populate the
		// full extent of the FE representation.
		Mesh::BSExtent<int> UnsExtent(flat_uns_extent);
		Mesh::BSExtent<int> GlobalExtent(flat_global_extent);
		Mesh::BSExtent<int> FDExtent(flat_fd_extent);
		
		// Get the indices of the original fd gridpoints wrt the fe mesh.
		std::vector<int> theflatindices;
		UnsExtent.GetFlatIndices(FDExtent,theflatindices);
		
		std::vector<Mesh::BSExtent<int> > AllFDExtents;
		std::vector<Mesh::BSExtent<int> > SharedExtents;
		std::vector<Mesh::BSExtent<int> > RecvExtents;
		std::vector<int> remote_ranks;
		//		std::vector<int>::iterator spi = shared_panes.begin();
		std::vector<std::vector<char> > SendBuffers;
		std::vector<std::vector<char> > RecvBuffers;
		std::vector<std::vector<int>  > FlatIndices;
		std::vector<int> recvranks;
		std::vector<int> sendranks;
		for(int i = 0;i < patch_nproc;i++)
		  AllFDExtents.push_back(Mesh::BSExtent<int>(&all_fd_extents[i*6]));
		
		int count = 0;
		int nsend = 0;
		int nrecv = 0;
		std::vector<std::string>::iterator ai = atts.begin();
		std::vector<int> asizes(atts.size(),0);
		size_t bsize = 0; // block size = all attribute sizes added up
		int attindex = 0;
		while(ai != atts.end()){
		  std::string att = *ai++;
		  int ncomp = 0;
		  char loc;
		  COM_Type type;
		  std::string unit;
		  if(ouf && !rank)
		    *ouf << "TRAIL: Processing attribute " << att << std::endl;
		  COM_get_dataitem((srcwin+"."+att).c_str(),&loc,&type,&ncomp,&unit);
		  asizes[attindex] = COM_get_sizeof(type,ncomp);
		  bsize += asizes[attindex++];
		}
		if(ouf && !rank)
		  *ouf << "TRAIL: All attributes sized up." << std::endl;
		std::vector<int>::iterator spi = shared_panes.begin();
		while(spi != shared_panes.end()){
		  Mesh::BSExtent<int> SharedExtent(&flat_shared_extents[6*count++]);
		  SharedExtents.push_back(SharedExtent);
		  int rpane_id = *spi++;
		  // If the shared extent (FE rep) overlaps with
		  // the FDExtent, then those nodes are send nodes
		  int remote_rank = std::find(all_pane_ids.begin(),all_pane_ids.end(),rpane_id) - all_pane_ids.begin();
		  remote_ranks.push_back(remote_rank);
		  Mesh::BSExtent<int> CommExtent;
		  FDExtent.Overlap(SharedExtent,CommExtent);
		  if(!CommExtent.empty()){ // Then we need to send data to the remote processor
		    // For each attribute to send, pack and send the buffer (which needs to be persistent)
		    sendranks.push_back(remote_rank);
		    int nnodes = CommExtent.NNodes();
		    std::vector<std::string>::iterator ai = atts.begin();
		    std::vector<int> indices;
		    FDExtent.GetFlatIndices(CommExtent,indices);
		    //		    SendBuffers.push_back(std::vector<char>(nnodes*bsize,0));
		    std::vector<char> sendbufv(nnodes*bsize,0);
		    SendBuffers.push_back(sendbufv);
		    std::vector<int>::iterator ii = indices.begin();
		    size_t byte_offset = 0;
		    while(ii != indices.end()){
		      int attindex = 0;
		      int index = *ii++ - 1;
		      ai = atts.begin();
		      while(ai != atts.end()){
			std::string att = *ai++;
			char *data_ptr = NULL;
			COM_get_array((srcwin+"."+att).c_str(),pane_id,&data_ptr);
			//			memcpy(&SendBuffers[nsend][byte_offset],&data_ptr[(index-1)*asizes[attindex]],asizes[attindex]);
			memcpy(&SendBuffers[nsend][byte_offset],&data_ptr[index*asizes[attindex]],asizes[attindex]);
			
			byte_offset += asizes[attindex++];
		      }
		    }
		    //	      PatchComm.ASend(SendBuffers[nsend++],remote_rank);
		    nsend++;
		  }
		  else { // Then we need to receive data from the remote processor
		    // First, use the overlap function to ensure that we build a receive
		    // buffer for only the *real* nodes on the remote processor.  Each
		    // processor sends only real nodes as seen in the above send block
		    AllFDExtents[remote_rank].Overlap(SharedExtent,CommExtent);
		    if(!CommExtent.empty()){ // if it's empty, then there were only ghosts, no receive necessary
		      RecvExtents.push_back(CommExtent);
		      // For each attribute to receive build a buffer into which to receive
		      int nnodes = CommExtent.NNodes();
		      std::vector<char> recvbuf(nnodes*bsize,0);
		      std::vector<std::string>::iterator ai = atts.begin();
		      //		      RecvBuffers.push_back(std::vector<char>(nnodes*bsize,0));
		      RecvBuffers.push_back(recvbuf);
		      if(ouf && !rank)
			*ouf << "TRAIL: Receive buffer size for " << nnodes << " nodes is "
			     << nnodes*bsize << "." << std::endl;
		      recvranks.push_back(remote_rank);
		      //		      PatchComm.ARecv(RecvBuffers[nrecv++],remote_rank);
		    }
		  }
		}
		PatchComm.Barrier();
		std::vector<std::vector<char> >::iterator rbi = RecvBuffers.begin();
		std::vector<int>::iterator rri = recvranks.begin();
		//		 std::vector<std::vector<char> > testbufr(recvranks.size());
		while(rbi != RecvBuffers.end()){
		  //		   int iiindex = rbi - RecvBuffers.begin();
		  //		   testbufr[iiindex].resize(10000);
		  //		   PatchComm.ARecv(testbufr[iiindex],*rri);
		  if(ouf)
		    *ouf << "Receiving from " << *rri << std::endl;
		  PatchComm.ARecv(*rbi,*rri);
		  rbi++;
		  rri++;
		}
		if(ouf)
		  *ouf << "SendBuffers size = " << SendBuffers.size();
		std::vector<std::vector<char> >::iterator sbi = SendBuffers.begin();
		rri = sendranks.begin();
		//		 std::vector<std::vector<char> > testbufs(sendranks.size());
		while(sbi != SendBuffers.end()){
		  //		   int iiindex = sbi - SendBuffers.begin();
		  //		   testbufs[iiindex].resize(10000);
		  //		   PatchComm.ASend(testbufs[iiindex],*rri);
		  if(ouf)
		    *ouf << "Sending to " << *rri << std::endl;
		  PatchComm.ASend(*sbi,*rri);
		  sbi++;
		  rri++;
		}
		PatchComm.WaitAll();
		PatchComm.Barrier();
		//		 if(ouf){
		//		   *ouf << "TRAIL: All communication are completed." << std::endl << std::flush;
		//		 }	  
		
		if(ouf){
		  *ouf << "TRAIL: All communication are initiated." << std::endl;
		}
		// Now all the communcation is initiated. We'll do some
		// work while the communication completes.
		// Populate the FE representation with the local values
		//		{
		if(true){
		  // Determine the indices of the FD vertices in the FE
		  // representation
		  //		  std::vector<int> flatindices;
		  //		  UnsExtent.GetFlatIndices(FDExtent,flatindices);
		  PatchComm.Barrier();
		  if(ouf){
		    *ouf << "TRAIL: Getting indices of fd nodes in unstructured representation." << std::endl;
		    std::vector<int> uflat;
		    std::vector<int> dflat;
		    UnsExtent.Flatten(uflat);
		    FDExtent.Flatten(dflat);
		    *ouf << "TRAIL: finding fd_extent: (["
			 << dflat[0] << "," << dflat[1] << "],[" 
			 << dflat[2] << "," << dflat[3] << "],["
			 << dflat[4] << "," << dflat[5] << "])"
			 << std::endl;
		    *ouf << "TRAIL: in uns_extent: (["
			 << uflat[0] << "," << uflat[1] << "],[" 
			 << uflat[2] << "," << uflat[3] << "],["
			 << uflat[4] << "," << uflat[5] << "])"
			 << std::endl;
		  }
		  //	    std::vector<int> flatindices;
		  //	    UnsExtent.GetFlatIndices(FDExtent,flatindices);
		  //	    if(ouf){
		  //	      *ouf << "TRAIL: Made it past finding indices." << std::endl
		  //		   << "TRAIL: found overlap: ";
		  //	      IRAD::Util::DumpContents(*ouf,theflatindices," ");
		  //	      *ouf << std::endl << std::flush;
		  //	    }
		  
		  PatchComm.Barrier();
		  if(ouf){
		    *ouf << "TRAIL: Ready to copy local FD data into FE arrays." << std::endl << std::flush;
		  }
		  
		  // For each attribute, copy the values from the FD vertices
		  // to the FE vertices
		  int attindex = 0;
		  ai = atts.begin();
		  while(ai != atts.end()){
		    std::string att = *ai++;
		    std::string tatt = destwin+"."+att;
		    std::string satt = srcwin+"."+att;
		    // First, make sure the destination attribute exists, if not then create and size it
		    //		    {
		    if(COM_get_dataitem_handle(tatt) <= 0){
		      if(ouf && !rank)
			*ouf << "TRAIL: " << tatt << " did not exist." << std::endl;
		      int ncomp = 0;
		      char loc;
		      COM_Type type;
		      std::string unit;
		      COM_get_dataitem(satt,&loc,&type,&ncomp,&unit);
		      COM_new_dataitem(tatt,loc,type,ncomp,unit.c_str());
		      //		  COM_set_size(tatt,pane_id,UnsExtent.NNodes());
		      COM_resize_array(tatt,pane_id);
		      if(ouf && !rank)
			*ouf << "TRAIL: Initializing buffer to 0." << std::endl;
		      char *trg_ptr = NULL;
		      COM_get_array(tatt.c_str(),pane_id,&trg_ptr);
		      if(trg_ptr){
			int cnt = 0;
			for(int acount = 0;acount < UnsExtent.NNodes();acount++)
			  for(int ncomp = 0; ncomp < asizes[attindex];ncomp++)
			    trg_ptr[cnt++] = 0;
		      }
		      //		      COM_window_init_done(destwin,0);
		    }
		    else {
		      if(ouf && !rank)
			*ouf << "TRAIL: Attribute appeared to exist on target" << std::endl;
		      COM_resize_array(tatt,pane_id);
		      if(ouf && !rank)
			*ouf << "TRAIL: Initializing buffer to 0." << std::endl;
		      char *trg_ptr = NULL;
		      COM_get_array(tatt.c_str(),pane_id,&trg_ptr);
		      assert(trg_ptr != NULL);
		      if(trg_ptr){
			int cnt = 0;
			for(int acount = 0;acount < UnsExtent.NNodes();acount++)
			  for(int ncomp = 0; ncomp < asizes[attindex];ncomp++)
			    trg_ptr[cnt++] = 0;
		      }
		      //		      COM_window_init_done(destwin,0);
		    }
		    char *src_ptr = NULL;
		    char *trg_ptr = NULL;
		    COM_get_array(tatt.c_str(),pane_id,&trg_ptr);
		    assert(trg_ptr != NULL);
		    COM_get_array(satt.c_str(),pane_id,&src_ptr);
		    assert(src_ptr != NULL);
		    std::vector<int>::iterator fii = theflatindices.begin();
		    int count = 0;
		    while(fii != theflatindices.end()){
		      int index = *fii++;
		      memcpy(&trg_ptr[(index-1)*asizes[attindex]],&src_ptr[count++*asizes[attindex]],asizes[attindex]);
		    }
		    attindex++;
		  }
		  //		}
		  PatchComm.Barrier();
		}
		// Wait for pending communication
		//		PatchComm.WaitAll();
		//	  if(ouf && !rank)
		//	    *ouf << "TRAIL: All communication completed." << std::endl;
		PatchComm.Barrier();
		if(ouf)
		  *ouf << "TRAIL: Base copy completed." << std::endl << std::flush;
		// Populate the FE representation from the received info
		if(true){
		  // Copying out of the receive buffers is dead easy
		  std::vector<Mesh::BSExtent<int> >::iterator rei = RecvExtents.begin();
		  int recvindex = 0;
		  while(rei != RecvExtents.end()){
		    std::vector<int> flatindices;
		    unsigned int nnodes = rei->NNodes();
		    //		    UnsExtent.GetFlatIndices(*rei++,flatindices);
		    UnsExtent.GetFlatIndices(*rei,flatindices);
		    std::vector<int>::iterator fii = flatindices.begin();
		    assert(nnodes == flatindices.size());
		    for(unsigned int i = 0;i < nnodes;i++){
		      int index = *fii++;
		      int attindex = 0;
		      ai = atts.begin();
		      int byte_offset = 0;
		      while(ai != atts.end()){
			std::string tatt = destwin+"."+*ai;
			char *trg_ptr = NULL;
			COM_get_array(tatt.c_str(),pane_id,&trg_ptr);
			assert(trg_ptr != NULL);
			//			memcpy(&trg_ptr[(index-1)*asizes[attindex]],&RecvBuffers[recvindex][i*bsize+byte_offset],asizes[attindex]);
			memcpy(&trg_ptr[(index-1)*asizes[attindex]],
			       &RecvBuffers[recvindex][i*bsize+byte_offset],asizes[attindex]);
			byte_offset += asizes[attindex++];
			ai++;
		      }
		    }
		    recvindex++;
		    rei++;
		  }
		}
	      } // if(patch_color == 1)
	    } // if(true) [debugging]
	    PatchComm.Barrier();
	  } // loop through patches on this block
	} // if(block_color == 1)
      } // if(true) [debugging]
      BlockComm.Barrier();
    } // loop through global blocks
  } // if(true) [debugging]
  BaseComm.Barrier();
  if(ouf && !rank)
    *ouf << "TRAIL: All processors done with transfer." << std::endl;
  COM_window_init_done(destwin.c_str());
  BaseComm.Barrier();
  return(0);
}

/// Takes as input a block structured FD grid.  An FD grid is one in which
/// the partitioning is *node* based instead of element based as in an FE
/// mesh - and produces a FE representation of the mesh, including a
/// description of the shared vertices.   Attributes can then be moved
/// back and forth between the two representations easily.
///
/// Needs the following data for each patch/pane:
/// Block ID: block_id
/// Patch ID: patch_id
/// Local Patch Extent: local_patch_extent [istart,iend,jstart,jend,kstart,kend] (global indices for local patch extent)
/// Global Patch Extent: global_patch_extent [istart,iend,jstart,jend,kstart,kend]
///
/// The output window has *not* transferred any of the non-mesh attributes.  (Currently needs NC transferred)
/// Produces the following data for each patch/pane:
/// Block ID: block_id
/// Patch ID: patch_id
/// Local Extent: local_extent
/// Global Extent: global_extent
/// Shared Extent: shared_extent[local_extent1 local_extent2 .... local_extent_n]
/// Shared Panes: shared_panes[pane1 pane2 .... pane_n]
int
TRAIL_FD2FE_WinCreate2(const std::string &wname,const std::string &outwname,std::ostream *ouf)
{
  // Get the window communicator
  MPI_Comm communicator = MPI_COMM_NULL;
  COM_get_communicator(wname.c_str(),&communicator);
  IRAD::Comm::CommunicatorObject BaseComm(communicator);
  int nprocs = BaseComm.Size();
  int rank   = BaseComm.Rank();
  if(ouf)
    *ouf << "TRAIL_AddBlockStructuredGhostZone::Nprocs = " << nprocs << std::endl;
  std::vector<int> pane_ids;
  COM_get_panes(wname.c_str(),pane_ids);
  if(ouf)
    *ouf << "Number of panes: " << pane_ids.size() << std::endl;
  // Form a list of unique block id's across all processors
  std::vector<int> local_block_ids;
  TRAIL_GetPanelAttribute(wname,"block_id",local_block_ids);
  if(ouf){
    *ouf << "Local Block_ids: ";
    IRAD::Util::DumpContents(*ouf,local_block_ids," ");
    *ouf << std::endl;
  }
  std::vector<int> global_blocks;
  TRAIL_UniqueAcrossProcs(local_block_ids,global_blocks,BaseComm.World());
  if(ouf){
    *ouf << "Global Block_ids: ";
    IRAD::Util::DumpContents(*ouf,global_blocks," ");
    *ouf << std::endl;
  }
  // Now all_block_ids is an array of all the unique block ids across all procs
  // For each block
  //  std::vector<int>::iterator bii = global_blocks.begin();
  std::string wname2 = outwname; // wname+"_uns";
  COM_new_window(wname2.c_str());
  COM_new_dataitem((wname2+".local_extent").c_str(),  'p',COM_INTEGER,6,"");
  COM_new_dataitem((wname2+".global_extent").c_str(), 'p',COM_INTEGER,6,"");
  COM_new_dataitem((wname2+".shared_extents").c_str(),'p',COM_INTEGER,6,"");
  COM_new_dataitem((wname2+".shared_panes").c_str(),  'p',COM_INTEGER,1,"");
  COM_new_dataitem((wname2+".bcflag").c_str(),        'p',COM_INTEGER,1,"");
  //  COM_new_dataitem((wname2+".send_extent").c_str(),  'p',COM_INTEGER,6,"");
  //  COM_new_dataitem((wname2+".send_panes").c_str(),   'p',COM_INTEGER,1,"");
  //  COM_new_dataitem((wname2+".recv_panes").c_str(),   'p',COM_INTEGER,1,"");
  //  COM_new_dataitem((wname2+".recv_extent").c_str(),  'p',COM_INTEGER,6,"");
  COM_new_dataitem((wname2+".block_id").c_str(),     'p',COM_INTEGER,1,"");
  COM_new_dataitem((wname2+".patch_id").c_str(),     'p',COM_INTEGER,1,"");
  //  int bcflag = 1;
  //  COM_set_array((wname2+".bcflag").c_str()),
  //  COM_new_dataitem((wname+".send_extent").c_str(),   'p',COM_INTEGER,6,"");
  //  COM_new_dataitem((wname+".send_panes").c_str(),    'p',COM_INTEGER,1,"");
  //  COM_new_dataitem((wname+".recv_extent").c_str(),   'p',COM_INTEGER,6,"");
  //  COM_new_dataitem((wname+".recv_panes").c_str(),    'p',COM_INTEGER,1,"");

  // Loop through the global block ids on all processors
  std::vector<int>::iterator bii = global_blocks.begin();
  while(bii != global_blocks.end()){
    int block_id = *bii++;
    std::vector<int>::iterator fi = std::find(local_block_ids.begin(),local_block_ids.end(),block_id);
    int block_color = 0;
    if(fi != local_block_ids.end()) // Then this processor has data for the given block
      block_color = 1;
    // Split the communicator into haves and have nots for this block
    if(ouf)
      *ouf << "Splitting communicator for blocks." << std::endl;
    IRAD::Comm::CommunicatorObject BlockComm;
    BaseComm.Split(block_color,rank,BlockComm);
    if(block_color == 1){ // all the guys with data on this block
      int block_nproc = BlockComm.Size();
      int block_rank  = BlockComm.Rank();
      if(ouf){
	*ouf << "Processor " << rank << " has block rank "
	     << block_rank << "/" << block_nproc << std::endl;
      }
      // Each block can have multiple "patches". Each patch should be
      // contained in its own pane in the source window. confusing naming:
      // block = grid, patch = block, new word patch
      // Form a list of unique patch id's across all processors
      std::vector<int> local_patch_ids;
      std::vector<int> panes;
      TRAIL_ExtractPanes(wname,"block_id",block_id,panes); // get pane id's of panes with block id = block_id
      if(ouf){
	*ouf << "Found " << panes.size() << " local panes: ";
	IRAD::Util::DumpContents(*ouf,panes," ");
	*ouf << std::endl;
      }
      std::vector<int>::iterator pi = panes.begin();
      while(pi != panes.end()){
	int *patch_id;
	COM_get_array((wname+".patch_id").c_str(),*pi,&patch_id);
	if(patch_id)
	  local_patch_ids.push_back(*patch_id);
	pi++;
      }
      if(ouf){
	*ouf << "Found " << local_patch_ids.size() << " local patch_ids: ";
	IRAD::Util::DumpContents(*ouf,local_patch_ids," ");
	*ouf << std::endl;
      }
      std::vector<int> all_patches(block_nproc);
      TRAIL_UniqueAcrossProcs(local_patch_ids,all_patches,BlockComm.World());
      // Now all_patches is an array of all the unique patch id across all procs having the current block
      if(ouf){
	*ouf << "Found " << all_patches.size() << " global patch_ids: ";
	IRAD::Util::DumpContents(*ouf,all_patches," ");
	*ouf << std::endl;
      }
      std::vector<int>::iterator pii = all_patches.begin();
      while(pii != all_patches.end()){ // For each patch
	int patch_id = *pii++;
	// Determine if local processor owns part of the patch
	// Split the communicator into haves and have nots
	int patch_color = 0;
	IRAD::Comm::CommunicatorObject PatchComm;
	std::vector<int>::iterator fp = std::find(local_patch_ids.begin(),local_patch_ids.end(),patch_id);
	if(fp != local_patch_ids.end())
	  patch_color = 1;
	if(ouf)
	  *ouf << "Splitting communicator for patches." << std::endl;
	BlockComm.Split(patch_color,block_rank,PatchComm);
	if(patch_color == 1) { // all of us that have data on the given block/patch
	  int patch_nproc = PatchComm.Size();
	  int patch_rank  = PatchComm.Rank();
	  std::vector<int> patch_pane;
	  TRAIL_ExtractPanes(wname,"patch_id",patch_id,patch_pane); // get the pane for this patch (hopefully only 1)
	  int *global_patch_extent_ptr = NULL;
	  COM_get_array((wname+".global_extent").c_str(),patch_pane[0],&global_patch_extent_ptr);
	  if(!global_patch_extent_ptr){
	    if(ouf)
	      *ouf << "ERROR: Window " << wname << " has no global_patch_extent attribute." << std::endl;
	    std::cerr << "ERROR: Window " << wname << " has no global_patch_extent attribute." << std::endl;
	    assert(0);
	  }
	  if(ouf){
	    *ouf << "global_patch_extent: (["
		 << global_patch_extent_ptr[0] << "," << global_patch_extent_ptr[1] << "],[" 
		 << global_patch_extent_ptr[2] << "," << global_patch_extent_ptr[3] << "],["
		 << global_patch_extent_ptr[4] << "," << global_patch_extent_ptr[5] << "])"
		 << std::endl;
	  }
	  int *local_patch_extent_ptr = NULL;
	  COM_get_array((wname+".local_extent").c_str(),patch_pane[0],&local_patch_extent_ptr);
	  if(!local_patch_extent_ptr){
	    if(ouf)
	      *ouf << "ERROR: Window " << wname << " has no local_patch_extent attribute." << std::endl;
	    std::cerr << "ERROR: Window " << wname << " has no local_patch_extent attribute." << std::endl;
	    assert(0);
	  }
	  if(ouf){
	    *ouf << "local_patch_extent: (["
		 << local_patch_extent_ptr[0] << "," << local_patch_extent_ptr[1] << "],[" 
		 << local_patch_extent_ptr[2] << "," << local_patch_extent_ptr[3] << "],["
		 << local_patch_extent_ptr[4] << "," << local_patch_extent_ptr[5] << "])"
		 << std::endl;
	  }
	  std::vector<int> local_patch_extent;
	  std::vector<int> flat_local_patch_extents(6*patch_nproc);
	  IRAD::Util::CopyIntoContainer(local_patch_extent,local_patch_extent_ptr,6);
	  Mesh::BSExtent<int> LocalPatchExtent(local_patch_extent);
	  std::vector<int> global_patch_extent;
	  IRAD::Util::CopyIntoContainer(global_patch_extent,global_patch_extent_ptr,6);
	  Mesh::BSExtent<int> GlobalPatchExtent(global_patch_extent);

	  // Here's a rather simple bit of code that extends the patch to
	  // to the right if needed.
	  Mesh::BSExtent<int> ExtendedPatchExtent(LocalPatchExtent);
	  Mesh::BSExtent<int>::iterator lpei = LocalPatchExtent.begin();
	  Mesh::BSExtent<int>::iterator gpei = GlobalPatchExtent.begin();
	  Mesh::BSExtent<int>::iterator epei = ExtendedPatchExtent.begin();
	  while(lpei != LocalPatchExtent.end()){
	    if((*lpei)[1] < (*gpei)[1]){
	      (*epei)[1]++;
	      
	    }
	    lpei++;
	    gpei++;
	    epei++;
	  }
	  

	  // Communicate the extended patch extents to all procs having
	  // nodes in this patch
	  std::vector<int> extended_patch_extent;
	  ExtendedPatchExtent.Flatten(extended_patch_extent);
	  std::vector<int> all_extended_patch_extents(6*patch_nproc);
	  //	  PatchComm.AllGather<std::vector<int>,int>(extended_patch_extent,all_extended_patch_extents,6,6);
	  PatchComm.AllGather(extended_patch_extent,all_extended_patch_extents,6,6);
	  //	  PatchComm.AllGather<std::vector<int>,int>(local_patch_extent,all_local_patch_extents,6,6);
	  std::vector<int> all_pane_ids(patch_nproc);
	  //	  PatchComm.AllGather<std::vector<int>,int>(patch_pane[0],all_pane_ids);
	  PatchComm.AllGather(patch_pane[0],all_pane_ids);
	  if(ouf)
	    *ouf << "All extended patch extents communicated." << std::endl;

	  // Fill the search pool, these are the grid extents of all
	  // remote panes.
	  std::vector<Mesh::BSExtent<int> > ExtentPool;
	  for(int pindex = 0;pindex < patch_nproc;pindex++){
	    if(pindex != patch_rank){
	      Mesh::BSExtent<int> swimmer(&all_extended_patch_extents[6*pindex]);
	      ExtentPool.push_back(swimmer);
	    }
	  }

	  if(ouf)
	    *ouf << "Search pool filled." << std::endl;
	  // Find the shared nodes
	  std::vector<Mesh::BSExtent<int> > shared_extents;
	  std::vector<int> neighborranks;
	  std::vector<int> neighborpanes;
	  if(true){  // Scope temporary vars
	    std::vector<int> neighborindices;
	    if(ouf)
	      *ouf << "Finding shared nodes." << std::endl;
	    ExtendedPatchExtent.FindSharedNodes(ExtentPool,shared_extents,neighborindices);
	    if(ouf){
	      *ouf << "Shared nodes found." << std::endl
		   << "Neighbor indices: ";
	      IRAD::Util::DumpContents(*ouf,neighborindices," ");
	      *ouf << std::endl << "Shared extents: " << std::endl;
	      std::vector<Mesh::BSExtent<int> >::iterator bsmi = shared_extents.begin();
	      while(bsmi != shared_extents.end()){
		unsigned int nid = bsmi - shared_extents.begin();
		*ouf << "Shared extent " << nid << ":" << std::endl;
		std::vector<int> tempflat_extent;
		bsmi++->Flatten(tempflat_extent);
		*ouf << "shared_extent: (["
		     << tempflat_extent[0] << "," << tempflat_extent[1] << "],[" 
		     << tempflat_extent[2] << "," << tempflat_extent[3] << "],["
		     << tempflat_extent[4] << "," << tempflat_extent[5] << "])"
		     << std::endl;
	      }
	    }
	    std::vector<int>::iterator nii = neighborindices.begin();
	    while(nii != neighborindices.end()){
	      int nrank = *nii++;
	      if(nrank >= patch_rank)
		nrank++;
	      neighborranks.push_back(nrank);
	      neighborpanes.push_back(all_pane_ids[nrank]);
	    }
	  }
	  if(ouf){
	    *ouf << "Neighbor ranks: ";
	    IRAD::Util::DumpContents(*ouf,neighborranks," ");
	    *ouf << std::endl << "Neighborpanes: ";
	    IRAD::Util::DumpContents(*ouf,neighborpanes," ");
	    *ouf << std::endl;
	  }

	  if(ouf)
	    *ouf << "Copying attributes." << std::endl;
	  TRAIL_Copy2Attribute((wname2+".local_extent"),extended_patch_extent,patch_pane[0]);
	  TRAIL_Copy2Attribute((wname2+".global_extent"),global_patch_extent,patch_pane[0]);
	  TRAIL_SetAttribute((wname2+".block_id"),patch_pane[0],block_id);
	  TRAIL_SetAttribute((wname2+".patch_id"),patch_pane[0],patch_id);
	  Mesh::Connectivity conn;
	  ExtendedPatchExtent.Sync();
	  if(ouf)
	    *ouf << "Creating unstructured mesh." << std::endl;
	  ExtendedPatchExtent.CreateUnstructuredMesh(conn);
	  conn.Sync();
	  conn.SyncSizes();
	  bool flip = false;
	  int flip_hndl = COM_get_dataitem_handle((wname+".flip").c_str());
	  if(flip_hndl > 0){
	    int *flip_ptr = NULL;
	    COM_get_array((wname+".flip").c_str(),patch_pane[0],&flip_ptr);
	    if(flip_ptr != NULL)
	      flip = (*flip_ptr > 0);
	  }
	  if(flip){
	    if(ouf)
	      *ouf << "Flipping unstructured surface." << std::endl;
	    Mesh::Connectivity::iterator connit = conn.begin();
	    while(connit != conn.end())
	      {
		Mesh::GenericCell_2 ge(connit->size());
		ge.ReOrient(*connit++);
	      }
	  }
// 	  if(ouf){
// 	    std::vector<int> tempout;
// 	    ExtendedPatchExtent.Flatten(tempout);
// 	    *ouf << "Unstructured extent: ";
// 	    IRAD::Util::DumpContents(*ouf,tempout," ");
// 	    *ouf << std::endl;
// 	    *ouf << "Connectivity: " << std::endl
// 		 << conn << std::endl;
// 	  }
	  if(ouf)
	    *ouf << "Setting window data." << std::endl;
	  std::vector<int> flatcon;
	  unsigned int nnodes_new = ExtendedPatchExtent.NNodes();
	  conn.Flatten(flatcon);
	  COM_set_size((wname2+".nc").c_str(),patch_pane[0],nnodes_new);
	  COM_resize_array((wname2+".nc").c_str(),patch_pane[0]);
	  COM_set_size((wname2+".:q4:").c_str(),patch_pane[0],conn.Nelem());
	  COM_resize_array((wname2+".:q4:").c_str(),patch_pane[0]);
	  COM_set_size((wname2+".bcflag").c_str(),patch_pane[0],1);
	  COM_resize_array((wname2+".bcflag").c_str(),patch_pane[0]);
	  int *bcflag;
	  COM_get_array((wname2+".bcflag").c_str(),patch_pane[0],&bcflag);
	  *bcflag = 1;

	  int *conndat;
	  COM_get_array((wname2+".:q4:").c_str(),patch_pane[0],&conndat);
	  memcpy(conndat,&flatcon[0],sizeof(int)*4*conn.Nelem());
	  std::vector<std::vector<int> > flattened_shared_extents;
	  std::vector<Mesh::BSExtent<int> >::iterator sei = shared_extents.begin();
	  while(sei != shared_extents.end()){
	    std::vector<int> flatextent;
	    sei++->Flatten(flatextent);
	    flattened_shared_extents.push_back(flatextent);
	  }
	  TRAIL_Copy2Attribute((wname2+".shared_panes"),neighborpanes,patch_pane[0],neighborpanes.size());
	  TRAIL_Copy2Attribute((wname2+".shared_extents"),flattened_shared_extents,patch_pane[0]);
	  if(ouf)
	    *ouf << "Building pconn." << std::endl;
	  std::vector<int> pconn;
	  pconn.push_back(neighborpanes.size());
	  std::vector<int> sortedpanes(neighborpanes);
	  std::sort(sortedpanes.begin(),sortedpanes.end());
	  std::vector<int>::iterator spit = sortedpanes.begin();
	  while(spit != sortedpanes.end()){
	    int remote_pane_id = *spit++;
	    pconn.push_back(remote_pane_id);
	    int remote_pane_index = std::find(neighborpanes.begin(),neighborpanes.end(),remote_pane_id) - neighborpanes.begin();
	    std::vector<int> sharednodes;
	    ExtendedPatchExtent.GetFlatIndices(shared_extents[remote_pane_index],sharednodes);
	    pconn.push_back(sharednodes.size());
	    std::vector<int>::iterator sni = sharednodes.begin();
	    while(sni != sharednodes.end())
	      pconn.push_back(*sni++);
	  }
	  TRAIL_Copy2Attribute((wname2+".pconn"),pconn,patch_pane[0],pconn.size());
	  if(ouf){
	    *ouf << "PConn: ";
	    IRAD::Util::DumpContents(*ouf,pconn," ");
	    *ouf << std::endl;
	  }
	}
	PatchComm.Barrier();
      }
    }
    BlockComm.Barrier();
  }
  if(ouf)
    *ouf << "window done, roccom finalizing it." << std::endl;
  COM_window_init_done(wname2.c_str());
  if(ouf)
    *ouf << "window finalized." << std::endl;
  if(ouf)
    *ouf << "end barrier" << std::endl;
  BaseComm.Barrier();
  if(ouf)
    *ouf << "all processors exiting." << std::endl;
  return(0);
}

void TRAIL_Window2UnstructuredMesh(const std::string &wname,std::vector<Mesh::UnstructuredMesh> &meshes,
				  std::vector<SolnMetaData> &smdv,
				  std::vector<std::vector<std::vector<double> > > &soln_data,
				  int verblevel, bool no_ghost)
{
  std::vector<int> pane_ids;
  std::vector<std::vector<std::vector<double> > > temp_soln;
  TRAIL_GetWindowSolnMetaData(wname,smdv,verblevel>1);
  TRAIL_GetWindowSolnData(wname,temp_soln,smdv,verblevel>1);
  soln_data.resize(temp_soln.size());
  meshes.resize(temp_soln.size());
  COM_get_panes(wname.c_str(),pane_ids);
  if(verblevel > 1)
    std::cout << "TRAIL_Window2UnstructuredMesh::Processing " << pane_ids.size()
	      << " panes." << std::endl;
  std::vector<int>::iterator pii = pane_ids.begin();
  while(pii != pane_ids.end()){
    int pane_index = pii - pane_ids.begin();
    int pane_id = *pii++;
    const int *conn_ptr = NULL;
    double *coordinates_ptr = NULL;
    COM_get_array_const((wname+".nc").c_str(),pane_id,(const double **)&coordinates_ptr);
    int nreal = 0;
    int nghost = 0;
    COM_get_array_const((wname+".:st3:").c_str(),pane_id,&conn_ptr);
    int nnodes = 0;
    if(conn_ptr){
      COM_get_size((wname+".:st3:").c_str(),pane_id,&nreal,&nghost);
      if(verblevel)
	std::cout << "Block Structured pane " << pane_id << ": " << std::endl
		  << "Coordinate array size (i,j,k): (" << conn_ptr[0]
		  << "," << conn_ptr[1] << "," << conn_ptr[2] << ")" << std::endl
		  << "Ghost zone is " << nghost << " cells wide."
		  << std::endl;


      nnodes = conn_ptr[0]*conn_ptr[1]*conn_ptr[2];
      meshes[pane_index].nc.init_copy(nnodes,coordinates_ptr);
      std::vector<Mesh::IndexType> extent;
      for(int i = 0;i < 3;i++){
	extent.push_back(1);
	extent.push_back(conn_ptr[i]);
      }
      // resize soln_data to the number of attributes
      soln_data[pane_index].resize(temp_soln[pane_index].size());

      Mesh::BSExtent<Mesh::IndexType> bsextent(extent);
      if(no_ghost){
	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Stripping ghosts from block structured Window"
		    << std::endl;
	std::vector<Mesh::IndexType> real_extent;
	for(int i = 0;i < 3;i++){
	  real_extent.push_back(1+nghost);
	  real_extent.push_back(conn_ptr[i]-nghost);
	}
	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Getting real node indices"
		    << std::endl;
	Mesh::BSExtent<Mesh::IndexType> realextent(real_extent);
	std::vector<Mesh::IndexType> real_indices;
	bsextent.GetFlatIndices(realextent,real_indices);
	unsigned int nreal_nodes = real_indices.size();
	std::vector<double> real_coordinates(nreal_nodes*3);
	std::vector<Mesh::IndexType>::iterator rii = real_indices.begin();
	Mesh::IndexType rni = 0;
	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Extracting real coordinates"
		    << std::endl;
	while(rii != real_indices.end()){
	  Mesh::IndexType cindex = *rii++;
	  double *crdptr = meshes[pane_index].nc[cindex];
	  real_coordinates[rni++] = *crdptr++;
	  real_coordinates[rni++] = *crdptr++;
	  real_coordinates[rni++] = *crdptr;
	}
	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Replacing coordinates with real"
		    << std::endl;
	meshes[pane_index].nc.init_copy(nreal_nodes,&real_coordinates[0]);
	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Real nodal coordinates extracted."
		    << std::endl;
	for(int i = 0;i < 3;i++){
	  real_extent[i*2]   = 1;
	  real_extent[i*2+1] -= nghost;
	}
	Mesh::BSExtent<Mesh::IndexType> rext(real_extent);
	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Creating mesh with real vertices/cells"
		    << std::endl;
	rext.CreateUnstructuredMesh(meshes[pane_index].con);
	meshes[pane_index].con.Sync();
	meshes[pane_index].con.SyncSizes();
	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Nodes: " << meshes[pane_index].nc.size()
		    << ", Cells: " << meshes[pane_index].con.Nelem() << std::endl
		    << "TRAIL_Window2UnstructuredMesh::Extracting real cell indices."
		    << std::endl;

	std::vector<Mesh::IndexType> cell_extent;
	for(int i = 0;i < 3;i++){
	  cell_extent.push_back(1);
	  cell_extent.push_back(conn_ptr[i]-1);
	}
	Mesh::BSExtent<Mesh::IndexType> cellextent(cell_extent);
	std::vector<Mesh::IndexType> real_cell_extent;
	for(int i = 0;i < 3;i++){
	  real_cell_extent.push_back(1+nghost);
	  real_cell_extent.push_back(conn_ptr[i]-1-nghost);
	}
	Mesh::BSExtent<Mesh::IndexType> realcellextent(real_cell_extent);
	std::vector<Mesh::IndexType> real_cell_indices;
	cellextent.GetFlatIndices(realcellextent,real_cell_indices);

	if(verblevel > 1)
	  std::cout << "TRAIL_Window2UnstructuredMesh::Extracting real solution data"
		    << std::endl;
	std::vector<std::vector<double> >::iterator tsai = temp_soln[pane_index].begin();
	std::vector<std::vector<double> >::iterator sdai = soln_data[pane_index].begin();
	while(tsai != temp_soln[pane_index].end()){
	  unsigned int att_index = tsai - temp_soln[pane_index].begin();
	  unsigned int ncomp = smdv[att_index].ncomp;
	  if(!tsai->empty()){ // make sure there's data for this attribute
	    if(smdv[att_index].loc == 'N' || smdv[att_index].loc == 'n'){
	      // It's a nodal attribute
	      sdai->resize(real_indices.size()*ncomp);
	      if(verblevel > 1)
		std::cout << "TRAIL_Window2UnstructuredMesh::Processing "
			  << ncomp << " nodal components of "
			  << smdv[att_index].name << std::endl;
	      std::vector<double>::iterator sddi = sdai->begin();
	      rii = real_indices.begin();
	      while(rii != real_indices.end()){
		Mesh::IndexType real_index = *rii++ - 1;
		real_index *= ncomp;
		for(unsigned int ncc = 0;ncc < ncomp;ncc++)
		  *sddi++ = temp_soln[pane_index][att_index][real_index++];

	      }
	    }
	    else if(smdv[att_index].loc == 'E' || smdv[att_index].loc == 'e'){
	      // It's a cellular attribute
	      sdai->resize(real_cell_indices.size()*ncomp);
	      std::vector<double>::iterator sddi = sdai->begin();
	      if(verblevel > 1)
		std::cout << "TRAIL_Window2UnstructuredMesh::Processing "
			  << ncomp << " cellular components of "
			  << smdv[att_index].name << std::endl;
	      rii = real_cell_indices.begin();
	      while(rii != real_cell_indices.end()){
		Mesh::IndexType real_index = *rii++ - 1;
		real_index *= ncomp;
		for(unsigned int ncc = 0;ncc < ncomp;ncc++)
		  *sddi++ = temp_soln[pane_index][att_index][real_index++];

	      }
	    }
	  }
	  else
	    sdai->resize(0);
	  tsai->resize(0);
	  sdai++;
	  tsai++;
	}
	temp_soln[pane_index].resize(0);
      }
      else{
	bsextent.CreateUnstructuredMesh(meshes[pane_index].con);
 	std::vector<std::vector<double> >::iterator tsai = temp_soln[pane_index].begin();
 	std::vector<std::vector<double> >::iterator sdai = soln_data[pane_index].begin();
 	while(tsai != temp_soln[pane_index].end()){
	  // 	  unsigned int att_index = tsai - temp_soln[pane_index].begin();
	  // 	  unsigned int ncomp = smdv[att_index].ncomp;
 	  sdai->resize(tsai->size());
 	  std::vector<double>::iterator sddi = sdai->begin();
 	  std::vector<double>::iterator tddi = tsai->begin();
 	  while(tddi != tsai->end())
 	    *sddi++ = *tddi++;
	  tsai->resize(0);
 	  tsai++;
	  sdai++;
 	}
	temp_soln[pane_index].resize(0);
      }
      //      meshes.push_back(mesh);
    }
    else{
      if(verblevel > 1)
	std::cout << "TRAIL_Window2UnstructuredMesh::Processing unstructured window" << std::endl;
      // its unstructured - easier to deal with (i.e. it's already
      // unstructured!
    }
    //    COM_get_array_cont((
	//		<< "ST3 for geopane " << pane_id << std::endl
	//		<< geo_array_ptr[0] << "," << geo_array_ptr[1]
	//		<< "," << geo_array_ptr[2] << std::endl;
  }
//   std::vector<Mesh::UnstructuredMesh>::iterator mi = meshes.begin();
//   std::ofstream Ouf;
//   Ouf.open("test_coords");
//   while(mi != meshes.end()){
//     GeoPrim::CBox mesh_box;
//     GeoPrim::CBox small_box;
//     GeoPrim::CBox large_box;
//     Ouf << mi->nc << std::endl;
//     Mesh::GetMeshBoxes(mi->nc,mi->con,mesh_box,small_box,large_box);
//     std::cout << "Mesh BOXES:" << std::endl
// 	      << "Mesh:  " << mesh_box << std::endl
// 	      << "Small: " << small_box << std::endl
// 	      << "Large: " << large_box << std::endl;
//     mi++;
//   }
//   Ouf.close();
}
//#endif

// Serial Function:
// Read the HDF file fname into the window wname
void TRAIL_HDF2Window( const std::string &fname, const std::string &wname,int verb) {
  COM_new_window( wname.c_str());
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");
  int IN_read;
  IN_read = COM_get_function_handle( "IN.read_window");
  MPI_Comm comm_null = MPI_COMM_NULL;
  std::string bufwin("bufwin");
  if(verb)
    std::cout << "Reading file " << fname << "..." << std::endl;
  COM_call_function( IN_read, fname.c_str(), bufwin.c_str(), &comm_null);
  if(verb)
    std::cout << "Done reading file " << fname << "." << std::endl;
  int IN_obtain = COM_get_function_handle( "IN.obtain_attribute");
  int buf_all = COM_get_dataitem_handle((bufwin+".all").c_str());
  COM_call_function( IN_obtain, &buf_all, &buf_all);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC( SimIN, "IN");
  if(verb)
    std::cout << "Obtained temp window from file " << fname
	      << ", cloning.." << std::endl;
  COM_window_init_done(bufwin.c_str());
  // Change the memory layout to contiguous
  COM_clone_dataitem( (wname+".all").c_str(), (bufwin+".all").c_str(), 1);
  COM_window_init_done(wname.c_str());
  COM_delete_dataitem((bufwin+".data").c_str());
  COM_delete_window( bufwin.c_str());
  if(verb)
    std::cout << "Window " << wname << " created." << std::endl;

}

void
TRAIL_GetWindowSolnMetaData(const std::string &wname,std::vector<SolnMetaData> &smdv,int verblevel)
{
  int na = 0;
  std::string atts;
  COM_get_dataitems(wname.c_str(),&na,atts);
  int count = 0;
  if(na > 0){
    std::istringstream Istr(atts);
    std::string aname;
    while(Istr >> aname) {
      std::string rstring(wname+"."+aname);
      char loc;
      int ncomp;
      COM_Type ta;
      std::string unit;
      COM_get_dataitem(rstring.c_str(),&loc,&ta,&ncomp,&unit);
      if( (loc == 'E' || loc == 'e' || loc == 'n' || loc == 'N') && ta == COM_DOUBLE)
      	count++;
    }
  }
  smdv.resize(count);
  count = 0;
  if(na > 0) {
    std::istringstream Istr(atts);
    std::string aname;
    while(Istr >> aname){
      std::string rstring(wname+"."+aname);
      char loc;
      int ncomp;
      COM_Type ta;
      std::string unit;
      COM_get_dataitem(rstring.c_str(),&loc,&ta,&ncomp,&unit);
      if( (loc == 'E' || loc == 'e' || loc == 'n' || loc == 'N') && ta == COM_DOUBLE) {
					//	SolnMetaData smd;
					smdv[count].loc   = loc;
					smdv[count].name  = aname;
					smdv[count].unit  = unit;
					smdv[count].ncomp = ncomp;
					//	smdv[count].push_back(smd);
					count++;
      }
    }
  }
}

void TRAIL_GetWindowSolnData(const std::string &wname,
			    std::vector<std::vector<std::vector<double> > > &soln_data,
			    std::vector<SolnMetaData> &smdv,int verblevel)
{


	int nAttrs;
	char *attrStr;
	COM_get_dataitems(wname.c_str(), &nAttrs, &attrStr);
	std::cout << "Number of attributes: " << nAttrs  << std::endl;
	std::cout << "Attribute string: " 		<< attrStr << std::endl;

  std::vector<int> pane_ids;
  COM_get_panes(wname.c_str(),pane_ids);

  if( soln_data.size( ) != pane_ids.size( ) )
  	soln_data.resize(pane_ids.size() );

  std::vector<int>::iterator pii = pane_ids.begin();
  while(pii != pane_ids.end()) {

    int pane_index = pii - pane_ids.begin();
    int pane_id = *pii++;

    if(verblevel > 0) {
      std::cout << "TRAIL_GetWindowSolnData::Pane id: " << pane_id << std::endl
								<< "TRAIL_GetWindowSolnData::Number of atts: "
								<< smdv.size() << std::endl;
    }

    soln_data[pane_index].resize(smdv.size());

    //BugFix: update smdi iterator for each pane (gzagaris)
    std::vector<SolnMetaData>::iterator smdi = smdv.begin();
    while(smdi != smdv.end()) {

    	std::string name 		= wname + "." + smdi ->name;
//    	std::string name	   = wname + ".siVel";
      unsigned int aindex = smdi - smdv.begin();
      int array_size = 0;

      COM_get_size( name.c_str(), pane_id, &array_size);
      if(verblevel)
      	std::cout << smdi->name << " size = " << array_size << std::endl;

      const double *solver_array = NULL;

      COM_get_array_const( name.c_str( ), pane_id, &solver_array);

//      if( array_size > 0 )
//      {
//      	assert( solver_array != NULL );
//      }

      if( solver_array ) {

      	if(verblevel)
					std::cout << "TRAIL_GetWindowSolnData::Copying " << array_size*smdi->ncomp << " doubles into " << smdi->name << " array." << std::endl;

				soln_data[pane_index][aindex].resize(array_size*smdi->ncomp);
				std::memcpy(&soln_data[pane_index][aindex][0],solver_array,sizeof(double)*array_size*smdi->ncomp);
				std::vector<double>::iterator ai = soln_data[pane_index][aindex].begin();
      }
      else {

      	if(verblevel)
					std::cout << "TRAIL_GetWindowSolnData::No " << smdi->name << " array found." << std::endl;

				soln_data[pane_index][aindex].resize(0);
      }

      smdi++;

    } // while smdi != smdv.end()

  } // For all panes

}

// Read the Rocin control file fname into he window wname on
// processors in comm (or MPI_COMM_NULL for 1 proc) - read only
// panes with matching bcflags (or all panes if bcflag.empty(),
// and if(apply_disp) then apply the uhat displacement field,
// get all attributes if(all) or only the mesh if(!all),
// and read ghost info or not.
void
TRAIL_File2Window( const std::string &fname,
		  const std::string &wname,
		  std::vector<int>  &bcflags,
		  MPI_Comm           comm,
		  bool               apply_disp,
		  bool               all,
		  bool               with_ghost)
{

  COM_new_window( wname.c_str());
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimIN, "TRAILIN");
  int IN_read;
  IN_read = COM_get_function_handle( "TRAILIN.read_by_control_file");


  std::string bufwin("bufwin");
  COM_call_function( IN_read, fname.c_str(), bufwin.c_str(), &comm);
  int IN_obtain = COM_get_function_handle( "TRAILIN.obtain_dataitem");

  // If bcflags specified, keep only panes with matching bcflags
  if(!bcflags.empty()){
    int bcflag = COM_get_dataitem_handle((bufwin+".bcflag").c_str());
    if (bcflag > 0) {
      COM_call_function( IN_obtain, &bcflag, &bcflag);
      int npanes, *pane_ids;
      COM_get_panes( bufwin.c_str(), &npanes, &pane_ids);

      // remove panes with bcflag not found in bcflags
      for ( int i=0; i<npanes; ++i) {
	int *flag;
	COM_get_array( (bufwin+".bcflag").c_str(), pane_ids[i], &flag);
	if ( flag==NULL )
	  COM_delete_pane( bufwin.c_str(),pane_ids[i]);
	bool delite = true;
	std::vector<int>::iterator bcfi = bcflags.begin();
	while(bcfi != bcflags.end() && delite)
	  if(*bcfi++ == *flag)
	    delite = false;
	if(delite)
	  COM_delete_pane( bufwin.c_str(), pane_ids[i]);
      }
      COM_free_buffer( &pane_ids);
    }
  }
  if(apply_disp){
    // This is NOT correct for problems with regression.
    int disp_hndl = COM_get_dataitem_handle((bufwin+".uhat").c_str());
    if(disp_hndl > 0){
      COM_call_function( IN_obtain, &disp_hndl, &disp_hndl);
      COM_LOAD_MODULE_STATIC_DYNAMIC( Simpal, "BLAS");
      int add = COM_get_function_handle( "BLAS.add");
      COM_call_function(IN_obtain,&disp_hndl,&disp_hndl);
      int nc_hndl = COM_get_dataitem_handle( bufwin + ".nc");
      COM_call_function( add, &disp_hndl, &nc_hndl, &nc_hndl);
      COM_UNLOAD_MODULE_STATIC_DYNAMIC(Simpal,"BLAS");
    }
  }
  // Read in the mesh.
  int buf_atts;
  if(all)
    buf_atts = COM_get_dataitem_handle((bufwin+".all").c_str());
  else
    buf_atts = COM_get_dataitem_handle((bufwin+".mesh").c_str());
  COM_call_function( IN_obtain, &buf_atts, &buf_atts);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC( SimIN, "TRAILIN");

  if(!all)
    // Remove all attributes except for the mesh
    COM_delete_dataitem(  (bufwin+".data").c_str());

  // Change the memory layout to contiguous.
  if(all)
    COM_clone_dataitem( (wname+".all").c_str(),
			 (bufwin+".all").c_str(), (with_ghost ? 1 : 0));
  else
    COM_clone_dataitem( (wname+".mesh").c_str(),
			 (bufwin+".mesh").c_str(), (with_ghost ? 1 : 0));

  COM_delete_window( bufwin.c_str());
  COM_window_init_done(wname);
}


// DO NOT CALL IN PARALLEL: This function processes files for every processor
void
TRAIL_MergeRocinFiles(const std::string &srcname,
		     const std::string &trgname,
		     const std::string &path,
		     double t,unsigned int np,
		     std::ostream *ouf)
{
  std::string homedir(TRAIL_CWD());
  std::string timestring(TRAIL_TimeString(t));
  std::string filepre(srcname+"_"+timestring+"_");
  std::ostringstream Ofstr;
  unsigned int id = 0;
  if(ouf)
    *ouf << "TRAIL_MergeRocinFiles: Entry" << std::endl;
  TRAIL_CD(path,ouf);
  if(ouf)
    *ouf << "Searching for files with prefix: " << filepre << std::endl;
  while(id <= np){
    std::ifstream Inf;
    std::ostringstream Ostr;
    Ostr << filepre << std::setw(5) << std::setfill('0')
	 << id++;
    std::string filename(Ostr.str()+"_in.txt");
    Inf.open(filename.c_str());
    if(Inf){
      Ofstr << Inf.rdbuf() << std::endl;
      Inf.close();
      unlink(filename.c_str());
    }
  }
  std::ofstream Ouf;
  Ouf.open((trgname+"_in_"+timestring+".txt").c_str());
  Ouf << Ofstr.str();
  Ouf.close();
  TRAIL_CD(homedir,ouf);
  if(ouf)
    *ouf << "TRAIL_MergeRocinFiles: Exit" << std::endl;
}

void
TRAIL_CreateRobustFC(const std::string &wname,const std::string &path)
{
  std::ofstream Ouf;
  Ouf.open((path+"/"+wname+".fc").c_str());
  Ouf << "0.5        0.6        3 0.17365" << std::endl
      << "1.3962634 0.314159265 3" << std::endl
      << "0.0        0.5        3" << std::endl
      << "6 1 1 0" << std::endl
      << "2" << std::endl;
  Ouf.close();
}
void
TRAIL_CreateRobustFC_old(const std::string &wname,const std::string &path)
{
  std::ofstream Ouf;
  Ouf.open((path+"/"+wname+".fc").c_str());
  Ouf << "0.76604444  0.98480775 3 0.98480775" << std::endl
      << "1.3962634 0.314159265  3" << std::endl
      << "0.17364818  0.96592583 3" << std::endl
      << "6 1 1 0" << std::endl
      << "2" << std::endl;
  Ouf.close();
}

// serial function:  **do not call on multiple procs**
// Reads inputs for the source and target windows, whos names are
// spec'd by src,trg from the srcpath and trgpath, and writes the common
// refinement to the destpath.
void
TRAIL_AutoSurfer(const std::string &src, const std::string &trg,
		const std::string &srcpath, const std::string &trgpath,
		const std::string &destpath, double t, MPI_Comm comm,
		std::ostream *ouf)
{
  std::string timestring(TRAIL_TimeString(t));
  std::string srcfile(src + "_in_" + timestring + ".txt");
  std::string trgfile(trg + "_in_" + timestring + ".txt");
  std::string trailwin(src+"_trail");
  std::string homedir(TRAIL_CWD());
  std::string format("HDF");
  std::string newpath;
  COM_LOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");
  int RFC_overlay = COM_get_function_handle( "RFC.overlay");
  int RFC_write   = COM_get_function_handle( "RFC.write_overlay");
  TRAIL_CreateRobustFC(trg,destpath);
  TRAIL_CreateRobustFC(trailwin,destpath);
  COM_set_default_communicator(MPI_COMM_NULL);
  std::vector<int> bcflags(3);
  bcflags[0] = 0;
  bcflags[1] = 1;
  bcflags[2] = 2;
  TRAIL_CD(srcpath,ouf);
  if(ouf)
    *ouf << "TRAIL_AutoSurfer:: homedir = " << homedir << std::endl
	 << "TRAIL_AutoSurfer:: CWD = " << TRAIL_CWD() << std::endl
	 << "TRAIL_AutoSurfer:: Creating common refinement." << std::endl
	 << "TRAIL_AutoSurfer:: Reading in source surface from " << srcpath
	 << std::endl;
  if(ouf)
    std::cout << "Roctrail> Reading in source surface from "
	      << srcpath << "(test of patience)" << std::endl;
  TRAIL_File2Window(srcfile,trailwin,bcflags,MPI_COMM_NULL,false,false,false);
  newpath.assign(homedir+"/"+trgpath);
  TRAIL_CD(newpath,ouf);
  if(ouf)
    *ouf << "TRAIL_AutoSurfer: Reading target surface from " << newpath
	 << std::endl;
  if(ouf)
    std::cout << "Roctrail> Reading target surface from " << newpath
	      << std::endl;
  TRAIL_File2Window(srcfile,trg,bcflags,MPI_COMM_NULL,false,false,false);
  newpath.assign(homedir+"/"+destpath);
  TRAIL_CD(newpath,ouf);
  int src_mesh    = COM_get_dataitem_handle( (trailwin+".mesh").c_str());
  int trg_mesh    = COM_get_dataitem_handle( (trg+".mesh").c_str());
  if(ouf){
    *ouf << "TRAIL_AutoSurfer: Calling Rocface Overlay...."
	 << std::endl;
    std::cout << "Roctrail> Rocface performing overlay..." << std::endl;
  }
  COM_call_function( RFC_overlay, &src_mesh, &trg_mesh);
  if(ouf){
    std::cout << "Roctrail> Overlay complete, writing overlay to " << newpath
	      << std::endl;
    *ouf << "TRAIL_AutoSurfer: Writing overlay to " << newpath << "..."
	 << std::endl;
  }
  COM_call_function( RFC_write, &src_mesh, &trg_mesh,
		     trailwin.c_str(), trg.c_str(), format.c_str());
  TRAIL_CD(homedir,ouf);
  if(ouf)
    *ouf << "TRAIL_AutoSurfer: Done. CWD = " << TRAIL_CWD() << std::endl;
  COM_delete_window(trg);
  COM_delete_window(trailwin);
  COM_set_default_communicator(comm);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");
}

void TRAIL_FixRocstarData(const std::string &wname,std::ostream *ouf = NULL);
void TRAIL_ExtractSurf0(const std::string &srcwin,
		       const std::string &trgwin,
		       std::ostream *ouf = NULL);

bool
TRAIL_TransferSurfDataFILE
(
 const std::string &src,
 const std::string &trg,
 const std::string &dest,
 const std::string &srcpath,
 const std::string &trgpath,
 const std::string &destpath,
 const std::string &crpath,
 double t,unsigned int id,
 MPI_Comm comm,
 std::ostream *ouf)
{
  std::string timestring(TRAIL_TimeString(t));
  std::string suffix("_in_"+timestring+".txt");
  std::string srcfile(src + suffix);
  std::string trgfile(trg + suffix);
  std::string r_trg(dest);
  std::string trailwin(src+"_trail");
  std::string homedir(TRAIL_CWD());
  std::string format("HDF");
  std::string newpath;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::vector<int> bcflags(3);
  bcflags[0] = 0;
  bcflags[1] = 1;
  bcflags[2] = 2;
  COM_LOAD_MODULE_STATIC_DYNAMIC( SurfX, "RFC");
  int RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
  int RFC_read = COM_get_function_handle( "RFC.read_overlay");
  if(ouf)
    *ouf << "Reading source window (" << trailwin << ") from "
	 << srcpath << "/" << srcfile << "...";
  TRAIL_CD(srcpath,ouf);
  TRAIL_File2Window(srcfile,trailwin,bcflags,comm,false,true,false);
  newpath.assign(homedir+"/"+trgpath);
  MPI_Barrier(comm);
  if(ouf)
    *ouf << "done" << std::endl
	 << "Reading target window from " << newpath
	 << "/" << trgfile << "...";
  TRAIL_CD(newpath,ouf);
  TRAIL_File2Window(trgfile,trg,bcflags,comm,false,true,false);
  MPI_Barrier(comm);
  if(ouf)
    *ouf << "done" << std::endl
	 << "Reading destination window " << r_trg << " from "
	 << trgfile << "...";
  TRAIL_File2Window(trgfile,r_trg,bcflags,comm,false,true,true);
  MPI_Barrier(comm);
  newpath.assign(homedir+"/"+crpath);
  if(ouf)
    *ouf << "done" << std::endl
	 << "Reading mesh overlay from  " << newpath << "..." << std::endl;
  TRAIL_CD(newpath,ouf);
  MPI_Barrier(comm);
  int srcmesh = COM_get_dataitem_handle( (trailwin+".mesh").c_str());
  int trgmesh = COM_get_dataitem_handle( (trg+".mesh").c_str());
  std::vector<int> pane_id;
  if(ouf)
    *ouf << "TRAIL_AutoSurfer: Reading mesh overlay for all surfaces."
	 << "TRAIL_AutoSurfer: CR DIR: " << TRAIL_CWD() << std::endl;
  COM_call_function( RFC_read, &srcmesh, &trgmesh, &comm,
		     trailwin.c_str(),trg.c_str(),format.c_str());
  if(ouf)
    *ouf << "TRAIL_AutoSurfer: Beginning transfer for all surfaces."
	 << std::endl;
  MPI_Barrier(comm);
  if(ouf)
    *ouf << "done" << std::endl
	 << "Transferring data ..." << std::endl;
  int num_attributes;
  std::string names;
  COM_get_dataitems( trailwin.c_str(),&num_attributes,names);
  char loc;
  COM_Type comtype;
  int ncomp;
  std::string unit;
  std::istringstream Istr(names);
  for(int i = 0;i < num_attributes;i++){
    std::string aname;
    Istr >> aname;
    COM_get_dataitem(trailwin+"."+aname,&loc,&comtype,&ncomp,&unit);
    if((loc == 'e' || loc == 'n') && comtype == COM_DOUBLE){
      if(!rank)
	if(ouf)
	  std::cout << "Roctrail> Transferring attribute: " << aname << " on "
		    << (loc == 'e' ? "elements" : "nodes") << "."
		    << std::endl;
      COM_resize_array((trailwin+"."+aname).c_str());
      COM_new_dataitem((trg+"."+aname).c_str(),(char)loc,
			COM_DOUBLE,(int)ncomp,unit.c_str());
      COM_new_dataitem((r_trg+"."+aname).c_str(),(char)loc,
			COM_DOUBLE,(int)ncomp,unit.c_str());
      COM_resize_array((r_trg+"."+aname).c_str());
      COM_resize_array((trg+"."+aname).c_str());
      int src_ahdl  = COM_get_dataitem_handle((trailwin+"."+aname).c_str());
      int trg_ahdl  = COM_get_dataitem_handle((trg+"."+aname).c_str());
      COM_call_function( RFC_transfer, &src_ahdl, &trg_ahdl);
      int *srcpane_ids;
      int npanes;
      COM_get_panes( trg.c_str(), &npanes, &srcpane_ids);
      pane_id.resize(npanes);
      for(int i = 0;i < npanes;i++)
	pane_id[i] = srcpane_ids[i];
      // These are no longer necessary as we've duped the info into
      // a locally allocated array
      COM_free_buffer( &srcpane_ids);
      for(int p = 0;p < npanes;p++){
	void *src_ptr = NULL;
	int src_std = 0;
	int src_cap = 0;
	void *trg_ptr = NULL;
	int trg_std = 0;
	int trg_cap = 0;
	COM_get_array((trg+"."+aname).c_str(),pane_id[p],
		      &src_ptr,&src_std,&src_cap);
	COM_get_array((r_trg+"."+aname).c_str(),pane_id[p],
		      &trg_ptr,&trg_std,&trg_cap);
	if(src_ptr && trg_ptr && (trg_std*trg_cap >= src_std*src_cap)){
	  if(ouf)
	    *ouf << "TRAIL_AutoSurfer: Transferred " << aname << "(" << src_std
		 << "," << src_cap << ") to " << aname << "(" << trg_std
		 << "," << trg_cap << ")" << std::endl;
	  memcpy(trg_ptr,src_ptr,sizeof(double)*src_std*src_cap);
	}
	else
	  if(ouf)
	    *ouf << "TRAIL_AutoSurfer: WARNING: non matching sizes for "
		 << aname << " on pane " << pane_id[p] << "."
		 << std::endl
		 << "TRAIL_AutoSurfer: src(" << src_std << "," << src_cap
		 << ") trg(" << trg_std << "," << trg_cap << ")"
		 << std::endl;
      }
    }
  }
  MPI_Barrier(comm);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC(SurfX,"RFC");
  newpath.assign(homedir+"/"+destpath);
  TRAIL_CD(newpath,ouf);
  MPI_Barrier(comm);
  if(ouf){
    *ouf << "Transfer is complete, massaging and writing results to "
	 << newpath << "/" << r_trg << std::endl;
    if(!rank)
      std::cout << "Roctrail> Transfer complete, writing new surface"
		<< std::endl;
  }
  COM_delete_window(trailwin);
  COM_delete_window(trg);

  // FIXME - Temporarily done here as a quick fix
  TRAIL_FixRocstarData(r_trg,ouf);
  TRAIL_ExtractSurf0(r_trg,"surf0",ouf);
  TRAIL_WriteWindow("surf0",".","surf0",".",t,id,comm);
  // FIXME

  TRAIL_WriteWindow(r_trg,".",r_trg,".",t,id,comm);
  TRAIL_CD(homedir,ouf);
  MPI_Barrier(comm);
  COM_delete_window(r_trg);
  if(ouf){
    *ouf << "Results written." << std::endl
	 << "Surface data transfer is complete." << std::endl;
    if(!rank)
      std::cout << "Roctrail> New surface written" << std::endl;
  }
  MPI_Barrier(comm);
  return(true);
}

void
TRAIL_WriteRocinControl(std::vector<int> &pane_id,const std::string &pre,
		       int rank)
{
  std::ofstream Ouf;
  std::string controlfilename(pre + "_in.txt");
  Ouf.open(controlfilename.c_str());
  Ouf << "@Proc: " << rank << std::endl
      << "@Files: " << pre << ".hdf" << std::endl;
  Ouf.clear();
  Ouf << "@Panes: ";
  std::vector<int>::iterator pii = pane_id.begin();
  while(pii != pane_id.end())
    Ouf << *pii++ << " ";
  Ouf << std::endl;
  Ouf.close();
}

bool
TRAIL_WriteWindow(const std::string &wname,const std::string &winpath,
		 const std::string &cntl_name,const std::string &cntl_path,
		 double t,unsigned int id,MPI_Comm comm,std::ostream *ouf)
{

  std::string timestring(TRAIL_TimeString(t));
  std::string homedir(TRAIL_CWD());
  int rank = 0;
  int nproc = 1;
  if(comm != MPI_COMM_NULL){
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nproc);
  }
  int *srcpane_ids;
  int npanes;
  std::vector<int> pane_id;
  if(ouf)
    *ouf << "TRAIL_WriteWindow: Entry" << std::endl;
  COM_get_panes( wname.c_str(), &npanes, &srcpane_ids);
  pane_id.resize(npanes);
  for(int i = 0;i < npanes;i++)
    pane_id[i] = srcpane_ids[i];
  COM_free_buffer( &srcpane_ids);


  COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT,"Rocout");
  int OUT_set_option = COM_get_function_handle( "Rocout.set_option");
  std::string rankstr("0");
  COM_call_function( OUT_set_option, "rankwidth", rankstr.c_str());
  int whand = COM_get_function_handle("Rocout.write_dataitem");


  int all = COM_get_dataitem_handle((wname+".all"));
  std::ostringstream Ostr;
  Ostr << wname << "_" << timestring << "_" << std::setw(5)
       << std::setfill('0') << id;
  TRAIL_CD(winpath,ouf);
  COM_call_function(whand,Ostr.str().c_str(),&all,
		    wname.c_str(),timestring.c_str());
  COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimOUT,"Rocout");

  TRAIL_WriteRocinControl(pane_id,Ostr.str(),rank);
  if(ouf)
    *ouf << "TRAIL_WriteWindow: Wrote window " << wname << " id("
	 << id << ") to " <<  Ostr.str() << ".hdf" << std::endl
	 << "TRAIL_WriteWindow: Merging Rocin control files."
	 << std::endl;
  if(comm != MPI_COMM_NULL)
    MPI_Barrier(comm);
  // Merge Rocin control files
  if(!rank)
    TRAIL_MergeRocinFiles(wname,wname,".",t,nproc,ouf);
  TRAIL_CD(homedir,ouf);
  if(comm != MPI_COMM_NULL)
    MPI_Barrier(comm);
  if(ouf)
    *ouf << "TRAIL_WriteWindow: Exit" << std::endl;
  return(true);
}

// Specify Rocout.remesh, time, and solver name
double
TRAIL_FindSourceTime(const std::string &dirpre,
		    double t,
		    const std::string &relpath)
{
  double targ_time = -1;
  Directory sd(relpath);
  Directory::iterator dfi = sd.begin();
  while(dfi != sd.end()){
    std::string cdir(*dfi++);
    std::string::size_type x = cdir.find(dirpre);
    if(x != std::string::npos){
      double tt = TRAIL_TimeString(cdir.substr(dirpre.size(),9));
      if(tt < t && tt > targ_time)
	targ_time = tt;

    }
  }
  return(targ_time);
}

void
TRAIL_FixRocstarData(const std::string &wname,
		    std::ostream *ouf)
{
  int *srcpane_ids;
  int npanes;
  std::vector<int> pane_id;
  COM_get_panes( wname.c_str(), &npanes, &srcpane_ids);
  pane_id.resize(npanes);
  for(int i = 0;i < npanes;i++)
    pane_id[i] = srcpane_ids[i];
  COM_free_buffer( &srcpane_ids);
  for(int i = 0;i < npanes;i++){
    double *nc_t0 = NULL;
    double *nc    = NULL;
    double *mdot  = NULL;
    int    *bflag = NULL;
    int    stride1 = 0;
    int    stride2 = 0;
    int    cap1    = 0;
    int    cap2    = 0;
    COM_get_array((wname+".nc_t0").c_str(),pane_id[i],
		  &nc_t0,&stride1,&cap1);
    COM_get_array((wname+".nc").c_str(),pane_id[i],
		  &nc,&stride2,&cap2);
    if(nc_t0){
      if(ouf)
	*ouf << "TRAIL_FixRocstarData: Fixing nc_t0 for pane(" << pane_id[i]
	     << ")" << std::endl;
      for(int c = 0;c < cap1;c++){
	unsigned int index = c*stride1;
	if(nc_t0[index]   == 0.0 &&
	   nc_t0[index+1] == 0.0 &&
	   nc_t0[index+2] == 0.0){
	  nc_t0[index]   = nc[index];
	  nc_t0[index+1] = nc[index+1];
	  nc_t0[index+2] = nc[index+2];
	}
      }
    }

    COM_get_array((wname+".bflag").c_str(),pane_id[i],
		  &bflag,&stride1,&cap1);
    COM_get_array((wname+".mdot").c_str(),pane_id[i],
		  &mdot,&stride2,&cap2);
    if(!bflag && mdot){
      COM_resize_array((wname+".bflag").c_str());
      COM_get_array((wname+".bflag").c_str(),pane_id[i],
		    &bflag,&stride1,&cap1);
    }
    if(bflag && mdot){
      if(ouf)
	*ouf << "TRAIL_FixRocstarData: Fixing bflag for pane(" << pane_id[i]
	     << ")" << std::endl;
      for(int c = 0;c < cap2;c++){
	if(mdot[c] > 0.0)
	  bflag[c] = 1;
      }
    }
  }
  COM_window_init_done(wname);
}

void
TRAIL_ExtractSurf0(const std::string &srcwin,
		  const std::string &trgwin,
		  std::ostream *ouf)
{
  COM_new_window(trgwin);
  COM_clone_dataitem( (trgwin+".mesh").c_str(),
		       (srcwin+".mesh").c_str(),1);
  COM_clone_dataitem( (trgwin+".bcflag").c_str(),
		       (srcwin+".bcflag").c_str(),1);
  int *srcpane_ids;
  int npanes;
  std::vector<int> pane_id;
  COM_get_panes( trgwin.c_str(), &npanes, &srcpane_ids);
  pane_id.resize(npanes);
  for(int i = 0;i < npanes;i++)
    pane_id[i] = srcpane_ids[i];
  COM_free_buffer( &srcpane_ids);
  for(int i = 0;i < npanes;i++){
    double *nc_t0 = NULL;
    double *nc    = NULL;
    int stride1   = 0;
    int stride2   = 0;
    int cap1      = 0;
    int cap2      = 0;
    COM_get_array((srcwin+".nc_t0").c_str(),pane_id[i],
		  &nc_t0,&stride1,&cap1);
    COM_get_array((trgwin+".nc").c_str(),pane_id[i],
		  &nc,&stride2,&cap2);
    if(nc && nc_t0)
      memcpy(nc,nc_t0,sizeof(double)*stride2*cap2);
  }
  COM_window_init_done(trgwin);
}

void
TRAIL_RocmopSmooth(GEM_Partition &gp, unsigned int niter)
{
  std::string wname("smooth_win");
  unsigned int nnodes = gp._nc.size()/3;
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RocmopSmooth: This partition has " << nnodes
	     << " nodes." << std::endl;
  MPI_Barrier(gp._comm);
  std::vector<double> disp(nnodes*3,0.0);
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RocmopSmooth: Creating window." << std::endl;
  MPI_Barrier(gp._comm);
  COM_new_window(wname);
  gp.PopulateVolumeWindow(wname);
  COM_new_dataitem((wname+".disp"),'n',COM_DOUBLE,3,"m");
  COM_set_array((wname+".disp"),gp.pane_id,&disp[0]);
  COM_window_init_done(wname);
  int meshhandle = COM_get_dataitem_handle((wname+".pmesh"));
  int disphandle = COM_get_dataitem_handle((wname+".disp"));
  MPI_Barrier(gp._comm);
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RocmopSmooth: Loading Rocmop." << std::endl;
  MPI_Barrier(gp._comm);
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocmop,"MOP");
  int smoothhandle = COM_get_function_handle("MOP.smooth");
  int handleOption = COM_get_function_handle("MOP.set_value");
  MPI_Barrier(gp._comm);
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RocmopSmooth: Setting Rocmop options." << std::endl;
  MPI_Barrier(gp._comm);
  int inverted = 1;
  COM_call_function(handleOption,2,"inverted",&inverted) ;
  MPI_Barrier(gp._comm);
  unsigned int count = 0;
  while(count++ < niter){
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RocmopSmooth: Calling Rocmop's smoothing function.."
	       << std::endl;
    MPI_Barrier(gp._comm);
    COM_call_function(smoothhandle,2,&meshhandle,&disphandle);
    MPI_Barrier(gp._comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RocmopSmooth: Smoothing done. Applying displacements."
	       << std::endl;
    unsigned int n = 0;
    while(n <3*nnodes){
      gp._nc[n] += disp[n];
      disp[n] = 0.0;
      n++;
    }
    MPI_Barrier(gp._comm);
    if(gp._debug && gp._out)
      *gp._out << "TRAIL_RocmopSmooth: Applying displacements done."
	       << std::endl;
    MPI_Barrier(gp._comm);
  }
  COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocmop,"MOP");
  COM_delete_window(wname);
  MPI_Barrier(gp._comm);
  if(gp._debug && gp._out)
    *gp._out << "TRAIL_RocmopSmooth: All done." << std::endl;
  MPI_Barrier(gp._comm);
}

void
TRAIL_RocpropSmoothSurf(double *nc,unsigned int nnodes,
		       unsigned int *ec,unsigned int nel,
		       unsigned int *cnstr_type,
		       unsigned int niter)
{
  std::string wname("smooth_win");
  COM_new_window(wname,MPI_COMM_SELF);
  COM_new_dataitem((wname+".speed"),'e',COM_DOUBLE,1,"m/s");
  COM_new_dataitem((wname+".offsets"),'n',COM_DOUBLE,3,"m");
  COM_new_dataitem((wname+".cnstr_type"),'e',COM_INT,1,"");


  for(unsigned int i = 0; i < nel*3;i++) ec[i]++;
  std::vector<double> offsets;
  offsets.resize(nnodes*3);
  std::vector<double> speed;
  speed.resize(nel,0.0);

  COM_set_size((wname+".nc"),102,nnodes);
  COM_set_array((wname+".nc"),102,nc,3);
  COM_set_size((wname+".:t3:real"),102,nel);
  COM_set_array((wname+".:t3:real"),102,ec,3);
  //  COM_set_size((wname+".speed"),102,nel);
  COM_set_array((wname+".speed"),102,&speed[0],3);
  //  COM_set_size((wname+".offsets"),102,nnodes,3);
  COM_set_array((wname+".offsets"),102,&offsets[0],3);
  //  COM_set_size((wname+".cnstr_type"),102,nel);
  COM_set_array((wname+".cnstr_type"),102,&cnstr_type[0],1);

  int cnstrtype  = COM_get_dataitem_handle((wname+".cnstr_type"));
  int offset_hdl = COM_get_dataitem_handle((wname+".offsets"));
  int speed_hdl  = COM_get_dataitem_handle((wname+".speed"));

  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocprop, "PROP");
  int setopt = COM_get_function_handle( "PROP.set_option");
  COM_call_function( setopt, "method", "fo");
  std::ostringstream Ostr;
  Ostr << niter << std::ends;
  COM_call_function( setopt, "rediter", Ostr.str().c_str());
  COM_call_function( setopt, "fangle", "35");
  int init = COM_get_function_handle("PROP.initialize");
  int prop = COM_get_function_handle("PROP.propagate");
  int setcnstr = COM_get_function_handle("PROP.set_constraints");

  int pmesh_hdl = COM_get_dataitem_handle((wname+".mesh"));
  COM_call_function(init, &pmesh_hdl);

  COM_call_function(setcnstr,&cnstrtype);

  double dt = 0.0;

  COM_call_function(prop,&pmesh_hdl,&speed_hdl,&dt,&offset_hdl);

  for(unsigned int i = 0; i < nnodes*3; i++)
    nc[i] += offsets[i];

  COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocprop,"PROP");
  COM_delete_window(wname);
  for(unsigned int i = 0; i < nel*3;i++) ec[i]--;
}


