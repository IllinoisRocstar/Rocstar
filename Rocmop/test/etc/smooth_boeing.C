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
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <bitset>
#include "../../Roccom/External/cgnslib_2.4/cgnslib.h"
#include "roccom.h"
#include "mapbasic.h"
#include "Rocblas.h"
#include "Rocmop_1_1.h"
#include "mopbasic.h"

using namespace std;

COM_EXTERN_MODULE( Rocmop);
COM_EXTERN_MODULE( Rocout);
COM_EXTERN_MODULE( Rocblas);

int main(int argc, char*argv[])
{

  MPI_Init( &argc, &argv);

  if (argc < 3){
    cout << "Usage: " << argv[0]
	 << " <fname><niter>" << endl;
    exit(-1);
  }

  COM_init( &argc, &argv);

  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocmop, "MOP");
  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocout, "OUT");

  // Get function handles
  int OUT_write = COM_get_function_handle( "OUT.write_attribute");
  int MOP_sb    = COM_get_function_handle( "MOP.smooth_boeing");
  int MOP_set   = COM_get_function_handle( "MOP.set_value");
  int MOP_angles= COM_get_function_handle( "MOP.obtain_extremal_dihedrals");

  string fname(argv[1]);
  int niter = atoi(argv[2]);

  const int nameSize = 33;
  int file_hdl;
  int nbases;
  int nzones;
  int bCellDim;
  int bPhysDim;
  char bName[33];  
  char zName[33];
  int sizes[3];
  SimulationType_t sType;

  const char* _3DUnstructInfo[3] = 
    {"NVertex", "NCell3D", "NBoundVertex"};

  std::cout << "Opening file: " << fname.c_str() << "\n";
  
  // open CGNS file
  if(cg_open(fname.c_str(), MODE_MODIFY, &file_hdl))
    cg_error_exit();

  // read the number of bases
  if(cg_nbases(file_hdl,&nbases))
    cg_error_exit();
  
  if(nbases != 1){
    cout << "Error: .cgns files with multiple bases not supported.\n";
    assert(0);
  }

  // get info about each base
  std::cout << "File opened successfully.\n";

  int i =1;

  if(cg_base_read(file_hdl, i, bName,&bCellDim, &bPhysDim))
    cg_error_exit();

  std::cout << "Checking mesh format.\n";
  if(bCellDim != 3 ||
     bPhysDim != 3){
    cout << "Error: Only volume meshes are supported by this utility.\n";
    exit(-1);
  }

  // Examine zone data:
  if(cg_nzones(file_hdl, i, &nzones))
    cg_error_exit();
  
  ZoneType_t zType;
    
  if( nzones != 1 ){
    cout << "Error: .cgns files with multiple zones not supported by this utility.\n";
    assert(0);
  }

  int j = 1;

  if(cg_zone_type(file_hdl, i, j, &zType))
    cg_error_exit();

  if(zType != Unstructured){
    cout << "Error: only structured meshes supported by this utility.\n";
    assert(0);
  }

  if(cg_zone_read(file_hdl, i, j, zName, &sizes[0]))
    cg_error_exit();

  for(int k = 0; k < 3; ++k)
    cout << _3DUnstructInfo[k] << " = " <<  sizes[k] << "\n";

  if(zType == Structured){
    cout << "Error: this utility does not supoort structured meshes.\n";
    exit(-1);
  }

  // Check for Grids associated w/ this node,zone:
  int nGrids;
  if(cg_ngrids(file_hdl, i, j, &nGrids))
    cg_error_exit();
      
  if( nGrids != 1 ){
    cout << "Erorr: this utility only supports files with 1 grid\n";
    exit(-1);
  }

  int ncoords = 0;
  if(cg_ncoords(file_hdl,i,j,&ncoords))
    cg_error_exit();

  if(ncoords != 3){
    cout << "Error: this utility only supports meshes with 3 coordinate components\n";
    exit(-1);
  }

  int nnodes = sizes[0];
  std::vector<double> coord_data(3*nnodes);

  for(int k=1; k<=ncoords; ++k){
	
    DataType_t dataType;
    char coordName[nameSize];
    if(cg_coord_info(file_hdl, i, j, k, &dataType, &coordName[0]))
      cg_error_exit();

    if(dataType != RealDouble){
      std::cout << "Error: this utility only supports mesh coordinates of CGNS type RealDouble\n";
      exit(-1);
    }

    int ranges[2] = {1,nnodes};
    if(cg_coord_read(file_hdl,i,j,&coordName[0],dataType,&ranges[0],
		     &ranges[1],&coord_data[(k-1)*nnodes]))
      cg_error_exit();
  }

  // Check for Conn tables associated w/ this node,zone:
  int nConns;
  if(cg_nsections(file_hdl, i, j, &nConns))
    cg_error_exit();
      
  char cName[nameSize];
  ElementType_t eType;
  int rmin_elems, rmax_elems, nbnrdy, parent_flag;
  int nodes_pe = 0; // nodes per element

  // Ok, we've spit out enough info, lets actually go smooth
  for(int k=1; k<nConns; ++k){

    if(cg_section_read(file_hdl, i, j, k, &cName[0], &eType, 
		       &rmin_elems, &rmax_elems, &nbnrdy, &parent_flag))
      cg_error_exit();

      /*
       The .cgns documentation says that parent_data "contains
       information on the cell(s) and cell faces(s) sharing the element"

       There appears to be no way to determine the size of the parent_data
       a priori, and I see no cg_ call to determine this information. 
       Since the Boeing meshes I've examined so far do not contain this
       piece of information, I've decided not to write the utility to support
       its presence.
      */
    if(parent_flag){
      cout << "Error: the parent flag is set, and this utility does not "
	   << "support meshes with parent data.\n";
      exit(-1);
    }

    // For the time being, we're only optimizing tet elements
    if(eType == TETRA_4){

      /* 
	 We will create node index mappings in both directions:
	 .cgns -> temp. Roccom window,
	 temp.Roccom window -> .cgns
	 This is overkill, but doesn't change the big O order of the
	 program. I've chosen this method because it is very clear and
	 should make future maintenance very easy.
      */	      

      // Determine the number of elements in the mesh
      int nElems;
      int nodes_pe = 4;
      if(cg_ElementDataSize(file_hdl, i, j, k, &nElems))
	cg_error_exit();
      nElems /= nodes_pe;
      
      // Allocate buffers for the tetrahedral connectivity tables for
      // both the .cgns file and the temp Roccom window
      std::vector<int> conn_data(nElems*nodes_pe,-1);	
      std::vector<int> temp_conn(nElems*nodes_pe,-2);
      int parent_data[1];

      // Allocate a buffer for the mapping the ids in the .cgns file
      // to those in the Roccom window we are going to create
      std::vector<int> cgns_to_roccom(nnodes,-1);

      // Read in the connectivity table to our buffer
      if(cg_elements_read(file_hdl, i, j, k, &conn_data[0], &parent_data[0]))
	cg_error_exit();

      cout << "Found a tetrahedral zone with " << nElems << " elements\n";
      // Make a pass over the connectivity table
      //   1) Determine the number of nodes used by tet elements.
      //   2) Assign new ids to these nodes
      //   3) Fill in the connectivity table for the Roccom window
      int used_nodes = 0;      
      for(int ii=0; ii< nElems; ++ii){
	for(int jj=0; jj< nodes_pe; ++jj){
	  int old_node_ind = conn_data[ii*nodes_pe + jj]-1;
	  if(cgns_to_roccom[old_node_ind] == -1)
	    cgns_to_roccom[old_node_ind] = used_nodes++;
	  temp_conn[ii*nodes_pe + jj] = cgns_to_roccom[old_node_ind]+1;
	}
      }

      // We're finished w/ the connectivity data now, so minimize its footprint
      std::vector<int>().swap(conn_data);

      // Now that we know how many nodes will be in the Roccom window,
      // allocate buffers for its nodal coordinates and for the 
      // mapping of node ids from the Roccom window to the .cgns file
      std::vector<double> temp_coords(used_nodes*3,-3);
      std::vector<int> roccom_to_cgns(used_nodes,-4);

      // Pass over the .cgns coordinate array
      //   1) Fill in the coordinate array for the Roccom window
      //   2) Fill in the 2nd node id mapping
      for(int ii =0; ii<nnodes; ++ii){
	if(cgns_to_roccom[ii] != -1){
	  int node_ind = cgns_to_roccom[ii];
	  temp_coords[node_ind*3] = coord_data[ii];
	  temp_coords[node_ind*3+1] = coord_data[ii+nnodes];
	  temp_coords[node_ind*3+2] = coord_data[ii+2*nnodes];
	}
      }

      // close CGNS file temporarily to decrease mem footprint
      if(cg_close(file_hdl))
	cg_error_exit();

      // Create a new window
      cout << "Initializing smoothing data structures\n";
      COM_new_window("unstr_temp");
      
      // Register nodal coordinates and connectivity table
      COM_set_size("unstr_temp.nc", 1, used_nodes);
      COM_set_array("unstr_temp.nc",1, &temp_coords[0],3);
      COM_set_size("unstr_temp.:T4:", 1, nElems);
      COM_set_array("unstr_temp.:T4:",1, &temp_conn[0],4);

      COM_window_init_done("unstr_temp");    

      // Obtain attribute handles
      int HDL_mesh = COM_get_attribute_handle("unstr_temp.mesh");

      #if 0

      // Write out the mesh in its initial configuration
      cout << "Writing initial mesh configuration\n";
      COM_call_function(OUT_write, "unstr_temp", &HDL_mesh,
			"unstr_temp", "001");
      
      #endif

      double min_angle = 190, max_angle=0;
      COM_call_function(MOP_angles, &HDL_mesh, &min_angle, &max_angle);
      cout << "Initial dihedral angle range = ( " << min_angle << " , " << max_angle << " )\n";
      cout << "Beginning mesh optimizaton\n";

      cout << setiosflags(std::ios::left);
      cout << std::setw(15) << "iteration" << std::setw(20) << "min dihedral " 
	   << std::setw(20) << "max dihedral " << endl;
      cout << std::setw(55) << std::setfill('-') << "" << endl;
      cout << std::setfill(' ');

      // Use a conjugate gradient solver to decrease the required
      // memory footprint
      int wrapper_cg = 1;
      COM_call_function(MOP_set, "wrapper", &wrapper_cg);

      for(int ii =0; ii< niter; ++ii){

	int one = 1;
	COM_call_function( MOP_sb, &HDL_mesh, &one);
	COM_call_function(MOP_angles, &HDL_mesh, &min_angle, &max_angle);

	cout << std::setw(15) << ii+1 << std::setw(20) << min_angle 
	     << std::setw(20) << max_angle << "\n";
      }


      #if 0
      // Writing out the smoothed mesh configuration
      cout << "Writing out smoothed configuration\n";
      COM_call_function(OUT_write, "unstr_temp_smoothed", &HDL_mesh,
			"unstr_temp", "001");
      #endif
      
      // Pass over the .cgns coordinate array
      //   1) Fill in the coordinate array for the Roccom window
      //   2) Fill in the 2nd node id mapping
      for(int ii =0; ii<nnodes; ++ii){
	if(cgns_to_roccom[ii] != -1){
	  int node_ind = cgns_to_roccom[ii];
	  coord_data[ii] = temp_coords[node_ind*3];
	  coord_data[ii+nnodes] = temp_coords[node_ind*3+1];
	  coord_data[ii+2*nnodes] = temp_coords[node_ind*3+2];
	}
      }

      // reopen CGNS file
      if(cg_open(fname.c_str(), MODE_MODIFY, &file_hdl))
	cg_error_exit();
      
      cout << "Writing updated coordinates to .cgns file\n";
      for(int ii= 1; ii<=3; ++ii){
	char coordName[33];
	int coord_ind;
	DataType_t dataType;
	
	if(cg_coord_info(file_hdl, i, j, ii, &dataType, &coordName[0]))
	  cg_error_exit();

	if(cg_coord_write(file_hdl,i,j,RealDouble,coordName, 
			  &coord_data[(ii-1)*nnodes], &coord_ind))
	  cg_error_exit();

	assert(coord_ind == ii);
      }
      
    } // if(eType == TETRA_4)
  }      

  // close CGNS file
  if(cg_close(file_hdl))
    cg_error_exit();
}






