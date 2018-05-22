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
// $Id: build_meshes.C,v 1.4 2008/12/06 08:45:25 mtcampbe Exp $

#include "com.h"
#include "mapbasic.h"
#include "Rocblas.h"
#include <iostream>
#include <sstream>

using namespace std;

COM_EXTERN_MODULE( SurfMap);
COM_EXTERN_MODULE( SimOUT);
COM_EXTERN_MODULE( Rocmop);
COM_EXTERN_MODULE( Simpal);

// build an unstructured tet mesh with 2 partitions
void build_unstr_tet_2();

// build a serial unstructured hex mesh consisting
// of a 2x2x2 stack of elements
void build_unstr_hex();

void build_unstr_pyr();

void build_unstr_prism();

void build_unstr_prism_tet();

void build_unstr_prism_tet_2();

template<class T>
void print_array(T* data,int size, int space){
  for(int i =0; i < size; ++i){
    cout << data[i] << " ";
    if((i+1)%space ==0)
      cout << endl;
  }
}

int main(int argc, char *argv[]) {

  // Initialize Roccom
  COM_init( &argc, &argv);

  // Set Rocom's verbose level pretty high, and
  // turn on profiling.
  COM_set_verbose(111);
  COM_set_profiling(1);

  // Load the required Roccom modules.
  COM_LOAD_MODULE_STATIC_DYNAMIC( SimOUT, "OUT");
  COM_LOAD_MODULE_STATIC_DYNAMIC( Rocmop, "MOP");
  COM_LOAD_MODULE_STATIC_DYNAMIC( Simpal, "BLAS");

  build_unstr_tet_2();
  build_unstr_hex();
  build_unstr_pyr();
  build_unstr_prism();
  build_unstr_prism_tet();
  build_unstr_prism_tet_2();

  COM_finalize();
}

// Predeclaration of functions needed to build multiple
// meshes

// Predeclaration of functions needed by build_unstr_tet_2

// Build 6 tets from the hex consisting of nodes
// a,b,c,d and a+9,b+9,c+9,d+9
void init_tets_from_hex(int *elmt, int a, int b, int c, int d);
// Build a tet mesh by subdividing a 2x2x4 hex mesh
void init_tet_unstr_2_elmts(int* elmts);
// Initialize the mesh coordinates for pane1
void init_crds1_tet_unstr_2(double* crds);
// Initialize the mesh coordinates for pane2
void init_crds2_tet_unstr_2(double* crds);

void build_unstr_tet_2(){

  const int unstr_tet_2_num_nodes = 27;
  const int unstr_tet_2_num_elmts = 48;
  double unstr_tet_2_coors1[unstr_tet_2_num_nodes*3];
  double unstr_tet_2_coors2[unstr_tet_2_num_nodes*3];
  int    unstr_tet_2_elmts1[unstr_tet_2_num_elmts*4];
  int    unstr_tet_2_elmts2[unstr_tet_2_num_elmts*4];
  int    unstr_tet_2_pconn1[12] = {1,2,9,3,6,9,12,15,18,21,24,27};
  int    unstr_tet_2_pconn2[12] = {1,1,9,1,4,7,10,13,16,19,22,25};

  // Create a new window named "unstr"
  COM_new_window("unstr_tet_2");

  // Regist nodal coordinates, connectivity tables,
  // and the pconn for pane 1.
  COM_set_size( "unstr_tet_2.nc", 1, unstr_tet_2_num_nodes);
  COM_set_array( "unstr_tet_2.nc", 1, &unstr_tet_2_coors1[0],3);
  COM_set_size( "unstr_tet_2.:T4:", 1, unstr_tet_2_num_elmts);
  COM_set_array( "unstr_tet_2.:T4:", 1, &unstr_tet_2_elmts1[0],4);
  COM_set_size( "unstr_tet_2.pconn", 1, 12);
  COM_set_array( "unstr_tet_2.pconn", 1, &unstr_tet_2_pconn1[0]);

  // Register nodal coordinates, connectivity tables,
  // and the pconn for pane 2.
  COM_set_size( "unstr_tet_2.nc", 2, unstr_tet_2_num_nodes);
  COM_set_array( "unstr_tet_2.nc", 2, &unstr_tet_2_coors2[0],3);
  COM_set_size( "unstr_tet_2.:T4:", 2, unstr_tet_2_num_elmts);
  COM_set_array( "unstr_tet_2.:T4:", 2, &unstr_tet_2_elmts2[0],4);
  COM_set_size( "unstr_tet_2.pconn", 2, 12);
  COM_set_array( "unstr_tet_2.pconn", 2, &unstr_tet_2_pconn2[0]);

  // Let Roccom know we have finished adding dataitems
  // to the window "unstr_tet_2"
  COM_window_init_done("unstr_tet_2");

  // Calculate coordinates and build connectivity arrays.
  init_crds1_tet_unstr_2(unstr_tet_2_coors1);
  init_crds2_tet_unstr_2(unstr_tet_2_coors2);
  init_tet_unstr_2_elmts(unstr_tet_2_elmts1);
  init_tet_unstr_2_elmts(unstr_tet_2_elmts2);

  // Write the mesh to file.
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int HDL_all = COM_get_dataitem_handle("unstr_tet_2.all");
  COM_call_function( OUT_write, "unstr_tet_2", &HDL_all, 
		     "unstr_tet_2", "001");

}

void init_crds1_tet_unstr_2(double* crds){
  for(int z = 0; z < 3; ++z){
    for(int y = 0; y < 3; ++y){
      for(int x = 0; x < 3; ++x){
	int index = z*9 + y*3 + x;
	crds[index*3]=x;
	crds[index*3+1]=2-y;
	crds[index*3+2]=2-z;
	//	if((index == 13)||(index == 14))
	if(index == 14)
	  crds[index*3] += .33;
      }
    }
  }
}

void init_crds2_tet_unstr_2(double* crds){
  init_crds1_tet_unstr_2(crds);
  for(int i =0; i < 27; ++i){
    crds[3*i] += 2;
    if(i == 12)
      crds[3*i] += .33;
    else if (i == 14)
      crds[3*i] -= .33;
  }
}

void init_tet_unstr_2_elmts(int* elmts){
  init_tets_from_hex(elmts+0,  1,2,4,5);
  init_tets_from_hex(elmts+24, 2,3,5,6);
  init_tets_from_hex(elmts+48, 4,5,7,8);
  init_tets_from_hex(elmts+72, 5,6,8,9);
  init_tets_from_hex(elmts+96, 10,11,13,14);
  init_tets_from_hex(elmts+120,11,12,14,15);
  init_tets_from_hex(elmts+144,13,14,16,17);
  init_tets_from_hex(elmts+168,14,15,17,18);  
}

void init_tets_from_hex(int *elmt, int a, int b, int c, int d){
  int e = a+9, f = b+9, g = c+9, h = d+9;
  elmt[0]  = a; elmt[1]  =b; elmt[2]  = c; elmt[3]  = e;
  elmt[4]  = b; elmt[5]  =d; elmt[6]  = c; elmt[7]  = e;
  elmt[8]  = c; elmt[9]  =e; elmt[10] = d; elmt[11] = g;
  elmt[12] = b; elmt[13] =d; elmt[14] = e; elmt[15] = f;
  elmt[16] = d; elmt[17] =e; elmt[18] = f; elmt[19] = g;
  elmt[20] = d; elmt[21] =g; elmt[22] = f; elmt[23] = h;
}

// predeclaration of functions used in build_unstr_hex
void init_hex(int *elmt, int a, int b, int c, int d);
void init_unstr_hex_elmts(int* elmts);
void init_unstr_hex_coords(double *coords);

// build a serial unstructured hex mesh consisting
// of a 2x2x2 stack of elements
void build_unstr_hex(){

  const int num_nodes = 27;
  const int num_elmts = 8;
  double crds[num_nodes*3];
  int    elmts[num_elmts*8];

  // Create a new window named "unstr"
  COM_new_window("unstr_hex");

  // Regist nodal coordinates, connectivity tables,
  // and the pconn for pane 1.
  COM_set_size( "unstr_hex.nc", 1, num_nodes);
  COM_set_array( "unstr_hex.nc", 1, &crds[0],3);
  COM_set_size( "unstr_hex.:H8:", 1, num_elmts);
  COM_set_array( "unstr_hex.:H8:", 1, &elmts[0],8);

  COM_window_init_done("unstr_hex");

  init_hex(elmts,   1,2,5,4); init_hex(elmts+8, 2,3,6,5);  
  init_hex(elmts+16,4,5,8,7); init_hex(elmts+24,5,6,9,8);  

  init_hex(elmts+32,10,11,14,13); init_hex(elmts+40,11,12,15,14);  
  init_hex(elmts+48,13,14,17,16); init_hex(elmts+56,14,15,18,17);

  init_unstr_hex_coords(&crds[0]);
  
  // Write the mesh to file.
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int HDL_all = COM_get_dataitem_handle("unstr_hex.all");
  COM_call_function( OUT_write, "unstr_hex", &HDL_all, 
		     "unstr_hex", "001");

}

void init_hex(int *elmt, int a, int b, int c, int d){
  elmt[0] = a;   elmt[1] = b;   elmt[2] = c;   elmt[3] = d;
  elmt[4] = a+9; elmt[5] = b+9; elmt[6] = c+9; elmt[7] = d+9;
}

void init_unstr_hex_coords(double *coords){

  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      for(int k=0;k<3;++k){
	coords[3*(9*k+3*j+i)] = (double)i;
	coords[3*(9*k+3*j+i)+1] = 1.0*(double)j;
	coords[3*(9*k+3*j+i)+2] = (double)k;
      }
    }
  }
  
  coords[39] = .73;
  coords[40] = .87;
  coords[41] = 1.27;
}

// predeclaration of functions used in build_unstr_hex
void init_pyr(int *elmt, int a, int b, int c, int d, int e);
void init_unstr_pyr_elmts(int* elmts);
void init_unstr_pyr_coords(double *coords);

// build a serial unstructured hex mesh consisting
// of a 2x2x2 stack of elements
void build_unstr_pyr(){

  const int num_nodes = 9;
  const int num_elmts = 6;
  double crds[num_nodes*3];
  int    elmts[num_elmts*5];

  // Create a new window named "pyr"
  COM_new_window("unstr_pyr");

  // Regist nodal coordinates, connectivity tables,
  // and the pconn for pane 1.
  COM_set_size( "unstr_pyr.nc", 1, num_nodes);
  COM_set_array( "unstr_pyr.nc", 1, &crds[0],3);
  COM_set_size( "unstr_pyr.:P5:", 1, num_elmts);
  COM_set_array( "unstr_pyr.:P5:", 1, &elmts[0],5);

  COM_window_init_done("unstr_pyr");

  init_pyr(elmts,    1,2,3,4,5);
  init_pyr(elmts+5,  1,6,7,2,5);
  init_pyr(elmts+10, 2,7,8,3,5);
  init_pyr(elmts+15, 3,8,9,4,5);
  init_pyr(elmts+20, 4,9,6,1,5);
  init_pyr(elmts+25, 6,9,8,7,5);

  init_unstr_pyr_coords(&crds[0]);
  
  // Write the mesh to file.
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int HDL_all = COM_get_dataitem_handle("unstr_pyr.all");
  COM_call_function( OUT_write, "unstr_pyr", &HDL_all, 
		     "unstr_pyr", "001");
}

void init_pyr(int *elmt, int a, int b, int c, int d, int e){
  elmt[0] = a; elmt[1] = b; elmt[2] = c; elmt[3] = d; elmt[4] = e;
}

void init_unstr_pyr_coords(double *coords){
  coords[0]  = 0.0; coords[1]  = 0.0; coords[2]  = 0.0;
  coords[3]  = 1.0; coords[4]  = 0.0; coords[5]  = 0.0;
  coords[6]  = 1.0; coords[7]  =1.0; coords[8]  = 0.0;
  coords[9]  = 0.0; coords[10] =1.0; coords[11] = 0.0;
  coords[12] = 0.5; coords[13] = 0.5; coords[14] = 0.7;  
  coords[15] = 0.0; coords[16] =  0.0; coords[17] = 1.0;  
  coords[18] = 1.0; coords[19] =  0.0; coords[20] = 1.0;  
  coords[21] = 1.0; coords[22] = 1.0; coords[23] = 1.0;  
  coords[24] = 0.0; coords[25] = 1.0; coords[26] = 1.0;  
}

// predeclaration of functions used in build_unstr_prism
void init_prism(int *elmt, int a, int b, int c, int d, int e, int f);
void init_unstr_prism_elmts(int* elmts);
void init_unstr_prism_coords(double *coords);

void build_unstr_prism(){

  const int num_nodes = 15;
  const int num_elmts = 8;
  double crds[num_nodes*3];
  int    elmts[num_elmts*6];

  // Create a new window named "pyr"
  COM_new_window("unstr_prism");

  // Regist nodal coordinates, connectivity tables,
  // and the pconn for pane 1.
  COM_set_size( "unstr_prism.nc", 1, num_nodes);
  COM_set_array( "unstr_prism.nc", 1, &crds[0],3);
  COM_set_size( "unstr_prism.:W6:", 1, num_elmts);
  COM_set_array( "unstr_prism.:W6:", 1, &elmts[0],6);

  COM_window_init_done("unstr_prism");

  init_prism(elmts, 1,5,4,6,10,9);
  init_prism(elmts+6, 1,2,5,6,7,10);
  init_prism(elmts+12, 2,3,5,7,8,10);
  init_prism(elmts+18, 3,4,5,8,9,10);
  init_prism(elmts+24, 6,10,9,11,15,14);
  init_prism(elmts+30, 6,7,10,11,12,15);
  init_prism(elmts+36, 7,8,10,12,13,15);
  init_prism(elmts+42, 8,9,10,13,14,15);

  init_unstr_prism_coords(&crds[0]);
  
  // Write the mesh to file.
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int HDL_all = COM_get_dataitem_handle("unstr_prism.all");
  COM_call_function( OUT_write, "unstr_prism", &HDL_all, 
		     "unstr_prism", "001");
}

void init_prism(int *elmt, int a, int b, int c, int d, int e, int f){
  elmt[0] = a;   elmt[1] = b;   elmt[2] = c;
  elmt[3] = d; elmt[4] = e; elmt[5] = f;
}

void init_unstr_prism_coords(double *coords){
  coords[0]  = 0.0; coords[1]  = 0.0; coords[2]  = 0.0;
  coords[3]  = 1.0; coords[4]  = 0.0; coords[5]  = 0.0;
  coords[6]  = 1.0; coords[7]  = 1.0; coords[8]  = 0.0;
  coords[9]  = 0.0; coords[10] = 1.0; coords[11] = 0.0;
  coords[12] = 0.5; coords[13] = 0.5; coords[14] = 0.0;  
  coords[15] = 0.0; coords[16] = 0.0; coords[17] = 1.0;  
  coords[18] = 1.0; coords[19] = 0.0; coords[20] = 1.0;  
  coords[21] = 1.0; coords[22] = 1.0; coords[23] = 1.0;  
  coords[24] = 0.0; coords[25] = 1.0; coords[26] = 1.0;  
  coords[27] = 0.34; coords[28] = 0.55; coords[29] = 1.23;  
  coords[30] = 0.0; coords[31] = 0.0; coords[32] = 2.0;
  coords[33] = 1.0; coords[34] = 0.0; coords[35] = 2.0;  
  coords[36] = 1.0; coords[37] = 1.0; coords[38] = 2.0;  
  coords[39] = 0.0; coords[40] = 1.0; coords[41] = 2.0;  
  coords[42] = 0.5; coords[43] = 0.5; coords[44] = 2.0;  
}

// predeclaration of functions used in build_unstr_prism_tet
void init_unstr_prism_tet_elmts(int* elmts);
void init_unstr_prism_tet_coords(double *coords);
void init_tet(int *elmt, int a, int b, int c, int d);

void build_unstr_prism_tet(){

  const int num_nodes = 16;
  const int num_prisms = 8;
  const int num_tets = 4;
  double crds[num_nodes*3];
  int    prism_conn[num_prisms*6];
  int    tet_conn[num_tets*4];

  // Create a new window
  COM_new_window("unstr_prism_tet");

  // Regist nodal coordinates, connectivity tables,
  // and the pconn for pane 1.

  COM_set_size( "unstr_prism_tet.nc", 1, num_nodes);
  COM_set_array( "unstr_prism_tet.nc", 1, &crds[0],3);
  COM_set_size( "unstr_prism_tet.:W6:", 1, num_prisms);
  COM_set_array( "unstr_prism_tet.:W6:", 1, &prism_conn[0],6);
  COM_set_size( "unstr_prism_tet.:T4:", 1, num_tets);
  COM_set_array( "unstr_prism_tet.:T4:", 1, &tet_conn[0],4);

  COM_window_init_done("unstr_prism_tet");

  init_prism(prism_conn, 1,5,4,6,10,9);
  init_prism(prism_conn+6, 1,2,5,6,7,10);
  init_prism(prism_conn+12, 2,3,5,7,8,10);
  init_prism(prism_conn+18, 3,4,5,8,9,10);
  init_prism(prism_conn+24, 6,10,9,11,15,14);
  init_prism(prism_conn+30, 6,7,10,11,12,15);
  init_prism(prism_conn+36, 7,8,10,12,13,15);
  init_prism(prism_conn+42, 8,9,10,13,14,15);

  init_tet(tet_conn,   11,12,15,16);
  init_tet(tet_conn+4, 12,13,15,16);
  init_tet(tet_conn+8, 13,14,15,16);
  init_tet(tet_conn+12,14,11,15,16);

  init_unstr_prism_tet_coords(&crds[0]);
  
  // Write the mesh to file.
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int HDL_all = COM_get_dataitem_handle("unstr_prism_tet.all");
  COM_call_function( OUT_write, "unstr_prism_tet", &HDL_all, 
		     "unstr_prism_tet", "001");
}

void init_tet(int *elmt, int a, int b, int c, int d){
  elmt[0] = a;   elmt[1] = b;   elmt[2] = c;   elmt[3] = d;
}

void init_unstr_prism_tet_coords(double *coords){
  coords[0]  = 0.0; coords[1]  = 0.0; coords[2]  = 0.0; // 1
  coords[3]  = 1.0; coords[4]  = 0.0; coords[5]  = 0.0; // 2
  coords[6]  = 1.0; coords[7]  = 1.0; coords[8]  = 0.0; // 3
  coords[9]  = 0.0; coords[10] = 1.0; coords[11] = 0.0; // 4
  coords[12] = 0.5; coords[13] = 0.5; coords[14] = 0.0; // 5 

  coords[15] = 0.0; coords[16] = 0.0; coords[17] = 1.0; // 6 
  coords[18] = 1.0; coords[19] = 0.0; coords[20] = 1.0; // 7 
  coords[21] = 1.0; coords[22] = 1.0; coords[23] = 1.0; // 8 
  coords[24] = 0.0; coords[25] = 1.0; coords[26] = 1.0; // 9 
  coords[27] = 0.65; coords[28] = 0.35; coords[29] = 1.0; // 10
  
  coords[30] = 0.0; coords[31] = 0.0; coords[32] = 2.0; // 11
  coords[33] = 1.0; coords[34] = 0.0; coords[35] = 2.0; // 12
  coords[36] = 1.0; coords[37] = 1.0; coords[38] = 2.0; // 13
  coords[39] = 0.0; coords[40] = 1.0; coords[41] = 2.0; // 14
  coords[42] = 0.27; coords[43] = 0.29; coords[44] = 2.24; // 15
  
  coords[45] = 0.5; coords[46] = 0.5; coords[47] = 2.5; // 16
}

void init_prism_tet_unit(int* prism_conn, int* tet_conn,
			 int n1, int n2, int n3, int n4,
			 int n5, int n6, int n7, int n8,
			 int n9, int n10, int n11, int n12,
			 int n13, int n14, int n15, int n16){
 
  init_prism(prism_conn, n1,n5,n4,n6,n10,n9);
  init_prism(prism_conn+6, n1,n2,n5,n6,n7,n10);
  init_prism(prism_conn+12, n2,n3,n5,n7,n8,n10);
  init_prism(prism_conn+18, n3,n4,n5,n8,n9,n10);
  init_prism(prism_conn+24, n6,n10,n9,n11,n15,n14);
  init_prism(prism_conn+30, n6,n7,n10,n11,n12,n15);
  init_prism(prism_conn+36, n7,n8,n10,n12,n13,n15);
  init_prism(prism_conn+42, n8,n9,n10,n13,n14,n15);

  init_tet(tet_conn,   n11,n12,n15,n16);
  init_tet(tet_conn+4, n12,n13,n15,n16);
  init_tet(tet_conn+8, n13,n14,n15,n16);
  init_tet(tet_conn+12,n14,n11,n15,n16);
}

void init_prism_tet_unit_coords(double *coords, double x, double y,
			 int n1, int n2, int n3, int n4,
			 int n5, int n6, int n7, int n8,
			 int n9, int n10, int n11, int n12,
			 int n13, int n14, int n15, int n16){

  coords[3*n1-3]  = 0.0+x; coords[3*n1-2]  = 0.0+y; coords[3*n1-1]  = 0.0; // 1
  coords[3*n2-3]  = 1.0+x; coords[3*n2-2]  = 0.0+y; coords[3*n2-1]  = 0.0; // 2
  coords[3*n3-3]  = 1.0+x; coords[3*n3-2]  = 1.0+y; coords[3*n3-1]  = 0.0; // 3
  coords[3*n4-3]  = 0.0+x; coords[3*n4-2]  = 1.0+y; coords[3*n4-1]  = 0.0; // 4
  coords[3*n5-3]  = 0.5+x; coords[3*n5-2]  = 0.5+y; coords[3*n5-1]  = 0.0; // 5 

  coords[3*n6-3]  = 0.0+x; coords[3*n6-2]  = 0.0+y; coords[3*n6-1]  = 1.0; // 6 
  coords[3*n7-3]  = 1.0+x; coords[3*n7-2]  = 0.0+y; coords[3*n7-1]  = 1.0; // 7 
  coords[3*n8-3]  = 1.0+x; coords[3*n8-2]  = 1.0+y; coords[3*n8-1]  = 1.0; // 8 
  coords[3*n9-3]  = 0.0+x; coords[3*n9-2]  = 1.0+y; coords[3*n9-1]  = 1.0; // 9 
  coords[3*n10-3] = .65+x; coords[3*n10-2] = .35+y; coords[3*n10-1] = 1.0; // 10
  
  coords[3*n11-3] = 0.0+x; coords[3*n11-2] = 0.0+y; coords[3*n11-1] = 2.0; // 11
  coords[3*n12-3] = 1.0+x; coords[3*n12-2] = 0.0+y; coords[3*n12-1] = 2.0; // 12
  coords[3*n13-3] = 1.0+x; coords[3*n13-2] = 1.0+y; coords[3*n13-1] = 2.0; // 13
  coords[3*n14-3] = 0.0+x; coords[3*n14-2] = 1.0+y; coords[3*n14-1] = 2.0; // 14
  coords[3*n15-3] = .27+x; coords[3*n15-2] = .29+y; coords[3*n15-1] = 2.24; // 15
  
  coords[3*n16-3] = 0.5+x; coords[3*n16-2] = 0.5+y; coords[3*n16-1] = 2.5; // 16
}

void build_unstr_prism_tet_2(){

  const int num_nodes   = 46;
  const int num_prisms1 = 32;
  const int num_prisms2 = 32;
  const int num_tets1   = 27;
  const int num_tets2   = 27;
  double coors1[num_nodes*3];
  double coors2[num_nodes*3];
  int    prisms1[num_prisms1*6];
  int    prisms2[num_prisms2*6];
  int    tets1[num_tets1*4];
  int    tets2[num_tets2*4];
  const int pconn_size = 15;
  int    pconn1[pconn_size] = {1,2,12,3,6,9,12,15,18,21,24,27,46,45,43};
  int    pconn2[pconn_size] = {1,1,12,1,4,7,10,13,16,19,22,25,40,44,46};
  

  // Create a new window named "unstr"
  COM_new_window("unstr_prism_tet_2");

  // Regist nodal coordinates, connectivity tables,
  // and the pconn for pane 1.
  COM_set_size( "unstr_prism_tet_2.nc", 1, num_nodes);
  COM_set_array( "unstr_prism_tet_2.nc", 1, &coors1[0],3);
  COM_set_size( "unstr_prism_tet_2.:W6:", 1, num_prisms1);
  COM_set_array( "unstr_prism_tet_2.:W6:", 1, &prisms1[0],6);
  COM_set_size( "unstr_prism_tet_2.:T4:", 1, num_tets1);
  COM_set_array( "unstr_prism_tet_2.:T4:", 1, &tets1[0],4);
  //  COM_set_size( "unstr_prism_tet_2.pconn", 1, pconn_size);
  //COM_set_array( "unstr_prism_tet_2.pconn", 1, &pconn1[0]);

  // Register nodal coordinates, connectivity tables,
  // and the pconn for pane 2.

  COM_set_size( "unstr_prism_tet_2.nc", 2, num_nodes);
  COM_set_array( "unstr_prism_tet_2.nc", 2, &coors2[0],3);
  COM_set_size( "unstr_prism_tet_2.:W6:", 2, num_prisms2);
  COM_set_array( "unstr_prism_tet_2.:W6:", 2, &prisms2[0],6);
  COM_set_size( "unstr_prism_tet_2.:T4:", 2, num_tets2);
  COM_set_array( "unstr_prism_tet_2.:T4:", 2, &tets2[0],4);
  //  COM_set_size( "unstr_prism_tet_2.pconn", 2, pconn_size);
  //COM_set_array( "unstr_prism_tet_2.pconn", 2, &pconn2[0]);

  COM_window_init_done("unstr_prism_tet_2");

  init_prism_tet_unit(prisms1, tets1,
		      1,2,5,4,28,10,11,14,13,32,
		      19,20,23,22,36,40);
  
  init_prism_tet_unit(prisms1+8*6,tets1+4*4,
		      2,3,6,5,29,11,12,15,14,33,
		      20,21,24,23,37,41);
  
  init_prism_tet_unit(prisms1+16*6,tets1+8*4,
		      4,5,8,7,30,13,14,17,16,34,
		      22,23,26,25,38,42);
  
  init_prism_tet_unit(prisms1+24*6,tets1+12*4,
		      5,6,9,8,31,14,15,18,17,35,
		      23,24,27,26,39,43);
  
  init_prism_tet_unit_coords(coors1, 0.0, 0.0,
			     1,2,5,4,28,10,11,14,13,32,
			     19,20,23,22,36,40);
  
  init_prism_tet_unit_coords(coors1, 1.0, 0.0,
			     2,3,6,5,29,11,12,15,14,33,
			     20,21,24,23,37,41);
  
  init_prism_tet_unit_coords(coors1,0.0,1.0,
			     4,5,8,7,30,13,14,17,16,34,
			     22,23,26,25,38,42);
  
  init_prism_tet_unit_coords(coors1,1.0,1.0,
			     5,6,9,8,31,14,15,18,17,35,
			     23,24,27,26,39,43);

  coors1[129] = 1.0; coors1[130] = 1.0; coors1[131] = 2.5;
  coors1[132] = 2.0; coors1[133] = 1.0; coors1[134] = 2.5;
  coors1[135] = 2.5; coors1[136] = 0.5; coors1[137] = 2.5;
  
  init_tet(tets1+64,44,41,40,23); init_tet(tets1+68,44,43,41,23);
  init_tet(tets1+72,44,42,43,23); init_tet(tets1+76,44,40,42,23);
  init_tet(tets1+80,40,41,20,23); init_tet(tets1+84,41,43,24,23);
  init_tet(tets1+88,43,42,26,23); init_tet(tets1+92,42,40,22,23);
  init_tet(tets1+96,41,46,21,24); init_tet(tets1+100,41,45,46,24);
  init_tet(tets1+104,41,43,45,24);

  init_prism_tet_unit(prisms2, tets2,
		      1,2,5,4,28,10,11,14,13,32,
		      19,20,23,22,36,40);
  
  init_prism_tet_unit(prisms2+8*6,tets2+4*4,
		      2,3,6,5,29,11,12,15,14,33,
		      20,21,24,23,37,41);
  
  init_prism_tet_unit(prisms2+16*6,tets2+8*4,
		      4,5,8,7,30,13,14,17,16,34,
		      22,23,26,25,38,42);
  
  init_prism_tet_unit(prisms2+24*6,tets2+12*4,
		      5,6,9,8,31,14,15,18,17,35,
		      23,24,27,26,39,43);
  
  init_prism_tet_unit_coords(coors2, 2.0, 0.0,
			     1,2,5,4,28,10,11,14,13,32,
			     19,20,23,22,36,40);
  
  init_prism_tet_unit_coords(coors2, 3.0, 0.0,
			     2,3,6,5,29,11,12,15,14,33,
			     20,21,24,23,37,41);
  
  init_prism_tet_unit_coords(coors2,2.0,1.0,
			     4,5,8,7,30,13,14,17,16,34,
			     22,23,26,25,38,42);
  
  init_prism_tet_unit_coords(coors2,3.0,1.0,
			     5,6,9,8,31,14,15,18,17,35,
			     23,24,27,26,39,43);

  coors2[129] = 2.0; coors2[130] = 1.0; coors2[131] = 2.5;
  coors2[132] = 3.0; coors2[133] = 1.0; coors2[134] = 2.5;  
  coors2[135] = 1.5; coors2[136] = 1.5; coors2[137] = 2.5;

  init_tet(tets2+64,45,41,40,23); init_tet(tets2+68,45,43,41,23);
  init_tet(tets2+72,45,42,43,23); init_tet(tets2+76,45,40,42,23);
  init_tet(tets2+80,40,41,20,23); init_tet(tets2+84,41,43,24,23);
  init_tet(tets2+88,43,42,26,23); init_tet(tets2+92,42,40,22,23);
  init_tet(tets2+96,40,44,42,22); init_tet(tets2+100,44,46,42,22);
  init_tet(tets2+104,46,25,42,22);

  // Write the mesh to file.
  int OUT_write = COM_get_function_handle( "OUT.write_dataitem");
  int HDL_all = COM_get_dataitem_handle("unstr_prism_tet_2.all");
  COM_call_function( OUT_write, "unstr_prism_tet_2", &HDL_all, 
		     "unstr_prism_tet_2", "001");
}






