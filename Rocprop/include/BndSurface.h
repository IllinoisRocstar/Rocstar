#ifndef BNDSURFACE_INCLUDED 
#define BNDSURFACE_INCLUDED

#include <iostream>
#include <vector>
#include <algorithm>
#include <ext/hash_map>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <queue>
#include "triangle.h"
#include "NVec.h"
#include "KNN_Grid.h"

//****************************************************************************************************
// A C++ class for a bounding surface described as a triangulated mesh
//---------------------------------
//
// An application should:
//                  a) create a new MeshBndSurf object
//                  b) call initialize()
//                  c) call intersect() as many times as needed
//                  d) delete the MeshBndSurf object
//
// External interface:
//  
//------------------------------------------------------------------------------------------------------
//  void initialize(double * nodes, int * triangles, unsigned int nv, unsigned int nt, int divisions = 50)
//
//  Purpose: Creates a data structure for accelerates intersection tests against a  bounding surface 
//           A MeshBndSurf object requires the bounding surface be expressed as a triangulated surface mesh
//           in 3D space. It uses a simple cartesian grid to accelerate tests. It has no special persistance
//           so if it is a local variable it will be destroyed when it goes out of scope. Most uses
//           will probably allocate a MeshBndSurf using new to avoid this. 
//
//  NOTE: the data in the "nodes" and "triangles" arrays will be copied. These arrays should be deleted
//        after initialize() returns to conserve storage space
//
//
//  Input:
//         double * nodes: the coordinates of 3D nodes of the bounding mesh in an array of size 3*nv
//                         node i has it's x coordinate at nodes[3*i], y coordinate at nodes[3*i+1]
//                                         z coordinate at nodes[3*i+2]
//
//         unsigned int nv: number of nodes in the input mesh
//
//         int * triangles: the triangles of the bounding mesh in an array of size 3*nt
//                          triangle i has the indices of it's nodes stored at triangles[3*i], 
//                                 triangles[3*i+1], and triangles[3*i+2]
//
//         unsigned int nt: number of triangles in the input mesh
//         
//         int divisions: the number of grid divisions along the longest axis
//                        MUST BE BETWEEN 1 and 1023
//                        If perfomance is slow due too a densely tesselated  mesh, increase divisions
//       
//------------------------------------------------------------------------------------------------------ 
//  bool intersection(double * p, double * d, double * x)
//
//  Purpose: determines if displacing a point p by offset d will cross the bounding surface
//           if so, the point of intersection is the output
//
//  Input:
//          double * p: a point in 3D space exressed as an array of 3 doubles
//                      x coordinate is in p[0], y in p[1] and z in p[2] 
//
//          double * d: a displacment of p in 3D space expressed as an array of doubles
//
//  Ouput: double * x: the calling code provides storage for a 3D point passed as this parameter
//                     if the function finds an intersection, x will hold the 3D coordinates of the point of
//                     intersection between the bounding mesh and the line segment (p,p+d)  
//
//  Returns: true if line segment (p,p+d) intersects the bounding mesh. Otherwise false
//------------------------------------------------------------------------------------------------------
//
//

#define GRID_NO 100
#define MAX_GRID_DIV 1023

using namespace std;
using namespace nvc;

class MeshBndSurf
{
  
 private:
  KNN_Grid< Tri > * gSpace;
  Point3D bbmin, bbmax;
  Vec3D * node_array;
  Tri * tri_array;
  unsigned int numv;
  unsigned int numt;

  void compute_bbox(Point3D &bbmin, Point3D &bbmax)
  {
	if (numv == 0) 
	  { 
	    bbmin = bbmax = 0.0;
	  }
	else
	  {
	    bbmin = bbmax = node_array[0];
	  }
 
	for(unsigned int i=1; i<numv; i++)
	  {
	    Vertex3D & v = node_array[i];
	    
	    if(v[0] < bbmin[0])  bbmin[0]=v[0];
	    if(v[1] < bbmin[1])  bbmin[1]=v[1];
	    if(v[2] < bbmin[2])  bbmin[2]=v[2];
	    
	    if(v[0] > bbmax[0])  bbmax[0]=v[0];
	    if(v[1] > bbmax[1])  bbmax[1]=v[1];
	    if(v[2] > bbmax[2])  bbmax[2]=v[2];
	  }
  }
	
 public:

 MeshBndSurf():gSpace(NULL),node_array(NULL),tri_array(NULL),numv(0),numt(0){}

  MeshBndSurf(Vec3D * n, Tri * t,unsigned int nv,
	      unsigned int nt):node_array(n),tri_array(t),numv(nv),numt(nt)
  {
	compute_bbox(bbmin, bbmax);

	gSpace = new KNN_Grid< Tri > (bbmin, bbmax, 8);
	for(int i=0; i<nt; ++i){
	   Point3D p(n[t[i][0]][0], n[t[i][0]][1], n[t[i][0]][2]);
		Point3D q(n[t[i][1]][0], n[t[i][1]][1], n[t[i][1]][2]);
		Point3D r(n[t[i][2]][0], n[t[i][2]][1], n[t[i][2]][2]);
	   gSpace->insertRange(t[i], p,q,r);
	}
  }

  ~MeshBndSurf(){
	delete gSpace;
	delete [] node_array;
	delete [] tri_array;
  }
  
bool intersection(const Vec3D & p, const Vec3D & d, Vec3D & x)
  {
    double orig[3];
    double dir[3];
    double v1[3];
    double v2[3];
    double v3[3];
    double t;
	
	Point3D a,b;
	a = p;
	b = p+d; 
	vector<Tri> * hitTri = gSpace->get_cell_range(a,b);
	
	t = 0;
	
	for(int j=0;j<3;j++)
	  {
	    orig[j] = p[j];
	    dir[j]  = d[j]; 
	  }
	
	for(vector<Tri>::iterator it = hitTri->begin(); it !=hitTri->end(); it++)
	  {
	    Tri & cur = *it;
	    for(int j=0;j<3;j++)
	      {
		v1[j]=(node_array[cur[0]])[j];
		v2[j]=(node_array[cur[1]])[j];
		v3[j]=(node_array[cur[2]])[j];
	      }
	    if (rayIntersectsTriangle(orig, dir,v1,v2,v3,x))
	      {
		//Debugging info
		//cerr << "XXXX Ray " << p << " to " << d << " Tri " << endl;
		//cerr << node_array[cur[0]] << endl;
		//cerr << node_array[cur[1]] << endl;
		//cerr << node_array[cur[2]] << endl << endl;
		return true;
	      }
	    //else 
	    //{
	    //cerr<<"0000 Ray "<< p << " to" << d << " Tri " << endl;
	    //}
	  }
	return false;
  }      

 void initialize(double * nodes, int * triangles, unsigned int nv, unsigned int nt, int divisions = 50)
 {
   numv=nv;
    numt=nt;
    //make sure grid divisions are sane
    if (divisions > MAX_GRID_DIV)
      {
	cerr << "Grid divsisions on longest dimension must be less than or equal to " << MAX_GRID_DIV << endl;
	divisions = MAX_GRID_DIV;
      }

    //Copy nodal coordinates
    if (node_array != NULL)
      delete [] node_array;

    node_array = new Vec3D[numv];
    for(int i=0;i<numv;i++)
      {
	for(int j=0;j<3;j++)
	  {
	    (node_array[i])[j]=nodes[i*3+j]; 
	  }
      }

    //Copy triangle connectivity
    if (tri_array != NULL)
      delete [] tri_array;

    tri_array = new Tri[numt];
    for(int i=0;i<numt;i++)
      {
	for(int j=0;j<3;j++)
	  {
	    (tri_array[i])[j]=triangles[i*3+j]; 
	  }
      }

    compute_bbox(bbmin, bbmax);
    if (gSpace != NULL)
      delete gSpace;

    gSpace = new KNN_Grid< Tri > (bbmin, bbmax, divisions);
    for(int i=0; i<nt; ++i)
      {
	Point3D p(node_array[tri_array[i][0]][0], node_array[tri_array[i][0]][1], node_array[tri_array[i][0]][2]);
	Point3D q(node_array[tri_array[i][1]][0], node_array[tri_array[i][1]][1], node_array[tri_array[i][1]][2]);
	Point3D r(node_array[tri_array[i][2]][0], node_array[tri_array[i][2]][1], node_array[tri_array[i][2]][2]);
	gSpace->insertRange(tri_array[i], p,q,r);
      }   
  }

  bool intersection(double * p, double * d, double * x)
  {
    double orig[3];
    double dir[3];
    double v1[3];
    double v2[3];
    double v3[3];
    double t;
    Vec3D xint;	
    Point3D a,b;
    for(int j=0;j<3;j++)
      {
	a[j] = p[j];
	b[j] = p[j]+d[j]; 
      } 
    vector<Tri> * hitTri = gSpace->get_cell_range(a,b);
    
    t = 0;
	
    for(int j=0;j<3;j++)
      {
	orig[j] = p[j];
	dir[j]  = d[j]; 
      }
	
	for(vector<Tri>::iterator it = hitTri->begin(); it !=hitTri->end(); it++)
	  {
	    Tri & cur = *it;
	    for(int j=0;j<3;j++)
	      {
		v1[j]=(node_array[cur[0]])[j];
		v2[j]=(node_array[cur[1]])[j];
		v3[j]=(node_array[cur[2]])[j];
	      }
	    if (rayIntersectsTriangle(orig, dir,v1,v2,v3,xint))
	      {
		for(int i=0;i<3;i++)
		  x[i]=xint[i];
		return true;
	      }
	  }
	return false;
  }  

    
};

  
#endif
