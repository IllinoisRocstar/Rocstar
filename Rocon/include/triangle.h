#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <cmath>
#include "NVec.h"

using namespace nvc;

/* Ray-Triangle Intersection Test Routines          */
/* code from journals of graphics tools (JGT)       */
/* http://www.acm.org/jgt/                          */
/* by Tomas Moller, May 2000                        */

#define EPSILON 1e-24

#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 

int intersect_ray_triangle( const Vec3D &RP0, 
			    const Vec3D & dir,
			    Vec3D & TV0, 
			    Vec3D & TV1,
			    Vec3D & TV2, 
			    Vec3D & I );

/* the original jgt code */
int intersect_triangle(double orig[3], double dir[3],
		       double vert0[3], double vert1[3], double vert2[3],
		       double *t, double *u, double *v);

/* code rewritten to do tests on the sign of the determinant */
/* the division is at the end in the code                    */
int intersect_triangle1(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v);


/* code rewritten to do tests on the sign of the determinant */
/* the division is before the test of the sign of the det    */
int intersect_triangle2(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v);


/* code rewritten to do tests on the sign of the determinant */
/* the division is before the test of the sign of the det    */
/* and one CROSS has been moved out from the if-else if-else */
int intersect_triangle3(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v);

// ray-triangle intersection -- prateek
//int P_triRay(const RAYTRI *rt, double &t);

/* Wrapper for all other methods ? */
bool rayIntersectsTriangle(double *p, double *d, 
			   double *v0, double *v1, double *v2, Vec3D &I);

#endif
