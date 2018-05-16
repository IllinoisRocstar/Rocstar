/* Ray-Triangle Intersection Test Routines          */
/* Different optimizations of my and Ben Trumbore's */
/* code from journals of graphics tools (JGT)       */
/* http://www.acm.org/jgt/                          */
/* by Tomas Moller, May 2000                        */

using namespace std;

#include <cmath>

#include "NVec.h"
using namespace nvc;

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
#define PLUCKER(pline, pt1, pt2) \
		  pline[0] = pt1[0]*pt2[1] - pt2[0]*pt1[1]; \
		  pline[1] = pt1[0]*pt2[2] - pt2[0]*pt1[2]; \
		  pline[2] = pt1[0]        - pt2[0];        \
		  pline[3] = pt1[1]*pt2[2] - pt2[1]*pt1[2]; \
		  pline[4] = pt1[2]        - pt2[2];        \
		  pline[5] = pt2[1]        - pt1[1];
#define SIDE(p,q) (p[0]*q[4]+p[1]*q[5]+p[2]*q[3]+p[3]*q[2]+p[4]*q[0]+p[5]*q[1])
#define	RAYPOINT(result,start,end,dist) { result[0]=start[0]+d*(end[0]-start[0]); result[1]=start[1]+d*(end[1]-start[1]); result[2]=start[2]+d*(end[2]-start[2]); }

/* Sample implementation of the line segment—triangle intersection test
 * presented in the paper:
 *
 *   Fast 3D Line Segment—Triangle Intersection Test
 *   Nick Chirkov
 *   journal of graphics tools 10(3):13-18, 2005
 *
 */

struct RAYTRI
{
	double org[3];
	double end[3];
	double dir[3];
	double v0[3],v1[3],v2[3];

	struct PLANE
	{
		double x, y, z, d;
		enum MAIN_AXIS { X, Y, Z };
		MAIN_AXIS type;
	};
	PLANE plane;
};

//Plucker Co-ordinates
int P_triRay(double *p, double *d, double *v0, 
			   double *v1, double *v2, Vec3D & I)
{
	int flag = 0;
	double *q,*L,*e1,*e2,*e3, s1, s2,s3;
	q = new double[3];
	for(int i=0;i<3;++i)q[i] = p[i]+d[i];
	
	L = new double[6];
	e1 = new double[6];
	e2 = new double[6];
	e3 = new double[6];
	
	PLUCKER(L,p,q);
	PLUCKER(e1,v0,v1);
	PLUCKER(e2,v1,v2);
	PLUCKER(e3,v2,v0);
	
	s1 = SIDE(L,e1);
	s2 = SIDE(L,e2);
	s3 = SIDE(L,e3);
	
	if(s1==0.0 && s2==0.0 && s3==0.0)
	{   // co-planar
		int pts[4][3] = {{0,0,0},{0,1,0},{0,0,1},{1,0,0}};
		int edge[3] = {0,0,0};
		int vertex[3] = {0,0,0};
		int *x,i,intersect=0;
		double cs1,cs2,cs3;
		for(i=0;i<4;++i)
		{
			PLUCKER(e1,v0,v1);
			PLUCKER(e2,v1,pts[i]);
			PLUCKER(e3,pts[i],v0);
			cs1 = SIDE(L,e1);
			cs2 = SIDE(L,e2);
			cs3 = SIDE(L,e3);		
			if(!(cs1==0 && cs2==0 && cs3==0)) break;
		}
		
		x = pts[i];
		
		for(i=0;i<3;++i)
		{
			if(i==0){
				PLUCKER(e2,v1,x);
				PLUCKER(e3,x,v0);
				cs2 = SIDE(L,e2);
				cs3 = SIDE(L,e3);		
			} else if(i==1){
				PLUCKER(e2,v2,x);
				PLUCKER(e3,x,v1);
				cs2 = SIDE(L,e2);
				cs3 = SIDE(L,e3);						
			} else {
				PLUCKER(e2,v0,x);
				PLUCKER(e3,x,v2);
				cs2 = SIDE(L,e2);
				cs3 = SIDE(L,e3);		
			}
			
			if(cs2*cs3<0);
			else if(cs2*cs3>0)
				{ intersect++;
				  edge[i]++;
				}
			else if(cs2==0 && cs3==0)
				{ intersect+=2;
				  edge[i]+=2;
				  break;}
			else if (cs2==0)
				{intersect++;
				 vertex[(i+1)-3*((i+1)/3)]++;
				}
			else if (cs3==0)
				{intersect++;
				 vertex[i]++;
				}
		}
		 // implement line-line intersect 
		 // implement vertex-line intersect 
		if(intersect!=0)
		{
			if(intersect<=3){
				if(vertex[0]!=0 && !flag)
				{
					I[0] = v0[0]; I[1]=v0[1]; I[2]=v0[2];
					flag = 1;
				}
				if(vertex[1]!=0 && !flag)
				{
					I[0] = v1[0]; I[1]=v1[1]; I[2]=v1[2];
					flag = 1;
				}
				if(vertex[2]!=0 && !flag)
				{
					I[0] = v2[0]; I[1]=v2[1]; I[2]=v2[2];
					flag = 1;
				}
			}
			if(!flag && (edge[0]==1||edge[0]==2) )
			{
				I[0]=(v0[0]+v1[0])/2;I[1]=(v0[1]+v1[1])/2;I[2]=(v0[2]+v1[2])/2;
				flag = 1;
			}
			if(!flag && (edge[1]==1||edge[1]==2) )
			{
				I[0]=(v1[0]+v2[0])/2;I[1]=(v1[1]+v2[1])/2;I[2]=(v1[2]+v2[2])/2;
				flag = 1;
			}
			if(!flag && (edge[2]==1||edge[2]==2) )
			{
				I[0]=(v2[0]+v0[0])/2;I[1]=(v0[1]+v2[1])/2;I[2]=(v0[2]+v2[2])/2;
				flag = 1;
			}
		}
	} else if ((s1>0 && s2>0 && s3>0) || (s1<0 && s2<0 && s3<0)){
		double al,be,ga;
		al = s1/(s1+s2+s3);
		be = s2/(s1+s2+s3);
		ga = s3/(s1+s2+s3);
		for(int i=0;i<3;++i)
			I[i] = al*v2[i] + be*v0[i] + ga*v1[i];
		flag = 1;
	} else if ( s1==0 && s2*s3>0){
		I[0]=(v0[0]+v1[0])/2;I[1]=(v0[1]+v1[1])/2;I[2]=(v0[2]+v1[2])/2;
		flag = 1;
	} else if ( s2==0 && s1*s3>0){
		I[0]=(v1[0]+v2[0])/2;I[1]=(v1[1]+v2[1])/2;I[2]=(v1[2]+v2[2])/2;
		flag = 1;
	} else if ( s3==0 && s1*s2>0){
		I[0]=(v2[0]+v0[0])/2;I[1]=(v0[1]+v2[1])/2;I[2]=(v0[2]+v2[2])/2;
		flag = 1;
	} else if ( s1==0 && s2==0){
		I[0]=v1[0];I[1]=v1[1];I[2]=v1[2];
		flag = 1;
	} else if ( s1==0 && s3==0){
		I[0]=v0[0];I[1]=v0[1];I[2]=v0[2];
		flag = 1;
	} else if ( s2==0 && s3==0){
		I[0]=v2[0];I[1]=v2[1];I[2]=v2[2];
		flag = 1;
	} else
		flag = 0;
	
	delete [] q;
	delete [] L;
	delete [] e1;
	delete [] e2;
	delete [] e3;
	
	return flag;
}

int c2005(const RAYTRI* rt)
{
	double e0[3],e1[3],e2[3],norm[3],point[3],v[3],av[3],vb[3],vc[3];
	SUB(e0, rt->v1, rt->v0);
	SUB(e1, rt->v2, rt->v0);
	CROSS(norm,e0,e1);

	double pd = DOT(norm, rt->v0);

	double signSrc = DOT(norm, rt->org) -pd;
	double signDst = DOT(norm, rt->end) -pd;

	if(signSrc*signDst > 0.0) return 0;

	double d = signSrc/(signSrc - signDst);

	RAYPOINT(point, rt->org, rt->end,d);
	SUB(v, point, rt->v0);		
	CROSS(av,e0,v);
	CROSS(vb,v,e1);

	if(DOT(av,vb) >= 0.0)
	{
		SUB(e2, rt->v1, rt->v2);
		SUB(v, point, rt->v1);
		CROSS(vc,v,e2);
		if(DOT(av,vc) >= 0.0) return 1;
	}
	return 0;
}

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// Assume that classes are already given for the objects:
//    Point and Vector with
//        coordinates {double x, y, z;}
//        operators for:
//            == to test equality
//            != to test inequality
//            (Vector)0 = (0,0,0)         (null vector)
//            Point  = Point ± Vector
//            Vector = Point - Point
//            Vector = Scalar * Vector    (scalar product)
//            Vector = Vector * Vector    (cross product)
//    Line and Ray and Segment with defining points {Point P0, P1;}
//        (a Line is infinite, Rays and Segments start at P0)
//        (a Ray extends beyond P1, but a Segment ends at P1)
//    Plane with a point and a normal {Point V0; Vector n;}
//    Triangle with defining vertices {Point V0, V1, V2;}
//    Polyline and Polygon with n vertices {int n; Point *V;}
//        (a Polygon has V[n]=V[0])
//===================================================================

#define SMALL_NUM  1e-24 // anything that avoids division overflow
// dot product (3D) which allows vector operations in arguments

// intersect_RayTriangle(): intersect a ray with a 3D triangle
//    Input:  a ray R, and a triangle T
//    Output: *I = intersection point (when it exists)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 = disjoint (no intersect)
//             1 = intersect in unique point I1
//             2 = are in the same plane
int intersect_ray_triangle( const Vec3D &RP0, const Vec3D & dir,Vec3D & TV0, 
			    Vec3D & TV1,Vec3D & TV2, Vec3D & I )
{
    Vec3D    u, v, n;             // triangle vectors
    Vec3D     w0, w;          // ray vectors
    double    r, a, b;             // params to calc ray-plane intersect

    // get triangle edge vectors and plane normal
    u = TV1 - TV0;
    v = TV2 - TV0;
    n = cross(u,v);             // cross product
    if (n == Vec3D(0,0,0))            // triangle is degenerate
        return 0;                 // do not deal with this case

    //dir = RP1 - RP0;             // ray direction vector
    w0 = RP0 - TV0;
    a = -1.0*(n*w0);
    b = n*dir;
    if (fabs(b) < SMALL_NUM) {     // ray is parallel to triangle plane
      if (a != 0)                // ray lies in triangle plane
		//return 1;
      //else 
		return 0;             // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                   // ray goes away from triangle
      return 0;                  // => no intersect
    if (r > 1.0)                   // ray goes away from triangle
      return 0;					//change by prateek, 0 instead of 3
    // for a segment, also test if (r > 1.0) => no intersect

    I = RP0 + (r * dir);           // intersect point of ray and plane

    // is I inside T?
    double    uu, uv, vv, wu, wv, D;
    uu = (u*u);
    uv = (u*v);
    vv = (v*v);
    w =  I - TV0;
    wu = w*u;
    wv = w*v;
    D = (uv * uv) -( uu * vv);

    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // I is outside T
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
        return 0;

    return 1;                      // I is in T
}

/* 
	*the original jgt code 
	*modified by prateek -- checking if 0<=t<=1 for line segment
*/
int intersect_triangle(double orig[3], double dir[3],
		       double vert0[3], double vert1[3], double vert2[3],
		       double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
   double det,inv_det;

   /* find vectors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pvec, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pvec);

   if (det > -EPSILON && det < EPSILON)
     return 0;
   inv_det = 1.0 / det;

   /* calculate distance from vert0 to ray origin */
   SUB(tvec, orig, vert0);

   /* calculate U parameter and test bounds */
   *u = DOT(tvec, pvec) * inv_det;
   if (*u < 0.0 || *u > 1.0)
     return 0;

   /* prepare to test V parameter */
   CROSS(qvec, tvec, edge1);

   /* calculate V parameter and test bounds */
   *v = DOT(dir, qvec) * inv_det;
   if (*v < 0.0 || *u + *v > 1.0)
     return 0;

   /* calculate t, ray intersects triangle */
   *t = DOT(edge2, qvec) * inv_det;
   
   if(*t<0 || *t>1)
     return 0; //enforce that ray is line segment

   return 1;
}


/* code rewritten to do tests on the sign of the determinant */
/* the division is at the end in the code                    */
int intersect_triangle1(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
   double det,inv_det;

   /* find vectors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pvec, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pvec);

   if (det > EPSILON)
   {
      /* calculate distance from vert0 to ray origin */
      SUB(tvec, orig, vert0);
      
      /* calculate U parameter and test bounds */
      *u = DOT(tvec, pvec);
      if (*u < 0.0 || *u > det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qvec, tvec, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qvec);
      if (*v < 0.0 || *u + *v > det)
	 return 0;
      
   }
   else if(det < -EPSILON)
   {
      /* calculate distance from vert0 to ray origin */
      SUB(tvec, orig, vert0);
      
      /* calculate U parameter and test bounds */
      *u = DOT(tvec, pvec);
/*      printf("*u=%f\n",(double)*u); */
/*      printf("det=%f\n",det); */
      if (*u > 0.0 || *u < det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qvec, tvec, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qvec) ;
      if (*v > 0.0 || *u + *v < det)
	 return 0;
   }
   else return 0;  /* ray is parallell to the plane of the triangle */


   inv_det = 1.0 / det;

   /* calculate t, ray intersects triangle */
   *t = DOT(edge2, qvec) * inv_det;
   (*u) *= inv_det;
   (*v) *= inv_det;

   return 1;
}

/* code rewritten to do tests on the sign of the determinant */
/* the division is before the test of the sign of the det    */
int intersect_triangle2(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
   double det,inv_det;

   /* find vectors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pvec, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pvec);

   /* calculate distance from vert0 to ray origin */
   SUB(tvec, orig, vert0);
   inv_det = 1.0 / det;
   
   if (det > EPSILON)
   {
      /* calculate U parameter and test bounds */
      *u = DOT(tvec, pvec);
      if (*u < 0.0 || *u > det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qvec, tvec, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qvec);
      if (*v < 0.0 || *u + *v > det)
	 return 0;
      
   }
   else if(det < -EPSILON)
   {
      /* calculate U parameter and test bounds */
      *u = DOT(tvec, pvec);
      if (*u > 0.0 || *u < det)
	 return 0;
      
      /* prepare to test V parameter */
      CROSS(qvec, tvec, edge1);
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qvec) ;
      if (*v > 0.0 || *u + *v < det)
	 return 0;
   }
   else return 0;  /* ray is parallell to the plane of the triangle */

   /* calculate t, ray intersects triangle */
   *t = DOT(edge2, qvec) * inv_det;
   (*u) *= inv_det;
   (*v) *= inv_det;

   return 1;
}

/* code rewritten to do tests on the sign of the determinant */
/* the division is before the test of the sign of the det    */
/* and one CROSS has been moved out from the if-else if-else */
int intersect_triangle3(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v)
{
   double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
   double det,inv_det;

   /* find vectors for two edges sharing vert0 */
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   /* begin calculating determinant - also used to calculate U parameter */
   CROSS(pvec, dir, edge2);

   /* if determinant is near zero, ray lies in plane of triangle */
   det = DOT(edge1, pvec);

   /* calculate distance from vert0 to ray origin */
   SUB(tvec, orig, vert0);
   inv_det = 1.0 / det;
   
   CROSS(qvec, tvec, edge1);
      
   if (det > EPSILON)
   {
      *u = DOT(tvec, pvec);
      if (*u < 0.0 || *u > det)
	 return 0;
            
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qvec);
      if (*v < 0.0 || *u + *v > det)
	 return 0;
      
   }
   else if(det < -EPSILON)
   {
      /* calculate U parameter and test bounds */
      *u = DOT(tvec, pvec);
      if (*u > 0.0 || *u < det)
	 return 0;
      
      /* calculate V parameter and test bounds */
      *v = DOT(dir, qvec) ;
      if (*v > 0.0 || *u + *v < det)
	 return 0;
   }
   else return 0;  /* ray is parallell to the plane of the triangle */

   *t = DOT(edge2, qvec) * inv_det;
   (*u) *= inv_det;
   (*v) *= inv_det;

   return 1;
}

/* Wrapper for all other methods ? */
bool rayIntersectsTriangle(double *p, double *d, 
			   double *v0, double *v1, double *v2, Vec3D & I)
{
	/*Vec3D org(p[0],p[1],p[2]);
	Vec3D dir(d[0],d[1],d[2]);
	Vec3D ver0(v0[0],v0[1],v0[2]);
	Vec3D ver1(v1[0],v1[1],v1[2]);
	Vec3D ver2(v2[0],v2[1],v2[2]);
	int result = intersect_ray_triangle(org,dir,ver0,ver1,ver2,I);*/
	
	/*RAYTRI * input = new RAYTRI;
	for(int i =0;i<3;++i)
	{
		input->org[i] = p[i];
		input->dir[i] = d[i];
		input->v0[i] = v0[i];
		input->v1[i] = v1[i];
		input->v2[i] = v2[i];
		input->end[i] = p[i]+d[i];
	}
	int result = c2005(input);*/
	
	/*double *t,*u,*v;
	t = new double;
	u = new double;
	v = new double;
	*t = *u = *v = 0;
	int result = intersect_triangle(p,d,v0,v1,v2,t,u,v);
	for(int j=0;j<3;j++)
	{
	  I[j]=p[j]+ (*t)*d[j];//I[j]=v0[j]+ (*u)*v1[j] + (*v)*v2[j];
	}*/
	
	int result = P_triRay(p,d,v0,v1,v2,I);
	
	if(result)
	{	double *t,*u,*v;
		t = new double;
		u = new double;
		v = new double;
		*t = *u = *v = 0;
		result = intersect_triangle(p,d,v0,v1,v2,t,u,v);
		for(int j=0;j<3;j++)
		{
		  I[j]=p[j]+ (*t)*d[j];
		}
		delete t;
		delete u;
		delete v;
	}
	return result;
}
