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
// $Id: geometry.h,v 1.3 2008/12/06 08:45:24 mtcampbe Exp $

/*! \file geometry.h 
    \brief Geometric helper function header file.
*/

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "mopbasic.h"
#include <cmath>
#include <vector>

MOP_BEGIN_NAMESPACE

//! Find minimum  and maximum of 3 numbers.
void min_max3(double a1, double a2, double a3, double &min, 
	      double &max);

//! Find minimum  and maximum of 4 numbers.
void min_max4(double a1, double a2, double a3, double a4, 
	      double &min, double &max);

//! Find minimum  and maximum of 6 numbers.
void min_max6(double a1, double a2, double a3, double a4, 
	      double a5, double a6, double &min, double &max);

//! Find minimum  and maximum of 12 numbers.
void min_max12(double a1, double a2, double a3, double a4, 
	       double a5, double a6, double a7, double a8,
	       double a9, double a10, double a11, double a12, 
	       double &min, double &max);

//! Compute the angle between two vectors.
double angle(Vector_3<double> v1, 
	     Vector_3<double> v2);

//! Compute the angle between two faces.
double f_angle(Vector_3<double> v1, Vector_3<double>v2);

//! Compute the dot Product of two vectors.
inline double dotP(const Vector_3<double> & v1,
		   const Vector_3<double> & v2){
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

//! Compute the cross Product of two vectors.
inline void crossP(const Vector_3<double> & v1, 
	    const Vector_3<double> & v2, 
	    Vector_3<double> & v3){
  v3[0]= (v1[1]*v2[2]-v2[1]*v1[2]); v3[1]= (v2[0]*v1[2]-v1[0]*v2[2]);
  v3[2]= (v1[0]*v2[1]-v2[0]*v1[1]);  
}

//! Compute the unit vector normal to two input vectors.
void unit_normal(const Vector_3<double> & v1, 
		 const Vector_3<double> & v2,
		 Vector_3<double> & v3);

//! Compute the edge length of a vector.
double edge_length(const Vector_3<double> & v);

//! Compute the area of a triangle defined by its three sides as vectors
double tri_area(const Vector_3<double> & v1,
		const Vector_3<double> & v2,
		const Vector_3<double> & v3);


//! Compute the volume of a tetrahedron given three edges from a point as vectors.
double tet_vol(const Vector_3<double> & a, 
	       const Vector_3<double> & b, 
	       const Vector_3<double> & c);

//! Print a vector.
void printv(const Vector_3<double> & v);

MOP_END_NAMESPACE

#endif






