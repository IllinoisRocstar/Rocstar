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
/*! \file geometry.C 
    \brief Geometric functions.

    Implementation of geometric helper functions used by 2D and 3D geometric mesh
    quality measures.
*/

#include "geometry.h"
#include <cassert>

MOP_BEGIN_NAMESPACE
using namespace std;

void min_max3(double a1, double a2, double a3, double &min, double &max){
min = ( a1 < a2 )? a1 : a2;    max = ( a1 > a2 )? a1 : a2;
min = ( min < a3 ) ? min : a3; max = ( max > a3 ) ? max : a3;
}

void min_max4(double a1, double a2, double a3, double a4, double &min, double &max){
  min_max3(a1,a2,a3,min,max);
  min = (min < a4)? min : a4; max = (max > a4)? max : a4;
}

void min_max6(double a1, double a2, double a3, double a4, 
	      double a5, double a6, double &min, double &max){
  double min2, max2;
  min_max3(a1,a2,a3,min,max);     min_max3(a4,a5,a6,min2,max2);
  min = (min < min2)? min : min2; max = (max > max2)? max : max2;
}

void min_max12(double a1, double a2, double a3, double a4, 
	      double a5, double a6, double a7, double a8,
	      double a9, double a10, double a11, double a12, double &min, double &max){
  double min2, max2;
  min_max6(a1,a2,a3,a4,a5,a6,min,max); min_max6(a7,a8,a9,a10,a11,a12,min2,max2);
  min = (min < min2)? min : min2;      max = (max > max2)? max : max2;
}

double angle(Vector_3<double> v1, 
	     Vector_3<double> v2){
  double dot_product, norm1, norm2, angle;
  dot_product = dotP(v1,v2);
  norm1 = sqrt(dotP(v1,v1)); norm2 = sqrt(dotP(v2,v2));
  angle = acos(dot_product/(norm1*norm2));
  assert(dot_product/(norm1*norm2) <=1 && dot_product/(norm1*norm2)>= -1);
  angle *= (double)180; angle  = angle/M_PI;
  return angle;
}

double f_angle(Vector_3<double> v1, Vector_3<double>v2){
  return double(180) - angle(v1,v2);
}

void printv(const Vector_3<double> & v){
  cout << v[0] <<"  " << v[1] << "  " << v[2] << endl;
}

void unit_normal(const Vector_3<double> & v1, 
		 const Vector_3<double> & v2,
		 Vector_3<double> &v3){
  crossP(v1,v2,v3);
  double normalizer = sqrt( dotP (v3,v3) );
  if (normalizer == 0) return;
  v3[0] = v3[0]/normalizer; v3[1] = v3[1]/normalizer; v3[2] = v3[2]/normalizer;
} 

double edge_length(const Vector_3<double> & v){
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double tri_area(const Vector_3<double> & v1,
		const Vector_3<double> & v2,
		const Vector_3<double> & v3){
  double a = edge_length(v1);
  double b = edge_length(v2);
  double c = edge_length(v3);
  double s = .5 * (a+b+c);
  return sqrt( s * (s-a) * (s-b) * (s-c) ); 
}

double tet_vol(const Vector_3<double> & a, 
	       const Vector_3<double> & b,
	       const Vector_3<double> & c){
  Vector_3<double> temp;
  crossP(b,c,temp);
  return  (1/(double)6) * abs (dotP(a,temp)  ) ;
}

MOP_END_NAMESPACE






