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
// $Id: Geometric_Metrics_2.C,v 1.5 2008/12/06 08:45:24 mtcampbe Exp $

/*! \file Geometric_Metrics_2.C 
    \brief Implementation of 2D geometric metric quality measures.

    Implementation of geometric metrics:
    minimum and maximum angles,
    circumradius over inradius (R/r), and
    circumradius over shortest edge length (R/l)
    The latter two measures are defined only for simplicial elements.
*/

#include "geometry.h"
#include "Geometric_Metrics_2.h"
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

MOP_BEGIN_NAMESPACE

using namespace std;

void Geo_Metric_Base_2::initialize(Element_node_enumerator &ene){
  type_ = ene.type();
  Element_node_vectors_k_const<double> n;
  n.set(ene.pane()->dataitem("nc"),ene);
  if(type_ == COM::Connectivity::TRI3){
    v[0][0]=n(2,0)-n(1,0); v[0][1]=n(2,1)-n(1,1); v[0][2]=n(2,2)-n(1,2);
    v[1][0]=n(1,0)-n(2,0); v[1][1]=n(1,1)-n(2,1); v[1][2]=n(1,2)-n(2,2);    
    v[2][0]=n(2,0)-n(0,0); v[2][1]=n(2,1)-n(0,1); v[2][2]=n(2,2)-n(0,2);
    v[3][0]=n(1,0)-n(0,0); v[3][1]=n(1,1)-n(0,1); v[3][2]=n(1,2)-n(0,2);
  }
  else if (type_ == COM::Connectivity::QUAD4){
    v[0][0]=n(1,0)-n(0,0); v[0][1]=n(1,1)-n(0,1); v[0][2]=n(1,2)-n(0,2);
    v[1][0]=n(3,0)-n(0,0); v[1][1]=n(3,1)-n(0,1); v[1][2]=n(3,2)-n(0,2);    
    v[2][0]=n(1,0)-n(2,0); v[2][1]=n(1,1)-n(2,1); v[2][2]=n(1,2)-n(2,2);
    v[3][0]=n(3,0)-n(2,0); v[3][1]=n(3,1)-n(2,1); v[3][2]=n(3,2)-n(2,2);
  }
  else COM_assertion_msg(0,"Element type not supported for 2D geometric metrics");
}

void Geo_Metric_Base_2::initialize(Vector_3<double> n[],
				   int type){
  type_ = type;
  if (type_ == COM::Connectivity::TRI3){
    v[0]= n[2]-n[1]; v[1]= n[1]-n[2]; v[2]= n[2]-n[0]; v[3]= n[1]-n[0];
  }
  else if (type_ == COM::Connectivity::QUAD4){
    v[0]= n[1]-n[0]; v[1]= n[3]-n[0]; v[2]= n[1]-n[2]; v[3]= n[3] - n[2];
  }
}

void Geo_Metric_Base_2::compute_angles(double& min, double& max) const{
  if (type_ == COM::Connectivity::TRI3){
    double a1,a2,a3; 
    a1= angle(v[0],v[2]); a2= angle(v[1],v[3]); a3= angle(v[2],v[3]);
    min_max3(a1,a2,a3,min,max);
  }
  if (type_ == COM::Connectivity::QUAD4){
    double a1,a2,a3,a4;
    a1 = angle (v[0],v[1]); a2 = angle (v[0],v[2]); 
    a3 = angle (v[2],v[3]); a4 = angle (v[1],v[3]);
    min_max4(a1,a2,a3,a4,min,max); 
  }
}

void Geo_Metric_Base_2::compute_aspects(double& R, double& r, double& l) const {
  // Defined for Triangles only
  assert (type_ == COM::Connectivity::TRI3);

  double a = edge_length(v[0]);
  double b = edge_length(v[2]);
  double c = edge_length(v[3]);

  // Use formulas from mathworld to find circumradius (R),
  // inradius (r), and shortest edge length (l)

  l = ( a < b )? a : b;
  l = ( l < c ) ? l : c;
  R = (a*b*c) / sqrt( (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c) );
  r = .5*sqrt( ((b+c-a)*(c+a-b)*(a+b-c))/(a+b+c) );

}

void Angle_Metric_2::compute(double atts[]) const {
  compute_angles(atts[0],atts[1]);
}

double Angle_Metric_2::maxValue() const { 
return 180.0; 
}

double Angle_Metric_2::minValue() const { 
return 0.0; 
}

void Aspect_Metric_2::compute(double atts[]) const {
  if (type_ == COM::Connectivity::QUAD4) {
    atts[0] = -1;
    atts[1] = -1;
  }
  else{
    double R, r, l;
    compute_aspects(R,r,l);
    atts[0] = (r/R)*2.0;
    atts[1] = (l/R)*.5773502692;  
  }
}

double Aspect_Metric_2::maxValue() const { 
return 1.0; 
}

double Aspect_Metric_2::minValue() const { 
return 0.0; 
}

MOP_END_NAMESPACE






