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
// $Id: metric_demo.C,v 1.5 2008/12/06 08:45:25 mtcampbe Exp $

#include "Algebraic_Metrics_2.h"
#include "Algebraic_Metrics_3.h"
#include "Geometric_Metrics_2.h"
#include "Geometric_Metrics_3.h"
#include "Connectivity.hpp"
#include <iostream>

using namespace std;
USE_MOP_NAMESPACE;

int main (){
  // Build the reference tetrahedron and triangle
  Vector_3<double> v0, v1, v2 , v3;
  v0[0] = 0.0; v0[1] = 0.0;     v0[2] = 0.0;
  v1[0] = 1.0; v1[1] = 0.0;     v1[2] = 0.0;
  v2[0] = 0.5; v2[1] = .866025; v2[2] = 0.0;
  v3[0] = 0.5; v3[1] = .288675; v3[2] = .816497;
  Vector_3<double> ref_tet[4] = {v0,v1,v2,v3};
  Vector_3<double> ref_tri[3] = {v0,v1,v2};

  // Build the logical tetrahedron and triangle (nodes are at position 1 or 0 along each axis)
  Vector_3<double> p0, p1, p2 , p3;
  p0[0] = 0.0; p0[1] = 0.0;     p0[2] = 0.0;
  p1[0] = 1.0; p1[1] = 0.0;     p1[2] = 0.0;
  p2[0] = 0.0; p2[1] = 1.0;     p2[2] = 0.0;
  p3[0] = 0.0; p3[1] = 0.0;     p3[2] = 1.0;
  Vector_3<double> log_tet[4] = {p0,p1,p2,p3};
  Vector_3<double> log_tri[3] = {p0,p1,p2,};

  double tet_shape[1];
  double tet_size[1];
  double tet_skew[1];
  double tet_angles[2];
  double tet_aspects[2];

  double tri_shape[1];
  double tri_size[1];
  double tri_skew[1];
  double tri_angles[2];
  double tri_aspects[2];

  // For the reference volume I used the volume of a regular tetrahedron with
  // edges of length 1.  

  Shape_Metric_3 Shape3;
  Size_Metric_3  Size3(.1178511);
  Skew_Metric_3  Skew3;
  Angle_Metric_3 Angle3;
  Aspect_Metric_3 Aspect3;

  Shape_Metric_2 Shape2;
  Size_Metric_2  Size2(.4330127);
  Skew_Metric_2  Skew2;
  Angle_Metric_2 Angle2;
  Aspect_Metric_2 Aspect2;

  // Initialize the metrics

  Shape3.initialize(ref_tet,COM::Connectivity::TET4);
  Size3.initialize(ref_tet,COM::Connectivity::TET4);
  Skew3.initialize(ref_tet,COM::Connectivity::TET4);
  Angle3.initialize(ref_tet,COM::Connectivity::TET4);
  Aspect3.initialize(ref_tet,COM::Connectivity::TET4);

  Shape2.initialize(ref_tet,COM::Connectivity::TRI3);
  Size2.initialize(ref_tet,COM::Connectivity::TRI3);
  Skew2.initialize(ref_tet,COM::Connectivity::TRI3);
  Angle2.initialize(ref_tri,COM::Connectivity::TRI3);
  Aspect2.initialize(ref_tri,COM::Connectivity::TRI3);

  // Compute the metric values

  Shape3.compute(tet_shape);
  Size3.compute(tet_size);
  Skew3.compute(tet_skew);
  Angle3.compute(tet_angles);
  Aspect3.compute(tet_aspects);

  Shape2.compute(tri_shape);
  Size2.compute(tri_size);
  Skew2.compute(tri_skew);
  Angle2.compute(tri_angles);
  Aspect2.compute(tri_aspects);

  // Spit out some junk
  cout << "For the regular tetrahedron, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << tet_shape[0]   << endl;
  cout << "    size        = " << tet_size[0]    << endl;
  cout << "    skew        = " << tet_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << tet_angles[0]  << endl;
  cout << "    max angle   = " << tet_angles[1]  << endl;
  cout << "    R/r         = " << tet_aspects[0] << endl;
  cout << "    R/l         = " << tet_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape3.minValue() << "," << Shape3.maxValue() << "]" << endl;
  cout << "  size      = [" << Size3.minValue() << "," << Size3.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew3.minValue() << "," << Skew3.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle3.minValue() << "," << Angle3.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect3.minValue() << "," << Aspect3.maxValue()   << "]" << endl << endl;
 
  cout << "For the equilateral triangle, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << tri_shape[0]   << endl;
  cout << "    size        = " << tri_size[0]    << endl;
  cout << "    skew        = " << tri_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << tri_angles[0]  << endl;
  cout << "    max angle   = " << tri_angles[1]  << endl;
  cout << "    R/r         = " << tri_aspects[0] << endl;
  cout << "    R/l         = " << tri_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape3.minValue() << "," << Shape3.maxValue() << "]" << endl;
  cout << "  size      = [" << Size3.minValue() << "," << Size3.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew3.minValue() << "," << Skew3.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle3.minValue() << "," << Angle3.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect3.minValue() << "," << Aspect3.maxValue()   << "]" << endl << endl;

  // Initialize the metrics

  Shape3.initialize(log_tet,COM::Connectivity::TET4);
  Size3.initialize(log_tet,COM::Connectivity::TET4);
  Skew3.initialize(log_tet,COM::Connectivity::TET4);
  Angle3.initialize(log_tet,COM::Connectivity::TET4);
  Aspect3.initialize(log_tet,COM::Connectivity::TET4);

  Shape2.initialize(log_tet,COM::Connectivity::TRI3);
  Size2.initialize(log_tet,COM::Connectivity::TRI3);
  Skew2.initialize(log_tet,COM::Connectivity::TRI3);
  Angle2.initialize(log_tri,COM::Connectivity::TRI3);
  Aspect2.initialize(log_tri,COM::Connectivity::TRI3);

  // Compute the metric values

  Shape3.compute(tet_shape);
  Size3.compute(tet_size);
  Skew3.compute(tet_skew);
  Angle3.compute(tet_angles);
  Aspect3.compute(tet_aspects);

  Shape2.compute(tri_shape);
  Size2.compute(tri_size);
  Skew2.compute(tri_skew);
  Angle2.compute(tri_angles);
  Aspect2.compute(tri_aspects);

  // Spit out more junk
  cout << "For the logical tetrahedron, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << tet_shape[0]   << endl;
  cout << "    size        = " << tet_size[0]    << endl;
  cout << "    skew        = " << tet_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << tet_angles[0]  << endl;
  cout << "    max angle   = " << tet_angles[1]  << endl;
  cout << "    R/r         = " << tet_aspects[0] << endl;
  cout << "    R/l         = " << tet_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape3.minValue() << "," << Shape3.maxValue() << "]" << endl;
  cout << "  size      = [" << Size3.minValue() << "," << Size3.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew3.minValue() << "," << Skew3.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle3.minValue() << "," << Angle3.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect3.minValue() << "," << Aspect3.maxValue()   << "]" << endl << endl;
 
  cout << "For the logical triangle, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << tri_shape[0]   << endl;
  cout << "    size        = " << tri_size[0]    << endl;
  cout << "    skew        = " << tri_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << tri_angles[0]  << endl;
  cout << "    max angle   = " << tri_angles[1]  << endl;
  cout << "    R/r         = " << tri_aspects[0] << endl;
  cout << "    R/l         = " << tri_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape3.minValue() << "," << Shape3.maxValue() << "]" << endl;
  cout << "  size      = [" << Size3.minValue() << "," << Size3.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew3.minValue() << "," << Skew3.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle3.minValue() << "," << Angle3.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect3.minValue() << "," << Aspect3.maxValue()   << "]" << endl << endl;


  // Build the reference/logical hexahedron=
  Vector_3<double> v2_0, v2_1, v2_2 ,v2_3, v2_4, v2_5, v2_6, v2_7;
  v2_0[0] = 0.0; v2_0[1] = 0.0;     v2_0[2] = 0.0;
  v2_1[0] = 1.0; v2_1[1] = 0.0;     v2_1[2] = 0.0;
  v2_2[0] = 1.0; v2_2[1] = 1.0;     v2_2[2] = 0.0;
  v2_3[0] = 0.0; v2_3[1] = 1.0;     v2_3[2] = 0.0;
  v2_4[0] = 0.0; v2_4[1] = 0.0;     v2_4[2] = 1.0;
  v2_5[0] = 1.0; v2_5[1] = 0.0;     v2_5[2] = 1.0;
  v2_6[0] = 1.0; v2_6[1] = 1.0;     v2_6[2] = 1.0;
  v2_7[0] = 0.0; v2_7[1] = 1.0;     v2_7[2] = 1.0;

  Vector_3<double> ref_hex[8] = {v2_0,v2_1,v2_2,v2_3,v2_4,v2_5,v2_6,v2_7};
  Vector_3<double> ref_quad[4] = {v2_0,v2_1,v2_2,v2_3};

  Size_Metric_3  Size_3(1.0);
  Size_Metric_2  Size_2(1.0);

  double hex_shape[1];
  double hex_size[1];
  double hex_skew[1];
  double hex_angles[2];
  double hex_aspects[2];

  double quad_shape[1];
  double quad_size[1];
  double quad_skew[1];
  double quad_angles[2];
  double quad_aspects[2];

  // Initialize the metrics
  Shape3.initialize(ref_hex,COM::Connectivity::HEX8);
  Size_3.initialize(ref_hex,COM::Connectivity::HEX8);
  Skew3.initialize(ref_hex,COM::Connectivity::HEX8);
  Angle3.initialize(ref_hex,COM::Connectivity::HEX8);
  Aspect3.initialize(ref_hex,COM::Connectivity::HEX8);

  Shape2.initialize(ref_quad,COM::Connectivity::QUAD4);
  Size_2.initialize(ref_quad,COM::Connectivity::QUAD4);
  Skew2.initialize(ref_quad,COM::Connectivity::QUAD4);
  Angle2.initialize(ref_quad,COM::Connectivity::QUAD4);
  Aspect2.initialize(ref_quad,COM::Connectivity::QUAD4);

  // Compute values
  Shape3.compute(hex_shape);
  Size_3.compute(hex_size);
  Skew3.compute(hex_skew);
  Angle3.compute(hex_angles);
  Aspect3.compute(hex_aspects);

  Shape2.compute(quad_shape);
  Size_2.compute(quad_size);
  Skew2.compute(quad_skew);
  Angle2.compute(quad_angles);
  Aspect2.compute(quad_aspects);

  // Spit out junk
  cout << "For the regular hexahedron, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << hex_shape[0]   << endl;
  cout << "    size        = " << hex_size[0]    << endl;
  cout << "    skew        = " << hex_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << hex_angles[0]  << endl;
  cout << "    max angle   = " << hex_angles[1]  << endl;
  cout << "    R/r         = " << hex_aspects[0] << endl;
  cout << "    R/l         = " << hex_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape3.minValue() << "," << Shape3.maxValue() << "]" << endl;
  cout << "  size      = [" << Size_3.minValue() << "," << Size_3.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew3.minValue() << "," << Skew3.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle3.minValue() << "," << Angle3.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect3.minValue() << "," << Aspect3.maxValue()   << "]" << endl << endl;
 
  cout << "For the square, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << quad_shape[0]   << endl;
  cout << "    size        = " << quad_size[0]    << endl;
  cout << "    skew        = " << quad_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << quad_angles[0]  << endl;
  cout << "    max angle   = " << quad_angles[1]  << endl;
  cout << "    R/r         = " << quad_aspects[0] << endl;
  cout << "    r/l         = " << quad_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape2.minValue() << "," << Shape2.maxValue() << "]" << endl;
  cout << "  size      = [" << Size2.minValue() << "," << Size2.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew2.minValue() << "," << Skew2.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle2.minValue() << "," << Angle2.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect2.minValue() << "," << Aspect2.maxValue()   << "]" << endl << endl;

// Build a non-ideal hexahedron and square
  Vector_3<double> p2_0, p2_1, p2_2 ,p2_3, p2_4, p2_5, p2_6, p2_7;
  p2_0[0] = 0.0; p2_0[1] = 0.0;     p2_0[2] = 0.0;
  p2_1[0] = .9; p2_1[1] = 0.0;     p2_1[2] = 0.1;
  p2_2[0] = 1.0; p2_2[1] = 1.1;     p2_2[2] = 0.0;
  p2_3[0] = 0.2; p2_3[1] = .85;     p2_3[2] = 0.0;
  p2_4[0] = 0.0; p2_4[1] = 0.0;     p2_4[2] = 1.1;
  p2_5[0] = 1.0; p2_5[1] = 0.0;     p2_5[2] = 1.0;
  p2_6[0] = 1.2; p2_6[1] = .9;     p2_6[2] = 1.0;
  p2_7[0] = 0.2; p2_7[1] = 1.0;     p2_7[2] = 1.2;

  Vector_3<double> misc_hex[8] = {p2_0,p2_1,p2_2,p2_3,p2_4,p2_5,p2_6,p2_7};
  Vector_3<double> misc_quad[4] = {p2_0,p2_1,p2_2,p2_3};

  // And so on, and so forth...
  Shape3.initialize(misc_hex,COM::Connectivity::HEX8);
  Size_3.initialize(misc_hex,COM::Connectivity::HEX8);
  Skew3.initialize(misc_hex,COM::Connectivity::HEX8);
  Angle3.initialize(misc_hex,COM::Connectivity::HEX8);
  Aspect3.initialize(misc_hex,COM::Connectivity::HEX8);

  Shape2.initialize(misc_quad,COM::Connectivity::QUAD4);
  Size_2.initialize(misc_quad,COM::Connectivity::QUAD4);
  Skew2.initialize(misc_quad,COM::Connectivity::QUAD4);
  Angle2.initialize(misc_quad,COM::Connectivity::QUAD4);
  Aspect2.initialize(misc_quad,COM::Connectivity::QUAD4);

  Shape3.compute(hex_shape);
  Size_3.compute(hex_size);
  Skew3.compute(hex_skew);
  Angle3.compute(hex_angles);
  Aspect3.compute(hex_aspects);

  Shape2.compute(quad_shape);
  Size_2.compute(quad_size);
  Skew2.compute(quad_skew);
  Angle2.compute(quad_angles);
  Aspect2.compute(quad_aspects);

  cout << "For the non-regular hexahedron, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << hex_shape[0]   << endl;
  cout << "    size        = " << hex_size[0]    << endl;
  cout << "    skew        = " << hex_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << hex_angles[0]  << endl;
  cout << "    max angle   = " << hex_angles[1]  << endl;
  cout << "    R/r         = " << hex_aspects[0] << endl;
  cout << "    R/l         = " << hex_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape3.minValue() << "," << Shape3.maxValue() << "]" << endl;
  cout << "  size      = [" << Size_3.minValue() << "," << Size_3.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew3.minValue() << "," << Skew3.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle3.minValue() << "," << Angle3.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect3.minValue() << "," << Aspect3.maxValue()   << "]" << endl << endl;
 
  cout << "For the non-regular quadrilateral, metric values are: " << endl;
  cout << "  Algebraic Metrics: " << endl;
  cout << "    shape       = " << quad_shape[0]   << endl;
  cout << "    size        = " << quad_size[0]    << endl;
  cout << "    skew        = " << quad_skew[0]    << endl;
  cout << "  Geometric Metrics: " << endl;
  cout << "    min angle   = " << quad_angles[0]  << endl;
  cout << "    max angle   = " << quad_angles[1]  << endl;
  cout << "    R/r         = " << quad_aspects[0] << endl;
  cout << "    r/l         = " << quad_aspects[1] << endl;
  cout << "The range of values for the metrics are: " << endl;
  cout << "  shape     = [" << Shape2.minValue() << "," << Shape2.maxValue() << "]" << endl;
  cout << "  size      = [" << Size2.minValue() << "," << Size2.maxValue()   << "]" << endl;
  cout << "  skew      = [" << Skew2.minValue() << "," << Skew2.maxValue()   << "]" << endl; 
  cout << "  angles    = [" << Angle2.minValue() << "," << Angle2.maxValue()   << "]" << endl;
  cout << "  a. ratios = [" << Aspect2.minValue() << "," << Aspect2.maxValue()   << "]" << endl << endl;

  return 0;
}






