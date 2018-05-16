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
/*! \file Algebraic_Metrics_3.C 
    \brief Implementation of 3D metric quality measures.

    Implementation of the algebraic mesh quality metrics outlined by Patrick M. 
    Knupp in "Algebraic mesh quality metrics for unstructured initial meshes."
    Finite Elements in Analysis and Design 39 (2003) 217-241.
*/

#include "Algebraic_Metrics_3.h"
#include <vector>
#include <cassert>
#include <iostream>

MOP_BEGIN_NAMESPACE

using namespace std;

void Alg_Metric_Base_3::initialize(Element_node_enumerator &ene){
  type_ = ene.type();
  Element_node_vectors_k_const<double> n;
  n.set(ene.pane()->dataitem("nc"),ene);
  vector<Vector_3<double> > v(2);
  if (type_ == COM::Connectivity::TET4){
    double premult = 1.0;
    for ( int i = 0; i < 4; i ++){
      for (int j = 0; j < 3; j ++){
	v[0][j] = premult * (n((i+1)%4,j) - n(i,j));
	v[1][j] = premult * (n((i+2)%4,j) - n(i,j));
	v[2][j] = premult * (n((i+3)%4,j) - n(i,j));
      }
      A[i] = J_Matrix(&v[0],3);
      L[i] ^= A[i];
      alpha[i] = A[i].det();
      premult *= -1.0;
    }
  }
  else if (type_ == COM::Connectivity::HEX8){
    double premult = 1.0;
    for ( int i = 0; i < 8; i ++){

      if ( (i==0) || (i==3) || (i==5) || (i==6) )      { premult = 1.0;  }
      else if ( (i==1) || (i==2) || (i==4) || (i==7) ) { premult = -1.0; }
      int a,b,c;
      if ( (i == 0) || (i == 4) )   { a=1; b=3; c=4; }
      else if ( (i ==3) || (i==7) ) { a=4;b=5;c=7;   }
      else {a=1;b=4;c=7;}

      for (int j = 0; j < 3; j ++){
	v[0][j] = premult*( n((i+a)%8,j) - n(i,j));
	v[1][j] = premult*( n((i+b)%8,j) - n(i,j));
	v[2][j] = premult*( n((i+c)%8,j) - n(i,j));
      }
      A[i] = J_Matrix(&v[0],3);
      L[i] ^= A[i];
      alpha[i] = A[i].det();
    }
  }
}

void Alg_Metric_Base_3::initialize(Vector_3<double> n[], int type){
  type_ = type;
  vector<Vector_3<double> > v(3);
  if (type_ == COM::Connectivity::TET4){
    double premult = 1.0;
    for ( int i = 0; i < 4; i ++){
      for (int j = 0; j < 3; j ++){
	v[0][j] = premult * (n[(i+1)%4][j] - n[i][j]);
	v[1][j] = premult * (n[(i+2)%4][j] - n[i][j]);
	v[2][j] = premult * (n[(i+3)%4][j] - n[i][j]);
      }
      A[i] = J_Matrix(&v[0],3);
      L[i] ^= A[i];
      alpha[i] = A[i].det();
      premult *= -1.0;
    }
  }
  else if (type_ == COM::Connectivity::HEX8){
    double premult = 1.0;
    for ( int i = 0; i < 8; i ++){

      if ( (i==0) || (i==3) || (i==5) || (i==6) )      { premult = 1.0;  }
      else if ( (i==1) || (i==2) || (i==4) || (i==7) ) { premult = -1.0; }
      int a,b,c;
      if ( (i == 0) || (i == 4) )   { a=1; b=3; c=4; }
      else if ( (i ==3) || (i==7) ) { a=4;b=5;c=7;   }
      else {a=1;b=4;c=7;}

      for (int j = 0; j < 3; j ++){
	v[0][j] = premult*( n[(i+a)%8][j] - n[i][j]);
	v[1][j] = premult*( n[(i+b)%8][j] - n[i][j]);
	v[2][j] = premult*( n[(i+c)%8][j] - n[i][j]);
      }
      A[i] = J_Matrix(&v[0],3);
      L[i] ^= A[i];
      alpha[i] = A[i].det();
    }
  }
}

double Alg_Metric_Base_3::compute_size( double ref_vol) const {
  double size = 0.0;
  if (type_ == COM::Connectivity::TET4){
    size = (alpha[0]/6.0)/ref_vol;
    size = (size < (1/size))? size : (1/size);
  }
  if (type_ == COM::Connectivity::HEX8){
    for (int i = 0; i < 8; i++){
      size += alpha[i];
    }
    if (size  < 0) {size = 0;}
    else {
      size = size / ( 8.0 * ref_vol );
      size = (size < (1/size))? size: (1/size);
    }
  }
  return size;
}

double Alg_Metric_Base_3::compute_shape() const {
  if (type_ == COM::Connectivity::TET4){
    double num = (3.0 * pow(alpha[0]*sqrt(2.0),(2.0/3.0)) );
    double denom = 1.5*(L[0](0,0)+L[0](1,1)+L[0](2,2))-
                       (L[0](0,1)+L[0](1,2)+L[0](0,2));
    return num/denom; 
  }
  else{
    assert (type_ == COM::Connectivity::HEX8);
    double denom = 0;
    for (int i = 0; i < 8; i ++){
      denom += (L[i](0,0) + L[i](1,1) + L[i](2,2))/(pow(abs(alpha[i]),(2.0/3.0)));
    }
    return  24.0/denom;
  }
}

double Alg_Metric_Base_3::compute_skew() const {
  // Not defined for Tets, so we return 1
  if (type_ == COM::Connectivity::TET4) {
    return 1.0;
  }

  // For quads, set flag in case of divide by zero (degenerate element)
  else{
    assert (type_ == COM::Connectivity::HEX8);
    double denom = 0.0;
    int flag = 0;
    for (int i = 0; i < 8; i ++){
      double root = sqrt(L[i](0,0)*L[i](1,1)*L[i](2,2));
      if (root == 0) { flag = 1;}
      denom += pow( abs(root/alpha[i]), (2.0/3.0));
    }
    if (flag) { return  0.0;}
    else { return 8.0/denom; }
  }
}

void Shape_Metric_3::compute(double atts[]) const {
  atts[0] = compute_shape();
}

void Size_Metric_3::compute(double atts[]) const {
  atts[0] = compute_size(ref_vol);
}

void Size_Shape_Metric_3::compute(double atts[]) const {
  atts[0] = compute_size(ref_vol) * compute_shape();
}

void Skew_Metric_3::compute(double atts[]) const {
  atts[0] = compute_skew();
}

void Size_Skew_Metric_3::compute(double atts[]) const {
  atts[0] = compute_size(ref_vol) * compute_skew();
}

ostream & operator<<(ostream &output, const Alg_Metric_Base_3 &m){
  for(int i =0; i < (m.type_-100); i ++){
    output << "A[" << i << "]:" << endl << m.A[i];
    output << "L[" << i << "]:" << endl << m.L[i];
    output << "alpha"<<i << ": " << m.alpha[i] << endl << endl;
  }
  return output;
}  

MOP_END_NAMESPACE






