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
/*! \file Algebraic_Metrics_2.C 
    \brief Implementation of 2D metric quality measures.

    Implementation of the algebraic mesh quality metrics outlined by Patrick M. 
    Knupp in "Algebraic mesh quality metrics for unstructured initial meshes."
    Finite Elements in Analysis and Design 39 (2003) 217-241.
*/

#include "Algebraic_Metrics_2.h"
#include <vector>
#include <cassert>
#include <iostream>

using namespace std;

MOP_BEGIN_NAMESPACE

void Alg_Metric_Base_2::initialize(Element_node_enumerator &ene){
  type_ = ene.type();
  Element_node_vectors_k_const<double> n;
  n.set(ene.pane()->dataitem("nc"),ene);
  vector<Vector_3<double> > v(2);
  if (type_ == COM::Connectivity::TRI3){
    for ( int i = 0; i < 3; i ++){
      int i1 = (i+1)%3, i2 = (i+2)%3;

      v[0][0] = n(i1,0)-n(i,0); 
      v[0][1] = n(i1,1)-n(i,1); 
      v[0][2] = n(i1,2)-n(i,2); 

      v[1][0] = n(i2,0)-n(i,0); 
      v[1][1] = n(i2,1)-n(i,1); 
      v[1][2] = n(i2,2)-n(i,2); 

      A[i] = J_Matrix(&v[0],2);
      L[i] ^= A[i];
      alpha[i] = A[i].det();
    }
  }
  else if (type_ == COM::Connectivity::QUAD4){
    for ( int i = 0; i < 4; i ++){
      for (int j = 0; j < 3; j ++){
	int i1 = (i+1)%4, i3 = (i+3)%4;

	v[0][j] = n(i1,j)-n(i,j); 	
	v[1][j] = n(i3,0)-n(i,j); 
      }
	A[i] = J_Matrix(&v[0],2);
	L[i] ^= A[i];
	alpha[i] = A[i].det();
    }
  }
  else COM_assertion_msg(0,"Element type not supported for 2D Algebraic Metrics");
}


void Alg_Metric_Base_2::initialize(Vector_3<double> n[], int type){
  type_ = type;
  vector<Vector_3<double> > v(type_);
  if (type_ == COM::Connectivity::TRI3){
    for ( int i = 0; i < 3; i ++){
      v[0] = n[(i+1)%3] - n[i];
      v[1] = n[(i+2)%3] - n[i];
      A[i] = J_Matrix(&v[0],2);
      L[i] ^= A[i];
      alpha[i] = A[i].det();
    }
  }
  else if (type_ == COM::Connectivity::QUAD4){
    for ( int i = 0; i < 4; i ++){
      for (int j = 0; j < 3; j ++){
	v[0][j] = (n[(i+1)%4][j] - n[i][j]);
	v[1][j] = (n[(i+3)%4][j] - n[i][j]);
      }
	A[i] = J_Matrix(&v[0],2);
	L[i] ^= A[i];
	alpha[i] = A[i].det();
    }
  }
}

double Alg_Metric_Base_2::compute_size( double ref_area) const {
  double tau;
  if (type_ == COM::Connectivity::TRI3){
    tau = alpha[0]/(2.0*ref_area);
    return (tau < (1/tau))? tau : (1/tau);
  }
  else {
    assert (type_ == COM::Connectivity::QUAD4);
    tau = (alpha[0] + alpha[2])/(2*ref_area);
    return (tau < (1/tau))? tau : (1/tau);
  }
}

double Alg_Metric_Base_2::compute_shape() const {
  // Triangles
  if (type_ == COM::Connectivity::TRI3){
    return sqrt(3.0) * alpha[0] / 
      ( L[0](0,0) + L[0](1,1) - L[0](0,1));
  }

  // Quads
  else{
    assert (type_ == COM::Connectivity::QUAD4);
    double denom = 0;
    for (int i = 0; i < 4; i ++){
      denom += ( (L[i](0,0) + L[i](1,1))/alpha[i] );
    }
    return 8.0/denom;
  }
}

double Alg_Metric_Base_2::compute_skew() const {
  // only defined for quads, so return 1 if triangle
  if (type_ == COM::Connectivity::TRI3) {
    return 1;
  }
  // flag checks for divide by 0
  else {
    double denom = 0;
    int flag = 0;
    for (int i = 0; i < 4 ; i++) {
      if (alpha[i]==0) { flag = 1;}
      else {
	denom += ( sqrt ( L[i](0,0) * L[i](1,1)) / alpha[i] );
      }
    }
    if (flag){ return 0; }
    else     { return  4.0/denom; } 
  }
}

void Shape_Metric_2::compute(double atts[]) const {
  atts[0] = compute_shape();
}

void Size_Metric_2::compute(double atts[]) const {
  atts[0] = compute_size(ref_area);
}

void Size_Shape_Metric_2::compute(double atts[]) const {
  atts[0] = compute_size(ref_area) * compute_shape();
}

void Skew_Metric_2::compute(double atts[]) const {
  atts[0] = compute_skew();
}

void Size_Skew_Metric_2::compute(double atts[]) const {
  atts[0] = compute_size(ref_area) * compute_skew();
}

ostream & operator<<(ostream &output, const Alg_Metric_Base_2 &m){
  for(int i =0; i < m.type_; i ++){
    output << "A[" << i << "]:" << endl << m.A[i];
    output << "L[" << i << "]:" << endl << m.L[i];
    output << "alpha"<<i << ": " << m.alpha[i] << endl << endl;
  }
  return output;
}  

MOP_END_NAMESPACE






