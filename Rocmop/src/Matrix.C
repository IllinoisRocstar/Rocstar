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
/*! \file Matrix.C 
    \brief Implementation of a limited matrix class.

    Supports 3 by 2, and 3 by 3 matrices.
*/

#include "Matrix.h"
#include <cmath>
#include <cassert>
#include <iostream>

MOP_BEGIN_NAMESPACE

using namespace std;

// For debugging purposes
void Matrix::info(){
  cout << "Hi, I'm a matrix with " << rows_ << " rows and " << cols_ << " cols_" << endl;
  cout << "    DATA = "; 
    for (int i = 0; i < rows_ * cols_; i ++){
      cout <<data_[i] << " ";
    }
  cout << endl;
}

double &Matrix::operator() (int r, int c){
  assert ( 0 <= r && r <= rows_ && 0 <= c && c <= cols_);
  return data_[r*cols_ + c];
}

double Matrix::operator() (int r, int c) const {
  assert ( 0 <= r && r <= rows_ && 0 <= c && c <= cols_);
  return data_[r*cols_ + c];
}

const Matrix &Matrix::operator = (const Matrix &right){
  if (&right != this) {
    rows_ = right.rows_;
    cols_ = right.cols_;

    for (int r = 0; r < rows_; r ++){
      for (int c = 0; c < cols_; c ++ ){
	data_[r*cols_+c] = right.data_[r*cols_+c];
      }
    }
  }
  return *this;
}

// An operator used for the creation of the tensor matrix
// L ^= A sets L to A transpose times A
const Matrix &Matrix::operator ^= (const Matrix &right){
  if (&right != this) {
    rows_ = right.cols_;
    cols_ = right.cols_;

    for (int r = 0; r < rows_; r ++){
      for (int c = 0; c < rows_; c ++ ){
	data_[r*cols_+c] = 0;
	for (int j = 0; j < right.rows_; j ++){
	  data_[r*rows_+c] +=  right(j,r) * right(j,c) ;
	}
      }
    }
  }
  return *this;
}
  
// For debugging purposes
ostream &operator <<(ostream &output, const Matrix &m){
  for (int r = 0; r <m.rows_; r++){
    for (int c = 0; c < m.cols_; c++){
      output << m(r,c) << "  ";
    }
    output << endl;
  }
  return output;
}

// The Jacobian matrix class, for 3 by 2 and 3 by 3 matrices 
//J_Matrix::J_Matrix()
//  :Matrix()
//{
//}

J_Matrix::J_Matrix(Vector_3<double> pnts[],int cols)
{
  //cols is the number of columns in the resulting jacobian matrix
  //cols == 2 for triangles and quads, 3 for tets, or hexes
  rows_ = 3;
  cols_ = cols;

  for (int r = 0; r < cols_; r++){
    for (int c = 0; c < 3; c++){
      data_[c*cols_+r] = pnts[r][c];
    }
  }
}

double J_Matrix::det (){
//       in the case for size == 2 (jacobian matrix not square)
//       the determinant should be twice the area, so just caclulate
//       the area and double it.
//
//                                     2                     2                    2
//                        | y_0 z_0 1 |         | x_0 x_0 1 |        | x_0 y_0 1 |
//   area = 1/2 * sqrt (  | y_1 z_1 1 |     +   | z_1 x_1 1 |    +   | x_1 y_1 1 |  )
//                        | y_2 z_2 1 |         | z_2 x_2 1 |        | x_2 y_2 1 |
//
// see http://mathworld.wolfram.com/TriangleArea.html for more details
  assert( (rows_ == 3) && (cols_ == 2 || cols_ ==3) );
  if (cols_ == 2){
    //return (*this)(0,0)*(*this)(1,1) - (*this)(0,1) * (*this)(1,0);
    double s1 = (*this)(1,0)*(*this)(2,1) - (*this)(1,1)*(*this)(2,0);
    double s2 = (*this)(2,0)*(*this)(0,1) - (*this)(2,1)*(*this)(0,0);
    double s3 = (*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0);
    s1 *= s1;
    s2 *= s2;
    s3 *= s3;
    return sqrt ( s1 + s2 + s3 );
  }
  else
    return 
	 (*this)(0,0) * ((*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1)) -
	 (*this)(0,1) * ((*this)(1,0)*(*this)(2,2)-(*this)(1,2)*(*this)(2,0)) +
	 (*this)(0,2) * ((*this)(1,0)*(*this)(2,1)-(*this)(1,1)*(*this)(2,0));
}

MOP_END_NAMESPACE






