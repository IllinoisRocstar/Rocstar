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
// $Id: Matrix.h,v 1.4 2008/12/06 08:45:23 mtcampbe Exp $

/*! \file Matrix.h 
    \brief Declaration of Matrix class.
*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "mopbasic.h"

MOP_BEGIN_NAMESPACE


//! A simple Matrix class.
/**
 * This class provides only those functions required to implement Knupp's
 * Algebriac quality metrics.
 */
class Matrix{
  friend std::ostream &operator <<(std::ostream &, const Matrix &);
public:
  
  //! Constructor
  Matrix() {}

  //! Destructor
  ~Matrix() {}

  //! Matrix data access by reference
  double&  operator() (int, int);

  //! Matrix data access by value
  double   operator() (int, int) const;

  //! Equal operator
  const Matrix &operator = (const Matrix &);

  //! Shortcut for initializing a metric tensor.
  /**
   *  If T is a Matrix, and A is a J_Matrix, then to
   *  initialize \f$T\f$ to \f$A^T\f$, use 
   *  \code T ^= A;
   */
  const Matrix &operator ^= (const Matrix &);

  //! Print out matrix information for debugging.
  void info ();
protected:
  int     rows_;
  int     cols_;
  double  data_[9];
};

//! Jacobian Matrix Class
/**
 * A Jacobian Matrix Class supporting 3 by 3 matrices for tetrahedrons 
 * and hexahedrons and 3 by 2 matrices for triangles and quadrilaterals.
 */ 
class J_Matrix : public Matrix {

public:
  
  // Constructor
  J_Matrix() {}

  // Destructor
  ~J_Matrix() {}

  // Constructor from noal coords and number of columns
  //  
  //  \param pnts[] A Vector_3 array of ordered nodal points.
  //  \param cols   Columns in the Jacobian.  Should be 2 for triangles or quadrilaterals
  //		    and 3 for tetrahedrons or hexahedrons.
  J_Matrix(Vector_3<double> pnts[],int cols);

  //! Obtain the determinant
  /** 
   * Only defined for 3 by 3, and 3 by 2 matrices.  In the case of a 3 by 2 jacobian 
   * matrix, returns twice the area of the triangle defined by the matrix.
   */
  double det();
};

MOP_END_NAMESPACE

#endif






