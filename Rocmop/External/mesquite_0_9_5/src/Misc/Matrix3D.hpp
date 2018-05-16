/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 18-Dec-02 at 11:08:22
//  LAST-MOD: 27-May-04 at 14:48:56 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file Matrix3D.hpp

3*3 Matric class, row-oriented, 0-based [i][j] indexing.

 \author Thomas Leurent
 
*/
// DESCRIP-END.
//



#ifndef Matrix3D_hpp
#define Matrix3D_hpp

#ifndef MSQ_USE_OLD_IO_HEADERS
#include <iostream>
#include <sstream>
#else
#include <iostream.h>
#include <strstream.h>
#endif

#ifndef MSQ_USE_OLD_C_HEADERS
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include "Mesquite.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{

  /*! \class Matrix3D
      \brief 3*3 Matric class, row-oriented, 0-based [i][j] indexing.

      Since the size of the object is fixed at compile time, the Matrix3D
      object is as fast as a double[9] array.
  */
  class Matrix3D 
  {
  protected:
    double v_[9];                  

   
    void copy(const double*  v)
    { memcpy(v_, v, 9*sizeof(double)); }

    void set(const double& val)
    {
      v_[0]=val;  v_[1]=val;  v_[2]=val;
      v_[3]=val;  v_[4]=val;  v_[5]=val;
      v_[6]=val;  v_[7]=val;  v_[8]=val;
    }

    void set_values(const char *s)
    {
#ifdef MSQ_USE_OLD_IO_HEADERS
      ::istrstream ins(s);
#else
      std::istringstream ins(s);
#endif
      ins>>v_[0];  ins>>v_[1];  ins>>v_[2]; 
      ins>>v_[3];  ins>>v_[4];  ins>>v_[5]; 
      ins>>v_[6];  ins>>v_[7];  ins>>v_[8]; 
    }
    
  public:

    // constructors
    //! Default constructor sets all entries to 0. 
    Matrix3D()
    {
      zero();
    }
    
    Matrix3D(const Matrix3D &A)
    {
      copy(A.v_);
    }

    //! sets all entries of the matrix to value.
    Matrix3D(const double& value)
    {
      set(value);
    }

    //! sets matrix entries to values in array.
    //! \param v is an array of 9 doubles. 
    Matrix3D(const double* v)
    {
      copy(v);
    }

    //! for test purposes, matrices can be instantiated as
    //! \code Matrix3D A("3 2 1  4 5 6  9 8 7"); \endcode
    Matrix3D(const char *s)
    {
      set_values(s);
    }

    // destructor
    ~Matrix3D() { }

    // assignments
    Matrix3D& operator=(const Matrix3D &A)
    {
      if (v_ == A.v_)
        return *this;
      copy(A.v_);
      return *this;
    }
        
    Matrix3D& operator=(const double& scalar)
    { 
      set(scalar); 
      return *this;
    }

    //! for test purposes, matrices can be assigned as follows
    //! \code A = "3 2 1  4 5 6  9 8 7"; \endcode
    Matrix3D& operator=(const char* s)
    { 
      set_values(s); 
      return *this;
    }

    //! Sets all entries to zero (more efficient than assignement).
    void zero()
    {
      v_[0]=0.;  v_[1]=0.;  v_[2]=0.;
      v_[3]=0.;  v_[4]=0.;  v_[5]=0.;
      v_[6]=0.;  v_[7]=0.;  v_[8]=0.;
    }
     
    //! Sets column j (0, 1 or 2) to Vector3D c.
    void set_column(int j, const Vector3D& c)
    {
      v_[0+j]=c[0];
      v_[3+j]=c[1];
      v_[6+j]=c[2];
    }
    
    //! returns the column length -- i is 0-based. 
    double column_length(int i) const 
    { return sqrt( v_[0+i]*v_[0+i] + v_[3+i]*v_[3+i] + v_[6+i]*v_[6+i] ); }

    // Matrix Operators
    friend bool operator==(const Matrix3D &lhs, const Matrix3D &rhs);
    friend bool operator!=(const Matrix3D &lhs, const Matrix3D &rhs);
    friend double Frobenius_2(const Matrix3D &A);
    friend Matrix3D transpose(const Matrix3D &A);
    friend const Matrix3D operator+(const Matrix3D &A, const Matrix3D &B);
    friend const Matrix3D operator-(const Matrix3D &A, const Matrix3D &B) ;
    friend const Matrix3D operator*(const Matrix3D &A, const Matrix3D &B);
    friend const Matrix3D mult_element(const Matrix3D &A, const Matrix3D &B);
    friend int matmult(Matrix3D& C, const Matrix3D  &A, const Matrix3D &B);
    friend const Vector3D operator*(const Matrix3D  &A, const Vector3D &x);
    friend const Vector3D operator*(const Vector3D &x, const Matrix3D  &A);
    const Matrix3D operator*(const double &s) const;
    friend const Matrix3D operator*(const double &s, const Matrix3D &A);    
    void operator+=(const Matrix3D &rhs);
    void operator-=(const Matrix3D &rhs);
    void operator*=(const double &s);
    Matrix3D plus_transpose(const Matrix3D &B) const;
    void plus_transpose_equal(const Matrix3D &B);
    Matrix3D& outer_product(const Vector3D &v1, const Vector3D &v2);
    void fill_lower_triangle();

    //! \f$ v = A*x \f$
    friend void eqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    //! \f$ v += A*x \f$
    friend void plusEqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    //! \f$ v += A^T*x \f$
    friend void plusEqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
     
    //! \f$ B += a*A \f$
    friend void plusEqaA(Matrix3D& B, const double a, const Matrix3D &A);

    //! determinant of matrix A, det(A).
    friend double det(const Matrix3D &A);

    //! \f$ B = A^{-1} \f$
    friend void inv(Matrix3D& B, const Matrix3D &A);

    //! \f$ B *= A^{-1} \f$
    friend void timesInvA(Matrix3D& B, const Matrix3D &A);

    //! \f$ Q*R = A \f$
    friend void QR(Matrix3D &Q, Matrix3D &R, const Matrix3D &A);

    size_t num_rows() const { return 3; }
    size_t num_cols() const { return 3; }

    //! returns a pointer to a row.
    inline double* operator[](unsigned i)
    {
      return v_ + 3*i;
    }

    //! returns a pointer to a row.
    inline const double* operator[](unsigned i) const
    {
      return v_ + 3*i;
    }

  };


  /* ***********  I/O  **************/

  inline msq_stdio::ostream& operator<<(msq_stdio::ostream &s, const Matrix3D &A)
  {
    for (size_t i=0; i<3; ++i)
      {
        for (size_t j=0; j<3; ++j)
          s << A[i][j] << " ";
        s << "\n";
      }
    return s;
  }

  inline msq_stdio::istream& operator>>(msq_stdio::istream &s, Matrix3D &A)
  {
    for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
        {
          s >>  A[i][j];
        }
    return s;
  }

  // *********** matrix operators *******************

  // comparison functions
  inline bool operator==(const Matrix3D &lhs, const Matrix3D &rhs)
  {
    return (memcmp(lhs.v_, rhs.v_, 9*sizeof(double)) == 0);
  }
  inline bool operator!=(const Matrix3D &lhs, const Matrix3D &rhs)
  {
    return (memcmp(lhs.v_, rhs.v_, 9*sizeof(double)) != 0);
  }

  //! \return A+B
  inline const Matrix3D operator+(const Matrix3D &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp;
    size_t i;
    for (i=0; i<3; ++i) {
      tmp[i][0] = A[i][0] + B[i][0];
      tmp[i][1] = A[i][1] + B[i][1];
      tmp[i][2] = A[i][2] + B[i][2];
    }
    return tmp;
  }

  //! \return A-B
  inline const Matrix3D operator-(const Matrix3D &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp;
    size_t i;
    for (i=0; i<3; ++i) {
      tmp[i][0] = A[i][0] - B[i][0];
      tmp[i][1] = A[i][1] - B[i][1];
      tmp[i][2] = A[i][2] - B[i][2];
    }
    return tmp;
  }

    //! Multiplies entry by entry. This is NOT a matrix multiplication. 
  inline const Matrix3D mult_element(const Matrix3D &A, 
                               const Matrix3D &B)
  {
    Matrix3D tmp;
    size_t i;
    for (i=0; i<3; ++i) {
      tmp[i][0] = A[i][0] * B[i][0];
      tmp[i][1] = A[i][1] * B[i][1];
      tmp[i][2] = A[i][2] * B[i][2];
    }
    return tmp;
  }

  //! Return the square of the Frobenius norm of A, i.e. sum (diag (A' * A))
  inline double Frobenius_2(const Matrix3D &A)
  {
    double fro=0.;
    for (int i=0; i<3; ++i) {
        fro += A[0][i]*A[0][i] + A[1][i]*A[1][i] + A[2][i]*A[2][i] ;
    }
    return fro;
  }

  inline Matrix3D transpose(const Matrix3D &A)
  {
    Matrix3D S;
    size_t i;
    for (i=0; i<3; ++i) {
        S[size_t(0)][i] = A[i][0];
        S[size_t(1)][i] = A[i][1];
        S[size_t(2)][i] = A[i][2];
    }
    return S;
  }

  inline void Matrix3D::operator+=(const Matrix3D &rhs)
  {
      v_[0] += rhs.v_[0]; v_[1] += rhs.v_[1]; v_[2] += rhs.v_[2];
      v_[3] += rhs.v_[3]; v_[4] += rhs.v_[4]; v_[5] += rhs.v_[5];
      v_[6] += rhs.v_[6]; v_[7] += rhs.v_[7]; v_[8] += rhs.v_[8];
  }

  inline void Matrix3D::operator-=(const Matrix3D &rhs)
  {
      v_[0] -= rhs.v_[0]; v_[1] -= rhs.v_[1]; v_[2] -= rhs.v_[2];
      v_[3] -= rhs.v_[3]; v_[4] -= rhs.v_[4]; v_[5] -= rhs.v_[5];
      v_[6] -= rhs.v_[6]; v_[7] -= rhs.v_[7]; v_[8] -= rhs.v_[8];
  }

  //! multiplies each entry by the scalar s
  inline void Matrix3D::operator*=(const double &s)
  {
      v_[0] *= s; v_[1] *= s; v_[2] *= s;
      v_[3] *= s; v_[4] *= s; v_[5] *= s;
      v_[6] *= s; v_[7] *= s; v_[8] *= s;
  }

  //! \f$ + B^T  \f$
  inline Matrix3D Matrix3D::plus_transpose(const Matrix3D &B) const
  {
    Matrix3D tmp;

    tmp.v_[0] = v_[0] + B.v_[0];
    tmp.v_[1] = v_[1] + B.v_[3];
    tmp.v_[2] = v_[2] + B.v_[6];
    
    tmp.v_[3] = v_[3] + B.v_[1];
    tmp.v_[4] = v_[4] + B.v_[4];
    tmp.v_[5] = v_[5] + B.v_[7];
    
    tmp.v_[6] = v_[6] + B.v_[2];
    tmp.v_[7] = v_[7] + B.v_[5];
    tmp.v_[8] = v_[8] + B.v_[8];

    return tmp;
  }

  //! \f$ += B^T  \f$
  inline void Matrix3D::plus_transpose_equal(const Matrix3D &B)
  {
    v_[0] += B.v_[0];
    v_[1] += B.v_[3];
    v_[2] += B.v_[6];
    
    v_[3] += B.v_[1];
    v_[4] += B.v_[4];
    v_[5] += B.v_[7];
    
    v_[6] += B.v_[2];
    v_[7] += B.v_[5];
    v_[8] += B.v_[8];
  }

  //! Computes \f$ A = v_1 v_2^T \f$
  inline Matrix3D& Matrix3D::outer_product(const Vector3D  &v1, const Vector3D &v2)
  {
    // remember, matrix entries are v_[0] to v_[8].
    
    // diagonal
    v_[0] = v1[0]*v2[0];
    v_[4] = v1[1]*v2[1];
    v_[8] = v1[2]*v2[2];

    // upper triangular part
    v_[1] = v1[0]*v2[1];
    v_[2] = v1[0]*v2[2];
    v_[5] = v1[1]*v2[2];

    // lower triangular part
    v_[3] = v2[0]*v1[1];
    v_[6] = v2[0]*v1[2];
    v_[7] = v2[1]*v1[2];

    return *this;
  }

  inline void Matrix3D::fill_lower_triangle()
  {
    v_[3] = v_[1];
    v_[6] = v_[2];
    v_[7] = v_[5];
  } 

  //! \return A*B
  inline const Matrix3D operator*(const Matrix3D  &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp;
    double sum;
    for (size_t i=0; i<3; ++i)
      for (size_t k=0; k<3; ++k)
        {
          sum = 0;
          for (size_t j=0; j<3; j++)
            sum = sum +  A[i][j] * B[j][k];
          tmp[i][k] = sum; 
        }
    return tmp;
  }
   
   //! multiplies each entry by the scalar s
  inline const Matrix3D Matrix3D::operator*(const double &s) const
  {
    Matrix3D temp;
    temp[0][0]=v_[0] * s; temp[0][1]=v_[1] * s; temp[0][2]=v_[2] * s;
    temp[1][0]=v_[3] * s; temp[1][1]=v_[4] * s; temp[1][2]=v_[5] * s;
    temp[2][0]=v_[6] * s; temp[2][1]=v_[7] * s; temp[2][2]=v_[8] * s;
    return temp;
  }
     //!friend function to allow for commutatative property of
     //! scalar mulitplication.
   inline const Matrix3D operator*(const double &s, const Matrix3D &A)
   {
     return (A.operator*(s));
   }
   
   
  //! \f$ C = A \times B \f$
  inline int matmult(Matrix3D& C, const Matrix3D  &A, 
                     const Matrix3D &B)
  {
    double sum;
    const double* row_i;
    const double* col_k;
    for (size_t i=0; i<3; ++i)
      for (size_t k=0; k<3; ++k)
        {
          row_i  = &(A[i][0]);
          col_k  = &(B[0][k]);
          sum = 0;
          for (size_t j=0; j<3; ++j)
            {
              sum  += *row_i * *col_k;
              row_i++;
              col_k += 3;
            }
          C[i][k] = sum; 
        }
    return 0;
  }

  /*! \brief Computes \f$ A v \f$ . */
  inline const Vector3D operator*(const Matrix3D  &A, const Vector3D &x)
  {
    Vector3D tmp; // initializes to 0
    for (size_t i=0; i<3; ++i)
      {
        const double* rowi = A[i];
        tmp[i] = rowi[0]*x[0] + rowi[1]*x[1] + rowi[2]*x[2];
      }
    return tmp;
  }

  /*! \brief Computes \f$ v^T A \f$ .
    
      This function implicitly considers the transpose of vector x times
      the matrix A and it is implicit that the returned vector must be
      transposed. */
  inline const Vector3D operator*(const Vector3D &x, const Matrix3D  &A)
  {
    Vector3D res(0., 0., 0.);
    for (size_t i=0; i<3; ++i)
      {
        const double* rowi = A[i];
        for (size_t j=0; j<3; ++j)
          res[j] += rowi[j] * x[i];
      }
    return res;
  }
   
  inline void eqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x)
  {
     v.mCoords[0] = A.v_[0]*x[0] + A.v_[1]*x.mCoords[1] + A.v_[2]*x.mCoords[2];
     v.mCoords[1] = A.v_[3]*x[0] + A.v_[4]*x.mCoords[1] + A.v_[5]*x.mCoords[2];
     v.mCoords[2] = A.v_[6]*x[0] + A.v_[7]*x.mCoords[1] + A.v_[8]*x.mCoords[2];
  }
   
  inline void plusEqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x)
  {
     v.mCoords[0] += A.v_[0]*x[0] + A.v_[1]*x.mCoords[1] + A.v_[2]*x.mCoords[2];
     v.mCoords[1] += A.v_[3]*x[0] + A.v_[4]*x.mCoords[1] + A.v_[5]*x.mCoords[2];
     v.mCoords[2] += A.v_[6]*x[0] + A.v_[7]*x.mCoords[1] + A.v_[8]*x.mCoords[2];
  }
   
  inline void plusEqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x)
  {
     v.mCoords[0] += A.v_[0]*x.mCoords[0] + A.v_[3]*x.mCoords[1] + A.v_[6]*x.mCoords[2];
     v.mCoords[1] += A.v_[1]*x.mCoords[0] + A.v_[4]*x.mCoords[1] + A.v_[7]*x.mCoords[2];
     v.mCoords[2] += A.v_[2]*x.mCoords[0] + A.v_[5]*x.mCoords[1] + A.v_[8]*x.mCoords[2];
  }
   
  inline void plusEqaA(Matrix3D& B, const double a, const Matrix3D &A) {
    B.v_[0] += a*A.v_[0]; B.v_[1] += a*A.v_[1]; B.v_[2] += a*A.v_[2]; 
    B.v_[3] += a*A.v_[3]; B.v_[4] += a*A.v_[4]; B.v_[5] += a*A.v_[5];
    B.v_[6] += a*A.v_[6]; B.v_[7] += a*A.v_[7]; B.v_[8] += a*A.v_[8];
  }

  inline double det(const Matrix3D &A) {
    return (  A.v_[0]*(A.v_[4]*A.v_[8]-A.v_[7]*A.v_[5])
            -A.v_[1]*(A.v_[3]*A.v_[8]-A.v_[6]*A.v_[5])
            +A.v_[2]*(A.v_[3]*A.v_[7]-A.v_[6]*A.v_[4]) );
  }

  inline void inv(Matrix3D &Ainv, const Matrix3D &A) {
    double inv_detA = 1 / (det(A));

    Ainv[0][0] = inv_detA*( A.v_[4]*A.v_[8]-A.v_[5]*A.v_[7] );
    Ainv[0][1] = inv_detA*( A.v_[2]*A.v_[7]-A.v_[8]*A.v_[1] );
    Ainv[0][2] = inv_detA*( A.v_[1]*A.v_[5]-A.v_[4]*A.v_[2] );

    Ainv[1][0] = inv_detA*( A.v_[5]*A.v_[6]-A.v_[8]*A.v_[3] );
    Ainv[1][1] = inv_detA*( A.v_[0]*A.v_[8]-A.v_[6]*A.v_[2] );
    Ainv[1][2] = inv_detA*( A.v_[2]*A.v_[3]-A.v_[5]*A.v_[0] );

    Ainv[2][0] = inv_detA*( A.v_[3]*A.v_[7]-A.v_[6]*A.v_[4] );
    Ainv[2][1] = inv_detA*( A.v_[1]*A.v_[6]-A.v_[7]*A.v_[0] );
    Ainv[2][2] = inv_detA*( A.v_[0]*A.v_[4]-A.v_[3]*A.v_[1] );
    return;
  }

  inline void timesInvA(Matrix3D& B, const Matrix3D &A) {

    Matrix3D Ainv;

    double inv_detA = 1 / ( det(A) );

    Ainv[0][0] = inv_detA*( A.v_[4]*A.v_[8]-A.v_[5]*A.v_[7] );
    Ainv[0][1] = inv_detA*( A.v_[2]*A.v_[7]-A.v_[8]*A.v_[1] );
    Ainv[0][2] = inv_detA*( A.v_[1]*A.v_[5]-A.v_[4]*A.v_[2] );

    Ainv[1][0] = inv_detA*( A.v_[5]*A.v_[6]-A.v_[8]*A.v_[3] );
    Ainv[1][1] = inv_detA*( A.v_[0]*A.v_[8]-A.v_[6]*A.v_[2] );
    Ainv[1][2] = inv_detA*( A.v_[2]*A.v_[3]-A.v_[5]*A.v_[0] );

    Ainv[2][0] = inv_detA*( A.v_[3]*A.v_[7]-A.v_[6]*A.v_[4] );
    Ainv[2][1] = inv_detA*( A.v_[1]*A.v_[6]-A.v_[7]*A.v_[0] );
    Ainv[2][2] = inv_detA*( A.v_[0]*A.v_[4]-A.v_[3]*A.v_[1] );

    B = B*Ainv;
  }

  inline void QR(Matrix3D &Q, Matrix3D &R, const Matrix3D &A) {
    // Compute the QR factorization of A.  This code uses the
    // Modified Gram-Schmidt method for computing the factorization.
    // The Householder version is more stable, but costs twice as many
    // floating point operations.

    Q = A;

    R[0][0] = sqrt(Q[0][0]*Q[0][0] + Q[1][0]*Q[1][0] + Q[2][0]*Q[2][0]);
    R[1][0] = 0.0L;
    R[2][0] = 0.0L;
    Q[0][0] /= R[0][0];
    Q[1][0] /= R[0][0];
    Q[2][0] /= R[0][0];

    R[0][1]  = Q[0][0]*Q[0][1] + Q[1][0]*Q[1][1] + Q[2][0]*Q[2][1];
    Q[0][1] -= Q[0][0]*R[0][1];
    Q[1][1] -= Q[1][0]*R[0][1];
    Q[2][1] -= Q[2][0]*R[0][1];

    R[0][2]  = Q[0][0]*Q[0][2] + Q[1][0]*Q[1][2] + Q[2][0]*Q[2][2];
    Q[0][2] -= Q[0][0]*R[0][2];
    Q[1][2] -= Q[1][0]*R[0][2];
    Q[2][2] -= Q[2][0]*R[0][2];

    R[1][1] = sqrt(Q[0][1]*Q[0][1] + Q[1][1]*Q[1][1] + Q[2][1]*Q[2][1]);
    R[2][1] = 0.0L;
    Q[0][1] /= R[1][1];
    Q[1][1] /= R[1][1];
    Q[2][1] /= R[1][1];

    R[1][2]  = Q[0][1]*Q[0][2] + Q[1][1]*Q[1][2] + Q[2][1]*Q[2][2];
    Q[0][2] -= Q[0][1]*R[1][2];
    Q[1][2] -= Q[1][1]*R[1][2];
    Q[2][2] -= Q[2][1]*R[1][2];
  
    R[2][2] = sqrt(Q[0][2]*Q[0][2] + Q[1][2]*Q[1][2] + Q[2][2]*Q[2][2]);
    Q[0][2] /= R[2][2];
    Q[1][2] /= R[2][2];
    Q[2][2] /= R[2][2];
    return;
  }

} // namespace Mesquite

#endif // Matrix3D_hpp
