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

/*! \file TargetMatrix.hpp

Header file for the Mesquite::TargetMatrix class

  \author Thomas Leurent
  \date   2004-09-29
 */


#ifndef TargetMatrix_hpp
#define TargetMatrix_hpp

#include "Mesquite.hpp"
#include "Matrix3D.hpp"

namespace Mesquite
{
  
  /*! \class TargetMatrix
    \brief Class containing the target corner matrices for the context based smoothing. 
  */
  class TargetMatrix : public Matrix3D
  {
  public:

    TargetMatrix() : Matrix3D()
    {
      cK = 1.;
    }

    TargetMatrix(const Matrix3D &A) : Matrix3D(A)
    {
      cK = 1.;
    }

    TargetMatrix& operator=(const Matrix3D &A)
    {
      if (&A != this)
        Matrix3D::operator=(A);
      cK = 1.;
      return *this;
    }

    TargetMatrix& operator=(const TargetMatrix &A)
    {
      if (&A != this)
        Matrix3D::operator=(A);
      cK = A.cK;
      return *this;
    }

    double get_cK() const { return cK; }
    
  protected:
 
  private:
    double cK;
    
  };
  
} //namespace


#endif // TargetMatrix_hpp
