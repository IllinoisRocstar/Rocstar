!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
subroutine PlaneNorm ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn, yn, &
  zn,same )
!
!*******************************************************************************
!
!! PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
!
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, REAL*8 :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of three points that constitute a line.  These points should not
!    be identical, nor colinear.
!
!    Output, REAL*8 :: XN, YN, ZN, the coordinates of the unit normal
!    vector to the plane containing the three points.
!
  implicit none
!
  REAL*8 :: temp
  REAL*8 :: x1
  REAL*8 :: x2
  REAL*8 :: x3
  REAL*8 :: xn
  REAL*8 :: y1
  REAL*8 :: y2
  REAL*8 :: y3
  REAL*8 :: yn
  REAL*8 :: z1
  REAL*8 :: z2
  REAL*8 :: z3
  REAL*8 :: zn
  integer :: same
  REAL*8 :: tol = 1.0e-10

  same = 0

!
!  The cross product (P2-P1) x (P3-P1) is a vector normal to
!  (P2-P1) and (P3-P1).
!
  call cross0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn, yn, zn )

  temp = sqrt ( xn * xn + yn * yn + zn * zn )

  if ( temp .LE. tol) then
     xn = 0.d0
     yn = 0.d0
     zn = 0.d0
     same = 1
  else
    xn = xn / temp
    yn = yn / temp
    zn = zn / temp
  end if

  return
end
subroutine cross0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! CROSS0_3D computes the cross product of (P1-P0) and (P2-P0) in 3D.
!
!
!  Discussion:
!
!    The vectors are specified with respect to a basis point P0.
!    We are computing the normal to the triangle (P0,P1,P2).
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, the coordinates of
!    three points.  The basis point is (X0,Y0,Z0).
!
!    Output, real X3, Y3, Z3, the cross product (P1-P0) x (P2-P0), or
!    (X1-X0,Y1-Y0,Z1-Z0) x (X2-X0,Y2-Y0,Z2-Z0).
!
  implicit none
!
  REAL*8 :: x0
  REAL*8 :: x1
  REAL*8 :: x2
  REAL*8 :: x3
  REAL*8 :: y0
  REAL*8 :: y1
  REAL*8 :: y2
  REAL*8 :: y3
  REAL*8 :: z0
  REAL*8 :: z1
  REAL*8 :: z2
  REAL*8 :: z3
!
  x3 = ( y1 - y0 ) * ( z2 - z0 ) - ( z1 - z0 ) * ( y2 - y0 )
  y3 = ( z1 - z0 ) * ( x2 - x0 ) - ( x1 - x0 ) * ( z2 - z0 )
  z3 = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )

  return
end

function CalcOffset ( x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! DOT_3D computes the dot product of a pair of vectors in 3D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real DOT_3D, the dot product.
!
  implicit none
!
  REAL*8 :: CalcOffset
  REAL*8 :: x1
  REAL*8 :: x2
  REAL*8 :: y1
  REAL*8 :: y2
  REAL*8 :: z1
  REAL*8 :: z2
!
  CalcOffset = x1 * x2 + y1 * y2 + z1 * z2

  return
end

function Offset ( vNorm,  x2, y2, z2 )
!
!*******************************************************************************
!
!! DOT_3D computes the dot product of a pair of vectors in 3D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real DOT_3D, the dot product.
!
  implicit none
!
  REAL*8 :: Offset
  real*8, dimension(1:3) :: Vnorm
  REAL*8 :: x2
  REAL*8 :: y2
  REAL*8 :: z2
!
  Offset = Vnorm(1) * x2 + Vnorm(2) * y2 + Vnorm(3) * z2

  return
end

