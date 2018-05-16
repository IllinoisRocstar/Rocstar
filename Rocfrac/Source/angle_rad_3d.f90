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

FUNCTION angle_rad_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )

!!****f* Rocfrac/Source/angle_rad_3d.f90
!!
!!  NAME
!!     angle_rad_3d
!!
!!  FUNCTION
!!     returns the angle in radians between two rays in 3D.
!!
!!  NOTES
!!
!!    The routine always computes the SMALLER of the two angles between
!!    two rays.  Thus, if the rays make an (exterior) angle of 
!!    1.5 radians, the (interior) angle of 0.5 radians will be reported.
!!
!!  Formula:
!!
!!    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )

!!  INPUTS
!!    real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are three points
!!    which define the rays.  The rays are:
!!    ( X1-X2, Y1-Y2, Z1-Z2 ) and ( X3-X2, Y3-Y2, Z3-Z2 ).
!!  
!!  OUTPUT
!!    real ANGLE_RAD_3D, the angle between the two rays, in radians.
!!    This value will always be between 0 and PI.  If either ray has 
!!    zero length, then the angle is returned as zero.
!!
!!***

  REAL*8 angle_rad_3d
  REAL*8 dot
  REAL*8 dot0_3d
  REAL*8 enorm0_3d
  REAL*8 v1norm
  REAL*8 v2norm
  REAL*8 x1
  REAL*8 x2
  REAL*8 x3
  REAL*8 y1
  REAL*8 y2
  REAL*8 y3
  REAL*8 z1
  REAL*8 z2
  REAL*8 z3

! computes the dot product of (P1-P0) and (P2-P0) in 3D.

  dot =  ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 ) + &
       ( z1 - z2 ) * ( z3 - z2 )
!
!  computes the Euclidean norm of (Point1-Point0) in 3D.

  v1norm = SQRT ( ( x2 - x1 )**2 + ( y2 - y1 )**2 + ( z2 - z1 )**2 )
  v2norm = SQRT ( ( x2 - x3 )**2 + ( y2 - y3 )**2 + ( z2 - z3 )**2 )
  
  IF ( v1norm == 0.d0 .OR. v2norm == 0.d0 ) THEN
     angle_rad_3d = 0.d0
  ELSE
     angle_rad_3d = ACOS ( dot / ( v1norm * v2norm ) )
  END IF
  
  RETURN
END FUNCTION angle_rad_3d

