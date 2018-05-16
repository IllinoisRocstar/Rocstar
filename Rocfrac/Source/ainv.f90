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
SUBROUTINE ainv(ajac,ajacin,det,ndim)

  IMPLICIT NONE

!!****f* Rocfrac/Rocfrac/Source/ainv.f90
!!
!!  NAME
!!    ainv
!!
!!  FUNCTION
!!    Computes the det and inverse of a (3x3) matrix
!!
!!  USED BY
!!    shcalc, shcalc_3d10
!!
!!  INPUTS
!!    ndim -- size of input array (must be 3)
!!    ajac -- Input array (ndim x ndim)
!!
!!  OUTPUT
!!    ajacin -- inverse of ajac
!!    det -- determinate of ajac
!!
!!***

  INTEGER :: ndim
  REAL*8, DIMENSION(ndim,ndim) :: ajac, ajacin
  REAL*8 :: det, asmall

  asmall = 1.d-15

  det = ajac(1,1)*( ajac(2,2)*ajac(3,3) - ajac(2,3)* &
       ajac(3,2) ) &
       - ajac(1,2)*( ajac(2,1)*ajac(3,3) - ajac(2,3)* &
       ajac(3,1) ) &
       + ajac(1,3)*( ajac(2,1)*ajac(3,2) - ajac(2,2)* &
       ajac(3,1) )
      
  IF(dabs(det).LT.asmall) THEN
     PRINT*, 'cannot invert matrix in ainv'
     PRINT*,'ABS(det) =', dabs(det)
     STOP
  END IF


  det = 1.d0/det

!     construct inverse

  ajacin(1,1) = det*( ajac(2,2)*ajac(3,3) &
       - ajac(2,3)*ajac(3,2) )
  ajacin(1,2) = -det*( ajac(1,2)*ajac(3,3) &
       - ajac(1,3)*ajac(3,2) )
  ajacin(1,3) = det*( ajac(1,2)*ajac(2,3) &
       - ajac(1,3)*ajac(2,2) )
  ajacin(2,1) = -det*( ajac(2,1)*ajac(3,3) &
       - ajac(2,3)*ajac(3,1) )
  ajacin(2,2) = det*( ajac(1,1)*ajac(3,3) &
       - ajac(1,3)*ajac(3,1) )
  ajacin(2,3) = -det*( ajac(1,1)*ajac(2,3) &
       - ajac(1,3)*ajac(2,1) )
  ajacin(3,1) = det*( ajac(2,1)*ajac(3,2) &
       - ajac(2,2)*ajac(3,1) )
  ajacin(3,2) = -det*( ajac(1,1)*ajac(3,2) &
       - ajac(1,2)*ajac(3,1) )
  ajacin(3,3) = det*( ajac(1,1)*ajac(2,2) &
       - ajac(1,2)*ajac(2,1) )
  
  det = 1.d0/det
  
  RETURN
END SUBROUTINE ainv

