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
SUBROUTINE VOL_ELEM_MAT_ORTHO( ci_full, ci, NumMatVol, NumMatOrtho, MatOrtho, E11, E22, E33, &
             xnu12, xnu13, xnu23, G12, G13, G23, vx1o, vy1o, vz1o, vx2o, vy2o, vz2o, &
             vx3o, vy3o, vz3o )

  IMPLICIT NONE

  ! Input/output vars
  INTEGER :: NumMatVol, NumMatOrtho
  REAL*8, DIMENSION(1:6,1:6,1:NumMatVol) :: ci_full
  REAL*8, DIMENSION(1:9,1:NumMatVol) :: ci
  INTEGER, DIMENSION(1:NumMatVol) :: MatOrtho
  REAL*8, DIMENSION(1:NumMatOrtho) :: E11, E22, E33
  REAL*8, DIMENSION(1:NumMatOrtho) :: xnu12, xnu13, xnu23
  REAL*8, DIMENSION(1:NumMatOrtho) :: G12, G13, G23
  REAL*8, DIMENSION(1:NumMatOrtho) :: vx1o, vy1o, vz1o
  REAL*8, DIMENSION(1:NumMatOrtho) :: vx2o, vy2o, vz2o
  REAL*8, DIMENSION(1:NumMatOrtho) :: vx3o, vy3o, vz3o

  ! Internal vars
  INTEGER :: i, ii, j, k
  REAL*8 :: a, b
  REAL*8, DIMENSION(1:6,1:6) :: temp0, temp1, temp2, temp3, temp4
  REAL*8, DIMENSION(1:6,1:6) :: T, Ti, R, Ri
  REAL*8, DIMENSION(1:3) :: x, y, z, v1, v2, v3
  REAL*8 :: m1, m2, m3, n1, n2, n3, p1, p2, p3

  DO i = 1, NumMatVol

     IF ( MatOrtho(i) == 0 ) THEN

        ! Put the isotropic ci into the full ci
        ci_full(:,:,i) = 0.0
        ci_full(1,1,i) = ci(1,i)
        ci_full(1,2,i) = ci(2,i)
        ci_full(1,3,i) = ci(4,i)
        ci_full(2,1,i) = ci(2,i)
        ci_full(2,2,i) = ci(3,i)
        ci_full(2,3,i) = ci(5,i)
        ci_full(3,1,i) = ci(4,i)
        ci_full(3,2,i) = ci(5,i)
        ci_full(3,3,i) = ci(6,i)
        ci_full(4,4,i) = ci(7,i)
        ci_full(5,5,i) = ci(8,i)
        ci_full(6,6,i) = ci(9,i)

     ELSE

        ii = MatOrtho(i)

        ! Compute the unrotated material compliance matrix

        temp0(:,:) = 0.0
        temp1(:,:) = 0.0
        temp2(:,:) = 0.0

        temp0(1,1) = 1.0/E11(ii)
        temp0(1,2) = -xnu12(ii)/E11(ii)
        temp0(1,3) = -xnu13(ii)/E11(ii)
        temp0(2,1) = -xnu12(ii)/E11(ii)
        temp0(2,2) = 1.0/E22(ii)
        temp0(2,3) = -xnu23(ii)/E11(ii)
        temp0(3,1) = -xnu13(ii)/E11(ii)
        temp0(3,2) = -xnu23(ii)/E11(ii)
        temp0(3,3) = 1.0/E33(ii)
        temp0(4,4) = 1.0/G23(ii)
        temp0(5,5) = 1.0/G13(ii)
        temp0(6,6) = 1.0/G12(ii)

        
        ! Compute the stiffness matrix

        CALL invert6x6(temp0,temp1)


        ! Compute the rotation matrices

        R(:,:) = 0.0
        R(1,1) = 1.0
        R(2,2) = 1.0
        R(3,3) = 1.0
        R(4,4) = 2.0
        R(5,5) = 2.0
        R(6,6) = 2.0

        Ri(:,:) = 0.0
        Ri(1,1) = 1.0
        Ri(2,2) = 1.0
        Ri(3,3) = 1.0
        Ri(4,4) = 0.5
        Ri(5,5) = 0.5
        Ri(6,6) = 0.5

        T(:,:) = 0.0
        Ti(:,:) = 0.0
        
        v1(1) = vx1o(ii)
        v1(2) = vy1o(ii)
        v1(3) = vz1o(ii)

        v2(1) = vx2o(ii)
        v2(2) = vy2o(ii)
        v2(3) = vz2o(ii)

        v3(1) = vx3o(ii)
        v3(2) = vy3o(ii)
        v3(3) = vz3o(ii)

        x(1) = 1.0
        x(2) = 0.0
        x(3) = 0.0

        y(1) = 0.0
        y(2) = 1.0
        y(3) = 0.0

        z(1) = 0.0
        z(2) = 0.0
        z(3) = 1.0

        m1 = DOT_PRODUCT(x,v1)
        m2 = DOT_PRODUCT(x,v2)
        m3 = DOT_PRODUCT(x,v3)

        n1 = DOT_PRODUCT(y,v1)
        n2 = DOT_PRODUCT(y,v2)
        n3 = DOT_PRODUCT(y,v3)

        p1 = DOT_PRODUCT(z,v1)
        p2 = DOT_PRODUCT(z,v2)
        p3 = DOT_PRODUCT(z,v3)

        T(1,1) = m1*m1
        T(1,2) = n1*n1
        T(1,3) = p1*p1
        T(2,1) = m2*m2
        T(2,2) = n2*n2
        T(2,3) = p2*p2
        T(3,1) = m3*m3
        T(3,2) = n3*n3
        T(3,3) = p3*p3

        T(1,4) = 2.0*n1*p1
        T(1,5) = 2.0*p1*m1
        T(1,6) = 2.0*m1*n1
        T(2,4) = 2.0*n2*p2
        T(2,5) = 2.0*p2*m2
        T(2,6) = 2.0*m2*n2
        T(3,4) = 2.0*n3*p3
        T(3,5) = 2.0*p3*m3
        T(3,6) = 2.0*m3*n3

        T(4,1) = m2*m3
        T(4,2) = n2*n3
        T(4,3) = p2*p3
        T(5,1) = m3*m1
        T(5,2) = n3*n1
        T(5,3) = p3*p1
        T(6,1) = m1*m2
        T(6,2) = n1*n2
        T(6,3) = p1*p2

        T(4,4) = n2*p3 + n3*p2
        T(4,5) = p2*m3 + p3*m2
        T(4,6) = m2*n3 + m3*n2
        T(5,4) = n3*p1 + n1*p3
        T(5,5) = p3*m1 + p1*m3
        T(5,6) = m3*n1 + m1*n3
        T(6,4) = n1*p2 + n2*p1
        T(6,5) = p1*m2 + p2*m1
        T(6,6) = m1*n2 + m2*n1
        
        ! Invert the rotation matrix
        CALL invert6x6(T,Ti)
        
        ! Rotate the material stiffness matrix
        temp2 = MATMUL(Ti,temp1)
        temp3 = MATMUL(temp2, R)
        temp4 = MATMUL(temp3, T)
        ci_full(:,:,i) = MATMUL(temp4, Ri)

     END IF

  END DO        

END SUBROUTINE VOL_ELEM_MAT_ORTHO

