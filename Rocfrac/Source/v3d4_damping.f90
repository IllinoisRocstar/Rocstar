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
SUBROUTINE v3d4_damping(vhalf,R_in,&
     numnp,&
     numcstet,lmcstet,coor,&
     nstart,nend, KappaDamp)

! Consistant Damping matrix [C]

  IMPLICIT NONE

  INTEGER :: i
  INTEGER :: numnp,numcstet
  INTEGER :: nstart, nend

  REAL*8 ::  Inv60=1.d0/60.d0, Inv120=1.d0/120.d0

  REAL*8 :: KappaDamp

! Tet connectivity table
  INTEGER, DIMENSION(1:4,1:numcstet) :: lmcstet
! nodal velocities
  REAL*8, DIMENSION(1:3*numnp) :: vhalf
!--   node numbers
  INTEGER :: n1,n2,n3,n4
!--   x, y, and z velocity of nodes
  REAL*8 :: u1,u2,u3,u4,v1,v2,v3,v4,w1,w2,w3,w4
!--   coordinate holding variable
  REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!-- Dummy
  REAL*8 :: c11, c21, c31

! mesh coordinates
  REAL*8, DIMENSION(1:3,1:numnp) :: coor
!--   6*volume, inverse(6*volume),  and the volume      
  REAL*8 :: Vx6, Vx6inv, vol,multiplier

  INTEGER :: k1n1,k1n2,k1n3,k1n4,k2n1,k2n2,k2n3,k2n4
  INTEGER :: k3n1,k3n2,k3n3,k3n4

! Contibution to the Net Force Vector
  REAL*8, DIMENSION(1:3*numnp) :: R_in

  DO i = nstart, nend

     n1 = lmcstet(1,i)
     n2 = lmcstet(2,i)
     n3 = lmcstet(3,i)
     n4 = lmcstet(4,i)
     
     k3n1 = 3*n1
     k3n2 = 3*n2
     k3n3 = 3*n3
     k3n4 = 3*n4
     
     k2n1 = k3n1 - 1
     k2n2 = k3n2 - 1
     k2n3 = k3n3 - 1
     k2n4 = k3n4 - 1
     
     k1n1 = k3n1 - 2
     k1n2 = k3n2 - 2
     k1n3 = k3n3 - 2
     k1n4 = k3n4 - 2

     ! k#n# dummy variables replaces:
     u1 = vhalf(k1n1)           ! (3*n1-2)
     u2 = vhalf(k1n2)           ! (3*n2-2)
     u3 = vhalf(k1n3)           ! (3*n3-2)
     u4 = vhalf(k1n4)           ! (3*n4-2)
     v1 = vhalf(k2n1)           ! (3*n1-1)
     v2 = vhalf(k2n2)           ! (3*n2-1)
     v3 = vhalf(k2n3)           ! (3*n3-1)
     v4 = vhalf(k2n4)           ! (3*n4-1)
     w1 = vhalf(k3n1)           ! (3*n1)
     w2 = vhalf(k3n2)           ! (3*n2)
     w3 = vhalf(k3n3)           ! (3*n3)
     w4 = vhalf(k3n4)           ! (3*n4)  

     x1 = coor(1,n1) ! Node 1, x-coor
     x2 = coor(1,n2) ! Node 2, x-coor
     x3 = coor(1,n3) ! Node 3, x-coor
     x4 = coor(1,n4) ! Node 4, x-coor
     y1 = coor(2,n1) ! Node 1, y-coor
     y2 = coor(2,n2) ! Node 2, y-coor
     y3 = coor(2,n3) ! Node 3, y-coor
     y4 = coor(2,n4) ! Node 4, y-coor
     z1 = coor(3,n1) ! Node 1, z-coor
     z2 = coor(3,n2) ! Node 2, z-coor
     z3 = coor(3,n3) ! Node 3, z-coor
     z4 = coor(3,n4) ! Node 4, z-coor
     
     x14 = x1 - x4
     x24 = x2 - x4
     x34 = x3 - x4
     y14 = y1 - y4
     y24 = y2 - y4
     y34 = y3 - y4
     z14 = z1 - z4
     z24 = z2 - z4
     z34 = z3 - z4
     c11 =    y24*z34 - z24*y34
     c21 = -( x24*z34 - z24*x34 )
     c31 =    x24*y34 - y24*x34
     
     Vx6 = -( x14*c11 + y14*c21 + z14*c31 )

     multiplier = KappaDamp*Vx6

! local node 1
     R_in(k1n1) = R_in(k1n1) - multiplier*(Inv60*u1+Inv120*(u2+u3+u4))
     R_in(k2n1) = R_in(k2n1) - multiplier*(Inv60*v1+Inv120*(v2+v3+v4))
     R_in(k3n1) = R_in(k3n1) - multiplier*(Inv60*w1+Inv120*(w2+w3+w4))
! local node 2 
     R_in(k1n2) = R_in(k1n2) - multiplier*(Inv60*u2+Inv120*(u1+u3+u4))
     R_in(k2n2) = R_in(k2n2) - multiplier*(Inv60*v2+Inv120*(v1+v3+v4))
     R_in(k3n2) = R_in(k3n2) - multiplier*(Inv60*w2+Inv120*(w1+w3+w4))
! local node 3 
     R_in(k1n3) = R_in(k1n3) - multiplier*(Inv60*u3+Inv120*(u1+u2+u4))
     R_in(k2n3) = R_in(k2n3) - multiplier*(Inv60*v3+Inv120*(v1+v2+v4))
     R_in(k3n3) = R_in(k3n3) - multiplier*(Inv60*w3+Inv120*(w1+w2+w4))
! local node 4
     R_in(k1n4) = R_in(k1n4) - multiplier*(Inv60*u4+Inv120*(u1+u2+u3))
     R_in(k2n4) = R_in(k2n4) - multiplier*(Inv60*v4+Inv120*(v1+v2+v3))
     R_in(k3n4) = R_in(k3n4) - multiplier*(Inv60*w4+Inv120*(w1+w2+w3))

  END DO
  
  RETURN
END SUBROUTINE v3d4_damping

