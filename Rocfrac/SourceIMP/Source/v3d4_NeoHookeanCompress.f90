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
SUBROUTINE V3D4_NeoHookeanCompress(coor,matcstet,lmcstet,R_in,d,ci, &
     S11,S22,S33,S12,S23,S13,numnp,nstart,nend,numcstet,numat_vol, &
     xmu,xlamda)

!________________________________________________________________________
!
!  V3D4 - Performs displacement based computations for Volumetric 3D 
!         4-node tetrahedra large deformation (i.e. nonlinear) elastic 
!         elements with linear interpolation functions. (constant strain 
!         tetrahedra). Returns the internal force vector R_in.
!
!  DATE: 04.2000                  AUTHOR: SCOT BREITENFELD
!
! Source:
!    "Nonlinear continuum mechanics for finite element analysis"
!     Javier Bonet and Richard D. Wood
!     ISBN 0-521-57272-X

!________________________________________________________________________

  IMPLICIT NONE
!---- Global variables
  INTEGER :: numnp          ! number of nodal points
  INTEGER :: numcstet       ! number of CSTet elements
  INTEGER :: numat_vol      ! number of volumetric materials
!--   coordinate array
  REAL*8, DIMENSION(1:3,1:numnp) :: coor
!--   elastic stiffness consts
  REAL*8, DIMENSION(1:9,1:numat_vol) :: ci
!--   internal force
  REAL*8, DIMENSION(1:3*numnp) :: R_in
!--  displacement vector
  REAL*8, DIMENSION(1:3*numnp) :: d
!--   CSTet stress
  REAL*8, DIMENSION(1:numcstet) :: S11, S22, S33, S12, S23, S13
!--  x, y and z displacements of nodes
  REAL*8 :: u1,u2,u3,u4,v1,v2,v3,v4,w1,w2,w3,w4
!--   coordinate holding variable
  REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  
  REAL*8 :: c11, c21, c31
!--   6*volume, inverse(6*volume),  and the volume      
  REAL*8 :: Vx6, Vx6inv, vol
!--   spacial derivatives
  REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12
!--   partial derivatives of the displacement 
  REAL*8 :: dudx,dvdy,dwdz,dudy,dvdx,dvdz,dwdy,dudz,dwdx
!--   strains
  REAL*8 :: E11,E22,E33,E12,E23,E13
! -- connectivity table for CSTet
  INTEGER, DIMENSION(1:4,1:numcstet) :: lmcstet
! -- mat number for CST element
  INTEGER, DIMENSION(1:numcstet) :: matcstet
! -- node numbers
  INTEGER :: n1,n2,n3,n4
! -- dummy and counters
  INTEGER :: i,j,nstart,nend
  INTEGER :: k1n1,k1n2,k1n3,k1n4,k2n1,k2n2,k2n3,k2n4
  INTEGER :: k3n1,k3n2,k3n3,k3n4
! -- 
  REAL*8 :: F11, F12, F13, F21, F22, F23, F31, F32, F33
  REAL*8 :: J1,J2,J2inv
  REAL*8 :: C11, C12, C13, C21, C22, C23, C31, C32, C33
  REAL*8 :: xmu, xlambda

  DO i = nstart, nend
     j = matcstet(i)
     
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
     u1 = d(k1n1)           ! (3*n1-2)
     u2 = d(k1n2)           ! (3*n2-2)
     u3 = d(k1n3)           ! (3*n3-2)
     u4 = d(k1n4)           ! (3*n4-2)
     v1 = d(k2n1)           ! (3*n1-1)
     v2 = d(k2n2)           ! (3*n2-1)
     v3 = d(k2n3)           ! (3*n3-1)
     v4 = d(k2n4)           ! (3*n4-1)
     w1 = d(k3n1)           ! (3*n1)
     w2 = d(k3n2)           ! (3*n2)
     w3 = d(k3n3)           ! (3*n3)
     w4 = d(k3n4)           ! (3*n4)
     
     x1 = coor(1,n1)
     x2 = coor(1,n2)
     x3 = coor(1,n3)
     x4 = coor(1,n4)
     y1 = coor(2,n1)
     y2 = coor(2,n2)
     y3 = coor(2,n3)
     y4 = coor(2,n4)
     z1 = coor(3,n1)
     z2 = coor(3,n2)
     z3 = coor(3,n3)
     z4 = coor(3,n4)
     
     
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

! See the maple worksheet 'V3D4.mws' for the derivation of Vx6

!         Vx6 = x2*y3*z4 - x2*y4*z3 - y2*x3*z4 + y2*x4*z3 + z2*x3*y4 
!     $        - z2*x4*y3 - x1*y3*z4 + x1*y4*z3 + x1*y2*z4 - x1*y2*z3
!     $        - x1*z2*y4 + x1*z2*y3 + y1*x3*z4 - y1*x4*z3 - y1*x2*z4 
!     $        + y1*x2*z3 + y1*z2*x4 - y1*z2*x3 - z1*x3*y4 + z1*x4*y3
!     $        + z1*x2*y4 - z1*x2*y3 - z1*y2*x4 + z1*y2*x3

     Vx6inv = 1.d0 / Vx6

! See the maple worksheet 'V3D4.mws' for the derivation of [B]

     B1  = -(y3*z4 - y4*z3 - y2*z4 + y2*z3 + z2*y4 - z2*y3) * Vx6inv
     B2  =  (x3*z4 - x4*z3 - x2*z4 + x2*z3 + z2*x4 - z2*x3) * Vx6inv
     B3  = -(x3*y4 - x4*y3 - x2*y4 + x2*y3 + y2*x4 - y2*x3) * Vx6inv
     B4  =  (y3*z4 - y4*z3 - y1*z4 + y1*z3 + z1*y4 - z1*y3) * Vx6inv
     B5  = -(x3*z4 - x4*z3 - x1*z4 + x1*z3 + z1*x4 - z1*x3) * Vx6inv
     B6  =  (x3*y4 - x4*y3 - x1*y4 + x1*y3 + y1*x4 - y1*x3) * Vx6inv
     B7  = -(y2*z4 - z2*y4 - y1*z4 + y1*z2 + z1*y4 - z1*y2) * Vx6inv
     B8  =  (x2*z4 - z2*x4 - x1*z4 + x1*z2 + z1*x4 - z1*x2) * Vx6inv
     B9  = -(x2*y4 - y2*x4 - x1*y4 + x1*y2 + y1*x4 - y1*x2) * Vx6inv
     B10 =  (y2*z3 - z2*y3 - y1*z3 + y1*z2 + z1*y3 - z1*y2) * Vx6inv
     B11 = -(x2*z3 - z2*x3 - x1*z3 + x1*z2 + z1*x3 - z1*x2) * Vx6inv
     B12 =  (x2*y3 - y2*x3 - x1*y3 + x1*y2 + y1*x3 - y1*x2) * Vx6inv
!
! Calculate displacement gradient (H)
!
     dudx = B1*u1 + B4*u2 + B7*u3 + B10*u4
     dvdy = B2*v1 + B5*v2 + B8*v3 + B11*v4
     dwdz = B3*w1 + B6*w2 + B9*w3 + B12*w4
     dudy = B2*u1 + B5*u2 + B8*u3 + B11*u4
     dvdx = B1*v1 + B4*v2 + B7*v3 + B10*v4
     dvdz = B3*v1 + B6*v2 + B9*v3 + B12*v4
     dwdy = B2*w1 + B5*w2 + B8*w3 + B11*w4
     dudz = B3*u1 + B6*u2 + B9*u3 + B12*u4
     dwdx = B1*w1 + B4*w2 + B7*w3 + B10*w4
! 
! deformation gradients F
!
     F11 = 1.d0 + dudx
     F12 = dudy
     F13 = dudz
     F21 = dvdx
     F22 = 1.d0 + dvdy
     F23 = dvdz
     F31 = dwdx
     F32 = dwdy
     F33 = 1.d0 + dwdz

     J1 = F11*(F22*F33-F23*F32)+F12*(F31*F23-F21*F33)+F13*(-F31*F22+F21*F32)
     IF(J1.LE.0.d0) THEN
        WRITE(*,*)'area has become zero for element',i
        STOP
     ENDIF
     J2 = J1*J1

!          T
!     C = F  F = right Cauchy-Green deformation tensor
!
!     Divided by third invariant squared

     J2inv = 1.d0/J2
     
     C11 = (F11*F11+F21*F21+F31*F31)*J2inv
     C12 = (F11*F12+F21*F22+F31*F32)*J2inv
     C13 = (F11*F13+F21*F23+F31*F33)*J2inv
     C21 = C12
     C22 = (F12*F12+F22*F22+F32*F32)*J2inv
     C23 = (F12*F13+F22*F23+F32*F33)*J2inv
     C31 = C13
     C32 = C23
     C33 = (F13*F13+F23*F23+F33*F33)*J2inv       

!
! Second Piola-Kirchoff tensor
! Eq. (5.28), pg. 124

     S11(i) = xmu*(1.d0-(C22*C33-C23*C32)) + xlambda*LOG(J1)*(C22*C33-C23*C32)
     S22(i) = xmu*(1.d0-(C11*C33-C31*C13)) + xlambda*LOG(J1)*(C11*C33-C31*C13)
     S33(i) = xmu*(1.d0-(C11*C22-C12*C21)) + xlambda*LOG(J1)*(C11*C22-C12*C21)
     S12(i) = (-C12*C33+C13*C32)*(xmu - xlambda*LOG(J1))
     S23(i) = (-C11*C23+C13*C21)*(xmu - xlambda*LOG(J1))
     S13(i) = (C12*C23-C13*C22)*(xmu - xlambda*LOG(J1))
         

! calculate the volume

     vol = Vx6 / 6.d0

! ASSEMBLE THE INTERNAL FORCE VECTOR
!
! local node 1
     R_in(k1n1) = R_in(k1n1) - vol* &
          ( S11(i)*B1*(1.d0+dudx) + S22(i)*B2*dudy + S33(i)*B3*dudz &
          + S12(i)*( B2*(1.d0+dudx) + B1*dudy ) &
          + S23(i)*( B3*dudy + B2*dudz ) &
          + S13(i)*( B3*(1.d0+dudx) + B1*dudz ) )
     R_in(k2n1) = R_in(k2n1) - vol* &
          ( S11(i)*B1*dvdx + S22(i)*B2*(1.d0+dvdy) + S33(i)*B3*dvdz &
          + S12(i)*( B1*(1.d0+dvdy) + B2*dvdx ) &
          + S23(i)*( B3*(1.d0+dvdy) + B2*dvdz ) &
          + S13(i)*( B3*dvdx + B1*dvdz ) )
     R_in(k3n1)   = R_in(k3n1)   - vol* &
          ( S11(i)*B1*dwdx + S22(i)*B2*dwdy + S33(i)*B3*(1.d0+dwdz) &
          + S12(i)*( B2*dwdx + B1*dwdy ) &
          + S23(i)*( B3*dwdy + B2*(1.d0 + dwdz) ) &
          + S13(i)*( B3*dwdx + B1*(1.d0 + dwdz) ) )
! local node 2 
     R_in(k1n2) = R_in(k1n2) - vol* &
          ( S11(i)*B4*(1.d0+dudx) + S22(i)*B5*dudy + S33(i)*B6*dudz &
          + S12(i)*( B5*(1.d0+dudx) + B4*dudy ) &
          + S23(i)*( B6*dudy + B5*dudz ) &
          + S13(i)*( B6*(1.d0+dudx) + B4*dudz ) )
     R_in(k2n2) = R_in(k2n2) - vol* &
          ( S11(i)*B4*dvdx + S22(i)*B5*(1.d0+dvdy) + S33(i)*B6*dvdz &
          + S12(i)*( B4*(1.d0+dvdy) + B5*dvdx ) &
          + S23(i)*( B6*(1.d0+dvdy) + B5*dvdz ) &
          + S13(i)*( B6*dvdx + B4*dvdz ) )
     R_in(k3n2)   = R_in(k3n2)   - vol* &
          ( S11(i)*B4*dwdx + S22(i)*B5*dwdy + S33(i)*B6*(1.d0+dwdz) &
          + S12(i)*( B5*dwdx + B4*dwdy ) &
          + S23(i)*( B6*dwdy + B5*(1.d0 + dwdz) ) &
          + S13(i)*( B6*dwdx + B4*(1.d0 + dwdz) ) )
! local node 3 
     R_in(k1n3) = R_in(k1n3) - vol* &
          ( S11(i)*B7*(1.d0+dudx) + S22(i)*B8*dudy + S33(i)*B9*dudz &
          + S12(i)*( B8*(1.d0+dudx) + B7*dudy ) &
          + S23(i)*( B9*dudy + B8*dudz ) &
          + S13(i)*( B9*(1.d0+dudx) + B7*dudz ) ) 
     R_in(k2n3) = R_in(k2n3) - vol* &
          ( S11(i)*B7*dvdx + S22(i)*B8*(1.d0+dvdy) + S33(i)*B9*dvdz &
          + S12(i)*( B7*(1.d0+dvdy) + B8*dvdx ) &
          + S23(i)*( B9*(1.d0+dvdy) + B8*dvdz ) &
          + S13(i)*( B9*dvdx + B7*dvdz ) )
     R_in(k3n3)   = R_in(k3n3)   - vol* &
          ( S11(i)*B7*dwdx + S22(i)*B8*dwdy + S33(i)*B9*(1.d0+dwdz) &
          + S12(i)*( B8*dwdx + B7*dwdy ) &
          + S23(i)*( B9*dwdy + B8*(1.d0 + dwdz) ) &
          + S13(i)*( B9*dwdx + B7*(1.d0 + dwdz) ) )
! local node 4         
     R_in(k1n4) = R_in(k1n4) - vol* &
          ( S11(i)*B10*(1.d0+dudx) + S22(i)*B11*dudy+S33(i)*B12*dudz &
          + S12(i)*( B11*(1.d0+dudx) + B10*dudy ) &
          + S23(i)*( B12*dudy + B11*dudz ) &
          + S13(i)*( B12*(1.d0+dudx) + B10*dudz ) )
     R_in(k2n4) = R_in(k2n4) - vol* &
          ( S11(i)*B10*dvdx + S22(i)*B11*(1.d0+dvdy)+S33(i)*B12*dvdz &
          + S12(i)*( B10*(1.d0+dvdy) + B11*dvdx ) &
          + S23(i)*( B12*(1.d0+dvdy) + B11*dvdz ) &
          + S13(i)*( B12*dvdx + B10*dvdz ) )
     R_in(k3n4)   = R_in(k3n4)   - vol* &
          ( S11(i)*B10*dwdx + S22(i)*B11*dwdy+S33(i)*B12*(1.d0+dwdz) &
          + S12(i)*( B11*dwdx + B10*dwdy ) &
          + S23(i)*( B12*dwdy + B11*(1.d0 + dwdz) ) &
          + S13(i)*( B12*dwdx + B10*(1.d0 + dwdz) ) )

          
  ENDDO
  
  RETURN
END SUBROUTINE V3D4_NeoHookeanCompress

