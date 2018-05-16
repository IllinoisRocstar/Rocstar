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
SUBROUTINE V3D4_NeoHookeanInCompressPrin(coor,matcstet,lmcstet,R_in,d,ci, &
     S11,S22,S33,S12,S23,S13,numnp,nstart,nend,numcstet,numat_vol, &
     xmu,xkappa,xlambda)

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
!--   coordinate holding variable
  REAL*8 :: x1d,x2d,x3d,x4d,y1d,y2d,y3d,y4d,z1d,z2d,z3d,z4d
!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  Coordinate subtractions: These are to speed up B calculation
  REAL*8 :: x12, x13, y12, y13, z12, z13
!--   6*volume, inverse(6*volume),  and the volume      
!      REAL*8 :: Vx6, vol
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
  REAL*8 :: J1,JSqrd,JSqrdInv
  REAL*8 :: C11, C12, C13, C21, C22, C23, C31, C32, C33
  REAL*8 :: xmu, xlambda
  REAL*8, DIMENSION(1:4,1:4) :: eledb
  REAL*8,DIMENSION(1:3,1:3) :: xj
  REAL*8,DIMENSION(1:3,1:3) :: xjV0
  REAL*8 :: weight = 1.d0/6.d0
  INTEGER :: id, jd, in, ip
  REAL*8 :: xcoor
  REAL*8 :: detjb,detjbV0,evol,theta,xkappa,press
  REAL*8 :: Ic, IIIc,Vx6old0,Vx6old
  REAL*8 :: onethird = 1.d0/3.d0
  REAL*8 :: val11, val21, val31
  
  REAL*8 :: V0x6, V0x6Inv
  REAL*8 :: Vx6, Vx6Inv
  REAL*8 :: vol, vol0
  REAL*8 :: term1, term2, Cinv
  REAL*8, DIMENSION(1:3,1:3) :: btens,princ,sigma
  REAL*8, DIMENSION(1:3) :: stret,sprin
  REAL*8 :: detf
  REAL*8 :: SH1, SH2, SH3, SH4, SH5, SH6, SH7, SH8, SH9, SH10, SH11, SH12
  REAL*8 :: voldef,Vx6invdef,Vx6def

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
     
     x12 = x1 - x2          ! not used in vol. calc
     x13 = x1 - x3          ! not used in vol. calc
     x14 = x1 - x4
     x24 = x2 - x4
     x34 = x3 - x4
     y12 = y1 - y2          ! not used in vol. calc
     y13 = y1 - y3          ! not used in vol. calc
     y14 = y1 - y4
     y24 = y2 - y4
     y34 = y3 - y4
     z12 = z1 - z2          ! not used in vol. calc
     z13 = z1 - z3          ! not used in vol. calc
     z14 = z1 - z4
     z24 = z2 - z4
     z34 = z3 - z4
     
     val11 =    y24*z34 - z24*y34
     val21 = -( x24*z34 - z24*x34 )
     val31 =    x24*y34 - y24*x34
     
     V0x6 = -( x14*val11 + y14*val21 + z14*val31 )
         
     Vx6Inv = 1.d0/V0x6

! See the maple worksheet 'V3D4.mws' for the derivation of Vx6

! See the maple worksheet 'V3D4.mws' for the derivation of [B]
! NOTE: Factored for a more equivalent/compact form then maple's
     
     B1  = (y34*z24 - y24*z34) * Vx6inv
     B2  = (z34*x24 - z24*x34) * Vx6inv
     B3  = (x34*y24 - x24*y34) * Vx6inv
     B4  = (y13*z14 - y14*z13) * Vx6inv
     B5  = (z13*x14 - z14*x13) * Vx6inv
     B6  = (x13*y14 - x14*y13) * Vx6inv
     B7  = (y14*z12 - y12*z14) * Vx6inv
     B8  = (z14*x12 - z12*x14) * Vx6inv
     B9  = (x14*y12 - x12*y14) * Vx6inv
     B10 = (y12*z13 - y13*z12) * Vx6inv
     B11 = (z12*x13 - z13*x12) * Vx6inv
     B12 = (x12*y13 - x13*y12) * Vx6inv
! 
! deformation gradients F
!
     F11 = 1.d0 + ( B1*u1 + B4*u2 + B7*u3 + B10*u4 ) ! 1 + ( dudx )
     F22 = 1.d0 + ( B2*v1 + B5*v2 + B8*v3 + B11*v4 ) ! 1 + ( dvdy )
     F33 = 1.d0 + ( B3*w1 + B6*w2 + B9*w3 + B12*w4 ) ! 1 + ( dwdz )
     F12 = B2*u1 + B5*u2 + B8*u3 + B11*u4 ! dudy
     F21 = B1*v1 + B4*v2 + B7*v3 + B10*v4 ! dvdx
     F23 = B3*v1 + B6*v2 + B9*v3 + B12*v4 ! dvdz
     F32 = B2*w1 + B5*w2 + B8*w3 + B11*w4 ! dwdy
     F13 = B3*u1 + B6*u2 + B9*u3 + B12*u4 ! dudz
     F31 = B1*w1 + B4*w2 + B7*w3 + B10*w4 ! dwdx

     detf = F11*(F22*F33-F23*F32)+F12*(F31*F23-F21*F33)+F13*(-F31*F22+F21*F32)
     IF(detf.LE.0.d0) THEN
        WRITE(*,*)'Jacobian has become zero for element',i
        STOP
     ENDIF

!                                     
!     obtains the tensor b, left cauchy-green tensor using equation ???
!
     btens(1,1) = F11**2+F12**2+F13**2
     btens(1,2) = F21*F11+F12*F22+F13*F23
     btens(1,3) = F31*F11+F12*F32+F13*F33
     btens(2,1) = btens(1,2)
     btens(2,2) = F21**2+F22**2+F23**2
     btens(2,3) = F21*F31+F32*F22+F23*F33
     btens(3,1) = btens(1,3)
     btens(3,2) = btens(2,3)
     btens(3,3) = F31**2+F32**2+F33**2

     CALL jacobi(btens,stret,princ)
     
     CALL CauchyStressPrinc(3,xmu,xlambda,detf,stret,princ,sigma,sprin)

     S11(i) = sigma(1,1)
     S12(i) = sigma(1,2)
     S13(i) = sigma(1,3)
     S22(i) = sigma(2,2)
     S23(i) = sigma(2,3)
     S33(i) = sigma(3,3)

     x1 = x1 + u1
     x2 = x2 + u2
     x3 = x3 + u3
     x4 = x4 + u4
     y1 = y1 + v1
     y2 = y2 + v2
     y3 = y3 + v3
     y4 = y4 + v4
     z1 = z1 + w1
     z2 = z2 + w2
     z3 = z3 + w3
     z4 = z4 + w4
     
     x12 = x1 - x2          ! not used in vol. calc
     x13 = x1 - x3          ! not used in vol. calc
     x14 = x1 - x4
     x24 = x2 - x4
     x34 = x3 - x4
     y12 = y1 - y2          ! not used in vol. calc
     y13 = y1 - y3          ! not used in vol. calc
     y14 = y1 - y4
     y24 = y2 - y4
     y34 = y3 - y4
     z12 = z1 - z2          ! not used in vol. calc
     z13 = z1 - z3          ! not used in vol. calc
     z14 = z1 - z4
     z24 = z2 - z4
     z34 = z3 - z4
     
     val11 =    y24*z34 - z24*y34
     val21 = -( x24*z34 - z24*x34 )
     val31 =    x24*y34 - y24*x34
     
     Vx6def = -( x14*val11 + y14*val21 + z14*val31 )
     
     Vx6invdef = 1.d0 / Vx6def 

! calculate the volume
     voldef = Vx6def/6.d0

     IF(Vx6def.LE.0.d0) THEN
        WRITE(*,100) i
        STOP
     ENDIF
     
     SH1  = (y34*z24 - y24*z34) * Vx6invdef
     SH2  = (z34*x24 - z24*x34) * Vx6invdef
     SH3  = (x34*y24 - x24*y34) * Vx6invdef
     SH4  = (y13*z14 - y14*z13) * Vx6invdef
     SH5  = (z13*x14 - z14*x13) * Vx6invdef
     SH6  = (x13*y14 - x14*y13) * Vx6invdef
     SH7  = (y14*z12 - y12*z14) * Vx6invdef
     SH8  = (z14*x12 - z12*x14) * Vx6invdef
     SH9  = (x14*y12 - x12*y14) * Vx6invdef
     SH10 = (y12*z13 - y13*z12) * Vx6invdef
     SH11 = (z12*x13 - z13*x12) * Vx6invdef
     SH12 = (x12*y13 - x13*y12) * Vx6invdef

! ASSEMBLE THE INTERNAL FORCE VECTOR (uses Cauchy Stress)
!
! local node 1
     R_in(k1n1) = R_in(k1n1) - voldef* &
          ( S11(i)*SH1 + S12(i)*SH2 + S13(i)*SH3 )
     R_in(k2n1) = R_in(k2n1) - voldef* &
          ( S12(i)*SH1 + S22(i)*SH2 + S23(i)*SH3 )
     R_in(k3n1) = R_in(k3n1) - voldef* &
          ( S13(i)*SH1 + S23(i)*SH2 + S33(i)*SH3 )
! local node 2 
     R_in(k1n2) = R_in(k1n2) - voldef* &
          ( S11(i)*SH4 + S12(i)*SH5 + S13(i)*SH6 )
     R_in(k2n2) = R_in(k2n2) - voldef* &
          ( S12(i)*SH4 + S22(i)*SH5 + S23(i)*SH6 )
     R_in(k3n2) = R_in(k3n2) - voldef* &
          ( S13(i)*SH4 + S23(i)*SH5 + S33(i)*SH6 )
! local node 3 
     R_in(k1n3) = R_in(k1n3) - voldef* &
          ( S11(i)*SH7 + S12(i)*SH8 + S13(i)*SH9 )
     R_in(k2n3) = R_in(k2n3) - voldef* &
          ( S12(i)*SH7 + S22(i)*SH8 + S23(i)*SH9 )
     R_in(k3n3) = R_in(k3n3) - voldef* &
          ( S13(i)*SH7 + S23(i)*SH8 + S33(i)*SH9 )
! local node 4  
     R_in(k1n4) = R_in(k1n4) - voldef* &
          ( S11(i)*SH10 + S12(i)*SH11 + S13(i)*SH12 )
     R_in(k2n4) = R_in(k2n4) - voldef* &
          ( S12(i)*SH10 + S22(i)*SH11 + S23(i)*SH12 )
     R_in(k3n4) = R_in(k3n4) - voldef* &
          ( S13(i)*SH10 + S23(i)*SH11 + S33(i)*SH12 )

  ENDDO

  RETURN
      
100 FORMAT(' Negative Jacobian for element: ',i10)

END SUBROUTINE V3D4_NeoHookeanInCompressPrin

