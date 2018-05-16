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
SUBROUTINE V3D4_NeoHookeanInCompressDef(coor,matcstet,lmcstet,R_in,d,ci, &
     S11,S22,S33,S12,S23,S13,numnp,nstart,nend,numcstet,numat_vol, &
     xmu,xkappa)
  
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
  REAL*8, DIMENSION(1:numat_vol) :: xmu, xkappa
  REAL*8, DIMENSION(1:4,1:4) :: eledb
  REAL*8,DIMENSION(1:3,1:3) :: xj
  REAL*8,DIMENSION(1:3,1:3) :: xjV0
  REAL*8 :: weight = 1.d0/6.d0
  INTEGER :: id, jd, in, ip
  REAL*8 :: xcoor
  REAL*8 :: detjb,detjbV0,evol,theta,press
  REAL*8 :: Ic, IIIc,Vx6old0,Vx6old
  REAL*8 :: onethird = 1.d0/3.d0
  REAL*8 :: val11, val21, val31
  
  REAL*8 :: V0x6, V0x6Inv
  REAL*8 :: Vx6, Vx6Inv
  REAL*8 :: vol, vol0
  REAL*8 :: term1, term2, Cinv
  REAL*8 :: SH1, SH2, SH3, SH4, SH5, SH6, SH7, SH8, SH9, SH10, SH11, SH12
!--  Coordinate subtractions
  REAL*8 :: x14d, x24d, x34d, y14d, y24d, y34d, z14d, z24d, z34d
!--  Coordinate subtractions: These are to speed up B calculation
  REAL*8 :: x12d, x13d, y12d, y13d, z12d, z13d
  REAL*8 :: sigma(3,3),F(3,3),FT(3,3),SPK2(3,3), btens(3,3)
  REAL*8 :: rmat3_det, detF,Jinv
  
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
     
     x1d = x1 + u1
     x2d = x2 + u2
     x3d = x3 + u3
     x4d = x4 + u4
     y1d = y1 + v1
     y2d = y2 + v2
     y3d = y3 + v3
     y4d = y4 + v4
     z1d = z1 + w1
     z2d = z2 + w2
     z3d = z3 + w3
     z4d = z4 + w4
!
!     creates the Jacobian matrix using equation ???
!     Using undeformed coordinates
     
     x14 = x1 - x4
     x24 = x2 - x4
     x34 = x3 - x4
     y14 = y1 - y4
     y24 = y2 - y4
     y34 = y3 - y4
     z14 = z1 - z4
     z24 = z2 - z4
     z34 = z3 - z4
     
     val11 =    y24*z34 - z24*y34
     val21 = -( x24*z34 - z24*x34 )
     val31 =    x24*y34 - y24*x34
     
     V0x6 = -( x14*val11 + y14*val21 + z14*val31 )
!
!     creates the Jacobian matrix using equation ???
!     Determines the determinant of xj and checks for zero values
!
!     Using deformed coordinates

     x12d = x1d - x2d          ! not used in vol. calc
     x13d = x1d - x3d          ! not used in vol. calc
     x14d = x1d - x4d
     x24d = x2d - x4d
     x34d = x3d - x4d
     y12d = y1d - y2d          ! not used in vol. calc
     y13d = y1d - y3d          ! not used in vol. calc
     y14d = y1d - y4d
     y24d = y2d - y4d
     y34d = y3d - y4d
     z12d = z1d - z2d          ! not used in vol. calc
     z13d = z1d - z3d          ! not used in vol. calc
     z14d = z1d - z4d
     z24d = z2d - z4d
     z34d = z3d - z4d
     
     val11 =    y24d*z34d - z24d*y34d
     val21 = -( x24d*z34d - z24d*x34d )
     val31 =    x24d*y34d - y24d*x34d
     
     Vx6 = -( x14d*val11 + y14d*val21 + z14d*val31 )
     
     IF(Vx6.LE.0.d0) THEN
        WRITE(*,100) i
        STOP
     ENDIF
     

! calculate the volume
      
     vol0 = V0x6/6.d0
     vol  = Vx6/6.d0
     
     V0x6Inv = 1.d0/V0x6
     Vx6inv = 1.d0/Vx6
     
! See the maple worksheet 'V3D4.mws' for the derivation of Vx6

! See the maple worksheet 'V3D4.mws' for the derivation of [B]
! NOTE: Factored for a more equivalent/compact form then maple's

     B1  = ( (y3-y4)*(z2-z4) - (y2-y4)*(z3-z4) ) * V0x6inv
     B2  = ( (z3-z4)*(x2-x4) - (z2-z4)*(x3-x4) ) * V0x6inv
     B3  = ( (x3-x4)*(y2-y4) - (x2-x4)*(y3-y4) ) * V0x6inv
     B4  = ( (y1-y3)*(z1-z4) - (y1-y4)*(z1-z3) ) * V0x6inv
     B5  = ( (z1-z3)*(x1-x4) - (z1-z4)*(x1-x3) ) * V0x6inv
     B6  = ( (x1-x3)*(y1-y4) - (x1-x4)*(y1-y3) ) * V0x6inv
     B7  = ( (y1-y4)*(z1-z2) - (y1-y2)*(z1-z4) ) * V0x6inv
     B8  = ( (z1-z4)*(x1-x2) - (z1-z2)*(x1-x4) ) * V0x6inv
     B9  = ( (x1-x4)*(y1-y2) - (x1-x2)*(y1-y4) ) * V0x6inv
     B10 = ( (y1-y2)*(z1-z3) - (y1-y3)*(z1-z2) ) * V0x6inv
     B11 = ( (z1-z2)*(x1-x3) - (z1-z3)*(x1-x2) ) * V0x6inv
     B12 = ( (x1-x2)*(y1-y3) - (x1-x3)*(y1-y2) ) * V0x6inv
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
     
     theta = Vol/Vol0
     
     press = xkappa(j)*(theta-1.d0)
     
     J1 = F11*(F22*F33-F23*F32)+F12*(F31*F23-F21*F33)+F13*(-F31*F22+F21*F32)
     IF(J1.LE.0.d0) THEN
        WRITE(*,*)'Volume has become zero for element',i
        STOP
     ENDIF
     JSqrd = J1*J1
     
! first invariant of C
     
     Ic = F11**2+F21**2+F31**2+F12**2+F22**2+F32**2+F13**2+F23**2+F33**2
     
! Third invariant of C
     
     IIIc = JSqrd        
     
!          T
!     C = F  F = right Cauchy-Green deformation tensor
!
     C11 = (F11*F11+F21*F21+F31*F31)
     C12 = (F11*F12+F21*F22+F31*F32)
     C13 = (F11*F13+F21*F23+F31*F33)
     C21 = C12 
     C22 = (F12*F12+F22*F22+F32*F32)
     C23 = (F12*F13+F22*F23+F32*F33) 
     C31 = C13
     C32 = C23
     C33 = (F13*F13+F23*F23+F33*F33)
!
!           T
!    b = F F  = left Cauchy-Green deformation tensor
!

     btens(1,1) = F11**2+F12**2+F13**2
     btens(1,2) = F11*F21+F12*F22+F13*F23
     btens(1,3) = F11*F31+F12*F32+F13*F33
     btens(2,1) = btens(1,2)
     btens(2,2) = F21**2+F22**2+F23**2
     btens(2,3) = F21*F31+F22*F32+F23*F33
     btens(3,1) = btens(1,3)
     btens(3,2) = btens(2,3)
     btens(3,3) = F31**2+F32**2+F33**2
!
!    Cauchy stress

     call NeoIncCauchyStress(3,xmu(j),J1,btens,sigma)
     call addpres(3,press,sigma)

!
! Second Piola-Kirchoff tensor
! Eq. (5.28), pg. 124 

     JSqrdInv = 1.d0/JSqrd
     
     term1 = xmu(j)*IIIc**(-onethird)
     term2 = press*J1
     
     Cinv = (C22*C33-C23*C32)*JSqrdInv
     S11(i) = term1*(1.d0-onethird*Ic*Cinv) + term2*Cinv
     Cinv = (C11*C33-C31*C13)*JSqrdInv
     S22(i) = term1*(1.d0-onethird*Ic*Cinv) + term2*Cinv
     Cinv = (C11*C22-C12*C21)*JSqrdInv
     S33(i) = term1*(1.d0-onethird*Ic*Cinv) + term2*Cinv
     Cinv = (-C12*C33+C13*C32)*JSqrdInv
     S12(i) = Cinv*(-term1*onethird*Ic + term2)
     Cinv = (-C11*C23+C13*C21)*JSqrdInv
     S23(i) = Cinv*(-term1*onethird*Ic + term2)
     Cinv = ( C12*C23-C13*C22)*JSqrdInv
     S13(i) = Cinv*(-term1*onethird*Ic + term2)

     SPK2(1,1) = S11(i)
     SPK2(1,2) = S12(i)
     SPK2(1,3) = S13(i)
     SPK2(2,1) = SPK2(1,2)
     SPK2(2,2) = S22(i)
     SPK2(2,3) = S23(i)
     SPK2(3,1) = SPK2(1,3)
     SPK2(3,2) = SPK2(2,3)
     SPK2(3,3) = S33(i)

     F(1,1) = F11
     F(1,2) = F12
     F(1,3) = F13
     F(2,1) = F21
     F(2,2) = F22
     F(2,3) = F23
     F(3,1) = F31
     F(3,2) = F32
     F(3,3) = F33

     FT(:,:) = TRANSPOSE(F(:,:))

     detF = rmat3_det(F)

     Jinv=1.d0/detF

!     print*,'lkjk',sigma(1,1),sigma(2,2), sigma(3,3)
!     SIGMA(:,:) = Jinv*F(:,:)*SPK2(:,:)*FT(:,:)
!     print*,sigma(1,1),sigma(2,2), sigma(3,3)

     SH1  = (y34d*z24d - y24d*z34d) * Vx6inv
     SH2  = (z34d*x24d - z24d*x34d) * Vx6inv
     SH3  = (x34d*y24d - x24d*y34d) * Vx6inv
     SH4  = (y13d*z14d - y14d*z13d) * Vx6inv
     SH5  = (z13d*x14d - z14d*x13d) * Vx6inv
     SH6  = (x13d*y14d - x14d*y13d) * Vx6inv
     SH7  = (y14d*z12d - y12d*z14d) * Vx6inv
     SH8  = (z14d*x12d - z12d*x14d) * Vx6inv
     SH9  = (x14d*y12d - x12d*y14d) * Vx6inv
     SH10 = (y12d*z13d - y13d*z12d) * Vx6inv
     SH11 = (z12d*x13d - z13d*x12d) * Vx6inv
     SH12 = (x12d*y13d - x13d*y12d) * Vx6inv

! ASSEMBLE THE INTERNAL FORCE VECTOR
!
! local node 1
     R_in(k1n1) = R_in(k1n1) - vol* &
          ( sigma(1,1)*SH1 + sigma(1,2)*SH2 + sigma(1,3)*SH3 )
     R_in(k2n1) = R_in(k2n1) - vol* &
          ( sigma(1,2)*SH1 + sigma(2,2)*SH2 + sigma(2,3)*SH3 )
     R_in(k3n1) = R_in(k3n1) - vol* &
          ( sigma(1,3)*SH1 + sigma(2,3)*SH2 + sigma(3,3)*SH3 )
! local node 2 
     R_in(k1n2) = R_in(k1n2) - vol* &
          ( sigma(1,1)*SH4 + sigma(1,2)*SH5 + sigma(1,3)*SH6 )
     R_in(k2n2) = R_in(k2n2) - vol* &
          ( sigma(1,2)*SH4 + sigma(2,2)*SH5 + sigma(2,3)*SH6 )
     R_in(k3n2) = R_in(k3n2) - vol* &
          ( sigma(1,3)*SH4 + sigma(2,3)*SH5 + sigma(3,3)*SH6 )
! local node 3 
     R_in(k1n3) = R_in(k1n3) - vol* &
          ( sigma(1,1)*SH7 + sigma(1,2)*SH8 + sigma(1,3)*SH9 )
     R_in(k2n3) = R_in(k2n3) - vol* &
          ( sigma(1,2)*SH7 + sigma(2,2)*SH8 + sigma(2,3)*SH9 )
     R_in(k3n3) = R_in(k3n3) - vol* &
          ( sigma(1,3)*SH7 + sigma(2,3)*SH8 + sigma(3,3)*SH9 )
! local node 4  
     R_in(k1n4) = R_in(k1n4) - vol* &
          ( sigma(1,1)*SH10 + sigma(1,2)*SH11 + sigma(1,3)*SH12 )
     R_in(k2n4) = R_in(k2n4) - vol* &
          ( sigma(1,2)*SH10 + sigma(2,2)*SH11 + sigma(2,3)*SH12 )
     R_in(k3n4) = R_in(k3n4) - vol* &
          ( sigma(1,3)*SH10 + sigma(2,3)*SH11 + sigma(3,3)*SH12 )

     
  ENDDO
  
  RETURN
      
100 FORMAT(' Negative Jacobian for element: ',i10)
END SUBROUTINE V3D4_NeoHookeanInCompressDef


FUNCTION rmat3_det ( a )
!
!*******************************************************************************
!
!! RMAT3_DET computes the determinant of a 3 by 3 matrix.
!
!
!  The determinant is the volume of the shape spanned by the vectors
!  making up the rows or columns of the matrix.
!
!  Formula:
!
!    det = a11 * a22 * a33 - a11 * a23 * a32
!        + a12 * a23 * a31 - a12 * a21 * a33
!        + a13 * a21 * a32 - a13 * a22 * a31
!
!  Parameters:
!
!    Input, real A(3,3), the matrix whose determinant is desired.
!
!    Output, real RMAT3_DET, the determinant of the matrix.
!
        IMPLICIT NONE
!
        REAL*8 a(3,3)
        REAL*8 rmat3_det
!
        rmat3_det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
             + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
             + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

        RETURN
      END FUNCTION rmat3_det

      SUBROUTINE NeoIncCauchyStress(ndime,xmu,detf,btens,sigma)
!
!-----------------------------------------------------------------------
!
!     Determines the Cauchy stress tensor for nearly incompressible 
!     neo-Hookean materials
!     
!
!     ndime  -->  number of dimensions
!     xmu    -->  mu coefficient
!     detf   -->  determinant of F
!     btens  -->  right Cauchy-Green tensor
!     sigma  -->  Cauchy stress tensor
!
!-----------------------------------------------------------------------
!
        IMPLICIT NONE
        
        REAL*8, DIMENSION(1:3,1:3) :: sigma, btens
        REAL*8 :: detf, xmu
        INTEGER :: ndime
        
        REAL*8 :: trace, xm
        INTEGER :: id, jd
!
!     Obtains the Cauchy stress using equation ???
!     First finds the trace
!
      trace = 0.d0
      DO id = 1, 3
         trace = trace + btens(id,id)
      ENDDO

      xm = xmu/(detf**(5./3.))
      DO id = 1, ndime      
         DO jd = 1, ndime
            sigma(id,jd) = xm*btens(id,jd)
         ENDDO
         sigma(id,id)=sigma(id,id)-xm*trace/3.
      ENDDO

      RETURN
    END SUBROUTINE NeoIncCauchyStress

    subroutine addpres(ndime,press,sigma)

!-----------------------------------------------------------------------
!
!     Adds the pressure to the deviatoric stresses
!
!     ndime  -->  number of dimensions
!     press  -->  pressure
!     sigma  -->  Cauchy stress tensor
!
!-----------------------------------------------------------------------
!

      implicit none

      integer :: ndime
      real*8, DIMENSION(1:3,1:3) :: sigma
      real*8 :: press

      integer :: id
!
!     Obtains the Cauchy stress from deviatoric and pressure
!     components using equation ???
!
      do id = 1,ndime      
         sigma(id,id) = sigma(id,id) + press
      enddo
      return
      end

