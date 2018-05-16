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
SUBROUTINE V3D10_NL_ARRUDA_BOYCE_DAMPING(coor,matcstet,lmcstet,R_in,d, &
     S11,S22,S33,S12,S23,S13, &
     numnp,nstart,nend,numcstet,numat_vol, &
     mu,kappa, rho, cd_fastest, DetFold, velo, dt)

!________________________________________________________________________
!
! V3D4 - Performs displacement based computations for Volumetric 3D, 
!        10-node Quadratic Linear Elastic Tetrahedron with quadratic 
!        interpolation functions. (linear strain tetrahedra). 
!        Large Deformation. Returns the internal force vector R_in.
!
!  DATE: 10.2001                  AUTHOR: SCOT BREITENFELD
!________________________________________________________________________

  IMPLICIT NONE
!-----Global variables
  INTEGER :: numnp          ! number of nodes
  INTEGER :: numat_vol      ! number of volumetric materials
  INTEGER :: numcstet       ! number of LSTets
!--   coordinate array
  REAL*8, DIMENSION(1:3,1:numnp) :: coor
!--   internal force
  REAL*8, DIMENSION(1:3*numnp) :: R_in
!--   displacement vector
  REAL*8, DIMENSION(1:3*numnp) :: d
!--   CSTet stress
  REAL*8, DIMENSION(1:4,1:numcstet) :: S11, S22, S33, S12, S23, S13
!--   connectivity table for CSTet  
  INTEGER, DIMENSION(1:10,1:numcstet) :: lmcstet
!--   mat number for CSTet element
  INTEGER, DIMENSION(1:numcstet) :: matcstet
!---- Local variables
!--   node numbers
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
!--   x, y, and z displacements of nodes  
  REAL*8 :: u1,u2,u3,u4,u5,u6,u7,u8,u9,u10
  REAL*8 :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
  REAL*8 :: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10
!--   6*volume and the volume      
  REAL*8 :: Vx6,Vx6inv
!--   spacial derivatives
  REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10
  REAL*8 :: B11,B12,B13,B14,B15,B16,B17,B18,B19,B20
  REAL*8 :: B21,B22,B23,B24,B25,B26,B27,B28,B29,B30
!--   partial derivatives of the displacement 
  REAL*8 :: dudx,dvdy,dwdz,dudy,dvdx,dvdz,dwdy,dudz,dwdx
!--   strains
  REAL*8 :: E11,E22,E33,E12,E23,E13
!--   coordinate holding variable
  REAL*8 :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
  REAL*8 :: y1,y2,y3,y4,y5,y6,y7,y8,y9,y10
  REAL*8 :: z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
!--   dummy and counters
  INTEGER :: i,j,nstart,nend
  REAL*8 :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8,aux9,aux10,aux11,aux12
!--   partial internal force
  REAL*8 :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18
  REAL*8 :: r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30
  REAL*8 :: g1, g2, g3, g4
  REAL*8 :: xN1, xN2, xN3, xN4

  INTEGER :: k1n1,k1n2,k1n3,k1n4,k1n5,k1n6,k1n7,k1n8,k1n9,k1n10
  INTEGER :: k2n1,k2n2,k2n3,k2n4,k2n5,k2n6,k2n7,k2n8,k2n9,k2n10
  INTEGER :: k3n1,k3n2,k3n3,k3n4,k3n5,k3n6,k3n7,k3n8,k3n9,k3n10
! -- 
  REAL*8 :: F11, F12, F13, F21, F22, F23, F31, F32, F33

!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  Coordinate subtractions: These are to speed up B calculation
  REAL*8 :: x12, x13, y12, y13, z12, z13
  REAL*8 :: val11, val21, val31
  REAL*8, DIMENSION(1:numat_vol) :: mu, kappa, rho
  real*8 :: cd_fastest
  real*8 , DIMENSION(1:4, 1:numcstet) :: DetFold
!--   velocity vector
  REAL*8, DIMENSION(1:3*numnp) :: velo
  real*8 :: dt

  REAL*8 :: StrssVisco(3,3)

!--   x, y, and z displacements of nodes  
  REAL*8 :: vel_u1,vel_u2,vel_u3,vel_u4,vel_u5,vel_u6,vel_u7,vel_u8,vel_u9,vel_u10
  REAL*8 :: vel_v1,vel_v2,vel_v3,vel_v4,vel_v5,vel_v6,vel_v7,vel_v8,vel_v9,vel_v10
  REAL*8 :: vel_w1,vel_w2,vel_w3,vel_w4,vel_w5,vel_w6,vel_w7,vel_w8,vel_w9,vel_w10

  REAL*8 :: F(3,3), Fdot(3,3)

  REAL*8,DIMENSION(1:4,1:4) :: GaussIntPt = RESHAPE( &
       (/0.58541020d0,0.13819660d0,0.13819660d0,0.13819660d0, &
       0.13819660d0,0.58541020d0,0.13819660d0,0.13819660d0, &
       0.13819660d0,0.13819660d0,0.58541020d0,0.13819660d0, &
       0.13819660d0,0.13819660d0,0.13819660d0,0.58541020d0/),(/4,4/) )
  
  INTEGER :: igpt

  REAL*8, DIMENSION(1:3,1:3) :: Fij, Cij
  INTEGER :: kk,ll,mm

  DO i = nstart, nend
     
     j   = matcstet(i)
     
     n1  = lmcstet(1,i)
     n2  = lmcstet(2,i)
     n3  = lmcstet(3,i)
     n4  = lmcstet(4,i)
     n5  = lmcstet(5,i)
     n6  = lmcstet(6,i)
     n7  = lmcstet(7,i)
     n8  = lmcstet(8,i)
     n9  = lmcstet(9,i)
     n10 = lmcstet(10,i)
     
     k3n1  = 3*n1
     k3n2  = 3*n2
     k3n3  = 3*n3
     k3n4  = 3*n4
     k3n5  = 3*n5
     k3n6  = 3*n6
     k3n7  = 3*n7
     k3n8  = 3*n8
     k3n9  = 3*n9
     k3n10 = 3*n10
     
     k2n1  = k3n1  - 1
     k2n2  = k3n2  - 1
     k2n3  = k3n3  - 1
     k2n4  = k3n4  - 1
     k2n5  = k3n5  - 1
     k2n6  = k3n6  - 1
     k2n7  = k3n7  - 1
     k2n8  = k3n8  - 1
     k2n9  = k3n9  - 1
     k2n10 = k3n10 - 1

     k1n1  = k3n1  - 2
     k1n2  = k3n2  - 2
     k1n3  = k3n3  - 2
     k1n4  = k3n4  - 2 
     k1n5  = k3n5  - 2
     k1n6  = k3n6  - 2
     k1n7  = k3n7  - 2
     k1n8  = k3n8  - 2 
     k1n9  = k3n9  - 2
     k1n10 = k3n10 - 2
                                ! k#n# dummy variables replaces:
     u1  = d(k1n1)          ! (3*n1 -2)
     u2  = d(k1n2)          ! (3*n2 -2)
     u3  = d(k1n3)          ! (3*n3 -2)
     u4  = d(k1n4)          ! (3*n4 -2)
     u5  = d(k1n5)          ! (3*n5 -2)
     u6  = d(k1n6)          ! (3*n6 -2)
     u7  = d(k1n7)          ! (3*n7 -2)
     u8  = d(k1n8)          ! (3*n8 -2)
     u9  = d(k1n9)          ! (3*n9 -2)
     u10 = d(k1n10)         ! (3*n10-2)         
     v1  = d(k2n1)          ! (3*n1 -1)
     v2  = d(k2n2)          ! (3*n2 -1)
     v3  = d(k2n3)          ! (3*n3 -1)
     v4  = d(k2n4)          ! (3*n4 -1)
     v5  = d(k2n5)          ! (3*n5 -1)
     v6  = d(k2n6)          ! (3*n6 -1)
     v7  = d(k2n7)          ! (3*n7 -1)
     v8  = d(k2n8)          ! (3*n8 -1)
     v9  = d(k2n9)          ! (3*n9 -1)
     v10 = d(k2n10)         ! (3*n10-1)
     w1  = d(k3n1)          ! (3*n1)
     w2  = d(k3n2)          ! (3*n2)
     w3  = d(k3n3)          ! (3*n3)
     w4  = d(k3n4)          ! (3*n4)
     w5  = d(k3n5)          ! (3*n5)
     w6  = d(k3n6)          ! (3*n6)
     w7  = d(k3n7)          ! (3*n7)
     w8  = d(k3n8)          ! (3*n8)
     w9  = d(k3n9)          ! (3*n9)
     w10 = d(k3n10)         ! (3*n10)

     vel_u1  = velo(k1n1)          ! (3*n1 -2)
     vel_u2  = velo(k1n2)          ! (3*n2 -2)
     vel_u3  = velo(k1n3)          ! (3*n3 -2)
     vel_u4  = velo(k1n4)          ! (3*n4 -2)
     vel_u5  = velo(k1n5)          ! (3*n5 -2)
     vel_u6  = velo(k1n6)          ! (3*n6 -2)
     vel_u7  = velo(k1n7)          ! (3*n7 -2)
     vel_u8  = velo(k1n8)          ! (3*n8 -2)
     vel_u9  = velo(k1n9)          ! (3*n9 -2)
     vel_u10 = velo(k1n10)         ! (3*n10-2)         
     vel_v1  = velo(k2n1)          ! (3*n1 -1)
     vel_v2  = velo(k2n2)          ! (3*n2 -1)
     vel_v3  = velo(k2n3)          ! (3*n3 -1)
     vel_v4  = velo(k2n4)          ! (3*n4 -1)
     vel_v5  = velo(k2n5)          ! (3*n5 -1)
     vel_v6  = velo(k2n6)          ! (3*n6 -1)
     vel_v7  = velo(k2n7)          ! (3*n7 -1)
     vel_v8  = velo(k2n8)          ! (3*n8 -1)
     vel_v9  = velo(k2n9)          ! (3*n9 -1)
     vel_v10 = velo(k2n10)         ! (3*n10-1)
     vel_w1  = velo(k3n1)          ! (3*n1)
     vel_w2  = velo(k3n2)          ! (3*n2)
     vel_w3  = velo(k3n3)          ! (3*n3)
     vel_w4  = velo(k3n4)          ! (3*n4)
     vel_w5  = velo(k3n5)          ! (3*n5)
     vel_w6  = velo(k3n6)          ! (3*n6)
     vel_w7  = velo(k3n7)          ! (3*n7)
     vel_w8  = velo(k3n8)          ! (3*n8)
     vel_w9  = velo(k3n9)          ! (3*n9)
     vel_w10 = velo(k3n10)         ! (3*n10)

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

     Vx6 = -( x14*val11 + y14*val21 + z14*val31 )

     Vx6inv = 1.d0 / Vx6

     aux1 = -(y3*z4 - y4*z3 - y2*z4 + y2*z3 + z2*y4 - z2*y3)
     aux2 =  (x3*z4 - x4*z3 - x2*z4 + x2*z3 + z2*x4 - z2*x3)
     aux3 = -(x3*y4 - x4*y3 - x2*y4 + x2*y3 + y2*x4 - y2*x3)
     aux4 =  (y3*z4 - y4*z3 - y1*z4 + y1*z3 + z1*y4 - z1*y3)
     aux5 = -(x3*z4 - x4*z3 - x1*z4 + x1*z3 + z1*x4 - z1*x3)
     aux6 =  (x3*y4 - x4*y3 - x1*y4 + x1*y3 + y1*x4 - y1*x3)
     aux7 = -(y2*z4 - z2*y4 - y1*z4 + y1*z2 + z1*y4 - z1*y2)
     aux8 =  (x2*z4 - z2*x4 - x1*z4 + x1*z2 + z1*x4 - z1*x2)
     aux9 = -(x2*y4 - y2*x4 - x1*y4 + x1*y2 + y1*x4 - y1*x2)
     aux10 = (y2*z3 - z2*y3 - y1*z3 + y1*z2 + z1*y3 - z1*y2)
     aux11 =-(x2*z3 - z2*x3 - x1*z3 + x1*z2 + z1*x3 - z1*x2)
     aux12 = (x2*y3 - y2*x3 - x1*y3 + x1*y2 + y1*x3 - y1*x2)

     r1 = 0.d0
     r2 = 0.d0
     r3 = 0.d0
     r4 = 0.d0
     r5 = 0.d0
     r6 = 0.d0
     r7 = 0.d0
     r8 = 0.d0
     r9 = 0.d0
     r10 = 0.d0
     r11 = 0.d0
     r12 = 0.d0
     r13 = 0.d0
     r14 = 0.d0
     r15 = 0.d0
     r16 = 0.d0
     r17 = 0.d0
     r18 = 0.d0
     r19 = 0.d0
     r20 = 0.d0
     r21 = 0.d0
     r22 = 0.d0
     r23 = 0.d0
     r24 = 0.d0
     r25 = 0.d0
     r26 = 0.d0
     r27 = 0.d0
     r28 = 0.d0
     r29 = 0.d0
     r30 = 0.d0

     DO igpt = 1, 4

        g1 = GaussIntPt(igpt,1)
        g2 = GaussIntPt(igpt,2)
        g3 = GaussIntPt(igpt,3)
        g4 = GaussIntPt(igpt,4)

     
        xN1 = (4.d0*g1-1.d0)   ! derivative of shape function
        xN2 = (4.d0*g2-1.d0)   ! dN_i/dzeta_i
        xN3 = (4.d0*g3-1.d0)
        xN4 = (4.d0*g4-1.d0)
!        xN5 = 4.d0*g1*g2
!        xN6 = 4.d0*g2*g3
!        xN7 = 4.d0*g3*g1
!        xN8 = 4.d0*g1*g4
!        xN9 = 4.d0*g2*g4
!        xN10= 4.d0*g3*g4

        B1 = aux1*xN1
        B2 = aux2*xN1
        B3 = aux3*xN1
        B4 = aux4*xN2
        B5 = aux5*xN2
        B6 = aux6*xN2
        B7 = aux7*xN3
        B8 = aux8*xN3
        B9 = aux9*xN3
        B10 =  aux10*xN4
        B11 =  aux11*xN4
        B12 =  aux12*xN4
        
        B13 =  4.d0*(g2*aux1 + g1*aux4)
        B14 =  4.d0*(g2*aux2 + g1*aux5)
        B15 =  4.d0*(g2*aux3 + g1*aux6)
        
        B16 =  4.d0*(g3*aux4 + g2*aux7)
        B17 =  4.d0*(g3*aux5 + g2*aux8)
        B18 =  4.d0*(g3*aux6 + g2*aux9)
        
        B19 =  4.d0*(g1*aux7 + g3*aux1)
        B20 =  4.d0*(g1*aux8 + g3*aux2)
        B21 =  4.d0*(g1*aux9 + g3*aux3)
        
        B22 =  4.d0*(g4*aux1 + g1*aux10)
        B23 =  4.d0*(g4*aux2 + g1*aux11)
        B24 =  4.d0*(g4*aux3 + g1*aux12)
        
        B25 =  4.d0*(g4*aux4 + g2*aux10)
        B26 =  4.d0*(g4*aux5 + g2*aux11)
        B27 =  4.d0*(g4*aux6 + g2*aux12)
        
        B28 =  4.d0*(g4*aux7 + g3*aux10)
        B29 =  4.d0*(g4*aux8 + g3*aux11)
        B30 =  4.d0*(g4*aux9 + g3*aux12)
        
 
!-----Calculate displacement gradient (H)
        dudx = (B1*u1 + B4*u2 + B7*u3 + B10*u4 + B13*u5 + B16*u6 + B19*u7 + B22*u8 + B25*u9 + B28*u10)*Vx6inv
        dvdy = (B2*v1 + B5*v2 + B8*v3 + B11*v4 + B14*v5 + B17*v6 + B20*v7 + B23*v8 + B26*v9 + B29*v10)*Vx6inv
        dwdz = (B3*w1 + B6*w2 + B9*w3 + B12*w4 + B15*w5 + B18*w6 + B21*w7 + B24*w8 + B27*w9 + B30*w10)*Vx6inv
        dudy = (B2*u1 + B5*u2 + B8*u3 + B11*u4 + B14*u5 + B17*u6 + B20*u7 + B23*u8 + B26*u9 + B29*u10)*Vx6inv
        dvdx = (B1*v1 + B4*v2 + B7*v3 + B10*v4 + B13*v5 + B16*v6 + B19*v7 + B22*v8 + B25*v9 + B28*v10)*Vx6inv
        dvdz = (B3*v1 + B6*v2 + B9*v3 + B12*v4 + B15*v5 + B18*v6 + B21*v7 + B24*v8 + B27*v9 + B30*v10)*Vx6inv
        dwdy = (B2*w1 + B5*w2 + B8*w3 + B11*w4 + B14*w5 + B17*w6 + B20*w7 + B23*w8 + B26*w9 + B29*w10)*Vx6inv
        dudz = (B3*u1 + B6*u2 + B9*u3 + B12*u4 + B15*u5 + B18*u6 + B21*u7 + B24*u8 + B27*u9 + B30*u10)*Vx6inv
        dwdx = (B1*w1 + B4*w2 + B7*w3 + B10*w4 + B13*w5 + B16*w6 + B19*w7 + B22*w8 + B25*w9 + B28*w10)*Vx6inv

! 
! deformation gradients F
!
        F11 = 1.d0 + ( dudx )
        F22 = 1.d0 + ( dvdy )
        F33 = 1.d0 + ( dwdz )
        F12 = dudy
        F21 = dvdx
        F23 = dvdz
        F32 = dwdy
        F13 = dudz
        F31 = dwdx

        F(1,1) = F11
        F(2,2) = F22
        F(3,3) = F33
        F(1,2) = F12
        F(1,3) = F13
        F(2,1) = F21
        F(2,3) = F23
        F(3,1) = F31
        F(3,2) = F32

        Fdot(1,1) = 1.d0 + &
             (B1*vel_u1 + B4*vel_u2 + B7*vel_u3 + B10*vel_u4 + B13*vel_u5 + &
             B16*vel_u6 + B19*vel_u7 + B22*vel_u8 + B25*vel_u9 + B28*vel_u10)*Vx6inv
        Fdot(2,2) = 1.d0 + &
             (B2*vel_v1 + B5*vel_v2 + B8*vel_v3 + B11*vel_v4 + B14*vel_v5 + &
             B17*vel_v6 + B20*vel_v7 + B23*vel_v8 + B26*vel_v9 + B29*vel_v10)*Vx6inv
        Fdot(3,3) = 1.d0 + &
             (B3*vel_w1 + B6*vel_w2 + B9*vel_w3 + B12*vel_w4 + B15*vel_w5 + &
             B18*vel_w6 + B21*vel_w7 + B24*vel_w8 + B27*vel_w9 + B30*vel_w10)*Vx6inv
        Fdot(1,2) = (B2*vel_u1 + B5*vel_u2 + B8*vel_u3 + B11*vel_u4 + B14*vel_u5 + &
             B17*vel_u6 + B20*vel_u7 + B23*vel_u8 + B26*vel_u9 + B29*vel_u10)*Vx6inv
        Fdot(2,1) = (B1*vel_v1 + B4*vel_v2 + B7*vel_v3 + B10*vel_v4 + B13*vel_v5 + &
             B16*vel_v6 + B19*vel_v7 + B22*vel_v8 + B25*vel_v9 + B28*vel_v10)*Vx6inv
        Fdot(2,3) = (B3*vel_v1 + B6*vel_v2 + B9*vel_v3 + B12*vel_v4 + B15*vel_v5 + &
             B18*vel_v6 + B21*vel_v7 + B24*vel_v8 + B27*vel_v9 + B30*vel_v10)*Vx6inv
        Fdot(3,2) = (B2*vel_w1 + B5*vel_w2 + B8*vel_w3 + B11*vel_w4 + B14*vel_w5 + &
             B17*vel_w6 + B20*vel_w7 + B23*vel_w8 + B26*vel_w9 + B29*vel_w10)*Vx6inv
        Fdot(1,3) = (B3*vel_u1 + B6*vel_u2 + B9*vel_u3 + B12*vel_u4 + B15*vel_u5 + &
             B18*vel_u6 + B21*vel_u7 + B24*vel_u8 + B27*vel_u9 + B30*vel_u10)*Vx6inv
        Fdot(3,1) = (B1*vel_w1 + B4*vel_w2 + B7*vel_w3 + B10*vel_w4 + B13*vel_w5 + &
             B16*vel_w6 + B19*vel_w7 + B22*vel_w8 + B25*vel_w9 + B28*vel_w10)*Vx6inv


        CALL ArtificialDamping(numcstet, numnp, &
             rho(j), cd_fastest, DetFold(igpt,i), dt, F, Fdot, Vx6, StrssVisco)

!
! Arruda-Boyce Nonlinear Elasticity Model 
!         
! -- NOTE: ci(7,j) : shear modulus

!--  (2) right Cauchy-Green deformation Tensor Cij = Fmi Fmj 
!             T
!        C = F  F
 
        DO kk=1,3
           DO ll=1,3
              Cij(kk,ll) = 0.d0
              DO mm=1,3
                 Cij(kk,ll)=Cij(kk,ll)+F(mm,kk)*F(mm,ll)
              ENDDO
           ENDDO
        ENDDO

        CALL ARRUDA_BOYCE(Cij, &
             S11(igpt,i),S22(igpt,i),S33(igpt,i),S12(igpt,i),S23(igpt,i),S13(igpt,i),i, &
             mu(j),kappa(j)) 

        s11(igpt,i) = s11(igpt,i) + StrssVisco(1,1)
        s22(igpt,i) = s22(igpt,i) + StrssVisco(2,2)
        s33(igpt,i) = s33(igpt,i) + StrssVisco(3,3)
        s12(igpt,i) = s12(igpt,i) + StrssVisco(1,2)
        s13(igpt,i) = s13(igpt,i) + StrssVisco(1,3)
        s23(igpt,i) = s23(igpt,i) + StrssVisco(2,3)

        r1 = r1 + &
             ( S11(igpt,i)*B1*(1.d0+dudx) + S22(igpt,i)*B2*dudy + S33(igpt,i)*B3*dudz &
             + S12(igpt,i)*( B2*(1.d0+dudx) + B1*dudy ) &
             + S23(igpt,i)*( B3*dudy + B2*dudz ) &
             + S13(igpt,i)*( B3*(1.d0+dudx) + B1*dudz ) )
        r2 = r2 +&
             ( S11(igpt,i)*B1*dvdx + S22(igpt,i)*B2*(1.d0+dvdy) + S33(igpt,i)*B3*dvdz &
             + S12(igpt,i)*( B1*(1.d0+dvdy) + B2*dvdx ) &
             + S23(igpt,i)*( B3*(1.d0+dvdy) + B2*dvdz ) &
             + S13(igpt,i)*( B3*dvdx + B1*dvdz ) )
        r3 = r3 +  &
             ( S11(igpt,i)*B1*dwdx + S22(igpt,i)*B2*dwdy + S33(igpt,i)*B3*(1.d0+dwdz) &
             + S12(igpt,i)*( B2*dwdx + B1*dwdy ) &
             + S23(igpt,i)*( B3*dwdy + B2*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B3*dwdx + B1*(1.d0 + dwdz) ) )
        
        r4 = r4 +  &
             ( S11(igpt,i)*B4*(1.d0+dudx) + S22(igpt,i)*B5*dudy + S33(igpt,i)*B6*dudz &
             + S12(igpt,i)*( B5*(1.d0+dudx) + B4*dudy ) &
             + S23(igpt,i)*( B6*dudy + B5*dudz ) &
             + S13(igpt,i)*( B6*(1.d0+dudx) + B4*dudz ) )
        r5 = r5 +   &
             ( S11(igpt,i)*B4*dvdx + S22(igpt,i)*B5*(1.d0+dvdy) + S33(igpt,i)*B6*dvdz &
             + S12(igpt,i)*( B4*(1.d0+dvdy) + B5*dvdx ) &
             + S23(igpt,i)*( B6*(1.d0+dvdy) + B5*dvdz ) &
             + S13(igpt,i)*( B6*dvdx + B4*dvdz ) )
        r6 = r6 +  &
             ( S11(igpt,i)*B4*dwdx + S22(igpt,i)*B5*dwdy + S33(igpt,i)*B6*(1.d0+dwdz) &
             + S12(igpt,i)*( B5*dwdx + B4*dwdy ) &
             + S23(igpt,i)*( B6*dwdy + B5*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B6*dwdx + B4*(1.d0 + dwdz) ) )
        
        r7 = r7 + &
             ( S11(igpt,i)*B7*(1.d0+dudx) + S22(igpt,i)*B8*dudy + S33(igpt,i)*B9*dudz &
             + S12(igpt,i)*( B8*(1.d0+dudx) + B7*dudy ) &
             + S23(igpt,i)*( B9*dudy + B8*dudz ) &
             + S13(igpt,i)*( B9*(1.d0+dudx) + B7*dudz ) )
        r8 = r8 + &
             ( S11(igpt,i)*B7*dvdx + S22(igpt,i)*B8*(1.d0+dvdy) + S33(igpt,i)*B9*dvdz &
             + S12(igpt,i)*( B7*(1.d0+dvdy) + B8*dvdx ) &
             + S23(igpt,i)*( B9*(1.d0+dvdy) + B8*dvdz ) &
             + S13(igpt,i)*( B9*dvdx + B7*dvdz ) )
        r9 = r9 + &
             ( S11(igpt,i)*B7*dwdx + S22(igpt,i)*B8*dwdy + S33(igpt,i)*B9*(1.d0+dwdz) &
             + S12(igpt,i)*( B8*dwdx + B7*dwdy ) &
             + S23(igpt,i)*( B9*dwdy + B8*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B9*dwdx + B7*(1.d0 + dwdz) ) )
        
        r10 = r10 + &
             ( S11(igpt,i)*B10*(1.d0+dudx) + S22(igpt,i)*B11*dudy+S33(igpt,i)*B12*dudz &
             + S12(igpt,i)*( B11*(1.d0+dudx) + B10*dudy ) &
             + S23(igpt,i)*( B12*dudy + B11*dudz ) &
             + S13(igpt,i)*( B12*(1.d0+dudx) + B10*dudz ) )
        r11 = r11 + & 
             ( S11(igpt,i)*B10*dvdx + S22(igpt,i)*B11*(1.d0+dvdy)+S33(igpt,i)*B12*dvdz &
             + S12(igpt,i)*( B10*(1.d0+dvdy) + B11*dvdx ) &
             + S23(igpt,i)*( B12*(1.d0+dvdy) + B11*dvdz ) &
             + S13(igpt,i)*( B12*dvdx + B10*dvdz ) )
        r12 = r12 +   &
             ( S11(igpt,i)*B10*dwdx + S22(igpt,i)*B11*dwdy+S33(igpt,i)*B12*(1.d0+dwdz) &
             + S12(igpt,i)*( B11*dwdx + B10*dwdy ) &
             + S23(igpt,i)*( B12*dwdy + B11*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B12*dwdx + B10*(1.d0 + dwdz) ) )
        
        r13 = r13 +   &
             ( S11(igpt,i)*B13*(1.d0+dudx) + S22(igpt,i)*B14*dudy + S33(igpt,i)*B15*dudz &
             + S12(igpt,i)*( B14*(1.d0+dudx) + B13*dudy ) &
             + S23(igpt,i)*( B15*dudy + B14*dudz ) &
             + S13(igpt,i)*( B15*(1.d0+dudx) + B13*dudz ) )
        r14 = r14 +   &
             ( S11(igpt,i)*B13*dvdx + S22(igpt,i)*B14*(1.d0+dvdy) + S33(igpt,i)*B15*dvdz &
             + S12(igpt,i)*( B13*(1.d0+dvdy) + B14*dvdx ) &
             + S23(igpt,i)*( B15*(1.d0+dvdy) + B14*dvdz ) &
             + S13(igpt,i)*( B15*dvdx + B13*dvdz ) )
        r15 = r15 +   &
             ( S11(igpt,i)*B13*dwdx + S22(igpt,i)*B14*dwdy + S33(igpt,i)*B15*(1.d0+dwdz) &
             + S12(igpt,i)*( B14*dwdx + B13*dwdy ) &
             + S23(igpt,i)*( B15*dwdy + B14*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B15*dwdx + B13*(1.d0 + dwdz) ) )
        
        r16 = r16 +  &
             ( S11(igpt,i)*B16*(1.d0+dudx) + S22(igpt,i)*B17*dudy + S33(igpt,i)*B18*dudz &
             + S12(igpt,i)*( B17*(1.d0+dudx) + B16*dudy ) &
             + S23(igpt,i)*( B18*dudy + B17*dudz ) &
             + S13(igpt,i)*( B18*(1.d0+dudx) + B16*dudz ) )
        r17 = r17 +  &
             ( S11(igpt,i)*B16*dvdx + S22(igpt,i)*B17*(1.d0+dvdy) + S33(igpt,i)*B18*dvdz &
             + S12(igpt,i)*( B16*(1.d0+dvdy) + B17*dvdx ) &
             + S23(igpt,i)*( B18*(1.d0+dvdy) + B17*dvdz ) &
             + S13(igpt,i)*( B18*dvdx + B16*dvdz ) )
        r18  = r18 +   &
             ( S11(igpt,i)*B16*dwdx + S22(igpt,i)*B17*dwdy + S33(igpt,i)*B18*(1.d0+dwdz) &
             + S12(igpt,i)*( B17*dwdx + B16*dwdy ) &
             + S23(igpt,i)*( B18*dwdy + B17*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B18*dwdx + B16*(1.d0 + dwdz) ) )
        
        r19 =  r19 +   &
             ( S11(igpt,i)*B19*(1.d0+dudx) + S22(igpt,i)*B20*dudy + S33(igpt,i)*B21*dudz &
             + S12(igpt,i)*( B20*(1.d0+dudx) + B19*dudy ) &
             + S23(igpt,i)*( B21*dudy + B20*dudz ) &
             + S13(igpt,i)*( B21*(1.d0+dudx) + B19*dudz ) )
        r20 =  r20 +  &
             ( S11(igpt,i)*B19*dvdx + S22(igpt,i)*B20*(1.d0+dvdy) + S33(igpt,i)*B21*dvdz &
             + S12(igpt,i)*( B19*(1.d0+dvdy) + B20*dvdx ) &
             + S23(igpt,i)*( B21*(1.d0+dvdy) + B20*dvdz ) &
             + S13(igpt,i)*( B21*dvdx + B19*dvdz ) )
        r21 =  r21 +  &
             ( S11(igpt,i)*B19*dwdx + S22(igpt,i)*B20*dwdy + S33(igpt,i)*B21*(1.d0+dwdz) &
             + S12(igpt,i)*( B20*dwdx + B19*dwdy ) &
             + S23(igpt,i)*( B21*dwdy + B20*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B21*dwdx + B19*(1.d0 + dwdz) ) )
        
        r22 = r22 +  &
             ( S11(igpt,i)*B22*(1.d0+dudx) + S22(igpt,i)*B23*dudy+S33(igpt,i)*B24*dudz &
             + S12(igpt,i)*( B23*(1.d0+dudx) + B22*dudy ) &
             + S23(igpt,i)*( B24*dudy + B23*dudz ) &
             + S13(igpt,i)*( B24*(1.d0+dudx) + B22*dudz ) )
        r23 = r23 +  &
             ( S11(igpt,i)*B22*dvdx + S22(igpt,i)*B23*(1.d0+dvdy)+S33(igpt,i)*B24*dvdz &
             + S12(igpt,i)*( B22*(1.d0+dvdy) + B23*dvdx ) &
             + S23(igpt,i)*( B24*(1.d0+dvdy) + B23*dvdz ) &
             + S13(igpt,i)*( B24*dvdx + B22*dvdz ) )
        r24 = r24 +  &
             ( S11(igpt,i)*B22*dwdx + S22(igpt,i)*B23*dwdy+S33(igpt,i)*B24*(1.d0+dwdz) &
             + S12(igpt,i)*( B23*dwdx + B22*dwdy ) &
             + S23(igpt,i)*( B24*dwdy + B23*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B24*dwdx + B22*(1.d0 + dwdz) ) )
        
        r25 = r25 +  &
             ( S11(igpt,i)*B25*(1.d0+dudx) + S22(igpt,i)*B26*dudy+S33(igpt,i)*B27*dudz &
             + S12(igpt,i)*( B26*(1.d0+dudx) + B25*dudy ) &
             + S23(igpt,i)*( B27*dudy + B26*dudz ) &
             + S13(igpt,i)*( B27*(1.d0+dudx) + B25*dudz ) )
        r26 = r26 +  &
             ( S11(igpt,i)*B25*dvdx + S22(igpt,i)*B26*(1.d0+dvdy)+S33(igpt,i)*B27*dvdz &
             + S12(igpt,i)*( B25*(1.d0+dvdy) + B26*dvdx ) &
             + S23(igpt,i)*( B27*(1.d0+dvdy) + B26*dvdz ) &
             + S13(igpt,i)*( B27*dvdx + B25*dvdz ) )
        r27 = r27 +  &
             ( S11(igpt,i)*B25*dwdx + S22(igpt,i)*B26*dwdy+S33(igpt,i)*B27*(1.d0+dwdz) &
             + S12(igpt,i)*( B26*dwdx + B25*dwdy ) &
             + S23(igpt,i)*( B27*dwdy + B26*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B27*dwdx + B25*(1.d0 + dwdz) ) )
        
        r28 = r28 +   &
             ( S11(igpt,i)*B28*(1.d0+dudx) + S22(igpt,i)*B29*dudy+S33(igpt,i)*B30*dudz &
             + S12(igpt,i)*( B29*(1.d0+dudx) + B28*dudy ) &
             + S23(igpt,i)*( B30*dudy + B29*dudz ) &
             + S13(igpt,i)*( B30*(1.d0+dudx) + B28*dudz ) )
        r29 = r29 +  &
             ( S11(igpt,i)*B28*dvdx + S22(igpt,i)*B29*(1.d0+dvdy)+S33(igpt,i)*B30*dvdz &
             + S12(igpt,i)*( B28*(1.d0+dvdy) + B29*dvdx ) &
             + S23(igpt,i)*( B30*(1.d0+dvdy) + B29*dvdz ) &
             + S13(igpt,i)*( B30*dvdx + B28*dvdz ) )
        r30 = r30 +  &
             ( S11(igpt,i)*B28*dwdx + S22(igpt,i)*B29*dwdy+S33(igpt,i)*B30*(1.d0+dwdz) &
             + S12(igpt,i)*( B29*dwdx + B28*dwdy ) &
             + S23(igpt,i)*( B30*dwdy + B29*(1.d0 + dwdz) ) &
             + S13(igpt,i)*( B30*dwdx + B28*(1.d0 + dwdz) ) )
     
     ENDDO
     
!Wi (i.e. weight) for 4 guass point integration is 1/4         

! Wi * 1/6 because the volume of a reference tetrahedra in 
! volume coordinates is 1/6

! ASSEMBLE THE INTERNAL FORCE VECTOR
!
! local node 1
     R_in(k1n1)  = R_in(k1n1)  - r1*0.04166666666666667d0
     R_in(k2n1)  = R_in(k2n1)  - r2*0.04166666666666667d0
     R_in(k3n1)  = R_in(k3n1)  - r3*0.04166666666666667d0
! local node 2
     R_in(k1n2)  = R_in(k1n2)  - r4*0.04166666666666667d0
     R_in(k2n2)  = R_in(k2n2)  - r5*0.04166666666666667d0
     R_in(k3n2)  = R_in(k3n2)  - r6*0.04166666666666667d0
! local node 3
     R_in(k1n3)  = R_in(k1n3)  - r7*0.04166666666666667d0
     R_in(k2n3)  = R_in(k2n3)  - r8*0.04166666666666667d0
     R_in(k3n3)  = R_in(k3n3)  - r9*0.04166666666666667d0
! local node 4
     R_in(k1n4)  = R_in(k1n4)  - r10*0.04166666666666667d0
     R_in(k2n4)  = R_in(k2n4)  - r11*0.04166666666666667d0
     R_in(k3n4)  = R_in(k3n4)  - r12*0.04166666666666667d0
! local node 5
     R_in(k1n5)  = R_in(k1n5)  - r13*0.04166666666666667d0
     R_in(k2n5)  = R_in(k2n5)  - r14*0.04166666666666667d0
     R_in(k3n5)  = R_in(k3n5)  - r15*0.04166666666666667d0
! local node 6
     R_in(k1n6)  = R_in(k1n6)  - r16*0.04166666666666667d0
     R_in(k2n6)  = R_in(k2n6)  - r17*0.04166666666666667d0
     R_in(k3n6)  = R_in(k3n6)  - r18*0.04166666666666667d0
! local node 7
     R_in(k1n7)  = R_in(k1n7)  - r19*0.04166666666666667d0
     R_in(k2n7)  = R_in(k2n7)  - r20*0.04166666666666667d0
     R_in(k3n7)  = R_in(k3n7)  - r21*0.04166666666666667d0
! local node 8
     R_in(k1n8)  = R_in(k1n8)  - r22*0.04166666666666667d0
     R_in(k2n8)  = R_in(k2n8)  - r23*0.04166666666666667d0
     R_in(k3n8)  = R_in(k3n8)  - r24*0.04166666666666667d0
! local node 9
     R_in(k1n9)  = R_in(k1n9)  - r25*0.04166666666666667d0
     R_in(k2n9)  = R_in(k2n9)  - r26*0.04166666666666667d0
     R_in(k3n9)  = R_in(k3n9)  - r27*0.04166666666666667d0
! local node 10
     R_in(k1n10) = R_in(k1n10) - r28*0.04166666666666667d0
     R_in(k2n10) = R_in(k2n10) - r29*0.04166666666666667d0
     R_in(k3n10) = R_in(k3n10) - r30*0.04166666666666667d0
     
  ENDDO
  RETURN
END SUBROUTINE V3D10_NL_ARRUDA_BOYCE_DAMPING

