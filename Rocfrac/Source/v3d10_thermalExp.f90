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
SUBROUTINE v3d10_thermalExp(coor,matcstet,lmcstet,R_in,ci, &
     S11l,S22l,S33l,&
     numnp,nstart,nend,numcstet,numat_vol, CoeffExp,Temperature,Temperature0)

!!****f* Rocfrac/Rocfrac/Source/v3d10_thermalExp.f90
!!
!!  NAME
!!    v3d10_thermalExp
!!
!!  FUNCTION
!!
!!    Computes the internal force for a small deformation, 
!!    linear elastic 10-node tetrahedral with thermal expansion
!!
!!  INPUTS
!!
!!   NumNP    -- Number of nodes
!!   numcstet -- Number of elements
!!   Coor     -- number of coordinates
!!   Matcstet  -- Material id
!!   d -- nodal displacement
!!   ci -- compliance matrix
!! 
!!       stored in the follow format:
!!
!!            [ ci(1) ci(2) ci(4)                   ]
!!            [       ci(3) ci(5)                   ]
!!            [             ci(6)                   ]
!!            [                   ci(7)             ]
!!            [                         ci(8)       ]
!!            [                               ci(9) ]
!!
!!   CoeffExp -- coeffecient of thermal expansion
!!   Temperature - nodal temperature
!!   Temperature0 - reference temperature
!!    
!!
!!   lmcstet  -- Nodal connectivity
!!   nstart, nend -- element beginning and end loop counter
!!   numat_vol -- number of materials
!!
!!  OUTPUT
!!
!!    Rin -- internal force vector
!!    S11l,S22l,S33l,S12l,S23l,S13l -- Stresses at gauss points
!!
!!****

!----- 4 Guass points 
  IMPLICIT NONE
!-----Global variables
  INTEGER :: numnp          ! number of nodes
  INTEGER :: numat_vol      ! number of volumetric materials
  INTEGER :: numcstet       ! number of LSTets
!--   coordinate array
  REAL*8, DIMENSION(1:3,1:numnp) :: coor
!--   elastic stiffness consts
  REAL*8, DIMENSION(1:9,1:numat_vol) :: ci

  REAL*8, DIMENSION(1:numat_vol) :: CoeffExp

!--   internal force
  REAL*8, DIMENSION(1:3*numnp) :: R_in
!--   CSTet stress
  REAL*8, DIMENSION(1:4,1:numcstet) :: S11l, S22l, S33l

  REAL*8, DIMENSION(1:4) :: StrssTh11, StrssTh22, StrssTh33


!--   connectivity table for CSTet  
  INTEGER, DIMENSION(1:10,1:numcstet) :: lmcstet
!--   mat number for CSTet element
  INTEGER, DIMENSION(1:numcstet) :: matcstet
!---- Local variables
!--   node numbers
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
!--   6*volume and inverse of 6*volume      
  REAL*8 :: Vx6, Vx6inv
!--   spacial derivatives
  REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10
  REAL*8 :: B11,B12,B13,B14,B15,B16,B17,B18,B19,B20
  REAL*8 :: B21,B22,B23,B24,B25,B26,B27,B28,B29,B30
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
!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  
  REAL*8 :: c11, c21, c31
  INTEGER :: k1n1,k1n2,k1n3,k1n4,k1n5,k1n6,k1n7,k1n8,k1n9,k1n10
  INTEGER :: k2n1,k2n2,k2n3,k2n4,k2n5,k2n6,k2n7,k2n8,k2n9,k2n10
  INTEGER :: k3n1,k3n2,k3n3,k3n4,k3n5,k3n6,k3n7,k3n8,k3n9,k3n10

  REAL*8, DIMENSION(1:NumNP) :: Temperature
  REAL*8 :: TemperatureGauss,Temperature0


  REAL*8,DIMENSION(1:4,1:4) :: GaussIntPt = RESHAPE( &
       (/0.58541020d0,0.13819660d0,0.13819660d0,0.13819660d0, &
       0.13819660d0,0.58541020d0,0.13819660d0,0.13819660d0, &
       0.13819660d0,0.13819660d0,0.58541020d0,0.13819660d0, &
       0.13819660d0,0.13819660d0,0.13819660d0,0.58541020d0/),(/4,4/) )
  
  INTEGER :: igpt

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

!        Vx6 = x2*y3*z4 - x2*y4*z3 - y2*x3*z4 + y2*x4*z3 + z2*x3*y4 
!     $        - z2*x4*y3 - x1*y3*z4 + x1*y4*z3 + x1*y2*z4 - x1*y2*z3
!     $        - x1*z2*y4 + x1*z2*y3 + y1*x3*z4 - y1*x4*z3 - y1*x2*z4 
!     $        + y1*x2*z3 + y1*z2*x4 - y1*z2*x3 - z1*x3*y4 + z1*x4*y3
!     $        + z1*x2*y4 - z1*x2*y3 - z1*y2*x4 + z1*y2*x3

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

     r1 = 0.
     r2 = 0.
     r3 = 0.
     r4 = 0.
     r5 = 0.
     r6 = 0.
     r7 = 0.
     r8 = 0.
     r9 = 0.
     r10 = 0.
     r11 = 0.
     r12 = 0.
     r13 = 0.
     r14 = 0.
     r15 = 0.
     r16 = 0.
     r17 = 0.
     r18 = 0.
     r19 = 0.
     r20 = 0.
     r21 = 0.
     r22 = 0.
     r23 = 0.
     r24 = 0.
     r25 = 0.
     r26 = 0.
     r27 = 0.
     r28 = 0.
     r29 = 0.
     r30 = 0.

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
        

! interpolate Temperature from nodes to Gauss Point

        TemperatureGauss = g1*(2.d0*g1-1.d0)*Temperature(n1) +                 &
             g2*(2.d0*g2-1.d0)*Temperature(n2) +                               &
             g3*(2.d0*g3-1.d0)*Temperature(n3) +                               &
             g4*(2.d0*g4-1.d0)*Temperature(n4) +                               &
             4.d0*g1*g2*Temperature(n5) + 4.d0*g2*g3*Temperature(n6) +         &
             4.d0*g3*g1*Temperature(n7) + 4.d0*g1*g4*Temperature(n8) +         &
             4.d0*g2*g4*Temperature(n9) + 4.d0*g3*g4*Temperature(n10)

        StrssTh11(igpt) = -ci(1,j)*CoeffExp(j)*( TemperatureGauss - Temperature0 ) !* Vx6inv
        StrssTh22(igpt) = -ci(3,j)*CoeffExp(j)*( TemperatureGauss - Temperature0 ) !* Vx6inv
        StrssTh33(igpt) = -ci(6,j)*CoeffExp(j)*( TemperatureGauss - Temperature0 ) !* Vx6inv


        r1 = r1 + StrssTh11(igpt)*B1
        r2 = r2 + StrssTh22(igpt)*B2
        r3 = r3 + StrssTh33(igpt)*B3
     
        r4 = r4 + StrssTh11(igpt)*B4
        r5 = r5 + StrssTh22(igpt)*B5
        r6 = r6 + StrssTh33(igpt)*B6
        
        r7 = r7 + StrssTh11(igpt)*B7 
        r8 = r8 + StrssTh22(igpt)*B8 
        r9 = r9 + StrssTh33(igpt)*B9 
        
        r10 = r10 + StrssTh11(igpt)*B10 
        r11 = r11 + StrssTh22(igpt)*B11 
        r12 = r12 + StrssTh33(igpt)*B12 
        
        r13 = r13 + StrssTh11(igpt)*B13 
        r14 = r14 + StrssTh22(igpt)*B14 
        r15 = r15 + StrssTh33(igpt)*B15 
        
        r16 = r16 + StrssTh11(igpt)*B16 
        r17 = r17 + StrssTh22(igpt)*B17 
        r18 = r18 + StrssTh33(igpt)*B18 
        
        r19 = r19 + StrssTh11(igpt)*B19 
        r20 = r20 + StrssTh22(igpt)*B20 
        r21 = r21 + StrssTh33(igpt)*B21 
        
        r22 = r22 + StrssTh11(igpt)*B22 
        r23 = r23 + StrssTh22(igpt)*B23 
        r24 = r24 + StrssTh33(igpt)*B24 
        
        r25 = r25 + StrssTh11(igpt)*B25 
        r26 = r26 + StrssTh22(igpt)*B26 
        r27 = r27 + StrssTh33(igpt)*B27 
        
        r28 = r28 + StrssTh11(igpt)*B28 
        r29 = r29 + StrssTh22(igpt)*B29 
        r30 = r30 + StrssTh33(igpt)*B30 

        S11l(igpt,i) = S11l(igpt,i) + StrssTh11(igpt)
        S22l(igpt,i) = S22l(igpt,i) + StrssTh22(igpt)
        S33l(igpt,i) = S33l(igpt,i) + StrssTh33(igpt)

     ENDDO

!     PRINT*,i, S11l(1:4,i)

! Wi (i.e. weight) for 4 guass integration is 1/4         

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
END SUBROUTINE v3d10_thermalExp




