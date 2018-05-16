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
SUBROUTINE v3d10_thermal(NumEl, NumNP, ElConnVol, Coor, Kappa,&
     Rnet, T, Rho,Cp, matcstet,numat_vol,MeshVel,ElemStart, ElemEnd)

  IMPLICIT NONE

!!****f* Rocfrac/Source/v3d10_thermal
!!
!!  NAME
!!    v3d10_thermal.f90
!!
!!  FUNCTION
!!    Calculates the internal force vector due to heat
!!    transfer for a 3D 10-node tet element.
!!
!!                        ALE
!!   {Rnet} = [K]{T} + [K]   {T}.
!!
!!               ALE
!!      where [K]    is the element matrices
!!            [K]    is the conductivity matrix
!!            {T}    is the temperature
!!
!!
!!  USED BY
!!     RocFracMain
!!
!!  INPUTS
!!     NumEl - Number of Elements
!!     NumNP - Number of Nodes
!!     ElConnVol - Element Connectivity Array
!!     Coor - Current Mesh coordinates
!!     Kappa - thermal conductivity
!!     Rnet - Accumulated "force" vector (RnetHT)
!!     T - Temperature
!!     RhoCp - Density * Cp
!!     MeshVel - Mesh motion velocity
!!
!!  OUTPUT
!!     
!!     Accumulated "force" vector with the addition of 
!!     the internal "force" vector (RnetHT) 
!!
!!***

  INTEGER :: NumNP, NumEl,numat_vol
  integer :: ElemStart, ElemEnd
  REAL*8, DIMENSION(1:numat_vol) :: Kappa
  REAL*8, DIMENSION(1:NumNP) :: T
  REAL*8, DIMENSION(1:NumNP) :: Rnet
  INTEGER , DIMENSION(1:10,1:NumEl) :: ElConnVol
  REAL*8, DIMENSION(1:3,1:NumNp) :: Coor
  REAL*8, DIMENSION(1:NumNp*3) :: MeshVel
  REAL*8, DIMENSION(1:numat_vol) :: Rho
  REAL*8, DIMENSION(1:numat_vol) :: Cp
  integer, dimension(1:NumEl) :: matcstet
  
  INTEGER :: i, imat, igpt,j
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10


  REAL*8, DIMENSION(1:10,1:10)  :: Kloc
  REAL*8, DIMENSION(1:10)  :: Tloc, Rinloc


!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  Coordinate subtractions: These are to speed up B calculation
  REAL*8 :: x12, x13, y12, y13, z12, z13

  REAL*8 :: c11, c21, c31
!--   6*volume and inverse of 6*volume      
  REAL*8 :: Vx6, Vx6inv, vol

!--   coordinate holding variable
  REAL*8 :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
  REAL*8 :: y1,y2,y3,y4,y5,y6,y7,y8,y9,y10
  REAL*8 :: z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
  REAL*8 :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8,aux9,aux10,aux11,aux12
  REAL*8 :: g1, g2, g3, g4
  REAL*8 :: xN1, xN2, xN3, xN4
  INTEGER :: k1n1,k1n2,k1n3,k1n4,k1n5,k1n6,k1n7,k1n8,k1n9,k1n10
  INTEGER :: k2n1,k2n2,k2n3,k2n4,k2n5,k2n6,k2n7,k2n8,k2n9,k2n10
  INTEGER :: k3n1,k3n2,k3n3,k3n4,k3n5,k3n6,k3n7,k3n8,k3n9,k3n10
!--   partial internal force
  REAL*8 :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18

  REAL*8,DIMENSION(1:4,1:4) :: GaussIntPt = RESHAPE( &
       (/0.58541020d0,0.13819660d0,0.13819660d0,0.13819660d0, &
         0.13819660d0,0.58541020d0,0.13819660d0,0.13819660d0, &
         0.13819660d0,0.13819660d0,0.58541020d0,0.13819660d0, &
         0.13819660d0,0.13819660d0,0.13819660d0,0.58541020d0/),(/4,4/) )

  real*8, dimension(1:3,1:10) :: B
  real*8, dimension(1:10,1:3) :: BT
  real*8, dimension(1:10) :: ShapFnct
  real*8, dimension(1:3) :: MeshVelNd 

  REAL*8, DIMENSION(1:3,1:3) :: KappaMatrx = RESHAPE( &
       (/1.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,1.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,1.000000000000000/),(/3,3/) )

  DO i = ElemStart, ElemEnd
     imat =  matcstet(i)
     
     n1  = ElConnVol(1,i)
     n2  = ElConnVol(2,i) 
     n3  = ElConnVol(3,i)
     n4  = ElConnVol(4,i)
     n5  = ElConnVol(5,i)
     n6  = ElConnVol(6,i) 
     n7  = ElConnVol(7,i)
     n8  = ElConnVol(8,i) 
     n9  = ElConnVol(9,i)
     n10 = ElConnVol(10,i)

    
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

     x1 = coor(1,n1)
     x2 = coor(1,n2)
     x3 = coor(1,n3)
     x4 = coor(1,n4) 
     x5 = coor(1,n5)
     x6 = coor(1,n6)
     x7 = coor(1,n7)
     x8 = coor(1,n8)
     x9  = coor(1,n9)
     x10 = coor(1,n10)

     y1 = coor(2,n1)
     y2 = coor(2,n2)
     y3 = coor(2,n3)
     y4 = coor(2,n4)
     y5 = coor(2,n5)
     y6 = coor(2,n6)
     y7 = coor(2,n7)
     y8 = coor(2,n8)
     y9 = coor(2,n9)
     y10 = coor(2,n10)

     z1 = coor(3,n1)
     z2 = coor(3,n2)
     z3 = coor(3,n3)
     z4 = coor(3,n4)
     z5 = coor(3,n5)
     z6 = coor(3,n6)
     z7 = coor(3,n7)
     z8 = coor(3,n8)
     z9 = coor(3,n9)
     z10 = coor(3,n10)

     x14 = x1 - x4
     x24 = x2 - x4
     x34 = x3 - x4
     y14 = y1 - y4
     y24 = y2 - y4
     y34 = y3 - y4
     z14 = z1 - z4
     z24 = z2 - z4
     z34 = z3 - z4

     x12 = x1 - x2          ! not used in vol. calc
     x13 = x1 - x3          ! not used in vol. calc
     y12 = y1 - y2          ! not used in vol. calc
     y13 = y1 - y3          ! not used in vol. calc
     z12 = z1 - z2          ! not used in vol. calc
     z13 = z1 - z3          ! not used in vol. calc

     c11 =    y24*z34 - z24*y34
     c21 = -( x24*z34 - z24*x34 )
     c31 =    x24*y34 - y24*x34

     Vx6 = -( x14*c11 + y14*c21 + z14*c31 )

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
  
     Tloc(1) = T(n1)
     Tloc(2) = T(n2)
     Tloc(3) = T(n3)
     Tloc(4) = T(n4)
     Tloc(5) = T(n5)
     Tloc(6) = T(n6)
     Tloc(7) = T(n7)
     Tloc(8) = T(n8)
     Tloc(9) = T(n9)
     Tloc(10) = T(n10)

     Kloc = 0.

     DO igpt = 1, 4

        g1 = GaussIntPt(igpt,1)
        g2 = GaussIntPt(igpt,2)
        g3 = GaussIntPt(igpt,3)
        g4 = GaussIntPt(igpt,4)

        ShapFnct(1) = g1*(2.*g1-1.)
        ShapFnct(2) = g2*(2.*g2-1.)
        ShapFnct(3) = g3*(2.*g3-1.)
        ShapFnct(4) = g4*(2.*g4-1.)
        ShapFnct(5) = 4.*g1*g2
        ShapFnct(6) = 4.*g2*g3
        ShapFnct(7) = 4.*g3*g1
        ShapFnct(8) = 4.*g1*g4
        ShapFnct(9) = 4.*g2*g4
        ShapFnct(10) = 4.*g3*g4


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

        B(1,1) = aux1*xN1
        B(2,1) = aux2*xN1
        B(3,1) = aux3*xN1
        B(1,2) = aux4*xN2
        B(2,2) = aux5*xN2
        B(3,2) = aux6*xN2
        B(1,3) = aux7*xN3
        B(2,3) = aux8*xN3
        B(3,3) = aux9*xN3
        B(1,4) =  aux10*xN4
        B(2,4) =  aux11*xN4
        B(3,4) =  aux12*xN4
        
        B(1,5) =  4.d0*(g2*aux1 + g1*aux4)
        B(2,5) =  4.d0*(g2*aux2 + g1*aux5)
        B(3,5) =  4.d0*(g2*aux3 + g1*aux6)
        
        B(1,6) =  4.d0*(g3*aux4 + g2*aux7)
        B(2,6) =  4.d0*(g3*aux5 + g2*aux8)
        B(3,6) =  4.d0*(g3*aux6 + g2*aux9)
        
        B(1,7) =  4.d0*(g1*aux7 + g3*aux1)
        B(2,7) =  4.d0*(g1*aux8 + g3*aux2)
        B(3,7) =  4.d0*(g1*aux9 + g3*aux3)
        
        B(1,8) =  4.d0*(g4*aux1 + g1*aux10)
        B(2,8) =  4.d0*(g4*aux2 + g1*aux11)
        B(3,8) =  4.d0*(g4*aux3 + g1*aux12)
        
        B(1,9) =  4.d0*(g4*aux4 + g2*aux10)
        B(2,9) =  4.d0*(g4*aux5 + g2*aux11)
        B(3,9) =  4.d0*(g4*aux6 + g2*aux12)
        
        B(1,10) =  4.d0*(g4*aux7 + g3*aux10)
        B(2,10) =  4.d0*(g4*aux8 + g3*aux11)
        B(3,10) =  4.d0*(g4*aux9 + g3*aux12)

        B = Vx6inv*B

!        BT = TRANSPOSE(B)
! Compute the local stiffness matrix

!        B = MATMUL(KappaMatrx,B)

        vol = Vx6 / 6.d0

!        print*,'***',Kloc(9,10)
        Kloc = Kloc + Kappa(imat)*MATMUL(MATMUL(TRANSPOSE(B),KappaMatrx),B) &
             *0.25d0*vol

!            *0.04166666666666667d0*vol ! w(i) * is this correct fix check
! take the product K*T

!!$        RinLoc(:) = MATMUL(Kloc,Tloc)

!!$        r1 = r1 - RinLoc(1)
!!$        r2 = r2 - RinLoc(2)
!!$        r3 = r3 - RinLoc(3)
!!$        r4 = r4 - RinLoc(4)
!!$        r5 = r5 - RinLoc(5)
!!$        r6 = r6 - RinLoc(6)
!!$        r7 = r7 - RinLoc(7)
!!$        r8 = r8 - RinLoc(8)
!!$        r9 = r9 - RinLoc(9)
!!$        r10 = r10 - RinLoc(10)

! Mesh Motion term
        
        MeshVelNd(1) = ShapFnct(1)*MeshVel(k1n1) + &
             ShapFnct(2)*MeshVel(k1n2) + &
             ShapFnct(3)*MeshVel(k1n3) + &
             ShapFnct(4)*MeshVel(k1n4) + &
             ShapFnct(5)*MeshVel(k1n5) + &
             ShapFnct(6)*MeshVel(k1n6) + &
             ShapFnct(7)*MeshVel(k1n7) + &
             ShapFnct(8)*MeshVel(k1n8) + &
             ShapFnct(9)*MeshVel(k1n9) + &
             ShapFnct(10)*MeshVel(k1n10)

        MeshVelNd(2) = ShapFnct(1)*MeshVel(k2n1) + &
             ShapFnct(2)*MeshVel(k2n2) + &
             ShapFnct(3)*MeshVel(k2n3) + &
             ShapFnct(4)*MeshVel(k2n4) + &
             ShapFnct(5)*MeshVel(k2n5) + &
             ShapFnct(6)*MeshVel(k2n6) + &
             ShapFnct(7)*MeshVel(k2n7) + &
             ShapFnct(8)*MeshVel(k2n8) + &
             ShapFnct(9)*MeshVel(k2n9) + &
             ShapFnct(10)*MeshVel(k2n10)

        MeshVelNd(3) = ShapFnct(1)*MeshVel(k3n1) + &
             ShapFnct(2)*MeshVel(k3n2) + &
             ShapFnct(3)*MeshVel(k3n3) + &
             ShapFnct(4)*MeshVel(k3n4) + &
             ShapFnct(5)*MeshVel(k3n5) + &
             ShapFnct(6)*MeshVel(k3n6) + &
             ShapFnct(7)*MeshVel(k3n7) + &
             ShapFnct(8)*MeshVel(k3n8) + &
             ShapFnct(9)*MeshVel(k3n9) + &
             ShapFnct(10)*MeshVel(k3n10)

        BT(1,1) = ShapFnct(1)*MeshVelNd(1)
        BT(1,2) = ShapFnct(1)*MeshVelNd(2)
        BT(1,3) = ShapFnct(1)*MeshVelNd(3)
        BT(2,1) = ShapFnct(2)*MeshVelNd(1)
        BT(2,2) = ShapFnct(2)*MeshVelNd(2)
        BT(2,3) = ShapFnct(2)*MeshVelNd(3)
        BT(3,1) = ShapFnct(3)*MeshVelNd(1)
        BT(3,2) = ShapFnct(3)*MeshVelNd(2)
        BT(3,3) = ShapFnct(3)*MeshVelNd(3)
        BT(4,1) = ShapFnct(4)*MeshVelNd(1)
        BT(4,2) = ShapFnct(4)*MeshVelNd(2)
        BT(4,3) = ShapFnct(4)*MeshVelNd(3)
        BT(5,1) = ShapFnct(5)*MeshVelNd(1)
        BT(5,2) = ShapFnct(5)*MeshVelNd(2)
        BT(5,3) = ShapFnct(5)*MeshVelNd(3)
        BT(6,1) = ShapFnct(6)*MeshVelNd(1)
        BT(6,2) = ShapFnct(6)*MeshVelNd(2)
        BT(6,3) = ShapFnct(6)*MeshVelNd(3)
        BT(7,1) = ShapFnct(7)*MeshVelNd(1)
        BT(7,2) = ShapFnct(7)*MeshVelNd(2)
        BT(7,3) = ShapFnct(7)*MeshVelNd(3)
        BT(8,1) = ShapFnct(8)*MeshVelNd(1)
        BT(8,2) = ShapFnct(8)*MeshVelNd(2)
        BT(8,3) = ShapFnct(8)*MeshVelNd(3)
        BT(9,1) = ShapFnct(9)*MeshVelNd(1)
        BT(9,2) = ShapFnct(9)*MeshVelNd(2)
        BT(9,3) = ShapFnct(9)*MeshVelNd(3)
        BT(10,1) = ShapFnct(10)*MeshVelNd(1)
        BT(10,2) = ShapFnct(10)*MeshVelNd(2)
        BT(10,3) = ShapFnct(10)*MeshVelNd(3)

        Kloc =  Kloc - Rho(imat)*Cp(imat)*MATMUL(BT,B)*0.25d0*vol

! take the product K*T

!!$!        RinLoc(:) = MATMUL(Kloc,Tloc)
!!$
!!$        r1 = r1 - RinLoc(1)
!!$        r2 = r2 - RinLoc(2)
!!$        r3 = r3 - RinLoc(3)
!!$        r4 = r4 - RinLoc(4)
!!$        r5 = r5 - RinLoc(5)
!!$        r6 = r6 - RinLoc(6)
!!$        r7 = r7 - RinLoc(7)
!!$        r8 = r8 - RinLoc(8)
!!$        r9 = r9 - RinLoc(9)
!!$        r10 = r10 - RinLoc(10)
     ENDDO


!!$     DO j = 1,10
!!$        WRITE(*,'(10(x,f10.1))')Kloc(j,1:10)
!!$     ENDDO
     RinLoc(:) = MATMUL(Kloc,Tloc)

!!$     DO j = 1,10
!!$        WRITE(*,'(10(x,f10.1))')Kloc(j,1:10)
!!$     ENDDO
!!$     PRINT*,'R'
!!$     WRITE(*,'(10(x,f10.1))') Tloc(1:10)
     

     r1 =  RinLoc(1)
     r2 =  RinLoc(2)
     r3 =  RinLoc(3)
     r4 =  RinLoc(4)
     r5 =  RinLoc(5)
     r6 =  RinLoc(6)
     r7 =  RinLoc(7)
     r8 =  RinLoc(8)
     r9 =  RinLoc(9)
     r10 =  RinLoc(10)

     Rnet(n1)  = Rnet(n1)  - r1 !*0.04166666666666667d0
     Rnet(n2)  = Rnet(n2)  - r2 !*0.04166666666666667d0
     Rnet(n3)  = Rnet(n3)  - r3 !*0.04166666666666667d0
     Rnet(n4)  = Rnet(n4)  - r4 !*0.04166666666666667d0
     Rnet(n5)  = Rnet(n5)  - r5 !*0.04166666666666667d0
     Rnet(n6)  = Rnet(n6)  - r6 !*0.04166666666666667d0
     Rnet(n7)  = Rnet(n7)  - r7 !*0.04166666666666667d0
     Rnet(n8)  = Rnet(n8)  - r8 !*0.04166666666666667d0
     Rnet(n9)  = Rnet(n9)  - r9 !*0.04166666666666667d0
     Rnet(n10)  = Rnet(n10) - r10 !*0.04166666666666667d0

     
  END DO

  RETURN

END SUBROUTINE v3d10_thermal

