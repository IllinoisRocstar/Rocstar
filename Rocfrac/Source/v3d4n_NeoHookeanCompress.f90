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
SUBROUTINE V3D4n_NeoHookeanCompress(numnp, numel, coor, disp, nodes, Rnet,&
     NumElNeigh, ElConn, alpha, Ahat,NumMatVol, &
     xmu,xlambda)
  
!     nprocs,TotNumNdComm,TotNumNeighProcs,NeighProcList,NumNdComm,neigh_lst,&
!     MPI_STATUS_SIZE,MPI_COMM_ROCFRAC,MPI_DOUBLE_PRECISION, &
!     ReqRcv,ReqSnd,StatRcv, StatSnd)


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

  USE ROCSTAR_RocFracComm 

  IMPLICIT NONE
  
  integer :: NumMatVol
  INTEGER :: numnp, numel
  REAL*8,dimension(1:NumMatVol) :: xmu, xlambda
  REAL*8, DIMENSION(1:numnp) :: Ahat
  REAL*8, DIMENSION(1:3*numnp) :: disp, Rnet
  REAL*8, DIMENSION(1:3,1:numnp) :: coor
  INTEGER, DIMENSION(1:numnp) ::NumElNeigh
  INTEGER, DIMENSION(1:numnp,1:40) :: ElConn  ! fix 40 should not be hard coded
  INTEGER,DIMENSION(1:4,1:numel) :: nodes
  REAL*8, DIMENSION(1:4,1:numel) :: alpha
  REAL*8, DIMENSION(1:numnp) :: S11n, S22n, S33n, S12n, S13n,S23n

  INTEGER :: i,k
!--  x, y and z displacements of nodes
  REAL*8 :: u1,u2,u3,u4,v1,v2,v3,v4,w1,w2,w3,w4
!--   coordinate holding variable
  REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
!--  Coordinate subtractions

  REAL*8 ::  z12, z13
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
  REAL*8 :: x12, x13, y12, y13

!--   6*volume, inverse(6*volume),  and the volume      
  REAL*8 :: Vx6, Vx6inv, vol
!--   spacial derivatives
  REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12
!--   partial derivatives of the displacement 
  REAL*8 :: dudx,dvdy,dwdz,dudy,dvdx,dvdz,dwdy,dudz,dwdx
!--   strains
  REAL*8 :: VolNd,VaInv
! -- node numbers
  INTEGER :: n1,n2,n3,n4
! -- dummy and counters
  INTEGER :: j,nstart,nend
  INTEGER :: k1n1,k1n2,k1n3,k1n4,k2n1,k2n2,k2n3,k2n4
  INTEGER :: k3n1,k3n2,k3n3,k3n4
! -- 
  REAL*8 :: F11, F12, F13, F21, F22, F23, F31, F32, F33
  REAL*8 :: J1,J2,J2inv
  REAL*8 :: C11, C12, C13, C21, C22, C23, C31, C32, C33
  INTEGER :: jj,iElNum
  REAL*8 :: VolNd0Inv,VolEl0
  REAL*8, DIMENSION(1:1) :: s11, s22, s33, s12, s23, s13
  INTEGER :: k3i,k2i,k1i
  INTEGER :: ix

! /* For each Node */

  DO i = 1, numnp ! for each node

! /* Undeformed Volume of node */

     VolNd0Inv = 1.d0/Ahat(i)

! /* Initialize to 0, Va, Va0, Fa */

     F11 = 0.d0
     F22 = 0.d0
     F33 = 0.d0
     F12 = 0.d0
     F13 = 0.d0
     F21 = 0.d0
     F23 = 0.d0
     F31 = 0.d0
     F32 = 0.d0
     VolNd = 0.d0
       
     DO j = 1, NumElNeigh(i)
        
        iElNum = ElConn(i,j)

        n1 = nodes(1,iElNum)
        n2 = nodes(2,iElNum)
        n3 = nodes(3,iElNum)
        n4 = nodes(4,iElNum)

        IF(i.EQ.n1) THEN
           ix = 1
        ELSE IF(i.EQ.n2)THEN
           ix = 2
        ELSE IF(i.EQ.n3)THEN
           ix = 3
        ELSE IF(i.EQ.n4)THEN
           ix = 4
        ENDIF

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
        u1 = disp(k1n1)           ! (3*n1-2)
        u2 = disp(k1n2)           ! (3*n2-2)
        u3 = disp(k1n3)           ! (3*n3-2)
        u4 = disp(k1n4)           ! (3*n4-2)
        v1 = disp(k2n1)           ! (3*n1-1)
        v2 = disp(k2n2)           ! (3*n2-1)
        v3 = disp(k2n3)           ! (3*n3-1)
        v4 = disp(k2n4)           ! (3*n4-1)
        w1 = disp(k3n1)           ! (3*n1)
        w2 = disp(k3n2)           ! (3*n2)
        w3 = disp(k3n3)           ! (3*n3)
        w4 = disp(k3n4)           ! (3*n4)
            
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

        x12 = x1 - x2 ! not used in vol. calc
        x13 = x1 - x3 ! not used in vol. calc
        x14 = x1 - x4
        x24 = x2 - x4
        x34 = x3 - x4
        y12 = y1 - y2 ! not used in vol. calc
        y13 = y1 - y3 ! not used in vol. calc
        y14 = y1 - y4
        y24 = y2 - y4
        y34 = y3 - y4
        z12 = z1 - z2 ! not used in vol. calc
        z13 = z1 - z3 ! not used in vol. calc
        z14 = z1 - z4
        z24 = z2 - z4
        z34 = z3 - z4

        c11 =    y24*z34 - z24*y34
        c21 = -( x24*z34 - z24*x34 )
        c31 =    x24*y34 - y24*x34

        Vx6 = -( x14*c11 + y14*c21 + z14*c31 )

        Vx6inv = 1.d0/Vx6

        VolEl0 = Vx6/6.d0 ! undeformed volume of element (Ve_O)
            

! See the maple worksheet 'V3D4.mws' for the derivation of Vx6


! See the maple worksheet 'V3D4.mws' for the derivation of [B]
! Compute the Shape functions
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
        F11 = F11 + alpha(ix,IElNum)*VolEl0*(1.d0 + dudx)
        F12 = F12 + alpha(ix,IElNum)*VolEl0*dudy
        F13 = F13 + alpha(ix,IElNum)*VolEl0*dudz
        F21 = F21 + alpha(ix,IElNum)*VolEl0*dvdx
        F22 = F22 + alpha(ix,IElNum)*VolEl0*(1.d0 + dvdy)
        F23 = F23 + alpha(ix,IElNum)*VolEl0*dvdz
        F31 = F31 + alpha(ix,IElNum)*VolEl0*dwdx
        F32 = F32 + alpha(ix,IElNum)*VolEl0*dwdy
        F33 = F33 + alpha(ix,IElNum)*VolEl0*(1.d0 + dwdz)

     ENDDO

! 2nd "for each node" in box starts here


  ! first bullet in box
     
     F11 = VolNd0Inv*F11
     F22 = VolNd0Inv*F22
     F33 = VolNd0Inv*F33
     F12 = VolNd0Inv*F12
     F21 = VolNd0Inv*F21
     F23 = VolNd0Inv*F23
     F32 = VolNd0Inv*F32
     F13 = VolNd0Inv*F13
     F31 = VolNd0Inv*F31

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

     S11n(i) = xmu(1)*(1.d0-(C22*C33-C23*C32)) + xlambda(1)*LOG(J1)*(C22*C33-C23*C32)
     S22n(i) = xmu(1)*(1.d0-(C11*C33-C31*C13)) + xlambda(1)*LOG(J1)*(C11*C33-C31*C13)
     S33n(i) = xmu(1)*(1.d0-(C11*C22-C12*C21)) + xlambda(1)*LOG(J1)*(C11*C22-C12*C21)
     S12n(i) = (-C12*C33+C13*C32)*(xmu(1) - xlambda(1)*LOG(J1))
     S23n(i) = (-C11*C23+C13*C21)*(xmu(1) - xlambda(1)*LOG(J1))
     S13n(i) = (C12*C23-C13*C22)*(xmu(1) - xlambda(1)*LOG(J1))
     
  ENDDO

! For each node

  DO i = 1, numnp

     k3i = 3*i
     k2i = 3*i-1
     k1i = 3*i-2

! For each element

     DO j = 1, NumElNeigh(i)
        
        iElNum = ElConn(i,j)
        
        n1 = nodes(1,iElNum)
        n2 = nodes(2,iElNum)
        n3 = nodes(3,iElNum)
        n4 = nodes(4,iElNum)  
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
        u1 = disp(k1n1)           ! (3*n1-2)
        u2 = disp(k1n2)           ! (3*n2-2)
        u3 = disp(k1n3)           ! (3*n3-2)
        u4 = disp(k1n4)           ! (3*n4-2)
        v1 = disp(k2n1)           ! (3*n1-1)
        v2 = disp(k2n2)           ! (3*n2-1)
        v3 = disp(k2n3)           ! (3*n3-1)
        v4 = disp(k2n4)           ! (3*n4-1)
        w1 = disp(k3n1)           ! (3*n1)
        w2 = disp(k3n2)           ! (3*n2)
        w3 = disp(k3n3)           ! (3*n3)
        w4 = disp(k3n4)           ! (3*n4)
        
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
        
        x12 = x1 - x2 ! not used in vol. calc
        x13 = x1 - x3 ! not used in vol. calc
        x14 = x1 - x4
        x24 = x2 - x4
        x34 = x3 - x4
        y12 = y1 - y2 ! not used in vol. calc
        y13 = y1 - y3 ! not used in vol. calc
        y14 = y1 - y4
        y24 = y2 - y4
        y34 = y3 - y4
        z12 = z1 - z2 ! not used in vol. calc
        z13 = z1 - z3 ! not used in vol. calc
        z14 = z1 - z4
        z24 = z2 - z4
        z34 = z3 - z4
        
        c11 =    y24*z34 - z24*y34
        c21 = -( x24*z34 - z24*x34 )
        c31 =    x24*y34 - y24*x34
        
        Vx6 = -( x14*c11 + y14*c21 + z14*c31 )
        VolEl0 = Vx6/6.d0

! See the maple worksheet 'V3D4.mws' for the derivation of Vx6

        Vx6inv = 1.d0 / Vx6

! See the maple worksheet 'V3D4.mws' for the derivation of [B]! Compute the Shape functions
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
        F11 = (1.d0 + dudx)
        F12 = dudy
        F13 = dudz
        F21 = dvdx
        F22 = (1.d0 + dvdy)
        F23 = dvdz
        F31 = dwdx
        F32 = dwdy
        F33 = (1.d0 + dwdz)

!             4
!            ----   ,
! Evaluate   \     P 
!            /      a
!            ----
!            a = 1

        S11(1) = 0.d0
        S22(1) = 0.d0
        S33(1) = 0.d0
        S12(1) = 0.d0
        S23(1) = 0.d0
        S13(1) = 0.d0
        
        DO k = 1, 4
           S11(1) = S11(1) + alpha(k,IElNum)*S11n(nodes(k,iElNum))
           S22(1) = S22(1) + alpha(k,IElNum)*S22n(nodes(k,iElNum))
           S33(1) = S33(1) + alpha(k,IElNum)*S33n(nodes(k,iElNum))
           S12(1) = S12(1) + alpha(k,IElNum)*S12n(nodes(k,iElNum))
           S23(1) = S23(1) + alpha(k,IElNum)*S23n(nodes(k,iElNum))
           S13(1) = S13(1) + alpha(k,IElNum)*S13n(nodes(k,iElNum)) 
        ENDDO

! calculate the volume

        vol = Vx6 / 6.d0

!!$
!!$! ASSEMBLE THE INTERNAL FORCE VECTOR
!!$!
!!$! local node 1
!!$
!!$         Rnet(k1i) = Rnet(k1i) - VolEl0* &
!!$              ( S11e*B1*(1.d0+dudx) + S22e*B2*dudy + S33e*B3*dudz &
!!$              + S12e*( B2*(1.d0+dudx) + B1*dudy ) &
!!$              + S23e*( B3*dudy + B2*dudz ) &
!!$              + S13e*( B3*(1.d0+dudx) + B1*dudz ) )
!!$            Rnet(k2i) = Rnet(k2i) - VolEl0* &
!!$                 ( S11e*B1*dvdx + S22e*B2*(1.d0+dvdy) + S33e*B3*dvdz &
!!$                 + S12e*( B1*(1.d0+dvdy) + B2*dvdx ) &
!!$                 + S23e*( B3*(1.d0+dvdy) + B2*dvdz ) &
!!$                 + S13e*( B3*dvdx + B1*dvdz ) )
!!$            Rnet(k3i)   = Rnet(k3i)   - VolEl0* &
!!$                 ( S11e*B1*dwdx + S22e*B2*dwdy + S33e*B3*(1.d0+dwdz) &
!!$                 + S12e*( B2*dwdx + B1*dwdy ) &
!!$                 + S23e*( B3*dwdy + B2*(1.d0 + dwdz) ) &
!!$                 + S13e*( B3*dwdx + B1*(1.d0 + dwdz) ) ) 
!!$! local node 2 
!!$            Rnet(k1i) = Rnet(k1i) - VolEl0* &
!!$                 ( S11e*B4*(1.d0+dudx) + S22e*B5*dudy + S33e*B6*dudz &
!!$                 + S12e*( B5*(1.d0+dudx) + B4*dudy ) &
!!$                 + S23e*( B6*dudy + B5*dudz ) &
!!$                 + S13e*( B6*(1.d0+dudx) + B4*dudz ) ) 
!!$            Rnet(k2i) = Rnet(k2i) - VolEl0* &
!!$                 ( S11e*B4*dvdx + S22e*B5*(1.d0+dvdy) + S33e*B6*dvdz &
!!$                 + S12e*( B4*(1.d0+dvdy) + B5*dvdx ) &
!!$                 + S23e*( B6*(1.d0+dvdy) + B5*dvdz ) &
!!$                 + S13e*( B6*dvdx + B4*dvdz ) )
!!$            Rnet(k3i)   = Rnet(k3i)   - VolEl0* &
!!$                 ( S11e*B4*dwdx + S22e*B5*dwdy + S33e*B6*(1.d0+dwdz) &
!!$                 + S12e*( B5*dwdx + B4*dwdy ) &
!!$                 + S23e*( B6*dwdy + B5*(1.d0 + dwdz) ) &
!!$                 + S13e*( B6*dwdx + B4*(1.d0 + dwdz) ) )
!!$! local node 3 
!!$            Rnet(k1i) = Rnet(k1i) - VolEl0* &
!!$                 ( S11e*B7*(1.d0+dudx) + S22e*B8*dudy + S33e*B9*dudz &
!!$                 + S12e*( B8*(1.d0+dudx) + B7*dudy ) &
!!$                 + S23e*( B9*dudy + B8*dudz ) &
!!$                 + S13e*( B9*(1.d0+dudx) + B7*dudz ) )
!!$            Rnet(k2i) = Rnet(k2i) - VolEl0* &
!!$                 ( S11e*B7*dvdx + S22e*B8*(1.d0+dvdy) + S33e*B9*dvdz &
!!$                 + S12e*( B7*(1.d0+dvdy) + B8*dvdx ) &
!!$                 + S23e*( B9*(1.d0+dvdy) + B8*dvdz ) &
!!$                 + S13e*( B9*dvdx + B7*dvdz ) )
!!$            Rnet(k3i)   = Rnet(k3i)   - VolEl0* &
!!$                 ( S11e*B7*dwdx + S22e*B8*dwdy + S33e*B9*(1.d0+dwdz) &
!!$                 + S12e*( B8*dwdx + B7*dwdy ) &
!!$                 + S23e*( B9*dwdy + B8*(1.d0 + dwdz) ) &
!!$                 + S13e*( B9*dwdx + B7*(1.d0 + dwdz) ) )
!!$! local node 4         
!!$            Rnet(k1i) = Rnet(k1i) - VolEl0* &
!!$                 ( S11e*B10*(1.d0+dudx) + S22e*B11*dudy+S33e*B12*dudz &
!!$                 + S12e*( B11*(1.d0+dudx) + B10*dudy ) &
!!$                 + S23e*( B12*dudy + B11*dudz ) &
!!$                 + S13e*( B12*(1.d0+dudx) + B10*dudz ) )
!!$            Rnet(k2i) = Rnet(k2i) - VolEl0* &
!!$                 ( S11e*B10*dvdx + S22e*B11*(1.d0+dvdy)+S33e*B12*dvdz &
!!$                 + S12e*( B10*(1.d0+dvdy) + B11*dvdx ) &
!!$                 + S23e*( B12*(1.d0+dvdy) + B11*dvdz ) &
!!$                 + S13e*( B12*dvdx + B10*dvdz ) ) 
!!$            Rnet(k3i)   = Rnet(k3i)   - VolEl0* &
!!$                 ( S11e*B10*dwdx + S22e*B11*dwdy+S33e*B12*(1.d0+dwdz) &
!!$                 + S12e*( B11*dwdx + B10*dwdy ) &
!!$                 + S23e*( B12*dwdy + B11*(1.d0 + dwdz) ) &
!!$                 + S13e*( B12*dwdx + B10*(1.d0 + dwdz) ) )

        Rnet(k1n1) = Rnet(k1n1) - VolEl0* &
             ( S11(1)*B1*F11 + S22(1)*B2*F12 + S33(1)*B3*F13 &
             + S12(1)*( B2*F11 + B1*F12 ) &
             + S23(1)*( B3*F12 + B2*F13 ) &
             + S13(1)*( B3*F11 + B1*F13 ) )
        Rnet(k2n1) = Rnet(k2n1) - VolEl0* &
             ( S11(1)*B1*F21 + S22(1)*B2*F22 + S33(1)*B3*F23 &
             + S12(1)*( B1*F22 + B2*F21 ) &
             + S23(1)*( B3*F22 + B2*F23 ) &
             + S13(1)*( B3*F21 + B1*F23 ) )
        Rnet(k3n1) = Rnet(k3n1) - VolEl0* &
             ( S11(1)*B1*F31 + S22(1)*B2*F32 + S33(1)*B3*F33 &
             + S12(1)*( B2*F31 + B1*F32 ) &
             + S23(1)*( B3*F32 + B2*F33 ) &
             + S13(1)*( B3*F31 + B1*F33 ) )
! local node 2 
        Rnet(k1n2) = Rnet(k1n2) - VolEl0* &
             ( S11(1)*B4*F11 + S22(1)*B5*F12 + S33(1)*B6*F13 &
             + S12(1)*( B5*F11 + B4*F12 ) &
             + S23(1)*( B6*F12 + B5*F13 ) &
             + S13(1)*( B6*F11 + B4*F13 ) )
        Rnet(k2n2) = Rnet(k2n2) - VolEl0* &
             ( S11(1)*B4*F21 + S22(1)*B5*F22 + S33(1)*B6*F23 &
             + S12(1)*( B4*F22 + B5*F21 ) &
             + S23(1)*( B6*F22 + B5*F23 ) &
             + S13(1)*( B6*F21 + B4*F23 ) )
        Rnet(k3n2) = Rnet(k3n2) - VolEl0* &
             ( S11(1)*B4*F31 + S22(1)*B5*F32 + S33(1)*B6*F33 &
             + S12(1)*( B5*F31 + B4*F32 ) &
             + S23(1)*( B6*F32 + B5*F33 ) &
             + S13(1)*( B6*F31 + B4*F33 ) )
! local node 3 
        Rnet(k1n3) = Rnet(k1n3) - VolEl0* &
             ( S11(1)*B7*F11 + S22(1)*B8*F12 + S33(1)*B9*F13 &
             + S12(1)*( B8*F11 + B7*F12 ) &
             + S23(1)*( B9*F12 + B8*F13 ) &
             + S13(1)*( B9*F11 + B7*F13 ) )
        Rnet(k2n3) = Rnet(k2n3) - VolEl0* &
             ( S11(1)*B7*F21 + S22(1)*B8*F22 + S33(1)*B9*F23 &
             + S12(1)*( B7*F22 + B8*F21 ) &
             + S23(1)*( B9*F22 + B8*F23 ) &
             + S13(1)*( B9*F21 + B7*F23 ) )
        Rnet(k3n3) = Rnet(k3n3) - VolEl0* &
             ( S11(1)*B7*F31 + S22(1)*B8*F32 + S33(1)*B9*F33 &
             + S12(1)*( B8*F31 + B7*F32 ) &
             + S23(1)*( B9*F32 + B8*F33 ) &
             + S13(1)*( B9*F31 + B7*F33 ) )
! local node 4         
        Rnet(k1n4) = Rnet(k1n4) - VolEl0* &
             ( S11(1)*B10*F11 + S22(1)*B11*F12+S33(1)*B12*F13 &
             + S12(1)*( B11*F11 + B10*F12 ) &
             + S23(1)*( B12*F12 + B11*F13 ) &
             + S13(1)*( B12*F11 + B10*F13 ) )
        Rnet(k2n4) = Rnet(k2n4) - VolEl0* &
             ( S11(1)*B10*F21 + S22(1)*B11*F22 + S33(1)*B12*F23 &
             + S12(1)*( B10*F22 + B11*F21 ) &
             + S23(1)*( B12*F22 + B11*F23 ) &
             + S13(1)*( B12*F21 + B10*F23 ) )
        Rnet(k3n4) = Rnet(k3n4) - VolEl0* &
             ( S11(1)*B10*F31 + S22(1)*B11*F32 + S33(1)*B12*F33 &
             + S12(1)*( B11*F31 + B10*F32 ) &
             + S23(1)*( B12*F32 + B11*F33 ) &
             + S13(1)*( B12*F31 + B10*F33 ) )
        
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE V3D4n_NeoHookeanCompress

