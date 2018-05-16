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
  SUBROUTINE V3D4_NL_ARRUDA_BOYCE(coor,matcstet,lmcstet,R_in,d, &
      S11,S22,S33,S12,S23,S13, &
      numnp,nstart,nend,numcstet,numat_vol, &
      mu,kapp)

!!****f* Rocfrac/Rocfrac/Source/v3d4_NL_ARRUDA_BOYCE.f90
!!
!!  NAME
!!    v3d4_NL_ARRUDA_BOYCE
!!
!!  FUNCTION
!!
!!    Computes the internal force for a large deformation, 
!!    hyperelastic arruda-boyce 4-node tetrahedral.
!!
!!  INPUTS
!!
!!   NumNP    -- Number of nodes
!!   numcstet -- Number of elements
!!   Coor     -- number of coordinates
!!   Matcstet  -- Material id
!!   d -- nodal displacement
!!   mu -- shear modulus
!!   kappa -- bulk modulus
!!
!!   lmcstet  -- Nodal connectivity
!!   nstart, nend -- element beginning and end loop counter
!!   numat_vol -- number of materials
!!
!!  OUTPUT
!!
!!    Rin -- internal force vector
!!    S11,S22,S33,S12,S23,S13 -- Stresses at gauss points
!!
!!****


      IMPLICIT NONE
!---- Global variables
      INTEGER :: numnp          ! number of nodal points
      INTEGER :: numcstet       ! number of CSTet elements
      INTEGER :: numat_vol      ! number of volumetric materials
!--   coordinate array
      REAL*8, DIMENSION(1:3,1:numnp) :: coor
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
!--  Coordinate subtractions
      REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  Coordinate subtractions: These are to speed up B calculation
      REAL*8 :: x12, x13, y12, y13, z12, z13
      REAL*8 :: val11, val21, val31
! -- 
      REAL*8 :: F11, F12, F13, F21, F22, F23, F31, F32, F33
      REAL*8, DIMENSION(1:numat_vol) :: mu, kapp
      REAL*8, DIMENSION(1:3,1:3) :: Fij, Cij

      INTEGER :: kk,ll,mm

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

         Vx6 = -( x14*val11 + y14*val21 + z14*val31 )

         Vx6inv = 1.d0 / Vx6

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
!
! Arruda-Boyce Nonlinear Elasticity Model 
!             
! -- NOTE: ci(7,j) : shear modulus 


!--   (1) Deformation Gradient Fij 
  
         Fij(1,1) = F11
         Fij(1,2) = F12
         Fij(1,3) = F13 
         Fij(2,1) = F21 
         Fij(2,2) = F22
         Fij(2,3) = F23 
         Fij(3,1) = F31
         Fij(3,2) = F32
         Fij(3,3) = F33

!--  (2) right Cauchy-Green deformation Tensor Cij = Fmi Fmj 
!             T
!        C = F  F
 
         DO kk=1,3
            DO ll=1,3
               Cij(kk,ll) = 0.d0
               DO mm=1,3
                  Cij(kk,ll)=Cij(kk,ll)+Fij(mm,kk)*Fij(mm,ll)
               ENDDO
            ENDDO
         ENDDO

         CALL ARRUDA_BOYCE(Cij, &
         S11(i),S22(i),S33(i),S12(i),S23(i),S13(i),i, &
         mu(j),kapp(j))
         
! calculate the volume

         vol = Vx6 / 6.d0

! ASSEMBLE THE INTERNAL FORCE VECTOR
!
! local node 1
         R_in(k1n1) = R_in(k1n1) - vol* &
             ( S11(i)*B1*F11 + S22(i)*B2*F12 + S33(i)*B3*F13 &
             + S12(i)*( B2*F11 + B1*F12 ) &
             + S23(i)*( B3*F12 + B2*F13 ) &
             + S13(i)*( B3*F11 + B1*F13 ) )
         R_in(k2n1) = R_in(k2n1) - vol* &
             ( S11(i)*B1*F21 + S22(i)*B2*F22 + S33(i)*B3*F23 &
             + S12(i)*( B1*F22 + B2*F21 ) &
             + S23(i)*( B3*F22 + B2*F23 ) &
             + S13(i)*( B3*F21 + B1*F23 ) )
         R_in(k3n1)   = R_in(k3n1)   - vol* &
             ( S11(i)*B1*F31 + S22(i)*B2*F32 + S33(i)*B3*F33 &
             + S12(i)*( B2*F31 + B1*F32 ) &
             + S23(i)*( B3*F32 + B2*F33 ) &
             + S13(i)*( B3*F31 + B1*F33 ) )
! local node 2 
         R_in(k1n2) = R_in(k1n2) - vol* &
             ( S11(i)*B4*F11 + S22(i)*B5*F12 + S33(i)*B6*F13 &
             + S12(i)*( B5*F11 + B4*F12 ) &
             + S23(i)*( B6*F12 + B5*F13 ) &
             + S13(i)*( B6*F11 + B4*F13 ) )
         R_in(k2n2) = R_in(k2n2) - vol* &
             ( S11(i)*B4*F21 + S22(i)*B5*F22 + S33(i)*B6*F23 &
             + S12(i)*( B4*F22 + B5*F21 ) &
             + S23(i)*( B6*F22 + B5*F23 ) &
             + S13(i)*( B6*F21 + B4*F23 ) )
         R_in(k3n2)   = R_in(k3n2)   - vol* &
             ( S11(i)*B4*F31 + S22(i)*B5*F32 + S33(i)*B6*F33 &
             + S12(i)*( B5*F31 + B4*F32 ) &
             + S23(i)*( B6*F32 + B5*F33 ) &
             + S13(i)*( B6*F31 + B4*F33 ) )
! local node 3 
         R_in(k1n3) = R_in(k1n3) - vol* &
             ( S11(i)*B7*F11 + S22(i)*B8*F12 + S33(i)*B9*F13 &
             + S12(i)*( B8*F11 + B7*F12 ) &
             + S23(i)*( B9*F12 + B8*F13 ) &
             + S13(i)*( B9*F11 + B7*F13 ) )
         R_in(k2n3) = R_in(k2n3) - vol* &
             ( S11(i)*B7*F21 + S22(i)*B8*F22 + S33(i)*B9*F23 &
             + S12(i)*( B7*F22 + B8*F21 ) &
             + S23(i)*( B9*F22 + B8*F23 ) &
             + S13(i)*( B9*F21 + B7*F23 ) )
         R_in(k3n3)   = R_in(k3n3)   - vol* &
             ( S11(i)*B7*F31 + S22(i)*B8*F32 + S33(i)*B9*F33 &
             + S12(i)*( B8*F31 + B7*F32 ) &
             + S23(i)*( B9*F32 + B8*F33 ) &
             + S13(i)*( B9*F31 + B7*F33 ) )
! local node 4         
         R_in(k1n4) = R_in(k1n4) - vol* &
             ( S11(i)*B10*F11 + S22(i)*B11*F12+S33(i)*B12*F13 &
             + S12(i)*( B11*F11 + B10*F12 ) &
             + S23(i)*( B12*F12 + B11*F13 ) &
             + S13(i)*( B12*F11 + B10*F13 ) )
         R_in(k2n4) = R_in(k2n4) - vol* &
             ( S11(i)*B10*F21 + S22(i)*B11*F22+S33(i)*B12*F23 &
             + S12(i)*( B10*F22 + B11*F21 ) &
             + S23(i)*( B12*F22 + B11*F23 ) &
             + S13(i)*( B12*F21 + B10*F23 ) )
         R_in(k3n4)   = R_in(k3n4)   - vol* &
             ( S11(i)*B10*F31 + S22(i)*B11*F32+S33(i)*B12*F33 &
             + S12(i)*( B11*F31 + B10*F32 ) &
             + S23(i)*( B12*F32 + B11*F33 ) &
             + S13(i)*( B12*F31 + B10*F33 ) )

      x1 = coor(1,n1) + u1
      x2 = coor(1,n2) + u2
      x3 = coor(1,n3) + u3
      x4 = coor(1,n4) + u4
      y1 = coor(2,n1) + v1
      y2 = coor(2,n2) + v2
      y3 = coor(2,n3) + v3
      y4 = coor(2,n4) + v4
      z1 = coor(3,n1) + w1
      z2 = coor(3,n2) + w2
      z3 = coor(3,n3) + w3
      z4 = coor(3,n4) + w4
     
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
     
      Vx6 = -( x14*val11 + y14*val21 + z14*val31 )

      IF(Vx6.LT.0.d0) THEN
         PRINT*,'ROCFRAC ::        ERROR * ERROR * ERROR'
         PRINT*,'ROCFRAC :: DEFORMED ELEMENT VOLUME TURNED NEGATIVE'
         STOP
      ENDIF
          
      ENDDO

      RETURN
      END

