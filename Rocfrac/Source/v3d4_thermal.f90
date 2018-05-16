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
SUBROUTINE V3D4_thermal(NumEl, NumNP, ElConnVol, Coor, Kappa, &
     Rnet, T, Rho,Cp, matcstet,numat_vol,MeshVelo,ElemStart, ElemEnd)

!!****f* Rocfrac/Source/v3d4_thermal
!!
!!  NAME
!!    v3d4_thermal.f90
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

  IMPLICIT NONE

!---- Global variables
  INTEGER :: numnp          ! number of nodal points
  INTEGER :: NumEl       ! number of CSTet elements
  INTEGER :: numat_vol      ! number of volumetric materials
  integer :: ElemStart, ElemEnd

  REAL*8, DIMENSION(1:numat_vol) :: Kappa
  REAL*8, DIMENSION(1:NumNP) :: T
  REAL*8, DIMENSION(1:numat_vol) :: Rho,Cp
!--   coordinate array
  REAL*8, DIMENSION(1:3,1:numnp) :: coor

  REAL*8, DIMENSION(1:3*numnp) :: MeshVelo
!--   internal force
  REAL*8, DIMENSION(1:numnp) :: Rnet
!--   connectivity table for CSTet  
  INTEGER, DIMENSION(1:4,1:NumEl) :: ElConnVol
!--   mat number for CSTet element
  INTEGER, DIMENSION(1:NumEl) :: matcstet
!---- Local variables
!--   node numbers
  INTEGER :: n1,n2,n3,n4
!--   x, y, and z displacements of nodes
  REAL*8 :: u1,u2,u3,u4,v1,v2,v3,v4,w1,w2,w3,w4
!--   6*volume, inverse(6*volume),  and the volume      
  REAL*8 :: Vx6, Vx6inv, vol
!--   spacial derivatives
  REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12
!--   strains
  REAL*8 :: E11,E22,E33,E12,E23,E13
!--   coordinate holding variable
  REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  Coordinate subtractions: These are to speed up B calculation
  REAL*8 :: x12, x13, y12, y13, z12, z13
!-- Dummy
  REAL*8 :: c11, c21, c31
!--   dummy and counters
  INTEGER :: i,j,nstart,nend
  INTEGER :: imat, igpt
  INTEGER :: k1n1,k1n2,k1n3,k1n4,k2n1,k2n2,k2n3,k2n4
  INTEGER :: k3n1,k3n2,k3n3,k3n4


  REAL*8, DIMENSION(1:4,1:4)  :: Kloc
  REAL*8, DIMENSION(1:4)  :: Tloc, Rinloc

  REAL*8, DIMENSION(1:3,1:4) :: B  
  REAL*8, DIMENSION(1:3,1:3) :: KappaMatrx = RESHAPE( &
       (/1.000000000000000,0.000000000000000,0.000000000000000, &
       0.000000000000000,1.000000000000000,0.000000000000000, &
       0.000000000000000,0.000000000000000,1.000000000000000/),(/3,3/) )

  DO i = ElemStart, ElemEnd
     imat =  matcstet(i)
     
     n1 = ElConnVol(1,i)
     n2 = ElConnVol(2,i)
     n3 = ElConnVol(3,i)
     n4 = ElConnVol(4,i)
     
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

     Tloc(1) = T(n1)
     Tloc(2) = T(n2)
     Tloc(3) = T(n3)
     Tloc(4) = T(n4)

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
     
     c11 =    y24*z34 - z24*y34
     c21 = -( x24*z34 - z24*x34 )
     c31 =    x24*y34 - y24*x34
     
     Vx6 = -( x14*c11 + y14*c21 + z14*c31 )

     Vx6inv = 1.d0 / Vx6

! See the maple worksheet 'V3D4.mws' for the derivation of [B]
! NOTE: Factored for a more equivalent/compact form then maple's

     B(1,1)  = (y34*z24 - y24*z34) * Vx6inv
     B(2,1)  = (z34*x24 - z24*x34) * Vx6inv
     B(3,1)  = (x34*y24 - x24*y34) * Vx6inv
     B(1,2)  = (y13*z14 - y14*z13) * Vx6inv
     B(2,2)  = (z13*x14 - z14*x13) * Vx6inv
     B(3,2)  = (x13*y14 - x14*y13) * Vx6inv
     B(1,3)  = (y14*z12 - y12*z14) * Vx6inv
     B(2,3)  = (z14*x12 - z12*x14) * Vx6inv
     B(3,3)  = (x14*y12 - x12*y14) * Vx6inv
     B(1,4) = (y12*z13 - y13*z12) * Vx6inv
     B(2,4) = (z12*x13 - z13*x12) * Vx6inv
     B(3,4) = (x12*y13 - x13*y12) * Vx6inv

     vol = Vx6 / 6.d0

     Kloc = Kappa(imat)*MATMUL(MATMUL(TRANSPOSE(B),KappaMatrx),B)*vol


     RinLoc(:) = MATMUL(Kloc,Tloc)

     Rnet(n1)  = Rnet(n1)  - RinLoc(1)
     Rnet(n2)  = Rnet(n2)  - RinLoc(2)
     Rnet(n3)  = Rnet(n3)  - RinLoc(3)
     Rnet(n4)  = Rnet(n4)  - RinLoc(4)

  ENDDO
  RETURN
END SUBROUTINE V3D4_THERMAL

