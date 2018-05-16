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

!!****
!!
!!  NAME
!!     implicit_v3d8_mass
!!
!!  FUNCTION
!!     This subroutine constructs the global mass matrix as a lumped matrix.
!!
!!  INPUTS
!!     NumEl -- The number of elements that this processor knows about
!!     NumNp -- The number of nodes that this processor knows about
!!     NumMat -- The global number of materials in the structure
!!     coor -- The coordinates of each of the nodes
!!     nodes -- The connectivity table
!!     MatType -- Mapping of material types to elements
!!     ri -- Gauss point positions within an element
!!     rho -- Density of the materials
!!     iElStart -- First element to compute the mass matrix for
!!     iElEnd -- Last element to compute the mass matrix for
!!
!!  OUTPUTS
!!     xm -- The inverse of the diagonal of the mass matrix
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE implicit_v3d8_mass( NumEl, NumNP, NumMat,&
     coor,nodes,MatType,ri,rho,xm,iElStart, iElEnd)


  USE Precision

  IMPLICIT NONE

  INTEGER :: iElStart, iElEnd

  INTEGER :: NumEl, NumNP, NumMat
  INTEGER, DIMENSION(1:NumEl) :: MatType
  INTEGER, DIMENSION(1:8,1:NumEl) :: nodes
  REAL(KIND=wp), DIMENSION(1:3,1:NumNP) :: coor
  REAL(KIND=wp), DIMENSION(1:3,1:8) :: ri
  REAL(KIND=wp), DIMENSION(1:NumMat) :: rho
  REAL(KIND=wp), DIMENSION(1:NumNP) :: xm

  REAL(kind=wp) :: nNn(8), dn(8,3)
  integer :: ngpts = 8
  INTEGER :: i,imat, igpt
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  REAL(kind=wp) :: element_volume
  REAL(kind=wp),DIMENSION(1:3,1:8) :: coord
  REAL(kind=wp),DIMENSION(1:3,1:3) :: jac, jacinv
  REAL(kind=wp) :: detj
  REAL(kind=wp) :: node_mass
  LOGICAL :: error

  INTEGER :: j, k

!!$  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: WT, HH, DH1, DH2, DH3
!!$  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: DET
!!$  REAL(KIND=wp), DIMENSION(1:3,1:280) :: DHG
!!$
!!$  REAL(KIND=wp), DIMENSION(1:8) :: NodalMass, Mass_ii
!!$  REAL(KIND=wp) :: Mass, SumMass_ii
!!$  REAL(KIND=wp) :: TotMass
!!$  INTEGER :: ii, j

!!$  ALLOCATE(WT(1:ngpts))
!!$  ALLOCATE(HH(1:ngpts*8))
!!$  ALLOCATE(DH1(1:ngpts*8),DH2(1:ngpts*8),DH3(1:ngpts*8))
!!$  ALLOCATE(DET(1:ngpts))

  error = .FALSE.

  

  DO i = iElStart, iElEnd

     imat = MatType(i)
     
     element_volume = 0._wp  ! Initialize

     n1 = nodes(1,i)
     n2 = nodes(2,i)
     n3 = nodes(3,i)
     n4 = nodes(4,i)
     n5 = nodes(5,i)
     n6 = nodes(6,i)
     n7 = nodes(7,i)
     n8 = nodes(8,i)

     coord(1,1) = coor(1,n1)
     coord(2,1) = coor(2,n1)
     coord(3,1) = coor(3,n1)

     coord(1,2) = coor(1,n2)
     coord(2,2) = coor(2,n2)
     coord(3,2) = coor(3,n2)

     coord(1,3) = coor(1,n3)
     coord(2,3) = coor(2,n3)
     coord(3,3) = coor(3,n3)

     coord(1,4) = coor(1,n4)
     coord(2,4) = coor(2,n4)
     coord(3,4) = coor(3,n4)

     coord(1,5) = coor(1,n5)
     coord(2,5) = coor(2,n5)
     coord(3,5) = coor(3,n5)

     coord(1,6) = coor(1,n6)
     coord(2,6) = coor(2,n6)
     coord(3,6) = coor(3,n6)

     coord(1,7) = coor(1,n7)
     coord(2,7) = coor(2,n7)
     coord(3,7) = coor(3,n7)

     coord(1,8) = coor(1,n8)
     coord(2,8) = coor(2,n8)
     coord(3,8) = coor(3,n8)

     DO igpt = 1, ngpts                 ! LOOP over gauss points

        CALL get_shape(ri,nNn,dn,igpt)  ! Get shape functions and derivatives

        jac(1:3,1:3) = 0.0
        jacinv(1:3,1:3) = 0.0
     
        CALL get_jacobien(coord,3,8,dn, &  ! Compute jacobien 
             jac,jacinv,detj,error)

        IF(error) THEN                           ! Check for singular jacobien
           WRITE(*,1000) igpt, i, detj
           STOP
        END IF

        element_volume = element_volume + detj ! * wi(igpt), weigth for hex element = 1 
     ENDDO


     node_mass = element_volume*rho(imat)*0.125_wp

     xm(n1  ) = xm(n1  ) + node_mass
     xm(n2  ) = xm(n2  ) + node_mass
     xm(n3  ) = xm(n3  ) + node_mass
     xm(n4  ) = xm(n4  ) + node_mass
     xm(n5  ) = xm(n5  ) + node_mass
     xm(n6  ) = xm(n6  ) + node_mass
     xm(n7  ) = xm(n7  ) + node_mass
     xm(n8  ) = xm(n8  ) + node_mass


! alternative

!!$     CALL DHH8 (ngpts,WT,HH,DH1,DH2,DH3)
!!$
!!$     CALL JAC3D (8,coord(1,:),coord(2,:),coord(3,:),ngpts, &
!!$             DH1,DH2,DH3,DH1,DH2,DH3,DHG,DET)
!!$
!!$     element_volume = 0.d0
!!$     Mass_ii(:) = 0.d0
!!$
!!$! form diagonal terms of consistent mass matrix
!!$
!!$     ! M   = int ( N  rho N  dV )
!!$     !  ii          i      i 
!!$
!!$
!!$     J = 0
!!$     DO igpt = 1, ngpts
!!$        element_volume = element_volume + det(igpt)* wt(igpt)
!!$        DO ii = 1, 8
!!$           Mass_ii(ii) = Mass_ii(ii) + HH(J+ii) *  rho(imat)* HH(J+ii) * det(igpt)* wt(igpt)
!!$        ENDDO
!!$        J = J + 8
!!$     ENDDO
!!$
!!$
!!$!
!!$!   /
!!$!   | rho Dv
!!$!   |
!!$!   /
!!$
!!$     Mass = rho(imat)*element_volume
!!$
!!$! Some the diagonal tems Mii of the consistent mass matix
!!$
!!$     SumMass_ii = 0.d0
!!$     DO ii = 1, 8
!!$        SumMass_ii = SumMass_ii + Mass_ii(ii)
!!$     ENDDO
!!$
!!$     DO ii = 1, 8
!!$        NodalMass(ii) = Mass_ii(ii) * Mass / SumMass_ii
!!$     ENDDO
!!$
!!$     PRINT*,NodalMass(:)

  ENDDO

  xm(:) = 1.0_wp/xm(:)

1000 FORMAT (/,2x,'>>> Fatal error!', &
          /,6x,'Jacobien at gauss point (',i2, &
          ') of element ',i6,' is singular with detj =',e14.8)

END SUBROUTINE implicit_v3d8_mass

