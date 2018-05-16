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
SUBROUTINE v3d8_mass( coor,nodes,MatType,rho,xm, &
     NumNP,NumEl,NumMat,nstart, nend)

!!****f* Rocfrac/Rocfrac/Source/v3d8_mass.f90
!!
!!  NAME
!!    v3d8_mass
!!
!!  FUNCTION
!!
!!    Computes the lumped mass matrix for a 8-node hexahedral.
!!
!!  INPUTS
!!
!!   NumNP   -- Number of nodes
!!   NumEL   -- Number of elements
!!   Coor    -- number of coordinates
!!   MatType -- Material id
!!   rho     -- Density
!!   nodes   -- Nodal connectivity
!!   nstart, nend -- element beginning and end loop counter
!!   NumMat -- number of materials
!!
!!  OUTPUT
!!
!!    xm -- lumped mass matrix
!!
!!****



  IMPLICIT NONE

  INTEGER :: NumEl, NumNP, NumMat
  INTEGER, DIMENSION(1:NumEl) :: MatType
  INTEGER, DIMENSION(1:8,1:NumEl) :: nodes
  REAL*8, DIMENSION(1:3,1:NumNP) :: coor
  REAL*8, DIMENSION(1:NumMat) :: rho
  REAL*8, DIMENSION(1:NumNP) :: xm

! root3 = 1./sqrt(3.)

  REAL*8, DIMENSION(1:3,1:8) :: ri = RESHAPE( &
       (/-0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626, 0.577350269189626, 0.577350269189626, &
         -0.577350269189626, 0.577350269189626, 0.577350269189626/),(/3,8/) )

  REAL*8 :: nNn(8), dn(8,3)
  integer :: ngpts = 8
  INTEGER :: i,imat, igpt
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  REAL*8 :: element_volume
  REAL*8,DIMENSION(1:3,1:8) :: coord
  REAL*8,DIMENSION(1:3,1:3) :: jac, jacinv
  REAL*8 :: detj
  REAL*8 :: node_mass
  LOGICAL :: error
  integer :: nstart, nend

  error = .FALSE.
  

  DO i = nstart, nend

     imat = MatType(i)
     
     element_volume = 0.  ! Initialize

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
     
        CALL get_jacobien(coord,3,8,dn, &  ! Compute jacobien 
             jac,jacinv,detj,error)

        IF(error) THEN                           ! Check for singular jacobien
           WRITE(*,1000) igpt, i, detj
           STOP
        END IF

        element_volume = element_volume + detj ! * wi(igpt), weigth for hex element = 1 
     ENDDO

     node_mass = element_volume*rho(imat)*0.125

     xm(n1  ) = xm(n1  ) + node_mass
     xm(n2  ) = xm(n2  ) + node_mass
     xm(n3  ) = xm(n3  ) + node_mass
     xm(n4  ) = xm(n4  ) + node_mass
     xm(n5  ) = xm(n5  ) + node_mass
     xm(n6  ) = xm(n6  ) + node_mass
     xm(n7  ) = xm(n7  ) + node_mass
     xm(n8  ) = xm(n8  ) + node_mass
      
  ENDDO

!  xm(:) = 1.d0/xm(:)

1000 FORMAT (/,2x,'>>> Fatal error!', &
          /,6x,'Jacobien at gauss point (',i2, &
          ') of element ',i6,' is singular with detj =',e15.8)

END SUBROUTINE v3d8_mass

