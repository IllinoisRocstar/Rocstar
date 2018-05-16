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
!!     implicit_v3d8_mass_consistent
!!
!!  FUNCTION
!!     This subroutine constructs the global mass matrix as a 
!!     consistent matrix.  The resulting matrix is put into
!!     the variables in comp_row_global.
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
!!     ElConnVol -- The connectivity table (redundant)
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE implicit_v3d8_mass_consistent( NumEl, NumNP, NumMat, &
     coor, nodes, MatType, ri, rho ,ElConnVol )

  USE Precision
  USE implicit_global
  USE comp_row_global

  IMPLICIT NONE

  ! Input variables
  INTEGER :: NumEl, NumNP, NumMat
  INTEGER, DIMENSION(1:NumEl) :: MatType
  INTEGER, DIMENSION(1:8,1:NumEl) :: nodes
  REAL(KIND=wp), DIMENSION(1:3,1:NumNP) :: coor
  REAL(KIND=wp), DIMENSION(1:3,1:8) :: ri
  REAL(KIND=wp), DIMENSION(1:NumMat) :: rho
  INTEGER, DIMENSION(1:8,1:NumEl) :: ElConnVol

  ! Internal variables
  REAL(kind=wp), DIMENSION(1:8) :: N
  REAL(kind=wp), DIMENSION(1:8,1:3) :: dN
  REAL(kind=wp), DIMENSION(1:8,1:8) :: NtN
  INTEGER :: igpt
  INTEGER :: ngpts = 8
  INTEGER :: ielem
  INTEGER :: imat
  INTEGER :: i, j
  INTEGER :: dof1, dof2, gdof1, gdof2, ldof1, ldof2
  REAL(kind=wp) :: tempval
  INTEGER :: nnzold
  REAL(kind=wp) :: element_volume
  REAL(kind=wp), DIMENSION(1:3,1:3) :: jac, jacinv
  REAL(kind=wp) :: detj
  LOGICAL :: error
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  REAL(kind=wp),DIMENSION(1:3,1:8) :: coord
  
  ! Output variables
  INTEGER, ALLOCATABLE, DIMENSION(:) :: rp, cval
  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: aval
  INTEGER :: nnz

  ! Initialize the output variables
  ALLOCATE(rp(1:3*GNumNp+1))
  ALLOCATE(cval(1:1))
  ALLOCATE(aval(1:1))
  nnz = 1
  rp(1) = 0
  rp(2:3*GNumNp+1) = 1
  cval(1) = 0
  aval(1) = 0.0

  ! Loop through the elements
  DO ielem = 1, NumEl

     ! Node positions
     n1 = nodes(1,ielem)
     n2 = nodes(2,ielem)
     n3 = nodes(3,ielem)
     n4 = nodes(4,ielem)
     n5 = nodes(5,ielem)
     n6 = nodes(6,ielem)
     n7 = nodes(7,ielem)
     n8 = nodes(8,ielem)
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

     ! Find which material this element is made of
     imat = MatType(ielem)

     ! Initialize stuff
     NtN(:,:) = 0.0
     element_volume = 0.0

     ! Loop throught the gauss points in this element
     DO igpt = 1, ngpts

        ! Find the shape functions
        N(:) = 0.0
        dN(:,:) = 0.0
        CALL get_shape(ri,N,dN,igpt)

        ! Find the volume of the element
        jac(1:3,1:3) = 0.0
        jacinv(1:3,1:3) = 0.0
        CALL get_jacobien(coord,3,8,dN,jac,jacinv,detj,error)
        element_volume = element_volume + detj ! * w

        ! Integrate and assemble into NtN
        DO i = 1, 8
           DO j = 1, 8
              NtN(i,j) = NtN(i,j) + rho(imat)*N(i)*N(j)*detj  ! * w
           ENDDO
        ENDDO

     ENDDO

     ! Put these parts of the mass matrix into the global mass matrix
     DO dof1 = 1, 3  ! dof1 loop
        DO i = 1, 8  ! node1 loop
           gdof1 = 3 * Local2Global(ElConnVol(i,ielem)) - 3 + dof1
           ldof1 = i
           DO j = 1, 8  ! node2 loop
              ldof2 = j
              gdof2 = 3 * Local2Global(ElConnVol(j,ielem)) - 3 + dof1
              
              ! Make sure the matrix is exactly symmetric
              IF ( ldof1 <= ldof2 ) THEN
                 tempval = NtN(ldof1,ldof2)
              ELSE
                 tempval = NtN(ldof2,ldof1)
              ENDIF

              IF ( tempval == 0.0 ) THEN
                 print*,'ZERO DETECTED IN LOCAL MASS MATRIX!',gdof1,gdof2,myid
              END IF
              
              ! Add the value to the global matrix
              nnzold = nnz
              CALL comp_row_addval(3*GNumNp,3*GNumNp,nnz,1,rp,cval,aval,gdof1,gdof2,tempval)
              nnz = nnz_temp
              IF (nnzold /= nnz) THEN
                 DEALLOCATE(cval)
                 DEALLOCATE(aval)
                 ALLOCATE(cval(1:nnz))
                 ALLOCATE(aval(1:nnz))
                 rp = rp_temp
                 cval = cval_temp
              ENDIF
              aval = aval_temp
              DEALLOCATE(rp_temp)
              DEALLOCATE(cval_temp)
              DEALLOCATE(aval_temp)
              
           ENDDO
        ENDDO
     ENDDO

  ENDDO

  ! Put the output variables into global variables
  nnz_temp = nnz
  ALLOCATE(rp_temp(1:3*GNumNp+1))
  ALLOCATE(cval_temp(1:nnz_temp))
  ALLOCATE(aval_temp(1:nnz_temp))
  rp_temp = rp
  cval_temp = cval
  aval_temp = aval
  DEALLOCATE(rp)
  DEALLOCATE(cval)
  DEALLOCATE(aval)


END SUBROUTINE implicit_v3d8_mass_consistent

