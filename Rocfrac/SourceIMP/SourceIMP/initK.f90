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
!!     numKnnz
!!
!!  FUNCTION
!!     This subroutine uses the node numbering combined with
!!     the connectivity table in order to calculate the number
!!     of nonzero entries in the stiffness matrix due only to
!!     the elements and nodes on this processor.
!!
!!  INPUTS
!!     NumNp -- Total number of nodes that this proc knows about
!!     NumEl -- Total number of elements that this proc knows about
!!     ElConnVol -- Connectivity table
!!
!!  OUTPUTS
!!     nnz -- The number of nonzeros in the K matrix
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE numKnnz(NumNp,NumEl,ElConnVol,nnz)

  USE Precision
  USE implicit_global

  IMPLICIT NONE

  ! Input variables
  INTEGER :: NumNp, NumEl
  INTEGER :: ElConnVol(1:8,1:NumEl)

  ! Output variables
  INTEGER :: nnz

  ! Internal variables
  INTEGER :: i, j, k, m, n, counter
  INTEGER :: innz(1:56), jnnz(1:26)


!
!   Count the nonzeros that should be in the K matrix
!

  ! Loop through nodes counting off-diagonal blocks of nonzeros
  nnz = 0
  DO i = 1, NumNp

     ! Initialize connectivity variables
     innz(:) = -1
     jnnz(:) = -1
     
     ! Construct nodal connectivity vector
     DO j = 1, NumEl
        DO k = 1, 8
           IF (ElConnVol(k,j) == i) THEN
              DO m = 1, 8
                 IF (k /= m) THEN
                    DO n = 1, 56
                       IF (innz(n) == -1) THEN
                          innz(n) = ElConnVol(m,j)
                          EXIT
                       ENDIF
                    ENDDO
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     
     ! Eliminate duplicates from nodal connectivity vector
     DO j = 1, 56
        IF (innz(j) /= -1) THEN
           DO k = 1, 26
              IF (jnnz(k) /= -1) THEN
                 IF (jnnz(k) == innz(j)) THEN
                    EXIT
                 ENDIF
              ELSE
                 jnnz(k) = innz(j)
                 nnz = nnz + 1
                 EXIT
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     
  ENDDO

  ! Take into account the blocks on the diagonal
  nnz = nnz + NumNp

  ! Convert number of blocks to number of nonzeros
  nnz = 9*nnz  ! since 3 dof

  
END SUBROUTINE numKnnz













!!****
!!
!!  NAME
!!     initK
!!
!!  FUNCTION
!!     This subroutine uses the node numbering combined with
!!     the connectivity table and number of nonzeros in the K
!!     matrix to assemble the pre-allocated compressed row
!!     storage arrays.  It does this by putting a definite
!!     (not assumed) zero at each point that can contain a
!!     nonzero value.  The actual values in the K matrix are
!!     added in the v3d8_me_K subroutine.
!!
!!  INPUTS
!!     NumNp -- Total number of nodes that this proc knows about
!!     NumEl -- Total number of elements that this proc knows about
!!     ElConnVol -- Connectivity table
!!     nnz -- The number of nonzeros in the K matrix
!!
!!  OUTPUTS
!!     rp -- The row mapping vector
!!     cval -- The collumn mapping vector
!!     aval -- The value vector (which will be returned as all zeros)
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE initK(NumNp,NumEl,ElConnVol,nnz,rp,cval,aval)


  USE Precision
  USE implicit_global

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  ! Input variables
  INTEGER :: NumNp, NumEl
  INTEGER :: ElConnVol(1:8,1:NumEl)
  INTEGER :: nnz

  ! Output variables
  INTEGER :: rp(1:3*GNumNp+1)
  INTEGER :: cval(1:nnz)
  INTEGER :: aval(1:nnz)

  ! Internal variables
  INTEGER :: i, j, k, m, n, counter, ii
  INTEGER :: innz(1:56), jnnz(1:27)


  ! Initialize CRS variables
  aval(:) = 0.0
  cval(:) = -1
  rp(:) = -1
  rp(1) = 0

  ! Loop through nodes
  counter = 0
  DO i = 1, GNumNp

     ! If this node is on this processor
     IF (Global2Local(i) /= -1) THEN

        ii = Global2Local(i)

        ! Initialize connectivity variables
        innz(:) = -1
        jnnz(:) = -1
        
        ! Make sure diagonal terms are included
        jnnz(27) = ii
        
        ! Construct nodal connectivity vector
        DO j = 1, NumEl
           DO k = 1, 8
              IF (ElConnVol(k,j) == ii) THEN
                 DO m = 1, 8
                    IF (k /= m) THEN
                       DO n = 1, 56
                          IF (innz(n) == -1) THEN
                             innz(n) = ElConnVol(m,j)
                             EXIT
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
        ENDDO

        ! Eliminate duplicates from nodal connectivity vector
        DO j = 1, 56
           IF (innz(j) /= -1) THEN
              DO k = 1, 26
                 IF (jnnz(k) /= -1) THEN
                    IF (jnnz(k) == innz(j)) THEN
                       EXIT
                    ENDIF
                 ELSE
                    jnnz(k) = innz(j)
                    EXIT
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
        
        ! Loop through rows with same nonzero blocks
        DO j = 1, 3

           ! Construct cval
           DO n = 1, GNumNp
              DO k = 1, 27
                 IF(jnnz(k) /= -1) THEN
                    m = Local2Global(jnnz(k))
                    IF (m == n) THEN
                       counter = counter + 3
                       cval(counter-2) = 3*m-3
                       cval(counter-1) = 3*m-2
                       cval(counter)   = 3*m-1
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO

           ! Construct rp
           rp(3*i+j-2) = counter

        ENDDO

     ! If this node is not on this processor
     ELSE

        DO j = 1, 3
           rp(3*i+j-2) = rp(3*i+j-3)
        ENDDO

     ENDIF

        
     
  ENDDO


END SUBROUTINE initK

