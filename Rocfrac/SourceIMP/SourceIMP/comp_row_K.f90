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
!!     comp_row_K
!!
!!  FUNCTION
!!     This subroutine takes the K matrix and puts it
!!     into compressed row storage.  It will always 
!!     assure that the matrix is symmetric. The resulting
!!     arrays are assigned to the global variables
!!     located in comp_row_global.
!!
!!  INPUTS
!!     ndim -- Dimension of one size of the square stiffness matrix
!!     A -- The part of the stiffness matrix assigned to this processor
!!     nstart -- The first row of the stiffness matrix assigned to this processor
!!     nrows -- The number of rows of the stiffness matrix assigned to this processor
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE comp_row_K(ndim,A,nstart,nrows)

  USE comp_row_global
  USE implicit_global
  USE Precision

  ! Input variables
  INTEGER :: ndim, nstart, nrows
  REAL(kind=wp), DIMENSION(1:ndim,1:ndim) :: A

  ! Internal variables
  INTEGER :: i, j, inode, jnode, idof, jdof
  INTEGER :: avali
  INTEGER :: counter

  ! Count the nonzeros in the matrix
  nnz_temp = 0
  DO i = nstart, nstart+nrows-1
     DO j = 1, 3*GNumNp
        inode = INT((i-0.5)/3.0)+1
        jnode = INT((j-0.5)/3.0)+1
        idof = i - 3*inode + 3
        jdof = j - 3*jnode + 3

         IF ((Global2Local(inode) /= -1))THEN
           IF ((Global2Local(jnode) /= -1))THEN
              !print*,3*Global2Local(inode)+idof-3,3*Global2Local(jnode)+jdof-3
              IF ( 3*Global2Local(inode)+idof-3 > 3*Global2Local(jnode)+jdof-3 ) THEN  ! something in lower triangle
                 !print*,A(3*Global2Local(jnode)+jdof-3,3*Global2Local(inode)+idof-3)
                 IF ( ABS(A(3*Global2Local(jnode)+jdof-3,3*Global2Local(inode)+idof-3)) /= 0.0 ) THEN
                    nnz_temp = nnz_temp + 1
                 ENDIF
              ELSE  ! something on diagonal or upper triangle
                 IF ( ABS(A(3*Global2Local(inode)+idof-3,3*Global2Local(jnode)+jdof-3)) /= 0.0 ) THEN
                    nnz_temp = nnz_temp + 1
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  ALLOCATE(rp_temp(1:nrows+1))
  ALLOCATE(cval_temp(1:nnz_temp))
  ALLOCATE(aval_temp(1:nnz_temp))
  
  !print*,myid,' nnz = ',nnz_temp

  ! Construct the matrix in compressed row format
  avali = 1
  rp_temp(1) = 0
  DO i = nstart, nstart+nrows-1
     DO j = 1, 3*GNumNp
        inode = INT((i-0.5)/3.0)+1
        jnode = INT((j-0.5)/3.0)+1
        idof = i - 3*inode + 3
        jdof = j - 3*jnode + 3
        IF ((Global2Local(inode) /= -1))THEN
           IF ((Global2Local(jnode) /= -1))THEN
              !print*,3*Global2Local(inode)+idof-3,3*Global2Local(jnode)+jdof-3
              IF ( 3*Global2Local(inode)+idof-3 > 3*Global2Local(jnode)+jdof-3 ) THEN  ! something in lower triangle
                 !print*,A(3*Global2Local(jnode)+jdof-3,3*Global2Local(inode)+idof-3)
                 IF ( ABS(A(3*Global2Local(jnode)+jdof-3,3*Global2Local(inode)+idof-3)) /= 0.0 ) THEN
                    aval_temp(avali) = A(3*Global2Local(jnode)+jdof-3,3*Global2Local(inode)+idof-3)
                    cval_temp(avali) = j-1
                    avali = avali + 1
                 ENDIF
              ELSE  ! something on diagonal or upper triangle
                 IF ( ABS(A(3*Global2Local(inode)+idof-3,3*Global2Local(jnode)+jdof-3)) /= 0.0 ) THEN
                    aval_temp(avali) = A(3*Global2Local(inode)+idof-3,3*Global2Local(jnode)+jdof-3)
                    cval_temp(avali) = j-1
                    avali = avali + 1
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     rp_temp(i-nstart+2) = avali-1
  ENDDO

END SUBROUTINE comp_row_K

