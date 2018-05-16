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
!!     comp_row_add
!!
!!  FUNCTION
!!     This subroutine adds 2 matrices stored in compressed
!!     row format and assigns them to the global variables
!!     located in comp_row_global.
!
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global square matrices
!!     gndim -- Currently the same as ndim
!!     nrows1 -- The number of rows in matrix1 assigned to this processor
!!     nrows2 -- The number of rows in matrix2 assigned to this processor (MUST be equal to nrows1 to work properly)
!!     nnz1 -- The number of nonzeros in this processors section of matrix1
!!     nnz2 -- The number of nonzeros in this processors section of matrix2
!!     nstart1 -- The global number of the first row in matrix1 assigned to this processor
!!     rp1 -- The row mapping vector for matrix1
!!     cval1 -- The collumn mapping vector for matrix1
!!     aval1 -- The nonzero values in matrix1
!!     nstart2 -- The global number of the first row in matrix2 assigned to this processor
!!     rp2 -- The row mapping vector for matrix2
!!     cval2 -- The collumn mapping vector for matrix2
!!     aval2 -- The nonzero values in matrix2
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE comp_row_add(ndim,gndim,nrows1,nrows2,nnz1,nnz2, &
     nstart1,rp1,cval1,aval1,nstart2,rp2,cval2,aval2)

  USE Precision
  USE comp_row_global
 
  IMPLICIT NONE

  ! Input variables
  INTEGER :: ndim, gndim, nnz1, nnz2, nstart1, nstart2, nrows1, nrows2
  REAL(kind=wp) :: a0
  REAL(kind=wp), DIMENSION(nnz1) :: aval1
  REAL(kind=wp), DIMENSION(nnz2) :: aval2
  INTEGER, DIMENSION(nnz1) :: cval1
  INTEGER, DIMENSION(nnz2) :: cval2
  INTEGER, DIMENSION(nrows1+1) :: rp1
  INTEGER, DIMENSION(nrows2+1) :: rp2

  ! Internal variables
  INTEGER :: i, j, counter1, counter2, counter3, nnzrow, nnz, foundnz
  REAL(kind=wp), DIMENSION(gndim) :: row1, row2, row3
  
  ! Count nonzeros
  counter1 = 0
  counter2 = 0
  counter3 = 0
  DO i = 1, nrows1 ! since nrows1 = nrows2
     ! Construct the ith row in matrix 1
     row1(1:gndim) = 0.0
     DO j = 1, nnz1  ! Loop through collumns
        !IF(i >= nstart1) THEN
           IF(j > rp1(i-nstart1+1)) THEN
              IF(j <= rp1(i-nstart1+2)) THEN
                 counter1 = counter1 + 1
                 row1(cval1(counter1)+1) = aval1(counter1)
              ENDIF
           ENDIF
        !ENDIF
     ENDDO
     ! Construct the ith row in matrix 2
     row2(1:gndim) = 0.0
     DO j = 1, nnz2  ! Loop through collumns
        !IF(i >= nstart2) THEN
           IF(j > rp2(i-nstart2+1)) THEN
              IF(j <= rp2(i-nstart2+2)) THEN
                 counter2 = counter2 + 1
                 row2(cval2(counter2)+1) = aval2(counter2)
              ENDIF
           ENDIF
        !ENDIF
     ENDDO
     ! Add the rows together
     DO j = 1, gndim
        row3(j) = row1(j) + row2(j)
     ENDDO
     ! Count the nonzeros in the new row
     DO j = 1, gndim
        IF(ABS(row3(j)) > 0.0 ) THEN
           counter3 = counter3 + 1
           !print*,'nonzero found at:  ',i,j
        ENDIF
     ENDDO 
     !WRITE(*,1000) (row3(j),j=1,ndim)
  ENDDO
  nnz_temp = counter3

  ! Allocate variables
  ALLOCATE(aval_temp(nnz_temp))
  ALLOCATE(cval_temp(nnz_temp))
  ALLOCATE(rp_temp(nrows1+1))  ! since nrows1 = nrows2

  ! Construct the new matrix
  counter1 = 0
  counter2 = 0
  counter3 = 0
  rp_temp(1) = 0
  DO i = 1, nrows1
     ! Construct the ith row in matrix 1
     row1(1:gndim) = 0.0
     DO j = 1, nnz1  ! Loop through collums
        !print*,j,rp1(i-nstart1+1),rp1(i-nstart1+2),counter1,aval1(counter1),cval1(counter1)
        !IF(i >= nstart1) THEN
           IF(j > rp1(i-nstart1+1)) THEN
              IF(j <= rp1(i-nstart1+2)) THEN
                 counter1 = counter1 + 1
                 row1(cval1(counter1)+1) = aval1(counter1)
                 !print*,'FOUND ONE:  row:  ',i,'  row1(',j,') = ',row1(j)
              ENDIF
           ENDIF
        !ENDIF
     ENDDO
     ! Construct the ith row in matrix 2
     row2(1:gndim) = 0.0
     DO j = 1, nnz2  ! Loop through collums
        !IF(i >= nstart2) THEN
           !print*,j,rp2(i-nstart2+1),rp2(i-nstart2+2),counter2,aval2(counter2)
           IF(j > rp2(i-nstart2+1)) THEN
              IF(j <= rp2(i-nstart2+2)) THEN
                 counter2 = counter2 + 1
                 row2(cval2(counter2)+1) = aval2(counter2)
                 !print*,'FOUND ONE:  row:  ',i,'  row2(',j,') = ',row2(j)
              ENDIF
           ENDIF
        !ENDIF
     ENDDO
     ! Add the rows together
     DO j = 1, gndim
        row3(j) = row1(j) + row2(j)
     ENDDO
     ! Put the new row into compressed row format
     nnzrow = 0
     DO j = 1, gndim
        IF(ABS(row3(j)) > 0.0 ) THEN
           counter3 = counter3 + 1
           nnzrow = nnzrow + 1
           cval_temp(counter3) = j - 1
           aval_temp(counter3) = row3(j)
        ENDIF
     ENDDO 
     rp_temp(i+1) = rp_temp(i) + nnzrow    

!!$     if(i==ndim)then
!!$        print*,'row1:'
!!$        WRITE(*,1000) (row1(j),j=1,ndim)
!!$        print*,'row2:'
!!$        WRITE(*,1000) (row2(j),j=1,ndim)
!!$        print*,'row3:'
!!$        WRITE(*,1000) (row3(j),j=1,ndim)
!!$     endif
  ENDDO

  1000 FORMAT(1000(f5.1,' '))

END SUBROUTINE comp_row_add



!!****
!!
!!  NAME
!!     comp_row_vecmult
!!
!!  FUNCTION
!!     This subroutine multiplies a matrix in compressed
!!     row format with an explicitly defined vector and
!!     returns the vector created by this multiplication.
!!
!!  INPUTS
!!     gndim -- The size of one dimension of the global square matrices
!!     nrows1 -- The number of rows in matrix1 assigned to this processor
!!     nnz1 -- The number of nonzeros in this processors section of matrix1
!!     nstart1 -- The global number of the first row in matrix1 assigned to this processor
!!     rp1 -- The row mapping vector for matrix1
!!     cval1 -- The collumn mapping vector for matrix1
!!     aval1 -- The nonzero values in matrix1
!!     vec -- The vector that the matrix will multiply
!!
!!  OUTPUTS
!!     ans -- The solution vector
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE comp_row_vecmult(gndim,nrows1,nnz1, &
     nstart1,rp1,cval1,aval1,vec,ans)

  USE Precision
 
  IMPLICIT NONE

  ! Input variables
  INTEGER :: ndim, gndim, nnz1, nstart1, nrows1
  REAL(kind=wp) :: a0
  REAL(kind=wp), DIMENSION(nnz1) :: aval1
  INTEGER, DIMENSION(nnz1) :: cval1
  INTEGER, DIMENSION(nrows1+1) :: rp1
  REAL(kind=wp) :: vec(gndim)

  ! Output variables
  REAL(kind=wp) :: ans(nrows1)

  ! Internal variables
  INTEGER :: i, j, counter1
  REAL(kind=wp), DIMENSION(gndim) :: row1
  REAL(kind=wp) :: tempval

  ! Construct the new vector
  counter1 = 0
  ans(1:nrows1) = 0.0
  DO i = 1, gndim
     IF ((i >= nstart1) .AND. (i < nstart1+nrows1)) THEN
        row1(1:gndim) = 0.0
        ! Construct the ith row in matrix 1
        DO j = 1, nnz1  ! Loop through collums
           IF(j > rp1(i-nstart1+1)) THEN
              IF(j <= rp1(i-nstart1+2)) THEN
                 counter1 = counter1 + 1
                 row1(cval1(counter1)+1) = aval1(counter1)
              ENDIF
           ENDIF
        ENDDO
        ! Take the dot product of the row with the vector
        DO j = 1, gndim
           ans(i-nstart1+1) = ans(i-nstart1+1) + row1(j) * vec(j) 
        ENDDO
     ENDIF
  END DO

  1000 FORMAT(1000(f5.1,' '))

END SUBROUTINE comp_row_vecmult




!!****
!!
!!  NAME
!!     comp_row_mat
!!
!!  FUNCTION
!!     This subroutine takes a matrix and puts it
!!     into compressed row storage.  The resulting
!!     arrays are assigned to the global variables
!!     located in comp_row_global.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global square matrices
!!     A -- The matrix
!!     nstart -- The first row in the matrix that is assigned to this processor
!!     nrows -- The number of rows in the matrix assigned to this processor
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE comp_row_mat(ndim,A,nstart,nrows)

  USE comp_row_global
  USE Precision

  ! Input variables
  INTEGER :: ndim, nstart, nrows
  REAL(kind=wp), DIMENSION(1:ndim,1:ndim) :: A

  ! Internal variables
  INTEGER :: i, j
  INTEGER :: avali
  INTEGER :: counter

  ! Count the nonzeros in the matrix
  nnz_temp = 0
  DO i = nstart, nstart+nrows-1
     DO j = 1, ndim
        IF ( ABS(A(i,j)) > 0.0 ) THEN
           nnz_temp = nnz_temp + 1
        ENDIF
     ENDDO
  ENDDO
  ALLOCATE(rp_temp(1:nrows+1))
  ALLOCATE(cval_temp(1:nnz_temp))
  ALLOCATE(aval_temp(1:nnz_temp))

  ! Construct the matrix in compressed row format
  avali = 1
  rp_temp(1) = 0
  DO i = nstart, nstart+nrows-1
     DO j = 1, ndim
        IF( ABS(A(i,j)) > 0.0 ) THEN
           aval_temp(avali) = A(i,j)
           cval_temp(avali) = j-1
           avali = avali + 1
        ENDIF
     ENDDO
     rp_temp(i-nstart+2) = avali-1
  ENDDO

END SUBROUTINE comp_row_mat





!!****
!!
!!  NAME
!!     comp_row_addval
!!
!!  FUNCTION
!!     This subroutine adds a value to a matrix stored in
!!     compressed row format.  If a value exists in the
!!     position for the new value to be added, they are 
!!     summed.  The modified matrix is then stored in the 
!!     global variables located in implicit_global.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global square matrices
!!     nrows -- The number of rows in matrix1 assigned to this processor
!!     nnz -- The number of nonzeros in this processors section of matrix1
!!     nstart -- The global number of the first row in matrix1 assigned to this processor
!!     rp -- The row mapping vector for matrix1
!!     cval -- The collumn mapping vector for matrix1
!!     aval -- The nonzero values in matrix1
!!     ipos -- The global i index of the position to add the value within the matrix
!!     jpos -- The global j index of the position to add the value within the matrix
!!     val -- The value to add to the matrix
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE comp_row_addval(ndim,nrows,nnz,nstart,rp,cval,aval,ipos,jpos,val)

  USE Precision
  USE comp_row_global
 
  IMPLICIT NONE

  ! Input variables
  INTEGER :: ndim, nnz, nstart, nrows, ipos, jpos
  REAL(kind=wp) :: val, tempval, tempval2
  REAL(kind=wp), DIMENSION(nnz) :: aval
  INTEGER, DIMENSION(nnz) :: cval
  INTEGER, DIMENSION(nrows+1) :: rp

  ! Internal variables
  INTEGER :: i, j, counter, counter1, nnzrow
  REAL(kind=wp), DIMENSION(ndim) :: row

  ! Count nonzeros in new matrix
  counter = 0
  DO i = rp(ipos-nstart+1)+1, rp(ipos-nstart+2)
     IF (cval(i) == jpos-1) THEN  ! value already there
        counter = 1
        EXIT
     ENDIF
  ENDDO
  IF (counter == 1) THEN
     nnz_temp = nnz  ! there was a value there already
  ELSE
     nnz_temp = nnz + 1  ! there was not a value there already
  ENDIF
  
  ! Allocate variables
  ALLOCATE(aval_temp(nnz_temp))
  ALLOCATE(cval_temp(nnz_temp))
  ALLOCATE(rp_temp(nrows+1))
  
  ! Construct the new matrix
  IF (nnz_temp == nnz) THEN
     rp_temp(1:nrows+1) = rp(1:nrows+1)
     cval_temp(1:nnz_temp) = cval(1:nnz)
     aval_temp(1:nnz_temp) = aval(1:nnz)
     DO i = rp(ipos-nstart+1)+1, rp(ipos-nstart+2)
        IF (cval(i) == jpos-1) THEN
           aval_temp(i) = aval_temp(i) + val
        ENDIF
     ENDDO
  ELSE
     DO i = 1, rp(ipos-nstart+2)
        cval_temp(i) = cval(i)
        aval_temp(i) = aval(i)
     ENDDO
     cval_temp(rp(ipos-nstart+2)+1) = jpos-1
     aval_temp(rp(ipos-nstart+2)+1) = val
     DO i = rp(ipos-nstart+2)+1, rp(nrows+1)
        cval_temp(i+1) = cval(i)
        aval_temp(i+1) = aval(i)
     ENDDO
     DO i = 1, ipos-nstart+1
        rp_temp(i) = rp(i)
     ENDDO
     DO i = ipos-nstart+2, nrows+1
        rp_temp(i) = rp(i) + 1
     ENDDO
  ENDIF


END SUBROUTINE comp_row_addval




!!****
!!
!!  NAME
!!     comp_row_getval
!!
!!  FUNCTION
!!     This subroutine returns a value within a matrix
!!     stored in compressed row format when given the
!!     i and j indices of the value.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global square matrices
!!     nrows -- The number of rows in matrix1 assigned to this processor
!!     nnz -- The number of nonzeros in this processors section of matrix1
!!     nstart -- The global number of the first row in matrix1 assigned to this processor
!!     rp -- The row mapping vector for matrix1
!!     cval -- The collumn mapping vector for matrix1
!!     aval -- The nonzero values in matrix1
!!     ipos -- The global i index of the value
!!     jpos -- The global j index of the value
!!
!!  OUTPUTS
!!     val -- The value at (ipos,jpos) in the matrix
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE comp_row_getval(ndim,nrows,nnz,nstart,rp,cval,aval,ipos,jpos,val)

  USE Precision
 
  IMPLICIT NONE

  ! Input variables
  INTEGER :: ndim, nnz, nstart, nrows, ipos, jpos
  REAL(kind=wp) :: val
  REAL(kind=wp), DIMENSION(nnz) :: aval
  INTEGER, DIMENSION(nnz) :: cval
  INTEGER, DIMENSION(nrows+1) :: rp

  ! Internal variables
  INTEGER :: i, j, counter
  REAL(kind=wp), DIMENSION(ndim) :: row

!!$  ! Count nonzeros in new matrix
!!$  counter = 0
!!$  DO i = 1, nrows
!!$     row(1:ndim) = 0.0
!!$     DO j = 1, nnz 
!!$        IF(j > rp(i)) THEN
!!$           IF(j <= rp(i+1)) THEN
!!$              counter = counter + 1
!!$              row(cval(counter)+1) = aval(counter)
!!$           ENDIF
!!$        ENDIF
!!$     ENDDO
!!$     IF (i == ipos-nstart+1) THEN
!!$        val = row(jpos)
!!$        EXIT
!!$     ENDIF
!!$  ENDDO

  val = 0.0
  DO i = rp(ipos-nstart+1)+1, rp(ipos-nstart+2)
     IF (cval(i) == jpos-1) THEN  ! value already there
        val = aval(i)
        EXIT
     ENDIF
  ENDDO

  !IF(val /= 0.0) print*,'COMP_ROW_GETVAL returning nonzero at ',ipos,jpos,' of value ',val

END SUBROUTINE comp_row_getval





!!****
!!
!!  NAME
!!     comp_row_resize
!!
!!  FUNCTION
!!     This subroutine resizes the rows of a matrix
!!     in compressed row storage and returns the new
!!     matrix in the global variables in
!!     comp_row_global
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global square matrices
!!     nrows -- The number of rows in matrix1 assigned to this processor
!!     nnz -- The number of nonzeros in this processors section of matrix1
!!     nstart -- The global number of the first row in matrix1 assigned to this processor
!!     rp -- The row mapping vector for matrix1
!!     cval -- The collumn mapping vector for matrix1
!!     aval -- The nonzero values in matrix1
!!     newnstart -- The new starting row for this processor
!!     newnrows -- The new number of rows for this processor
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE comp_row_resize(ndim,nrows,nnz,nstart,rp,cval,aval,newnstart,newnrows)


  USE Precision
  USE comp_row_global
 
  IMPLICIT NONE

  ! Input variables
  INTEGER :: ndim, nnz, nstart, nrows, ipos, jpos, newnstart, newnrows
  REAL(kind=wp) :: val
  REAL(kind=wp), DIMENSION(nnz) :: aval
  INTEGER, DIMENSION(nnz) :: cval
  INTEGER, DIMENSION(nrows+1) :: rp

  ! Internal variables
  INTEGER :: i, j, counter, counter1, counter2, nnzrow
  REAL(kind=wp), DIMENSION(ndim) :: row

  ! Count nonzeros in new matrix
  counter = 0
  nnz_temp = 0
  DO i = 1, nrows
     row(1:ndim) = 0.0
     DO j = 1, nnz
        IF(j > rp(i)) THEN
           IF(j <= rp(i+1)) THEN
              counter = counter + 1
              IF(i >= newnstart-nstart+1) THEN
                 IF(i < newnstart -nstart + 1 + newnrows) THEN
                    nnz_temp = nnz_temp + 1
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  
  !print*,'nnz = ',nnz_temp

  ! Allocate variables
  ALLOCATE(aval_temp(nnz_temp))
  ALLOCATE(cval_temp(nnz_temp))
  ALLOCATE(rp_temp(newnrows+1))

  ! Construct the new matrix
  counter = 0
  counter1 = 0
  counter2 = 0
  rp_temp(1) = 0
  DO i = 1, nrows
     ! Construct the ith row in matrix 1
     row(1:ndim) = 0.0
     DO j = 1, nnz  ! Loop through collums
        IF(j > rp(i)) THEN
           IF(j <= rp(i+1)) THEN
              counter = counter + 1
              IF(i >= newnstart - nstart + 1) THEN
                 IF(i < newnstart + newnrows - nstart + 1) THEN
                    row(cval(counter)+1) = aval(counter)
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     ! Put the new row into compressed row format
     nnzrow = 0
     DO j = 1, ndim
        IF(row(j) /= 0.0) THEN
           IF(i >= newnstart - nstart + 1) THEN
              IF(i < newnstart + newnrows - nstart + 1) THEN
                 counter1 = counter1 + 1
                 nnzrow = nnzrow + 1
                 cval_temp(counter1) = j - 1
                 aval_temp(counter1) = row(j)
              ENDIF
           ENDIF
        ENDIF
     ENDDO 
     IF(i >= newnstart - nstart + 1) THEN
        IF(i < newnstart + newnrows - nstart + 1) THEN
           counter2 = counter2 + 1
           rp_temp(counter2+1) = rp_temp(counter2) + nnzrow
        ENDIF
     ENDIF
  ENDDO


END SUBROUTINE comp_row_resize




