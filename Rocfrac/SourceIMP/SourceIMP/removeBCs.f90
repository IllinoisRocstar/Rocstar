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
!!     removeBCs_meff
!!
!!  FUNCTION
!!     This subroutine removes rows and collumns from the effective
!!     mass matrix that are associated with perscribed boundary
!!     conditions.  The new effective mass matrix is then put into
!!     the global variables in comp_row_global.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global Meff matrix
!!     nrows -- The number of rows assigned to this processor
!!     nnz -- The number of nonzeros in section of the Meff matrix on this processor
!!     nstart -- The global index of the first row assigned to this processor
!!     rp1 -- The row mapping vector
!!     cval -- The collumn mapping vector
!!     aval -- The nonzero value vector
!!
!!  OUTPUTS
!!     newnrows -- The number of rows assigned to this proc after BCs have been removed
!!     newnstart -- The global index of the first row assigned to this proc after BCs have been removed
!!     newndim -- The size of the global Meff matrix after BCs have been removed
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE removeBCs_meff(ndim,nrows,nnz,nstart,rp1,cval,aval,newnrows,newnstart,newndim,global)

  USE removeBCs_global
  USE comp_row_global
  USE implicit_global
  USE Precision
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL) :: global

  INTEGER :: i, j, m, counter1, counter2, counter3, counter4, counter5

  INTEGER :: ndim, nnz, nstart, nrows
  REAL(kind=wp), DIMENSION(nnz) :: aval
  INTEGER, DIMENSION(nnz) :: cval
  INTEGER, DIMENSION(nrows+1) :: rp1

  INTEGER :: numdisp  ! Local number of displacement BCs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dispbc  ! Local displacement BCs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: tempintv
  INTEGER :: newnrows
  INTEGER :: newnstart
  INTEGER :: newndim


  ! Count the number of displacement BCs locally
  numdisp = 0
  DO i = 1, global%NumNp
     IF (NodeProc(i) == myid) THEN
        DO j = 1, 3
           IF(node_flag(i,j) == 8) THEN  ! Imposed constant nodal displacement
              numdisp = numdisp + 1
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !print*,myid,' number of local disp BCs = ',numdisp

  ! Sum the local disp BCs, then broadcast to all procs
  CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
  CALL MPI_REDUCE(numdisp, gnumdisp, 1, MPI_INTEGER, MPI_SUM, 0 , ROCSTAR_COMMUNICATOR, ierr)
  CALL MPI_BCAST(gnumdisp, 1, MPI_INTEGER, 0, ROCSTAR_COMMUNICATOR, ierr)
  newndim = 3*GNumNp - gnumdisp
  !print*,myid,' number of global disp BCs = ',gnumdisp

  ! Figure out which DOF have displacements perscribed
  IF(numdisp > 0) THEN
     ALLOCATE(dispbc(1:numdisp))
     dispbc(:) = 0
     counter2 = 0
     DO m = 1, GNumNp
        DO i = 1, global%NumNp
           IF (Local2Global(i) == m) THEN
              IF (NodeProc(i) == myid) THEN
                 DO j = 1, 3
                    IF(node_flag(i,j) == 8) THEN  ! Imposed constant nodal displacement
                       counter2 = counter2 + 1
                       dispbc(counter2) = 3*m - 3 + j
                    ENDIF
                 ENDDO
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !print*,myid,' local disp bcs = ',dispbc(:)

  ! Have each proc in turn broadcast its number of displacement BCs to the other procs
  ALLOCATE(numdispproc(1:nprocs))
  numdispproc(:) = 0
  numdispproc(myid+1) = numdisp
  DO i = 1, nprocs
     CALL MPI_BCAST(numdispproc(i), 1, MPI_INTEGER, i-1, ROCSTAR_COMMUNICATOR, ierr)
  ENDDO
  newnstart = nstart
  DO i = 1, myid
     newnstart = newnstart - numdispproc(i)
  ENDDO
  !print*,myid,' number of disp bcs on other procs = ',numdispproc(:)

  ! Have each proc in turn broadcast the numbers of its displacement BCs
  ALLOCATE(gdispbc(1:gnumdisp))
  gdispbc(:) = 0
  counter1 = 0
  DO i = 1, nprocs
     IF (numdispproc(i) > 0) THEN
        ALLOCATE(tempintv(1:numdispproc(i)))
        tempintv(:) = 0
        if (i-1 == myid) tempintv(:) = dispbc(:)
        CALL MPI_BCAST(tempintv(1), numdispproc(i), MPI_INTEGER, i-1, ROCSTAR_COMMUNICATOR, ierr)
        DO j = 1, numdispproc(i)
           counter1 = counter1 + 1
           gdispbc(counter1) = tempintv(j)
        ENDDO
        DEALLOCATE(tempintv)
     ENDIF
  ENDDO
  !print*,myid,' global disp bcs = ',gdispbc(:)

  

  ! Count nonzeros in new matrix
  nnz_temp = 0
  counter2 = 0
  DO i = 1, nrows
     DO j = rp1(i)+1, rp1(i+1)
        counter2 = counter2 + 1
        counter1 = 0
        DO m = 1, gnumdisp
           IF((gdispbc(m) == i + nstart - 1).OR.(gdispbc(m) == cval(counter2)+1)) THEN
              counter1 = 1
           ENDIF
        ENDDO
        !if(counter1==1) print*,myid,' non-zero removed at ',i+nstart-1,cval(counter2)+1
        IF(counter1 == 0) nnz_temp = nnz_temp + 1
     ENDDO
  ENDDO

  ! Count rows in new matrix
  newnrows = 0
  counter2 = 0
  DO i = 1, nrows
     counter1 = 0
     DO m = 1, gnumdisp
        IF(gdispbc(m) == i + nstart - 1) THEN
           counter1 = 1
        ENDIF
     ENDDO
     !if(counter1==1) print*,myid,' row removed at ',i+nstart-1
     IF(counter1 == 0) newnrows = newnrows + 1
  ENDDO

  ! Allocate variables for new matrix
  ALLOCATE(rp_temp(1:newnrows+1))
  ALLOCATE(cval_temp(1:nnz_temp))
  ALLOCATE(aval_temp(1:nnz_temp))

  ! Construct new matrix
  counter2 = 0
  counter3 = 0
  counter4 = 0
  rp_temp(1) = 0
  DO i = 1, nrows
     DO j = rp1(i)+1, rp1(i+1)
        counter2 = counter2 + 1
        counter1 = 0
        counter5 = 0
        DO m = 1, gnumdisp
           IF(gdispbc(m) < cval(counter2)+1) THEN
              counter5 = counter5 + 1
           ENDIF
           IF((gdispbc(m) == i + nstart - 1).OR.(gdispbc(m) == cval(counter2)+1)) THEN
              counter1 = 1
           ENDIF
        ENDDO
        IF(counter1 == 0) THEN
           counter3 = counter3 + 1
           aval_temp(counter3) = aval(counter2)
           cval_temp(counter3) = cval(counter2) - counter5
        ENDIF
     ENDDO
     counter1 = 0
     DO m = 1, gnumdisp
        IF(gdispbc(m) == i + nstart - 1) THEN
           counter1 = 1
        ENDIF
     ENDDO
     IF(counter1 == 0) THEN
        counter4 = counter4 + 1
        rp_temp(counter4+1) = counter3
     ENDIF
  ENDDO
  

END SUBROUTINE removeBCs_meff

















!!****
!!
!!  NAME
!!     removeBCs_pbar
!!
!!  FUNCTION
!!     This subroutine removes values from the load vector
!!     that are associated with perscribed boundary
!!     conditions.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global Meff matrix
!!     nstart -- The global index of the first row assigned to this processor
!!     pbar -- The part of the global load vector that is assigned to this processor
!!     newndim -- The size of the global Meff matrix after BCs have been removed
!!
!!  OUTPUTS
!!     newpbar -- The part of the global load vector that is assigned to this processor after BCs have been removed
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE removeBCs_pbar(nstart,ndim,pbar,newndim,newpbar)

  USE removeBCs_global
  USE comp_row_global
  USE implicit_global
  USE Precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: i, m, counter1, counter2

  INTEGER :: ndim, newndim, nstart
  REAL(kind=wp), DIMENSION(ndim) :: pbar
  REAL(kind=wp), DIMENSION(newndim) :: newpbar
  

  ! Construct new vector
  counter2 = 0
  DO i = 1, ndim
     counter1 = 0
     DO m = 1, gnumdisp
        IF(gdispbc(m) == i + nstart - 1) THEN
           counter1 = 1
        ENDIF
     ENDDO
     IF(counter1 == 0) THEN
        counter2 = counter2 + 1
        newpbar(counter2) = pbar(i)
     ENDIF
  ENDDO
  

END SUBROUTINE removeBCs_pbar









!!****
!!
!!  NAME
!!     removeBCs_newa
!!
!!  FUNCTION
!!     This subroutine takes the acceleration vector that
!!     does not include boundary conditions and expands it
!!     back into the acceleration vector with rows associated
!!     with boundary conditions in it.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global Meff matrix
!!     nstart -- The global index of the first row assigned to this processor
!!     newa -- The part of the global acceleration vector assigned to this proc without rows associated with boundary conditions
!!     newndim -- The size of the global Meff matrix after BCs have been removed
!!
!!  OUTPUTS
!!     a -- The part of the global acceleration vector assigned to this proc with rows associated with boundary conditions
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE removeBCs_newa(nstart,ndim,a,newndim,newa)

  USE removeBCs_global
  USE comp_row_global
  USE implicit_global
  USE Precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: i, m, counter1, counter2, counter3

  INTEGER :: ndim, newndim, nstart
  REAL(kind=wp), DIMENSION(ndim) :: a
  REAL(kind=wp), DIMENSION(newndim) :: newa
  

  ! Construct new vector
  counter2 = 0
  DO i = 1, newndim
     counter3 = 0
     DO m = 1, gnumdisp
        IF((gdispbc(m) <= i + counter3 + nstart - 1).AND.(gdispbc(m) >= nstart)) THEN
           counter3 = counter3 + 1
        ENDIF
     ENDDO
     a(i+counter3) = newa(i)
     !if(myid==2) print*,myid,i,counter3,i+counter3
  ENDDO
  

END SUBROUTINE removeBCs_newa

