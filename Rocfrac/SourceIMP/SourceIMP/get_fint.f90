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
!!     get_fint
!!
!!  FUNCTION
!!     This subroutine constructs the internal force vector for this
!!     processor by multiplying its part of the stiffness matrix by
!!     the global displacment vector.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global stiffness matrix
!!     nrows -- The number of rows that are assigned to this processor
!!     nnz -- The number of nonzeros in this processor's part of the stiffness matrix
!!     nstart -- The global index of the first row assigned to this processor
!!     rp_k -- The row mapping vector
!!     cval_k -- The collumn mapping vector
!!     aval_k -- The nonzero values in this processor's part of the stiffness matrix
!!     dispin -- The part of the global displacement vector thats assigned to this processor
!!
!!  OUTPUTS
!!     fint -- The part of the global internal force vector thats assigned to this processor
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE get_fint(ndim,nrows,nnz,nstart,rp_k,cval_k,aval_k,dispin,fint)

  USE Precision
  USE implicit_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  ! Input variables
  INTEGER :: ndim, nrows, nnz, nstart
  INTEGER, DIMENSION(nrows+1) :: rp_k
  INTEGER, DIMENSION(nnz) :: cval_k
  REAL(kind=wp), DIMENSION(nnz) :: aval_k
  REAL(kind=wp), DIMENSION(nrows) :: dispin
  
  ! Output variables
  REAL(kind=wp), DIMENSION(nrows) :: fint

  ! Internal variables
  REAL(kind=wp), DIMENSION(ndim) :: disptemp
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumDispFrom
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd
  INTEGER :: i, j, k, m, counter1, counter2
  REAL(kind=wp) :: tempval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: req_rcv, req_snd
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: stat_rcv, stat_snd

  !
  ! Communicate how many pieces of disp to send to each other proc
  !

  IF (nprocs > 1) THEN
  	
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(req_rcv(1:nprocs))
     ALLOCATE(req_snd(1:nprocs))
     ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:nprocs))
     ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:nprocs))
     ALLOCATE(NumDispFrom(1:nprocs))
     DO i = 1, nprocs
        req_rcv(i) = 0
        req_snd(i) = 0
        DO j = 1, MPI_STATUS_SIZE
           stat_rcv(j,i) = 0
           stat_snd(j,i) = 0
        ENDDO
     ENDDO
     DO i = 1, nprocs
     	IF (i-1 /= myid) THEN
           CALL MPI_IRECV(NumDispFrom(i),1, &
              MPI_INTEGER,i-1,10,ROCSTAR_COMMUNICATOR, &
              req_rcv(i),ierr)
        ENDIF
     ENDDO
     DO i = 1, nprocs
        IF (i-1 /= myid) THEN
           !CALL MPI_ISEND(nrows,1,MPI_INTEGER, &
           !   i-1,10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
           CALL MPI_SEND(nrows,1,MPI_INTEGER, &
              i-1,10,ROCSTAR_COMMUNICATOR,ierr)
        ENDIF
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)
     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)
        
     
     
     !
     ! Communicate displacement to all the other procs
     !
          
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(req_rcv(1:nprocs))
     ALLOCATE(req_snd(1:nprocs))
     ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:nprocs))
     ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:nprocs))
     DO i = 1, nprocs
        req_rcv(i) = 0
        req_snd(i) = 0
        DO j = 1, MPI_STATUS_SIZE
           stat_rcv(j,i) = 0
           stat_snd(j,i) = 0
        ENDDO
     ENDDO
     ALLOCATE(frmproc(1:nprocs))
     DO i = 1, nprocs
     	IF (i-1 /= myid) THEN
           ALLOCATE(frmproc(i)%rcvbuf(1:4*INT(NumDispFrom(i)/3)))
           CALL MPI_IRECV(frmproc(i)%rcvbuf(1),4*INT(NumDispFrom(i)/3), &
                MPI_DOUBLE_PRECISION,i-1,10,ROCSTAR_COMMUNICATOR, &
                req_rcv(i),ierr)
        ENDIF
     ENDDO
     DO i = 1, nprocs
     	IF (i-1 /= myid) THEN
           ALLOCATE(bufsnd(1:4*LNumNp))
           counter1 = 0
           DO j = 1, GNumNp
              IF ((Global2Local(j) /= -1) .AND. (NodeProc(Global2Local(j)) == myid)) THEN
                 counter1 = counter1 + 1
                 bufsnd(4*counter1-3) = j
                 bufsnd(4*counter1-2) = dispin(3*counter1-2)
                 bufsnd(4*counter1-1) = dispin(3*counter1-1)
                 bufsnd(4*counter1)   = dispin(3*counter1)
              ENDIF
           ENDDO
           !CALL MPI_ISEND(bufsnd,4*INT(nrows/3),MPI_DOUBLE_PRECISION, &
           !     i-1,10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
           CALL MPI_SEND(bufsnd,4*INT(nrows/3),MPI_DOUBLE_PRECISION, &
                i-1,10,ROCSTAR_COMMUNICATOR,ierr)
           DEALLOCATE(bufsnd)
        ENDIF
     ENDDO
     !CALL MPI_WAITALL(NumCommProcsFrom2,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(NumCommProcs2,req_snd,stat_snd,ierr)
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)
     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)

  ENDIF
  
  ! Create a global displacement vector
  DO i = 1, GNumNp
     IF ((Global2Local(i) /= -1) .AND. (NodeProc(Global2Local(i)) == myid)) THEN
        disptemp( 3*i - 2 ) = dispin( 3*(i - (nstart-1)/3)-2 )
        disptemp( 3*i - 1 ) = dispin( 3*(i - (nstart-1)/3)-1 )
        disptemp( 3*i )     = dispin( 3*(i - (nstart-1)/3)   )
     ELSE
        DO j = 1, nprocs
           IF (j-1 /= myid) THEN
              DO m = 1, INT(NumDispFrom(j)/3)
                 IF (INT(frmproc(j)%rcvbuf(4*m-3)) == i) THEN
                    disptemp( 3*i - 2 ) = frmproc(j)%rcvbuf(4*m-2)
                    disptemp( 3*i - 1 ) = frmproc(j)%rcvbuf(4*m-1)
                    disptemp( 3*i )   = frmproc(j)%rcvbuf(4*m)
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  IF (nprocs > 1) DEALLOCATE(frmproc)


  !
  ! Multiply the K matrix by the displacement to find fint
  !
  CALL comp_row_vecmult(3*GNumNp,nrows,nnz,nstart,rp_k,cval_k,aval_k,disptemp,fint)
  
  
  
END SUBROUTINE get_fint

