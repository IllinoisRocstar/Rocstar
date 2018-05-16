
!!****
!!
!!  NAME
!!     IntLoad
!!
!!  FUNCTION
!!     This subroutine constructs the internal thermal load vector for this
!!     processor by multiplying its part of the stiffness or capacitance matrix by
!!     the global temperature vector.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global stiffness or capacitance matrix
!!     nrows -- The number of rows that are assigned to this processor
!!     nnz -- The number of nonzeros in this processor's part of the stiffness or capacitance matrix
!!     nstart -- The global index of the first row assigned to this processor
!!     rp -- The row mapping vector
!!     cval -- The collumn mapping vector
!!     aval -- The nonzero values in this processor's part of the stiffness or capacitance matrix
!!     tempin -- The part of the global temperature vector that's assigned to this processor
!!
!!  OUTPUTS
!!     rint -- The part of the global internal thermal load vector that's assigned to this processor
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE IntLoad(ndim,nrows,nnz,nstart,rp,cval,aval,tempin,rint)

  USE Precision
  USE implicit_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  ! ... Input variables
  INTEGER :: ndim, nrows, nnz, nstart
  INTEGER, DIMENSION(nrows+1) :: rp
  INTEGER, DIMENSION(nnz) :: cval
  REAL(kind=wp), DIMENSION(nnz) :: aval
  REAL(kind=wp), DIMENSION(nrows) :: tempin
  
  ! ... Output variables
  REAL(kind=wp), DIMENSION(nrows) :: rint

  ! ... Internal variables
  REAL(kind=wp), DIMENSION(ndim) :: temptemp
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumTempFrom
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd
  INTEGER :: i, j, k, m, counter1, counter2
  REAL(kind=wp) :: tempval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: req_rcv, req_snd
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: stat_rcv, stat_snd


  ! ... Communicate how many pieces of temp to send to each other proc
  IF ( nprocs > 1) THEN

     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(req_rcv(1:nprocs))
     ALLOCATE(req_snd(1:nprocs))
     ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:nprocs))
     ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:nprocs))
     ALLOCATE(NumTempFrom(1:nprocs))
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
           CALL MPI_IRECV(NumTempFrom(i),1, &
                MPI_INTEGER,i-1,10,ROCSTAR_COMMUNICATOR, &
                req_rcv(i),ierr)
        ENDIF
     ENDDO
     DO i = 1, nprocs
        IF (i-1 /= myid) THEN

           CALL MPI_SEND(nrows,1,MPI_INTEGER, &
                i-1,10,ROCSTAR_COMMUNICATOR,ierr)
        ENDIF
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)

     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)




     ! ... Communicate temperature to all the other procs
     ! ... receive
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
        ! ... if processor i is not this processor (i=1=processor 0)
     	IF (i-1 /= myid) THEN
         ! ... Allocate a buffer of size 2 times the number of nodes from
         ! ... processor j.  The buffer is orgnized as {node #, temp, node #, temp...}
           ALLOCATE(frmproc(i)%rcvbuf(1:2*INT(NumTempFrom(i))))
           CALL MPI_IRECV(frmproc(i)%rcvbuf(1),2*INT(NumTempFrom(i)), &
                MPI_DOUBLE_PRECISION,i-1,10,ROCSTAR_COMMUNICATOR, &
                req_rcv(i),ierr)
        ENDIF
     ENDDO

     ! ... send
     DO i = 1, nprocs
        IF (i-1 /= myid) THEN
           ! ... if processor i is not this processor
           ! ... allocate a buffer 2 times the number of the nodes 
           ! ... on this processor to send to processor i
           ALLOCATE(bufsnd(1:2*LNumNp))
           counter1 = 0
           DO j = 1, GNumNp
              ! ... if global node j belongs to this processor
              ! ... then the corresponding temperature is copied
              ! ... to a buffer to go to processor i
              IF ((Global2Local(j) /= -1) .AND. (NodeProc(Global2Local(j)) == myid)) THEN
                 ! ... if global node j belongs to this processor
                 ! ... then counter1 determines the index in tempin where
                 ! ... the corresponding temperature is stored
                 counter1 = counter1 + 1
                 bufsnd(2*counter1-1) = j
                 bufsnd(2*counter1)   = tempin(counter1)
              ENDIF
           ENDDO
           ! ... send the buffer of length 2*nodes on this processor
           ! ... to processor (i-1).
           CALL MPI_SEND(bufsnd,2*nrows,MPI_DOUBLE_PRECISION, &
                i-1,10,ROCSTAR_COMMUNICATOR,ierr)
           DEALLOCATE(bufsnd)
        ENDIF
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)

     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)
     
  ENDIF
  
  ! ... create a global temperature vector
  DO i = 1, GNumNp
     ! ... if global node i belongs to this processor
     ! ... then the corresponding temperature is copied
     ! ... in from the input array tempin
     IF ((Global2Local(i) /= -1) .AND. (NodeProc(Global2Local(i)) == myid)) THEN
        temptemp( i ) = tempin( i - (nstart-1))
     ELSE
        ! ... if global node i does not belong to this processor
        ! ... then find which processor it does belong to
        ! ... and get the temperature from the buffer received from
        ! ... that processor
        DO j = 1, nprocs
           IF (j-1 /= myid) THEN
              ! ... loop through the number of temperature values from processor j
              DO m = 1, NumTempFrom(j)
                 ! ... frmproc(j)%rcvbuf(2*m-1) is where the node number
                 ! ... is stored for the temperature value stored at 
                 ! ... frmproc(j)%rcvbuf(2*m)
                 IF (INT(frmproc(j)%rcvbuf(2*m-1)) == i) THEN
                    temptemp( i )   = frmproc(j)%rcvbuf(2*m)
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  IF (nprocs > 1) DEALLOCATE(frmproc)



  ! ... Multiply the input matrix by the temperature to find internal load
  CALL comp_row_vecmult(GNumNp,nrows,nnz,nstart,rp,cval,aval,temptemp,rint)
  
END SUBROUTINE IntLoad
