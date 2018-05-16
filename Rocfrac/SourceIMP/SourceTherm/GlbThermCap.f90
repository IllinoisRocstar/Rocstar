!!****
!!
!!  NAME
!!     GlbThermCap
!!
!!  FUNCTION
!!     This subroutine constructs the global consistent capacitance matrix,
!!     then puts it into CRS format.  The matrix information is stored in
!!     global variables: rp_c, aval_c, cval_c, nnz_c, nrows_ktc, nstart_ktc
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE GlbThermCap(global)

  USE implicit_global
  USE comp_row_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL) :: global

  INTEGER :: i, j, m, p, inode, jnode, counter
  
  REAL(kind=wp) :: tempCval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: AmountToReceive, AmountToSend
  
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd
  
  ! ... Gauss point coordinates
  REAL(kind=wp), DIMENSION(1:3,1:8) :: ri = RESHAPE( &
       (/-0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626, 0.577350269189626, 0.577350269189626, &
         -0.577350269189626, 0.577350269189626, 0.577350269189626/),(/3,8/) )



  ! ... Construct the global capacitance matrix
  IF(myid == 0) PRINT*,'CONSTRUCTING THE CAPACITANCE MATRIX'
  CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)

  IF ( myid == 0) PRINT*,'RocTherm, GlbThermCap :: before LocThermCap'

  ! ... Construct the capacitance matrix as a consistent matrix
  CALL LocThermCap_v3d8(global%NumElVol,global%NumNp,global%NumMatVol,global%MeshCoor,global%ElConnVol,global%MatIdVol,ri,global%rho,global%Cp,global%ElConnVol)
  
  IF(myid == 0) PRINT*,'FINISHED CONSTRUCTING THE CAPACITANCE MATRIX'
  CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
 

  ! ... number of non-zeros and CRS vectors ' '_temp computed by LocThermCap_v3d8
  nnz_c = nnz_temp
  ALLOCATE(cval_c(1:nnz_c))
  ALLOCATE(aval_c(1:nnz_c))
  ALLOCATE(rp_c(1:GNumNp+1))
  rp_c = rp_temp
  cval_c = cval_temp
  aval_c = aval_temp
  DEALLOCATE(rp_temp)
  DEALLOCATE(cval_temp)
  DEALLOCATE(aval_temp)


  IF(myid == 0) PRINT*,'COMMUNICATING CAPACITANCE MATRIX TO ALL PROCESSES'

  ! ... Communicate mass to the other procs
  ! ... boundary nodes on this process are communicated to the processes
  ! ... to which they are assigned
  IF (nprocs > 1) THEN
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ! ... initialize request and status variables for MPI sends and receives
     ALLOCATE(req_rcv(1:nprocs))
     ALLOCATE(req_snd(1:nprocs))
     ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:nprocs))
     ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:nprocs))
     DO i = 1, nprocs
        req_rcv(i) = 0
        req_snd(i) = 0
        DO j = 1, MPI_STATUS_SIZE
           stat_rcv(1,i) = 0
           stat_snd(j,i) = 0
        ENDDO
     ENDDO
     
     ! ... number of nodes each process will send to this process
     ! ... (nodes owned by this process but in the partition of other processes)
     ALLOCATE(AmountToReceive(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        CALL MPI_IRECV(AmountToReceive(i),1, &
             MPI_INTEGER,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
             req_rcv(i),ierr)
     ENDDO
print*,'RocTherm, GlbCap :: after MPI_IRECV',myid
     ! ... number of nodes each process will be sent from this process
     ! ... (nodes owned by other processes but in the partition of this one)
     ALLOCATE(AmountToSend(1:NumCommProcs1))
     DO i = 1, NumCommProcs1
        counter = 0
        ! ... loop through the nodes going to be sent to process i
        DO j = 1, NumCommNodes1(i)
           ! ... inode is the global node number corresponding to the jth node going to process i
           ! ... from this process
           inode = Local2Global(CommNodes1(i,j))
           DO jnode = 1, GNumNp
              tempCval = 0.0
              ! ... return the value (tempCval) stored in CRS at explicit matrix
              ! ... location (inode, jnode)
              CALL comp_row_getval(GNumNp,GNumNP,nnz_c,1,rp_c,cval_c,aval_c, &
                   inode,jnode,tempCval)



              ! ... if there is a value at this node, then add 1 to the ammount to send
              ! ... to process i.
              IF (tempCval /= 0.0) THEN
                 counter = counter + 1
              ENDIF
           ENDDO
        ENDDO

        ! ... amount to send process i is the sum of all non-zeros in the loop above 
        AmountToSend(i) = counter

        CALL MPI_SEND(AmountToSend(i),1,MPI_INTEGER, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
     ENDDO

     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)


     ! ... Communicate parts of the capcitance matrix to other procs

     ! ... receive information from other processes
     ! ... (nodes owned by this process but in the partition of other processes)
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(frmproc(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(frmproc(i)%rcvbuf(1:3*AmountToReceive(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),3*AmountToReceive(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
             req_rcv(i),ierr)
     ENDDO
     
     ! ... send information to other processes
     ! ... (nodes owned by other processes but in the partition of this one)    
     DO i = 1, NumCommProcs1
        ALLOCATE(bufsnd(3*AmountToSend(i)))
        bufsnd(:) = 0.0
        counter = 0
        ! ... loop through the nodes going to be sent to process i
        DO j = 1, NumCommNodes1(i)
           ! ... inode is the global node number corresponding to the jth node going to process i
           ! ... from this process
           inode = Local2Global(CommNodes1(i,j))
           DO jnode = 1, GNumNp
              tempCval = 0.0
              ! ... return the value (tempCval) stored in CRS at explicit matrix
              ! ... location (inode, jnode)
              CALL comp_row_getval(GNumNp,GNumNP,nnz_c,1,rp_c,cval_c,aval_c, &
                   inode,jnode,tempCval)

              ! ... if there is a value at this node, put it in the location '3*counter'
              ! ... (preceeded by the row and column locations (fortran)) in the buffer to
              ! ... be sent to process i.
              IF (tempCval /= 0.0) THEN
                 counter = counter + 1
                 bufsnd(3*counter-2) = inode
                 bufsnd(3*counter-1) = jnode
                 bufsnd(3*counter) = tempCval
              ENDIF
           ENDDO
        ENDDO

        CALL MPI_SEND(bufsnd,3*AmountToSend(i),MPI_DOUBLE_PRECISION, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)

        DEALLOCATE(bufsnd)
     ENDDO

     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)

     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)
  ENDIF

  IF(myid == 0) PRINT*,'FINISHED COMMUNICATING CAPACITANCE MATRIX TO ALL PROCESSES'

  ! ... add values from other processes into this processes capacitance matrix
  IF (nprocs > 1) THEN
     DO i = 1, NumCommProcsFrom1
        DO j = 1, AmountToReceive(i)
           CALL comp_row_addval(GNumNp,GNumNp,nnz_c,1,rp_c,cval_c,aval_c, &
                INT(frmproc(i)%rcvbuf(3*j-2)),INT(frmproc(i)%rcvbuf(3*j-1)),frmproc(i)%rcvbuf(3*j))

           ! ... if after the addition, the number of non-zeroes has changed
           ! ... reallocate the aval and cval vectors to accomodate
           IF (nnz_temp /= nnz_c) THEN
              DEALLOCATE(cval_c,aval_c)
              ALLOCATE(cval_c(nnz_temp),aval_c(nnz_temp))
              nnz_c = nnz_temp
              rp_c = rp_temp
              cval_c = cval_temp
           ENDIF
           aval_c = aval_temp
           DEALLOCATE(rp_temp,cval_temp,aval_temp)
        ENDDO
     ENDDO
     DEALLOCATE(frmproc)
     DEALLOCATE(AmountToReceive)
     DEALLOCATE(AmountToSend)
  ENDIF

  ! ... Resize the matrix to the size needed on this proc
  CALL comp_row_resize(GNumNp,GNumNp,nnz_c,1,rp_c,cval_c,aval_c,nstart_ktc,nrows_ktc)
  DEALLOCATE(rp_c,cval_c,aval_c)
  ALLOCATE(rp_c(1:nrows_ktc+1))
  ALLOCATE(cval_c(1:nnz_temp))
  ALLOCATE(aval_c(1:nnz_temp))
  nnz_c = nnz_temp
  rp_c = rp_temp
  cval_c = cval_temp
  aval_c = aval_temp
  DEALLOCATE(rp_temp,cval_temp,aval_temp)
  

END SUBROUTINE GlbThermCap
