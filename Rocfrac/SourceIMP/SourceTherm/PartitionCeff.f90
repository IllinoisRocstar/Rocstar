
!!****
!!
!!  NAME
!!     Partition_Ceff
!!
!!  FUNCTION
!!     This subroutine partitions the full effective capacitance matrix
!!     into the blocks corresponding to the prescribed-prescribed, 
!!     prescribed-free and free-free DOFs.  This is accomplished by converting
!!     Ceff into its explicit form, partitioning, and converting the free-free
!!     part back to CRS format.  The prescribed-prescribed and free-prescribed
!!     blocks can stay in explicit form because they are not going to be input
!!     to BlockSolve95.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global Ceff matrix
!!     nrows -- The number of rows assigned to this processor
!!     nnz -- The number of nonzeros in section of the Ceff matrix on this processor
!!     nstart -- The global index of the first row assigned to this processor
!!     rp1 -- The row mapping vector
!!     cval -- The collumn mapping vector
!!     aval -- The nonzero value vector
!!     ProcTemp -- The part of the global temperature vector that's assigned to this processor
!!
!!  OUTPUTS
!!     newnrows -- The number of rows assigned to this proc after BCs have been removed
!!     newnstart -- The global index of the first row assigned to this proc after BCs have been removed
!!     newndim -- The size of the global CeffFF matrix after BCs have been removed
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE Partition_Ceff(ndim,nrows,nnz,nstart,rp1,cval,aval,newnrows,newnstart,newndim,ProcTemp,global)

  USE removeBCs_global
  USE comp_row_global
  USE implicit_global
  USE Precision
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL) :: global

  ! ... Arguments
  INTEGER :: ndim, nnz, nstart, nrows
  REAL(kind=wp), DIMENSION(nnz) :: aval
  INTEGER, DIMENSION(nnz) :: cval
  INTEGER, DIMENSION(nrows+1) :: rp1
  REAL(kind=wp), DIMENSION(LNumNp) :: ProcTemp

  INTEGER :: newnrows
  INTEGER :: newnstart
  INTEGER :: newndim

  ! ... local variables
  INTEGER :: numdisp  ! Local number of displacement BCs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dispbc  ! Local displacement BCs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: tempintv
  
  INTEGER :: i, j, n, m, counter1, counter2, counter3, counter4, counter5

 
  INTEGER :: ncols, iDOF, jDOF

  REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) :: CeffExplicit, CeffFF, CeffFP, CeffPF, CeffPP


  ! ... Deallocate removeBCs_global variables
  IF(ALLOCATED(GTempBC)) DEALLOCATE(GTempBC)
  IF(ALLOCATED(NumTempProc)) DEALLOCATE(NumTempProc)

  ! ... Count the number of temperature boundary conditions in this 
  ! ... processes partition.
  numdisp = 0
  DO i = 1, global%NumNp
     IF (NodeProc(i) == myid) THEN
        ! ... If the value of the flag at node i on in this processes
        ! ... is 8, add 1 to the count of temperature BCs.
        IF(node_flag(i,1) == 8) THEN
           numdisp = numdisp + 1
        ENDIF

     ENDIF
  ENDDO
  print*,myid,' number of local disp BCs = ',numdisp

  ! Sum the local disp BCs, then broadcast to all procs
  CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)

  ! ... MPI_REDUCE sums values on all processes to one value. (SUM(numdisp)=GNumTemp)
  CALL MPI_REDUCE(numdisp, GNumTemp, 1, MPI_INTEGER, MPI_SUM, 0 , ROCSTAR_COMMUNICATOR, ierr)
  ! ... The number GNumTemp is then communicated to all processes, all know the global
  ! ... number of temperature boundary condtions.
  CALL MPI_BCAST(GNumTemp, 1, MPI_INTEGER, 0, ROCSTAR_COMMUNICATOR, ierr)

  ! ... newndim is the size of the Global capacitance matrix after the temperature BC's are removed.
  newndim = GNumNp - GNumTemp
  print*,myid,' number of global disp BCs = ',GNumTemp

  ! ... Figure out which nodes have perscribed temperatures
  ! ... If this processes partition contains prescribed temperatures, then
  ! ... allocate an array to store the global index of that temperature BC.
  IF(numdisp > 0) THEN
     ALLOCATE(dispbc(1:numdisp))
     dispbc(:) = 0
     counter2 = 0

     DO m = 1, GNumNp
        DO i = 1, global%NumNp
           ! ... If global node m is owned by this process and there is a temperature
           ! ... BC at that node, record the node number as the next entry in array dispbc.
           IF (Local2Global(i) == m) THEN
              IF (NodeProc(i) == myid) THEN
                 IF(node_flag(i,1) == 8) THEN  ! Imposed constant nodal displacement
                    counter2 = counter2 + 1
                    dispbc(counter2) = m
                  ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  If(numdisp >0) print*,myid,' local disp bcs = ',dispbc(:)

  ! ... Have each proc in turn broadcast its number of temperature BCs to the other procs
  ALLOCATE(NumTempProc(1:nprocs))
  NumTempProc(:) = 0
  NumTempProc(myid+1) = numdisp
  DO i = 1, nprocs
     CALL MPI_BCAST(NumTempProc(i), 1, MPI_INTEGER, i-1, ROCSTAR_COMMUNICATOR, ierr)
  ENDDO

  ! ... The global row index of the first row on this process is the old index
  ! ... minus the number of temperature BCs on global nodes in the range from 
  ! ... node 1 to node nstart(the old starting row index).
  newnstart = nstart
  DO i = 1, myid
     newnstart = newnstart - NumTempProc(i)
  ENDDO
  print*,myid,' number of disp bcs on other procs = ',NumTempProc(:)

  ! ... Have each proc in turn broadcast the node numbers that it owns which have temperature 
  ! ... boundary conditions.  The global list of these node numbers will be stored in the array
  ! ... GTempBC.
  ALLOCATE(GTempBC(1:GNumTemp))
  GTempBC(:) = 0
  counter1 = 0
  ! ... Loop through the processes
  DO i = 1, nprocs
     ! ... If there are temperature boundary conditions on process i then allocate the temporary array
     ! ... tempintv to the number of temperature BCs on process i.
     IF (NumTempProc(i) > 0) THEN
        ALLOCATE(tempintv(1:NumTempProc(i)))
        tempintv(:) = 0
        ! ... It process i is this process, fill the array tempintv with the node numbers of the 
        ! ... temperature boundary conditions on this process.
        if (i-1 == myid) tempintv(:) = dispbc(:)
        ! ... The process i broadcasts the array tempintv to all other processes.
        CALL MPI_BCAST(tempintv(1), NumTempProc(i), MPI_INTEGER, i-1, ROCSTAR_COMMUNICATOR, ierr)
        ! ... Loop through the number of temperature BC nodes on process i and record them on 
        ! ... array GTempBC.
        DO j = 1, NumTempProc(i)
           counter1 = counter1 + 1
           GTempBC(counter1) = tempintv(j)
        ENDDO
        DEALLOCATE(tempintv)
     ENDIF
  ENDDO
  print*,myid,' global disp bcs = ',GTempBC(:)

  ! ... convert Ceff to explicit form, CeffExplicit

  ! ... Allocate explicit matrix
  ALLOCATE(CeffExplicit(nrows,ndim))

  ! ... initialize CeffExplicit
  CeffExplicit = 0.0d0

  ! ... n is the CRS index for the aval and cval vectors
  n = 1

  ! ... First loop through the rows i
  DO i = 1, nrows
     ! ... loop through the columns j
     DO j = 1, ndim
        ! ... is n in the index range for row i?
        IF ((rp1(i) < n).AND.(rp1(i+1) >= n)) THEN
           ! ... if this is the correct column corresponding to index n
           IF (cval(n) == j - 1) THEN
              CeffExplicit(i,j) = aval(n)
              n = n + 1
           ENDIF
        ENDIF
        ! ... if end of CRS vectors, no need to keep looping
        IF (n > nnz) EXIT
     ENDDO
     ! ... if end of CRS vectors, no need to keep looping
     IF (n > nnz) EXIT
  ENDDO



  ! ... Partition matrix CeffExplicit into explicit forms of
  ! ... - CeffFF: Free-free DOFs
  ! ... - CeffPP: Prescribed-prescribed DOFs
  ! ... - CeffPF: Prescribed-free DOFs
  ! ... - CeffFP: Free-prescribed DOFs

  ! ... Allocate the partitioned matrices
  ALLOCATE(CeffFF(nrows - numdisp, ndim - GNumTemp))
  ALLOCATE(CeffPP(numdisp,GNumTemp))
  ALLOCATE(CeffPF(numdisp, ndim - GNumTemp))
  ALLOCATE(CeffFP(nrows - numdisp, GNumTemp))
  ALLOCATE(Ceff_fpTp(nrows - numdisp))

  ! ... Begin partitioning
  counter1 = 0
  counter2 = 0

  ! ... outer loop through rows i
  DO i = nstart, nstart + nrows - 1
     ! ... iDOF stores the DOF number of the current essential BC
     iDOF = 0

     ! ... reset the nmber of prescribed BSs counted in a row (counter 2)
     counter2 = 0

     ! ... check the DOFs of all local prescribed BCs to see if they are in
     ! ... this row
     DO n = 1, numdisp

        IF(i == dispBC(n)) THEN
           ! ... if i is a prescribed DOF, then increase row BC counter1
           counter1 = counter1 + 1


           ! ... if i is a prescribed DOF, then record the value of the DOF
           iDOF = i

           ! ... stop searching through dispBC, i can only equal one of them
           EXIT
        ENDIF
     ENDDO
     jDOF = 0
     ! ... if row i corresponds to a prescribed DOF, then fill CeffPP and CeffPF
     IF(iDOF > 0) THEN
        ! ... loop through the columns in prescribed row i
        DO j = 1, ndim
           ! ... check the DOFs of all local prescribed BCs to see if they are in
           ! ... this column, j
           DO n = 1, GNumTemp
              IF(j == GTempBC(n)) THEN
                 ! ... if i is a prescribed DOF, then increase column BC counter2
                 counter2 = counter2 + 1
                 ! ... if i is a prescribed DOF, then record the value of the dof
                 jDOF = j
                 ! ... stop searching through dispBC, j can only equal one of them
                 EXIT
              ENDIF
           ENDDO
           IF(jDOF > 0) THEN
              ! ... if row j corresponds to a prescribed DOF, then fill CeffPP

              CeffPP(counter1,counter2) = CeffExplicit(iDOF-nstart+1,jDOF)

              ! ... reset jDOF for next loop
              jDOF = 0
           ELSE

              ! ... if row j corresponds to a free DOF, then fill CeffPF
              CeffPF(counter1,j-counter2) = CeffExplicit(iDOF-nstart+1,j)
           ENDIF
        ENDDO
     ELSE
        ! ... loop through the columns in free row i
        DO j = 1, ndim

           ! ... check the DOFs of all local prescribed BCs to see if they are in
           ! ... this column, j
           DO n = 1, GNumTemp

              IF(j == GTempBC(n)) THEN
                 ! ... if i is a prescribed DOF, then increase column BC counter2
                 counter2 = counter2 + 1
                 ! ... if i is a prescribed DOF, then record the value of the dof
                 jDOF = j
                 ! ... stop searching through dispBC, j can only equal one of them
                 EXIT
              ENDIF
           ENDDO


           IF(jDOF > 0) THEN
              ! ... if row j corresponds to a prescribed DOF
              CeffFP(i-counter1-nstart+1,counter2) = CeffExplicit(i-nstart+1,jDOF)

              ! ... reset jDOF for next loop
              jDOF = 0
           ELSE
              ! ... if row j corresponds to a free DOF, then fill CeffFF
              CeffFF(i-counter1-nstart+1,j-counter2) = CeffExplicit(i-nstart+1,j)
           ENDIF
        ENDDO
     ENDIF
  ENDDO

  ! ... Now convert CeffFF back to CRS format

  ! ... Count nonzeros in new matrix
  nnz_temp = 0
  counter2 = 0
  ! ... Loop throught the number of rows on this process
  DO i = 1, nrows
     ! ... Loop from the beginning to the end of each row, keeping track for the column
     ! ... index using counter2.
     DO j = rp1(i)+1, rp1(i+1)
        counter2 = counter2 + 1
        counter1 = 0
        DO m = 1, GNumTemp
           ! ... If the global index of a temperature BC matches either the row or column
           ! ... indices, then do not increment the number of non-zeros.
           IF((GTempBC(m) == i + nstart - 1).OR.(GTempBC(m) == cval(counter2)+1)) THEN
              counter1 = 1
           ENDIF
        ENDDO

        IF(counter1 == 0) nnz_temp = nnz_temp + 1
     ENDDO
  ENDDO

  ! ... Count rows in new matrix
  newnrows = 0
  counter2 = 0
  ! ... Loop throught the number of rows on this process
  DO i = 1, nrows
     counter1 = 0
     DO m = 1, GNumTemp
        ! ... If the global index of a termperature BC matches the row index, then
        ! ... do not increment the new number of rows on this process.
        IF(GTempBC(m) == i + nstart - 1) THEN
           counter1 = 1
        ENDIF
     ENDDO
     if(counter1==1) print*,myid,' row removed at ',i+nstart-1
     IF(counter1 == 0) newnrows = newnrows + 1
  ENDDO

  ! ... Allocate variables for new matrix
  ALLOCATE(rp_temp(1:newnrows+1))
  ALLOCATE(cval_temp(1:nnz_temp))
  ALLOCATE(aval_temp(1:nnz_temp))

  ! ... n is the CRS index for the vectors aval_temp and cval_temp
  n = 0
  ! ... loop through the rows in the matrix CeffFF
  DO i = 1, newnrows
     rp_temp(i) = n
     ! ... loop through the columns in the matrix CeffFF
     DO j = 1, newndim
        ! ... check for a non-zero at i,j in CeffFF
        IF (CeffFF(i,j) > 0.0d0) THEN
           n = n + 1
           ! ... if there is a non-zero, it is entered into position n of 
           ! ... aval_temp and the column j(-1) is entered into positon n 
           ! ... of cval_temp
           aval_temp(n) = CeffFF(i,j)
           cval_temp(n) = j - 1
        ENDIF
     ENDDO
  ENDDO
  ! ... fill in the last entry of rp_temp

  rp_temp(newnrows+1) = nnz_temp



  ! ... multiply CeffFP with the global prescribed temperature vector
  CALL PrescribedLoad(nrows,nstart,numdisp,dispBC,CeffFP,ProcTemp)



END SUBROUTINE Partition_Ceff


!!****
!!
!!  NAME
!!     PrescribedLoad
!!
!!  FUNCTION
!!     This subroutine multiplies the Free-Prescribed partition of the effective thermal capacitance
!!     matrix with the list of global prescribed temperatures.
!!
!!  INPUTS
!!     ncols -- The number of columns in Ceff_fp, same as number of prescribed temperatures
!!     nrows -- The number of rows in Ceff_fp
!!     nstart -- The global index of the first row assigned to this processor
!!     GNumTemp -- The global number of prescribed temperatures
!!     GTempBC -- Vector containing global nodes of the prescribed temperatures 
!!     NumpTemp -- The local number of prescribed temperatures
!!     TempBC -- Vector containing local nodes of the prescribed temperatures
!!     Ceff_fp -- The free-perscribed partition of the effective thermal capacitance matrix
!!     ProcTemp -- The part of the global temperature vector that's assigned to this processor 
!!
!!  OUTPUTS
!!     Ceff_fpTp -- The product of Ceff_fp and the global prescribed temperature vector
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE PrescribedLoad(nrows,nstart,NumTemp,TempBC,Ceff_fp,ProcTemp)

  USE Precision
  USE implicit_global
  USE comp_row_global
  USE ROCSTAR_RocFracComm
  USE removeBCs_global

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  ! ... Input variables
  ! ... the Ceff_fp matrix
  INTEGER :: nrows, NumTemp
  REAL(kind=wp), DIMENSION(nrows-NumTemp,GNumTemp) :: Ceff_fp
  ! ... the local temperature boundary condions
  INTEGER, DIMENSION(NumTemp) :: TempBC
  ! ... first row of this processor
  INTEGER :: nstart
  ! ... all temperatures on this processor
  REAL(kind=wp), DIMENSION(LNumNp) :: ProcTemp

  ! ... Internal variables
  ! ... the prescribed temperature vector
  REAL(kind=wp), DIMENSION(GNumTemp) :: GTempBCValue

  ! ... communication
  REAL(kind=wp), DIMENSION(GNumNp) :: temptemp
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumTempFrom
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd
  INTEGER :: i, j, k, m, counter1, counter2
  REAL(kind=wp) :: tempval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: req_rcv, req_snd
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: stat_rcv, stat_snd

  !
  ! Communicate how many pieces of temp to send to each other proc
  !

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
     ! Communicate temperature to all the other procs
     !
     
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
                 ! ... then counter1 determines the index in ProcTemp where
                 ! ... the corresponding temperature is stored
                 counter1 = counter1 + 1
                 bufsnd(2*counter1-1) = j
                 bufsnd(2*counter1)   = ProcTemp(counter1)
              ENDIF
           ENDDO
           !CALL MPI_ISEND(bufsnd,4*INT(nrows/3),MPI_DOUBLE_PRECISION, &
           !     i-1,10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)

           ! ... send the buffer of length 2*nodes on this processor
           ! ... to processor (i-1).
           CALL MPI_SEND(bufsnd,2*nrows,MPI_DOUBLE_PRECISION, &
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

  ! ... create a global temperature vector
  DO i = 1, GNumNp
     ! ... if global node i belongs to this processor
     ! ... then the corresponding temperature is copied
     ! ... in from the input array ProcTemp
     IF ((Global2Local(i) /= -1) .AND. (NodeProc(Global2Local(i)) == myid)) THEN
        temptemp( i ) = ProcTemp( i - (nstart-1))
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


  ! ... Get only the prescribed temperatures
  ! ... put in GTempBCValues
  DO i = 1, GNumNp
     ! ... loop through prescribed temperature nodes
     DO j = 1, GNumTemp
        IF (i == GTempBC(j)) THEN
           GTempBCValue(j)= temptemp(i)
        ENDIF
     ENDDO
  ENDDO

Ceff_fpTp = 0.0d0

  ! ... multiply Ceff_fp with GTempBCValue
  DO i = 1, nrows-Numtemp
     DO j = 1, GNumTemp
        Ceff_fpTp(i) = Ceff_fpTP(i) + Ceff_fp(i,j) * GTempBCValue(j)
     ENDDO
  ENDDO

  
END SUBROUTINE PrescribedLoad




