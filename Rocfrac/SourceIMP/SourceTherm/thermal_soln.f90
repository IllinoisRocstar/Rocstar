!!****
!!
!!  NAME
!!     thermal_soln
!!
!!  FUNCTION
!!     This subroutine computes the thermal solution
!!     at t=t+dt using the beta method.
!!
!!  INPUTS
!!     delt -- Length of one timestep
!!     t -- Current timestep
!!     rext_in -- External loads
!!     global -- Holds all global quantities
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     MPI
!!
!!****


SUBROUTINE thermal_soln(delt,t,rext_in,global,istep)

  USE implicit_global
  USE comp_row_global
  USE Precision
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  ! ... variables
  TYPE(ROCFRAC_GLOBAL) :: global
  REAL(kind=wp) :: delt, t
  REAL(kind=wp), DIMENSION(1:global%NumNp) :: rext_in
  INTEGER :: istep
  
  ! ... local variables
  INTEGER :: BS95debug
  INTEGER :: i, j, k, m, n, p, counter, ii
  INTEGER :: newnrows_ceff, newnstart_ceff, newndim
  REAL(kind=wp) :: contol, beta
  REAL(kind=wp), DIMENSION(1:LNumNp) :: rint, rint2, rext_new, Temp, reff
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: newreff, newTemp
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd, ftemp
  REAL(wp) :: per1, delay

  
  ! ... Initialize
  rint(:) = 0.0
  rint2(:) = 0.0
  Temp(:) = 0.0
  contol = 1.0e-09
  beta = .5
  
  IF ( global%debug_state ) THEN
     BS95debug = 1
  ELSE
     BS95debug = 0
  END IF

  ! ... Assemble the temperature and load vectors
  ! ... loop through all node points
  DO i = 1, global%NumNp
     ! ... if a node belongs to the current process
     IF ( NodeProc(i) == myid ) THEN
        ! ... k is the row index for the given node
        k = Local2Global(i) - nstart_ktc + 1

        ! ... rext_new are external thermal loads for the timestep currently being calculated
        rext_new( k ) = rext_in( i )
        
        ! ... temperature at each node at last timestep
        Temp( k ) = global%Temperature( i )

     END IF
  END DO

  ! ... Output list of nodal temperatures
!  if (real(istep/10) == real(istep)/10) CALL RockOut_8node(global,istep)
!  if (istep < 31) CALL RockOut_8node(global,istep)
  If (ALLOCATED(r_old) .eqv. .false.) then
     ALLOCATE(rext_imp(LNumNp),r_old(LNumNp))
     rext_imp = 0.0
     r_old = 0.0
  endif



  ! ... if it is the first timestep, then compute the effective capacitance matrix and 
  ! ... apply the boundary conditions.
  ! ... this will only be done once unless
  ! ... (1) the imposed temperature boundary conditions change 
  ! ... (2) the timestep changes
  ! ... (3) the material properties change
  IF (istep == 1) THEN
     
     IF(ALLOCATED(aval_ceff).EQV. .TRUE.) DEALLOCATE(aval_ceff)
     IF(ALLOCATED(cval_ceff).EQV. .TRUE.) DEALLOCATE(cval_ceff)
     IF(ALLOCATED(rp_ceff).EQV. .TRUE.) DEALLOCATE(rp_ceff)
     IF(ALLOCATED(Ceff_fpTp).EQV. .TRUE.) DEALLOCATE(Ceff_fpTp)



     ! ... Create the effective capacitance matrix
     !
     !
     ! m_eff = (1/dt * C + beta * K)
     !
     !

     PRINT*,'FORMING EFFECTIVE CAPACITANCE MATRIX',myid

     ! ... For some reason, the first GNumNp is not used at all in 'comp_row_add'
     ! ... The ouput from 'comp_row_add are stored in the global variables ..._temp
     CALL comp_row_add(GNumNp,GNumNp,nrows_ktc,nrows_ktc,nnz_kt,nnz_c,1, &
          rp_kt,cval_kt,beta*aval_kt,1,rp_c,cval_c,(1/delt)*aval_c)

     ALLOCATE(rp_ceff(nrows_ktc+1),cval_ceff(nnz_temp),aval_ceff(nnz_temp))

     ! ... compressed row vectors for the effective capacitance matrix
     ! ... the vectors are stored under the global variables ..._ceff (implicit_global.f90)
     nnz_ceff = nnz_temp
     nstart_ceff = nstart_ktc
     nrows_ceff = nrows_ktc
     rp_ceff = rp_temp
     cval_ceff = cval_temp
     aval_ceff = aval_temp
     DEALLOCATE(rp_temp,cval_temp,aval_temp)

     ! ... Enforce imposed temperature boundary conditions.
     CALL EnforceThermalBC(global%NumNp,LNumNp,Temp,node_flag,boundary_value,0,myid)

     ! ... partition effective capacitance matrix into CeffFF, CeffPP, CeffFP.
     ! ... multiply CeffFP by the prescribed temperatures
     CALL  Partition_Ceff(GNumNp,nrows_ceff,nnz_ceff,nstart_ceff,rp_ceff,cval_ceff,aval_ceff, &
          newnrows_ceff,newnstart_ceff,newndim,Temp,global)

     ! ... deallocate effective capacitance matrix, dimensions have changed after BCs removed

     DEALLOCATE(rp_ceff,cval_ceff,aval_ceff)

     ! ... position and size of local effective capacitance matrix have changed due to BC's applied 
     ! ... to global matrix
     nrows_ceff = newnrows_ceff
     nstart_ceff = newnstart_ceff

     ! ... reallocate local effective capacitance matrix to new dimensions, and assign new values 
     ALLOCATE(rp_ceff(nrows_ceff+1),cval_ceff(nnz_temp),aval_ceff(nnz_temp))
     nnz_ceff = nnz_temp
     rp_ceff = rp_temp
     cval_ceff = cval_temp
     aval_ceff = aval_temp
     DEALLOCATE(rp_temp,cval_temp,aval_temp)

     CALL BS95SETUP(GNumNp,nnz_ceff,nstart_ceff-1,nrows_ceff,rp_ceff,cval_ceff,aval_ceff,1,BS95debug)
  ENDIF

  ! ... initialize thermal load vectors
  rint(:) = 0.0
  rint2(:) = 0.0
  reff(:) = 0.0


  ! ... Calculate internal thermal loads at time = t
  !
  !                   .
  !     (1/dt * C - (1-beta) * K) * Temp 
  !                                  k 
  !
  !     where C: thermal capacitance, K: thermal stiffness, dt: timestep,
  !           Temp : temperature at timestep k
  !               k
  !

  ! ... calculate rint = [-(1-beta) * K] * Temp 
  CALL IntLoad(GNumNp,nrows_ktc,nnz_kt,nstart_ktc,rp_kt,cval_kt,(-1.0*(1.0-beta)*aval_kt),Temp,rint)

  ! ... calculate rint2 = [1/dt * C] * Temp
  CALL IntLoad(GNumNp,nrows_ktc,nnz_kt,nstart_ktc,rp_c,cval_c,1/delt*aval_c,Temp,rint2)


  ! ... Apply boundary conditions
  counter=0
  DO m = 1, GNumNp
     DO i = 1, global%NumNp
        IF (Local2Global(i) == m) THEN
           IF (NodeProc(i) == myid) THEN
              counter = counter + 1
              IF(node_flag(i,1) == 7) THEN
                 ! ... imposed heat flux in wall normal
                 ! ... direction flowing into the solid
                 rext_imp(counter) = boundary_value(i,1)
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO



  ! ... Construct the effective load vector
  !
  !
  !   reff =  (1/dt * C - (1-beta) * K) * Temp  + (1-beta) * R + beta * R
  !                                           k               k          k+1
  !
  !
  !
  !   where R   = rext_new (time varying load) + rext_imp (load from input file)
  !          k+1
  !

  reff = rint + rint2 +(1.0 - beta) * r_old + beta * (rext_new + rext_imp)

  ! ... the current load vector is stored in r_old for the calculation at
  ! ... the next timestep.
  r_old = (rext_new + rext_imp)

  ! Allocate shortened arrays
  ALLOCATE(newreff(1:nrows_ceff))
  ALLOCATE(newTemp(1:nrows_ceff))

  ! ... Initialize newTemp
  newTemp(:) = 0.0d0

  ! ... Remove rows from the effective load vector that are constrained by displacement BC's
  CALL RemoveBCHT_reff(nstart_ktc,LNumNp,reff,nrows_ceff,newreff)

  newreff(:) = newreff(:) - Ceff_fpTp(:)


  
  !
  !     Solve for nodal temperatures at time t + Dt
  !     
  !     C Temp    = R
  !           k+1    k+1
  !

  ! ... Solve using BlockSolve95 functions



  CALL BS95SOLVE(nrows_ceff,newreff,newTemp,contol,BS95debug)

!  CALL BS95FREE(BS95debug)


  ! ... Put newTemp back into Temp, the full temperature vector for this process
  ! ... without BCs missing
  CALL RemoveBCHT_newTemp(nstart_ktc,LNumNp,Temp,nrows_ceff,newTemp)

  ! Deallocate shortened arrays

  DEALLOCATE(newreff)
  DEALLOCATE(newTemp)




  IF(myid==0) PRINT*,'ASSEMBLING NEW TEMPERATURE VECTOR'


  
  ! ... Put the vectors back into their Rocfrac equivalent vectors

  DO i = 1, global%NumNp
     IF ( NodeProc(i) == myid ) THEN
        k = Local2Global(i) - nstart_ktc + 1
        global%Temperature( i ) = Temp( k )
     END IF
  END DO

  ! ... Communicate variables to other processors and assemble

  IF (nprocs > 1) THEN

     ! clear the status and request vars
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
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

     ! ... receive from other procs
     ! ... Allocate a buffer of size 2 times the number of nodes from
     ! ... processor i.  The buffer is orgnized as {node #, temp, node #, temp...}
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(frmproc(1:NumCommProcsFrom2))
     DO i = 1, NumCommProcsFrom2
        ALLOCATE(frmproc(i)%rcvbuf(1:2*NumCommNodesFrom2(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),2*NumCommNodesFrom2(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom2(i),10,ROCSTAR_COMMUNICATOR, &
             req_rcv(i),ierr)
     ENDDO

     ! ... Send to the other processes, the number of which is NumCommProcs2
     ! ... Allocate a buffer of size 2 times the number of nodes going to
     ! ... process i.  The buffer is orgnized as {node #, temp, node #, temp...}
     DO i = 1, NumCommProcs2
        ALLOCATE(bufsnd(1:2*NumCommNodes2(i)))
        counter = 0
        DO j = 1, NumCommNodes2(i)
           counter = counter + 1
           bufsnd(counter) = 1.0 * CommNodes2(i,j)
           ! ... k is the local node number that corresponds to the global node number
           ! ... that this process is sending to process i.  This processor has j nodes
           ! ... to send to process i.
           k = Global2Local(CommNodes2(i,j))
           counter = counter + 1
           bufsnd(counter) = global%Temperature( k )
        ENDDO
        CALL MPI_SEND(bufsnd(:),2*NumCommNodes2(i),MPI_DOUBLE_PRECISION, &
             CommProcs2(i),10,ROCSTAR_COMMUNICATOR,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)

     ! ... deallocate status and request vars
     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)

     ! ... assemble into Rocfrac nodal variable vectors
     ! ... loop through all processes from which this process has received information
     DO i = 1, NumCommProcsFrom2
        counter = 0
        ! ... loop through the number of nodes that were sent to this process
        ! ... process i
        DO j = 1, NumCommNodesFrom2(i)
           counter = counter + 1
           ! ... frmproc(i)%rcvbuf(1+2*(j-1)) is the global node number of the informaition
           ! ... frmproc(i)%rcvbuf(2*j)
           k = Global2Local( INT( frmproc(i)%rcvbuf(counter) ) )  ! local index of node
           counter = counter + 1
           global%Temperature( k ) = frmproc(i)%rcvbuf(counter)
        END DO
     END DO

     ! ... clear the receive buffer
     DEALLOCATE(frmproc)
     
     ! ... Deallocate effective capacitance matrices for the current 
     ! ... timestep.


  ENDIF


!  if (real(istep/10) == real(istep)/10 ) CALL RockOut_8node(global,istep)

END SUBROUTINE thermal_soln
