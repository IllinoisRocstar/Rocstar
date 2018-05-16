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
!!     fractography3d
!!
!!  FUNCTION
!!     This program inputs a model of a solids structure and
!!     performs an implicit time-dependant analysis of its
!!     motion based on the Newmark-Wilson scheme.
!!
!!  INPUTS
!!     File: fractography3d.inp
!!     File: timestep.in
!!     Files: <casename>/<casename>.<proc#>.inp
!!
!!  OUTPUTS
!!     Files: output_<proc#>.dat
!!
!!  USES
!!     MPI
!!     readpat
!!     feminp
!!     ReadMeshtran
!!     InitComm1
!!     InitComm2
!!     get_mat_stiffness
!!     implicit_v3d8_mass_consistent
!!     implicit_v3d8_mass
!!     implicit_v3d8_me_K
!!     comp_row_getval
!!     comp_row_addval
!!     comp_row_resize
!!     comp_row_add
!!     implicit_bc_enforce
!!     vol_elem_mat
!!     get_fint
!!     BS95SETUP
!!     BS95SOLVE
!!     BS95FREE
!!     BS95FINALIZE
!!     removebcs_meff
!!     removebcs_pbar
!!     removebcs_newa
!!     
!!
!!****

PROGRAM fractography3D

  USE Precision

  USE comp_row_global
  USE implicit_global
  
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  ! Declare types that you'll need
  TYPE int_buf
     INTEGER, DIMENSION(:), POINTER :: rcvbuf
  END TYPE int_buf

  ! Declare the variables that you'll need
  INTEGER :: ios  ! error control in opening file
  INTEGER :: n, nstep  ! current timestep & total timesteps
  INTEGER :: BS95debug = 0  ! debug statements in BlockSolve95
  REAL(kind=wp) :: delt  ! time step size
  REAL(kind=wp) :: thetaimp, alphaimp, deltaimp  ! implicit direct integration constants
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: fext  ! vector containing external forces
  CHARACTER*180 :: word180  ! dummy string for reading in input files
  INTEGER :: ifreq  ! frequency of output
  INTEGER :: nnz_k  ! number of nonzeros in the K matrix
  INTEGER :: nnz_m  ! number of nonzeros in the M matrix
  INTEGER :: i, j, m, p, counter  ! counters
  INTEGER :: ngpts = 8  ! number of gauss points per element
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: tempmg, mg  ! Global mass matrix (lumped)
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_m, aval_k  ! M and K matrices in compressed row storage
  INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_m, cval_m, rp_k, cval_k  ! M and K matrices in compressed row storage
  INTEGER :: nnz_m  ! Number of nonzeros in the mass matrix
  INTEGER :: nstart_km, nrows_km  ! Dimensions of parts of M  and K matrices
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_meff  ! Effective M matrix in compressed row storage
  INTEGER,ALLOCATABLE,DIMENSION(:) :: cval_meff, rp_meff  ! Effective M matrix in compressed row storage
  INTEGER :: nnz_meff, nstart_meff, nrows_meff  ! Dimensions of parts of effective M matrix
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: pbar  ! Effective load vector
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: newpbar  ! Pbar after removing DOF related to displacement BCs
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: fint  ! Internal force vector
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:) :: S  ! Stress used in calculation of internal force vector
  REAL(kind=wp), DIMENSION(1:3,1:8) :: ri = RESHAPE( &
       (/-0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626, 0.577350269189626, 0.577350269189626, &
         -0.577350269189626, 0.577350269189626, 0.577350269189626/),(/3,8/) )
  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: v  ! nodal velocity
  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: a  ! nodal acceleration
  REAL(kind=wp) :: t  ! current time within iteration
  REAL(kind=wp) :: maxdisp, gmaxdisp  ! maximum nodal displacement
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: ftemp  ! temporary vector used to calculate initial acceleration
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd  ! MPI send buffer
  INTEGER, ALLOCATABLE, DIMENSION(:) :: AmountToReceive, AmountToSend  ! Amount of stuff that will be received via MPI
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: finttemp, disptemp  ! Temporary arrays used to store info about all NumNp nodes
  INTEGER :: maxdispnode  ! Node that has the maximum displacement
  REAL(kind=wp) :: tempKval  ! temporary value from the K matrix
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_ktemp  ! Temporaray K matrix in compressed row storage
  INTEGER,ALLOCATABLE,DIMENSION(:) :: cval_ktemp, rp_ktemp  ! Temporaray K matrix in compressed row storage
  INTEGER :: nnz_ktemp  ! Number of nonzeros in the temporary K matrix
  REAL(kind=wp) :: contol  ! Convergence tolerance for BS95
  REAL(kind=wp), DIMENSION(5,9,9) :: dmat  ! Material stiffness matrix
  REAL(kind=wp), DIMENSION(2) :: props  ! Material properties needed to find material stiffness matrix
  REAL(kind=wp), DIMENSION(8) :: xi, eta, zeta  ! Gauss point info
  REAL(kind=wp), DIMENSION(8) :: xiE, etaE, zetaE  ! More Gauss point info
  REAL(kind=wp) :: one = 1.0, three = 3.0
  REAL(KIND=wp), DIMENSION(1:8,1:9,1:12) :: mixed_map  ! Mixed map tensor
  REAL(KIND=wp), DIMENSION(1:8,1:9,1:9)  :: enhanced_map  ! Enhanced map tensor
  INTEGER :: igpt  ! Gauss point counter
  REAL(kind=wp) :: tempval  ! value from the K matrix to be sent via MPI
  INTEGER :: idof, jdof, inode, jnode  ! Temporary numbering of nodes during K matrix communication
  REAL(kind=wp) :: per1  ! Period of the first mode
  INTEGER :: newnrows_meff  ! Number of rows in the Meff matrix after removing BC's
  INTEGER :: newndim  ! Number of DOF after removing DOF associated with imposed displacements
  INTEGER :: newnstart_meff  ! Starting row of the Meff matrix after removing BC's
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: newa  ! Acceleration vector after removing rows pertaining to displacement BCs
  TYPE(int_buf), ALLOCATABLE, DIMENSION(:) :: Ki, Kj  ! i and j indices of nonzeros in K matrix  
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_m_temp ! temporary portion of the M matrix using during construction
  INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_m_temp, cval_m_temp  ! temporary portion of the M matrix using during construction
  REAL(kind=8) :: elapsedtime


  ! Initialize MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  WRITE(myid_chr,'(i4.4)') myid

  ! Start timing of program
  elapsedtime = -1.0*MPI_WTIME()
  
!
!     Read geometry, material properties, and FE information
!

  ! Read in timestep information file
  OPEN(20,FILE='timestep.in' , FORM='formatted') !'
  READ(20,'(a180)') word180
  READ(20,*) nstep,delt,ifreq,contol  ! read in timestep info
  DO i=1,4
     READ(20,'(a180)') word180
  ENDDO
  READ(20,*) thetaimp,alphaimp,deltaimp  ! read in implicit direct integration constants
  CLOSE(20)

  !
  ! -- Open Analysis Deck File
  ScaleLoad = 1.0
  OPEN(io_input,FILE='fractography3d.inp',STATUS='old',IOSTAT=ios)
  CLOSE(io_input)
  IF(ios.NE.0)THEN
     IF(myid.EQ.0) PRINT*, '**WARNING: Unable to find Deck File -- fractography3d.inp'
     IF(myid.EQ.0) PRINT*, '        - Assuming Patran Neutral File Input Format'
     CALL ReadPat()

     LNumNp = NumNp
     GNumNp = NumNp
     ALLOCATE(Local2Global(1:NumNp))
     ALLOCATE(Global2Local(1:NumNp))
     ALLOCATE(NodeProc(1:NumNp))
     DO i = 1, NumNp
        Local2Global(i) = i
        Global2Local(i) = i
        NodeProc(i) = myid
     ENDDO

  ELSE
  !
  ! -- Read the control deck file
     IF(myid.EQ.0) PRINT*,'READING INPUT DECK'
     CALL feminp()
     CALL ReadMeshTran()

  ENDIF


!
! Initialize MPI communications
!
  IF(myid==0) PRINT*,'INITIALIZING COMMUNICATIONS'
  IF (nprocs > 1) THEN
     CALL InitComm1()
     CALL InitComm2()
  ENDIF

!
! Construct the material stiffness matrix for each material
!
  DO i = 1, NumMat
     props(1) = E(i)
     props(2) = xnu(i)
     CALL get_mat_stiffness(props,dmat(i,:,:),i)
  ENDDO


!
! Get Gauss point info
!
  xi   = (/ one/SQRT(three), -one/SQRT(three), -one/SQRT(three), &
       one/SQRT(three),  one/SQRT(three), -one/SQRT(three), &
       -one/SQRT(three),  one/SQRT(three)     /)
  
  eta  = (/ one/SQRT(three),  one/SQRT(three), -one/SQRT(three), &
       -one/SQRT(three),  one/SQRT(three),  one/SQRT(three), &
       -one/SQRT(three), -one/SQRT(three)     /)
  
  zeta = (/ one/SQRT(three),  one/SQRT(three),  one/SQRT(three), &
       one/SQRT(three), -one/SQRT(three), -one/SQRT(three), &  
       -one/SQRT(three), -one/SQRT(three)     /)

  xiE   = (/ SQRT(three), -SQRT(three), -SQRT(three), &
       SQRT(three),  SQRT(three), -SQRT(three), &
       -SQRT(three),  SQRT(three)     /)
  
  etaE  = (/ SQRT(three),  SQRT(three), -SQRT(three), &
       -SQRT(three),  SQRT(three),  SQRT(three), &
       -SQRT(three), -SQRT(three)     /)
  
  zetaE = (/ SQRT(three),  SQRT(three),  SQRT(three), &
       SQRT(three), -SQRT(three), -SQRT(three), &  
       -SQRT(three), -SQRT(three)     /)
  mixed_map(:,:,:) = 0.0
  enhanced_map(:,:,:) = 0.0
  DO igpt = 1, 8
     mixed_map(igpt,1,1)  = eta(igpt)
     mixed_map(igpt,1,2)  = zeta(igpt)
     mixed_map(igpt,1,3)  = eta(igpt) * zeta(igpt)
     mixed_map(igpt,2,10) = zeta(igpt)
     mixed_map(igpt,3,12) = eta(igpt)
     mixed_map(igpt,4,10) = zeta(igpt)
     mixed_map(igpt,5,4)  = xi(igpt)        
     mixed_map(igpt,5,5)  = zeta(igpt)
     mixed_map(igpt,5,6)  = xi(igpt) * zeta(igpt)
     mixed_map(igpt,6,11) = xi(igpt)
     mixed_map(igpt,7,12) = eta(igpt)
     mixed_map(igpt,8,11) = xi(igpt)
     mixed_map(igpt,9,7)  = xi(igpt)        
     mixed_map(igpt,9,8)  = eta(igpt)
     mixed_map(igpt,9,9)  = xi(igpt) * eta(igpt)
     enhanced_map(igpt,1,1)  = xi(igpt)    
     enhanced_map(igpt,1,2)  = xi(igpt) * eta(igpt)
     enhanced_map(igpt,1,3)  = xi(igpt) * zeta(igpt)
     enhanced_map(igpt,5,4)  = eta(igpt)
     enhanced_map(igpt,5,5)  = eta(igpt) * zeta(igpt)
     enhanced_map(igpt,5,6)  = eta(igpt) * xi(igpt)
     enhanced_map(igpt,9,7)  = zeta(igpt)        
     enhanced_map(igpt,9,8)  = zeta(igpt) * eta(igpt)
     enhanced_map(igpt,9,9)  = zeta(igpt) * xi(igpt)
  ENDDO
  ALLOCATE(Aenh(1:9,1:NumElv))
  ALLOCATE(stress(1:ngpts,1:9,1:NumElv))
  Aenh(:,:) = 0.0
  stress(:,:,:) = 0.0
  
  ! Set up some constants
  nstart_km = GNumNp
  DO m = 1, GNumNp
     DO i = 1, NumNp
        IF (Local2Global(i) == m) THEN
           IF (NodeProc(i) == myid) THEN
              nstart_km = MIN(nstart_km,Local2Global(i))
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  nstart_km = 3 * (nstart_km - 1) + 1
  nrows_km = 3*LNumNp

!
! 1. Form the stiffness, and mass matrices
!


  ! Construct the global mass matrix
  IF(myid==0) PRINT*,'CONSTRUCTING THE MASS MATRIX'

  ! Construct the mass matrix as a consistent matrix
  CALL implicit_v3d8_mass_consistent(NumElv,NumNp,NumMat,coor,ElConnVol,MatType,ri,rho,ElConnVol)
  nnz_m = nnz_temp
  ALLOCATE(cval_m(1:nnz_m))
  ALLOCATE(aval_m(1:nnz_m))
  ALLOCATE(rp_m(1:3*GNumNp+1))
  rp_m = rp_temp
  cval_m = cval_temp
  aval_m = aval_temp
  DEALLOCATE(rp_temp)
  DEALLOCATE(cval_temp)
  DEALLOCATE(aval_temp)

  ! Communicate mass to the other procs
  IF (nprocs > 1) THEN
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
     ALLOCATE(AmountToReceive(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        CALL MPI_IRECV(AmountToReceive(i),1, &
             MPI_INTEGER,CommProcsFrom1(i),10,MPI_COMM_WORLD, &
             req_rcv(i),ierr)
     ENDDO
     ALLOCATE(AmountToSend(1:NumCommProcs1))
     DO i = 1, NumCommProcs1
        counter = 0
        DO j = 1, NumCommNodes1(i)
           DO p = 1, 3
              DO m = 1, 3*GNumNp
                 tempKval = 0.0
                 inode = Local2Global(CommNodes1(i,j))
                 jnode = INT((m-0.5)/3)+1
                 idof = p
                 jdof = m - 3*jnode + 3
                 CALL comp_row_getval(3*GNumNp,3*GNumNp,nnz_m,1,rp_m,cval_m,aval_m, &
                      inode*3-3+idof,jnode*3-3+jdof,tempKval)
                 IF (tempKval /= 0.0) THEN
                    counter = counter + 1
                    !print*,myid,Global2Local(inode)*3-3+idof,Global2Local(jnode)*3-3+jdof,tempKval
                 ENDIF
              ENDDO
           ENDDO
           !print*,myid,' sending stuff about node ',Local2Global(CommNodes1(i,j)),' to ',CommProcs1(i),' : ',counter
        ENDDO
        AmountToSend(i) = counter
       !CALL MPI_ISEND(AmountToSend(i),1,MPI_INTEGER, &
       !      CommProcs1(i),10,MPI_COMM_WORLD,req_snd(i),ierr)
        CALL MPI_SEND(AmountToSend(i),1,MPI_INTEGER, &
             CommProcs1(i),10,MPI_COMM_WORLD,ierr)
     ENDDO
     !CALL MPI_WAITALL(NumCommProcsFrom1,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(NumCommProcs1,req_snd,stat_snd,ierr)
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

     ! Communicate parts of the M matrix to other procs
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     ALLOCATE(frmproc(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(frmproc(i)%rcvbuf(1:3*AmountToReceive(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),3*AmountToReceive(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom1(i),10,MPI_COMM_WORLD, &
             req_rcv(i),ierr)
     ENDDO
     DO i = 1, NumCommProcs1
        ALLOCATE(bufsnd(1:3*AmountToSend(i)))
        bufsnd(:) = 0.0
        counter = 0
        DO j = 1, NumCommNodes1(i)
           DO p = 1, 3
              DO m = 1, 3*GNumNp
                 tempKval = 0.0
                 inode = Local2Global(CommNodes1(i,j))
                 jnode = INT((m-0.5)/3)+1
                 idof = p
                 jdof = m - 3*jnode + 3
                 CALL comp_row_getval(3*GNumNp,3*GNumNp,nnz_m,1,rp_m,cval_m,aval_m, &
                      inode*3-3+idof,jnode*3-3+jdof,tempKval)
                 IF (tempKval /= 0.0) THEN
                    counter = counter + 1
                    bufsnd(3*counter-2) = Local2Global(CommNodes1(i,j))*3-3+p
                    bufsnd(3*counter-1) = m
                    bufsnd(3*counter) = tempKval
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !print*,myid,' sending to ',CommProcs1(i),' : ',bufsnd(:)
        !CALL MPI_ISEND(bufsnd,3*AmountToSend(i),MPI_DOUBLE_PRECISION, &
        !     CommProcs1(i),10,MPI_COMM_WORLD,req_snd(i),ierr)
        CALL MPI_SEND(bufsnd,3*AmountToSend(i),MPI_DOUBLE_PRECISION, &
             CommProcs1(i),10,MPI_COMM_WORLD,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     !CALL MPI_WAITALL(NumCommProcsFrom1,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(NumCommProcs1,req_snd,stat_snd,ierr)
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)
     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)
  ENDIF

  IF (nprocs > 1) THEN
     DO i = 1, NumCommProcsFrom1
        !print*,myid,' received from ',CommProcsFrom1(i), ' : ',frmproc(i)%rcvbuf(:)
        DO j = 1, AmountToReceive(i)
           CALL comp_row_addval(3*GNumNp,3*GNumNp,nnz_m,1,rp_m,cval_m,aval_m, &
                INT(frmproc(i)%rcvbuf(3*j-2)),INT(frmproc(i)%rcvbuf(3*j-1)),frmproc(i)%rcvbuf(3*j))
           IF (nnz_temp /= nnz_m) THEN
              DEALLOCATE(cval_m,aval_m)
              ALLOCATE(cval_m(nnz_temp),aval_m(nnz_temp))
              nnz_m = nnz_temp
              rp_m = rp_temp
              cval_m = cval_temp
           ENDIF
           aval_m = aval_temp
           DEALLOCATE(rp_temp,cval_temp,aval_temp)
        ENDDO
     ENDDO
     DEALLOCATE(frmproc)
     DEALLOCATE(AmountToReceive)
     DEALLOCATE(AmountToSend)
  ENDIF

  ! Resize the matrix to the size needed on this proc
  CALL comp_row_resize(3*GNumNp,3*GNumNp,nnz_m,1,rp_m,cval_m,aval_m,nstart_km,nrows_km)
  DEALLOCATE(rp_m,cval_m,aval_m)
  ALLOCATE(rp_m(1:nrows_km+1))
  ALLOCATE(cval_m(1:nnz_temp))
  ALLOCATE(aval_m(1:nnz_temp))
  nnz_m = nnz_temp
  rp_m = rp_temp
  cval_m = cval_temp
  aval_m = aval_temp
  DEALLOCATE(rp_temp,cval_temp,aval_temp)

  ! Construct the mass matrix as lumped
  ALLOCATE(tempmg(1:NumNp))
  ALLOCATE(mg(1:3*LNumNp))
  ALLOCATE(cval_m_temp(1:3*LNumNp))
  ALLOCATE(aval_m_temp(1:3*LNumNp))
  ALLOCATE(rp_m_temp(1:3*LNumNp+1))
  tempmg(1:NumNp) = 0.0
  CALL implicit_v3d8_mass(NumElv,NumNp,NumMat,coor,ElConnVol,MatType,ri,rho,tempmg,1,NumElv)
  DO i = 1, NumNp
     tempmg(i) = 1.0 / tempmg(i)  ! Get the mass matrix, not its inverse
  ENDDO

  IF (nprocs > 1) THEN
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !ALLOCATE(req_rcv(1:NumCommProcsFrom1))
     !ALLOCATE(req_snd(1:NumCommProcs1))
     !ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:NumCommProcs1))
     !ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:NumCommProcsFrom1))
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
     ALLOCATE(frmproc(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(frmproc(i)%rcvbuf(1:2*NumCommNodesFrom1(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),2*NumCommNodesFrom1(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom1(i),10,MPI_COMM_WORLD, &
             req_rcv(i),ierr)
     ENDDO
     DO i = 1, NumCommProcs1
        ALLOCATE(bufsnd(2*NumCommNodes1(i)))
        DO j = 1, NumCommNodes1(i)
           bufsnd(2*j-1) = Local2Global(CommNodes1(i,j))
           bufsnd(2*j) = tempmg(CommNodes1(i,j))
        ENDDO
        !CALL MPI_ISEND(bufsnd,2*NumCommNodes1(i),MPI_DOUBLE_PRECISION, &
        !     CommProcs1(i),10,MPI_COMM_WORLD,req_snd(i),ierr)
        CALL MPI_SEND(bufsnd,2*NumCommNodes1(i),MPI_DOUBLE_PRECISION, &
             CommProcs1(i),10,MPI_COMM_WORLD,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     !CALL MPI_WAITALL(NumCommProcsFrom1,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(NumCommProcs1,req_snd,stat_snd,ierr)
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)
     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)

     ! Add the communicated mass to the current mass
     DO i = 1, NumCommProcsFrom1
        DO j = 1, NumCommNodesFrom1(i)
           tempmg(Global2Local(INT(frmproc(i)%rcvbuf(2*j-1)))) = &
                tempmg(Global2Local(INT(frmproc(i)%rcvbuf(2*j-1)))) + frmproc(i)%rcvbuf(2*j)
        ENDDO
     ENDDO

     DEALLOCATE(frmproc)
  ENDIF

  ! Put the mass into compressed row format
  counter = 0
  DO m = 1, GNumNp
     DO i = 1, NumNp
        IF (Local2Global(i) == m) THEN
           IF (NodeProc(i) == myid) THEN
              DO j = 1, 3
                 counter = counter + 1
                 mg(counter) = tempmg(i)
              ENDDO
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  DO i = 1, 3*LNumNp
     rp_m_temp(i) = i-1
     cval_m_temp(i) = i-1 + (nstart_km - 1)
     aval_m_temp(i) = mg(i)
  ENDDO
  rp_m_temp(3*LNumNp+1) = 3*LNumNp
  DEALLOCATE(tempmg)
  DEALLOCATE(mg)

  ! Average the consistent and lumped mass matrices to get a higher-order mass matrix
  CALL comp_row_add(3*LNumNp,3*GNumNp,nrows_km,nrows_km,nnz_m,3*LNumNp,1, &
       rp_m,cval_m,0.5*aval_m,1,rp_m_temp,cval_m_temp,0.5*aval_m_temp)
  DEALLOCATE(rp_m,cval_m,aval_m)
  DEALLOCATE(rp_m_temp,cval_m_temp,aval_m_temp)
  ALLOCATE(rp_m(nrows_km+1),cval_m(nnz_temp),aval_m(nnz_temp))
  nnz_m = nnz_temp
  rp_m = rp_temp
  cval_m = cval_temp
  aval_m = aval_temp
  DEALLOCATE(rp_temp,cval_temp,aval_temp)

  !! Output the matrix to matlab files for analysis
  !OPEN(30,FILE='m_'//myid_chr//'.m',FORM='formatted')
  !counter = 0
  !DO i = 1, nrows_km
  !   DO j = rp_m(i)+1, rp_m(i+1)
  !      counter = counter + 1
  !      m = cval_m(counter)+1
  !      WRITE(30,*) 'A(', i+nstart_km-1, ',', m, ') = ',aval_m(counter), ';'
  !   ENDDO
  !ENDDO
  !CLOSE(30)

     
  ! Construct the global stiffness matrix
  IF(myid==0) PRINT*,'CONSTRUCTING THE STIFFNESS MATRIX'
  CALL implicit_v3d8_me_K(coor, ElConnVol, dmat, NumNp, 1, NumElv, NumElv, &
       MatType, NumMat, enhanced_map, mixed_map)
  ALLOCATE(rp_k(1:3*GNumNp+1))
  ALLOCATE(cval_k(1:nnz_temp))
  ALLOCATE(aval_k(1:nnz_temp))
  nnz_k = nnz_temp
  rp_k = rp_temp
  cval_k = cval_temp
  aval_k = aval_temp
  DEALLOCATE(rp_temp,cval_temp,aval_temp)

  ! Communicate the K matrix to other procs
  IF (nprocs > 1) THEN
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
     ALLOCATE(AmountToReceive(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        CALL MPI_IRECV(AmountToReceive(i),1, &
             MPI_INTEGER,CommProcsFrom1(i),10,MPI_COMM_WORLD, &
             req_rcv(i),ierr)
     ENDDO
     ALLOCATE(AmountToSend(1:NumCommProcs1))
     DO i = 1, NumCommProcs1
        counter = 0
        DO j = 1, NumCommNodes1(i)
           DO p = 1, 3
              DO m = 1, 3*GNumNp
                 tempKval = 0.0
                 inode = Local2Global(CommNodes1(i,j))
                 jnode = INT((m-0.5)/3)+1
                 idof = p
                 jdof = m - 3*jnode + 3
                 CALL comp_row_getval(3*GNumNp,3*GNumNp,nnz_k,1,rp_k,cval_k,aval_k, &
                      inode*3-3+idof,jnode*3-3+jdof,tempKval)
                 IF (tempKval /= 0.0) THEN
                    counter = counter + 1
                    !print*,myid,Global2Local(inode)*3-3+idof,Global2Local(jnode)*3-3+jdof,tempKval
                 ENDIF
              ENDDO
           ENDDO
           !print*,myid,' sending stuff about node ',Local2Global(CommNodes1(i,j)), &
           !' to ',CommProcs1(i),' : ',counter
        ENDDO
        AmountToSend(i) = counter
        !CALL MPI_ISEND(AmountToSend(i),1,MPI_INTEGER, &
        !     CommProcs1(i),10,MPI_COMM_WORLD,req_snd(i),ierr)
        CALL MPI_SEND(AmountToSend(i),1,MPI_INTEGER, &
             CommProcs1(i),10,MPI_COMM_WORLD,ierr)
     ENDDO
     !CALL MPI_WAITALL(NumCommProcsFrom1,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(NumCommProcs1,req_snd,stat_snd,ierr)
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)


     ! Communicate i locations of nonzeros in K matrix to other procs
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     ALLOCATE(Ki(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(Ki(i)%rcvbuf(1:AmountToReceive(i)))
        CALL MPI_IRECV(Ki(i)%rcvbuf(1),AmountToReceive(i), &
             MPI_INTEGER,CommProcsFrom1(i),10,MPI_COMM_WORLD, &
             req_rcv(i),ierr)
     ENDDO
     DO i = 1, NumCommProcs1
        ALLOCATE(bufsnd(1:AmountToSend(i)))
        counter = 0
        DO j = 1, NumCommNodes1(i)
           DO p = 1, 3
              DO m = 1, 3*GNumNp
                 tempKval = 0.0
                 inode = Local2Global(CommNodes1(i,j))
                 jnode = INT((m-0.5)/3)+1
                 idof = p
                 jdof = m - 3*jnode + 3
                 CALL comp_row_getval(3*GNumNp,3*GNumNp,nnz_k,1,rp_k,cval_k,aval_k, &
                      inode*3-3+idof,jnode*3-3+jdof,tempKval)
                 IF (tempKval /= 0.0) THEN
                    counter = counter + 1
                    bufsnd(counter) = Local2Global(CommNodes1(i,j))*3-3+p
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !CALL MPI_ISEND(INT(bufsnd(:)),AmountToSend(i),MPI_INTEGER, &
        !     CommProcs1(i),10,MPI_COMM_WORLD,req_snd(i),ierr)
        CALL MPI_SEND(INT(bufsnd(:)),AmountToSend(i),MPI_INTEGER, &
             CommProcs1(i),10,MPI_COMM_WORLD,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

     ! Communicate j locations of nonzeros in K matrix to other procs
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     ALLOCATE(Kj(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(Kj(i)%rcvbuf(1:AmountToReceive(i)))
        CALL MPI_IRECV(Kj(i)%rcvbuf(1),AmountToReceive(i), &
             MPI_INTEGER,CommProcsFrom1(i),10,MPI_COMM_WORLD, &
             req_rcv(i),ierr)
     ENDDO
     DO i = 1, NumCommProcs1
        ALLOCATE(bufsnd(1:AmountToSend(i)))
        counter = 0
        DO j = 1, NumCommNodes1(i)
           DO p = 1, 3
              DO m = 1, 3*GNumNp
                 tempKval = 0.0
                 inode = Local2Global(CommNodes1(i,j))
                 jnode = INT((m-0.5)/3)+1
                 idof = p
                 jdof = m - 3*jnode + 3
                 CALL comp_row_getval(3*GNumNp,3*GNumNp,nnz_k,1,rp_k,cval_k,aval_k, &
                      inode*3-3+idof,jnode*3-3+jdof,tempKval)
                 IF (tempKval /= 0.0) THEN
                    counter = counter + 1
                    bufsnd(counter) = m
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !CALL MPI_ISEND(INT(bufsnd(:)),AmountToSend(i),MPI_INTEGER, &
        !     CommProcs1(i),10,MPI_COMM_WORLD,req_snd(i),ierr)
        CALL MPI_SEND(INT(bufsnd(:)),AmountToSend(i),MPI_INTEGER, &
             CommProcs1(i),10,MPI_COMM_WORLD,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

     ! Communicate values of nonzeros in K matrix to other procs
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     ALLOCATE(frmproc(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(frmproc(i)%rcvbuf(1:AmountToReceive(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),AmountToReceive(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom1(i),10,MPI_COMM_WORLD, &
             req_rcv(i),ierr)
     ENDDO
     DO i = 1, NumCommProcs1
        ALLOCATE(bufsnd(1:AmountToSend(i)))
        counter = 0
        DO j = 1, NumCommNodes1(i)
           DO p = 1, 3
              DO m = 1, 3*GNumNp
                 tempKval = 0.0
                 inode = Local2Global(CommNodes1(i,j))
                 jnode = INT((m-0.5)/3)+1
                 idof = p
                 jdof = m - 3*jnode + 3
                 CALL comp_row_getval(3*GNumNp,3*GNumNp,nnz_k,1,rp_k,cval_k,aval_k, &
                      inode*3-3+idof,jnode*3-3+jdof,tempKval)
                 IF (tempKval /= 0.0) THEN
                    counter = counter + 1
                    bufsnd(counter) = tempKval
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !CALL MPI_ISEND(bufsnd,AmountToSend(i),MPI_DOUBLE_PRECISION, &
        !     CommProcs1(i),10,MPI_COMM_WORLD,req_snd(i),ierr)
        CALL MPI_SEND(bufsnd,AmountToSend(i),MPI_DOUBLE_PRECISION, &
             CommProcs1(i),10,MPI_COMM_WORLD,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)

  ENDIF
  
  CALL comp_row_resize(3*GNumNp,3*GNumNp,nnz_k,1,rp_k,cval_k,aval_k,nstart_km,nrows_km)
  DEALLOCATE(rp_k,cval_k,aval_k)
  ALLOCATE(rp_k(1:nrows_km+1))
  ALLOCATE(cval_k(1:nnz_temp))
  ALLOCATE(aval_k(1:nnz_temp))
  nnz_k = nnz_temp
  rp_k = rp_temp
  cval_k = cval_temp
  aval_k = aval_temp
  DEALLOCATE(rp_temp,cval_temp,aval_temp)


  IF (nprocs > 1) THEN
     DO i = 1, NumCommProcsFrom1
        !print*,myid,' received from ',CommProcsFrom1(i), ' : ',frmproc(i)%rcvbuf(:)
        DO j = 1, AmountToReceive(i)
           CALL comp_row_addval(3*GNumNp,nrows_km,nnz_k,nstart_km,rp_k,cval_k,aval_k, &
              Ki(i)%rcvbuf(j),Kj(i)%rcvbuf(j),frmproc(i)%rcvbuf(j))
           IF (nnz_temp /= nnz_k) THEN
              DEALLOCATE(cval_k,aval_k)
              ALLOCATE(cval_k(nnz_temp),aval_k(nnz_temp))
              nnz_k = nnz_temp
              rp_k = rp_temp
           ENDIF
           cval_k = cval_temp
           aval_k = aval_temp
           DEALLOCATE(rp_temp,cval_temp,aval_temp)
        ENDDO
     ENDDO
     DEALLOCATE(frmproc)
     DEALLOCATE(Ki)
     DEALLOCATE(Kj)
  ENDIF 


!                                                   .
! 2. Initialize displacement and velocity array: Uo,Uo
!
  IF(myid==0) PRINT*,'INITIALIZING DISPLACEMENT AND VELOCITY'
  ALLOCATE(disp(1:3*LNumNp),v(1:3*LNumNp),a(1:3*LNumNp))
  disp(1:3*LNumNp) = 0.0
  v(1:3*LNumNp) = 0.0
  IF (ALLOCATED(node_flag)) THEN
  ELSE
     ALLOCATE(node_flag(1:NumNp,1:3))
     ALLOCATE(boundary_value(1:NumNp,1:3))
     node_flag(1:NumNp,1:3) = 0
     DO i = 1, numbound
        DO j = 1, 3
!!$           IF (BC_ID(j+1,i) == 0) THEN ! displacement BC
!!$              node_flag(BC_ID(1,i),j) = 8
!!$              boundary_value(BC_ID(1,i),j) = BC_VAL(j,i)
           IF ( global%BCFlag(j+1,i) == 0 ) THEN
              node_flag(global%BCFlag(1,i),j) = 8
              boundary_value(global%BCFlag(1,i),j) = global%BCvalue(j,i)
           ENDIF
!!$           IF (BC_ID(j+1,i) == 1) THEN ! force BC
!!$              node_flag(BC_ID(1,i),j) = 7
!!$              boundary_value(BC_ID(1,i),j) = BC_VAL(j,i)
           IF ( global%BCFlag(j+1,i) == 1 ) THEN
              node_flag(global%BCFlag(1,i),j) = 7
              boundary_value(global%BCFlag(1,i),j) = global%BCvalue(j,i)
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  CALL implicit_bc_enforce(NumNp,LNumNp,disp,a,v,node_flag,boundary_value,0,myid)


  !
  !     Form the material compliance matrix
  !

  IF(ALLOCATED(c)) DEALLOCATE(c)
  IF(ALLOCATED(ci)) DEALLOCATE(ci)
  ALLOCATE(c(1:9,1:NumMat),ci(1:9,1:NumMat))
  c(1:9,1:NumMat) = 0.0
  ci(1:9,1:NumMat) = 0.0
  CALL vol_elem_mat(E,xnu,c,ci,NumMat)


!                                              ..
! 3. Calculate the initial acceleration array: Uo
!


  ! Use applied forces and mass matrix
  ! i.e. solve {F} - [K]*{disp} = [M]*{a} for {a}

  IF(myid==0) PRINT*,'CALCULATING INITIAL ACCELERATION'

  ALLOCATE(fint(1:3*LNumNp))
  CALL get_fint(3*GNumNp,nrows_km,nnz_k,nstart_km,rp_k,cval_k,aval_k,disp,fint)
  
  ! Construct the external force vector
  ALLOCATE(fext(1:3*LNumNp))
  fext(1:3*LNumNp) = 0.0
  counter = 0
  DO m = 1, GNumNp
     DO i = 1, NumNp
        IF (Local2Global(i) == m) THEN
           DO j = 1, 3
              IF (NodeProc(i) == myid) THEN
                 counter = counter + 1
                 IF(node_flag(i,j) == 7) THEN  ! Imposed force
                    fext(counter) = boundary_value(i,j)
                    !print*,myid,' external force on global node ',Local2Global(i),' in direction ',j,' of value ',fext(counter)
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  ALLOCATE(ftemp(1:3*LNumNp))
  ftemp(1:3*LNumNp) = 0.0


  ! Construct the effective force vector
  DO i = 1, 3*LNumNp
     ftemp(i) = fext(i) - fint(i)
     !print*,myid,i,ftemp(i)
  ENDDO
  DEALLOCATE(fint)

  ! Remove rows that are constrained by displacement BC's
  counter = 0
  DO m = 1, GNumNp
     DO i = 1, NumNp
        IF (Local2Global(i) == m) THEN
           IF (NodeProc(i)==myid) THEN
              DO j = 1, 3
                 counter = counter + 1
                 IF(node_flag(i,j) == 8) THEN  ! Imposed constant nodal displacement
                    DO n = rp_m(counter)+1,rp_m(counter+1)
                       IF (cval_m(n)+1 /= counter+nstart_km-1) THEN
                          aval_m(n) = 0.0
                       ELSE
                          aval_m(n) = 1.0
                       ENDIF
                    ENDDO
                 ENDIF
              ENDDO
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  counter = 0
  DO m = 1, GNumNp
     DO i = 1, NumNp
        IF (Local2Global(i) == m) THEN
           IF (NodeProc(i)==myid) THEN
              DO j = 1, 3
                 counter = counter + 1
                 IF(node_flag(i,j) == 8) THEN  ! Imposed constant nodal displacement
                    ftemp(counter) = 0.0  ! If constant displacement, acceleration = 0
                 ENDIF
              ENDDO
           ENDIF
        ENDIF
     ENDDO
  ENDDO

  
!
!     3.1 Solve for acceleration array
!
  ! Solve using BlockSolve95 functions
  a(:) = 0.0  ! Initial guess
  CALL BS95SETUP(3*GNumNp,nnz_m,nstart_km-1,nrows_km,rp_m,cval_m,aval_m,1,BS95debug)
  CALL BS95SOLVE(3*LNumNp,ftemp,a,contol,BS95debug)
  CALL BS95FREE(BS95debug)
  DEALLOCATE(ftemp)


!
!     Make sure initial conditions conform to BC's
!
  CALL implicit_bc_enforce(NumNp,LNumNp,disp,a,v,node_flag,boundary_value,0,myid)


!                                _  _
!     Form effective mass matrix M: M = M + gamma*Dt*C + beta*Dt^2*K
!

  IF(myid==0) PRINT*,'FORMING EFFECTIVE MASS MATRIX'
  CALL comp_row_add(3*LNumNp,3*GNumNp,nrows_km,nrows_km,nnz_k,nnz_m,1, &
       rp_k,cval_k,alphaimp*delt*delt*aval_k,1,rp_m,cval_m,aval_m)
  ALLOCATE(rp_meff(nrows_km+1),cval_meff(nnz_temp),aval_meff(nnz_temp))
  nnz_meff = nnz_temp
  nstart_meff = nstart_km
  nrows_meff = nrows_km
  rp_meff = rp_temp
  cval_meff = cval_temp
  aval_meff = aval_temp
  DEALLOCATE(rp_temp,cval_temp,aval_temp)

  ! Remove rows that are constrained by displacement BC's
  CALL removebcs_meff(3*GNumNp,nrows_meff,nnz_meff,nstart_meff,rp_meff,cval_meff,aval_meff, &
       newnrows_meff,newnstart_meff,newndim)
  DEALLOCATE(rp_meff,cval_meff,aval_meff)
  nrows_meff = newnrows_meff
  nstart_meff = newnstart_meff
  ALLOCATE(rp_meff(nrows_meff+1),cval_meff(nnz_temp),aval_meff(nnz_temp))
  nnz_meff = nnz_temp
  rp_meff = rp_temp
  cval_meff = cval_temp
  aval_meff = aval_temp
  DEALLOCATE(rp_temp,cval_temp,aval_temp)

  !! Output matrix to matlab files for analysis
  !OPEN(30,FILE='meff_'//myid_chr//'.m',FORM='formatted')
  !counter = 0
  !DO i = 1, nrows_meff
  !   DO j = rp_meff(i)+1, rp_meff(i+1)
  !      counter = counter + 1
  !      m = cval_meff(counter)+1
  !      WRITE(30,*) 'A(', i+nstart_meff-1, ',', m, ') = ',aval_meff(counter), ';'
  !   ENDDO
  !ENDDO
  !CLOSE(30)
  !OPEN(30,FILE='k_'//myid_chr//'.m',FORM='formatted')
  !counter = 0
  !DO i = 1, nrows_km
  !   DO j = rp_k(i)+1, rp_k(i+1)
  !      counter = counter + 1
  !      m = cval_k(counter)+1
  !      WRITE(30,*) 'A(', i+nstart_km-1, ',', m, ') = ',aval_k(counter), ';'
  !   ENDDO
  !ENDDO
  !CLOSE(30)

  CALL BS95SETUP(newndim,nnz_meff,nstart_meff-1,nrows_meff,rp_meff,cval_meff,aval_meff,1,BS95debug)


!
!     Allocate variables for use within the time loop
!
  ALLOCATE(pbar(1:3*LNumNp))
  ALLOCATE(fint(1:3*LNumNp))
  ALLOCATE(S(1:6,1:ngpts,1:NumElv))
  ALLOCATE(newpbar(1:nrows_meff))
  ALLOCATE(newa(1:nrows_meff))
  newa(:) = 0.0  ! Initial guess


!
!     Open output files
!
   OPEN(30,FILE='output_'//myid_chr//'.dat',FORM='formatted')

   !DO i = 1, NumNp
   !   IF(NodeProc(i) == myid) THEN
   !      counter = counter + 1
   !      !IF(Local2Global(i) == 41)THEN
   !      !IF(Local2Global(i) == 88) THEN
   !      IF(Local2Global(i) == 718) THEN
   !         OPEN(40,FILE='waveprop_'//myid_chr//'.dat',FORM='formatted')
   !         EXIT
   !      ENDIF
   !   ENDIF
   !ENDDO
   
  

!
!     For each time step

   IF (myid==0) THEN
      PRINT*,'BEGINNING TIMESTEPPING'
      PRINT*,' '
      PRINT*,'  step         time         maxdisp'
      PRINT*,'------------------------------------'
   ENDIF
   t = 0.0
   DO n=1,nstep

      t = t + delt
      !IF(myid==0) print*,'timestep ',n,' of ',nstep,'        delt = ',delt,'        t = ',t


!     Calculate displacement and velocity predictors:
!     ~                 .                              ..
!     x    =  x  + Dt * x  + 1/2 * Dt^2 * (1-2*beta) * x
!      k+1     k         k                              k
!
!     ~.      .                     ..
!     x    =  x  + Dt * (1-gamma) * x
!      k+1     k                     k
!

     DO i = 1, 3*LNumNp
        disp(i) = disp(i) + delt * v(i) + 0.5*delt*delt * (1.0 - 2.0*alphaimp) * a(i)
        v(i) = v(i) + delt * (1.0 - deltaimp) * a(i)
     ENDDO

     
!     Calculate effective loads at time t + Dt:
!
!     _              .
!     P    = P   - C x - I 
!      k+1    k+1         k
!
!     where I is the internal force vector
!
     ! Get the internal force vector
     CALL get_fint(3*GNumNp,nrows_km,nnz_k,nstart_km,rp_k,cval_k,aval_k,disp,fint)

     ! Construct the load vector
     ! BEGIN TIME-DEPENDENT FORCE CODE MODIFICATIONS
     !per1 = 359.159245157598 
     DO i = 1,3*LNumNp
        !if(t < per1) then
        !   pbar(i) = fext(i)*SIN(3.14159*t/per1) - fint(i)
        !else
        !   pbar(i) = -fint(i)
        !endif
        pbar(i) = fext(i) - fint(i)
     ENDDO 
     ! END TIME-DEPENDENT FORCE CODE MODIFICATIONS

     ! Remove rows that are constrained by displacement BC's
     CALL removeBCs_pbar(nstart_km,3*LNumNp,pbar,nrows_meff,newpbar)
     
     
!
!     Solve for acceleration at time t + Dt
!     _        _
!     M a    = P
!        k+1    k+1
!
     ! Solve using BlockSolve95 functions
     !newa(:) = 0.0
     CALL BS95SOLVE(nrows_meff,newpbar,newa,contol,BS95debug)
!!$     CALL BS95SOLVE(nrows_meff,pbar,a,contol,BS95debug)


     ! Put newa back into a
     CALL removeBCs_newa(nstart_km,3*LNumNp,a,nrows_meff,newa)


!     Calculate displacement and velocity correctors:
!            ~                    ..
!     x    = x    + beta * Dt^2 * x
!      k+1    k+1                  k+1
!
!     .      ~.                  ..
!     x    = x    + gamma * Dt * x
!      k+1    k+1                 k+1
!
     DO i = 1, 3*LNumNp
        disp(i) = disp(i) + alphaimp*delt*delt*a(i)
        v(i) = v(i) + deltaimp*delt*a(i)
     ENDDO


!     
!     4. Apply boundary conditions and output results
!
     CALL implicit_bc_enforce(NumNp,LNumNp,disp,v,a,node_flag,boundary_value,t,myid)


!
!     Compute maximum nodal displacement
!
     maxdispnode = 0
     maxdisp = 0.0
     i = 0
     DO j = 1, NumNp
        IF (NodeProc(j) == myid) THEN
           i = i + 1
           IF ( maxdisp < SQRT(disp(3*(i-1)+1)**2.0 + disp(3*(i-1)+2)**2.0 + disp(3*(i-1)+3)**2.0) ) maxdispnode = j
           maxdisp = MAX(maxdisp,SQRT(disp(3*(i-1)+1)**2.0 + disp(3*(i-1)+2)**2.0 + disp(3*(i-1)+3)**2.0))
        ENDIF
     ENDDO   


!
!     Output results to a file
!
     IF(MOD(n,ifreq).EQ.0)THEN
        CALL MPI_REDUCE(maxdisp, gmaxdisp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
        IF(myid==0) WRITE(*,3000) n, t, gmaxdisp
        DO i=1,3*LNumNp
           WRITE(30,1001) Local2Global((i-1)/3+1),MOD(i-1,3)+1,disp(i),v(i),a(i),t
        ENDDO
        WRITE(30,'()')
        !counter = 0
        !DO i = 1, NumNp
        !   IF(NodeProc(i) == myid) THEN
        !      counter = counter + 1
        !      !IF(Local2Global(i) == 41)THEN
        !      !IF(Local2Global(i) == 88) THEN
        !      IF(Local2Global(i) == 718) THEN
        !         IF(t<per1)THEN
        !            WRITE(40,*) t, disp(3*counter-2),disp(3*counter-1),disp(3*counter), SIN(3.14159*t/per1)
        !         ELSE
        !            WRITE(40,*) t, disp(3*counter-2),disp(3*counter-1),disp(3*counter), 0.0
        !         ENDIF
        !         EXIT
        !      ENDIF
        !   ENDIF
        !ENDDO
     ENDIF

  ENDDO

  ! Stop timing of program
  elapsedtime = elapsedtime + MPI_WTIME()
  if(myid==0) print*,'Elapsed time of program:  ',elapsedtime,' seconds'

  ! Finalize BS95
  CALL BS95FREE(BS95debug)
  CALL BS95FINALIZE(BS95debug)

  ! Close output files
  CLOSE(30)
  !CLOSE(40)

  ! Deallocate arrays
  DEALLOCATE(cval_m)
  DEALLOCATE(aval_m)
  DEALLOCATE(rp_m)
  DEALLOCATE(cval_k)
  DEALLOCATE(aval_k)
  DEALLOCATE(rp_k)
  DEALLOCATE(cval_meff)
  DEALLOCATE(aval_meff)
  DEALLOCATE(rp_meff)
  DEALLOCATE(fext)
  DEALLOCATE(pbar)
  DEALLOCATE(fint)
  DEALLOCATE(disp)
  DEALLOCATE(v)
  DEALLOCATE(a)
  DEALLOCATE(S)
  DEALLOCATE(c)
  DEALLOCATE(ci)
  DEALLOCATE(Local2Global)
  DEALLOCATE(Global2Local)
  DEALLOCATE(NodeProc)
  IF (nprocs > 1) THEN
     DEALLOCATE(CommProcs1)
     DEALLOCATE(NumCommNodes1)
     DEALLOCATE(CommNodes1)
     DEALLOCATE(CommProcsFrom1)
     DEALLOCATE(NumCommNodesFrom1)
     DEALLOCATE(CommProcs2)
     DEALLOCATE(NumCommNodes2)
     DEALLOCATE(CommNodes2)
     DEALLOCATE(CommProcsFrom2)
     DEALLOCATE(NumCommNodesFrom2)
  ENDIF

  ! Formats
  1001 FORMAT(i5,' ',i1,' ',4(f15.6,' '))
  2000 FORMAT(100(f6.2,' '))
  3000 FORMAT('  ',i5,'     ',f12.3,'     ',f10.6)

END PROGRAM fractography3D

