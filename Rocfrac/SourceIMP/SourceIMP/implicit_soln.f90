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
SUBROUTINE implicit_soln(delt,t,fext_in,global)

  USE implicit_global
  USE comp_row_global
  USE Precision
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

! Input
  TYPE(ROCFRAC_GLOBAL) :: global
  REAL(kind=wp) :: delt, t
  REAL(kind=wp), DIMENSION(1:3*global%NumNp) :: fext_in
! Internal
  INTEGER :: BS95debug
  INTEGER :: i, j, k, m, n, p, counter, ii
  INTEGER :: newnrows_meff, newnstart_meff, newndim
  REAL(kind=wp) :: contol, alphaimp, deltaimp
  REAL(kind=wp), DIMENSION(1:3*LNumNp) :: pbar, fint, fext, v, a, d
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: newpbar, newa
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd, ftemp
  REAL(wp) :: per1, delay


! Initialize
  fint(:) = 0.0
  contol = 1.0e-09
  alphaimp = 0.25
  deltaimp = 0.5
  IF ( global%debug_state ) THEN
     BS95debug = 1
  ELSE
     BS95debug = 0
  END IF


! Assemble the vectors
  DO i = 1, global%NumNp
     IF ( NodeProc(i) == myid ) THEN
        DO j = 1, 3
           k = 3*Local2Global(i) - nstart_km + 1 + j - 3
           m = 3*i + j - 3
           fext( k ) = fext_in( m )
           d( k ) = global%Disp( m )
           v( k ) = global%VeloHalf( m )
           a( k ) = global%Accel( m )
        END DO
     END IF
  END DO


!
! If this is the first timestep, initialize the acceleration array
!
  IF ( initAccel .EQV. .TRUE. ) THEN

     initAccel = .FALSE.

     ! Use applied forces and mass matrix
     ! i.e. solve {F} - [K]*{disp} = [M]*{a} for {a}

     IF(myid==0) PRINT*,'CALCULATING INITIAL ACCELERATION'

     ! Construct the imposed external force vector
     ALLOCATE(fext_imp(1:3*LNumNP))
     fext_imp(:) = 0.0
     counter = 0
     DO m = 1, GNumNp
        DO i = 1, global%NumNp
           IF (Local2Global(i) == m) THEN
              DO j = 1, 3
                 IF (NodeProc(i) == myid) THEN
                    counter = counter + 1
                    IF(node_flag(i,j) == 7) THEN  ! Imposed force
                       fext_imp(counter) = boundary_value(i,j)
                           
                    ENDIF
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
     ENDDO

     ! Construct the internal force vector
     CALL get_fint(3*GNumNp,nrows_km,nnz_k,nstart_km,rp_k,cval_k,aval_k,d,fint)

     ALLOCATE(ftemp(1:3*LNumNp))
     ftemp(1:3*LNumNp) = 0.0
     
     ! Construct the effective force vector
     DO i = 1, 3*LNumNp
        ftemp(i) = fext(i) + fext_imp(i) - fint(i)
     ENDDO
     
     ! Remove rows that are constrained by displacement BC's
     counter = 0
     DO m = 1, GNumNp
        DO i = 1, global%NumNp
           IF (Local2Global(i) == m) THEN
              IF (NodeProc(i) == myid) THEN
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
        DO i = 1, global%NumNp
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
     
     ! Solve using BlockSolve95 functions
     a(:) = 0.0  ! Initial guess
     CALL BS95SETUP(3*GNumNp,nnz_m,nstart_km-1,nrows_km,rp_m,cval_m,aval_m,1,BS95debug)
     CALL BS95SOLVE(3*LNumNp,ftemp,a,contol,BS95debug)
     CALL BS95FREE(BS95debug)
     DEALLOCATE(ftemp)

     ! Make sure initial conditions conform to BC's
     CALL implicit_bc_enforce(global%NumNp,LNumNp,d,a,v,node_flag,boundary_value,0,myid)

     ! Create the effective mass matrix
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
          newnrows_meff,newnstart_meff,newndim,global)
     DEALLOCATE(rp_meff,cval_meff,aval_meff)
     nrows_meff = newnrows_meff
     nstart_meff = newnstart_meff
     ALLOCATE(rp_meff(nrows_meff+1),cval_meff(nnz_temp),aval_meff(nnz_temp))
     nnz_meff = nnz_temp
     rp_meff = rp_temp
     cval_meff = cval_temp
     aval_meff = aval_temp
     DEALLOCATE(rp_temp,cval_temp,aval_temp)

     ! Set up the new system to be solved
     CALL BS95SETUP(newndim,nnz_meff,nstart_meff-1,nrows_meff,rp_meff,cval_meff,aval_meff,1, &
          BS95debug)

  END IF


!!$! DEBUG
!!$
!!$  IF ( t == 0.0 ) THEN
!!$
!!$     OPEN(30,FILE='m_'//global%MyIdChr//'.m',FORM='formatted')
!!$     counter = 0
!!$     DO i = 1, nrows_km
!!$        DO j = rp_m(i)+1, rp_m(i+1)
!!$           counter = counter + 1
!!$           m = cval_m(counter)+1
!!$           WRITE(30,*) 'M(', i+nstart_km-1, ',', m, ') = ',aval_m(counter), ';'
!!$        ENDDO
!!$     ENDDO
!!$     CLOSE(30)
!!$     
!!$     OPEN(30,FILE='k_'//global%MyIdChr//'.m',FORM='formatted')
!!$     counter = 0
!!$     DO i = 1, nrows_km
!!$        DO j = rp_k(i)+1, rp_k(i+1)
!!$           counter = counter + 1
!!$           m = cval_k(counter)+1
!!$           WRITE(30,*) 'K(', i+nstart_km-1, ',', m, ') = ',aval_k(counter), ';'
!!$        ENDDO
!!$     ENDDO
!!$     CLOSE(30)
!!$
!!$     OPEN(30,FILE='meff_'//global%MyIdChr//'.m',FORM='formatted')
!!$     counter = 0
!!$     DO i = 1, nrows_meff
!!$        DO j = rp_meff(i)+1, rp_meff(i+1)
!!$           counter = counter + 1
!!$           m = cval_meff(counter)+1
!!$           WRITE(30,*) 'Meff(', i+nstart_meff-1, ',', m, ') = ',aval_meff(counter), ';'
!!$        ENDDO
!!$     ENDDO
!!$     CLOSE(30)
!!$
!!$  END IF
!!$     
!!$! END DEBUG



  IF(myid==0) PRINT*,'CALCULATING NEW DISPLACEMENT VECTOR'


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
     d(i) = d(i) + delt * v(i) + 0.5*delt*delt * (1.0 - 2.0*alphaimp) * a(i)
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
  fint(:) = 0.0
  CALL get_fint(3*GNumNp,nrows_km,nnz_k,nstart_km,rp_k,cval_k,aval_k,d,fint)

  ! Construct the load vector
  !per1 = 359.1592  ! beambending2
  !delay = 0.0
  !print*,'constructing load vector for beambending2 model'
  per1 = 0.03237  ! agard (Yates SOLID 1 flutter)
  delay = 0.0 * per1
  print*,'constructing load vector for YATES SOLID 1 model in flutter'
  !per1 = 100.0
  !delay = 0.0
  !print*,'constructing load vector for orthotropic test cases'
  DO i = 1,3*LNumNp
 !    IF ( ( t > delay ) .AND. ( t < per1 + delay ) ) THEN
 !       pbar(i) =  fext_imp(i)*SIN(3.14159*(t-delay)/per1) + fext(i) - fint(i)
 !       !pbar(i) = fext(i) + (0.5 - 0.5*COS(3.14159/per1*(t-delay)))*fext_imp(i) - fint(i)
 !    ELSE
 !       pbar(i) =  fext(i) - fint(i)
 !       !pbar(i) = fext(i) + fext_imp(i) - fint(i)
 !    END IF
     pbar(i) = fext(i) + fext_imp(i) - fint(i)
  ENDDO
  
  print*,myid,'MAXIMUM FORCING TERM:',MAXVAL(pbar(:) - fext(:) + fint(:))

  ! Allocate shortened arrays
  ALLOCATE(newpbar(1:nrows_meff))
  ALLOCATE(newa(1:nrows_meff))
  newpbar(:) = 0.0
  newa(:) = 0.0

  ! Remove rows that are constrained by displacement BC's
  CALL removeBCs_pbar(nstart_km,3*LNumNp,pbar,nrows_meff,newpbar)
     
     
!
!     Solve for acceleration at time t + Dt
!     _        _
!     M a    = P
!        k+1    k+1
!

  ! Solve using BlockSolve95 functions
  CALL BS95SOLVE(nrows_meff,newpbar,newa,contol,BS95debug)

  ! Put newa back into a
  CALL removeBCs_newa(nstart_km,3*LNumNp,a,nrows_meff,newa)

  ! Deallocate shortened arrays
  DEALLOCATE(newpbar)
  DEALLOCATE(newa)


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
     d(i) = d(i) + alphaimp*delt*delt*a(i)
     v(i) = v(i) + deltaimp*delt*a(i)
  ENDDO


!     
!     4. Apply boundary conditions
!
  CALL implicit_bc_enforce(global%NumNp,LNumNp,d,v,a,node_flag,boundary_value,t,myid)


  IF(myid==0) PRINT*,'ASSEMBLING NEW DISPLACEMENT VECTOR'


!
! Put the vectors back into their Rocfrac equivalent vectors
!

  DO i = 1, global%NumNp
     IF ( NodeProc(i) == myid ) THEN
        DO j = 1, 3
           k = 3*Local2Global(i) - nstart_km + 1 + j - 3
           m = 3*i + j - 3
           global%Disp( m ) = d( k )
           global%VeloHalf( m ) = v( k )
           global%Accel( m ) = a( k )
        END DO
     END IF
  END DO


!
! Communicate variables to other processors and assemble
!

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

! receive from other procs
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(frmproc(1:NumCommProcsFrom2))
     DO i = 1, NumCommProcsFrom2
        ALLOCATE(frmproc(i)%rcvbuf(1:10*NumCommNodesFrom2(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),10*NumCommNodesFrom2(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom2(i),10,ROCSTAR_COMMUNICATOR, &
             req_rcv(i),ierr)
     ENDDO

! send to other procs
     DO i = 1, NumCommProcs2
        ALLOCATE(bufsnd(1:10*NumCommNodes2(i)))
        counter = 0
        DO j = 1, NumCommNodes2(i)
           counter = counter + 1
           bufsnd(counter) = 1.0 * CommNodes2(i,j)
           k = Global2Local(CommNodes2(i,j))
           DO p = 1, 3
              counter = counter + 1
              bufsnd(counter) = global%Disp( k * 3 - 3 + p )
           ENDDO
           DO p = 1, 3
              counter = counter + 1
              bufsnd(counter) = global%VeloHalf( k * 3 - 3 + p )
           ENDDO
           DO p = 1, 3
              counter = counter + 1
              bufsnd(counter) = global%Accel( k * 3 - 3 + p )
           ENDDO
        ENDDO
        CALL MPI_SEND(bufsnd(:),10*NumCommNodes2(i),MPI_DOUBLE_PRECISION, &
             CommProcs2(i),10,ROCSTAR_COMMUNICATOR,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)

! deallocate status and request vars
     DEALLOCATE(req_rcv)
     DEALLOCATE(req_snd)
     DEALLOCATE(stat_snd)
     DEALLOCATE(stat_rcv)

! assemble into Rocfrac nodal variable vectors
     DO i = 1, NumCommProcsFrom2
        counter = 0
        DO j = 1, NumCommNodesFrom2(i)
           counter = counter + 1
           k = Global2Local( INT( frmproc(i)%rcvbuf(counter) ) )  ! local index of node
           DO p = 1, 3
              counter = counter + 1
              global%Disp( 3*k - 3 + p ) = frmproc(i)%rcvbuf(counter)
           END DO
           DO p = 1, 3
              counter = counter + 1
              global%VeloHalf( 3*k - 3 + p ) = frmproc(i)%rcvbuf(counter)
           END DO
           DO p = 1, 3
              counter = counter + 1
              global%Accel( 3*k - 3 + p ) = frmproc(i)%rcvbuf(counter)
           END DO
        END DO
     END DO

! clear the receive buffer
     DEALLOCATE(frmproc)

  ENDIF


END SUBROUTINE implicit_soln

