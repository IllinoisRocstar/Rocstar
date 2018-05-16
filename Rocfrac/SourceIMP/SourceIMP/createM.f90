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
SUBROUTINE createM(global)

  USE implicit_global
  USE comp_row_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL) :: global

  INTEGER :: i, j, m, p, inode, jnode, counter, idof, jdof

  REAL(kind=wp) :: tempKval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: AmountToReceive, AmountToSend
  
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd
  
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: tempmg, mg
  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_m_temp
  INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_m_temp, cval_m_temp

  REAL(kind=wp), DIMENSION(1:3,1:8) :: ri = RESHAPE( &
       (/-0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626, 0.577350269189626, 0.577350269189626, &
         -0.577350269189626, 0.577350269189626, 0.577350269189626/),(/3,8/) )



  ! Construct the global mass matrix
  IF(myid==0) PRINT*,'CONSTRUCTING THE MASS MATRIX'

  ! Construct the mass matrix as a consistent matrix
  CALL implicit_v3d8_mass_consistent(global%NumElVol,global%NumNp,global%NumMatVol,global%MeshCoor,global%ElConnVol,global%MatIdVol,ri,global%rho,global%ElConnVol)
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
     ALLOCATE(AmountToReceive(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        CALL MPI_IRECV(AmountToReceive(i),1, &
             MPI_INTEGER,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
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
       !      CommProcs1(i),10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
        CALL MPI_SEND(AmountToSend(i),1,MPI_INTEGER, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
     ENDDO
     !CALL MPI_WAITALL(NumCommProcsFrom1,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(NumCommProcs1,req_snd,stat_snd,ierr)
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

     ! Communicate parts of the M matrix to other procs
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(frmproc(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(frmproc(i)%rcvbuf(1:3*AmountToReceive(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),3*AmountToReceive(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
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
        !     CommProcs1(i),10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
        CALL MPI_SEND(bufsnd,3*AmountToSend(i),MPI_DOUBLE_PRECISION, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
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
  ALLOCATE(tempmg(1:global%NumNp))
  ALLOCATE(mg(1:3*LNumNp))
  ALLOCATE(cval_m_temp(1:3*LNumNp))
  ALLOCATE(aval_m_temp(1:3*LNumNp))
  ALLOCATE(rp_m_temp(1:3*LNumNp+1))
  tempmg(1:global%NumNp) = 0.0
  CALL implicit_v3d8_mass(global%NumElVol,global%NumNp,global%NumMatVol,global%MeshCoor,global%ElConnVol,global%MatIdVol,ri,global%rho,tempmg,1,global%NumElVol)
  DO i = 1, global%NumNp
     tempmg(i) = 1.0 / tempmg(i)  ! Get the mass matrix, not its inverse
  ENDDO

  IF (nprocs > 1) THEN
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
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
             MPI_DOUBLE_PRECISION,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
             req_rcv(i),ierr)
     ENDDO
     DO i = 1, NumCommProcs1
        ALLOCATE(bufsnd(2*NumCommNodes1(i)))
        DO j = 1, NumCommNodes1(i)
           bufsnd(2*j-1) = Local2Global(CommNodes1(i,j))
           bufsnd(2*j) = tempmg(CommNodes1(i,j))
        ENDDO
        !CALL MPI_ISEND(bufsnd,2*NumCommNodes1(i),MPI_DOUBLE_PRECISION, &
        !     CommProcs1(i),10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
        CALL MPI_SEND(bufsnd,2*NumCommNodes1(i),MPI_DOUBLE_PRECISION, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
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
     DO i = 1, global%NumNp
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

END SUBROUTINE createM

