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
SUBROUTINE createK(global)

  USE implicit_global
  USE comp_row_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE int_buf
     INTEGER, DIMENSION(:), POINTER :: rcvbuf
  END TYPE int_buf

  TYPE(ROCFRAC_GLOBAL) :: global
  INTEGER :: i, j, m, p, idof, jdof, inode, jnode, counter
  REAL(kind=wp) :: tempKval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: AmountToReceive, AmountToSend
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd
  TYPE(int_buf), ALLOCATABLE, DIMENSION(:) :: Ki, Kj  


  ! Construct the global stiffness matrix
  IF(myid==0) PRINT*,'CONSTRUCTING THE STIFFNESS MATRIX'
!  CALL implicit_v3d8_me_K(global%MeshCoor, global%ElConnVol, global%ci, global%NumNp, 1, global%NumElVol, &
!       global%NumElVol, global%MatIdVol, global%NumMatVol, global%enhanced_map, global%mixed_map, &
!       global%NumMatVol)
  CALL implicit_v3d8_me_K(global%MeshCoor, global%ElConnVol, global%ci_full, global%NumNp, 1, global%NumElVol, &
       global%NumElVol, global%MatIdVol, global%NumMatVol, global%enhanced_map, global%mixed_map, &
       global%NumMatVol)
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
        !     CommProcs1(i),10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
        CALL MPI_SEND(AmountToSend(i),1,MPI_INTEGER, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
     ENDDO
     !CALL MPI_WAITALL(NumCommProcsFrom1,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(NumCommProcs1,req_snd,stat_snd,ierr)
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)


     ! Communicate i locations of nonzeros in K matrix to other procs
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(Ki(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(Ki(i)%rcvbuf(1:AmountToReceive(i)))
        CALL MPI_IRECV(Ki(i)%rcvbuf(1),AmountToReceive(i), &
             MPI_INTEGER,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
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
        !     CommProcs1(i),10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
        CALL MPI_SEND(INT(bufsnd(:)),AmountToSend(i),MPI_INTEGER, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

     ! Communicate j locations of nonzeros in K matrix to other procs
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(Kj(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(Kj(i)%rcvbuf(1:AmountToReceive(i)))
        CALL MPI_IRECV(Kj(i)%rcvbuf(1),AmountToReceive(i), &
             MPI_INTEGER,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
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
        !     CommProcs1(i),10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
        CALL MPI_SEND(INT(bufsnd(:)),AmountToSend(i),MPI_INTEGER, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
        DEALLOCATE(bufsnd)
     ENDDO
     CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
     !CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

     ! Communicate values of nonzeros in K matrix to other procs
     CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
     ALLOCATE(frmproc(1:NumCommProcsFrom1))
     DO i = 1, NumCommProcsFrom1
        ALLOCATE(frmproc(i)%rcvbuf(1:AmountToReceive(i)))
        CALL MPI_IRECV(frmproc(i)%rcvbuf(1),AmountToReceive(i), &
             MPI_DOUBLE_PRECISION,CommProcsFrom1(i),10,ROCSTAR_COMMUNICATOR, &
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
        !     CommProcs1(i),10,ROCSTAR_COMMUNICATOR,req_snd(i),ierr)
        CALL MPI_SEND(bufsnd,AmountToSend(i),MPI_DOUBLE_PRECISION, &
             CommProcs1(i),10,ROCSTAR_COMMUNICATOR,ierr)
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

END SUBROUTINE createK


! LocalWords:  RocFracComm
