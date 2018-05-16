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
!!     InitComm1
!!
!!  FUNCTION
!!     Initializes global communication variables for communicating
!!     boundary data from this proc to other procs
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE InitComm1(global)

  USE Precision
  USE implicit_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL)  :: global
  INTEGER :: i, j, counter  ! counters
  INTEGER :: tempnodes
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumCommNodesFromAll  ! Number of nodes to be recieved from all processors



!
! Determine communications to be sent to other procs
!


  ! Count the processors to be communicated with
  NumCommProcs1 = 0
  DO i = 0, nprocs
     IF (i /= myid) THEN
        counter = 0
        DO j = 1, global%NumNp
           IF (NodeProc(j) == i) THEN
              counter = 1
              EXIT
           ENDIF
        ENDDO
        IF (counter == 1) THEN  ! New proc found
           NumCommProcs1 = NumCommProcs1 + 1
        ENDIF
     ENDIF
  ENDDO
  
  ! Number the processors to be communicated with
  ALLOCATE(CommProcs1(1:NumCommProcs1))
  counter = 0
  DO i = 0, nprocs
     IF (i /= myid) THEN
        DO j = 1, global%NumNp
           IF (NodeProc(j) == i) THEN
              counter = counter + 1
              CommProcs1(counter) = i
              EXIT
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  
  ! Count the nodes to be sent to each processor
  ALLOCATE(NumCommNodes1(1:NumCommProcs1))
  DO i = 1, NumCommProcs1
     NumCommNodes1(i) = 0
     DO j = 1, global%NumNp
        IF (NodeProc(j) == CommProcs1(i)) THEN
           NumCommNodes1(i) = NumCommNodes1(i) + 1
        ENDIF
     ENDDO
  ENDDO
  
  ! Number the nodes to be sent to each processor
  MaxNumCommNodes1 = MAXVAL(NumCommNodes1)
  ALLOCATE(CommNodes1(1:NumCommProcs1,1:MaxNumCommNodes1))
  DO i = 1, NumCommProcs1
     CommNodes1(i,1:MaxNumCommNodes1) = -1  ! Flag everything as unused
     counter = 0
     DO j = 1, global%NumNp
        IF (NodeProc(j) == CommProcs1(i)) THEN
           counter = counter + 1
           CommNodes1(i,counter) = j
        ENDIF
     ENDDO
  ENDDO


!
! Communicate how many nodes will be sent to each processor
!
  
  ! Allocate MPI arrays
  ALLOCATE(frmproc(0:nprocs-1))
  ALLOCATE(req_rcv(1:nprocs))
  ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:nprocs))
  ALLOCATE(req_snd(1:nprocs))
  ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:nprocs))
  DO i = 1, nprocs
     req_rcv(i) = 0
     req_snd(i) = 0
     DO j = 1, MPI_STATUS_SIZE
        stat_rcv(1,i) = 0
        stat_snd(j,i) = 0
     ENDDO
  ENDDO

  ! Perform the communication
  ALLOCATE(NumCommNodesFromAll(0:nprocs-1))
  DO i = 0, nprocs-1
     CALL MPI_IRECV(NumCommNodesFromAll(i),1, &
          MPI_INTEGER,i,10,ROCSTAR_COMMUNICATOR, &
          req_rcv(i+1),ierr)
  ENDDO
  DO i = 0, nprocs-1
     tempnodes = 0
     DO j = 1, NumCommProcs1
        IF (CommProcs1(j) == i) THEN
           tempnodes = NumCommNodes1(j)
        ENDIF
     ENDDO
     CALL MPI_ISEND(tempnodes,1,MPI_INTEGER, &
     i,10,ROCSTAR_COMMUNICATOR,req_snd(i+1),ierr)
  ENDDO
  CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
  CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

  ! Deallocate MPI arrays
  DEALLOCATE(frmproc)
  DEALLOCATE(req_rcv)
  DEALLOCATE(stat_rcv)
  DEALLOCATE(req_snd)
  DEALLOCATE(stat_snd)

  ! Count the number of procs that will send me data
  NumCommProcsFrom1 = 0
  DO i = 0, nprocs-1
     IF (NumCommNodesFromAll(i) > 0) THEN
        NumCommProcsFrom1 = NumCommProcsFrom1 + 1
     ENDIF
  ENDDO

  ! Number the procs that will send me data
  ALLOCATE(CommProcsFrom1(1:NumCommProcsFrom1))
  ALLOCATE(NumCommNodesFrom1(1:NumCommProcsFrom1))
  counter = 0
  DO i = 0, nprocs-1
     IF (NumCommNodesFromAll(i) > 0) THEN
        counter = counter + 1
        CommProcsFrom1(counter) = i
        NumCommNodesFrom1(counter) = NumCommNodesFromAll(i)
     ENDIF
  ENDDO


  ! Deallocate stuff
  DEALLOCATE(NumCommNodesFromAll)



END SUBROUTINE InitComm1


















!!****
!!
!!  NAME
!!     InitComm2
!!
!!  FUNCTION
!!     Initializes global communication variables for communicating
!!     boundary data from other procs to this proc
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE InitComm2(global)

  USE Precision
  USE implicit_global

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL)  :: global
  INTEGER :: i, j, counter  ! counters
  INTEGER :: tempnodes
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumCommNodesFromAll  ! Number of nodes to be recieved from all processors
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CommNodesFrom2
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: bufsnd  ! MPI send buffer


!
! Determine communications to be sent to other procs
!


  ! Count the processors that will send me data
  NumCommProcsFrom2 = 0
  DO i = 0, nprocs
     IF (i /= myid) THEN
        counter = 0
        DO j = 1, global%NumNp
           IF (NodeProc(j) == i) THEN
              counter = 1
              EXIT
           ENDIF
        ENDDO
        IF (counter == 1) THEN  ! New proc found
           NumCommProcsFrom2 = NumCommProcsFrom2 + 1
        ENDIF
     ENDIF
  ENDDO
  
  ! Number the processors that will send me data
  ALLOCATE(CommProcsFrom2(1:NumCommProcsFrom2))
  counter = 0
  DO i = 0, nprocs
     IF (i /= myid) THEN
        DO j = 1, global%NumNp
           IF (NodeProc(j) == i) THEN
              counter = counter + 1
              CommProcsFrom2(counter) = i
              EXIT
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  
  ! Count the nodes to be received from each processor
  ALLOCATE(NumCommNodesFrom2(1:NumCommProcsFrom2))
  DO i = 1, NumCommProcsFrom2
     NumCommNodesFrom2(i) = 0
     DO j = 1, global%NumNp
        IF (NodeProc(j) == CommProcsFrom2(i)) THEN
           NumCommNodesFrom2(i) = NumCommNodesFrom2(i) + 1
        ENDIF
     ENDDO
  ENDDO


!
! Communicate how many nodes are needed to each processor
!
  
  ! Allocate MPI arrays
  ALLOCATE(frmproc(0:nprocs-1))
  ALLOCATE(req_rcv(1:nprocs))
  ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:nprocs))
  ALLOCATE(req_snd(1:nprocs))
  ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:nprocs))
  DO i = 1, nprocs
     req_rcv(i) = 0
     req_snd(i) = 0
     DO j = 1, MPI_STATUS_SIZE
        stat_rcv(1,i) = 0
        stat_snd(j,i) = 0
     ENDDO
  ENDDO

  ! Perform the communication
  ALLOCATE(NumCommNodesFromAll(0:nprocs-1))
  DO i = 0, nprocs-1
     CALL MPI_IRECV(NumCommNodesFromAll(i),1, &
          MPI_INTEGER,i,10,ROCSTAR_COMMUNICATOR, &
          req_rcv(i+1),ierr)
  ENDDO
  DO i = 0, nprocs-1
     tempnodes = 0
     DO j = 1, NumCommProcsFrom2
        IF (CommProcsFrom2(j) == i) THEN
           tempnodes = NumCommNodesFrom2(j)
        ENDIF
     ENDDO
     CALL MPI_ISEND(tempnodes,1,MPI_INTEGER, &
     i,10,ROCSTAR_COMMUNICATOR,req_snd(i+1),ierr)
  ENDDO
  CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
  CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

  ! Deallocate MPI arrays
  DEALLOCATE(frmproc)
  DEALLOCATE(req_rcv)
  DEALLOCATE(stat_rcv)
  DEALLOCATE(req_snd)
  DEALLOCATE(stat_snd)

  ! Count the number of procs that I will send data to
  NumCommProcs2 = 0
  DO i = 0, nprocs-1
     IF (NumCommNodesFromAll(i) > 0) THEN
        NumCommProcs2 = NumCommProcs2 + 1
     ENDIF
  ENDDO

  ! Number the procs that I will send data to
  ALLOCATE(CommProcs2(1:NumCommProcs2))
  ALLOCATE(NumCommNodes2(1:NumCommProcs2))
  counter = 0
  DO i = 0, nprocs-1
     IF (NumCommNodesFromAll(i) > 0) THEN
        counter = counter + 1
        CommProcs2(counter) = i
        NumCommNodes2(counter) = NumCommNodesFromAll(i)
     ENDIF
  ENDDO


  ! Deallocate stuff
  DEALLOCATE(NumCommNodesFromAll)

  ! Number the nodes to be received from each processor
  MaxNumCommNodes2 = MAXVAL(NumCommNodesFrom2)
  ALLOCATE(CommNodesFrom2(1:NumCommProcsFrom2,1:MaxNumCommNodes2))
  DO i = 1, NumCommProcsFrom2
     CommNodesFrom2(i,1:MaxNumCommNodes2) = -1  ! Flag everything as unused
     counter = 0
     DO j = 1, global%NumNp
        IF (NodeProc(j) == CommProcsFrom2(i)) THEN
           counter = counter + 1
           CommNodesFrom2(i,counter) = j
        ENDIF
     ENDDO
  ENDDO

  !DO i = 1, NumCommProcs
  !   print*,myid,' has to receive info about ',NumCommNodesFrom(i),' nodes from proc ',CommProcsFrom(i)
  !   print*,myid,' CommNodesFrom(i,:) = ',Local2Global(CommNodesFrom(i,:))
  !ENDDO


!
! Communicate which nodes are needed to each processor
!
  
  ! Allocate MPI arrays
  !ALLOCATE(req_rcv(1:NumCommProcs2))
  !ALLOCATE(req_snd(1:NumCommProcsFrom2))
  !ALLOCATE(stat_snd(1:MPI_STATUS_SIZE,1:NumCommProcsFrom2))
  !ALLOCATE(stat_rcv(1:MPI_STATUS_SIZE,1:NumCommProcs2))
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

  ! Perform the communication
  ALLOCATE(frmproc(1:NumCommProcs2))
  DO i = 1, NumCommProcs2
     ALLOCATE(frmproc(i)%rcvbuf(1:NumCommNodes2(i)))
     CALL MPI_IRECV(frmproc(i)%rcvbuf(1),NumCommNodes2(i), &
          MPI_DOUBLE_PRECISION,CommProcs2(i),10,ROCSTAR_COMMUNICATOR, &
          req_rcv(i+1),ierr)
  ENDDO
  DO i = 1, NumCommProcsFrom2
     ALLOCATE(bufsnd(1:NumCommNodesFrom2(i)))
     DO j = 1, NumCommNodesFrom2(i)
        bufsnd(j) = Local2Global(CommNodesFrom2(i,j))
     ENDDO
     !print*,myid,' bufsnd = ',bufsnd(:)
     CALL MPI_ISEND(bufsnd,NumCommNodesFrom2(i),MPI_DOUBLE_PRECISION, &
     CommProcsFrom2(i),10,ROCSTAR_COMMUNICATOR,req_snd(i+1),ierr)
     DEALLOCATE(bufsnd)
  ENDDO
  !CALL MPI_WAITALL(NumCommProcs2,req_rcv,stat_rcv,ierr)
  !CALL MPI_WAITALL(NumCommProcsFrom2,req_snd,stat_snd,ierr)
  CALL MPI_WAITALL(nprocs,req_rcv,stat_rcv,ierr)
  CALL MPI_WAITALL(nprocs,req_snd,stat_snd,ierr)

  ! Deallocate MPI arrays
  DEALLOCATE(req_rcv)
  DEALLOCATE(stat_rcv)
  DEALLOCATE(req_snd)
  DEALLOCATE(stat_snd)

  !DO i = 1, NumCommProcs
  !   DO j = 1, NumCommNodes(i)
  !      print*,myid,' has to send info about global node ',frmproc(i)%rcvbuf(j),' to proc ',CommProcs(i)
  !   ENDDO
  !ENDDO

  ! Number the nodes to be sent to each processor
  MaxNumCommNodes2 = MAXVAL(NumCommNodes2)
  ALLOCATE(CommNodes2(1:NumCommProcs2,1:MaxNumCommNodes2))
  DO i = 1, NumCommProcs2
     CommNodes2(i,1:MaxNumCommNodes2) = -1  ! Flag everything as unused
     DO j = 1, NumCommNodes2(i)
        CommNodes2(i,j) = INT(frmproc(i)%rcvbuf(j))
     ENDDO
  ENDDO

  ! Deallocate stuff
  DEALLOCATE(frmproc)
  DEALLOCATE(CommNodesFrom2)

  

END SUBROUTINE InitComm2

