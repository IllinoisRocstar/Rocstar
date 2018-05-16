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
SUBROUTINE NewCommList(NumEl, NumNP, NumProcs, NumVertex, &
     NodeConn,ElPart,NumProcPerNd,ProcNdList,MaxNumberOfProcsToShareNode,&
     NumElPerProc,NumNdPerProc)

  USE CommGlobal
  USE Linked_List2


  IMPLICIT NONE

  INTEGER :: i, j, k, m, npart,icnt

  INTEGER :: NumEl
  INTEGER :: NumNP 
  INTEGER :: NumProcs
  INTEGER :: NumVertex
  integer :: MaxNumberOfProcsToShareNode

  INTEGER, DIMENSION(1:NumVertex,1:NumEl) :: NodeConn
  INTEGER, DIMENSION(1:NumEl) :: ElPart
  INTEGER, DIMENSION(1:NumNP) ::  NumProcPerNd ! NumProcPerNd
  
  INTEGER, DIMENSION(1:NumNP,1:MaxNumberOfProcsToShareNode) :: ProcNdList
  INTEGER, DIMENSION(1:NumProcs) :: NumElperProc, NumNdperProc
  INTEGER :: LocNodeNum, ElID
  LOGICAL :: debug = .TRUE.
  TYPE proclist
     INTEGER, DIMENSION(:), POINTER :: proc_list
  END TYPE proclist
  TYPE(proclist), ALLOCATABLE, DIMENSION(:) :: node


  INTEGER, ALLOCATABLE, DIMENSION(:) :: nproc_neigh_lst

  TYPE(ProcElemList_data_ptr), pointer :: ptr

  ALLOCATE(ID_sendto(1:NumProcs,1:NumProcs))

  ID_sendto(1:NumProcs,1:NumProcs)%num_border_comm = 0


!         Here's a quick way to build communication lists in any
! FEM-type program.  There are  e  elements,  p  processors,
! and  n  nodes.  When I write that a table maps foo->bar,
! I mean table[foo] .eq. bar

! 1.) Partition the elements (using, e.g. metis).  This will
!     give you a table mapping element->processor.

! 2.) Build a list of the elements on each processor, by going
! through the element->processor table once, adding each element
! to the processor's local element list.  Time is O(e).

! Initialize list

  allocate(ProcElemList(1:NumProcs))
  DO i = 1, NumProcs
     NULLIFY(ProcElemList(i)%head,ProcElemList(i)%tail)
  enddo

! Build up list (add to stack)
  DO i = 1, NumEl
     npart = ElPart(i)
     allocate(ProcElem_Item)
     ProcElem_Item%GlbElNum = i
     CALL add_ProcElemList(ProcElem_Item,ProcElemList(npart)%head,ProcElemList(npart)%tail)
  ENDDO


! 3.) Build a list of the nodes on each processor.  The way
! you do this is critical:
!         -First, build an emtpy node->processors table, to
! store a list of processors that touch each node.  Time is O(n).
!         -For each processor, for each local element,
! add this processor to all the element's nodes'
! processor lists (if it's not already in the list).
! The nodes' processor lists will be very short
! (almost always 0 or 1), so this step will be O(n).
!         -Go through the node->processors table, and
! add each node to all its adjacent processors.
! Nearly all node->processors lists have length 1, so
! this step is also O(n).

  NumProcPerNd(:) = 0

  DO i = 1, NumProcs ! Number of processors
     NumElPerProc(i) = 0 ! Get_Len_ProcElemList(ProcElemList_head(i)) ! length of linked list, i.e. number of elements per processor

     ptr => ProcElemList(i)%head
 
     DO WHILE(ASSOCIATED(ptr))

        NumElPerProc(i) =  NumElPerProc(i) + 1

        ElID = ptr%GlbElNum

        DO k = 1, NumVertex
           LocNodeNum = NodeConn(k,ElID)
           DO m = 1, NumProcPerNd(LocNodeNum)
              IF(ProcNdList(LocNodeNum,m).EQ.i)THEN
                 GOTO 2
              ENDIF
           ENDDO
           NumProcPerNd(LocNodeNum) = NumProcPerNd(LocNodeNum) + 1
           IF(NumProcPerNd(LocNodeNum).GT.MaxNumberOfProcsToShareNode)THEN
              PRINT*,'Error in NewCommlist'
              PRINT*,'Greater then ',MaxNumberOfProcsToShareNode, &
                   ' processors share a node on the boundary'
              PRINT*,'Increase MaxNumberOfProcsToShareNode'
              PRINT*,'MaxNumberOfProcsToShareNode =', MaxNumberOfProcsToShareNode
              PRINT*,' Needed Size= ', NumProcPerNd(LocNodeNum)
              PRINT*,'Stopping'
              STOP
           endif
           ProcNdList(LocNodeNum,NumProcPerNd(LocNodeNum)) = i
           NumNdPerProc(i) = NumNdPerProc(i) + 1
2          CONTINUE
        ENDDO
        
        ptr => ptr%next

     ENDDO
  ENDDO

! 4.) Build the communication lists.  We can re-use the node->processors
! table we built before-- go through the table, and for each node:
!         -If the processor list for this node has length 1 (usual case),
! the node is local to the processor and never needs to be sent anywhere.
!         -If the length is greater than 1, the node is shared
! and you add it to the comm. lists for each listed processor.
! Time for this step is O(n).
  DO i = 1, NumProcs
     DO j = 1, NumProcs
        NULLIFY(ID_sendto(i,j)%comm_head )
        NULLIFY(ID_sendto(i,j)%comm_tail )
     enddo
  enddo

  DO i = 1, NumNp
     IF(NumProcPerNd(i).NE.1)THEN
        DO j = 1, NumProcPerNd(i)
           DO k = 1, NumProcPerNd(i)
              IF(k.NE.j) CALL addcommNd(ID_sendto(ProcNdList(i,j),ProcNdList(i,k)),i)
           ENDDO
        ENDDO
     ENDIF
  ENDDO

!!$  ALLOCATE(nproc_neigh_lst(1:NumProcs))
!!$  nproc_neigh_lst(1:NumProcs) = 0
!!$  DO i = 1, NumProcs
!!$!
!!$!     Determine the neighbor of processors 'i' is communicating with.
!!$!     
!!$     DO j = 1,NumProcs
!!$        IF(ID_sendto(i,j)%num_border_comm.NE.0) &
!!$             nproc_neigh_lst(i) = nproc_neigh_lst(i) + 1
!!$     ENDDO
!!$
!!$     WRITE(4000+i-1,*) nproc_neigh_lst(i)
!!$
!!$
!!$     
!!$     DO j = 1, NumProcs            ! receiving processor
!!$        IF(ID_sendto(j,i)%num_border_comm.NE.0)THEN
!!$!     Number of nodes that need to be communicated for R_in calculation
!!$           WRITE(4000+i-1,*) j-1,ID_sendto(j,i)%num_border_comm ! common
!!$!     List of nodes that need to be communicated for R_in calculation
!!$           CALL print_comm_list(ID_sendto(j,i))
!!$        ENDIF
!!$     ENDDO
!!$
!!$  ENDDO

! That's it-- total time is O(n+e), as long as most nodes and
! elements are not shared.  Note that the node->processors
! table we build in 3 is critical to the performance of both
! 3 and 4.  You could use a simpler implementation for step 3,
! such as:
!         for each processor; for each element; add the
! nodes around the element to the processor's node list
! (if it's not already there).
! 

! This alternate step 3 actually runs slower, because
! you have to check the processor's long node list for
! duplicates [it's O(p*(n/p)^2)].  Worse yet, 4 becomes
! horrifically expensive without a node->processors table--
! because you have to check if every node is shared,
! by running through all other processor's nodes, it's O(n^2)!
!    
! For  n = 10 million nodes, n^2 is unbelievably painful--
! 100 trillion.  The node->processors table, by contrast, can easily be
! done (in C) with about 12 bytes of memory per node (a processor #,
! a processor-local node #, and a rarely-used next pointer), almost
! all allocated contiguously.
! The node->processors table is the way to go for neighbor lists.




END SUBROUTINE NewCommList

