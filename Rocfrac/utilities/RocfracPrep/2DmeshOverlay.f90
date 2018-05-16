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
SUBROUTINE mesh2dOverlay(nprocs,iProcs,ichr4)

  use meshdata
  USE Linked_List2

  IMPLICIT NONE

  INCLUDE 'comf90.h'

  CHARACTER*4 :: ichr4
  INTEGER :: NEW
  INTEGER :: i,k
  INTEGER :: icounter
  INTEGER :: iProcs
  INTEGER :: nprocs
  INTEGER :: NumNpNew
  INTEGER :: glbNdNum
  LOGICAL :: SurfaceElOnProc
  INTEGER :: TestNumSurfEl
  INTEGER, POINTER, DIMENSION(:,:) :: NodeFlag2D
  INTEGER :: iTestCnt

! User-defined list element
! The Link_Type field MUST be the FIRST in the user-defined list element 
! Note pointer to data so as to easily create sublists

  TYPE User_Type_SurfNdList
     TYPE(Link_Type) :: Link
     TYPE(User_Data_Type_ProcNodeList), POINTER :: Data  
  END TYPE User_Type_SurfNdList

  TYPE User_Data_Type_SurfNdList
     INTEGER :: GlbNdNum
  END TYPE User_Data_Type_SurfNdList

! Auxilliary data type required for the transfer function
  TYPE User_Ptr_Type_SurfNdList
     TYPE(User_Type_SurfNdList), POINTER :: P
  END TYPE User_Ptr_Type_SurfNdList

  TYPE(User_Ptr_Type_SurfNdList)  :: User_SurfNdList
  TYPE(List_Type) :: SurfNdList

  TYPE(Link_Ptr_Type)  :: Link

  TYPE(SurfMesh_tri3_ptr),POINTER :: ptrtri3
  TYPE(SurfMesh_tri6_ptr), POINTER :: ptrtri6
  TYPE(SurfMesh_hex8_ptr), POINTER :: ptrhex8

  INTEGER, POINTER, DIMENSION(:) :: BCSurfFlagZero,BCSurfFlagOne,BCSurfFlagTwo


  CHARACTER(*), PARAMETER :: surWin = "sfrac"

  REAL*8, POINTER, DIMENSION(:,:) :: MeshCoor

  INTEGER, POINTER, DIMENSION(:) :: ElFlag_List
  INTEGER, POINTER, DIMENSION(:) :: FaceOnCell

  INTEGER :: write_attr, set_option,  ierrFlg, sur_all

  INTEGER :: iNB, iBU, iNI

!  INTEGER :: NodeFlagOrg(1:numnp_prmry)

  iNI = iProcs*100 + 3

! obtain function handle ------------------------------------------------------

  write_attr = COM_get_function_handle( 'OUT.write_dataitem')
  set_option = COM_get_function_handle( 'OUT.set_option')



  PRINT*,'Overlay processor', iProcs - 1
  WRITE(456,*) iProcs - 1

!  NodeFlagOrg(1:numnp_prmry) = NodeFlag(1:numnp_prmry)


! --------------------------------------------------
! ------------------------------------------------
! (1)  Non Fluid/Structure Overlay
! ---------------------------------------
! -------------------------------------

! renumber elements nodes

  CALL COM_new_window( surWin)

!!$  IF(MeshType2D.EQ.6)THEN
!!$
!!$     NumNpNew = 0
!!$     NULLIFY(ptrtri6)
!!$! renumber elements nodes
!!$     ptrtri6 => SurfMesh_tri6_S_head
!!$
!!$     DO WHILE(ASSOCIATED(ptrtri6))
!!$! mark that changed surface mesh by changing it to a negative
!!$        DO i = 1, 6
!!$           glbNdNum = ptrtri6%ElemData(i)
!!$           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
!!$              NumNpNew = NumNpNew + 1
!!$              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
!!$           ENDIF
!!$        ENDDO
!!$        ptrtri6 => ptrtri6%next
!!$     ENDDO
!!$
!!$
!!$! check for special situation where there is just an edge on the surface
!!$
!!$     ptrtri6 => SurfMesh_tri6_S_head
!!$     do while(associated(ptrtri6))
!!$        DO i = 1, 6
!!$           glbNdNum = ptrtri6%ElemData(i)
!!$           DO k = 1,MaxNumberOfProcsToShareNode
!!$              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
!!$                 PRINT*,'Dectected surface sliver for pressure'
!!$                 NumNpNew = NumNpNew + 1
!!$                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
!!$              ENDIF
!!$           ENDDO
!!$        ENDDO
!!$        ptrtri6 => ptrtri6%next
!!$     ENDDO
!!$
!!$!     WRITE(4002,'(2i10)') NumNpNew,  NumEltet2D(2,iProcs)
!!$!     write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumEltet2D(2,iProcs),'F=FEPOINT, ET=TRIANGLE'
!!$
!!$     NodeFlag(:)= ABS( NodeFlag(:)) ! reset
!!$
!!$     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting
!!$
!!$     NumNpNew = 0
!!$! renumber elements nodes
!!$
!!$     ptrtri6 => SurfMesh_tri6_S_head
!!$
!!$     do while(associated(ptrtri6))
!!$! mark that changed surface mesh by changing it to a negative
!!$        DO i = 1, 6
!!$           glbNdNum = ptrtri6%ElemData(i)
!!$           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
!!$              NumNpNew = NumNpNew + 1
!!$!              WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!!$!              write(4003,*) coor(1:3,glbNdNum)
!!$              NodeFlag2D(1,NumNpNew) = glbNdNum
!!$              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
!!$              NodeFlag(ptrtri6%ElemData(i)) = - NumNpNew
!!$           endif
!!$        enddo
!!$        ptrtri6 => ptrtri6%next
!!$     ENDDO
!!$     ptrtri6 => SurfMesh_tri6_S_head
!!$     do while(associated(ptrtri6))
!!$        DO i = 1, 6
!!$           glbNdNum = ptrtri6%ElemData(i)
!!$           DO k = 1,MaxNumberOfProcsToShareNode
!!$              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
!!$                 PRINT*,'Dectected surface sliver for pressure'
!!$                 NumNpNew = NumNpNew + 1
!!$!                 WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!!$                 NodeFlag2D(1,NumNpNew) = glbNdNum
!!$                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
!!$                 NodeFlag(glbNdNum) = - NumNpNew
!!$              ENDIF
!!$           ENDDO
!!$        ENDDO
!!$        ptrtri6 => ptrtri6%next
!!$     ENDDO
!!$
!!$     
!!$     CALL COM_new_dataitem( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
!!$     CALL COM_set_size( surWin//'.nc', iNI, NumNpNew )
!!$     CALL COM_resize_array(surWin//'.nc', iNI, MeshCoor, 3)
!!$
!!$     DO i = 1, NumNpNew 
!!$        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
!!$     END DO
!!$
!!$
!!$     CALL COM_new_dataitem( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!!$!     CALL COM_set_size( surWin//'.bv', iNI, NumNpNew )
!!$     CALL COM_set_array( surWin//'.bv', iNI, NodeFlag2D(2,1), 2)
!!$
!!$
!!$! connectiveity
!!$     TestNumSurfEl = NumEltet2D(2,iProcs)
!!$     PRINT*,'Number of Non-Interacting nodes and elements',NumNpNew,TestNumSurfEl
!!$
!!$     CALL COM_set_size( surWin//'.:t6', iNI, TestNumSurfEl)
!!$     CALL COM_resize_array( surWin//'.:t6', iNI, ElConnTable, 6)
!!$
!!$
!!$     CALL COM_new_dataitem( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
!!$     CALL COM_resize_array(surWin//'.bf2c', iNI, ElFlag_List, 1)
!!$
!!$     CALL COM_new_dataitem( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
!!$     CALL COM_set_size( surWin//'.bcflag', iNI,  1)
!!$     CALL COM_resize_array(surWin//'.bcflag', iNI, BCSurfFlagTwo, 1)
!!$     
!!$     BCSurfFlagTwo(1) = 2
!!$
!!$     TestNumSurfEl = 0
!!$     ptrtri6 => SurfMesh_tri6_S_head
!!$     do while(associated(ptrtri6))
!!$        SurfaceElOnProc = .false.
!!$        IF(epart(ptrtri6%ElemData(7)).EQ.iProcs) SurfaceElOnProc = .TRUE.
!!$           TestNumSurfEl = TestNumSurfEl  + 1
!!$           ElConnTable(1:6,TestNumSurfEl) = ABS(NodeFlag(ptrtri6%ElemData(1:6)))
!!$        ENDIF
!!$        ptrtri6 => ptrtri6%next
!!$     ENDDO
!!$
!!$     IF(TestNumSurfEl.NE.NumEltet2D(2,iProcs))THEN
!!$        PRINT*,'ERROR, number of triangles from link list in Mesh2d'
!!$        PRINT*,'  different then that in read_patran'
!!$        PRINT*,TestNumSurfEl,NumEltet2D(2,iProcs)
!!$        STOP
!!$     ENDIF
!!$
!!$     NULLIFY(ptrtri6)     ! obtain function handle ------------------------------------------------------

!     write_attr = COM_get_function_handle( 'OUT.write_dataitem')
!     set_option = COM_get_function_handle( 'OUT.set_option')

!     IF(iProcs.EQ.1)THEN

!        CALL COM_call_function( set_option, 2, 'mode', 'w')

!     ELSE

!        CALL COM_call_function( set_option, 2, 'mode', 'a')
        
!     ENDIF

! do not append process rank -----------------

!     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

!     sur_all = Com_get_dataitem_handle( surWin//'.all')
!     CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_S.hdf', sur_all,&
!          "solid_surf","00.000000")

     


!!$  ELSE IF(MeshType2D.EQ.4)THEN
!!$
!!$
!!$     NumNpNew = 0
!!$     NULLIFY(ptrhex8)
!!$! renumber elements nodes
!!$
!!$     ptrhex8 => SurfMesh_hex8_S_head
!!$
!!$     DO WHILE(ASSOCIATED(ptrhex8))
!!$! mark that changed surface mesh by changing it to a negative
!!$        DO i = 1, 4
!!$           glbNdNum = ptrhex8%ElemData(i)
!!$           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
!!$              NumNpNew = NumNpNew + 1
!!$              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
!!$           ENDIF
!!$        ENDDO
!!$        ptrhex8 => ptrhex8%next
!!$     ENDDO
!!$
!!$
!!$     IF(NumElhex2D(2,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
!!$        PRINT*,'ERROR, found in surface mesh for Solid/Fluid interface'
!!$        PRINT*,'Found ',NumElhex2D(2,iProcs),'quads but no nodes'
!!$        PRINT*,'Stopping'
!!$        STOP
!!$     ENDIF
!!$     
!!$     IF(NumElhex2D(2,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
!!$        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid interface'
!!$        PRINT*,'Found ',NumNpNew,' nodes but no quads'
!!$     ENDIF
!!$
!!$! write the surface mesh files
!!$
!!$!     WRITE(4002,'(2i10)') NumNpNew,  NumElhex2D(2,iProcs)
!!$
!!$     NodeFlag(:)= ABS( NodeFlag(:)) ! reset
!!$
!!$     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting
!!$     NumNpNew = 0
!!$! renumber elements nodes
!!$
!!$     ptrhex8 => SurfMesh_hex8_S_head
!!$     do while(associated(ptrhex8))
!!$! mark that changed surface mesh by changing it to a negative
!!$        DO i = 1, 4
!!$           glbNdNum = ptrhex8%ElemData(i)
!!$           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
!!$              NumNpNew = NumNpNew + 1
!!$!              WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!!$              NodeFlag2D(1,NumNpNew) = glbNdNum
!!$              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
!!$              NodeFlag(glbNdNum) = - NumNpNew
!!$           endif
!!$        enddo
!!$        ptrhex8 => ptrhex8%next
!!$     ENDDO
!!$
!!$     CALL COM_new_dataitem( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
!!$     CALL COM_set_size( surWin//'.nc', iNI, NumNpNew )
!!$     CALL COM_resize_array(surWin//'.nc', iNI, MeshCoor, 3)
!!$
!!$     DO i = 1, NumNpNew 
!!$        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
!!$     END DO
!!$
!!$
!!$     CALL COM_new_dataitem( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!!$!     CALL COM_set_size( surWin//'.bv', iNI, NumNpNew )
!!$     CALL COM_set_array( surWin//'.bv', iNI, NodeFlag2D(2,1), 2)
!!$
!!$
!!$! connectivity
!!$     TestNumSurfEl = NumElhex2D(2,iProcs)
!!$     PRINT*,'Number of Non-Interacting nodes and elements',NumNpNew,TestNumSurfEl
!!$
!!$     CALL COM_set_size( surWin//'.:q4', iNI, TestNumSurfEl)
!!$     CALL COM_resize_array( surWin//'.:q4', iNI, ElConnTable, 4)
!!$
!!$
!!$     CALL COM_new_dataitem( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
!!$     CALL COM_resize_array(surWin//'.bf2c', iNI, ElFlag_List, 1)
!!$
!!$     CALL COM_new_dataitem( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
!!$     CALL COM_set_size( surWin//'.bcflag', iNI,  1)
!!$     CALL COM_resize_array(surWin//'.bcflag', iNI, BCSurfFlagTwo, 1)
!!$     
!!$
!!$     BCSurfFlagTwo(1) = 2
!!$! connectivity
!!$     NULLIFY(ptrhex8)
!!$     TestNumSurfEl = 0
!!$
!!$     ptrhex8 => SurfMesh_hex8_S_head
!!$     do while(associated(ptrhex8))
!!$        SurfaceElOnProc = .false.
!!$        IF(epart(ptrhex8%ElemData(5)).EQ.iProcs) SurfaceElOnProc = .TRUE.
!!$        IF(SurfaceElOnProc)THEN
!!$!           WRITE(4002,'(10i10)') ABS(NodeFlag(ptrhex8%ElemData(1:4))),&
!!$!             ElFlag(ptrhex8%ElemData(5)),1,1 ! last two are not really used
!!$           TestNumSurfEl = TestNumSurfEl + 1
!!$
!!$           ElConnTable(1:4,TestNumSurfEl) = ABS(NodeFlag(ptrhex8%ElemData(1:4)))
!!$        ENDIF
!!$        ptrhex8 => ptrhex8%next
!!$     ENDDO
!!$
!!$     IF(TestNumSurfEl.NE.NumElhex2D(2,iProcs))THEN
!!$        PRINT*,'Error: Number of Solid Surface hex elements in linked list'
!!$        PRINT*,' does not match that given by NumElhex2D(2,iProcs)'
!!$        PRINT*,' In linked list =',TestNumSurfEl
!!$        PRINT*,' NumElhex2D(2,iProcs) =',NumElhex2D(2,iProcs)
!!$        PRINT*,'Stopping'
!!$     ENDIF
!!$
!!$     NULLIFY(ptrhex8)
!!$
!!$     IF(associated(NodeFlag2D)) deallocate(NodeFlag2D)
!!$
!!$  ENDIF

  ! write surface window ------------------------
!!$  CALL COM_call_function( set_option, 2, 'mode', 'a')
!!$
!!$  sur_all = Com_get_dataitem_handle( SurWin//'.all')
!!$
!!$  CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh.'//ichr4, sur_all,&
!!$       "isolid","00.000000")

  IF(MeshType2D.EQ.3)THEN

     NumNpNew = 0
     NULLIFY(ptrtri3)
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_Ov1_head

     do while(associated(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO


! check for special situation where there is just an edge on the surface

     ptrtri3 => SurfMesh_tri3_Ov1_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO
!     WRITE(4002,'(2i10)') NumNpNew,  NumEltet2D(2,iProcs)
!     write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumEltet2D(2,iProcs),'F=FEPOINT, ET=TRIANGLE'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_Ov1_head

     do while(associated(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!              WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!              write(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(ptrtri3%ElemData(i)) = - NumNpNew
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO
     ptrtri3 => SurfMesh_tri3_Ov1_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
!                 WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
                 NodeFlag2D(1,NumNpNew) = glbNdNum
                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
                 NodeFlag(ptrtri3%ElemData(i)) = - NumNpNew
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO

     
     CALL COM_new_dataitem( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNI, NumNpNew )
     CALL COM_resize_array(surWin//'.nc', iNI, MeshCoor, 3)

     DO i = 1, NumNpNew 
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO

!!$     CALL COM_new_dataitem( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!!$     CALL COM_set_size( surWin//'.bv', iNI, NumNpNew )
!!$     CALL COM_set_array( surWin//'.bv', iNI, NodeFlag2D(2,1), 2)


! connectiveity
     TestNumSurfEl = NumEltet2D(4,iProcs) !0


     CALL COM_set_size( surWin//'.:t3', iNI, TestNumSurfEl)
     CALL COM_resize_array( surWin//'.:t3', iNI, ElConnTable, 3)

!!$     CALL COM_new_dataitem( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
!!$     CALL COM_resize_array(surWin//'.bf2c', iNI, ElFlag_List, 1)
!!$
!!$     CALL COM_new_dataitem( surWin//'.faceOnCell', 'e', COM_INTEGER, 1, '')
!!$     CALL COM_resize_array(surWin//'.faceOnCell', iNI, FaceOnCell, 1)

     WRITE(456,*) NumEltet2D(4,iProcs)

     TestNumSurfEl = 0
     ptrtri3 => SurfMesh_tri3_Ov1_head
     do while(associated(ptrtri3))
        SurfaceElOnProc = .false.
!!$        do i = 1, 3
!!$           IF(NodeFlag(ptrtri3%ElemData(i)).eq.0)THEN
!!$              SurfaceElOnProc = .false.
!!$              exit
!!$           endif
!!$        enddo
        IF(epart(ptrtri3%ElemData(4)).EQ.iProcs) SurfaceElOnProc = .true.
        IF(SurfaceElOnProc) THEN
 !          WRITE(4002,'(10i10)') ABS(NodeFlag(ptrtri3%ElemData(1:3))),&
 !            ElFlag(ptrtri3%ElemData(4)),2,1 ! last two are not really used
!!$        IF(SurfaceElOnProc) THEN
!!$           write(4003,*) ABS(NodeFlag(ptrtri3%ElemData(1:6)))
!!$        ENDIF
           TestNumSurfEl = TestNumSurfEl  + 1
           ElConnTable(1:3,TestNumSurfEl) = ABS(NodeFlag(ptrtri3%ElemData(1:3)))
           
           WRITE(456,*) ElFlag(ptrtri3%ElemData(4)),ptrtri3%ElemData(5)
!           ElFlag_List(TestNumSurfEl) = ElFlag(ptrtri3%ElemData(4))
!           FaceOnCell(TestNumSurfEl)  = ElFlag(ptrtri3%ElemData(5))
        ENDIF
        ptrtri3 => ptrtri3%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(4,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in Mesh2d'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumEltet2D(4,iProcs)
        STOP
     ENDIF



     IF(iProcs.EQ.1)THEN

        CALL COM_call_function( set_option, 2, 'mode', 'w')

     ELSE

        CALL COM_call_function( set_option, 2, 'mode', 'a')
        
     ENDIF


! do not append process rank -----------------

     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

     sur_all = Com_get_dataitem_handle( surWin//'.all')
     CALL COM_call_function( write_attr, 4, 'Rocin/A_SurfMesh_Ov.hdf', sur_all,&
          "solid_surf","00.000000")


     ! Correct the NodeFlag values that were over written

     NumNpNew = 0
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_Ov1_head

     DO WHILE(ASSOCIATED(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).LT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
              NodeFlag(ptrtri3%ElemData(i)) = NodeFlag2D(2,NumNpNew)
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO
     
! fix you need to but silever here if keeping that section

     IF(ASSOCIATED(NodeFlag2D)) DEALLOCATE(NodeFlag2D)
  ENDIF

  ! write surface window ------------------------
!!$  CALL COM_call_function( set_option, 2, 'mode', 'a')
!!$
!!$  sur_all = Com_get_dataitem_handle( SurWin//'.all')
!!$
!!$  CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_ov.'//ichr4, sur_all,&
!!$       "isolid","00.000000")

  CALL COM_delete_window( surWin)


  CALL COM_new_window( surWin)

!!$  DO i = 1, numnp_prmry
!!$     IF(NodeFlag(i).NE.NodeFlagOrg(i)) PRINT*,'lkjlk',NodeFlag(i), NodeFlagOrg(i)
!!$     stop
!!$  ENDDO

  IF(MeshType2D.EQ.3)THEN

     NumNpNew = 0
     NULLIFY(ptrtri3)
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_Ov2_head

     DO WHILE(ASSOCIATED(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           PRINT*,glbNdNum
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO

     PRINT*,'NumNpNew',NumNpNew


! check for special situation where there is just an edge on the surface

     ptrtri3 => SurfMesh_tri3_Ov2_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO

!     WRITE(4002,'(2i10)') NumNpNew,  NumEltet2D(2,iProcs)
!     write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumEltet2D(2,iProcs),'F=FEPOINT, ET=TRIANGLE'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset
     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_Ov2_head

     DO WHILE(ASSOCIATED(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!              WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!              write(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(ptrtri3%ElemData(i)) = - NumNpNew
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO
     ptrtri3 => SurfMesh_tri3_Ov2_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
!                 WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
                 NodeFlag2D(1,NumNpNew) = glbNdNum
                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
                 NodeFlag(ptrtri3%ElemData(i)) = - NumNpNew
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO


     
     CALL COM_new_dataitem( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNI, NumNpNew )
     CALL COM_resize_array(surWin//'.nc', iNI, MeshCoor, 3)

     DO i = 1, NumNpNew 
        PRINT*,coor(1:3,NodeFlag2D(1,i))
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO


! connectiveity
     TestNumSurfEl = NumEltet2D(5,iProcs) !0

     CALL COM_set_size( surWin//'.:t3', iNI, TestNumSurfEl)
     CALL COM_resize_array( surWin//'.:t3', iNI, ElConnTable, 3)

     WRITE(456,*) NumEltet2D(5,iProcs)

     TestNumSurfEl = 0
     ptrtri3 => SurfMesh_tri3_Ov2_head

     DO WHILE(ASSOCIATED(ptrtri3))
        SurfaceElOnProc = .FALSE.
        IF(epart(ptrtri3%ElemData(4)).EQ.iProcs) SurfaceElOnProc = .TRUE.
        IF(SurfaceElOnProc) THEN
           TestNumSurfEl = TestNumSurfEl  + 1
           ElConnTable(1:3,TestNumSurfEl) = ABS(NodeFlag(ptrtri3%ElemData(1:3)))
           WRITE(456,*) ElFlag(ptrtri3%ElemData(4)),ptrtri3%ElemData(5)
!!$         !  ElFlag_List(TestNumSurfEl) = ElFlag(ptrtri3%ElemData(4))
!!$         !  FaceOnCell(TestNumSurfEl)  = ElFlag(ptrtri3%ElemData(5))
        ENDIF
        ptrtri3 => ptrtri3%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(5,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in Mesh2d'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumEltet2D(5,iProcs)
        STOP
     ENDIF

     IF(iProcs.EQ.1)THEN

        CALL COM_call_function( set_option, 2, 'mode', 'w')

     ELSE

        CALL COM_call_function( set_option, 2, 'mode', 'a')
        
     ENDIF
     

! do not append process rank -----------------

     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

     sur_all = Com_get_dataitem_handle( surWin//'.all')


     CALL COM_call_function( write_attr, 4, 'Rocin/B_SurfMesh_Ov.hdf', sur_all,&
          "solid_surf","00.000000")


  ENDIF


  CALL COM_delete_window( surWin)


  IF(ASSOCIATED(NodeFlag2D)) DEALLOCATE(NodeFlag2D)

END SUBROUTINE mesh2dOverlay

