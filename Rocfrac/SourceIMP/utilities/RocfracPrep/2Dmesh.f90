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
SUBROUTINE mesh2d(nprocs,iProcs,ichr4)

  use meshdata
  USE Linked_List2

  IMPLICIT NONE

  INCLUDE 'roccomf90.h'

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

  INTEGER :: write_attr, set_option,  ierrFlg, sur_all
  INTEGER :: iNB, iBU, iNI

! obtain function handle ------------------------------------------------------

  write_attr = COM_get_function_handle( 'OUT.write_attribute')
  set_option = COM_get_function_handle( 'OUT.set_option')

! do not append process rank -----------------
  
  CALL COM_call_function( set_option, 2, 'rankwidth', '0')
!  CALL COM_call_function( set_option, 2, 'mode', 'w')

  CALL COM_new_window( surWin )

! --------------------------------------------------
! ------------------------------------------------
! (1)  Ignitable interface (BURNING)
! ---------------------------------------
! -------------------------------------
  
  iNB = iProcs*100 + 1
  iBU = iProcs*100 + 2
  iNI = iProcs*100 + 3
  
  NumNpNew = 0

  IF(MeshType2D.EQ.6)THEN

! renumber elements nodes

     ptrtri6 => SurfMesh_tri6_SF_head
     DO WHILE(ASSOCIATED(ptrtri6))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
           ENDIF
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO

! check for special situation where there is just an edge on the surface

     ptrtri6 => SurfMesh_tri6_SF_head
     DO WHILE(ASSOCIATED(ptrtri6))
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Detected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
              ENDIF
           ENDDO
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO

! write the surface mesh files

     IF(NumEltet2D(1,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
        PRINT*,'ERROR, found in surface mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumEltet2D(1,iProcs),'triangles but no nodes'
        PRINT*,'Stopping'
        STOP
     ENDIF
     
     IF(NumEltet2D(1,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumNpNew,' nodes but no triangles'
     ENDIF


!     IF(NumNPNew.NE.0.AND.NumEltet2D(1,iProcs).NE.0) &
!          write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumEltet2D(1,iProcs)*4,'F=FEPOINT, ET=TRIANGLE'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     ptrtri6 => SurfMesh_tri6_SF_head
     DO WHILE(ASSOCIATED(ptrtri6))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!!$              WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!!$              WRITE(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(glbNdNum) = - NumNpNew
           ENDIF
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO
     ptrtri6 => SurfMesh_tri6_SF_head
     DO WHILE(ASSOCIATED(ptrtri6))
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 NumNpNew = NumNpNew + 1
                 PRINT*,'Dectected surface sliver for pressure, writing'
!                 WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
                 NodeFlag2D(1,NumNpNew) = glbNdNum
                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
                 NodeFlag(glbNdNum) = - NumNpNew
              ENDIF
           ENDDO
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO
     
     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iBU, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iBU, MeshCoor, 3)

     DO i = 1, NumNpNew
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO

! nodeflag  to VolumeNode

     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, 'm')
     CALL COM_set_array( surWin//'.bv', iBU, NodeFlag2D(2,1), 2)


! connectiveity
     TestNumSurfEl = NumEltet2D(1,iProcs)

     CALL COM_set_size( surWin//'.:t6', iBU, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:t6', iBU, ElConnTable, 6)


! elemflag to VolumeElem

     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iBU, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iBU,  1)
     CALL COM_allocate_array(surWin//'.bcflag', iBU, BCSurfFlagOne, 1)
     
     BCSurfFlagOne(1) = 1

     TestNumSurfEl = 0
     ptrtri6 => SurfMesh_tri6_SF_head
     DO WHILE(ASSOCIATED(ptrtri6))
        SurfaceElOnProc = .FALSE.
!!$        do i = 1, 6
!!$           IF(NodeFlag(ptrtri6%ElemData(i)).eq.0) THEN
!!$              SurfaceElOnProc = .false.
!!$              exit
!!$           endif
!!$        enddo

        IF(epart(ptrtri6%ElemData(7)).EQ.iProcs) SurfaceElOnProc = .TRUE.
        IF(SurfaceElOnProc) THEN
!!$           WRITE(4001,'(10i10)') ABS(NodeFlag(ptrtri6%ElemData(1:6))),&
!!$                ElFlag(ptrtri6%ElemData(7)),1,1 ! last two are not really used

           TestNumSurfEl = TestNumSurfEl  + 1
           ELConnTable(1:6,TestNumSurfEl) =  ABS(NodeFlag(ptrtri6%ElemData(1:6)))
           ElFlag_List(TestNumSurfEl) = ElFlag(ptrtri6%ElemData(7))

!!$           write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(1))),ABS(NodeFlag(ptrtri6%ElemData(4))),ABS(NodeFlag(ptrtri6%ElemData(6)))
!!$           write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(4))),ABS(NodeFlag(ptrtri6%ElemData(2))),ABS(NodeFlag(ptrtri6%ElemData(5)))
!!$           write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(6))),ABS(NodeFlag(ptrtri6%ElemData(4))),ABS(NodeFlag(ptrtri6%ElemData(5)))
!!$           write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(3))),ABS(NodeFlag(ptrtri6%ElemData(6))),ABS(NodeFlag(ptrtri6%ElemData(5)))
        ENDIF
        ptrtri6 => ptrtri6%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(1,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in for ignitable surfaces'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumEltet2D(1,iProcs)
        STOP
     ENDIF

     NULLIFY(ptrtri6)


     ! obtain function handle ------------------------------------------------------

!     write_attr = COM_get_function_handle( 'OUT.write_attribute')
!     set_option = COM_get_function_handle( 'OUT.set_option')

!     IF(iProcs.EQ.1)THEN

        CALL COM_call_function( set_option, 2, 'mode', 'w')

!     ELSE
!
!        CALL COM_call_function( set_option, 2, 'mode', 'a')
!        
!     ENDIF

! do not append process rank -----------------

!     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

!     sur_all = Com_get_attribute_handle( surWin//'.all')
!     CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_FS_Ignt.hdf', sur_all,&
!          "solid_surf","00.000000")

!     CALL COM_delete_window( surWin)


  ELSE IF(MeshType2D.EQ.3)THEN

! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_SF_head
     do while(associated(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO

! check for special situation where there is just an edge on the surface

     ptrtri3 => SurfMesh_tri3_SF_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO

! write the surface mesh files

     IF(NumEltet2D(1,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
        PRINT*,'ERROR, found in surface mesh for Solid/Fluid interface'
        PRINT*,'Found ',NumEltet2D(1,iProcs),'triangles but no nodes'
        PRINT*,'Stopping'
        STOP
     ENDIF
     
     IF(NumEltet2D(1,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid interface'
        PRINT*,'Found ',NumNpNew,' nodes but no triangles'
     ENDIF

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_SF_head
     do while(associated(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!!$              WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!!$              write(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(glbNdNum) = - NumNpNew
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO
     ptrtri3 => SurfMesh_tri3_SF_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 NumNpNew = NumNpNew + 1
                 PRINT*,'Dectected surface sliver for pressure, writing'
!                 WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
                 NodeFlag2D(1,NumNpNew) = glbNdNum
                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
                 NodeFlag(glbNdNum) = - NumNpNew
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO

     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iBU, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iBU, MeshCoor, 3)

     DO i = 1, NumNpNew
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO

! nodeflag  to VolumeNode

     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, 'm')
     CALL COM_set_array( surWin//'.bv', iBU, NodeFlag2D(2,1), 2)


! connectiveity
     TestNumSurfEl = NumEltet2D(1,iProcs) ! 0

     CALL COM_set_size( surWin//'.:t3', iBU, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:t3', iBU, ElConnTable, 3)


! elemflag to VolumeElem

     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iBU, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iBU,  1)
     CALL COM_allocate_array(surWin//'.bcflag', iBU, BCSurfFlagOne, 1)
     
     BCSurfFlagOne(1) = 1

     TestNumSurfEl = 0
     ptrtri3 => SurfMesh_tri3_SF_head
     do while(associated(ptrtri3))
        SurfaceElOnProc = .false.
!!$        do i = 1, 6
!!$           IF(NodeFlag(ptrtri6%ElemData(i)).eq.0) THEN
!!$              SurfaceElOnProc = .false.
!!$              exit
!!$           endif
!!$        enddo

        IF(epart(ptrtri3%ElemData(4)).EQ.iProcs) SurfaceElOnProc = .true.
        IF(SurfaceElOnProc) THEN
!!$           WRITE(4001,'(10i10)') ABS(NodeFlag(ptrtri3%ElemData(1:3))),&
!!$                ElFlag(ptrtri3%ElemData(4)),1,1 ! last two are not really used
           
           TestNumSurfEl = TestNumSurfEl  + 1
           ELConnTable(1:3,TestNumSurfEl) =  ABS(NodeFlag(ptrtri3%ElemData(1:3)))
           ElFlag_List(TestNumSurfEl) = ElFlag(ptrtri3%ElemData(4))

!           write(4003,*) ABS(NodeFlag(ptrtri3%ElemData(1))),ABS(NodeFlag(ptrtri3%ElemData(2))),ABS(NodeFlag(ptrtri3%ElemData(3)))
        ENDIF
        ptrtri3 => ptrtri3%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(1,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in Ignitable surface mesh'
        PRINT*,'  different then that in read_patran'
        PRINT*,'  For Fluid Solid Inteface mesh'
        PRINT*,TestNumSurfEl,NumEltet2D(1,iProcs)
        STOP
     ENDIF

     NULLIFY(ptrtri3)
     ! obtain function handle ------------------------------------------------------

!     write_attr = COM_get_function_handle( 'OUT.write_attribute')
!     set_option = COM_get_function_handle( 'OUT.set_option')

!     IF(iProcs.EQ.1)THEN

     CALL COM_call_function( set_option, 2, 'mode', 'w')

!     ELSE
!
!        CALL COM_call_function( set_option, 2, 'mode', 'a')
!        
!     ENDIF

! do not append process rank -----------------

!     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

!     sur_all = Com_get_attribute_handle( surWin//'.all')
!     CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_FS_Ignt.hdf', sur_all,&
!          "solid_surf","00.000000")

!     CALL COM_delete_window( surWin)



  ELSE IF(MeshType2D.EQ.4)THEN ! quad surface mesh ! 8 node brick


! renumber elements nodes

     ptrhex8 => SurfMesh_hex8_SF_head
     DO WHILE(ASSOCIATED(ptrhex8))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 4
           glbNdNum = ptrhex8%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
           ENDIF
        ENDDO
        ptrhex8 => ptrhex8%next
     ENDDO

! write the surface mesh files

     IF(NumElhex2D(1,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
        PRINT*,'ERROR, found in surface mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumElhex2D(1,iProcs),'quads but no nodes'
        PRINT*,'Stopping'
        STOP
     ENDIF
     
     IF(NumElhex2D(1,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumNpNew,' nodes but no quads'
     ENDIF


! write the surface mesh files

 !    WRITE(4001,'(2i10)') NumNpNew,  NumElhex2D(1,iProcs)    

!!$     IF(NumNPNew.NE.0.AND.NumElhex2D(1,iProcs).NE.0) &
!!$          write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumElhex2D(1,iProcs),'F=FEPOINT, ET=QUADRILATERAL'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     TestNumSurfEl = 0

     ptrhex8 => SurfMesh_hex8_SF_head
     do while(associated(ptrhex8))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 4
           glbNdNum = ptrhex8%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!              WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!              write(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(glbNdNum) = - NumNpNew
           endif
        enddo
        ptrhex8 => ptrhex8%next
     ENDDO

     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iBU, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iBU, MeshCoor, 3)

     DO i = 1, NumNpNew
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO

! nodeflag  to VolumeNode

     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, 'm')
     CALL COM_set_array( surWin//'.bv', iBU, NodeFlag2D(2,1), 2)

! connectiveity

     TestNumSurfEl = NumElhex2D(1,iProcs)

     CALL COM_set_size( surWin//'.:q4', iBU, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:q4', iBU, ElConnTable, 4)


! elemflag to VolumeElem

     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iBU, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iBU,  1)
     CALL COM_allocate_array(surWin//'.bcflag', iBU, BCSurfFlagOne, 1)
     
     BCSurfFlagOne(1) = 1

     TestNumSurfEl = 0
     ptrhex8 => SurfMesh_hex8_SF_head
     DO WHILE(ASSOCIATED(ptrhex8))
        SurfaceElOnProc = .FALSE.
!!$        do i = 1, 4
!!$           IF(NodeFlag(ptrhex8%ElemData(i)).eq.0) THEN
!!$              SurfaceElOnProc = .false.
!!$              exit
!!$           endif
!!$        enddo
        IF(epart(ptrhex8%ElemData(5)).EQ.iProcs) SurfaceElOnProc = .TRUE.
        IF(SurfaceElOnProc)THEN
!           WRITE(4001,'(10i10)') ABS(NodeFlag(ptrhex8%ElemData(1:4))),&
!                ElFlag(ptrhex8%ElemData(5)),1,1 ! last two are not really used

           TestNumSurfEl = TestNumSurfEl  + 1
           ELConnTable(1:4,TestNumSurfEl) =  ABS(NodeFlag(ptrhex8%ElemData(1:4)))
           ElFlag_List(TestNumSurfEl) = ElFlag(ptrhex8%ElemData(5))
!           write(4003,*) ABS(NodeFlag(ptrhex8%ElemData(1:4)))
        ENDIF

        ptrhex8 => ptrhex8%next
     ENDDO

    IF(TestNumSurfEl.NE.NumElhex2D(1,iProcs))THEN
        PRINT*,'ERROR, number of quads from link list in Ignitable interface'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumElhex2D(1,iProcs)
        STOP
     ENDIF

     NULLIFY(ptrhex8)

     CALL COM_call_function( set_option, 2, 'mode', 'w')

  ENDIF
  CALL COM_call_function( set_option, 2, 'mode', 'w') ! why would this be called twice ? fix

  sur_all = Com_get_attribute_handle( SurWin//'.all')

  CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh.'//ichr4, sur_all,&
       "isolid","00.000000")

  
  CALL COM_delete_window( surWin)

! --------------------------------------------------
! ------------------------------------------------
! (2)  Non Ignitable interface (INTERACTING)
! ---------------------------------------
! -------------------------------------

  CALL COM_new_window( surWin )


! Reset NodeFlag for non solid/fluid interface

  IF(NumNpNew.NE.0)THEN
     
     DO i = 1, NumNpNew
        NodeFlag(NodeFlag2D(1,i)) = NodeFlag2D(2,i)
     enddo
  ENDIF
  IF(associated(NodeFlag2D)) deallocate(NodeFlag2D)
 
!!$  do i= 1, numnp_prmry
!!$     write(1001,*) i,NodeFlag(i)
!!$  enddo

  IF(MeshType2D.EQ.6)THEN

! renumber elements nodes

     ptrtri6 => SurfMesh_tri6_SF_NonIgnt_head
     DO WHILE(ASSOCIATED(ptrtri6))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
           ENDIF
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO

! check for special situation where there is just an edge on the surface

     ptrtri6 => SurfMesh_tri6_SF_NonIgnt_head
     DO WHILE(ASSOCIATED(ptrtri6))
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Detected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
              ENDIF
           ENDDO
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO

! write the surface mesh files

     IF(NumEltet2D(3,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
        PRINT*,'ERROR, found in surface mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumEltet2D(3,iProcs),'triangles but no nodes'
        PRINT*,'Stopping'
        STOP
     ENDIF
     
     IF(NumEltet2D(3,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumNpNew,' nodes but no triangles'
     ENDIF



!!$     IF(NumNPNew.NE.0.AND.NumEltet2D(3,iProcs).NE.0) &
!!$          write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumEltet2D(3,iProcs)*4,'F=FEPOINT, ET=TRIANGLE'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     ptrtri6 => SurfMesh_tri6_SF_NonIgnt_head
     DO WHILE(ASSOCIATED(ptrtri6))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
              !WRITE(4001,'(3(1X,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
              !WRITE(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(glbNdNum) = - NumNpNew
           ENDIF
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO
     ptrtri6 => SurfMesh_tri6_SF_NonIgnt_head
     DO WHILE(ASSOCIATED(ptrtri6))
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 NumNpNew = NumNpNew + 1
                 PRINT*,'Dectected surface sliver for pressure, writing'
                 !WRITE(4001,'(3(1X,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
                 NodeFlag2D(1,NumNpNew) = glbNdNum
                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
                 NodeFlag(glbNdNum) = - NumNpNew
              ENDIF
           ENDDO
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO
     
     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNB, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iNB, MeshCoor, 3)

     DO i = 1, NumNpNew 
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO


     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!     CALL COM_set_size( surWin//'.bv', iNB, NumNpNew )
     CALL COM_set_array(surWin//'.bv', iNB, NodeFlag2D(2,1), 2)

! connectiveity


     TestNumSurfEl = NumEltet2D(3,iProcs)

     CALL COM_set_size( surWin//'.:t6', iNB, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:t6', iNB, ElConnTable, 6)


     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iNB, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iNB, 1)
     CALL COM_allocate_array(surWin//'.bcflag', iNB, BCSurfFlagZero, 1)
     
     BCSurfFlagZero(1) = 0



     TestNumSurfEl = 0
     ptrtri6 => SurfMesh_tri6_SF_NonIgnt_head
     DO WHILE(ASSOCIATED(ptrtri6))
        SurfaceElOnProc = .FALSE.

        IF(epart(ptrtri6%ElemData(7)).EQ.iProcs) SurfaceElOnProc = .TRUE.
        IF(SurfaceElOnProc) THEN
!!$           WRITE(4001,'(10i10)') ABS(NodeFlag(ptrtri6%ElemData(1:6))),&
!!$                ElFlag(ptrtri6%ElemData(7)),1,1 ! last two are not really used

           TestNumSurfEl = TestNumSurfEl  + 1
           ELConnTable(1:6,TestNumSurfEl) =  ABS(NodeFlag(ptrtri6%ElemData(1:6)))
           ElFlag_List(TestNumSurfEl) = ElFlag(ptrtri6%ElemData(7))
           !write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(1))),ABS(NodeFlag(ptrtri6%ElemData(4))),ABS(NodeFlag(ptrtri6%ElemData(6)))
           !write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(4))),ABS(NodeFlag(ptrtri6%ElemData(2))),ABS(NodeFlag(ptrtri6%ElemData(5)))
           !write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(6))),ABS(NodeFlag(ptrtri6%ElemData(4))),ABS(NodeFlag(ptrtri6%ElemData(5)))
           !write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(3))),ABS(NodeFlag(ptrtri6%ElemData(6))),ABS(NodeFlag(ptrtri6%ElemData(5)))
        ENDIF
        ptrtri6 => ptrtri6%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(3,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in Mesh2d'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumEltet2D(3,iProcs)
        STOP
     ENDIF

     NULLIFY(ptrtri6)


     ! obtain function handle ------------------------------------------------------

!     write_attr = COM_get_function_handle( 'OUT.write_attribute')
!     set_option = COM_get_function_handle( 'OUT.set_option')

!     IF(iProcs.EQ.1)THEN

!        CALL COM_call_function( set_option, 2, 'mode', 'w')

!     ELSE

!        CALL COM_call_function( set_option, 2, 'mode', 'a')
        
!     ENDIF

! do not append process rank -----------------

!     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

!     sur_all = Com_get_attribute_handle( surWin//'.all')

!     CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_FS_NonIgnt.hdf', sur_all,&
!          "solid_surf","00.000000")

!     CALL COM_delete_window( surWin)


  ELSE IF(MeshType2D.EQ.3)THEN

! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_SF_NonIgnt_head
     do while(associated(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO

! check for special situation where there is just an edge on the surface

     ptrtri3 => SurfMesh_tri3_SF_NonIgnt_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
                 
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO

! write the surface mesh files


     !WRITE(4001,'(2i10)') NumNpNew,  NumEltet2D(3,iProcs)

     IF(NumEltet2D(3,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
        PRINT*,'ERROR, found in surface mesh for Solid/Fluid interface'
        PRINT*,'Found ',NumEltet2D(3,iProcs),'triangles but no nodes'
        PRINT*,'Stopping'
        STOP
     ENDIF
     
     IF(NumEltet2D(3,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid interface'
        PRINT*,'Found ',NumNpNew,' nodes but no triangles'
     ENDIF

!     IF(NumNPNew.NE.0.AND.NumEltet2D(3,iProcs).NE.0) &
!          write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumEltet2D(3,iProcs),'F=FEPOINT, ET=TRIANGLE'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_SF_NonIgnt_head
     do while(associated(ptrtri3))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!              WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!              write(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(glbNdNum) = - NumNpNew
           endif
        enddo
        ptrtri3 => ptrtri3%next
     ENDDO
     ptrtri3 => SurfMesh_tri3_SF_NonIgnt_head
     do while(associated(ptrtri3))
        DO i = 1, 3
           glbNdNum = ptrtri3%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 NumNpNew = NumNpNew + 1
                 PRINT*,'Dectected surface sliver for pressure, writing'
!                 WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
                 NodeFlag2D(1,NumNpNew) = glbNdNum
                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
                 NodeFlag(glbNdNum) = - NumNpNew
              ENDIF
           ENDDO
        ENDDO
        ptrtri3 => ptrtri3%next
     ENDDO

     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNB, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iNB, MeshCoor, 3)

     DO i = 1, NumNpNew 
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO


     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!     CALL COM_set_size( surWin//'.bv', iNB, NumNpNew )
     CALL COM_set_array(surWin//'.bv', iNB, NodeFlag2D(2,1), 2)
! connectiveity

     TestNumSurfEl = NumEltet2D(3,iProcs)

     CALL COM_set_size( surWin//'.:t3', iNB, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:t3', iNB, ElConnTable, 3)


     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iNB, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iNB, 1)
     CALL COM_allocate_array(surWin//'.bcflag', iNB, BCSurfFlagZero, 1)
     
     BCSurfFlagZero(1) = 0

     TestNumSurfEl = 0
     ptrtri3 => SurfMesh_tri3_SF_NonIgnt_head
     do while(associated(ptrtri3))
        SurfaceElOnProc = .false.
!!$        do i = 1, 6
!!$           IF(NodeFlag(ptrtri6%ElemData(i)).eq.0) THEN
!!$              SurfaceElOnProc = .false.
!!$              exit
!!$           endif
!!$        enddo

        IF(epart(ptrtri3%ElemData(4)).EQ.iProcs) SurfaceElOnProc = .true.
        IF(SurfaceElOnProc) THEN
!           WRITE(4001,'(10i10)') ABS(NodeFlag(ptrtri3%ElemData(1:3))),&
!                ElFlag(ptrtri3%ElemData(4)),1,1 ! last two are not really used
           TestNumSurfEl = TestNumSurfEl  + 1
           ELConnTable(1:3,TestNumSurfEl) =  ABS(NodeFlag(ptrtri3%ElemData(1:3)))
           ElFlag_List(TestNumSurfEl) = ElFlag(ptrtri3%ElemData(4))
        ENDIF
        ptrtri3 => ptrtri3%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(3,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in Mesh2d'
        PRINT*,'  different then that in read_patran'
        PRINT*,'  For Fluid Solid Inteface mesh'
        PRINT*,TestNumSurfEl,NumEltet2D(3,iProcs)
        STOP
     ENDIF

     NULLIFY(ptrtri3)

     ! obtain function handle ------------------------------------------------------

!     write_attr = COM_get_function_handle( 'OUT.write_attribute')
!     set_option = COM_get_function_handle( 'OUT.set_option')

!     IF(iProcs.EQ.1)THEN

!        CALL COM_call_function( set_option, 2, 'mode', 'w')

!     ELSE

!        CALL COM_call_function( set_option, 2, 'mode', 'a')
        
!     ENDIF

! do not append process rank -----------------

!     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

!     sur_all = Com_get_attribute_handle( surWin//'.all')

!     CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_FS_NonIgnt.hdf', sur_all,&
!          "solid_surf","00.000000")

!     CALL COM_delete_window( surWin)

  ELSE IF(MeshType2D.EQ.4)THEN ! quad surface mesh ! 8 node brick


! renumber elements nodes

     ptrhex8 => SurfMesh_hex8_SF_NonIgnt_head
     DO WHILE(ASSOCIATED(ptrhex8))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 4
           glbNdNum = ptrhex8%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
           ENDIF
        ENDDO
        ptrhex8 => ptrhex8%next
     ENDDO

! write the surface mesh files

     IF(NumElhex2D(3,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
        PRINT*,'ERROR, found in surface mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumElhex2D(3,iProcs),'quads but no nodes'
        PRINT*,'Stopping'
        STOP
     ENDIF
     
     IF(NumElhex2D(3,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid Ignitable interface'
        PRINT*,'Found ',NumNpNew,' nodes but no quads'
     ENDIF

!     WRITE(4001,'(2i10)') NumNpNew,  NumElhex2D(3,iProcs)    

!     IF(NumNPNew.NE.0.AND.NumElhex2D(3,iProcs).NE.0) &
!          write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumElhex2D(3,iProcs),'F=FEPOINT, ET=QUADRILATERAL'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes
     ptrhex8 => SurfMesh_hex8_SF_NonIgnt_head
     do while(associated(ptrhex8))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 4
           glbNdNum = ptrhex8%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!              WRITE(4001,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!              write(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(glbNdNum) = - NumNpNew
           endif
        enddo
        ptrhex8 => ptrhex8%next
     ENDDO

     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNB, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iNB, MeshCoor, 3)

     DO i = 1, NumNpNew 
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO

     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!     CALL COM_set_size( surWin//'.bv', iNB, NumNpNew )
     CALL COM_set_array(surWin//'.bv', iNB, NodeFlag2D(2,1), 2)

! connectiveity


     TestNumSurfEl = NumElhex2D(3,iProcs)

     CALL COM_set_size( surWin//'.:q4', iNB, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:q4', iNB, ElConnTable, 4)


     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iNB, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iNB, 1)
     CALL COM_allocate_array(surWin//'.bcflag', iNB, BCSurfFlagZero, 1)
     
     BCSurfFlagZero(1) = 0


     TestNumSurfEl = 0
     ptrhex8 => SurfMesh_hex8_SF_NonIgnt_head
     DO WHILE(ASSOCIATED(ptrhex8))
        SurfaceElOnProc = .FALSE.

        IF(epart(ptrhex8%ElemData(5)).EQ.iProcs) SurfaceElOnProc = .TRUE.
        IF(SurfaceElOnProc)THEN
!           WRITE(4001,'(10i10)') ABS(NodeFlag(ptrhex8%ElemData(1:4))),&
!                ElFlag(ptrhex8%ElemData(5)),1,1 ! last two are not really used

           TestNumSurfEl = TestNumSurfEl  + 1
           ELConnTable(1:4,TestNumSurfEl) =  ABS(NodeFlag(ptrhex8%ElemData(1:4)))
           ElFlag_List(TestNumSurfEl) = ElFlag(ptrhex8%ElemData(5))

!           write(4003,*) ABS(NodeFlag(ptrhex8%ElemData(1:4)))
        ENDIF
        ptrhex8 => ptrhex8%next
     ENDDO

    IF(TestNumSurfEl.NE.NumElhex2D(3,iProcs))THEN
        PRINT*,'ERROR, number of quads from link list in Non-Ignitable surface mesh'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumElhex2D(3,iProcs)
        STOP
     ENDIF

     NULLIFY(ptrhex8)

  ENDIF

! Reset NodeFlag for non solid/fluid interface

  IF(NumNpNew.NE.0)THEN
     
     DO i = 1, NumNpNew
        NodeFlag(NodeFlag2D(1,i)) = NodeFlag2D(2,i)
     enddo
  ENDIF

  ! write surface window ------------------------

  CALL COM_call_function( set_option, 2, 'mode', 'a')

  sur_all = Com_get_attribute_handle( SurWin//'.all')

  CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh.'//ichr4, sur_all,&
       "isolid","00.000000")


  CALL COM_delete_window( surWin)

  IF(associated(NodeFlag2D)) deallocate(NodeFlag2D)

! --------------------------------------------------
! ------------------------------------------------
! (3)  Non Fluid/Structure interface
! ---------------------------------------
! -------------------------------------

! renumber elements nodes

  CALL COM_new_window( surWin)

  IF(MeshType2D.EQ.6)THEN

     NumNpNew = 0
     NULLIFY(ptrtri6)
! renumber elements nodes
     ptrtri6 => SurfMesh_tri6_S_head

     DO WHILE(ASSOCIATED(ptrtri6))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
           ENDIF
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO


! check for special situation where there is just an edge on the surface

     ptrtri6 => SurfMesh_tri6_S_head
     do while(associated(ptrtri6))
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
                 NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)
              ENDIF
           ENDDO
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO

!     WRITE(4002,'(2i10)') NumNpNew,  NumEltet2D(2,iProcs)
!     write(4003,*) 'ZONE N=', NumNPNew, 'E=',NumEltet2D(2,iProcs),'F=FEPOINT, ET=TRIANGLE'

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting

     NumNpNew = 0
! renumber elements nodes

     ptrtri6 => SurfMesh_tri6_S_head

     do while(associated(ptrtri6))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!              WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
!              write(4003,*) coor(1:3,glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(ptrtri6%ElemData(i)) = - NumNpNew
           endif
        enddo
        ptrtri6 => ptrtri6%next
     ENDDO
     ptrtri6 => SurfMesh_tri6_S_head
     do while(associated(ptrtri6))
        DO i = 1, 6
           glbNdNum = ptrtri6%ElemData(i)
           DO k = 1,MaxNumberOfProcsToShareNode
              IF(ProcNdList(glbNdNum,k).EQ.iProcs.AND.NodeFlag(glbNdNum).GT.0)THEN
                 PRINT*,'Dectected surface sliver for pressure'
                 NumNpNew = NumNpNew + 1
!                 WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
                 NodeFlag2D(1,NumNpNew) = glbNdNum
                 NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
                 NodeFlag(glbNdNum) = - NumNpNew
              ENDIF
           ENDDO
        ENDDO
        ptrtri6 => ptrtri6%next
     ENDDO

     
     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNI, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iNI, MeshCoor, 3)

     DO i = 1, NumNpNew 
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO


     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!     CALL COM_set_size( surWin//'.bv', iNI, NumNpNew )
     CALL COM_set_array( surWin//'.bv', iNI, NodeFlag2D(2,1), 2)


! connectiveity
     TestNumSurfEl = NumEltet2D(2,iProcs)
     PRINT*,'Number of Non-Interacting nodes and elements',NumNpNew,TestNumSurfEl

     CALL COM_set_size( surWin//'.:t6', iNI, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:t6', iNI, ElConnTable, 6)


     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iNI, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iNI,  1)
     CALL COM_allocate_array(surWin//'.bcflag', iNI, BCSurfFlagTwo, 1)
     
     BCSurfFlagTwo(1) = 2

     TestNumSurfEl = 0
     ptrtri6 => SurfMesh_tri6_S_head
     do while(associated(ptrtri6))
        SurfaceElOnProc = .false.
!!$        do i = 1, 6
!!$           IF(NodeFlag(ptrtri6%ElemData(i)).eq.0)THEN
!!$              SurfaceElOnProc = .false.
!!$              exit
!!$           endif
!!$        enddo
        IF(epart(ptrtri6%ElemData(7)).EQ.iProcs) SurfaceElOnProc = .TRUE.
        IF(SurfaceElOnProc) THEN
 !          WRITE(4002,'(10i10)') ABS(NodeFlag(ptrtri6%ElemData(1:6))),&
!             ElFlag(ptrtri6%ElemData(7)),2,1 ! last two are not really used
!!$        IF(SurfaceElOnProc) THEN
!!$           write(4003,*) ABS(NodeFlag(ptrtri6%ElemData(1:6)))
!!$        ENDIF
           TestNumSurfEl = TestNumSurfEl  + 1
           ElConnTable(1:6,TestNumSurfEl) = ABS(NodeFlag(ptrtri6%ElemData(1:6)))
        ENDIF
        ptrtri6 => ptrtri6%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(2,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in Mesh2d'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumEltet2D(2,iProcs)
        STOP
     ENDIF

     NULLIFY(ptrtri6)     ! obtain function handle ------------------------------------------------------

!     write_attr = COM_get_function_handle( 'OUT.write_attribute')
!     set_option = COM_get_function_handle( 'OUT.set_option')

!     IF(iProcs.EQ.1)THEN

!        CALL COM_call_function( set_option, 2, 'mode', 'w')

!     ELSE

!        CALL COM_call_function( set_option, 2, 'mode', 'a')
        
!     ENDIF

! do not append process rank -----------------

!     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

!     sur_all = Com_get_attribute_handle( surWin//'.all')
!     CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_S.hdf', sur_all,&
!          "solid_surf","00.000000")

     

  ELSE IF(MeshType2D.EQ.3)THEN

     NumNpNew = 0
     NULLIFY(ptrtri3)
! renumber elements nodes

     ptrtri3 => SurfMesh_tri3_S_head

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

     ptrtri3 => SurfMesh_tri3_S_head
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

     ptrtri3 => SurfMesh_tri3_S_head

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
     ptrtri3 => SurfMesh_tri3_S_head
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

     
     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNI, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iNI, MeshCoor, 3)

     DO i = 1, NumNpNew 
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO

     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!     CALL COM_set_size( surWin//'.bv', iNI, NumNpNew )
     CALL COM_set_array( surWin//'.bv', iNI, NodeFlag2D(2,1), 2)


! connectiveity
     TestNumSurfEl = NumEltet2D(2,iProcs) !0

     CALL COM_set_size( surWin//'.:t3', iNI, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:t3', iNI, ElConnTable, 3)


     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iNI, ElFlag_List, 1)


     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iNI,1)
     CALL COM_allocate_array(surWin//'.bcflag', iNI, BCSurfFlagTwo, 1)
     
     BCSurfFlagTwo(1) = 2

     TestNumSurfEl = 0
     ptrtri3 => SurfMesh_tri3_S_head
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
           ElFlag_List(TestNumSurfEl) = ElFlag(ptrtri3%ElemData(4))
        ENDIF
        ptrtri3 => ptrtri3%next
     ENDDO

     IF(TestNumSurfEl.NE.NumEltet2D(2,iProcs))THEN
        PRINT*,'ERROR, number of triangles from link list in Mesh2d'
        PRINT*,'  different then that in read_patran'
        PRINT*,TestNumSurfEl,NumEltet2D(2,iProcs)
        STOP
     ENDIF

!     write_attr = COM_get_function_handle( 'OUT.write_attribute')
!     set_option = COM_get_function_handle( 'OUT.set_option')

!     IF(iProcs.EQ.1)THEN

!        CALL COM_call_function( set_option, 2, 'mode', 'w')

!     ELSE

!        CALL COM_call_function( set_option, 2, 'mode', 'a')
        
!     ENDIF

! do not append process rank -----------------

!     CALL COM_call_function( set_option, 2, 'rankwidth', '0')

! write surface window ------------------------

!     sur_all = Com_get_attribute_handle( surWin//'.all')
!     CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh_S.hdf', sur_all,&
!          "solid_surf","00.000000")

     

  ELSE IF(MeshType2D.EQ.4)THEN


     NumNpNew = 0
     NULLIFY(ptrhex8)
! renumber elements nodes

     ptrhex8 => SurfMesh_hex8_S_head

     DO WHILE(ASSOCIATED(ptrhex8))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 4
           glbNdNum = ptrhex8%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor, neg means already renumber node
              NumNpNew = NumNpNew + 1
              NodeFlag(glbNdNum) = -NodeFlag(glbNdNum)  ! overwrite nodeflag here but i need it for non S/F interface
           ENDIF
        ENDDO
        ptrhex8 => ptrhex8%next
     ENDDO


     IF(NumElhex2D(2,iProcs).NE.0.AND.NumNpNew.EQ.0)THEN
        PRINT*,'ERROR, found in surface mesh for Solid/Fluid interface'
        PRINT*,'Found ',NumElhex2D(2,iProcs),'quads but no nodes'
        PRINT*,'Stopping'
        STOP
     ENDIF
     
     IF(NumElhex2D(2,iProcs).EQ.0.AND.NumNpNew.NE.0)THEN
        PRINT*,'Warning: found in sufrace mesh for Solid/Fluid interface'
        PRINT*,'Found ',NumNpNew,' nodes but no quads'
     ENDIF

! write the surface mesh files

!     WRITE(4002,'(2i10)') NumNpNew,  NumElhex2D(2,iProcs)

     NodeFlag(:)= ABS( NodeFlag(:)) ! reset

     ALLOCATE(NodeFlag2D(1:2,1:NumNpNew)) ! for storing what I'm overwriting
     NumNpNew = 0
! renumber elements nodes

     ptrhex8 => SurfMesh_hex8_S_head
     do while(associated(ptrhex8))
! mark that changed surface mesh by changing it to a negative
        DO i = 1, 4
           glbNdNum = ptrhex8%ElemData(i)
           IF( NodeFlag(glbNdNum).GT.0) THEN ! if zero means not on this processor
              NumNpNew = NumNpNew + 1
!              WRITE(4002,'(3(x,e16.9),1i10)') coor(1:3,glbNdNum),NodeFlag(glbNdNum)
              NodeFlag2D(1,NumNpNew) = glbNdNum
              NodeFlag2D(2,NumNpNew) = NodeFlag(glbNdNum)
              NodeFlag(glbNdNum) = - NumNpNew
           endif
        enddo
        ptrhex8 => ptrhex8%next
     ENDDO

     CALL COM_new_attribute( surWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( surWin//'.nc', iNI, NumNpNew )
     CALL COM_allocate_array(surWin//'.nc', iNI, MeshCoor, 3)

     DO i = 1, NumNpNew 
        MeshCoor(1:3,i) = coor(1:3,NodeFlag2D(1,i))
     END DO


     CALL COM_new_attribute( surWin//'.bv', 'n', COM_INTEGER, 1, '')
!     CALL COM_set_size( surWin//'.bv', iNI, NumNpNew )
     CALL COM_set_array( surWin//'.bv', iNI, NodeFlag2D(2,1), 2)


! connectivity
     TestNumSurfEl = NumElhex2D(2,iProcs)
     PRINT*,'Number of Non-Interacting nodes and elements',NumNpNew,TestNumSurfEl

     CALL COM_set_size( surWin//'.:q4', iNI, TestNumSurfEl)
     CALL COM_allocate_array( surWin//'.:q4', iNI, ElConnTable, 4)


     CALL COM_new_attribute( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
     CALL COM_allocate_array(surWin//'.bf2c', iNI, ElFlag_List, 1)

     CALL COM_new_attribute( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( surWin//'.bcflag', iNI,  1)
     CALL COM_allocate_array(surWin//'.bcflag', iNI, BCSurfFlagTwo, 1)
     

     BCSurfFlagTwo(1) = 2
! connectivity
     NULLIFY(ptrhex8)
     TestNumSurfEl = 0

     ptrhex8 => SurfMesh_hex8_S_head
     do while(associated(ptrhex8))
        SurfaceElOnProc = .false.
        IF(epart(ptrhex8%ElemData(5)).EQ.iProcs) SurfaceElOnProc = .TRUE.
        IF(SurfaceElOnProc)THEN
!           WRITE(4002,'(10i10)') ABS(NodeFlag(ptrhex8%ElemData(1:4))),&
!             ElFlag(ptrhex8%ElemData(5)),1,1 ! last two are not really used
           TestNumSurfEl = TestNumSurfEl + 1

           ElConnTable(1:4,TestNumSurfEl) = ABS(NodeFlag(ptrhex8%ElemData(1:4)))
        ENDIF
        ptrhex8 => ptrhex8%next
     ENDDO

     IF(TestNumSurfEl.NE.NumElhex2D(2,iProcs))THEN
        PRINT*,'Error: Number of Solid Surface hex elements in linked list'
        PRINT*,' does not match that given by NumElhex2D(2,iProcs)'
        PRINT*,' In linked list =',TestNumSurfEl
        PRINT*,' NumElhex2D(2,iProcs) =',NumElhex2D(2,iProcs)
        PRINT*,'Stopping'
     ENDIF

     NULLIFY(ptrhex8)

     IF(associated(NodeFlag2D)) deallocate(NodeFlag2D)

  ENDIF

  ! write surface window ------------------------
  CALL COM_call_function( set_option, 2, 'mode', 'a')

  sur_all = Com_get_attribute_handle( SurWin//'.all')

  CALL COM_call_function( write_attr, 4, 'Rocin/SurfMesh.'//ichr4, sur_all,&
       "isolid","00.000000")


  CALL COM_delete_window( surWin)


  IF(associated(NodeFlag2D)) deallocate(NodeFlag2D)

  

end subroutine mesh2d

