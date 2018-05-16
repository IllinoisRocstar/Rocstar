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
MODULE linked_list
  IMPLICIT NONE
  PUBLIC :: coh_insert_head
  PUBLIC :: coh_insert_tail
  PUBLIC :: vol_insert_head
  PUBLIC :: vol_insert_tail
  TYPE, PUBLIC :: coh_info_type
     INTEGER :: mat_coh
     INTEGER, DIMENSION(1:6) :: lmcoh
     INTEGER :: clst_type
  END TYPE coh_info_type
  
  TYPE, PUBLIC :: coh_comm_info_type
     INTEGER, DIMENSION(1:3) :: lmcoh_comm
  ENDTYPE coh_comm_info_type
  
  TYPE, PUBLIC :: vol_info_type
     INTEGER :: mat_vol
     INTEGER, DIMENSION(1:10) :: lmvol
     INTEGER :: iface
     INTEGER :: press
  END TYPE vol_info_type
  
  TYPE, PUBLIC :: comm_info_type
     INTEGER :: NdId
  END TYPE comm_info_type
  
  TYPE :: vol_box_point
     TYPE(vol_node_type), POINTER :: ptr
  ENDTYPE vol_box_point
  
  TYPE :: coh_box_type
     TYPE(coh_node_type), POINTER :: coh_next_np ! => NULL()
  ENDTYPE coh_box_type
  TYPE :: coh_comm_box_type
     TYPE(coh_comm_node_type), POINTER :: coh_comm_next_np ! => NULL()
  ENDTYPE coh_comm_box_type
  
  TYPE :: vol_box_type
     TYPE(vol_node_type), POINTER :: vol_next_np ! => NULL()
  ENDTYPE vol_box_type
  TYPE :: comm_box_type
     TYPE(comm_node_type), POINTER :: comm_next_np ! => NULL()
  ENDTYPE comm_box_type
  
  TYPE :: coh_node_type
     TYPE(coh_info_type) :: coh_info
     TYPE(coh_box_type) :: coh_box
  END TYPE coh_node_type
  TYPE :: coh_comm_node_type
     TYPE(coh_comm_info_type) :: coh_comm_info
     TYPE(coh_comm_box_type) :: coh_comm_box
  END TYPE coh_comm_node_type
  
  TYPE :: vol_node_type
     TYPE(vol_info_type) :: vol_info
     TYPE(vol_box_type) :: vol_box
  END TYPE vol_node_type
  
  
  TYPE :: comm_node_type
     TYPE(comm_info_type) :: comm_info
     TYPE(comm_box_type) :: comm_box
  END TYPE comm_node_type
  
  TYPE :: coh_list_type
     TYPE(coh_node_type), POINTER :: coh_head ! => NULL()
     TYPE(coh_node_type), POINTER :: coh_tail ! => NULL()
  ENDTYPE coh_list_type
  TYPE :: coh_comm_list_type
     TYPE(coh_comm_node_type), POINTER :: coh_comm_head ! => NULL()
     TYPE(coh_comm_node_type), POINTER :: coh_comm_tail ! => NULL()
     INTEGER :: num_comm_nodes ! =  0
  ENDTYPE coh_comm_list_type
  TYPE :: vol_list_type
     TYPE(vol_node_type), POINTER :: vol_head ! => NULL()
     TYPE(vol_node_type), POINTER :: vol_tail ! => NULL()
     INTEGER :: num_border_vol ! = 0 ! Number of bordering volumetric. elements
  ENDTYPE vol_list_type
  
  TYPE :: comm_list_type
     TYPE(comm_node_type), POINTER :: comm_head  !=> NULL()
     TYPE(comm_node_type), POINTER :: comm_tail  !=> NULL()
     INTEGER :: num_border_comm ! = 0 ! Number of bordering volumetric. elements
  ENDTYPE comm_list_type
  
  TYPE :: ElemListType
     TYPE(comm_node_type), POINTER :: comm_head ! => NULL()
     TYPE(comm_node_type), POINTER :: comm_tail ! => NULL()
     
     INTEGER :: NumElperProc ! = 0 ! Number of bordering volumetric. elements
     TYPE(ElemListType), POINTER :: next
  ENDTYPE ElemListType
  
CONTAINS
  
  SUBROUTINE coh_insert_head(arg_b,coh_item)
    TYPE(coh_list_type), TARGET, INTENT(in out) :: arg_b
    TYPE(coh_info_type), intent(in) :: coh_item
    TYPE(coh_node_type), POINTER :: new_node
    
    ALLOCATE(new_node)
    new_node%coh_info = coh_item
    
    IF(.NOT.ASSOCIATED(arg_b%coh_head))THEN
       arg_b%coh_head => new_node
       arg_b%coh_tail => new_node
    ELSE
       new_node%coh_box%coh_next_np => arg_b%coh_head
       arg_b%coh_head => new_node
    ENDIF
    RETURN
  END SUBROUTINE coh_insert_head
  
  SUBROUTINE coh_insert_tail(arg_b,coh_item)
    TYPE(coh_list_type), TARGET, INTENT(in out) :: arg_b
    TYPE(coh_info_type), intent(in) :: coh_item
    TYPE(coh_node_type), POINTER :: new_node
    
    ALLOCATE(new_node)
    new_node%coh_info = coh_item
    
    IF(.NOT.ASSOCIATED(arg_b%coh_head))THEN
       arg_b%coh_head => new_node
       arg_b%coh_tail => new_node
    ELSE
       arg_b%coh_tail%coh_box%coh_next_np => new_node
       arg_b%coh_tail => new_node
    ENDIF
    
    RETURN
  END SUBROUTINE coh_insert_tail
  
  SUBROUTINE Print_coh_list(arg_b)
    TYPE(coh_list_type), TARGET, INTENT(in) :: arg_b
    TYPE(coh_node_type), POINTER :: current
    current => arg_b%coh_head
    DO WHILE (ASSOCIATED(current))
       CALL print_coh_target(current)
       current => current%coh_box%coh_next_np
    ENDDO
    RETURN
  CONTAINS
    SUBROUTINE print_coh_target(arg_np)
      USE meshdata
      TYPE(coh_node_type), POINTER :: arg_np
      IF(IOformat.EQ.0)THEN
         WRITE(4000) arg_np%coh_info
      ELSE
         WRITE(4000,'(8i9)') arg_np%coh_info
      ENDIF
      RETURN
    END SUBROUTINE print_coh_target
  END SUBROUTINE print_coh_list
  SUBROUTINE print_coh_comm(arg_b)
    
    USE meshdata
!
!     Number of communicated nodes needed for R_co calculation 
!     List of the nodes involved in communication for R_co calculation
!      (i.e. what nodes are going to be sent by the 'i' processor)
!
    TYPE(coh_comm_list_type), TARGET, INTENT(in) :: arg_b
    TYPE(coh_comm_node_type), POINTER :: current
    current => arg_b%coh_comm_head
    IF(IOformat.EQ.0)THEN
       WRITE(4000) arg_b%num_comm_nodes*3
    ELSE
       WRITE(4000,*) arg_b%num_comm_nodes*3
    ENDIF
    DO WHILE (ASSOCIATED(current))
       CALL print_coh_comm_target(current)
       current => current%coh_comm_box%coh_comm_next_np
    ENDDO
    RETURN
  CONTAINS
    SUBROUTINE print_coh_comm_target(arg_np)
      USE meshdata
      TYPE(coh_comm_node_type), POINTER :: arg_np
      IF(IOformat.EQ.0)THEN
         WRITE(4000) arg_np%coh_comm_info
      ELSE
         WRITE(4000,'(6i9)') arg_np%coh_comm_info
      ENDIF
      RETURN
    END SUBROUTINE print_coh_comm_target
  END SUBROUTINE print_coh_comm
  
  SUBROUTINE vol_insert_head(arg_b,vol_item)
    TYPE(vol_list_type), TARGET, INTENT(in out) :: arg_b
    TYPE(vol_info_type), intent(in) :: vol_item
    TYPE(vol_node_type), POINTER :: new_node
    
    ALLOCATE(new_node)
    new_node%vol_info = vol_item
    
    IF(.NOT.ASSOCIATED(arg_b%vol_head))THEN
       arg_b%vol_head => new_node
       arg_b%vol_tail => new_node
       arg_b%num_border_vol = 1
    ELSE
       new_node%vol_box%vol_next_np => arg_b%vol_head
       arg_b%vol_head => new_node
       arg_b%num_border_vol = arg_b%num_border_vol + 1
    ENDIF
    RETURN
  END SUBROUTINE vol_insert_head
  
  SUBROUTINE vol_insert_tail(arg_b,vol_item)
    TYPE(vol_list_type), TARGET, INTENT(in out) :: arg_b
    TYPE(vol_info_type), intent(in) :: vol_item
    TYPE(vol_node_type), POINTER :: new_node
    
    ALLOCATE(new_node)
    new_node%vol_info = vol_item
    
    IF(.NOT.ASSOCIATED(arg_b%vol_head))THEN
       arg_b%vol_head => new_node
       arg_b%vol_tail => new_node
    ELSE
       arg_b%vol_tail%vol_box%vol_next_np => new_node
       arg_b%vol_tail => new_node
    ENDIF
    
    RETURN
  END SUBROUTINE vol_insert_tail
  
  SUBROUTINE addcommNd(arg_b,comm_item)
    TYPE(comm_list_type), TARGET, INTENT(in out) :: arg_b
    integer, intent(in) :: comm_item
    TYPE(comm_node_type), POINTER :: new_node
    
    ALLOCATE(new_node)
    new_node%comm_info%NdId = comm_item
    
    IF(.NOT.ASSOCIATED(arg_b%comm_head))THEN
       arg_b%comm_head => new_node
       arg_b%comm_tail => new_node
       arg_b%num_border_comm = 1
    ELSE
       arg_b%comm_tail%comm_box%comm_next_np => new_node
       arg_b%comm_tail => new_node
       arg_b%num_border_comm = arg_b%num_border_comm + 1
    ENDIF
    
    RETURN
  END SUBROUTINE addcommNd
  
  SUBROUTINE AddElPart(arg_b,comm_item)
    TYPE(ElemListType), TARGET, INTENT(in out) :: arg_b
    integer, intent(in) :: comm_item
    TYPE(comm_node_type), POINTER :: new_node
    
    ALLOCATE(new_node)
    new_node%comm_info%NdId = comm_item
    
    IF(.NOT.ASSOCIATED(arg_b%comm_head))THEN
       arg_b%comm_head => new_node
       arg_b%comm_tail => new_node
       arg_b%NumElperProc = 1
    ELSE
       arg_b%comm_tail%comm_box%comm_next_np => new_node
       arg_b%comm_tail => new_node
       arg_b%NumElperProc = arg_b%NumElperProc + 1
    ENDIF
    
    RETURN
  END SUBROUTINE AddElPart
  
  
  SUBROUTINE vol_insert_ptr(arg_b,ik1_c4)
    TYPE(vol_list_type), TARGET, INTENT(in) :: arg_b
    TYPE(vol_box_point) :: ik1_c4
    
    ik1_c4%ptr => arg_b%vol_tail
    
    RETURN
  END SUBROUTINE vol_insert_ptr
  
  SUBROUTINE coh_comm(arg_b,coh_comm_item)
    TYPE(coh_comm_list_type), TARGET, INTENT(in out) :: arg_b
    TYPE(coh_comm_info_type), INTENT(in) :: coh_comm_item
    TYPE(coh_comm_node_type), POINTER :: new_node
    
    ALLOCATE(new_node)
    new_node%coh_comm_info = coh_comm_item
    
    IF(.NOT.ASSOCIATED(arg_b%coh_comm_head))THEN
       arg_b%coh_comm_head => new_node
       arg_b%coh_comm_tail => new_node
       arg_b%num_comm_nodes = 1 
    ELSE
       arg_b%coh_comm_tail%coh_comm_box%coh_comm_next_np => new_node
       arg_b%coh_comm_tail => new_node
       arg_b%num_comm_nodes = arg_b%num_comm_nodes + 1 
    ENDIF
    
    RETURN
  END SUBROUTINE coh_comm
  
  SUBROUTINE Print_vol_list(arg_b,imat, Inum, ElemCount)
    TYPE(vol_list_type), TARGET, INTENT(in) :: arg_b
    TYPE(vol_node_type), POINTER :: current
    integer :: imat, icnter, Inum, i
    Integer :: ElemCount
    current => arg_b%vol_head
!    icnter = 0

    ! ... ElemCount added to keep track of elements put into this processors 
    ! ... connectivity table (ElConnTable) regardless of the 
    ! ... material.  COstoich 10/27/09
    icnter = ElemCount

!   DO WHILE (ASSOCIATED(current))
    DO i = 1, Inum
       CALL print_vol_target(current,imat,icnter)
       current => current%vol_box%vol_next_np
    ENDDO

    RETURN
  CONTAINS
    SUBROUTINE Print_vol_Target(arg_np,imat,icnter)
      USE meshdata
      integer :: imat,icnter
      TYPE(vol_node_type), POINTER :: arg_np
      IF(IOformat.EQ.0)THEN
         WRITE(4000) arg_np%vol_info
      ELSE

!!$         WRITE(4000,'(15i9)') imat, &
!!$              arg_np%vol_info%lmvol(1:numvertx) 

         icnter = icnter + 1
         ElFlag(arg_np%vol_info%press) = icnter

         matType(icnter) = imat
         ElConnTable(1:numvertx,icnter) = arg_np%vol_info%lmvol(1:numvertx)


      ENDIF
      RETURN
    END SUBROUTINE print_vol_target
  END SUBROUTINE Print_vol_list
  
  
  SUBROUTINE Print_comm_list(arg_b,ip, iaux)
    
    integer :: ip
    integer :: iaux
    integer :: icnt
    TYPE(comm_list_type), TARGET, INTENT(in) :: arg_b
    TYPE(comm_node_type), POINTER :: current

    icnt = 0
    current => arg_b%comm_head
    DO WHILE (ASSOCIATED(current))
       icnt = icnt + 1
       CALL print_comm_target(current,ip,iaux)
       IF(icnt.EQ.arg_b%num_border_comm) EXIT ! should not need to do this, fix pointer
       current => current%comm_box%comm_next_np
    ENDDO
    RETURN
  CONTAINS
    SUBROUTINE Print_comm_Target(arg_np,ip,iaux)
      USE meshdata
      
      integer :: ip, iaux
      
      TYPE(comm_node_type), POINTER :: arg_np

!      NodesToCommunicate_cnt = NodesToCommunicate_cnt + 1
!      print*,'NodesToCommunicate_cnt',NodesToCommunicate_cnt
!      NodesToCommunicate(NodesToCommunicate_cnt) = NodeFlag(arg_np%comm_info%NdId) 
 
      iaux = iaux + 1

      Pconn_Comm(iaux) = NodeFlag(arg_np%comm_info%NdId)


      RETURN
    END SUBROUTINE print_comm_target
  END SUBROUTINE Print_comm_list
  
  SUBROUTINE LinkedList(arg_b,ip)
    
    integer :: ip
    
    TYPE(ElemListType), TARGET, INTENT(in) :: arg_b
    TYPE(comm_node_type), POINTER, SAVE :: current
    current => arg_b%comm_head
!      DO WHILE (ASSOCIATED(current))
    CALL print_comm_target2(current,ip)
    current => current%comm_box%comm_next_np
!      ENDDO
    RETURN
  CONTAINS
    SUBROUTINE Print_comm_Target2(arg_np,ip)
      
      integer :: ip
      
      TYPE(comm_node_type), POINTER :: arg_np
      
      
      ip = arg_np%comm_info%NdId
      
      RETURN
    END SUBROUTINE print_comm_target2
  END SUBROUTINE LinkedList
  
  
END MODULE linked_list

