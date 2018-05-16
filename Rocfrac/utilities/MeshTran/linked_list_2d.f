       MODULE linked_list_2d
         IMPLICIT NONE
         PUBLIC :: vol_insert_head
         PUBLIC :: vol_insert_tail
         TYPE, PUBLIC :: coh_info_type
            INTEGER :: mat_coh
            INTEGER, DIMENSION(1:6) :: lmcoh
            INTEGER :: clst_type
         END TYPE coh_info_type

         TYPE, PUBLIC :: coh_comm_info_type
            INTEGER, DIMENSION(1:3) :: lmcoh_comm
         ENDTYPE

         TYPE, PUBLIC :: lst_comm_info_type
            INTEGER :: lst_node_comm
         ENDTYPE

         TYPE, PUBLIC :: cohlst_comm_info_type
            INTEGER :: cohlst_node_comm
         ENDTYPE
C
C -- derived volumetric type
C
         TYPE, PUBLIC :: vol_info_type
            INTEGER, DIMENSION(1:6) :: lmvol    ! list of nodes
            INTEGER :: mapelem
            INTEGER :: mat_vol                  ! material type
            INTEGER :: flag

         END TYPE vol_info_type

         TYPE :: vol_box_point
           TYPE(vol_node_type), POINTER :: ptr
         ENDTYPE vol_box_point        

         TYPE :: coh_box_type
            TYPE(coh_node_type), POINTER :: coh_next_np ! => NULL()
         ENDTYPE coh_box_type
         TYPE :: coh_comm_box_type
            TYPE(coh_comm_node_type), POINTER :: coh_comm_next_np ! => NULL()
         ENDTYPE coh_comm_box_type

         TYPE :: lst_comm_box_type
            TYPE(lst_comm_node_type), POINTER :: lst_comm_next_np ! => NULL()
         ENDTYPE lst_comm_box_type

         TYPE :: cohlst_comm_box_type
            TYPE(cohlst_comm_node_type), POINTER :: cohlst_comm_next_np ! => NULL()
         ENDTYPE cohlst_comm_box_type

         TYPE :: vol_box_type
            TYPE(vol_node_type), POINTER :: vol_next_np ! => NULL()
         ENDTYPE vol_box_type

         TYPE :: coh_node_type
            TYPE(coh_info_type) :: coh_info
            TYPE(coh_box_type) :: coh_box
         END TYPE coh_node_type
         TYPE :: coh_comm_node_type
            TYPE(coh_comm_info_type) :: coh_comm_info
            TYPE(coh_comm_box_type) :: coh_comm_box
         END TYPE coh_comm_node_type
         TYPE :: lst_comm_node_type
            TYPE(lst_comm_info_type) :: lst_comm_info
            TYPE(lst_comm_box_type) :: lst_comm_box
         END TYPE lst_comm_node_type
         TYPE :: cohlst_comm_node_type
            TYPE(cohlst_comm_info_type) :: cohlst_comm_info
            TYPE(cohlst_comm_box_type) :: cohlst_comm_box
         END TYPE cohlst_comm_node_type

         TYPE :: vol_node_type
            TYPE(vol_info_type) :: vol_info
            TYPE(vol_box_type) :: vol_box
         END TYPE vol_node_type

         TYPE :: coh_list_type
            TYPE(coh_node_type), POINTER :: coh_head ! => NULL() 
            TYPE(coh_node_type), POINTER :: coh_tail ! => NULL()
         ENDTYPE
         TYPE :: coh_comm_list_type
            TYPE(coh_comm_node_type), POINTER :: coh_comm_head ! => NULL() 
            TYPE(coh_comm_node_type), POINTER :: coh_comm_tail ! => NULL()
            INTEGER :: num_comm_nodes ! =  0
         ENDTYPE
         TYPE :: lst_comm_list_type
            TYPE(lst_comm_node_type), POINTER :: lst_comm_head ! => NULL()
            TYPE(lst_comm_node_type), POINTER :: lst_comm_tail ! => NULL()
            INTEGER :: lst_num_comm_nodes ! =  0
         ENDTYPE lst_comm_list_type
         TYPE :: cohlst_comm_list_type
            TYPE(cohlst_comm_node_type), POINTER :: cohlst_comm_head ! => NULL()
            TYPE(cohlst_comm_node_type), POINTER :: cohlst_comm_tail ! => NULL()
            INTEGER :: cohlst_num_comm_nodes ! =  0
         ENDTYPE cohlst_comm_list_type
         TYPE :: vol_list_type
            TYPE(vol_node_type), POINTER :: vol_head ! => NULL(), 
            TYPE(vol_node_type), POINTER :: vol_tail ! => NULL()
            INTEGER :: num_border_vol ! = 0 ! Number of bordering volumetric. elements
         ENDTYPE
      CONTAINS

          SUBROUTINE print_lst_comm(arg_b)
C
C     1) Number of communicated nodes needed for R_in calculation 
C     2) List of nodes involved in communication for R_in calculation
C      (i.e. what nodes are going to be sent by the 'i' processor)
C
            TYPE(lst_comm_list_type), TARGET, INTENT(in) :: arg_b
            TYPE(lst_comm_node_type), POINTER :: current
            current => arg_b%lst_comm_head

            WRITE(4000,*) arg_b%lst_num_comm_nodes
            DO WHILE (ASSOCIATED(current))
               CALL print_lst_comm_target(current)
               current => current%lst_comm_box%lst_comm_next_np
            ENDDO
            RETURN
         CONTAINS
            SUBROUTINE print_lst_comm_target(arg_np)
               TYPE(lst_comm_node_type), POINTER :: arg_np
               WRITE(4000,*) arg_np%lst_comm_info
               RETURN
            END SUBROUTINE print_lst_comm_target
          END SUBROUTINE print_lst_comm

         SUBROUTINE vol_insert_head(arg_b,vol_item)
            ! Add the volumetric element to the beginning of the list

            ! The volumetric list 
            TYPE(vol_list_type), TARGET, INTENT(in out) :: arg_b
            ! Item to be added to the list
            TYPE(vol_info_type), intent(in) :: vol_item
            ! Dummy arguement
            TYPE(vol_node_type), POINTER :: new_node

            ALLOCATE(new_node)
            new_node%vol_info = vol_item

            IF(.NOT.ASSOCIATED(arg_b%vol_head))THEN
               ! List is empty
               arg_b%vol_head => new_node
               arg_b%vol_tail => new_node
c               arg_b%num_border_vol = 1
            ELSE
               ! List is not empty
               new_node%vol_box%vol_next_np => arg_b%vol_head
               arg_b%vol_head => new_node
c               arg_b%num_border_vol = arg_b%num_border_vol + 1
            ENDIF
            RETURN
         END SUBROUTINE vol_insert_head

         SUBROUTINE vol_insert_tail(arg_b,vol_item)
            ! Add the volumetric element to the end of the list
            TYPE(vol_list_type), TARGET, INTENT(in out) :: arg_b
            ! The volumetric list
            TYPE(vol_info_type), intent(in) :: vol_item
            !  Dummy arguement
            TYPE(vol_node_type), POINTER :: new_node

            ALLOCATE(new_node)
            new_node%vol_info = vol_item

            IF(.NOT.ASSOCIATED(arg_b%vol_head))THEN
               arg_b%vol_head => new_node
               arg_b%vol_tail => new_node
               arg_b%num_border_vol = 1
            ELSE
               arg_b%vol_tail%vol_box%vol_next_np => new_node
               arg_b%vol_tail => new_node
               arg_b%num_border_vol = arg_b%num_border_vol + 1
            ENDIF

            RETURN
         END SUBROUTINE vol_insert_tail

         SUBROUTINE vol_insert_ptr(arg_b,ik1_c4)
            TYPE(vol_list_type), TARGET, INTENT(in) :: arg_b
            TYPE(vol_box_point) :: ik1_c4

            ik1_c4%ptr => arg_b%vol_tail

            RETURN
         END SUBROUTINE vol_insert_ptr

         SUBROUTINE lst_comm(arg_b,lst_comm_item)
            ! Add to list of volumetric element communication
            TYPE(lst_comm_list_type), TARGET, INTENT(in out) :: arg_b
            TYPE(lst_comm_info_type), INTENT(in) :: lst_comm_item
            TYPE(lst_comm_node_type), POINTER :: new_node

            ALLOCATE(new_node)
            new_node%lst_comm_info = lst_comm_item
            
            IF(.NOT.ASSOCIATED(arg_b%lst_comm_head))THEN
               arg_b%lst_comm_head => new_node
               arg_b%lst_comm_tail => new_node
               arg_b%lst_num_comm_nodes = 1 
            ELSE
               arg_b%lst_comm_tail%lst_comm_box%lst_comm_next_np => new_node
               arg_b%lst_comm_tail => new_node
               arg_b%lst_num_comm_nodes = arg_b%lst_num_comm_nodes + 1 
            ENDIF

            RETURN
         END SUBROUTINE lst_comm

         SUBROUTINE Print_vol_list(arg_b,numvertx,iunit)
            TYPE(vol_list_type), TARGET, INTENT(in) :: arg_b
            TYPE(vol_node_type), POINTER :: current
            INTEGER :: numvertx,iunit
            current => arg_b%vol_head
            DO WHILE (ASSOCIATED(current))
               CALL print_vol_target(current,numvertx,iunit)
               current => current%vol_box%vol_next_np
            ENDDO
            RETURN
         CONTAINS
            SUBROUTINE Print_vol_Target(arg_np,numvertx,iunit)
               TYPE(vol_node_type), POINTER :: arg_np
               INTEGER :: numvertx,iunit
c               PRINT*,arg_np%vol_info
               WRITE(iunit,'(10i10)') arg_np%vol_info%lmvol(1:numvertx), 
     $              arg_np%vol_info%mapelem,
     $              arg_np%vol_info%flag, arg_np%vol_info%mat_vol
               RETURN
            END SUBROUTINE print_vol_target
          END SUBROUTINE Print_vol_list

         SUBROUTINE Print_vol_list_tec(arg_b)
            TYPE(vol_list_type), TARGET, INTENT(in) :: arg_b
            TYPE(vol_node_type), POINTER :: current
            current => arg_b%vol_head
            DO WHILE (ASSOCIATED(current))
               CALL print_vol_target(current)
               current => current%vol_box%vol_next_np
            ENDDO
            RETURN
         CONTAINS
            SUBROUTINE Print_vol_Target(arg_np)
               TYPE(vol_node_type), POINTER :: arg_np
               WRITE(2001,'(3i10)') arg_np%vol_info%lmvol(1:3)
               RETURN
            END SUBROUTINE print_vol_target
          END SUBROUTINE Print_vol_list_tec


      END MODULE linked_list_2d
