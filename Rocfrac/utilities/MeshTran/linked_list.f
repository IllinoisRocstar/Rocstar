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
         ENDTYPE

         TYPE, PUBLIC :: vol_info_type
            INTEGER :: mat_vol
            INTEGER, DIMENSION(1:10) :: lmvol
            INTEGER :: iface
            REAL*8 :: press
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
         ENDTYPE
         TYPE :: coh_comm_list_type
            TYPE(coh_comm_node_type), POINTER :: coh_comm_head ! => NULL()
            TYPE(coh_comm_node_type), POINTER :: coh_comm_tail ! => NULL()
            INTEGER :: num_comm_nodes ! =  0
         ENDTYPE
         TYPE :: vol_list_type
            TYPE(vol_node_type), POINTER :: vol_head ! => NULL()
            TYPE(vol_node_type), POINTER :: vol_tail ! => NULL()
            INTEGER :: num_border_vol ! = 0 ! Number of bordering volumetric. elements
         ENDTYPE

         TYPE :: comm_list_type
            TYPE(comm_node_type), POINTER :: comm_head ! => NULL()
            TYPE(comm_node_type), POINTER :: comm_tail ! => NULL()
            INTEGER :: num_border_comm ! = 0 ! Number of bordering volumetric. elements
         ENDTYPE
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
C
C     Number of communicated nodes needed for R_co calculation 
C     List of the nodes involved in communication for R_co calculation
C      (i.e. what nodes are going to be sent by the 'i' processor)
C
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

         SUBROUTINE Print_vol_list(arg_b,imat)
            TYPE(vol_list_type), TARGET, INTENT(in) :: arg_b
            TYPE(vol_node_type), POINTER :: current
            integer :: imat
            current => arg_b%vol_head
            DO WHILE (ASSOCIATED(current))
               CALL print_vol_target(current,imat)
               current => current%vol_box%vol_next_np
            ENDDO
            RETURN
         CONTAINS
            SUBROUTINE Print_vol_Target(arg_np,imat)
               USE meshdata
               integer :: imat
               TYPE(vol_node_type), POINTER :: arg_np
               IF(IOformat.EQ.0)THEN
                  WRITE(4000) arg_np%vol_info
               ELSE
                  
                  WRITE(4000,'(15i9)') imat,
     $                 arg_np%vol_info%lmvol(1:numvertx),
     $                 arg_np%vol_info%iface,0
               ENDIF
               RETURN
             END SUBROUTINE print_vol_target
      END SUBROUTINE Print_vol_list
      SUBROUTINE Print_comm_list(arg_b)
        TYPE(comm_list_type), TARGET, INTENT(in) :: arg_b
        TYPE(comm_node_type), POINTER :: current
        current => arg_b%comm_head
        DO WHILE (ASSOCIATED(current))
           CALL print_comm_target(current)
           current => current%comm_box%comm_next_np
        ENDDO
        RETURN
      CONTAINS
        SUBROUTINE Print_comm_Target(arg_np)
          USE meshdata
          TYPE(comm_node_type), POINTER :: arg_np
          IF(IOformat.EQ.0)THEN
             WRITE(4000) arg_np%comm_info%NdId
          ELSE
          !   print*,arg_np%comm_info%NdId
             WRITE(4000,'(15i9)') arg_np%comm_info%NdId
          ENDIF
          RETURN
        END SUBROUTINE print_comm_target
      END SUBROUTINE Print_comm_list

      END MODULE linked_list
