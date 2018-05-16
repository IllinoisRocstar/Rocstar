      SUBROUTINE coh_opposite_side(numlst,ni,node1,node2,node3,node4,
     $     gnode1,gnode2,gnode3,gnode4,neigh,j,i,coh_item)

      USE linked_list
      IMPLICIT NONE

      INTEGER :: ni,node1,node2,node3,node4,i,j,numlst
      INTEGER :: gnode1,gnode2,gnode3,gnode4
      INTEGER, DIMENSION(1:4,1:numlst) :: neigh
      TYPE(coh_info_type), intent(in out) :: coh_item

      IF(neigh(1,neigh(j,i)).EQ.i)THEN
         IF(gnode2.EQ.ni)THEN
            coh_item%lmcoh(4) = node2
            coh_item%lmcoh(5) = node3
            coh_item%lmcoh(6) = node4
         ELSEIF(gnode3.EQ.ni)THEN
            coh_item%lmcoh(4) = node3
            coh_item%lmcoh(5) = node4
            coh_item%lmcoh(6) = node2
         ELSEIF(gnode4.EQ.ni)THEN
            coh_item%lmcoh(4) = node4
            coh_item%lmcoh(5) = node2
            coh_item%lmcoh(6) = node3
         ENDIF
      ELSE IF(neigh(2,neigh(j,i)).EQ.i)THEN
         IF(gnode1.EQ.ni)THEN
            coh_item%lmcoh(4) = node1
            coh_item%lmcoh(5) = node4
            coh_item%lmcoh(6) = node3
         ELSEIF(gnode3.EQ.ni)THEN
            coh_item%lmcoh(4) = node3
            coh_item%lmcoh(5) = node1
            coh_item%lmcoh(6) = node4
         ELSEIF(gnode4.EQ.ni)THEN
            coh_item%lmcoh(4) = node4
            coh_item%lmcoh(5) = node3
            coh_item%lmcoh(6) = node1
         ENDIF
      ELSEIF(neigh(3,neigh(j,i)).EQ.i)THEN
         IF(gnode1.EQ.ni)THEN
            coh_item%lmcoh(4) = node1
            coh_item%lmcoh(5) = node2
            coh_item%lmcoh(6) = node4
         ELSEIF(gnode2.EQ.ni)THEN
            coh_item%lmcoh(4) = node2
            coh_item%lmcoh(5) = node4
            coh_item%lmcoh(6) = node1
         ELSEIF(gnode4.EQ.ni)THEN
            coh_item%lmcoh(4) = node4
            coh_item%lmcoh(5) = node1
            coh_item%lmcoh(6) = node2
         ENDIF
      ELSEIF(neigh(4,neigh(j,i)).EQ.i)THEN
         IF(gnode1.EQ.ni)THEN
            coh_item%lmcoh(4) = node1
            coh_item%lmcoh(5) = node3
            coh_item%lmcoh(6) = node2
         ELSEIF(gnode2.EQ.ni)THEN
            coh_item%lmcoh(4) = node2
            coh_item%lmcoh(5) = node1
            coh_item%lmcoh(6) = node3
         ELSEIF(gnode3.EQ.ni)THEN
            coh_item%lmcoh(4) = node3
            coh_item%lmcoh(5) = node2
            coh_item%lmcoh(6) = node1
         ENDIF
      ENDIF
      
      RETURN 
      END
