
!!****
!!
!!  NAME
!!     KtNumNnz
!!
!!  FUNCTION
!!     This subroutine uses the node numbering combined with
!!     the connectivity table in order to calculate the number
!!     of nonzero entries in the stiffness matrix due only to
!!     the elements and nodes on this processor.
!!
!!  INPUTS
!!     NumNp -- Total number of nodes that this proc knows about
!!     NumEl -- Total number of elements that this proc knows about
!!     ElConnVol -- Connectivity table
!!
!!  OUTPUTS
!!     nnz -- The number of nonzeros in the K matrix
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE KtNumNnz(NumNp,NumEl,ElConnVol,nnz)

  USE Precision
  USE implicit_global

  IMPLICIT NONE

  ! ... Input variables
  INTEGER :: NumNp, NumEl
  INTEGER :: ElConnVol(1:8,1:NumEl)

  ! ... Output variables
  INTEGER :: nnz

  ! ... local variables
  INTEGER :: i, j, k, m, n, counter
  INTEGER :: innz(1:56), jnnz(1:26)


!
!   Count the nonzeros that should be in the K matrix
!

  ! Loop through nodes counting off-diagonal blocks of nonzeros
  nnz = 0
  DO i = 1, NumNp

     ! Initialize connectivity variables
     innz(:) = -1
     jnnz(:) = -1
     
     ! ... Construct nodal connectivity vector
     ! ... there are a maximum of 8 elements with node i 
     ! ... innz caries the nodes of the elements associated with node i
     DO j = 1, NumEl
        DO k = 1, 8
           IF (ElConnVol(k,j) == i) THEN
              DO m = 1, 8
                 IF (k /= m) THEN
                    DO n = 1, 56
                       IF (innz(n) == -1) THEN
                          innz(n) = ElConnVol(m,j)
                          EXIT
                       ENDIF
                    ENDDO
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     
     ! ...Eliminate duplicates from nodal connectivity vector
     ! ... Of the 56 in innz, 3 * 12 - 6 are duplicates, (30)
     ! ... jnnz caries the unique entries 
     DO j = 1, 56
        IF (innz(j) /= -1) THEN
           DO k = 1, 26
              IF (jnnz(k) /= -1) THEN
                 IF (jnnz(k) == innz(j)) THEN
                    EXIT
                 ENDIF
              ELSE
                 jnnz(k) = innz(j)

                 ! ... the number of non-zeros is incremented each time
                 ! ... a unique entry is found in innz
                 nnz = nnz + 1
                 EXIT
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  ! ... Take into account the blocks on the diagonal
  ! ... these were the nodes that were looped over above
  ! ... and not added to the connectivity vectors
  nnz = nnz + NumNp

  ! ... Convert number of blocks to number of nonzeros
  ! ... three nodes corresponds to a 3x3 stiffness matrix
  ! ... multiplied by (1) because of 1 degree of freedom
 ! nnz = (1)*nnz 

  
END SUBROUTINE KtNumNnz




!!****
!!
!!  NAME
!!     KtInitCRS
!!
!!  FUNCTION
!!     This subroutine uses the node numbering combined with
!!     the connectivity table and number of nonzeros in the K
!!     matrix to assemble the pre-allocated compressed row
!!     storage arrays.  It does this by putting a definite
!!     (not assumed) zero at each point that can contain a
!!     nonzero value.  The actual values in the K matrix are
!!     added in the LocThermCap_v3d8 subroutine.
!!
!!  INPUTS
!!     NumNp -- Total number of nodes that this proc knows about
!!     NumEl -- Total number of elements that this proc knows about
!!     ElConnVol -- Connectivity table
!!     nnz -- The number of nonzeros in the K matrix
!!
!!  OUTPUTS
!!     rp -- The row mapping vector
!!     cval -- The collumn mapping vector
!!     aval -- The value vector (which will be returned as all zeros)
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE KtInitCRS(NumNp,NumEl,ElConnVol,nnz,rp,cval,aval)


  USE Precision
  USE implicit_global

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  ! ... Input variables
  INTEGER :: NumNp, NumEl
  INTEGER :: ElConnVol(1:8,1:NumEl)
  INTEGER :: nnz

  ! ... Output variables
  INTEGER :: rp(1:GNumNp+1)
  INTEGER :: cval(1:nnz)
  INTEGER :: aval(1:nnz)

  ! ... local variables
  INTEGER :: i, j, k, m, n, counter, ii
  INTEGER :: innz(1:56), jnnz(1:27)


  ! ... Initialize CRS variables
  aval(:) = 0.0
  cval(:) = -1
  rp(:) = -1
  rp(1) = 0

  ! ... Loop through nodes
  counter = 0
  DO i = 1, GNumNp

     ! ... If this node is on this processor
     IF (Global2Local(i) /= -1) THEN

        ii = Global2Local(i)

        ! ... Initialize connectivity variables
        innz(:) = -1
        jnnz(:) = -1
        
        ! ... Make sure diagonal terms are included
        jnnz(27) = ii
        
        ! ... Construct nodal connectivity vector
        ! ... there are a maximum of 8 elements that contain node ii 
        ! ... innz caries the nodes of the elements associated with node ii
        DO j = 1, NumEl
           DO k = 1, 8
              ! ... Find the element and element local node number corresponding
              ! ... to process local node number ii.  This is identified by 
              ! ... integer k.
              IF (ElConnVol(k,j) == ii) THEN
                 DO m = 1, 8
                    ! ... Find other 7 nodes (m) in element j.
                    IF (k /= m) THEN
                       DO n = 1, 56
                          ! ... Store the process local nodes belonging to element j in 
                          ! ... the vector innz.  There are a maximum of 56 nodes that
                          ! ... can belong to the same elements as node ii.
                          ! ... (8 elements) X (7 non-ii nodes per element) = 56 nodes
                          IF (innz(n) == -1) THEN
                             innz(n) = ElConnVol(m,j)
                             EXIT
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
        ENDDO


        ! ...Eliminate duplicates from nodal connectivity vector
        ! ... Of the 56 in innz, 3 * 12 - 6 are duplicates, (30)
        ! ... jnnz caries the 26 unique entries + 1 for the diagonal (ii)
        ! ... which is not counted
        DO j = 1, 56
           ! ... Look though the 56 entries in innz for a non -1.
           IF (innz(j) /= -1) THEN
              ! ... If so, look through the 27-1 entries in jnnz for a non -1.
              DO k = 1, 26
                 IF (jnnz(k) /= -1) THEN
                    IF (jnnz(k) == innz(j)) THEN
                       ! ... If the local node number found in innz already recorded
                       ! ... in jnnz, then skip it.
                       EXIT
                    ENDIF
                 ELSE
                    ! ... If the local node number is not already recorded in jnnz
                    ! ... then enter it into jnnz.
                    jnnz(k) = innz(j)
                    EXIT
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
        
        ! ... Loop through rows with same nonzero blocks
       
        ! ... To construct cval, find the column value of a global node that
        ! ... is equal to the node stored in jnnz.
        DO n = 1, GNumNp
           DO k = 1, 27
              IF(jnnz(k) /= -1) THEN
                 ! ... Check to see that global node n is the same node as local node jnnz(k)=m.
                 ! ... If so, then enter the node column number in the next entry of cval.
                 m = Local2Global(jnnz(k))
                 IF (m == n) THEN
                    counter = counter + 1
                    ! ... using C convention, column n is n - 1
                    cval(counter) = m - 1
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
           
        ! ... Construct rp, this value corresponding to the index number
        ! ... counter + 1 in C convention.
        rp(i+1) = counter
        
        
        ! ... If this node is not on this processor.  My guess, this is just a place holder ... COstoich.
     ELSE
        
        rp(i+1) = rp(i)
        
     ENDIF
     
     
     
  ENDDO
  
  
END SUBROUTINE KtInitCRS
