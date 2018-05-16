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
     SUBROUTINE READ_FRAC
       use data_declarations

       IMPLICIT NONE

       INTEGER :: numsurfnp,       &  ! number of surface nodes
                  numsurfelem,     &  ! number of surface elements (3-node)
                  numsurfsis,      &  ! Number of intended sister node matches
                  numsurfsisb,     &  ! for the backface of the geometry
                  iteration,       &  ! Current iteration number
                  max_iterations,  &  ! Maximum number of iterations
                  scount,          &  ! The number of nodes communicating with
                  step                ! The step value for redefining the tolerances

       INTEGER :: num_matches         ! The current number of sister node matches
       REAL*8  :: low_tol             ! The low tolerance used to compare coord. values
       REAL*8  :: hi_tol              ! The high tolerance used to compare coord. values

       INTEGER :: Nfront_elem         ! The number of front boundary elements
       INTEGER :: Nback_elem          ! The number of back boundary elements
       INTEGER :: NBoundryEL          ! The total number of boundary elements

       ! Data type used to hold information on the boundary elements
       TYPE bound_elem
          INTEGER :: ID
          INTEGER :: conn(3)
          INTEGER :: flag
          TYPE(bound_elem), POINTER :: next
       END TYPE bound_elem


       INTEGER :: jj      ! some variable scot put in...  counter

       TYPE(bound_elem), POINTER ::        &
                         front_elem,       &  ! A list of the front boundary elements
                         back_elem,        &  ! A list of the back boundary elements
                         first_back_elem,  &  
                         first_front_elem

!       INTEGER, ALLOCATABLE, DIMENSION(:) :: NBoundryEl

!       ALLOCATE(NBoundryEl(1:scale))
!       NBoundryEl(:) = 0


       ! Variable initializations
       low_tol = 99999.00000
       hi_tol = 100001.00000
       max_iterations = 5
       iteration = 1
       step = 5

       ! Number of boundary element counters are initialized to '0'
       NBoundryEL = 0
       Nfront_elem = 0
       Nback_elem = 0

       ! Number of surface sister nodes initialized to '0'
       numsurfsis = 0
       numsurfsisb = 0

!---------------------------------------------------------------------------------------------------!
!                           READ 'frac.im' SURFACE MESH DESCRIPTION FILE                            !
!---------------------------------------------------------------------------------------------------!
       OPEN(UNIT = 40, FILE="fracSF.im", STATUS = "OLD")

       WRITE(6,*) 
       WRITE(6,*) 
       WRITE(6,'(A,$)') "Reading 'fracSF.im' file... "

       READ(40,*)  ! INFORMATION NOT NEEDED
       READ(40,*)  ! INFORMATION NOT NEEDED

       ! Reading in number of surface node points and surface elements
       READ(40,*) numsurfnp, numsurfelem

       ! Allocating surface node array
       ALLOCATE(surfnode(numsurfnp))

        ! Setting up L-lists for front surface and back surface nodes
       ALLOCATE(frontsnode); first_frontsnode => frontsnode; NULLIFY(frontsnode%next)
       ALLOCATE(backsnode); first_backsnode => backsnode; NULLIFY(backsnode%next)      

       ! Loops through the number of surface nodes in the file
       DO l = 1, numsurfnp

          ! Reading the surface node coordinates and marker
          READ(40,*) surfnode(l)%coord(1:3), surfnode(l)%marker

          ! IF this node lies on the back edge...
          IF (surfnode(l)%coord(3) .LE. minZ) THEN

             ! ...THEN add this node to the back surface node L-list
             backsnode%ID = l
             backsnode%coord(1:3) = surfnode(l)%coord(1:3)
             backsnode%marker = surfnode(l)%marker
             ALLOCATE(backsnode%next); backsnode => backsnode%next; NULLIFY(backsnode%next)

             ! Updates the counter for number of surface sister nodes (edge nodes)
             numsurfsis = numsurfsis + 1

          ! ELSE IF this node is located on the front edge...
          ELSE IF (surfnode(l)%coord(3) .GE. maxZ) THEN

             ! ...THEN add this node to the front surface node L-list
             frontsnode%ID = l
             frontsnode%coord(1:3) = surfnode(l)%coord(1:3)
             frontsnode%marker = surfnode(l)%marker
             ALLOCATE(frontsnode%next); frontsnode => frontsnode%next; NULLIFY(frontsnode%next)

             ! Updates the counter for number of surface sister nodes (edge nodes)
             numsurfsisb = numsurfsisb + 1   ! These are BACK edge nodes
          END IF
          
       END DO

       elemlist => firstelem       ! Lists all of the surface elements

       ! Allocation of the L-lists that contain element info only for front and back boundary elements
       ALLOCATE(front_elem); first_front_elem => front_elem; NULLIFY(front_elem%next)
       ALLOCATE(back_elem); first_back_elem => back_elem; NULLIFY(back_elem%next)


       ALLOCATE(surface(numsurfelem))

       ! Loop that reads in element data and searches for boundary elements
       DO l = 1, numsurfelem

          ! Reads in the element data from the 'frac.im' file
          READ(40,*) surface(l)%conn(1:3), surface(l)%marker, m, n

          ! Loops through all nodes.  (used in searching for boundary elements)
          DO m = 1, numsurfnp

             ! If there is a node on the front boundary... (checks node z-coordinate)
             IF(surfnode(m)%coord(3) .GE. maxZ) THEN
                
                ! ...THEN check to see if the current element has this node in its conn array
                DO n = 1, 3

                   ! If the element conn array contains the node in question...
                   IF (m .EQ. surface(l)%conn(n)) THEN
 
                      ! ...THEN add this element and its information to 'front_elem' L-list
                      ! This list keeps track of only the front elements
                      front_elem%ID = l
                      front_elem%conn(1:3) = surface(l)%conn(1:3)
                      front_elem%flag = surface(l)%marker
                      ALLOCATE(front_elem%next); front_elem => front_elem%next; NULLIFY(front_elem%next)

                      ! These two variables keep count of the elements
                      Nfront_elem = Nfront_elem + 1    ! Keeps count of only the front boundary elements
                      NBoundryEL = NBoundryEL + 1      ! Keeps count of all boundary elements

                      ! Since the element was found to be on the boundary...
                      GOTO 50  ! ...exit these loops and move on to the next element
                   END IF
                END DO
             
             ! ELSE IF the node is on the back boundary... (checks node z-coordinate)
             ELSE IF (surfnode(m)%coord(3) .LE. minZ) THEN

                ! ...THEN check to see if the current element has this node in its conn array
                DO n = 1, 3

                   ! IF the element conn array contains the node in question...
                   IF (m .EQ. surface(l)%conn(n)) THEN

                      ! THEN add this element and its information to 'back_elem' L-list
                      ! This list keeps track of only the back elements
                      back_elem%ID = l
                      back_elem%conn(1:3) = surface(l)%conn(1:3)
                      back_elem%flag = surface(l)%marker
                      ALLOCATE(back_elem%next); back_elem => back_elem%next; NULLIFY(back_elem%next)

                      ! These two variables keep count of the elements
                      Nback_elem = Nback_elem + 1    ! Keeps count of only the back boundary elements 
                      NBoundryEL = NBoundryEL + 1    ! Keeps count of all boundary elements
                      
                      ! Since the element was found to be on the boundary...
                      GOTO 50  ! ...exit these loops and move on to the next element
                   END IF
                END DO
             END IF
          END DO

          ! Now that it has been determined - if the element lies on the boundary...
          ! ...we can add its information to the main element L-list
          elemlist%ID = l
          elemlist%conn(1:3) = surface(l)%conn(1:3)
          elemlist%conn(4) = 0  ! since this list was used for vol. Elem... we assign '0'
                                ! to the unneeded cell in the conn array
          elemlist%flag = surface(l)%marker
         ALLOCATE(elemlist%next); elemlist => elemlist%next; NULLIFY(elemlist%next)

50     END DO

      ! We are now finished reading the 'frac.im' file
      CLOSE(40)


!---------------------------------------------------------------------------------------------------!
!                            RUNNING SEARCH ON SISTER NODE MATCHES                                  !
!---------------------------------------------------------------------------------------------------!

      iteration = 0   ! initializes the number of iterations needed to find sister nodes

       WRITE(6,*) "DONE."
       WRITE(6,*) "Number of back edge nodes: ", numsurfsisb
       WRITE(6,*) "Number of front edge nodes: ", numsurfsis

       ! Just checks to see if there were the same number of nodes on both the front and back edges
       IF (numsurfsisb .EQ. numsurfsis) THEN
          WRITE(6,*) "...node quantities confirmed."
       ELSE
          ! If there wasn't a match, the program exits because nothing else can be done
          WRITE(6,*) "There was an error in calculating the number of matching edge nodes."
          STOP
       END IF

       scount = numsurfsis
100    frontsnode => first_frontsnode
       backsnode => first_backsnode
       iteration = iteration + 1
       num_matches = 0

       WRITE(6,*) "Beginning search on sister nodes..."

       ! Loop through all back edge nodes
       DO WHILE(associated(backsnode%next))

          ! Also loop through all front edge nodes
          DO WHILE(associated(frontsnode%next))

             ! If the 'X' coord is almost identical for both nodes THEN
             ! check the 'Y' coord value
             IF ((backsnode%coord(1) .EQ. frontsnode%coord(1)) .OR. &
                  ((backsnode%coord(1)/frontsnode%coord(1) .GT. low_tol/100000) .AND. &
                  (backsnode%coord(1)/frontsnode%coord(1) .LT. hi_tol/100000))) THEN

                ! IF the 'Y' coord is almost identical for both nodes THEN
                ! these nodes are sister nodes (a match is found)
                IF((backsnode%coord(2) .EQ. frontsnode%coord(2)) .OR. &
                     ((backsnode%coord(2)/frontsnode%coord(2) .GT. low_tol/100000) .AND. &
                     (backsnode%coord(2)/frontsnode%coord(2) .LT. hi_tol/100000))) THEN

                   backsnode%sister => frontsnode
                   frontsnode%sister => backsnode

                   ! updates the number of sister node matches found
                   num_matches = num_matches + 1

                ENDIF
             ENDIF

             ! proceed to the next node
             frontsnode => frontsnode%next
          ENDDO

          ! Once all the front nodes are cycled through, proceed to the next backsnode
          ! and begin the search again
          backsnode => backsnode%next
          frontsnode => first_frontsnode
       ENDDO

       WRITE(6,*) "Number of verified matches: ", num_matches
       WRITE(6,*) "Number of intended matches: ", numsurfsis


       WRITE(6,*) "Writing node match information to 'nodematch.dat' file..."



       OPEN(UNIT = 500, FILE = '2Dnodematch.dat', STATUS = 'UNKNOWN')

       WRITE(500,'(A8,3(6X,A1,5X),5X,A9,3(6X,A1,5X))') &
            "BACKNODE", "X", "Y", "Z", "FRONTNODE", "X", "Y", "Z"
       
       frontsnode => first_frontsnode
       backsnode => first_backsnode

       DO WHILE(associated(frontsnode%next) .AND. associated(backsnode%next))
          WRITE(500,'(I8,3F12.8,6X,I8,3F12.8)') backsnode%ID, backsnode%coord(1:3), &
                                      backsnode%sister%ID, backsnode%sister%coord(1:3)
          frontsnode => frontsnode%next
          backsnode => backsnode%next
       END DO

       CLOSE(500)



       IF (iteration .LE. max_iterations) THEN
          IF (numsurfsis .NE. num_matches) THEN
             WRITE(6,*) "There was a failure in matching up the desired number of nodes."
             WRITE(6,*) "Redefining tolerances and starting iteration number ", iteration
             low_tol = low_tol - step
             hi_tol = hi_tol + step
             GOTO 100
          ELSE IF (numscale_np/2 .EQ. num_matches) THEN
                WRITE(6,*) "Number of matches confirmed.  Proceeding to next step."
          ENDIF
       ENDIF


!---------------------------------------------------------------------------------------------------!
!                                    WRITING OUTPUT FILE                                            !
!---------------------------------------------------------------------------------------------------!

       OPEN(UNIT = 45, FILE = PREFIX(1:LENGTH(PREFIX))//'/fracSF.im', &
            STATUS = "UNKNOWN")
       WRITE(45,*) scale, 3

       ! Loops through all processors
       DO k = 1, scale

          WRITE(45,*) k  ! Writes the current Processor Number

          ! Writes the number of surface nodes and elements (This is the same for all Processors)
          WRITE(45,*) numsurfnp, numsurfelem

          ! Loops through all surface nodes
          DO l = 1, numsurfnp
             ! Writes the current node information to output file
             WRITE(45,*) surfnode(l)%coord(1:3), surfnode(l)%marker,0
          END DO

          ! Returns L-lists to the beginning
          elemlist => firstelem
          front_elem => first_front_elem
          back_elem => first_back_elem


          ! First it lists the back edge boundary elements...
          DO WHILE(associated(back_elem%next))
             WRITE(45,*) back_elem%conn(1:3), back_elem%flag, m, n
             back_elem => back_elem%next
             m = m + 1
          END DO

          ! ...Then it lists the front edge boundary elements...
          DO WHILE(associated(front_elem%next))
             WRITE(45,*) front_elem%conn(1:3), front_elem%flag, m, n
             front_elem => front_elem%next
             m = m + 1
          END DO

          ! Finally, we list all the internal nodes
          DO WHILE(associated(elemlist%next))
             WRITE(45,*) elemlist%conn(1:3), elemlist%flag, m, n
             elemlist => elemlist%next
          END DO

          ! Beginning translation
          DO l = 1, numsurfnp
             surfnode(l)%coord(3) = surfnode(l)%coord(3) + Ztrans
          END DO
       END DO

       CLOSE(45)

       OPEN(UNIT = 45, FILE = PREFIX(1:LENGTH(PREFIX))//'/fracS.im', &
            STATUS = "UNKNOWN")
       WRITE(45,*) scale, 3
       DO i = 1, scale
          WRITE(45,*) i ! processor id
          WRITE(45,*) 0, 0
       ENDDO

       CLOSE(45)


       ! Starting communication information
!!$       DO j = 0, scale - 1
!!$          
!!$          ! Writes out the current processor and total number of boundary nodes
!!$          IF (j .EQ. 0) THEN
!!$             WRITE(45,*) j, NBoundryEL  ! Nfront_elem for the number of front edge elements
!!$
!!$             ! Loop through each processor again
!!$             DO jj = 0, scale - 1
!!$
!!$                ! IF the current proc is the first one...
!!$                IF(jj .EQ. 0)THEN
!!$
!!$                   ! There is no communication information sent to itself
!!$                   WRITE(45,*) 0
!!$
!!$                ! ELSE IF the current proc is the very next one in the list...
!!$                ELSE IF(jj .EQ. 1)THEN
!!$
!!$                   ! Write the number of sister nodes communicating with
!!$                   WRITE(45,*) numsurfsis  ! Number of nodes sending
!!$                   backsnode => first_backsnode
!!$
!!$                   ! Then list all of these nodes
!!$                   DO WHILE(associated(backsnode%next))
!!$                      WRITE(45,*) backsnode%sister%ID
!!$                      backsnode => backsnode%next
!!$                   END DO
!!$
!!$                ! ELSE no other information is being sent to any other proc 
!!$                ELSE
!!$                   WRITE(45,*) 0 
!!$                ENDIF
!!$             ENDDO
!!$          
!!$
!!$          ! ELSE IF the current processor is somewhere in the middle...
!!$          ELSE IF ((j .GT. 0).AND.(j .LT. scale-1)) THEN  ! middle blocks
!!$             ! Write the number of boundary elements
!!$             WRITE(45,*) j, NBoundryEL
!!$
!!$             ! Loops through all processors again
!!$             DO jj = 0, scale - 1
!!$
!!$                ! IF the current processor is 1 less than the other looped proc #...
!!$                IF(jj .EQ. j+1)THEN
!!$                   
!!$                   WRITE(45,*) numsurfsis  ! Number of nodes sending
!!$                   backsnode => first_backsnode
!!$                   ! Then write the list of these nodes
!!$                   DO WHILE(associated(backsnode%next))
!!$                      WRITE(45,*) backsnode%sister%ID
!!$                      backsnode => backsnode%next
!!$                   END DO
!!$                   
!!$                ! IF the current processor is equal to the other looped proc #...
!!$                ELSE IF(jj .EQ. j)THEN
!!$                    WRITE(45,*) 0   ! Nothing is communicated with itself
!!$
!!$                ! ELSE IF the current proc is the other looped proc right after this one...
!!$                ELSE IF(jj .EQ. j-1)THEN
!!$
!!$                    WRITE(45,*) numsurfsisb  ! Number of nodes sending
!!$                    backsnode => first_backsnode
!!$                    ! Then write the list of these nodes
!!$                    DO WHILE(associated(backsnode%next))
!!$                       WRITE(45,*) backsnode%ID
!!$                       backsnode => backsnode%next
!!$                    END DO
!!$
!!$                 ! Else write '0' for the communication information
!!$                 ELSE
!!$                    WRITE(45,*) 0
!!$                 ENDIF
!!$              END DO
!!$
!!$          ! Else if the current processor is the last proc...
!!$          ELSE IF (j .EQ. scale-1) THEN
!!$             ! Write the number of boundary elements
!!$             WRITE(45,*) j, NBoundryEL  ! Nback_elem for just the back edge element total
!!$
!!$             ! Loops through all processors again             
!!$             DO jj = 0, scale - 1
!!$
!!$                ! IF the current processor is 1 less than the other looped proc #...
!!$                IF(jj .EQ. scale-2)THEN
!!$
!!$                   WRITE(45,*) scount  ! Number of nodes sending
!!$                   backsnode => first_backsnode
!!$                   ! Then write the list of these nodes
!!$                   DO WHILE(associated(backsnode%next))
!!$                      WRITE(45,*) backsnode%ID
!!$                      backsnode => backsnode%next
!!$                   END DO
!!$                ! Else write '0' for the communication information
!!$                ELSE
!!$                   WRITE(45,*) 0
!!$                ENDIF
!!$             ENDDO
!!$          END IF
!!$          
!!$          ! Writes out neighboring file information
!!$          IF (j .EQ. 0) THEN
!!$             WRITE(45,*) 1            ! The number of files communicating with
!!$             WRITE(45,*) 1            ! The fileID the current proc. is communicating with
!!$
!!$          ELSE IF (j .EQ. scale - 1) THEN
!!$             WRITE(45,*) 1            ! The number of files communicating with
!!$             WRITE(45,*) (scale - 2)    ! The fileID the current proc. is communicating with
!!$
!!$          ELSE
!!$             WRITE(45,*) 2            ! The number of files communicating with
!!$             WRITE(45,*) (j - 1), (j + 1) ! The fileID's the current proc. is communicating with
!!$          END IF
!!$       END DO
!!$       
     END SUBROUTINE READ_FRAC






