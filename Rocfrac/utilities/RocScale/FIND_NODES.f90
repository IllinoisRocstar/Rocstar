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
     SUBROUTINE FIND_NODES
       use data_declarations

       IMPLICIT NONE

       INTEGER :: num_matches    ! The number of nodes with matching XY coord values
                                 ! This should equal numscale_np and is used for verification
       INTEGER :: iteration      ! The number of iteration the subroutine has gone through
       INTEGER :: max_iterations ! The maximum number of iterations for the search algorithm
       INTEGER :: step           ! The step value for redifining the tolerances

       REAL*8 :: low_tol       ! The lower tolerance by which 'FIND_NODES' determines a match
       REAL*8 :: hi_tol        ! The upper tolerance by which 'FIND_NODES' determines a match

       ! Variable initializations
       low_tol = 99999.00000
       hi_tol = 100001.00000
       max_iterations = 5
       iteration = 1
       step = 5

       IF (numscale_ft .NE. numscale_bk) THEN
          WRITE(6,*) "ERROR:  The number of nodes on the front and back surfaces do not equal."
          WRITE(6,*) "Please check your mesh and try again."
          STOP
       ENDIF

       WRITE(6,'(//A)') "Beginning search on sister nodes... This could take a while"

100    frontnode => first_front_node
       backnode => first_back_node
       num_matches = 0
       iteration = iteration + 1

       DO WHILE(associated(backnode%next))
          DO WHILE(associated(frontnode%next))

!             WRITE(6,*) "COORDS"
!             WRITE(6,*) "backnode", backnode%ID, backnode%coord(1:3)
!             WRITE(6,*) "frontnode", frontnode%ID, frontnode%coord(1:3)

!             WRITE(6,*) "DIVISION"
!             WRITE(6,*) backnode%coord(1)/frontnode%coord(1)
!             WRITE(6,*) backnode%coord(2)/frontnode%coord(2)

             ! If the 'X' coord is almost identical for both nodes THEN
             ! check the 'Y' coord value
             IF ((backnode%coord(1) .EQ. frontnode%coord(1)) .OR. &
                  ((backnode%coord(1)/frontnode%coord(1) .GT. low_tol/100000) .AND. &
                  (backnode%coord(1)/frontnode%coord(1) .LT. hi_tol/100000))) THEN
                ! IF the 'Y' coord is almost identical for both nodes THEN
                ! these nodes are sister nodes (a match is found)
                IF((backnode%coord(2) .EQ. frontnode%coord(2)) .OR. &
                     ((backnode%coord(2)/frontnode%coord(2) .GT. low_tol/100000) .AND. &
                     (backnode%coord(2)/frontnode%coord(2) .LT. hi_tol/100000))) THEN
                   backnode%sister => frontnode
                   frontnode%sister => backnode
                   num_matches = num_matches + 1
                ENDIF
             ENDIF
             ! proceed to the next node
             frontnode => frontnode%next
          ENDDO
          ! Once all the front nodes are cycled through, proceed to the next backnode
          ! and begin the search again
          backnode => backnode%next
          frontnode => first_front_node
       ENDDO



!!$       WRITE(6,*) "Writing node match information to '3Dnodematch.dat' file..."
!!$
!!$       OPEN(UNIT = 400, FILE = '3Dnodematch.dat', STATUS = 'UNKNOWN')
!!$
!!$       frontnode => first_front_node
!!$       backnode => first_back_node
!!$
!!$       WRITE(400,'(A8,3(6X,A1,5X),5X,A9,3(6X,A1,5X))') &
!!$            "BACKNODE", "X", "Y", "Z", "FRONTNODE", "X", "Y", "Z"
!!$
!!$       DO WHILE(associated(frontnode%next) .AND. associated(backnode%next))
!!$          WRITE(400,'(I8,3F12.8,6X,I8,3F12.8)') backnode%ID, backnode%coord(1:3), &
!!$                                      backnode%sister%ID, backnode%sister%coord(1:3)
!!$          frontnode => frontnode%next
!!$          backnode => backnode%next
!!$       END DO
!!$
!!$       CLOSE(400)



       WRITE(6,*) "Number of verified matches: ", num_matches
       WRITE(6,*) "Number of intended matches: ", numscale_np/2

       IF (iteration .LE. max_iterations) THEN
          IF (numscale_np/2 .NE. num_matches) THEN
             WRITE(6,*) "There was a failure in matching up the desired number of nodes."
             WRITE(6,*) "Redefining tolerances and starting iteration number ", iteration
             low_tol = low_tol - step
             hi_tol = hi_tol + step
             GOTO 100
             ELSE IF (numscale_np/2 .EQ. num_matches) THEN
                WRITE(6,*) "Number of matches confirmed.  Proceeding to next step."
          ENDIF
       ENDIF

     END SUBROUTINE FIND_NODES









