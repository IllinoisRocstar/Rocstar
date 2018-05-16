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
     SUBROUTINE WRITE_OUTPUT
       use data_declarations

       IMPLICIT NONE

       ! If the current file number is less than 10, 
       ! then put a '.00' before the number in the file extension
       IF (i .LT. 10) THEN
          OPEN(35, FILE = PREFIX(1:LENGTH(PREFIX))//'.00'//scaleID(1:LENGTH(scaleID)), &
               STATUS = 'UNKNOWN')
       ! else if the # is greater than 10 but less than 100 put a '.0' before the # in the extension
       ELSE IF ((i .GE. 10).AND.(i .LT. 100)) THEN
          OPEN(35, FILE = PREFIX(1:LENGTH(PREFIX))//'.0'//scaleID(1:LENGTH(scaleID)), &
               STATUS = 'UNKNOWN')
       ! Else if the # is greater than or equal to 100, just add the '.' in front of the extension
       ELSE IF (i .GE. 100) THEN
          OPEN(35, FILE = PREFIX(1:LENGTH(PREFIX))//'.'//scaleID(1:LENGTH(scaleID)), &
               STATUS = 'UNKNOWN')
       END IF

       ! Writing first record:  <# of elements> <# of nodes> <Number of sister nodes>
       WRITE(35,*) numel, numnp, numscale_np

       ! Record #2:  <Node ID> < X > < Y > < Z >
       DO l = 1, numnp
          WRITE(35,'(I8,3F15.8)') l, node(l)%coord(1:3)
       ENDDO

       ! Record #3:  <# of Boundary nodes>
       WRITE(35,*) numbndnp

       ! Record #4:  <Boundary node ID> <marker>
       BNPlist => firstnode
       DO WHILE(associated(BNPlist%next))
          WRITE(35,*) BNPlist%ID, BNPlist%marker
          BNPlist => BNPlist%next
       END DO

       ! Record #5:  <# of Boundary faces (w/flags)>
       WRITE(35,*) numbndfaces

       ! Record #6:  <faceID> <face connectivity> <marker>
       facelist => firstface
       DO WHILE(associated(facelist%next))
          WRITE(35,*) facelist%ID, facelist%conn(1:3), facelist%marker
          facelist => facelist%next
       END DO

       ! Record #7:  <element ID> <element connectivity> <marker>
       DO l = 1, numel
          WRITE(35,*) l, element(l)%conn(1:4), element(l)%marker
       END DO

       frontnode => first_front_node
       backnode => first_back_node
       
100    FORMAT(I4,1X,$)  ! for the file ID write statement (used only once per list)
110    FORMAT(I8,1X,$)  ! for the nodelist (used numscale_np/2 times)
       ! If this is the first file to be written, Then only write information for one other file
       IF (count .EQ. 1) THEN
          frontfile = 2
          ! Record #8:  <current file #> <the mesh file next> <Number of invlolved sister nodes>
          WRITE(35,*) i, frontfile, numscale_np/2
          ! Record #9:  <current fileID> <1:numscale/2 nodes>
          WRITE(35,100) i
          DO WHILE(associated(frontnode%next))
             WRITE(35,110) frontnode%ID
             frontnode => frontnode%next
          END DO
          WRITE(35,*)  ! ends the current record being written

          frontnode => first_front_node
          ! Record #10:  <next fileID> <1:numscale/2 sister nodes>
          WRITE(35,100) frontfile
          DO WHILE(associated(frontnode%next))
             WRITE(35,110) frontnode%sister%ID
             frontnode => frontnode%next
          END DO
       
       ! Else if writing the last file, Then only write the current and last file information
       ELSE IF (count .EQ. scale) THEN
          backfile = scale - 1
          ! Record #8:  <current file #> <the previous fileID> <# of involved sister nodes>
          WRITE(35,*) i, backfile, numscale_np/2

          ! Record #9:  <current fileID> <1:numscale/2 nodes>
          WRITE(35,100) i
          DO WHILE(associated(backnode%next))
             WRITE(35,110) backnode%ID
             backnode => backnode%next
          END DO
          WRITE(35,*)  ! ends the current record being written

          backnode => first_back_node
          ! Record #10:  <next fileID> <1:numscale/2 sister nodes>
          WRITE(35,100) backfile
          DO WHILE(associated(backnode%next))
             WRITE(35,110) backnode%sister%ID
             backnode => backnode%next
          END DO

       ! Else if writing a middle file, Then write both back and front file information
       ELSE
         backfile = i - 1
         frontfile = i + 1
         ! Record #8:  <current file #> <previous fileID> <next fileID> <# of involved sister nodes>
         WRITE(35,*) i, backfile, frontfile, numscale_np/2

         backnode => first_back_node
         ! Record #9:  <current fileID> <1:numscale/2 nodes>
         WRITE(35,100) i
         DO WHILE(associated(backnode%next))
            WRITE(35,110) backnode%ID
            backnode => backnode%next
         END DO
         WRITE(35,*)  ! ends the current record being written
         
         backnode => first_back_node
         ! Record #10:  <next fileID> <1:numscale/2 sister nodes>
         WRITE(35,100) backfile
         DO WHILE(associated(backnode%next))
            WRITE(35,110) backnode%sister%ID
            backnode => backnode%next
         END DO
         WRITE(35,*)  ! ends the current record being written

         frontnode => first_front_node
         ! Record #11:  <current fileID> <1:numscale/2 nodes>
         WRITE(35,100) i
         DO WHILE(associated(frontnode%next))
            WRITE(35,110) frontnode%ID
            frontnode => frontnode%next
         END DO
         WRITE(35,*)  ! ends the current record being written

         frontnode => first_front_node
         ! Record #12:  <next fileID> <1:numscale/2 nodes>
         WRITE(35,100) frontfile
         DO WHILE(associated(frontnode%next))
            WRITE(35,110) frontnode%sister%ID
            frontnode => frontnode%next
         END DO
         
      END IF

      count = count + 1
      CLOSE(35)

    END SUBROUTINE WRITE_OUTPUT

