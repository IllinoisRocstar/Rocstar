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
     SUBROUTINE READ_FILES
       use data_declarations

       IMPLICIT NONE

       ! The following statements are because etype not user defined
       etype = 4   ! Element type (4-node tetrahedra)
       ftype = 3   ! Face type (3-node triangles)

!-----------------------------------------------------------------------------------------!
!                                1.0) .noboite READ SECTION                               !
!-----------------------------------------------------------------------------------------!

       OPEN(UNIT = 10, FILE = PREFIX(1:LENGTH(PREFIX))//'.noboite', STATUS = 'OLD', ERR = 1001)
       
write(6,*) "       ! Reading record #1 of .noboite file"
       READ(10,*) numel,    &  ! Number of elements
            numnp,          &  ! Number of nodes
            npfixed            ! Number of original nodes (from surface mesh)
       
       ALLOCATE(element(numel))
       ALLOCATE(node(numnp))

       DO i = 1, numel
          ALLOCATE(element(i)%conn(etype))
       ENDDO
       
       ! Reading record #2 of .noboite file (element connectivity array)
       WRITE(6,*) "Reading record #2 of .noboite file (element connectivity array)"
       READ(10,*) (element(i)%conn(1:etype), i = 1, numel)
       
       ! Reading record #3 of .noboite file (node coordinates)
       WRITE(6,*) "Reading record #3 of .noboite file (node coordinates)"
       READ(10,*) (node(i)%coord(1:3), i = 1, numnp)
       
       ! Reading the number of subregions designated by Tetmesh
       READ (10,*) numsub
       
       ! Reading the subregion triplets (faces located within the region)
       ! This information is not needed and is bypassed
       WRITE(6,*) "Reading the subregion triplets (faces located within the region)"
       READ (10,*)
       
       IF (numsub .EQ. 1) THEN
          DO i = 1, numel
             element(i)%marker = 1
          ENDDO
       ELSE IF (numsub .NE. 1) THEN
          ! Reading in the assigned subregion numbers given to each element
          WRITE(6,*) "Reading in the assigned subregion numbers given to each element"
          READ(10,*) (element(i)%marker, i = 1, numel)
          write(6,*) "finished reading subregion number information"
       ENDIF
       
       ! Closing .noboite file
       CLOSE(10)

!-----------------------------------------------------------------------------------------!
!                                1.2) .faces READ SECTION                                 !
!-----------------------------------------------------------------------------------------!

       OPEN(UNIT = 20, FILE = PREFIX(1:LENGTH(PREFIX))//'.faces', STATUS = 'OLD', ERR = 1002)

       READ(20,*) numfaces    ! Number of faces
       
       ! Initializes the face linked list
       ALLOCATE(facelist); NULLIFY(facelist%next)
       firstface => facelist; numbndfaces = 0
       ALLOCATE(facelist%conn(ftype))

       ! Reads through all surface faces and assigns corresponding connectivity and flags
       DO i = 1, numfaces
          READ(20,*) facelist%ID, facelist%conn(1:ftype), facelist%marker

          IF (facelist%marker .NE. 0) THEN
             ALLOCATE(facelist%next); facelist => facelist%next; NULLIFY(facelist%next)
             ALLOCATE(facelist%conn(ftype))
             numbndfaces = numbndfaces + 1
          ENDIF
       ENDDO

       ! Closing .faces file
       CLOSE(20)
       
!-----------------------------------------------------------------------------------------!
!                                1.3) .points READ SECTION                                !
!-----------------------------------------------------------------------------------------!

       OPEN(UNIT = 30, FILE = PREFIX(1:LENGTH(PREFIX))//'.points', STATUS = 'OLD', ERR = 1003)

       ALLOCATE(frontnode); NULLIFY(frontnode%next)
       first_front_node => frontnode; numscale_np = 0

       ALLOCATE(backnode); NULLIFY(backnode%next)
       first_back_node => backnode; numscale_np = 0
       
       ALLOCATE(BNPlist); NULLIFY(BNPlist%next)
       firstnode => BNPlist; numbndnp = 0

       numscale_bk = 0
       numscale_ft = 0
       ! Reading node marker information from .points file
       READ(30,*) ! npfixed... this information was already acquired in the .noboite file
       WRITE(6,*) "Reading node marker information from .points file"
       DO i = 1, npfixed
          READ(30,*) XYZ(1:3), mark
          
          ! If the marker is '1' then they are assigned to the backnode list
          IF (mark .EQ. 1) THEN
             backnode%ID = i
             backnode%coord(1:3) = XYZ(1:3)
             backnode%marker = mark
             ALLOCATE(backnode%next); backnode => backnode%next; NULLIFY(backnode%next)
             numscale_np = numscale_np + 1
             numscale_bk = numscale_bk + 1
          ! Else if the marker is '2' then they are assigned to the frontnode list
          ELSE IF (mark .EQ. 2) THEN
             frontnode%ID = i
             frontnode%coord(1:3) = XYZ(1:3)
             frontnode%marker = mark
             ALLOCATE(frontnode%next); frontnode => frontnode%next; NULLIFY(frontnode%next)
             numscale_np = numscale_np + 1
             numscale_ft = numscale_ft + 1
          ! Else if the marker is other than '0' then they are assigned to the boundary node list
          ELSE IF (mark .NE. 0) THEN
             numbndnp = numbndnp + 1
             BNPlist%ID = i
             BNPlist%marker = mark
             ALLOCATE(BNPlist%next); BNPlist => BNPlist%next; NULLIFY(BNPlist%next)
          END IF
          ! All nodes are added to the general list
          node(i)%coord(1:3) = XYZ(1:3)
          node(i)%marker = mark
       ENDDO

       ! Closing .points file
       CLOSE(30)

1000   GOTO 1010

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   ERROR STATEMENTS                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1001   WRITE(6,*) "AN ERROR OCCURRED WHILE OPENING ", PREFIX(1:LENGTH(PREFIX))//'.noboite'
       WRITE(6,*) "FILE DOES NOT EXIST."
       STOP
1002   WRITE(6,*) "AN ERROR OCCURRED WHILE OPENING ", PREFIX(1:LENGTH(PREFIX))//'.faces'
       WRITE(6,*) "FILE DOES NOT EXIST."
       STOP
1003   WRITE(6,*) "AN ERROR OCCURRED WHILE OPENING ", PREFIX(1:LENGTH(PREFIX))//'.points'
       WRITE(6,*) "FILE DOES NOT EXIST."
       STOP


1010   END SUBROUTINE READ_FILES

