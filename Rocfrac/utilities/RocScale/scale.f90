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
     PROGRAM SCALE_MESH
       use data_declarations

       IMPLICIT NONE


       WRITE(6,'(A24,$)') "Enter generic filename: "

       READ(5,*) PREFIX

       CALL READ_IO_FILES
       
       CALL FIND_NODES

       WRITE(6,'(A24,$)') "ENTER THE SCALE FACTOR: "
       READ(5,*) scale

       frontnode => first_front_node
       backnode => first_back_node

       ALLOCATE(FileID(0:scale-1))
       DO j = 0, scale-1
          IF (j .EQ. 0) THEN
             FileID(j)%num_sister_np = numscale_np/2
             FileID(j)%comm(1) = -1
             FileID(j)%comm(2) = 1
             FileID(j)%numboundnp = numbmeshfirst + numboundfirst
          ELSE IF (j .NE. (scale - 1)) THEN
             FileID(j)%num_sister_np = numscale_np
             FileID(j)%comm(1) = j-1
             FileID(j)%comm(2) = j+1
             FileID(j)%numboundnp = numbmeshmid + numboundmid
          ELSE IF (j .EQ. (scale - 1)) THEN
             FileID(j)%num_sister_np = numscale_np/2
             FileID(j)%comm(1) = j - 1
             FileID(j)%comm(2) = -1
             FileID(j)%numboundnp = numbmeshmid + numboundmid
          ENDIF
       ENDDO

       Ztrans = frontnode%coord(3) - backnode%coord(3)
       WRITE(6,*) Ztrans

       WRITE(6,*) "STARTING DUPLICATION OF MESH..."

       count = 1
       DO i = 0, scale-1
          OPEN(UNIT = 99, FILE = 'TEMP', STATUS = 'UNKNOWN')
          WRITE(99,*) i
          CLOSE(99)
          OPEN(UNIT = 99, FILE = 'TEMP', STATUS = 'UNKNOWN')
          READ(99,*) scaleID
          CLOSE(99)

          CALL WRITE_OUTPUT_2

          CALL TRANSLATE

          count = count + 1

       ENDDO

       PRINT*, '...Reading Surface Meshes'

       CALL READ_FRAC

       WRITE(6,*) "FINISHED WRITING ALL FILES"

     END PROGRAM SCALE_MESH




