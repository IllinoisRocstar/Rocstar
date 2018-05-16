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
     SUBROUTINE WRITE_OUTPUT_2
       use data_declarations

       IMPLICIT NONE

       CHARACTER*4 :: ichr4

       INTEGER :: ii

       WRITE(ichr4,'(i4.4)') i

       OPEN(35, FILE = PREFIX(1:LENGTH(PREFIX))//'/'//PREFIX(1:LENGTH(PREFIX))//'.'//ichr4//'.inp', &
            STATUS = 'UNKNOWN')

       ! Writing first record:  <# of elements> <# of nodes> <Number of sister nodes>
!       WRITE(35,*) numnp, numel
!
!       WRITE(35,*) delta
!       WRITE(35,*) numprop
!
!       DO l = 1, numprop
!         WRITE(35,*) ts_proportion(l),proportion(l)
!       ENDDO
!
       write(35,*) 1
       write(35,*) 2.5


       write(35,*) 2
       WRITE(35,*) numnp,0,0,0,0

       ! <Node ID> < X > < Y > < Z >
       DO l = 1,numnp
          WRITE(35,*) l , node(l)%coord(1:3),0  ! , 0 ! no corners
       ENDDO

       
100    FORMAT(I4,1X,$)  ! for the file ID write statement (used only once per list)
110    FORMAT(I8,1X,$)  ! for the nodelist (used numscale_np/2 times)
       ! If this is the first file to be written, Then only write information for one other file
       IF (count .EQ. 1) THEN
          write(35,*) 3
          WRITE(35,*) numboundfirst,0
          frontfile = 1

          DO l = 1, numbound
            IF((id(l)%NdBCflag .GE. 10) .AND. (id(l)%NdBCflag .LT. 100)) THEN
               WRITE(35,*) id(l)%ID, id(l)%NdBCflag/10,0
            ELSE IF (id(l)%NdBCflag .LT. 10) THEN
               WRITE(35,*) id(l)%ID, id(l)%NdBCflag,0
            END IF
          END DO

          write(35,*) 4
          WRITE(35,*) numbmeshfirst,0
          WRITE(6,*) numbmeshfirst

          DO l = 1, numboundmesh
            IF((idmesh(l)%NdBCflag .GE. 10) .AND. (idmesh(l)%NdBCflag .LT. 100)) THEN
               WRITE(35,*) idmesh(l)%ID, idmesh(l)%NdBCflag/10,0
            ELSE IF (idmesh(l)%NdBCflag .LT. 10) THEN
               WRITE(35,*) idmesh(l)%ID, idmesh(l)%NdBCflag,0
            END IF
          END DO

       ! Else if writing the last file, Then only write the current and last file information
       ELSE IF (count .EQ. scale) THEN
          backfile = scale - 2

          write(35,*) 3
          WRITE(35,*) numboundend,0

          DO l = 1, numbound
            IF (id(l)%NdBCflag .LT. 10) THEN
               WRITE(35,*) id(l)%ID, id(l)%NdBCflag,0
            ELSE IF (id(l)%NdBCflag .GE. 100) THEN
               WRITE(35,*) id(l)%ID, id(l)%NdBCflag/100,0
            END IF
          END DO

          ! Next Record
          write(35,*) 4
          WRITE(35,*) numbmeshend,0
          
         DO l = 1, numboundmesh
            IF(idmesh(l)%NdBCflag .LT. 10) THEN
               WRITE(35,*) idmesh(l)%ID, idmesh(l)%NdBCflag,0
            ELSE IF (idmesh(l)%NdBCflag .GE. 100) THEN
               WRITE(35,*) idmesh(l)%ID, idmesh(l)%NdBCflag/100,0
            END IF
         END DO

       ! Else if writing a middle file, Then write both back and front file information
       ELSE
          backfile = i - 1
          frontfile = i + 1
          write(35,*) 3
          WRITE(35,*) numboundmid,0
          
          DO l = 1, numbound
             IF (id(l)%NdBCflag .LT. 10) THEN
                WRITE(35,*) id(l)%ID, id(l)%NdBCflag,0
             END IF
          END DO
          
          ! Next Record
          write(35,*) 4
          WRITE(35,*) numbmeshmid
          
          DO l = 1, numboundmesh
             IF(idmesh(l)%NdBCflag .LT. 10) THEN
                WRITE(35,*) idmesh(l)%ID, idmesh(l)%NdBCflag,0
             END IF
          END DO
          
       END IF
       
       write(35,*) 5
       WRITE(35,*) numcstet, NBoundryEl3D,numcstet, NBoundryEl3D,4,0
       
       elemlist => firstelem
       eboundlist => firstebound

! fix:Scot -- Added a one for where the material id goes 
!  instead of the element id number

       ii = 0

       DO WHILE(associated(eboundlist%next))
          ii = ii + 1
          WRITE(35,*) 1, eboundlist%conn(1:4),0,0
          eboundlist => eboundlist%next
       ENDDO
       PRINT*,'Number of Boundary Elements =', ii

       DO WHILE(associated(elemlist%next))
          ii = ii + 1
          WRITE(35,*) 1, elemlist%conn(1:4),0,0
          elemlist => elemlist%next
       ENDDO
       PRINT*,'Total Number =',ii,numcstet
       
      ! These are all zero because cohesive elements are not yet supported
!       iaux = 0
!       DO l = 0, scale - 1
!          WRITE(35,*) iaux
!       ENDDO
!       WRITE(35,*) iaux

       write(35,*) 6

       j = i ! id processor


       IF (j .EQ. 0) THEN
          WRITE(35,*) 1            ! The number of files communicating with
!          WRITE(35,*) 1            ! The fileID the current proc. is communicating with
       ELSE IF (j .EQ. scale - 1) THEN
          WRITE(35,*) 1            ! The number of files communicating with
!          WRITE(35,*) scale - 2    ! The fileID the current proc. is communicating with
       ELSE
          WRITE(35,*) 2            ! The number of files communicating with
!          WRITE(35,*) j - 1, j + 1 ! The fileID's the current proc. is communicating with
       END IF

       DO k = 0, scale-1
          IF (k .EQ. j+1) THEN
             IF(j.EQ.0)THEN
                WRITE(35,*) 1,FileID(j)%num_sister_np
             ELSE IF(j.EQ.scale-1)THEN
                WRITE(35,*) scale - 2,FileID(j)%num_sister_np
             ELSE
                WRITE(35,*) j+1,FileID(j)%num_sister_np/2
             ENDIF
             backnode => first_back_node
             DO WHILE(associated(backnode%next))
                WRITE(35,*) backnode%sister%ID
                backnode => backnode%next
             END DO
       
          ELSE IF(k.EQ.j-1)THEN
             IF(j.EQ.0)THEN
                WRITE(35,*) 1,FileID(j)%num_sister_np
             elseif(j.EQ.scale-1)THEN
                WRITE(35,*) scale - 2,FileID(j)%num_sister_np
             ELSE
                WRITE(35,*) j-1,FileID(j)%num_sister_np/2
             ENDIF
             backnode => first_back_node
             DO WHILE(associated(backnode%next))
                WRITE(35,*) backnode%ID
                backnode => backnode%next
             END DO
!          ELSE
!             WRITE(35,*) 0
          ENDIF

       END DO
       
       ! Writes out neighboring file information

!!$       iaux = 0
!!$       DO l = 0, scale - 1
!!$          WRITE(35,*) iaux
!!$       ENDDO
!!$       WRITE(35,*) iaux

       write(35,*) 99
         
    CLOSE(35)
         
    END SUBROUTINE WRITE_OUTPUT_2






