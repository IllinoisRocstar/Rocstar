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

!!****
!!
!!  NAME
!!     implicit_bc_enforce
!!
!!  FUNCTION
!!     Enforces imposed displacement, velocity, and acceleration boundary
!!     conditions.  Code may be added for time-dependant boundary conditions.
!!
!!  INPUTS
!!     NumNp -- Total number of nodes that this proc knows about
!!     LocNumNp -- Total number of nodes assigned to this processor
!!     disp -- Local displacment vector
!!     v -- Local velocity vector
!!     a -- Local acceleration vector
!!     node_flag -- Flags for each dof at each node as to what kind of BC is imposed
!!     boundary_value -- The magnitudes of the imposed boundary conditions
!!     t -- Simulation time.  Used for time-dependant boundary conditions.
!!     myid -- The rank of this processor.  Mainly used for debugging purposes.
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****


SUBROUTINE implicit_bc_enforce(NumNp,LocNumNp,disp,v,a,node_flag,boundary_value,t,myid)

USE Precision
USE implicit_global

! Input variables
INTEGER :: NumNp
REAL(kind=wp) :: t
REAL(kind=wp),DIMENSION(1:3*LocNumNp) :: disp, v, a
INTEGER,DIMENSION(1:NumNp,1:3) :: node_flag
REAL(kind=wp),DIMENSION(1:NumNp,1:3) :: boundary_value
INTEGER :: myid

! Internal variables
INTEGER :: i,j,counter,m

! Impose displacement boundary conditions
counter = 0
DO m = 1, GNumNp
   DO i = 1, NumNp
      IF (Local2Global(i) == m) THEN
         IF (NodeProc(i)==myid) THEN
            DO j = 1, 3
               counter = counter + 1
               IF(node_flag(i,j) == 8) THEN  ! Imposed constant nodal displacement
                  disp(counter) = boundary_value(i,j)
                  v(counter) = 0.0
                  a(counter) = 0.0
               ENDIF
            ENDDO
         ENDIF
      ENDIF
   ENDDO
ENDDO


!
! Add code here for special conditions such as time dependant BC's
!

!!$print*,'NOTE:  APPLYING TIME-DEPENDANT BOUNDARY CONDITIONS'

!!$IF (t < 5) THEN
!!$   DO i = 1, 4
!!$      disp(3*(i-1)+3) = 0.0
!!$      v(3*(i-1)+3) = 0.0
!!$      a(3*(i-1)+3) = 0.0
!!$   ENDDO
!!$ELSEIF (t < 10) THEN
!!$   DO i = 1, 4
!!$      disp(3*(i-1)+3) = 0.1 * (t-5)
!!$      v(3*(i-1)+3) = 0.1
!!$      a(3*(i-1)+3) = 0.0
!!$   ENDDO
!!$ELSE
!!$   DO i = 1, 4
!!$      disp(3*(i-1)+3) = 0.5
!!$      v(3*(i-1)+3) = 0.0
!!$      a(3*(i-1)+3) = 0.0
!!$   ENDDO
!!$ENDIF

!!$IF (t < 5) THEN
!!$   DO i = 1, 4
!!$      disp(3*(i-1)+3) = 0.0
!!$      v(3*(i-1)+3) = 0.0
!!$      a(3*(i-1)+3) = 0.0
!!$   ENDDO
!!$ELSEIF (t < 10) THEN
!!$   DO i = 1, 4
!!$      disp(3*(i-1)+3) = 0.1 * (t-5) * (t-5)
!!$      v(3*(i-1)+3) = 0.1 * (t-5)
!!$      a(3*(i-1)+3) = 0.1
!!$   ENDDO
!!$ELSE
!!$   DO i = 1, 4
!!$      disp(3*(i-1)+3) = 0.5 * (t-5)
!!$      v(3*(i-1)+3) = 0.5
!!$      a(3*(i-1)+3) = 0.0
!!$   ENDDO
!!$ENDIF


END SUBROUTINE implicit_bc_enforce

