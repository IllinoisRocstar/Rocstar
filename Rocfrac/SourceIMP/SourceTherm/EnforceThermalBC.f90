
!!****
!!
!!  NAME
!!     EnforceThermalBC
!!
!!  FUNCTION
!!     Enforces imposed temperature boundary conditions.
!!     Code may be added for time-dependant boundary conditions.
!!
!!  INPUTS
!!     NumNp -- Total number of nodes that this proc knows about
!!     LocNumNp -- Total number of nodes assigned to this processor
!!     disp -- Local temperature vector
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


SUBROUTINE EnforceThermalBC(NumNp,LocNumNp,Temp,node_flag,boundary_value,t,myid)

USE Precision
USE implicit_global

! Input variables
INTEGER :: NumNp
REAL(kind=wp) :: t
REAL(kind=wp),DIMENSION(1:LocNumNp) :: Temp
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
            counter = counter + 1
            IF(node_flag(i,1) == 8) THEN  ! Imposed constant nodal temperature
               Temp(counter) = boundary_value(i,1)
            ENDIF
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


END SUBROUTINE EnforceThermalBC
