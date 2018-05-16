! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
!******************************************************************************
!
! Purpose: write thrust history into a file.
!
! Description: none.
!
! Input: global%thrustTotal = total thrust
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: WriteThrust.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteThrust( global )

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(CHRLEN+4) :: fname

#ifdef MPI
  REAL(RFREAL) :: localThrust(2), globalThrust(2)
#endif

!******************************************************************************

  CALL RegisterFunction( global,'WriteThrust',&
  'WriteThrust.F90' )

! sum up data from other processors

#ifdef MPI
  localThrust(1) = global%thrustMom
  localThrust(2) = global%thrustPress

  CALL MPI_Reduce( localThrust,globalThrust,2,MPI_RFREAL,MPI_SUM,MASTERPROC, &
                   global%mpiComm,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

  global%thrustMom   = globalThrust(1)
  global%thrustPress = globalThrust(2)
#endif
  global%thrustTotal = global%thrustMom + global%thrustPress

! steady flow

  IF (global%flowType==FLOW_STEADY .AND. global%myProcid==MASTERPROC) THEN
    WRITE(IF_THRUST,1000,err=10) global%currentIter,global%thrustMom, &
                                 global%thrustPress,global%thrustTotal

! unsteady flow

  ELSE IF (global%flowType==FLOW_UNSTEADY .AND. global%myProcid==MASTERPROC) THEN
    WRITE(IF_THRUST,2000,err=10) global%currentTime,global%thrustMom, &
                                 global%thrustPress,global%thrustTotal
  ENDIF

! close and open file (instead of fflush)

  IF (global%thrustOpenClose .AND. global%myProcid==MASTERPROC) THEN
    WRITE(fname,'(A)') TRIM(global%outDir)//TRIM(global%casename)//'.thr'
    CLOSE(IF_THRUST)
    OPEN(IF_THRUST,FILE=fname,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
  ENDIF

! finalize

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,'Thrust history file.' )

1000 FORMAT(I6,1PE13.4,2E13.4)
2000 FORMAT(1PE12.5,3E13.4)

999  CONTINUE

END SUBROUTINE WriteThrust

!******************************************************************************
!
! RCS Revision history:
!
! $Log: WriteThrust.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:52:27  haselbac
! Initial revision after changing case
!
! Revision 1.1  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
!******************************************************************************







