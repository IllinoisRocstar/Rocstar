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
! Purpose: open file for thrust history.
!
! Description: none.
!
! Input: global = case name, steady/unsteady flow.
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_OpenThrustFile.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_OpenThrustFile( global )

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

  INTEGER :: errorFlag

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_OpenThrustFile',&
  'RFLO_OpenThrustFile.F90' )

! open file

  IF (global%myProcid == MASTERPROC) THEN
    fname = TRIM(global%outDir)//TRIM(global%casename)//'.thr'

! - append to existing file (restart) or create new file

    IF ((global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL).OR.&
        (global%flowType==FLOW_STEADY   .AND. global%currentIter>1)) THEN
      OPEN(IF_THRUST,file=fname,form='formatted',status='old', &
                     position='append',iostat=errorFlag)
    ELSE
      OPEN(IF_THRUST,file=fname,form='formatted',status='unknown', &
                     iostat=errorFlag)
    ENDIF
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! - write header

    IF (global%thrustPlane == XCOORD) &
      WRITE(IF_THRUST,1000,iostat=errorFlag) 'x',global%thrustCoord
    IF (global%thrustPlane == YCOORD) &
      WRITE(IF_THRUST,1000,iostat=errorFlag) 'y',global%thrustCoord
    IF (global%thrustPlane == ZCOORD) &
      WRITE(IF_THRUST,1000,iostat=errorFlag) 'z',global%thrustCoord
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )
  ENDIF

! finalize

  CALL DeregisterFunction( global )

! format

1000 FORMAT('# thrust history (iteration/time, momentum, ', &
            ' pressure, total thrust [N])',/,'# ',A,' = ',1PE13.5)

END SUBROUTINE RFLO_OpenThrustFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_OpenThrustFile.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.1  2003/06/02 17:12:01  jblazek
! Added computation of thrust.
!
!******************************************************************************







