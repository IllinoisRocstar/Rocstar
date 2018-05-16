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
! Purpose: Set values derived from user input.
!
! Description: None.
!
! Input:
!   regions        Data for regions
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_DerivedInputValues.F90,v 1.4 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_DerivedInputValues(regions)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE PLAG_ModParameters

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg,nCont
#ifdef RFLO
  INTEGER :: iLev
#endif
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_DerivedInputValues.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'PLAG_DerivedInputValues',&
  'PLAG_DerivedInputValues.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
      'Setting derived input variables for particles...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set derived input values
! ******************************************************************************

#ifdef RFLO
  DO iReg = 1,global%nRegions

    regions(iReg)%mixtInput%computeTv = .TRUE.
    nCont = regions(iReg)%plagInput%nCont

    DO iLev = 1,regions(iReg)%nGridLevels

      IF (regions(iReg)%levels(iLev)%mixt%nTv < 2) THEN
        regions(iReg)%levels(iLev)%mixt%nTv = 2
      END IF ! nTv

      regions(iReg)%levels(iLev)%plag%nCv  = CV_PLAG_LAST + nCont
      regions(iReg)%levels(iLev)%plag%nDv  = DV_PLAG_LAST + nCont
      regions(iReg)%levels(iLev)%plag%nTv  = TV_PLAG_LAST
      regions(iReg)%levels(iLev)%plag%nAiv = AIV_PLAG_LAST
      regions(iReg)%levels(iLev)%plag%nArv = ARV_PLAG_LAST
      regions(iReg)%levels(iLev)%plag%nEv  = EV_PLAG_LAST + nCont
    END DO ! iLev
  END DO ! iReg
#endif

#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)

    regions(iReg)%mixtInput%computeTv = .TRUE.
    nCont = regions(iReg)%plagInput%nCont

    IF (regions(iReg)%mixtInput%nTv < 2) THEN
      regions(iReg)%mixtInput%nTv = 2
    END IF ! nTv

    regions(iReg)%plag%nCv  = CV_PLAG_LAST + nCont
    regions(iReg)%plag%nDv  = DV_PLAG_LAST + nCont
    regions(iReg)%plag%nTv  = TV_PLAG_LAST
    regions(iReg)%plag%nAiv = AIV_PLAG_LAST
    regions(iReg)%plag%nArv = ARV_PLAG_LAST
    regions(iReg)%plag%nEv  = EV_PLAG_LAST + nCont
  END DO ! iReg
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Setting derived input variables for particles done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DerivedInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_DerivedInputValues.F90,v $
! Revision 1.4  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/01/08 20:39:54  fnajjar
! Added nEv definition for Eulerian-based field
!
! Revision 1.1  2004/12/01 20:57:30  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/02/27 16:01:57  haselbac
! Fixed bug for RFLO
!
! Revision 1.1  2004/02/26 21:00:53  haselbac
! Initial revision
!
!******************************************************************************







