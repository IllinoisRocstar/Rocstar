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
! ******************************************************************************
!
! Purpose: Set restart time and flag.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!
! Output: None.
!
! Notes: 
!   1. Cannot be set in routine which reads restart information file because 
!      in coupled runs that file is not called.
!
! ******************************************************************************
!
! $Id: RFLU_SetRestartTimeFlag.F90,v 1.4 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetRestartTimeFlag(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModMPI
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetRestartTimeFlag.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'RFLU_SetRestartTimeFlag',&
  'RFLU_SetRestartTimeFlag.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Setting restart time and flag...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set restart time and flag
! ******************************************************************************

  IF ( global%flowType == FLOW_STEADY ) THEN ! steady flow
    global%restartIter = global%currentIter

    IF ( global%restartIter == 0 ) THEN 
      global%restartFromScratch = .TRUE.
    ELSE 
      global%restartFromScratch = .FALSE.
    END IF ! global%restartIter
  ELSE ! unsteady flow
    global%restartTime = global%currentTime

    IF ( global%restartTime == 0.0_RFREAL ) THEN 
      global%restartFromScratch = .TRUE.
    ELSE 
      global%restartFromScratch = .FALSE.
    END IF ! global%restartIter     
  END IF ! global%flowType

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Setting restart time and flag done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetRestartTimeFlag

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetRestartTimeFlag.F90,v $
! Revision 1.4  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2004/06/25 20:04:01  haselbac
! Initial revision
!
! ******************************************************************************







