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
! Purpose: Determine whether to print data to screen.
!
! Description: None.
!
! Input:
!   global                      Pointer to global data
!
! Output: 
!   RFLU_DecidePrint = .TRUE.   If should print to screen
!   RFLU_DecidePrint = .FALSE.  If should not print to screen
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_DecidePrint.F90,v 1.6 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

LOGICAL FUNCTION RFLU_DecidePrint(global)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
      
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_global), POINTER :: global

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: logical1,logical2,logical3

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_DecidePrint.F90,v $ $Revision: 1.6 $'

  CALL RegisterFunction(global,'RFLU_DecidePrint',&
  'RFLU_DecidePrint.F90')

! *****************************************************************************
! Initialize
! *****************************************************************************

  RFLU_DecidePrint = .FALSE.
  
! *****************************************************************************
! Determine whether should print to screen
! *****************************************************************************

! =============================================================================
! Unsteady flow
! =============================================================================

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    logical1 = ABS(global%timeSincePrint-global%printTime) & 
             < 0.1_RFREAL*global%dtMin             
    logical2 = (global%timeSincePrint > global%printTime)    
    logical3 = (global%iterSinceRestart == 1)
  
    IF ( logical1 .OR. logical2 .OR. logical3 ) THEN
      RFLU_DecidePrint = .TRUE.
    END IF ! logical1 

! =============================================================================
! Steady flow. NOTE must have RFLU_DecidePrint be TRUE if currentIter == 1 
! because only then will RFLU_ResidualNorm be called from RFLU_TimeStepping, 
! and at the first iteration, this is needed to get a residual for correct 
! normalization.
! =============================================================================

  ELSE    
    RFLU_DecidePrint = (MOD(global%currentIter,global%printIter) == 0)

    IF ( global%currentIter == 1 ) THEN 
      RFLU_DecidePrint = .TRUE. 
    END IF ! global%currentIter
  END IF ! global%flowType

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_DecidePrint

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DecidePrint.F90,v $
! Revision 1.6  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.3  2004/01/31 03:56:54  haselbac
! Added explanatory comment
!
! Revision 1.2  2004/01/29 22:56:28  haselbac
! Removed setting of timeSincePrint to make routine more usable
!
! Revision 1.1  2003/10/29 21:36:58  haselbac
! Initial revision
!
!******************************************************************************







