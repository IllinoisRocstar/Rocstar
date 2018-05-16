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
! Purpose: Determine whether to build geometry.
!
! Description: None.
!
! Input:
!   global                      Pointer to global data
!
! Output: 
!   RFLU_DecideBuildGeometry    TRUE or FALSE
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_DecideBuildGeometry.F90,v 1.7 2008/12/06 08:45:05 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_DecideBuildGeometry(global)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
      
  IMPLICIT NONE

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

  RCSIdentString = '$RCSfile: RFLU_DecideBuildGeometry.F90,v $ $Revision: 1.7 $'

! ******************************************************************************
! Initialize
! ******************************************************************************

  RFLU_DecideBuildGeometry = .FALSE.
  
! ******************************************************************************
! Determine whether should build geometry. NOTE need initFlowFlag here also 
! because this signals that errors can be computed...
! ******************************************************************************

  IF ( (global%initFlowFlag == INITFLOW_FROMHARDCODE) .OR. & 
       (global%postExtractFlag .EQV. .TRUE.) .OR. &
       (global%postInterpType == INTERP_TYPE_PROPER) .OR. & 
       (global%postDiscFlag .EQV. .TRUE.) .OR. & 
       (global%postVortFlag .EQV. .TRUE.) .OR. & 
       (global%postVortCoreFlag .EQV. .TRUE.) ) THEN
    RFLU_DecideBuildGeometry = .TRUE.
  END IF ! global%initFlowFlag 

! ******************************************************************************
! End
! ******************************************************************************

END FUNCTION RFLU_DecideBuildGeometry

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DecideBuildGeometry.F90,v $
! Revision 1.7  2008/12/06 08:45:05  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.4  2005/12/02 22:12:43  haselbac
! Bug fix: Added test for postVort{Core}Flag
!
! Revision 1.3  2005/06/02 17:44:06  haselbac
! Bug fix: Added test for postDiscFlag
!
! Revision 1.2  2005/04/22 15:23:54  haselbac
! Included extraction flag in making decision
!
! Revision 1.1  2004/07/21 14:58:55  haselbac
! Initial revision
!
! ******************************************************************************






