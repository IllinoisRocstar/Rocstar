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
! Purpose: Determine whether to build stencils and gradient weights.
!
! Description: None.
!
! Input:
!   global                              Pointer to global data
!
! Output: 
!   RFLU_DecideBuildStencilsWeights     TRUE or FALSE
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_DecideBuildStencilsWeights.F90,v 1.4 2008/12/06 08:45:05 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_DecideBuildStencilsWeights(global)

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

  RCSIdentString = '$RCSfile: RFLU_DecideBuildStencilsWeights.F90,v $ $Revision: 1.4 $'

! ******************************************************************************
! Initialize
! ******************************************************************************

  RFLU_DecideBuildStencilsWeights = .FALSE.
  
! ******************************************************************************
! Determine whether should build stencils. NOTE need initFlowFlag here also 
! because this signals that errors can be computed...
! ******************************************************************************

  IF ( (global%postDiscFlag .EQV. .TRUE.) .OR. & 
       (global%postGradFlag .EQV. .TRUE.) .OR. & 
       (global%postVortFlag .EQV. .TRUE.) .OR. & 
       (global%postVortCoreFlag .EQV. .TRUE.) ) THEN
    RFLU_DecideBuildStencilsWeights = .TRUE.
  END IF ! global%postDiscFlag

! ******************************************************************************
! End
! ******************************************************************************

END FUNCTION RFLU_DecideBuildStencilsWeights

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DecideBuildStencilsWeights.F90,v $
! Revision 1.4  2008/12/06 08:45:05  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.1  2006/01/06 22:02:58  haselbac
! Initial revision
!
! ******************************************************************************






