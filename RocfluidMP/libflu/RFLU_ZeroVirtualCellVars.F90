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
! Purpose: Zero variables in virtual cells.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ZeroVirtualCellVars.F90,v 1.5 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ZeroVirtualCellVars(pRegion,var)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY : t_region
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), POINTER :: var(:,:)
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ZeroVirtualCellVars.F90,v $ $Revision: 1.5 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ZeroVirtualCellVars',&
  'RFLU_ZeroVirtualCellVars.F90')

! ******************************************************************************
! Compute cell gradients
! ******************************************************************************

  DO icg = pRegion%grid%nCells+1,pRegion%grid%nCellsTot
    var(:,icg) = 0.0_RFREAL
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ZeroVirtualCellVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ZeroVirtualCellVars.F90,v $
! Revision 1.5  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.2  2005/05/17 21:49:41  haselbac
! Bug fix: Missing POINTER qualifier for pRegion
!
! Revision 1.1  2005/05/16 20:36:28  haselbac
! Initial revision
!
! ******************************************************************************







