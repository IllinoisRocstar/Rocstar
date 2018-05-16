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
! Purpose: Determine whether need boundary face gradients and related functions
!
! Description: None.
!
! Input:
!   pRegion			Pointer to region data
!   pPatch                      Pointer to Patch data
!
! Output: 
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_DecideNeedBGradFace.F90,v 1.3 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

LOGICAL FUNCTION RFLU_DecideNeedBGradFace(pRegion,pPatch)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
      
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  RCSIdentString = '$RCSfile: RFLU_DecideNeedBGradFace.F90,v $'

  CALL RegisterFunction(global,'RFLU_DecideNeedBGradFace',&
  'RFLU_DecideNeedBGradFace.F90')

! *****************************************************************************
! Initialize
! *****************************************************************************

  RFLU_DecideNeedBGradFace = .FALSE.
  
! *****************************************************************************
! Determine whether need boundary face gradients 
! *****************************************************************************

! =============================================================================
! Flow Model - true if Navier Stokes
! =============================================================================

  IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
    RFLU_DecideNeedBGradFace = .TRUE.
  END IF ! pRegion%mixtInput%flowModel 

! =============================================================================
! Boundary Condition Kind - true if NSCBC
! =============================================================================

  IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
    RFLU_DecideNeedBGradFace = .TRUE.
  END IF ! pPatch%bcKind

! =============================================================================
! Virtual boundary - always false
! =============================================================================

  IF ( pPatch%bcType == BC_VIRTUAL ) THEN
    RFLU_DecideNeedBGradFace = .FALSE.
  END IF ! pPatch%bcType

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_DecideNeedBGradFace

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DecideNeedBGradFace.F90,v $
! Revision 1.3  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/08/19 15:37:35  mparmar
! Initial revision
!
!******************************************************************************







