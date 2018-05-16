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
! Purpose: Allocate memory for species solution.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_AllocateMemorySol.F90,v 1.6 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_AllocateMemorySol(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iSpec,nSpecies
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_AllocateMemorySol.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_AllocateMemorySol',&
  'SPEC_RFLU_AllocateMemorySol.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid  => pRegion%grid

  nSpecies = pRegion%specInput%nSpecies

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! State vectors
! ==============================================================================

  ALLOCATE(pRegion%spec%cv(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%cv')
  END IF ! global%error

  ALLOCATE(pRegion%spec%cvInfo(nSpecies),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%cvInfo')
  END IF ! global%error

! ==============================================================================
! Transport variables
! ==============================================================================

  IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
    ALLOCATE(pRegion%spec%tv(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%tv')
    END IF ! global%error
  ELSE
    NULLIFY(pRegion%spec%tv)
  END IF ! pMixtInput%flowModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_AllocateMemorySol

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_AllocateMemorySol.F90,v $
! Revision 1.6  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/11/02 02:33:52  haselbac
! Removed init of cvInfo, cosmetics
!
! Revision 1.3  2004/03/03 23:54:41  jferry
! Corrected typo
!
! Revision 1.2  2004/01/29 22:59:37  haselbac
! Added allocation of cvInfo and tv
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







