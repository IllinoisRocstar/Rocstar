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
! Purpose: Allocate memory for Equilibrium Eulerian variables of species 
!   solution.
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
! $Id: SPEC_RFLU_AllocateMemoryEEv.F90,v 1.3 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_AllocateMemoryEEv(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  
  USE SPEC_ModParameters

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
  INTEGER :: errorFlag,icg,iSpecEE,nSpeciesEE
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_AllocateMemoryEEv.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_AllocateMemoryEEv',&
  'SPEC_RFLU_AllocateMemoryEEv.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid  => pRegion%grid

  nSpeciesEE = pRegion%specInput%nSpeciesEE

! ******************************************************************************
! Allocate and initialize memory
! ******************************************************************************

  IF ( nSpeciesEE > 0 ) THEN 
    ALLOCATE(pRegion%spec%eev(EEV_SPEC_XVEL:EEV_SPEC_TEMP,nSpeciesEE, & 
                              pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%eev')
    END IF ! global%error
    
    DO icg = 1,pGrid%nCellsTot
      DO iSpecEE = 1,nSpeciesEE
        pRegion%spec%eev(EEV_SPEC_XVEL,iSpecEE,icg) = 0.0_RFREAL
        pRegion%spec%eev(EEV_SPEC_YVEL,iSpecEE,icg) = 0.0_RFREAL
        pRegion%spec%eev(EEV_SPEC_ZVEL,iSpecEE,icg) = 0.0_RFREAL
        pRegion%spec%eev(EEV_SPEC_TEMP,iSpecEE,icg) = 0.0_RFREAL                        
      END DO ! iSpecEE
    END DO ! icg
  END IF ! nSpeciesEE

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_AllocateMemoryEEv

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_AllocateMemoryEEv.F90,v $
! Revision 1.3  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/11/27 01:47:26  haselbac
! Initial revision
!
! ******************************************************************************







