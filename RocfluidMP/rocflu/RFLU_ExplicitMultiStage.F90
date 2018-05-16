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
! Purpose: Calculate solution at a new time level/iteration.
!
! Description: None.
!
! Input: 
!   regions     Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ExplicitMultiStage.F90,v 1.7 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ExplicitMultiStage(regions)

  USE ModDataTypes
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt
  USE ModMPI

  USE RFLU_ModExplicitMultiStage
  USE RFLU_ModMPI
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformWrapper
  
  USE ModInterfaces, ONLY: RFLU_CheckPositivityWrapper, &
                           RFLU_CheckValidityWrapper, &
                           RFLU_SetVars, & 
                           RFLU_ZeroVirtualCellVars

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iReg,iStage
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt), POINTER :: pMixt
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_ExplicitMultiStage',&
  'RFLU_ExplicitMultiStage.F90')

! ******************************************************************************
! Loop over stages
! ******************************************************************************

  DO iStage = 1,global%nrkSteps
  
! ==============================================================================
!   Loop over regions, update and send variables
! ==============================================================================
 
    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
      pMixt   => pRegion%mixt
      
      pRegion%irkStep = iStage

! ------------------------------------------------------------------------------
!     Store previous solution
! ------------------------------------------------------------------------------

      IF ( iStage == 1 ) THEN
        CALL RFLU_EMS_SetCvOld(pRegion)
      END IF ! iStage

! ------------------------------------------------------------------------------
!     Compute residual and zero residuals in dummy cells 
! ------------------------------------------------------------------------------

      CALL RFLU_EMS_ComputeResidual(pRegion)
      CALL RFLU_ZeroVirtualCellVars(pRegion,pMixt%rhs)

! ------------------------------------------------------------------------------
!     Update conserved variables, check validity and positivity
! ------------------------------------------------------------------------------

      CALL RFLU_EMS_UpdateConservedVars(pRegion)
      CALL RFLU_CheckValidityWrapper(pRegion)
      CALL RFLU_CheckPositivityWrapper(pRegion)

! ------------------------------------------------------------------------------
!     Send variables to other regions
! ------------------------------------------------------------------------------

      CALL RFLU_MPI_ISendWrapper(pRegion)

! ------------------------------------------------------------------------------
!     Update other variables
! ------------------------------------------------------------------------------

      CALL RFLU_SetVars(pRegion,1,pGrid%nCells)
    END DO ! iReg

! ==============================================================================
!   Copy variables between regions on same process
! ==============================================================================

    CALL RFLU_MPI_CopyWrapper(regions)

! ==============================================================================
!   Loop over regions, receive variables and update 
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

! ------------------------------------------------------------------------------
!     Receive variables from other regions
! ------------------------------------------------------------------------------

      CALL RFLU_MPI_RecvWrapper(pRegion)    

! ------------------------------------------------------------------------------
!     Update other variables in virtual cells
! ------------------------------------------------------------------------------

      CALL RFLU_SetVars(pRegion,pGrid%nCells+1,pGrid%nCellsTot) 
      
! ------------------------------------------------------------------------------
!     Transform variables in virtual cells associated with related patches
! ------------------------------------------------------------------------------

      CALL RFLU_RELP_TransformWrapper(pRegion)            
    END DO ! iReg

! ==============================================================================
!   Loop over regions, clear requests
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg
  END DO ! iStage

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExplicitMultiStage

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ExplicitMultiStage.F90,v $
! Revision 1.7  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.4  2006/03/26 20:28:43  haselbac
! Changed to wrappers bcos of GL model
!
! Revision 1.3  2006/03/25 22:01:35  haselbac
! Added call to transforming data on related patches
!
! Revision 1.2  2005/12/03 19:45:20  haselbac
! Apparent bug fix: Separated call to RFLU_MPI_ClearRequestWrapper into separate loop
!
! Revision 1.1  2005/05/16 20:36:29  haselbac
! Initial revision
!
! ******************************************************************************







