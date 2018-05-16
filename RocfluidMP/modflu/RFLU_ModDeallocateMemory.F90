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
!*******************************************************************************
!
! Purpose: Suite of routines to deallocate memory.
!
! Description: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: RFLU_ModDeallocateMemory.F90,v 1.11 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModDeallocateMemory

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_DeallocateMemoryGSpeeds, &
            RFLU_DeallocateMemorySol, &
            RFLU_DeallocateMemorySolCv, &
            RFLU_DeallocateMemorySolDv, &
            RFLU_DeallocateMemorySolGv, &
            RFLU_DeallocateMemorySolTv, &
            RFLU_DeallocateMemoryTStep, & 
            RFLU_DeallocateMemoryTStep_C, & 
            RFLU_DeallocateMemoryTStep_I

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString = &
    '$RCSfile: RFLU_ModDeallocateMemory.F90,v $ $Revision: 1.11 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Deallocate memory for grid speeds.
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

SUBROUTINE RFLU_DeallocateMemoryGSpeeds(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryGSpeeds',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Grid motion active
! ==============================================================================

  IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN

! ------------------------------------------------------------------------------
!   Interior faces
! ------------------------------------------------------------------------------

    DEALLOCATE(pGrid%gs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%gs')
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Patch faces
! ------------------------------------------------------------------------------

    IF ( pGrid%nPatches > 0 ) THEN
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        DEALLOCATE(pPatch%gs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%gs')
        END IF ! global%error
      END DO ! iPatch
    END IF ! pGrid%nPatches

! ==============================================================================
! Grid motion not active
! ==============================================================================

  ELSE

! ------------------------------------------------------------------------------
!   Interior faces
! ------------------------------------------------------------------------------

    DEALLOCATE(pGrid%gs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%gs')
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Patch faces
! ------------------------------------------------------------------------------

    IF ( pGrid%nPatches > 0 ) THEN
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        DEALLOCATE(pPatch%gs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%gs')
        END IF ! global%error
      END DO ! iPatch
    END IF ! pGrid%nPatches
  END IF ! pMixtInput%moveGrid

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryGSpeeds






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture solution.
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

SUBROUTINE RFLU_DeallocateMemorySol(pRegion)

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
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySol',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  CALL RFLU_DeallocateMemorySolCv(pRegion)
  CALL RFLU_DeallocateMemorySolDv(pRegion)
  CALL RFLU_DeallocateMemorySolGv(pRegion)
  CALL RFLU_DeallocateMemorySolTv(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySol






! ******************************************************************************
!
! Purpose: Deallocate memory for conserved variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolCv(pRegion)

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
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolCv',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%cv,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cv')
  END IF ! global%error

  DEALLOCATE(pRegion%mixt%cvInfo,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvInfo')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolCv






! ******************************************************************************
!
! Purpose: Deallocate memory for dependent variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolDv(pRegion)

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
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolDv',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%dv,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%dv')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolDv






! ******************************************************************************
!
! Purpose: Deallocate memory for gas variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolGv(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolGv',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%gv,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gv')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolGv





! ******************************************************************************
!
! Purpose: Deallocate memory for transport variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolTv(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolTv',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  IF ( pMixtInput%computeTv .EQV. .TRUE. ) THEN
    DEALLOCATE(pRegion%mixt%tv,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tv')
    END IF ! global%error
  END IF ! pMixtInput%computeTv

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolTv






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture time-stepping.
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

SUBROUTINE RFLU_DeallocateMemoryTStep(pRegion)

  USE RFLU_ModOLES

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

  INTEGER :: errorFlag,iPatch,nBFaces,nBFacesTot
  TYPE(t_grid), POINTER :: pGrid,pGridOld,pGridOld2
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryTStep',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pGridOld   => pRegion%gridOld
  pGridOld2  => pRegion%gridOld2
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Old solutions
! ==============================================================================

  DEALLOCATE(pRegion%mixt%cvOld,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld')
  END IF ! global%error

  IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
    DEALLOCATE(pRegion%mixt%cvOld1,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld1')
    END IF ! global%error

    DEALLOCATE(pRegion%mixt%cvOld2,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld2')
    END IF ! global%error
  END IF ! global%solverType

! ==============================================================================
! Time step
! ==============================================================================

  DEALLOCATE(pRegion%dt,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%dt')
  END IF ! global%error

! ==============================================================================
! Residuals
! ==============================================================================

  DEALLOCATE(pRegion%mixt%rhs,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%rhs')
  END IF ! global%error

  DEALLOCATE(pRegion%mixt%diss,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%diss')
  END IF ! global%error

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    DEALLOCATE(pRegion%mixt%rhsSum,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%rhsSum')
    END IF ! global%error
  END IF ! global%flowType

! ==============================================================================
! Gradients
! ==============================================================================

! ------------------------------------------------------------------------------
! Cell gradients
! ------------------------------------------------------------------------------

  IF ( (pMixtInput%spaceDiscr == DISCR_UPW_ROE     ) .OR. &
       (pMixtInput%spaceDiscr == DISCR_UPW_HLLC    ) .OR. & 
       (pMixtInput%spaceDiscr == DISCR_UPW_AUSMPLUS) ) THEN
    IF ( pMixtInput%spaceOrder > 1 ) THEN
      DEALLOCATE(pRegion%mixt%gradCell,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCell')
      END IF ! global%error
    END IF ! pMixtInput%spaceOrder
  ELSE IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN
    DEALLOCATE(pRegion%mixt%gradCell,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCell')
    END IF ! global%error
  END IF ! pMixtInput%spaceDiscr

! ------------------------------------------------------------------------------
! Face gradients
! ------------------------------------------------------------------------------

  IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
    DEALLOCATE(pRegion%mixt%gradFace,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradFace')
    END IF ! global%error
  END IF ! pMixtInput%flowModel

  IF ( pGrid%nPatches > 0 ) THEN
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)
  
      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        DEALLOCATE(pPatch%mixt%gradFace,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mixt%gradFace')
        END IF ! global%error
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
  END IF ! pGrid%nPatches

! ==============================================================================
! Grid motion. NOTE grid speeds are allocated separately because they are
! written into grid file, and hence they need to be allocated in pre- and
! postprocessors also.
! ==============================================================================

  IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN

! ------------------------------------------------------------------------------
!   Residual
! ------------------------------------------------------------------------------

    DEALLOCATE(pGrid%rhs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%rhs')
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Displacement
! ------------------------------------------------------------------------------

    IF ( pMixtInput%moveGridType /= MOVEGRID_TYPE_XYZ ) THEN
      DEALLOCATE(pGrid%disp,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%disp')
      END IF ! global%error
    END IF ! pMixtInput%moveGridType

! ------------------------------------------------------------------------------
!   Old coordinates
! ------------------------------------------------------------------------------

    DEALLOCATE(pGridOld%xyz,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld%xyz')
    END IF ! global%error

    IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
      DEALLOCATE(pGridOld2%xyz,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld2%xyz')
      END IF ! global%error
    END IF ! global%solverType

! ------------------------------------------------------------------------------
!   Old volume
! ------------------------------------------------------------------------------

    DEALLOCATE(pGridOld%vol,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld%vol')
    END IF ! global%error

    IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
      DEALLOCATE(pGridOld2%vol,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld2%vol')
      END IF ! global%error
    END IF ! global%solverType

#ifndef GENX
! ------------------------------------------------------------------------------
!   Patch displacements. NOTE allocate here only if not running inside GENX,
!   because when running inside GENX allocate displacements also for virtual
!   vertices.
! ------------------------------------------------------------------------------

    IF ( pGrid%nPatches > 0 ) THEN
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        DEALLOCATE(pPatch%dXyz,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%dXyz')
        END IF ! global%error

      END DO ! iPatch
    END IF ! pGrid%nPatches
#endif
  END IF ! pMixtInput%moveGrid

#ifdef STATS
! ==============================================================================
! Time averaged statistics
! ==============================================================================

  IF ( (global%flowType == FLOW_UNSTEADY) .AND. &
       (global%doStat == ACTIVE) ) THEN
    DEALLOCATE(pRegion%mixt%tav,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tav')
    END IF ! global%error
  END IF ! global%flowType
#endif

! ==============================================================================
! Optimal LES
! ==============================================================================

  IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN
    CALL RFLU_DestroyStencilsWeightsOLES(pRegion)
  END IF ! pMixtInput%spaceDiscr

! ==============================================================================
! Substantial derivative
! ==============================================================================

  DEALLOCATE(pRegion%mixt%sd,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%sd')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryTStep






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture time-stepping for compressible fluid.
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

SUBROUTINE RFLU_DeallocateMemoryTStep_C(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryTStep_C',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Mass fluxes
! ==============================================================================

  DEALLOCATE(pRegion%mixt%mfMixt,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%mfMixt')
  END IF ! global%error

  IF ( pRegion%grid%nPatches > 0 ) THEN
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      DEALLOCATE(pPatch%mfMixt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mfMixt')
      END IF ! global%error
    END DO ! iPatch
  END IF ! pRegion%grid%nPatches

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryTStep_C






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture time-stepping for incompressible fluid.
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

SUBROUTINE RFLU_DeallocateMemoryTStep_I(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryTStep_I',&
  'RFLU_ModDeallocateMemory.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Face velocities
! ==============================================================================

  DEALLOCATE(pRegion%mixt%vfMixt,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%vfMixt')
  END IF ! global%error

  IF ( pRegion%grid%nPatches > 0 ) THEN
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      DEALLOCATE(pPatch%vfMixt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%vfMixt')
      END IF ! global%error
    END DO ! iPatch
  END IF ! pRegion%grid%nPatches

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryTStep_I





! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModDeallocateMemory


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModDeallocateMemory.F90,v $
! Revision 1.11  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/11/01 15:50:00  haselbac
! Changed so implicit-solver arrays dealt with properly
!
! Revision 1.8  2006/08/19 15:39:03  mparmar
! Moved region%mixt%bGradFace to patch%mixt%gradFace
!
! Revision 1.7  2006/02/08 21:03:47  hdewey2
! Added old2 quantities
!
! Revision 1.6  2005/09/22 17:11:04  hdewey2
! Added deallocation of cvOld1 and cvOld2 for transient implicit solver.
!
! Revision 1.5  2005/07/14 21:43:26  haselbac
! Added AUSM flux function to memory deallocation IF statement
!
! Revision 1.4  2004/12/19 15:47:08  haselbac
! Added memory deallocation for incompressible solver
!
! Revision 1.3  2004/10/19 19:38:48  haselbac
! Updated for GEN3
!
! Revision 1.2  2004/07/30 22:47:36  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.1  2004/03/19 21:15:19  haselbac
! Initial revision
!
! ******************************************************************************















