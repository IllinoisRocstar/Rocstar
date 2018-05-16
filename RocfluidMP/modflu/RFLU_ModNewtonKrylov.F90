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
! Purpose: Collection of routines related to Newton-Krylov schemes.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModNewtonKrylov.F90,v 1.16 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModNewtonKrylov

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch 
  USE ModMixture, ONLY: t_mixt,t_mixt_input
  USE ModMPI
  
  USE RFLU_ModDimensions
  USE RFLU_ModForcesMoments
  USE RFLU_ModGeometry
  USE RFLU_ModMPI
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModPETScNewtonKrylov  
  USE RFLU_ModProbes
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModWeights  
  
  USE ModInterfaces, ONLY: IntegrateSourceTermsMP, &
                           RFLU_ComputeGridSpeeds, &
                           RFLU_ComputeIntegralValues, & 
                           RFLU_CheckGridSpeeds, &
                           RFLU_DecideNeedBGradFace, &
                           RFLU_DecidePrint, & 
                           RFLU_DecideWrite, & 
                           RFLU_ExplicitMultiStage, &
                           RFLU_GetDeformationWrapper, & 
                           RFLU_MinimumTimeStep, &
                           RFLU_MoveGridWrapper, &
                           RFLU_PrintChangeInfo, &  
                           RFLU_PrintGridInfo, &
                           RFLU_PrintFlowInfoWrapper, &
                           RFLU_PrintFlowInfo, &
                           RFLU_PrintWriteConvergence, &
                           RFLU_ResidualNorm, &
                           RFLU_TimeStepInviscid, &
                           RFLU_TimeStepViscous, &
                           RFLU_WriteRestartInfo, &
                           RFLU_WriteStatsFileOLES, &
                           RungeKuttaMP, &
                           WriteProbe, &
                           WriteTotalMass, &
                           RFLU_CheckPositivity, &
                           RFLU_CheckValidity, &
                           RFLU_SetDependentVars

  IMPLICIT NONE

  INCLUDE 'mpif.h'

! DEBUG
!#ifdef GENX
!  INCLUDE 'rocmanf90.h'
!#endif
! END DEBUG

  PRIVATE
  PUBLIC :: RFLU_NK_TimeStepping
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModNewtonKrylov.F90,v $ $Revision: 1.16 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  


! ******************************************************************************
!
! Purpose: Integrate the governing equations in time.
!
! Description: None.
!
! Input: 
!   dTimeSystem         Total solution time (unsteady flow)
!   dIterSystem         Total number of iterations (steady flow)
!   regions             Data for all grid regions
!
! Output: 
!   regions             Data for all grid regions
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NK_TimeStepping(dTimeSystem,dIterSystem,regions)

  USE ModDataTypes
  USE ModParameters
  USE ModError  
  USE ModGlobal, ONLY: t_global    
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModMixture, ONLY: t_mixt_input    
  USE ModMPI
  
  USE RFLU_ModPETScNewtonKrylov
  
  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscsnes.h"

! ******************************************************************************
! INTERFACE blocks for PETSc SNES callback definition functions
! ******************************************************************************

  INTERFACE
    SUBROUTINE SNESSetFunction(snes,r,fun,ctx,ierr)
      USE ModDataStruct, ONLY: t_region
      SNES :: snes
      Vec :: r
      EXTERNAL :: fun
      TYPE(t_region),POINTER :: ctx
      INTEGER :: ierr
    END SUBROUTINE
    SUBROUTINE SNESSetJacobian(snes,A,pA,fun,ctx,ierr)
      USE ModDataStruct, ONLY: t_region
      SNES :: snes
      Mat :: A, pA
      EXTERNAL :: fun
      TYPE(t_region),POINTER :: ctx
      INTEGER :: ierr
    END SUBROUTINE
    SUBROUTINE MatFDColoringSetFunctionSNES(matfd,fun,ctx,ierr)
      USE ModDataStruct, ONLY: t_region
      MatFDColoring :: matfd
      EXTERNAL :: fun
      TYPE(t_region),POINTER :: ctx
      INTEGER :: ierr
    END SUBROUTINE
  END INTERFACE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL :: doPrint,doProbe,doWrite
  LOGICAL :: finished = .FALSE., moveGrid
  INTEGER :: dIterSystem
  REAL(RFREAL) :: dTimeSystem
  TYPE(t_region), POINTER :: pRegion,regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: gtemp,ic,ierr,iPatch,iReg,iter,iv,ltemp
  TYPE(t_global), POINTER :: global
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_patch), POINTER :: pPatch
  SNESConvergedReason reason
  PetscLogDouble time

! DEBUG
!#ifdef GENX
!  DOUBLE PRECISION, DIMENSION(MAN_INTEG_SIZE) :: integ
!#endif
! END DEBUG

! TEMPORARY
  INTEGER :: iSub,ivl,iPatch,nSub
  REAL(RFREAL) :: iNSub, iSubTerm
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_grid), POINTER :: pGrid, pGridOld, pGridOld2
! END TEMPORARY

! ******************************************************************************
! Start, set pointers and variables
! ******************************************************************************

  global  => regions(1)%global
  pRegion => regions(1)

  nSub  = 100
  iNSub = 1.0_RFREAL/nSub

  CALL RegisterFunction(global,'RFLU_NK_TimeStepping',&
  'RFLU_ModNewtonKrylov.F90')
  
! ==============================================================================
! Determine whether have moving grids
! ==============================================================================

  moveGrid = .FALSE.

  DO iReg = 1,global%nRegionsLocal
    IF ( regions(iReg)%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      moveGrid = .TRUE.
    END IF ! regions
  END DO ! iReg

! ==============================================================================
! Set SNES callback functions
! ==============================================================================

  CALL SNESSetFunction(pRegion%snes,pRegion%r, &
       RFLU_PETSC_FormResidual,pRegion,ierr)
  CALL SNESSetJacobian(pRegion%snes,pRegion%A,pRegion%preA, &
       RFLU_PETSC_FormJacobian,pRegion,ierr)
  CALL MatFDColoringSetFunctionSNES(pRegion%fdcolor, &
       RFLU_PETSC_FormResidualFirstOrder,pRegion,ierr)  

! ==============================================================================
! Initialize old variables as current variables
! ==============================================================================

  pRegion%mixt%cvOld(:,:)   = pRegion%mixt%cv(:,:)
  pRegion%mixt%cvOld1(:,:)  = pRegion%mixt%cv(:,:)
  pRegion%mixt%cvOld2(:,:)  = pRegion%mixt%cv(:,:)
  pRegion%gridOld%vol(:)    = pRegion%grid%vol(:)
  pRegion%gridOld2%vol(:)   = pRegion%grid%vol(:)
  pRegion%gridOld%xyz(:,:)  = pRegion%grid%xyz(:,:)
  pRegion%gridOld2%xyz(:,:) = pRegion%grid%xyz(:,:)

! ******************************************************************************
! Loop over iterations/time steps
! ****************************************************************************** 

! DEBUG
  IF ( global%myProcid == 0 ) THEN
    CALL PetscGetTime(time,ierr)
    print*,'PETSC TIME 1 : ',time
  END IF ! global%myProcid
! END DEBUG

  finished = .FALSE.

  DO

! ==============================================================================
!   Update iteration counter for steady flow.
! ==============================================================================

    IF ( global%flowType == FLOW_STEADY ) THEN
      global%currentIter      = global%currentIter      + 1
      global%iterSinceRestart = global%iterSinceRestart + 1
    END IF ! global%flowType

! ==============================================================================
!   Compute time step.
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal  
      pRegion => regions(iReg)

      IF ( pRegion%mixtInput%flowModel == FLOW_EULER ) THEN  
        CALL RFLU_TimeStepInviscid(pRegion)
      ELSE 
        CALL RFLU_TimeStepViscous(pRegion)
      END IF ! pRegion
    END DO ! iReg

    global%dtmin = global%dtImposed

! ==============================================================================
!   Move or generate new grid
! ==============================================================================

    IF ( ( global%flowType == FLOW_UNSTEADY ) .AND.  & 
         ( global%iterSinceRestart > 0 ) ) THEN
      IF ( moveGrid .EQV. .TRUE. ) THEN

! ------------------------------------------------------------------------------
!       Get the displacements from GENX (or other source).
! ------------------------------------------------------------------------------

        CALL RFLU_GetDeformationWrapper(regions)

! ------------------------------------------------------------------------------
!       Store the old grid.
! ------------------------------------------------------------------------------

        DO iReg = 1,global%nRegionsLocal 
          pGrid       => regions(iReg)%grid
          pGridOld    => regions(iReg)%gridOld
          pGridOld2   => regions(iReg)%gridOld2    
           
          DO ic = 1,pGrid%nCellsTot ! Explicit copy to avoid ASCI White problem
            pGridOld2%vol(ic) = pGridOld%vol(ic)
          END DO ! ic  
           
          DO iv = 1,pGrid%nVertTot ! Explicit copy to avoid ASCI White problem    
            pGridOld2%xyz(XCOORD,iv) = pGridOld%xyz(XCOORD,iv)
            pGridOld2%xyz(YCOORD,iv) = pGridOld%xyz(YCOORD,iv)
            pGridOld2%xyz(ZCOORD,iv) = pGridOld%xyz(ZCOORD,iv)            
          END DO ! iv 
        END DO ! iReg
        
        DO iReg = 1,global%nRegionsLocal 
          pGrid    => regions(iReg)%grid
          pGridOld => regions(iReg)%gridOld    
          
          DO ic = 1,pGrid%nCellsTot ! Explicit copy to avoid ASCI White problem
            pGridOld%vol(ic) = pGrid%vol(ic)
          END DO ! ic  
           
          DO iv = 1,pGrid%nVertTot ! Explicit copy to avoid ASCI White problem    
            pGridOld%xyz(XCOORD,iv) = pGrid%xyz(XCOORD,iv)
            pGridOld%xyz(YCOORD,iv) = pGrid%xyz(YCOORD,iv)
            pGridOld%xyz(ZCOORD,iv) = pGrid%xyz(ZCOORD,iv)            
          END DO ! iv 
        END DO ! iReg

! ------------------------------------------------------------------------------
!       We will be substepping Mesquite, so start with a smaller displacement.
! ------------------------------------------------------------------------------

        DO iReg = 1,global%nRegionsLocal 
          pGrid => regions(iReg)%grid    

          DO iPatch = 1,pGrid%nPatches
            pPatch => regions(iReg)%patches(iPatch)
             
            DO ivl = 1,pPatch%nBVert
              pPatch%dXyz(XCOORD,ivl) = iNSub*pPatch%dXyz(XCOORD,ivl)
              pPatch%dXyz(YCOORD,ivl) = iNSub*pPatch%dXyz(YCOORD,ivl)
              pPatch%dXyz(ZCOORD,ivl) = iNSub*pPatch%dXyz(ZCOORD,ivl)
            END DO ! ivl
          END DO ! iPatch         
        END DO ! iReg

! ------------------------------------------------------------------------------
!       Substep displacements and smoothing to avoid inverting elements.
! ------------------------------------------------------------------------------

        DO iSub = 1,nSub
          iSubTerm = (1.0_RFREAL*iSub)/(1.0_RFREAL*iSub-1.0_RFREAL)

! ------- Displace the mesh by steadily increasing displacements until the 
!         original displacement is reached.
           
          DO iReg = 1,global%nRegionsLocal 
            pGrid => regions(iReg)%grid    

            DO iPatch = 1,pGrid%nPatches
              pPatch => regions(iReg)%patches(iPatch)
                 
              DO ivl = 1,pPatch%nBVert 
                IF ( iSub > 1 ) THEN
                  pPatch%dXyz(XCOORD,ivl) = iSubTerm*pPatch%dXyz(XCOORD,ivl)
                  pPatch%dXyz(YCOORD,ivl) = iSubTerm*pPatch%dXyz(YCOORD,ivl)
                  pPatch%dXyz(ZCOORD,ivl) = iSubTerm*pPatch%dXyz(ZCOORD,ivl)
                END IF ! iSub
              END DO ! ivl
            END DO ! iPatch         
          END DO ! iReg

! ------- Have Mesquite smooth the grid at each substep. NOTE:  
!         This function will not store the old grid. That must be done before 
!         Mesquite substepping begins.

          CALL RFLU_MoveGridWrapper(regions)
        END DO ! iSub

! ------------------------------------------------------------------------------
!       Finish by generating geometry and grid speeds for the new grid.
! ------------------------------------------------------------------------------

        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_BuildGeometry(pRegion)
          CALL RFLU_ComputeGridSpeeds(pRegion)

          IF ( global%checkLevel == CHECK_HIGH ) THEN
            CALL RFLU_CheckGridSpeeds(pRegion)
          END IF ! global%checkLevel
        END DO ! iReg
      END IF ! moveGrid
    END IF ! global%flowType 

! ==============================================================================
!   Recompute weights
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( moveGrid .EQV. .TRUE. ) THEN
        DO iReg=1,global%nRegionsLocal
          pRegion => regions(iReg)
          pMixtInput => pRegion%mixtInput

          IF ( pMixtInput%spaceOrder > 1 ) THEN
            CALL RFLU_ComputeWtsC2CWrapper(pRegion,pMixtInput%spaceOrder-1)
          END IF ! pMixtInput%spaceOrder

          IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
            CALL RFLU_ComputeWtsF2CWrapper(pRegion,pMixtInput%spaceOrder-1)
          END IF ! pMixtInput%flowModel
            
          DO iPatch = 1,pRegion%grid%nPatches
            pPatch => pRegion%patches(iPatch)

            IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
              CALL RFLU_ComputeWtsBF2CWrapper(pRegion,pPatch, &
                                              pPatch%spaceOrder)
            END IF ! RFLU_DecideNeedBGradFace
          END DO ! iPatch  

        END DO ! iReg
      END IF ! moveGrid
    END IF ! global%flowType

! ==============================================================================
!   Relocate probes
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( moveGrid .EQV. .TRUE. ) THEN
        IF ( global%nProbes > 0 ) THEN
          DO iReg = 1,global%nRegionsLocal
            pRegion => regions(iReg)
            CALL RFLU_FindProbeCells(pRegion)
          END DO ! iReg

          CALL RFLU_PrintProbeInfo(global)
        END IF ! global
      END IF ! moveGrid
    END IF ! global%flowType

! ==============================================================================
!   Compute new solution.
! ==============================================================================

    IF ( global%flowType == FLOW_STEADY ) THEN

! ------------------------------------------------------------------------------
!   For steady-state problems, simply compute the new solution
! ------------------------------------------------------------------------------

      global%forceX  = 0.0_RFREAL
      global%forceY  = 0.0_RFREAL
      global%forceZ  = 0.0_RFREAL
      global%massIn  = 0.0_RFREAL
      global%massOut = 0.0_RFREAL
      pRegion => regions(1)
      pRegion%mixt%cvOld(:,:)  = pRegion%mixt%cv(:,:)
      pRegion%irkStep = 1
      CALL SNESSolve(pRegion%snes,PETSC_NULL_OBJECT,pRegion%x,ierr)
      CALL SNESGetIterationNumber(pRegion%snes,iter,ierr)
!      CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCells)

      CALL RFLU_MPI_ISendWrapper(pRegion)
!        CALL RFLU_MPI_CopyWrapper(regions)
      CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCells)
      CALL RFLU_MPI_RecvWrapper(pRegion)
      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
      CALL RFLU_SetDependentVars(pRegion,pRegion%grid%nCells+1,pRegion%grid%nCellsTot)

    ELSE

! ------------------------------------------------------------------------------
!   For transient problems, first update the iteration variables for the
!   pseudo steady-state problem.
! ------------------------------------------------------------------------------

      pRegion => regions(1)

      iter = global%iterSinceRestart
      global%currentIter      = 0
      global%iterSinceRestart = 0

      pRegion%mixt%cvOld2(:,:)  = pRegion%mixt%cvOld1(:,:)
      pRegion%mixt%cvOld1(:,:)  = pRegion%mixt%cv(:,:)

      DO

! ------------------------------------------------------------------------------
!   Repeatedly solve the pseudo steady-state problem until the residual drops
!   below the residual tolerance.
! ------------------------------------------------------------------------------

        global%currentIter      = global%currentIter      + 1
        global%iterSinceRestart = global%iterSinceRestart + 1

        DO iReg = 1,global%nRegionsLocal  
          pRegion => regions(iReg)
          IF ( pRegion%mixtInput%flowModel == FLOW_EULER ) THEN  
            CALL RFLU_TimeStepInviscid(pRegion)
          ELSE 
            CALL RFLU_TimeStepViscous(pRegion)
          END IF ! pRegion
        END DO ! iReg
        global%dtmin = global%dtImposed

        pRegion%mixt%cvOld(:,:)  = pRegion%mixt%cv(:,:)

        global%forceX  = 0.0_RFREAL
        global%forceY  = 0.0_RFREAL
        global%forceZ  = 0.0_RFREAL
        global%massIn  = 0.0_RFREAL
        global%massOut = 0.0_RFREAL

        pRegion%irkStep = 1

        CALL SNESSolve(pRegion%snes,PETSC_NULL_OBJECT,pRegion%x,ierr)
!        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCells)

        CALL RFLU_MPI_ISendWrapper(pRegion)
!        CALL RFLU_MPI_CopyWrapper(regions)
        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCells)
        CALL RFLU_MPI_RecvWrapper(pRegion)
        CALL RFLU_MPI_ClearRequestWrapper(pRegion)
        CALL RFLU_SetDependentVars(pRegion,pRegion%grid%nCells+1,pRegion%grid%nCellsTot)
        
        CALL RFLU_CheckValidity(pRegion)
        CALL RFLU_CheckPositivity(pRegion)

        CALL RFLU_ResidualNorm(regions)    

! DEBUG
        IF ( global%myProcid == 0 ) THEN
          print*,'pseudo-steady info:',global%currentIter, &
               global%residual/global%resInit,global%resTol
        END IF ! global%myProcid
! END DEBUG

        IF ( global%residual/global%resInit <= global%resTol ) THEN 
          EXIT
        END IF ! global%residual
      END DO
    END IF ! global%flowType

! ==============================================================================
!   Check for convergence problems
! ==============================================================================

   CALL SNESGetConvergedReason(pRegion%snes,reason,ierr)
   !print*,'SNES REASON = ',reason
   IF(reason == SNES_DIVERGED_FUNCTION_COUNT) THEN
     global%warnCounter = global%warnCounter + 1
     print*,'*** WARNING ***  SNES FUNCTION COUNT EXCEEDED.  ', &
          'COPYING cvOld INTO cv TO PRESERVE PARTIALLY-CONVERGED SOLUTION.'
     pRegion%mixt%cv(:,:) = pRegion%mixt%cvOld(:,:)
   ENDIF ! reason
   IF((reason /= SNES_DIVERGED_MAX_IT).AND.(reason /= SNES_CONVERGED_PNORM_RELATIVE).AND. &
        (reason /= SNES_DIVERGED_FUNCTION_COUNT)) THEN
     global%warnCounter = global%warnCounter + 1
     print*,'*** WARNING ***  SNES REASON DOES NOT INDICATE CONTINUATION OF SNES STEPPING: ',reason
   ENDIF

! DEBUG
!   CALL RFLU_PrintFlowInfo(pRegion)
! END DEBUG

! ==============================================================================
!   Check validity and positivity
! ==============================================================================

   CALL RFLU_CheckValidity(pRegion)
   CALL RFLU_CheckPositivity(pRegion)

! ==============================================================================
!   Get statistics
! ==============================================================================

#ifdef STATS
    CALL GetStatistics(regions)
#endif

! ==============================================================================
!   Reset global%timeSince* variables. NOTE this must be done right before
!   the update of the various time variables, otherwise calling the routines
!   RFLU_Decide* from within rungeKuttaMP or explicitMultiStage will not work
!   correctly.
! ==============================================================================

    IF ( RFLU_DecidePrint(global) .EQV. .TRUE. ) THEN 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN    
        IF ( global%iterSinceRestart > 1 ) THEN 
          global%timeSincePrint = 0.0_RFREAL
        END IF ! global%iterSinceRestart
      END IF ! global%flowType      
    END IF ! RFLU_DecidePrint

    IF ( RFLU_DecideWrite(global) .EQV. .TRUE. ) THEN 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN    
        global%timeSinceWrite = 0.0_RFREAL
      END IF ! global%flowType      
    END IF ! RFLU_DecideWrite

    IF ( RFLU_DecideWriteProbes(global) .EQV. .TRUE. ) THEN 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN    
        IF ( global%iterSinceRestart > 1 ) THEN 
          global%timeSinceProbe = 0.0_RFREAL
        END IF ! global%iterSinceRestart 
      END IF ! global%flowType      
    END IF ! RFLU_DecideWriteProbes  

! ==============================================================================
!   Update times for unsteady flow
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      global%currentTime      = global%currentTime      + global%dtImposed    
      global%timeSinceRestart = global%timeSinceRestart + global%dtImposed
      
      global%timeSincePrint = global%timeSincePrint + global%dtImposed
      global%timeSinceWrite = global%timeSinceWrite + global%dtImposed
      global%timeSinceProbe = global%timeSinceProbe + global%dtImposed  
      
      global%iterSinceRestart = iter + 1
    END IF ! global%flowType

! ==============================================================================
!   Decide whether to print/write convergence, data, and probes
! ==============================================================================

    doPrint = RFLU_DecidePrint(global)
    doWrite = RFLU_DecideWrite(global)
    doProbe = RFLU_DecideWriteProbes(global)

! ==============================================================================
!   Check for end of timestepping
! ==============================================================================
   
   IF ( global%flowType == FLOW_UNSTEADY ) THEN    
      IF ( global%timeSinceRestart >= dTimeSystem ) THEN 
        finished = .TRUE.
      END IF ! global%timeSinceRestart
    ELSE
      IF ( (doPrint .EQV. .TRUE.) .OR. (global%iterSinceRestart >= dIterSystem) ) THEN 
        CALL RFLU_ResidualNorm(regions)    
        IF ( (global%iterSinceRestart >= dIterSystem) .OR. &
             (global%residual/global%resInit <= global%resTol) ) THEN 
          finished = .TRUE.
        END IF ! global%iterSinceRestart      
      END IF ! doPrint
    END IF ! global%flowType

! ==============================================================================
!   Write convergence (file & screen) and total mass (file) history
! ==============================================================================

    IF ( (doPrint .EQV. .TRUE.) .OR. (finished .EQV. .TRUE.) ) THEN
      CALL RFLU_PrintWriteConvergence(global)
      
#ifndef GENX      
      IF ( moveGrid .EQV. .TRUE. ) THEN
        CALL RFLU_ComputeIntegralValues(regions) 
        CALL WriteTotalMass(regions)
      END IF ! moveGrid
#endif    
      
      DO iReg = 1,global%nRegionsLocal
        IF ( regions(iReg)%mixtInput%spaceDiscr == DISCR_OPT_LES ) THEN 
          CALL RFLU_WriteStatsFileOLES(global)
        END IF ! mixtInput
      END DO ! iReg
    END IF ! doPrint

! ==============================================================================
!   Compute forces and mass flow
! ==============================================================================

    IF ( global%forceFlag .EQV. .TRUE. ) THEN 
      IF ( (doWrite .EQV. .TRUE.) .OR. (finished .EQV. .TRUE.) ) THEN 
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_ComputeLocalForcesMoments(pRegion)     
        END DO ! iReg

        CALL RFLU_ComputeGlobalForcesMoments(regions)

        pRegion => regions(1)

        IF ( pRegion%global%myProcid == MASTERPROC ) THEN 
          CALL RFLU_PrintGlobalForcesMoments(pRegion)
          CALL RFLU_WriteGlobalForcesMoments(pRegion)
        END IF ! pRegion
      END IF ! doWrite
    END IF ! global%forceFlag

! ==============================================================================
!   Store probe data
! ==============================================================================

    IF ( global%nProbes > 0 ) THEN
      IF ( doProbe .EQV. .TRUE. ) THEN 
        DO iReg = 1,global%nRegionsLocal
          CALL WriteProbe(regions,iReg)
        END DO ! iReg
      END IF ! doProbe
    END IF ! global

! DEBUG
#ifdef GENX
!    CALL RFLU_WriteGridWrapper(pRegion)
!    CALL RFLU_WriteFlowWrapper(pRegion)
!    CALL RFLU_ComputeIntegralValues(regions,integ) 
!    print*,global%myProcid,global%currentTime,'777 777',integ(MAN_INTEG_MASS),integ(MAN_INTEG_VOL)
#endif
! END DEBUG

#ifndef GENX
! ==============================================================================
!   Store flow field (and grid if moving). Write restart info file after 
!   flow (and grid) files so that should those not be written due to reaching
!   the time limit, the restart file will not contain the time level of the 
!   incomplete flow (and grid) files.
! ==============================================================================

    IF ( (doWrite .EQV. .TRUE.) .AND. (finished .EQV. .FALSE.) ) THEN            
      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_WriteDimensionsWrapper(pRegion,WRITE_DIMENS_MODE_MAYBE)
                
        IF ( moveGrid .EQV. .TRUE. ) THEN        
          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteGridSpeedsWrapper(pRegion)
          
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_NONE ) THEN           
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%myProcid
        END IF ! moveGrid
                       
        CALL RFLU_WriteFlowWrapper(pRegion)
        CALL RFLU_WritePatchCoeffsWrapper(pRegion)
        
#ifdef PLAG
        CALL PLAG_WriteSurfStatsWrapper(pRegion)
#endif  
        
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN          
! TEMPORARY
!          CALL RFLU_PrintFlowInfoWrapper(pRegion)
! END TEMPORARY
          
          IF ( global%verbLevel > VERBOSE_LOW ) THEN 
            CALL RFLU_PrintChangeInfo(pRegion)
          END IF ! global%verbLevel
        END IF ! global%myProcid        
      END DO ! iReg
      
      CALL RFLU_WriteRestartInfo(global)      
    END IF ! doWrite

#ifdef STATS
! ==============================================================================
!   Output statistics
! ==============================================================================

    IF ( (doWrite .EQV. .TRUE.) .AND. (finished .EQV. .FALSE.) .AND. &
         (global%doStat == ACTIVE) ) THEN
      IF (global%myProcid==MASTERPROC .AND. &
          global%verbLevel/=VERBOSE_NONE) THEN
        WRITE(STDOUT,'(A)') SOLVER_NAME,'Saving statistics ...'
      ENDIF
      
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL RFLU_WriteStat(pRegion)
      END DO ! iReg
    END IF ! iReg
#endif
#endif

! ==============================================================================
!   If run finished, update GENX buffers and exit
! ==============================================================================

    IF ( finished .EQV. .TRUE. ) THEN 
#ifdef GENX
       DO iReg = 1,global%nRegionsLocal
         CALL RFLU_PutBoundaryValues(regions(iReg))
       END DO ! iReg
       
       global%timeStamp = global%currentTime       
#endif   
! DEBUG
       IF ( global%myProcid == 0 ) THEN
         CALL PetscGetTime(time,ierr)
         print*,'PETSC TIME 2 : ',time
       END IF ! global%myProcid
! END DEBUG
       EXIT
    END IF ! finished

! ==============================================================================
!   End of timestepping loop
! ==============================================================================

  END DO ! loop over time/iterations

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NK_TimeStepping






  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModNewtonKrylov


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModNewtonKrylov.F90,v $
! Revision 1.16  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2006/08/19 15:39:10  mparmar
! Added use of RFLU_DecideNeedBGradFace and pPatch%spaceOrder
!
! Revision 1.13  2006/04/07 16:04:02  haselbac
! Adapted to changes in bf2c wts computation
!
! Revision 1.12  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.11  2006/04/07 14:49:41  haselbac
! Adapted to changes in f and bf diff modules and new WENO module
!
! Revision 1.10  2006/03/25 21:54:01  haselbac
! Cosmetics only
!
! Revision 1.9  2006/03/22 14:58:28  hdewey2
! Fixed bug by adding additional virtual cell communication after SNESSolve calls.
!
! Revision 1.8  2006/03/09 14:07:42  haselbac
! Now call wrapper routine for F2C weights
!
! Revision 1.7  2006/02/08 21:23:39  hdewey2
! Added disp ramping and changed grid storage
!
! Revision 1.6  2006/01/06 22:12:44  haselbac
! Added call to cell grad wrapper routine
!
! Revision 1.5  2005/10/25 19:38:47  haselbac
! Reordered modules, added IF on forceFlag
!
! Revision 1.4  2005/09/22 17:15:26  hdewey2
! Added dual timestepping functionality to solve implicit time-dependent problems.
!
! Revision 1.3  2005/08/02 18:19:43  hdewey2
! Added actual procedures
!
! Revision 1.2  2005/05/19 18:19:00  haselbac
! Cosmetics only
!
! Revision 1.1  2005/05/16 21:12:45  haselbac
! Initial revision
!
! ******************************************************************************
  







