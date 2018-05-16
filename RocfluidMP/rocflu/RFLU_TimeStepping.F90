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
! Purpose: Integrate the governing equations in time; move/regenerate the grid.
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
!
! $Id: RFLU_TimeStepping.F90,v 1.73 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_TimeStepping(dTimeSystem,dIterSystem,regions)

  USE ModDataTypes
  USE ModError  
  USE ModParameters
  USE ModGlobal, ONLY: t_global  
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input    
  USE ModMPI
  
  USE RFLU_ModDimensions
  USE RFLU_ModForcesMoments
  USE RFLU_ModThrustSpecImpulse
  USE RFLU_ModGeometry
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModProbes
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModWeights

  USE RFLU_ModTimeZoom, ONLY: RFLU_UnZoomGridSpeeds,&
                              RFLU_ZoomGridSpeeds
  USE RFLU_ModBoundXvUtils, ONLY: RFLU_BXV_WriteVarsWrapper
  
#ifdef PLAG
  USE PLAG_ModDimensions, ONLY: PLAG_CalcNPclsGlobal, &
                                PLAG_PrintNPclsGlobal, &
                                PLAG_RFLU_WriteDimensions
  USE PLAG_ModSurfStats
#endif  
  
  USE ModInterfaces, ONLY: IntegrateSourceTermsMP, &
                           RFLU_ComputeGridSpeeds, &
                           RFLU_ComputeIntegralValues, & 
                           RFLU_CheckGridSpeeds, &
                           RFLU_DecidePrint, &
                           RFLU_DecideNeedBGradFace, & 
                           RFLU_DecideWrite, & 
                           RFLU_ExplicitMultiStage, &
                           RFLU_GetDeformationWrapper, & 
                           RFLU_MinimumTimeStep, &
                           RFLU_MoveGridWrapper, &
                           RFLU_PrintChangeInfo, &  
                           RFLU_PrintGridInfo, &
                           RFLU_PrintFlowInfoWrapper, &
                           RFLU_PrintWriteConvergence, &
                           RFLU_ResidualNorm, &
                           RFLU_TimeStepInviscid, &
                           RFLU_TimeStepViscous, &
                           RFLU_WriteRestartInfo, &
                           RFLU_WriteStatsFileOLES, &
                           RungeKuttaMP, &
                           WriteProbe, &
                           WriteTotalMass

#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_PutBoundaryValues
#endif
#ifdef STATS
  USE ModInterfaces, ONLY: RFLU_WriteStat
  USE ModStatsRoutines, ONLY: GetStatistics
#endif

  IMPLICIT NONE

#ifdef GENX
#include "comf90.h"
#endif

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: dIterSystem
  REAL(RFREAL) :: dTimeSystem
  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iPatch,iReg,stophandle,dtlimcount
  LOGICAL :: doPrint,doProbe,doWrite,finished,ftermNew,moveGrid,residFterm, &
             stopFileExists
  TYPE(t_global), POINTER :: global
  TYPE(t_mixt_input), POINTER :: pMixtInput   
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_TimeStepping.F90,v $ $Revision: 1.73 $'

! ******************************************************************************
! Set pointers and variables
! ****************************************************************************** 

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_TimeStepping',&
  'RFLU_TimeStepping.F90')

  finished = .FALSE. ! run not finished yet

#ifdef GENX
  global%timeSinceRestart = 0.0_RFREAL
#endif

! ==============================================================================
! No multigrid here
! ==============================================================================

  ftermNew   = .FALSE. ! no new forcing term
  residFterm = .FALSE. ! do not add forcing term to residual

! ==============================================================================
! Determine whether have moving grids
! ==============================================================================
  dtlimcount = 0
  moveGrid = .FALSE.

  DO iReg = 1,global%nRegionsLocal
    IF ( regions(iReg)%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      moveGrid = .TRUE.
    END IF ! regions
  END DO ! iReg

! ******************************************************************************
! Loop over iterations/time steps
! ****************************************************************************** 

  DO

! ==============================================================================
!   Check for non-zero system time or iteration step. This is needed because 
!   of changed restart mechanism. If restarting code when final time was 
!   reached, can get zero time step which gives indeterminate or infinite 
!   grid speeds for moving grid computations. 
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN 
      IF ( dTimeSystem == 0.0_RFREAL ) THEN
        global%warnCounter = global%warnCounter + 1       
      
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_NONE ) THEN 
          WRITE(STDOUT,'(A,2(1X,A))') SOLVER_NAME,'*** WARNING *** ', & 
                'Nothing to be done. Returning to calling procedure.'
        END IF ! global
        
        EXIT
      END IF ! dTimeSystem
    ELSE
      IF ( dIterSystem == 0 ) THEN
        global%warnCounter = global%warnCounter + 1       
      
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_NONE ) THEN 
          WRITE(STDOUT,'(A,2(1X,A))') SOLVER_NAME,'*** WARNING *** ', & 
                'Nothing to be done. Returning to calling procedure.'
        END IF ! global
        
        EXIT
      END IF ! dIterSystem
    END IF ! global%flowType

! ==============================================================================
!   Update iteration counter for steady flow. NOTE does not have to be done
!   here, but is more consistent, otherwise can get output for iteration 0 
!   from RFLU_DecidePrint, and again from iteration 1 because of how MOD 
!   function works. If update iteration counters here, only get output once.
! ==============================================================================

    IF ( global%flowType == FLOW_STEADY ) THEN
      global%currentIter      = global%currentIter      + 1
      global%iterSinceRestart = global%iterSinceRestart + 1
    END IF ! global%flowType

! ==============================================================================
!   Compute time step. For unsteady flow, compute minimum time step and make 
!   sure do not run over maximum time.
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal  
       pRegion => regions(iReg)
       
       IF ( pRegion%mixtInput%flowModel == FLOW_EULER ) THEN      
          CALL RFLU_TimeStepInviscid(pRegion)
       ELSE 
          CALL RFLU_TimeStepViscous(pRegion)
       END IF ! pRegion
    END DO ! iReg
    
    IF ( global%flowType == FLOW_UNSTEADY ) THEN         
       CALL RFLU_MinimumTimeStep(regions)
#ifdef GENX
       IF ( global%dtMin < global%dtImposed .AND. &
            global%dtMin < global%dtMinLimit ) THEN
          dtlimcount = dtlimcount + 1
          IF ( dtlimcount > 2) THEN
             dtlimcount = 0
             stophandle   = COM_get_function_handle('Rocman.interrupt')
             IF(stophandle <= 0) THEN
                WRITE(*,*) 'Could not get Rocman.stop function handle.'
               ELSE
                  CALL COM_call_function(stophandle,2,3,'Rocflu dt < limit, requesting remesh')
               ENDIF
            ENDIF
         ENDIF
#endif
         
         IF ( global%timeSinceRestart + global%dtMin > dTimeSystem ) THEN 
            global%dtMin = dTimeSystem - global%timeSinceRestart
            finished     = .TRUE.                     
         END IF ! global%timeSinceRestart            
      END IF ! global
      
! ==============================================================================
!   Move or generate new grid
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN    
      IF ( moveGrid .EQV. .TRUE. ) THEN     
        CALL RFLU_GetDeformationWrapper(regions)        
        CALL RFLU_MoveGridWrapper(regions)
            
        DO iReg=1,global%nRegionsLocal
          pRegion => regions(iReg)
                   
          CALL RFLU_BuildGeometry(pRegion)          
          CALL RFLU_ComputeGridSpeeds(pRegion)

          IF(global%zoomFactor > 1) THEN
             CALL RFLU_UnZoomGridSpeeds(pRegion)
          ENDIF
          
          IF ( global%checkLevel == CHECK_HIGH ) THEN           
            CALL RFLU_CheckGridSpeeds(pRegion)               
          END IF ! global%checkLevel
        END DO ! iReg
      END IF ! moveGrid     
    END IF ! global%flowType  

! ==============================================================================
!   Recompute weights
! ==============================================================================
#ifndef PROPONLY

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
              CALL RFLU_ComputeWtsBF2CWrapper(pRegion,pPatch,pPatch%spaceOrder)
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
!   Compute new solution
! ==============================================================================

    global%forceX  = 0.0_RFREAL
    global%forceY  = 0.0_RFREAL
    global%forceZ  = 0.0_RFREAL
    global%massIn  = 0.0_RFREAL
    global%massOut = 0.0_RFREAL  

    IF ( global%flowType == FLOW_UNSTEADY ) THEN    
      CALL RungeKuttaMP(regions)      
! TEMPORARY - At present, do not use operator-split integration of chemistry
!             source terms based on modifications by Luca. This means that 
!             call source term residual function directly in SourceTermsMP.F90
!      CALL IntegrateSourceTermsMP(regions)  
! END TEMPORARY
    ELSE
      CALL RFLU_ExplicitMultiStage(regions)
    END IF ! global%flowType
 
#ifdef STATS
! ==============================================================================
!   Get statistics
! ==============================================================================

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

! ENDIF PROPONLY
#endif

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      global%currentTime      = global%currentTime      + global%dtMin      
      global%timeSinceRestart = global%timeSinceRestart + global%dtMin
      
      global%timeSincePrint = global%timeSincePrint + global%dtMin
      global%timeSinceWrite = global%timeSinceWrite + global%dtMin
      global%timeSinceProbe = global%timeSinceProbe + global%dtMin      
      
      global%iterSinceRestart = global%iterSinceRestart + 1
    END IF ! global%flowType

! ==============================================================================
!   Decide whether to print/write convergence, data, and probes
! ==============================================================================

    doPrint = RFLU_DecidePrint(global)
    doWrite = RFLU_DecideWrite(global)
    doProbe = RFLU_DecideWriteProbes(global)

! ==============================================================================
!   Check for stop file
! ==============================================================================

    INQUIRE(FILE="STOP",EXIST=stopFileExists)
    IF ( stopFileExists .EQV. .TRUE. ) THEN 
      IF ( global%myProcid  == MASTERPROC .AND. &
           global%verbLevel /= VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A)')      SOLVER_NAME
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Stop file detected!'
        WRITE(STDOUT,'(A)')      SOLVER_NAME        
      END IF ! global

      finished = .TRUE.
    END IF ! stopFileExists
  
! ==============================================================================
!   Check for end of time stepping
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
      
#ifndef PROPONLY
      DO iReg = 1,global%nRegionsLocal
        IF ( regions(iReg)%mixtInput%spaceDiscr == DISCR_OPT_LES ) THEN 
          CALL RFLU_WriteStatsFileOLES(global)
        END IF ! mixtInput
      END DO ! iReg
#endif
    END IF ! doPrint

#ifndef PROPONLY
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

        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_TSI_ComputeGlobalThrustSI(pRegion)     
        END DO ! iReg

        pRegion => regions(1)

        IF ( pRegion%global%myProcid == MASTERPROC ) THEN 
          CALL RFLU_PrintGlobalForcesMoments(pRegion)
          CALL RFLU_WriteGlobalForcesMoments(pRegion)
          CALL RFLU_TSI_PrintGlobalVals(pRegion)
          CALL RFLU_TSI_WriteGlobalVals(pRegion)
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

#ifndef GENX
! ==============================================================================
!   Store flow field (and grid if moving). Write restart info file after 
!   flow (and grid) files so that should those not be written due to reaching
!   the time limit, the restart file will not contain the time level of the 
!   incomplete flow (and grid) files.
! ==============================================================================

    IF ( (doWrite .EQV. .TRUE.) .AND. (finished .EQV. .FALSE.) ) THEN            
#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL PLAG_CalcNPclsGlobal(regions)

        IF ( global%myProcid == MASTERPROC ) THEN 
          pRegionSerial => regions(0)

          CALL PLAG_RFLU_WriteDimensions(pRegionSerial)
        END IF ! global%myProcid
      END IF ! global%plagUsed
#endif

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
        CALL RFLU_BXV_WriteVarsWrapper(pRegion)

        IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
          CALL RFLU_WritePatchCoeffsWrapper(pRegion)
        END IF ! global%patchCoeffFlag
        
#ifdef PLAG
        IF ( global%plagUsed .EQV. .TRUE. ) THEN 
          CALL PLAG_WriteSurfStatsWrapper(pRegion)
        END IF ! global%plagUsed
#endif  
        
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN          
          CALL RFLU_PrintFlowInfoWrapper(pRegion)
          
          IF ( global%verbLevel > VERBOSE_LOW ) THEN 
            CALL RFLU_PrintChangeInfo(pRegion)
          END IF ! global%verbLevel

#ifdef PLAG
          IF ( global%plagUsed .EQV. .TRUE. ) THEN 
            CALL PLAG_PrintNPclsGlobal(pRegion)
          END IF ! global%plagUsed
#endif
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
      EXIT
    END IF ! finished

  END DO ! loop over time/iterations

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TimeStepping

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_TimeStepping.F90,v $
! Revision 1.73  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.72  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.71  2007/04/14 14:25:56  mtcampbe
! Updated for TZ
!
! Revision 1.70  2007/03/31 23:54:34  haselbac
! Bug fix: Writing of PLAG dims should only be called by serial region on master proc
!
! Revision 1.69  2007/03/27 00:41:17  haselbac
! Added calls to calculate, write, and print particle dimensions
!
! Revision 1.68  2007/02/18 03:17:56  mtcampbe
! Added proponly for Rocstar proponly runs
!
! Revision 1.67  2006/10/20 21:33:07  mparmar
! Added calls to write thrust/specific impulse and again added code for NSCBC implementation
!
! Revision 1.66  2006/09/12 14:58:44  mtcampbe
! Moved include of Roccom after IMPLICIT NONE
!
! Revision 1.65  2006/09/11 15:43:11  mtcampbe
! Added Rocman interrupt call to support automatic remeshing.
!
! Revision 1.64  2006/08/19 15:40:04  mparmar
! Changed because of NSCBC implementation
!
! Revision 1.63  2006/04/07 16:04:03  haselbac
! Adapted to changes in bf2c wts computation
!
! Revision 1.62  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.61  2006/04/07 14:55:06  haselbac
! Adapted to changes in bface stencil routine
!
! Revision 1.60  2006/03/09 14:10:31  haselbac
! Now call wrapper routine for F2C weights
!
! Revision 1.59  2006/01/06 22:16:13  haselbac
! Adapted to name changes, removed commented-out call to ExplMS routine
!
! Revision 1.58  2005/11/10 16:51:29  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.57  2005/10/25 19:39:23  haselbac
! Added IF on forceFlag
!
! Revision 1.56  2005/10/05 14:20:35  haselbac
! Added call to recompute bface wts for moving grids
!
! Revision 1.55  2005/08/09 00:59:56  haselbac
! Enclosed writing of patch coeffs within IF (patchCoeffFlag)
!
! Revision 1.54  2005/06/06 14:23:35  haselbac
! Adapted to Lucas changes
!
! Revision 1.53  2005/05/16 20:44:55  haselbac
! Now compute time step also for steady flow, call RFLU_ExplicitMultiStage
!
! Revision 1.52  2005/04/29 12:50:00  haselbac
! Added USE RFLU_ModProbes, removed interfaces for probe routines
!
! Revision 1.51  2005/04/20 14:44:04  haselbac
! Removed CHECK_UNIFLOW code section
!
! Revision 1.50  2005/04/15 15:07:26  haselbac
! Now use RFLU_PrintWriteConvergence, fixed bug in writing forces
!
! Revision 1.49  2004/12/28 20:28:20  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.48  2004/12/21 15:05:11  fnajjar
! Included calls for PLAG surface statistics
!
! Revision 1.47  2004/11/29 17:17:29  wasistho
! use ModInterfacesStatistics
!
! Revision 1.46  2004/11/03 17:05:53  haselbac
! Removed HACK_PERIODIC ifdef
!
! Revision 1.45  2004/10/19 19:29:36  haselbac
! Added writing of grid speeds, changed time step interfaces, cosmetics
!
! Revision 1.44  2004/07/06 15:14:53  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.43  2004/06/16 20:01:15  haselbac
! Added forces and moments stuff, cosmetics
!
! Revision 1.42  2004/04/14 02:10:07  haselbac
! Removed call to DescaleGridSpeeds
!
! Revision 1.41  2004/04/01 21:30:35  haselbac
! Added IntegrateSourceTermsMP, cosmetic changes
!
! Revision 1.40  2004/03/17 04:28:25  haselbac
! Adapted call to RFLU_WriteDimensionsWrapper
!
! Revision 1.39  2004/03/11 16:34:32  fnajjar
! ACH: Changed call to RFLU_WriteDimWrapper bcos of Rocpart
!
! Revision 1.38  2004/01/31 03:59:03  haselbac
! Improved updating of iteration counters
!
! Revision 1.37  2004/01/29 22:59:24  haselbac
! Added setting of timeSince* vars for improved useability of RFLU_Decide* 
! funcs
!
! Revision 1.36  2003/12/04 03:30:09  haselbac
! Adapted recomputation of weights
!
! Revision 1.35  2003/11/25 21:04:46  haselbac
! Added call to RFLU_PrintFlowInfoWrapper, cosmetic changes
!
! Revision 1.34  2003/11/04 01:35:27  haselbac
! Bug fix: added init for timeSinceRestart for GENX
!
! Revision 1.33  2003/10/29 21:39:57  haselbac
! Bug fix and clean-up: Writing of data to screen and files
!
! Revision 1.32  2003/10/15 02:44:21  haselbac
! Removed unnecessary initialization of dtMin
!
! Revision 1.31  2003/10/03 20:47:21  haselbac
! Changed from RungeKutta to RungeKuttaMP
!
! Revision 1.30  2003/08/13 20:28:49  haselbac
! Fixed bug with writing probe data within GENx
!
! Revision 1.29  2003/07/22 02:11:40  haselbac
! Added global%warnCounter
!
! Revision 1.28  2003/07/09 22:37:48  haselbac
! Added check to avoid wrong restart, find probes on moving grids
!
! Revision 1.27  2003/06/20 22:36:23  haselbac
! Added call to RFLU_WriteRestartInfo, fixed bug in steady conv check
!
! Revision 1.26  2003/04/24 15:43:35  haselbac
! Adapted interface to RFLU_PutBoundaryValues
!
! Revision 1.25  2003/04/07 14:27:44  haselbac
! Added interval writing of probe info
!
! Revision 1.24  2003/03/31 16:18:17  haselbac
! Replaced MoveGrid call by call to wrapper
!
! Revision 1.23  2003/03/15 18:59:56  haselbac
! Added calls to DescaleGridSpeeds and GetDeformationWrapper
!
! Revision 1.22  2003/02/25 21:47:39  haselbac
! Added call to DescaleGridSpeeds
!
! Revision 1.21  2003/01/28 14:52:08  haselbac
! Changes to be consistent with rewrite of RFLU_InitFlowSolver
!
! Revision 1.20  2002/12/20 23:21:33  haselbac
! Fixed output bug: no output for verbosity=0
!
! Revision 1.19  2002/11/15 21:26:49  haselbac
! Added RFLU_ComputeIntegralValues
!
! Revision 1.18  2002/11/08 21:36:25  haselbac
! Fixed bug in grid-speed comp, added RFLU_CheckGridSpeeds and WriteTotalMass
!
! Revision 1.17  2002/11/02 02:04:49  wasistho
! Added TURB statistics
!
! Revision 1.16  2002/10/27 19:19:31  haselbac
! Logical modifications of last part, cosmetic redesign
!
! Revision 1.15  2002/10/16 21:18:50  haselbac
! Added call to RFLU_NewGrid, some rearrangement of calls
!
! Revision 1.14  2002/10/12 15:00:44  haselbac
! Added statements for time check for GENX runs
!
! Revision 1.13  2002/10/05 19:35:19  haselbac
! Call wrapper routines for output, write probe data
!
! Revision 1.12  2002/09/09 15:51:56  haselbac
! global and mixtInput under regions, adapated calls to RFLU_CheckCalcGrad, 
! OLES routines removed, added viscous calls
!
! Revision 1.11  2002/07/25 14:23:40  haselbac
! Added OLES calls, gradient check, and HACK_PERIODIC segment
!
! Revision 1.10  2002/06/18 00:26:20  wasistho
! Added prefix SOLVER NAME to satistics STDOutput
!
! Revision 1.9  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.8  2002/06/14 22:26:08  wasistho
! update statistics
!
! Revision 1.7  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.6  2002/06/14 20:19:46  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal
!
! Revision 1.5  2002/06/10 21:35:50  haselbac
! Changed order of convergence and file writing, changed flag to CHECK_UNIFLOW
!
! Revision 1.4  2002/06/05 18:59:12  haselbac
! Miscellaneous small functionality changes
!
! Revision 1.3  2002/05/28 14:02:04  haselbac
! Activated unsteady routines
!
! Revision 1.2  2002/05/04 17:13:15  haselbac
! Many changes to enable flow solution
!
! Revision 1.1  2002/04/11 18:57:06  haselbac
! Initial revision, commented out RFLO stuff
!
! ******************************************************************************







