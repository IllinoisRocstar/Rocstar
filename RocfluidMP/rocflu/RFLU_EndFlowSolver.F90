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
! Purpose: Shut down Rocflu-MP.
!
! Description: None.
!
! Input: 
!   levels      Data associated with levels
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_EndFlowSolver.F90,v 1.62 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_EndFlowSolver(levels)

  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input   
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModMPI  
  
  USE RFLU_ModBFaceGradAccessList
  USE RFLU_ModBoundLists 
  USE RFLU_ModBoundXvUtils
  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
  USE RFLU_ModCellMapping 
  USE RFLU_ModCommLists 
  USE RFLU_ModDimensions
  USE RFLU_ModEdgeList   
  USE RFLU_ModFaceList  
  USE RFLU_ModForcesMoments
  USE RFLU_ModGeometry      
  USE RFLU_ModMPI
  USE RFLU_ModOLES
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModProbes
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds  
  USE RFLU_ModStencilsBFaces
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsFaces
  USE RFLU_ModWeights
  
#ifdef PETSC
  USE RFLU_ModPETScAdmin
  USE RFLU_ModPETScNewtonKrylov
#endif
  
#ifdef PLAG
  USE PLAG_ModDimensions, ONLY: PLAG_CalcNPclsGlobal, &
                                PLAG_PrintNPclsGlobal, &
                                PLAG_RFLU_WriteDimensions
  USE PLAG_ModSurfStats
#endif  
  
  USE ModInterfaces, ONLY: RFLU_DeallocateMemoryWrapper, &
                           RFLU_DecideNeedBGradFace, &
                           RFLU_DestroyGrid, & 
                           RFLU_PrintFlowInfoWrapper, & 
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintWarnInfo, &  
                           RFLU_WriteRestartInfo
    
#ifdef STATS
  USE ModInterfaces, ONLY: RFLU_WriteStat
#endif

  IMPLICIT NONE

! ******************************************************************************
! Arguments
! ******************************************************************************

  TYPE(t_level), POINTER :: levels(:)

! ******************************************************************************
! Locals
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: moveGrid  
  INTEGER :: errorFlag,iPatch,iProbe,iReg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput 
  TYPE(t_patch), POINTER :: pPatch   
  TYPE(t_region), POINTER :: pRegion,pRegionSerial  

! ******************************************************************************
! Start, initialize some variables
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_EndFlowSolver.F90,v $ $Revision: 1.62 $'

  moveGrid = .FALSE. 

! ******************************************************************************
! Set global pointer and initialize global type, register function
! ****************************************************************************** 

  global => levels(1)%regions(1)%global

  CALL RegisterFunction(global,'RFLU_EndFlowSolver',&
  'RFLU_EndFlowSolver.F90')

! ******************************************************************************
! Determine whether have moving grid
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg)
  
    IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      moveGrid = .TRUE.
      EXIT
    END IF ! regions
  END DO ! iReg

! ******************************************************************************
! Write grid file (if necessary) and flow file. Write restart info file after 
! flow (and grid) file so that incomplete flow (and grid) files due to
! exceeding time limit do not show up as iteration or time stamp in restart 
! info file.
! ******************************************************************************

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    CALL PLAG_CalcNPclsGlobal(levels(1)%regions)

    IF ( global%myProcid == MASTERPROC ) THEN
      pRegionSerial => levels(1)%regions(0)

      CALL PLAG_RFLU_WriteDimensions(pRegionSerial)
    END IF ! global%myProcid
  END IF ! global%plagUsed
#endif

#ifndef GENX
  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg) ! single-level grids for now
    
    CALL RFLU_WriteDimensionsWrapper(pRegion,WRITE_DIMENS_MODE_MAYBE)    

    IF ( moveGrid .EQV. .TRUE. ) THEN
      CALL RFLU_WriteGridWrapper(pRegion)
      CALL RFLU_WriteGridSpeedsWrapper(pRegion)
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
    
  END DO ! iReg

  CALL RFLU_WriteRestartInfo(global)
  
#ifdef STATS
    IF ( global%doStat == ACTIVE ) THEN
      DO iReg = 1,global%nRegionsLocal
        pRegion => levels(1)%regions(iReg)
        CALL RFLU_WriteStat(pRegion)
      END DO ! iReg
    END IF ! global%doStat
#endif
#endif

! ******************************************************************************
! Print (if necessary, grid, and) flow information 
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    IF ( moveGrid .EQV. .TRUE. ) THEN
      DO iReg = 1,global%nRegionsLocal
        pRegion => levels(1)%regions(iReg) ! single-level grids for now
        CALL RFLU_PrintGridInfo(pRegion)
      END DO ! iReg     
    END IF ! moveGrid       
       
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) ! single-level grids for now
      CALL RFLU_PrintFlowInfoWrapper(pRegion)

#ifdef PLAG                     
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL PLAG_PrintNPclsGlobal(pRegion)
      END IF ! global%plagUsed
#endif
    END DO ! iReg    
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory. NOTE must be done before calling deallocation wrapper. 
! NOTE weights must be deallocated before stencils.
! ******************************************************************************

! ==============================================================================
! Weights for cell and face gradients
! ==============================================================================

  DO iReg = 1,global%nRegionsLocal 
    pRegion => levels(1)%regions(iReg)
    pMixtInput => pRegion%mixtInput   

    IF ( pMixtInput%spaceOrder > 1 ) THEN 
      CALL RFLU_DestroyWtsC2CWrapper(pRegion)
    END IF ! pMixtInput%spaceOrder      

    IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN      
      CALL RFLU_DestroyWtsF2CWrapper(pRegion)
    END IF ! pMixtInput%flowModel 
      
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        CALL RFLU_DestroyWtsBF2CWrapper(pRegion,pPatch)      
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
  END DO ! iReg   
  
! ==============================================================================
! Weights for optimal LES approach
! ==============================================================================

  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg) 
    pMixtInput => pRegion%mixtInput 
         
    IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN 
      CALL RFLU_DestroyStencilsWeightsOLES(pRegion)
    END IF ! pMixtInput
  END DO ! iReg
 
! ******************************************************************************
! Deallocate memory for stencils
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => levels(1)%regions(iReg)
    pMixtInput => pRegion%mixtInput 
    
    IF ( pMixtInput%spaceOrder > 1 ) THEN 
      CALL RFLU_DestroyC2CStencilWrapper(pRegion)
      CALL RFLU_DestroyListCC2CStencil(pRegion)  
    END IF ! pMixtInput%spaceOrder
    
    IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN      
      CALL RFLU_DestroyF2CStencilWrapper(pRegion) 
      CALL RFLU_DestroyListCF2CStencil(pRegion)
    END IF ! pMixtInput%flowModel

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)
    
      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        CALL RFLU_DestroyBF2CStencilWrapper(pRegion,pPatch)
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch

  END DO ! iReg

! ******************************************************************************
! Destroy boundary-face gradient access list for viscous flows
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => levels(1)%regions(iReg) 
    
    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN          
      CALL RFLU_DestroyBFaceGradAccessList(pRegion)
    END IF ! pRegion%mixtInput
  END DO ! iReg  

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg) ! single-level grids for now
    pGrid   => pRegion%grid
        
    CALL RFLU_DestroyPatchCoeffs(pRegion)

    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_DestroyVarsCv(pRegion)
      CALL RFLU_BXV_DestroyVarsDv(pRegion)
      CALL RFLU_BXV_DestroyVarsTStep(pRegion)
    END IF ! RFLU_NSCBC_DecideHaveNSCBC(pRegion)   
 
    IF ( global%forceFlag .EQV. .TRUE. ) THEN 
      CALL RFLU_DestroyForcesMoments(pRegion)   
      CALL RFLU_DestroyGlobalThrustFlags(pRegion)   
    END IF ! global%forceFlag

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN 
      CALL PLAG_DestroySurfStats(pRegion)
    END IF ! global%plagUsed
#endif    
            
    IF ( pGrid%nBorders > 0 ) THEN 
      CALL RFLU_MPI_DestroyBuffersWrapper(pRegion)

#ifdef PLAG
      IF ( pRegion%global%plagUsed .EQV. .TRUE. ) THEN 
        CALL RFLU_MPI_DestroyBufferIPclSend(pRegion)
      END IF ! pRegion%global%plagUsed
#endif

      CALL RFLU_COMM_DestroyCommLists(pRegion)
      CALL RFLU_COMM_DestroyBorders(pRegion)
    END IF ! pGrid%nBorders 
        
    CALL RFLU_DeallocateMemoryWrapper(pRegion)

    CALL RFLU_DestroyCellMapping(pRegion)
    CALL RFLU_DestroyFaceList(pRegion)
    
#ifdef PLAG
    IF ( (pRegion%global%plagUsed .EQV. .TRUE.) .AND. & 
         (pRegion%grid%nFacesAV > 0) ) THEN
      CALL RFLU_DestroyAVFace2BorderList(pRegion)      
      CALL RFLU_DestroyAVFace2PatchList(pRegion)      
    END IF ! pRegion%global%plagUsed
#endif    
    
    IF ( pRegion%mixtInput%movegrid .EQV. .TRUE. ) THEN
      CALL RFLU_DestroyEdgeList(pRegion)
    END IF ! pRegion 
    
    CALL RFLU_DestroyGeometry(pRegion)
    CALL RFLU_DestroyGrid(pRegion)           
  END DO ! iReg

! ******************************************************************************
! Close convergence, mass, and probe files
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN
    CLOSE(IF_CONVER,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
    END IF ! global%error
  END IF ! global%myProcid

  IF ( (global%myProcid == MASTERPROC) .AND. (moveGrid .EQV. .TRUE.) ) THEN
    CLOSE(IF_MASS,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
    END IF ! global%error
  END IF ! global%myProcid  

  IF ( global%nProbes > 0 ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) ! single-level grids for now 
         
      CALL RFLU_CloseProbeFiles(pRegion)
    END DO ! iReg
  END IF ! global%nProbes

! ******************************************************************************
! Print info about warnings
! ******************************************************************************

  CALL RFLU_PrintWarnInfo(global)

! ******************************************************************************
! Deallocate PETSc memory & finalize PETSc
! ******************************************************************************

#ifdef PETSC
  IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
    CALL RFLU_PETSC_DestroyVectors(pRegion)
    CALL RFLU_PETSC_DestroyJacobian(pRegion)
    CALL RFLU_PETSC_Finalize(global)
  END IF ! global%solverType
#endif

#ifndef GENX
! ******************************************************************************
! Finalize Rocprof
! ****************************************************************************** 

#ifdef ROCPROF
  CALL Rocprof_Finalize(1)
#endif

! ******************************************************************************
! Finalize MPI
! ****************************************************************************** 

  CALL MPI_Finalize(errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel /= VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finalization done.'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Program finished.'
    WRITE(STDOUT,'(A)') SOLVER_NAME         
  END IF ! global%myProcid 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EndFlowSolver

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_EndFlowSolver.F90,v $
! Revision 1.62  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.61  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.60  2007/03/31 23:53:39  haselbac
! Added calls to determine, write, and print nPclsGlobal
!
! Revision 1.59  2006/10/20 21:32:08  mparmar
! Added call to RFLU_DestroyGlobalThrustFlags
!
! Revision 1.58  2006/08/19 15:46:41  mparmar
! Changed logic for NSCBC and added calls to deallocate patch arrays
!
! Revision 1.57  2006/08/18 14:04:02  haselbac
! Added call to destroy AVFace2Patch list
!
! Revision 1.56  2006/04/07 16:04:03  haselbac
! Adapted to changes in bf2c wts computation
!
! Revision 1.55  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.54  2006/04/07 14:53:11  haselbac
! Adapted to changes in bface stencil routines
!
! Revision 1.53  2006/03/09 14:09:49  haselbac
! Now call wrapper routines for stencils
!
! Revision 1.52  2006/01/06 22:15:07  haselbac
! Adapted to name changes
!
! Revision 1.51  2005/11/10 16:51:29  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.50  2005/10/27 19:20:06  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.49  2005/10/25 19:39:23  haselbac
! Added IF on forceFlag
!
! Revision 1.48  2005/10/05 14:18:41  haselbac
! Adapted to changes in stencil modules, added call to destroy bface wts
!
! Revision 1.47  2005/09/14 15:59:33  haselbac
! Minor clean-up
!
! Revision 1.46  2005/09/13 20:44:06  mtcampbe
! Added Rocprof finalization
!
! Revision 1.45  2005/08/09 00:59:42  haselbac
! Enclosed writing of patch coeffs within IF (patchCoeffFlag)
!
! Revision 1.44  2005/08/03 18:28:11  hdewey2
! Enclosed PETSc deallocation calls within IF
!
! Revision 1.43  2005/08/02 18:24:32  hdewey2
! Added PETSc support
!
! Revision 1.42  2005/05/18 22:12:55  fnajjar
! ACH: Added destruction of iPclSend buffers, now use nFacesAV
!
! Revision 1.41  2005/04/29 23:03:22  haselbac
! Added destruction of avf2b list
!
! Revision 1.40  2005/04/29 12:48:23  haselbac
! Updated closing of probe files to changes in probe handling
!
! Revision 1.39  2005/04/15 16:31:18  haselbac
! Removed calls to XyzEdge2RegionDegrList routines
!
! Revision 1.38  2005/04/15 15:07:15  haselbac
! Converted to MPI
!
! Revision 1.37  2005/01/18 15:18:19  haselbac
! Commented out COMM calls for now
!
! Revision 1.36  2005/01/14 21:33:56  haselbac
! Added calls to destroy comm lists and borders
!
! Revision 1.35  2005/01/03 16:14:57  haselbac
! Added call to destroy bface stencils
!
! Revision 1.34  2004/12/21 15:05:11  fnajjar
! Included calls for PLAG surface statistics
!
! Revision 1.33  2004/10/19 19:29:15  haselbac
! Adapted to GENX changes
!
! Revision 1.32  2004/07/06 15:14:50  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.31  2004/06/16 20:01:08  haselbac
! Added writing of patch coeffs, destruction of memory
!
! Revision 1.30  2004/03/17 04:28:26  haselbac
! Adapted call to RFLU_WriteDimensionsWrapper
!
! Revision 1.29  2004/03/11 16:33:22  fnajjar
! ACH: Moved call to RFLU_WriteDimWrapper outside if bcos of Rocpart
!
! Revision 1.28  2004/03/08 22:01:43  fnajjar
! ACH: Changed call to RFLU_WriteDimensionsWrapper so PLAG data also written
!
! Revision 1.27  2004/01/29 22:59:19  haselbac
! Removed hardcoded error computation for supersonic vortex
!
! Revision 1.26  2003/12/04 03:30:01  haselbac
! Added destruction calls for gradients, cleaned up
!
! Revision 1.25  2003/11/25 21:04:40  haselbac
! Added call to RFLU_PrintFlowInfoWrapper, cosmetic changes
!
! Revision 1.24  2003/11/03 03:51:17  haselbac
! Added call to destroy boundary-face gradient access list
!
! Revision 1.23  2003/08/07 15:33:30  haselbac
! Added call to RFLU_PrintWarnInfo
!
! Revision 1.22  2003/07/22 02:07:00  haselbac
! Added writing out of warnings
!
! Revision 1.21  2003/06/20 22:35:26  haselbac
! Added call to RFLU_WriteRestartInfo
!
! Revision 1.20  2003/01/28 14:36:53  haselbac
! Added deallocation, merged printing of grid and flow info, cosmetics
!
! Revision 1.19  2002/12/20 23:20:06  haselbac
! Fixed output bug: no output for verbosity=0
!
! Revision 1.18  2002/11/08 21:31:32  haselbac
! Added closing of total-mass file
!
! Revision 1.17  2002/11/02 02:03:20  wasistho
! Added TURB statistics
!
! Revision 1.16  2002/10/27 19:12:09  haselbac
! Added writing of grid for moving grid calcs
!
! Revision 1.15  2002/10/19 16:12:16  haselbac
! Cosmetic changes to output
!
! Revision 1.14  2002/10/16 21:16:06  haselbac
! Added writing of header when running
!
! Revision 1.13  2002/10/12 14:57:37  haselbac
! Enclosed RFLU_WriteFlowWrapper between ifndef GENX
!
! Revision 1.12  2002/10/08 15:49:29  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.11  2002/10/05 19:20:55  haselbac
! GENX integration, close probe file, use flow wrapper routine
!
! Revision 1.10  2002/09/09 15:49:26  haselbac
! global and mixtInput now under region
!
! Revision 1.9  2002/07/25 14:25:07  haselbac
! No longer called finalize for CHARM=1
!
! Revision 1.8  2002/06/27 15:27:11  haselbac
! Change name if CHARM defined, comment out destruction for now (crashes with CHARM)
!
! Revision 1.7  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.6  2002/06/14 22:26:08  wasistho
! update statistics
!
! Revision 1.5  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.4  2002/06/14 20:20:43  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal, added destroy flag
!
! Revision 1.3  2002/05/04 17:08:13  haselbac
! Close convergence file and print flow information
!
! Revision 1.2  2002/04/11 19:03:22  haselbac
! Added calls and cosmetic changes
!
! Revision 1.1  2002/03/14 19:11:00  haselbac
! Initial revision
!
! ******************************************************************************







