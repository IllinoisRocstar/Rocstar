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
! Purpose: Wrapper for postprocessing results with merging the regions.
!
! Description: None.
!
! Input: 
!   levels      Level data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_MergePostProcessRegions.F90,v 1.29 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_MergePostProcessRegions(levels)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModBndPatch, ONLY: t_patch  
  USE ModParameters
  
  USE RFLU_ModAllocateMemory
  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModDeallocateMemory
  USE RFLU_ModDimensions
  USE RFLU_ModENSIGHT
  USE RFLU_ModExtractFlowData
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModInterpolation
  USE RFLU_ModMergeRegions
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModPlottingVars
  USE RFLU_ModProbes
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModRenumberings
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsUtils
  USE RFLU_ModStencilsVert
  USE RFLU_ModVertexLists   
  USE RFLU_ModWeights  

#ifndef NO_TECPLOT  
  USE RFLU_ModTECPLOT
#endif

#ifdef PLAG
  USE PLAG_ModDimensions, ONLY: PLAG_PrintDimensions, &
                                PLAG_SetDimensions
#endif

  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, &
                           RFLU_AllocMemVertWrapper, & 
                           RFLU_ComputeExactFlowError, &
                           RFLU_ComputeExactFlowProbeError, & 
                           RFLU_CreateGrid, & 
                           RFLU_DeallocMemSolWrapper, &
                           RFLU_DeallocMemVertWrapper, & 
                           RFLU_DecideBuildGeometry, &
                           RFLU_DestroyGrid, &  
                           RFLU_InterpolateWrapper, &
                           RFLU_PrintFlowInfoWrapper, &
                           RFLU_PrintGridInfo, &
                           RFLU_PrintHeader, &
                           RFLU_SetPatchPlotFlags, & 
                           RFLU_SetVarInfoWrapper, &
                           RFLU_SetVars, & 
                           RFLU_SetVarsWrapper
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_level), POINTER :: levels(:)

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: casename,RCSIdentString
  INTEGER :: errorFlag,iPatchSerial,iReg
  TYPE(t_region), POINTER :: pRegion,pRegionSerial
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatchSerial  

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_MergePostProcessRegions.F90,v $ $Revision: 1.29 $'

  global => levels(1)%regions(1)%global

  CALL RegisterFunction(global,'RFLU_MergePostProcessRegions',&
  'RFLU_MergePostProcessRegions.F90')

#ifndef NO_TECPLOT
! ******************************************************************************
! Initialize TECPLOT interface
! ******************************************************************************

  IF ( global%postOutputFormat == POST_OUTPUT_FORMAT_TECPLOT ) THEN 
    CALL RFLU_TEC_Init(global)
  END IF ! global%postOutputFormat 
#endif

! ******************************************************************************
! Allocate memory for serial region
! ******************************************************************************

  pRegionSerial => levels(1)%regions(0)
  
  CALL RFLU_ReadDimensionsWrapper(pRegionSerial)
  CALL RFLU_CreateGrid(pRegionSerial) 

  IF ( pRegionSerial%grid%nPatches > 0 ) THEN      
    CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)    
  END IF ! pRegionSerial%grid%nPatches  

  CALL RFLU_CreateCellMapping(pRegionSerial)
  CALL RFLU_ReadLoc2GlobCellMapping(pRegionSerial)
  CALL RFLU_BuildGlob2LocCellMapping(pRegionSerial) 
    
  CALL RFLU_AllocMemSolWrapper(pRegionSerial) 
  CALL RFLU_AllocateMemoryGSpeeds(pRegionSerial)     
  CALL RFLU_SetVarInfoWrapper(pRegionSerial) 
  
  IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
    CALL RFLU_CreatePatchCoeffs(pRegionSerial)   
  END IF ! global%patchCoeffFlag

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN 
    CALL PLAG_SetDimensions(pRegionSerial,0)
  END IF ! global%plagUsed
#endif

! ******************************************************************************  
! Loop over regions and merge region data into serial region
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg)

    CALL RFLU_ReadDimensionsWrapper(pRegion)
    CALL RFLU_CreateGrid(pRegion) 

    IF ( pRegion%grid%nPatches > 0 ) THEN      
      CALL RFLU_ReadBCInputFileWrapper(pRegion)    
    END IF ! pRegion%grid%nPatches

    CALL RFLU_ReadGridWrapper(pRegion)  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      CALL RFLU_PrintGridInfo(pRegion)
    END IF ! global%verbLevel

    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion) 

    CALL RFLU_RNMB_CreatePC2SCMap(pRegion)
    CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
    CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)    
    CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)

    CALL RFLU_AllocMemSolWrapper(pRegion) 
    CALL RFLU_SetVarInfoWrapper(pRegion)
  
    IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN   
      CALL RFLU_CreatePatchCoeffs(pRegion)
    END IF ! global%patchCoeffFlag 

    IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
      CALL RFLU_ReadFlowWrapper(pRegion)
      
      IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)
      END IF ! global%patchCoeffFlag
      
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        CALL RFLU_PrintFlowInfoWrapper(pRegion)
      END IF ! global%verbLevel       
    END IF ! global%postPlotType
            
    CALL RFLU_MERG_MergeGrid(pRegion,pRegionSerial)
    
    IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
      CALL RFLU_MERG_MergeSolWrapper(pRegion,pRegionSerial)
      
      IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
        CALL RFLU_MERG_MergePatchCoeffs(pRegion,pRegionSerial)
      END IF ! global%patchCoeffFlag    
    END IF ! global%postPlotType

    CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)
    CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
    CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)

    IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
      CALL RFLU_DestroyPatchCoeffs(pRegion)
    END IF ! global%patchCoeffFlag  
      
    CALL RFLU_DeallocMemSolWrapper(pRegion)
    CALL RFLU_DestroyCellMapping(pRegion) 
    CALL RFLU_DestroyGrid(pRegion)
  END DO ! iReg    

! ******************************************************************************  
! Serial region is now available
! ******************************************************************************

  global%nRegions      = 1
  global%nRegionsLocal = 1

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    CALL RFLU_PrintGridInfo(pRegionSerial)
  END IF ! global%verbLevel
  
! ==============================================================================
! Define dependent variables of merged region
! ==============================================================================

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
    CALL RFLU_SetVarsWrapper(pRegionSerial,1,pRegionSerial%grid%nCellsTot)  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      CALL RFLU_PrintFlowInfoWrapper(pRegionSerial)
    END IF ! global%verbLevel       

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        CALL PLAG_PrintDimensions(pRegionSerial)
      END IF ! global%verbLevel
    END IF ! global%plagUsed
#endif
  END IF ! global%postPlotType

! ==============================================================================
! Write files for merged region
! ==============================================================================

  IF ( global%postWriteMergeFlag .EQV. .TRUE. ) THEN
    CALL RFLU_WriteDimensionsWrapper(pRegionSerial,WRITE_DIMENS_MODE_FORCE)
    CALL RFLU_WriteLoc2GlobCellMapping(pRegionSerial)
    CALL RFLU_WriteGridWrapper(pRegionSerial)

    IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
      CALL RFLU_WriteFlowWrapper(pRegionSerial)
    END IF ! global%postPlotType
  END IF ! global%postWriteMergeFlag

! ==============================================================================
! Build data structure and compute geometry of merged region
! ============================================================================== 

  CALL RFLU_CreateBVertexLists(pRegionSerial)
  CALL RFLU_BuildBVertexLists(pRegionSerial)              

  CALL RFLU_CreateFaceList(pRegionSerial)
  CALL RFLU_BuildFaceList(pRegionSerial)  
  CALL RFLU_RenumberBFaceLists(pRegionSerial)      

  IF ( RFLU_DecideBuildGeometry(pRegionSerial%global) .EQV. .TRUE. ) THEN        
    CALL RFLU_CreateGeometry(pRegionSerial)
    CALL RFLU_BuildGeometry(pRegionSerial)
  END IF ! RFLU_DecideBuildGeometry

! ******************************************************************************
! Carry out specific operations on solution. NOTE this can only be done if have
! read solution.
! ******************************************************************************

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN

! ============================================================================== 
!   Create and compute plotting variables. NOTE always use order 1 for cell
!   gradients.
! ============================================================================== 
        
    IF ( RFLU_DecideComputePlottingVars(pRegionSerial) .EQV. .TRUE. ) THEN
      CALL RFLU_CountPlottingVars(pRegionSerial)
      CALL RFLU_CreatePlottingVarMaps(pRegionSerial)
      CALL RFLU_BuildPlottingVarMaps(pRegionSerial)
      CALL RFLU_PrintPlottingVarsInfo(pRegionSerial)
      CALL RFLU_CreatePlottingVars(pRegionSerial)          

      IF ( (global%postDiscFlag .EQV. .TRUE.) .OR. & 
           (global%postVortFlag .EQV. .TRUE.) .OR. & 
           (global%postVortCoreFlag .EQV. .TRUE.) ) THEN
        CALL RFLU_CreateVert2CellList(pRegionSerial)
        CALL RFLU_BuildVert2CellList(pRegionSerial)

        CALL RFLU_SetInfoC2CStencilWrapper(pRegionSerial,1)
        CALL RFLU_CreateC2CStencilWrapper(pRegionSerial)
        CALL RFLU_BuildC2CStencilWrapper(pRegionSerial)

        CALL RFLU_DestroyVert2CellList(pRegionSerial)

        CALL RFLU_CreateWtsC2CWrapper(pRegionSerial,1)
        CALL RFLU_ComputeWtsC2CWrapper(pRegionSerial,1)
      END IF ! RFLU_DecideComputePlottingVars

      CALL RFLU_ComputePlottingVarsWrapper(pRegionSerial)       

      IF ( (global%postDiscFlag .EQV. .TRUE.) .OR. &
           (global%postVortFlag .EQV. .TRUE.) .OR. &
           (global%postVortCoreFlag .EQV. .TRUE.) ) THEN
        CALL RFLU_DestroyWtsC2CWrapper(pRegionSerial)   
        CALL RFLU_DestroyC2CStencilWrapper(pRegionSerial)        
      END IF ! global%postDiscFlag
    END IF ! pRegionSerial%mixtInput%fluidModel

! ============================================================================== 
!   Compute errors
! ============================================================================== 
        
    IF ( global%initFlowFlag == INITFLOW_FROMHARDCODE ) THEN          
      CALL RFLU_ComputeExactFlowError(pRegionSerial)  
    END IF ! global%initFlowFlag

! ==============================================================================
!   Compute errors at probe locations
! ==============================================================================

    IF ( (global%nProbes > 0) .AND. & 
         (global%postCompErrFlag .EQV. .TRUE.) .AND. & 
         (global%flowType == FLOW_UNSTEADY) .AND. &        
         (global%initFlowFlag == INITFLOW_FROMHARDCODE) ) THEN
      CALL RFLU_CreateCell2FaceList(pRegionSerial)
      CALL RFLU_BuildCell2FaceList(pRegionSerial)
      CALL RFLU_FindProbeCells(pRegionSerial) 
      CALL RFLU_DestroyCell2FaceList(pRegionSerial)
      CALL RFLU_OpenProbeFiles(pRegionSerial)
      CALL RFLU_ComputeExactFlowProbeError(pRegionSerial)
      CALL RFLU_CloseProbeFiles(pRegionSerial)
    END IF ! global%nProbes

! ==============================================================================
!   Extract data from solution
! ==============================================================================

    IF ( global%postExtractFlag .EQV. .TRUE. ) THEN 
      CALL RFLU_ExtractFlowData(pRegionSerial)
    END IF ! global%postExtractFlag
  END IF ! global%postPlotType

! ******************************************************************************
! Interpolate data from cell centers to vertices
! ******************************************************************************

  IF ( (global%postPlotType == PLOT_GRID_FLOW) .AND. & 
       (global%postInterpType /= INTERP_TYPE_NONE) ) THEN 
    CALL RFLU_CreateVert2CellList(pRegionSerial)
    CALL RFLU_BuildVert2CellList(pRegionSerial)

    IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN 
      CALL RFLU_SetInfoStencilVert2Cell(pRegionSerial,global%postInterpOrder)    
      CALL RFLU_CreateStencilVert2Cell(pRegionSerial)
      CALL RFLU_BuildStencilVert2Cell(pRegionSerial)
    END IF ! global%postInterpType
      
    CALL RFLU_AllocMemVertWrapper(pRegionSerial)
    CALL RFLU_InterpolateWrapper(pRegionSerial)
    
    IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN
      CALL RFLU_DestroyStencilVert2Cell(pRegionSerial) 
    END IF ! global%postInterpType
         
    CALL RFLU_DestroyVert2CellList(pRegionSerial)
  END IF ! global%postPlotType

! ******************************************************************************
! Write out files for postprocessing
! ******************************************************************************

  SELECT CASE ( global%postOutputFormat )   
#ifndef NO_TECPLOT
    CASE ( POST_OUTPUT_FORMAT_TECPLOT ) 
      CALL RFLU_TEC_OpenFileField(pRegionSerial)

      IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN    
        CALL RFLU_TEC_OpenFilePatch(pRegionSerial)  
      END IF ! global%patchCoeffFlag   

#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_OpenFilePnt(pRegion)
        CALL RFLU_TEC_OpenFilePatchStats(pRegion)
      END IF ! global%plagUsed
#endif

      IF ( global%postPlotVolFlag .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_BuildDataFieldVol(pRegionSerial)      
        CALL RFLU_TEC_WriteFileFieldVol(pRegionSerial)         
        CALL RFLU_TEC_DestroyDataFieldVol(pRegionSerial)
      END IF ! global%postPlotVolFlag

      CALL RFLU_SetPatchPlotFlags(pRegionSerial)

      DO iPatchSerial = 1,pRegionSerial%grid%nPatches
        pPatchSerial => pRegionSerial%patches(iPatchSerial)

        IF ( pPatchSerial%plotFlag .EQV. .TRUE. ) THEN 
          CALL RFLU_TEC_BuildDataFieldSurf(pRegionSerial,pPatchSerial)             
          CALL RFLU_TEC_WriteFileFieldSurf(pRegionSerial,pPatchSerial)
          CALL RFLU_TEC_DestroyDataFieldSurf(pRegionSerial,pPatchSerial) 
        END IF ! pPatchSerial%plotFlag
      END DO ! iPatchSerial

      IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN   
        DO iPatchSerial = 1,pRegionSerial%grid%nPatches
          pPatchSerial => pRegionSerial%patches(iPatchSerial)

          CALL RFLU_TEC_BuildDataPatch(pRegionSerial,pPatchSerial)  
          CALL RFLU_TEC_WriteFilePatch(pRegionSerial,pPatchSerial)
          CALL RFLU_TEC_DestroyDataPatch(pRegionSerial,pPatchSerial)
        END DO ! iPatchSerial
      END IF ! global%patchCoeffFlag

#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_WriteFilePnt(pRegionSerial)

        DO iPatchSerial = 1,pRegionSerial%grid%nPatches
          pPatchSerial => pRegionSerial%patches(iPatchSerial)

          IF ( pPatchSerial%plotStatsFlag .EQV. .TRUE. ) THEN
            CALL RFLU_TEC_BuildDataPatchStats(pRegionSerial,pPatchSerial)
            CALL RFLU_TEC_WriteFilePatch(pRegionSerial,pPatchSerial)
            CALL RFLU_TEC_DestroyDataPatch(pRegionSerial,pPatchSerial)
          END IF ! pPatchSerial%plotStatsFlag
        END DO ! iPatchSerial
      END IF ! global%plagUsed
#endif

      CALL RFLU_TEC_CloseFileField(global)

      IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN  
        CALL RFLU_TEC_CloseFilePatch(global)
      END IF ! global%patchCoeffFlag      

#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_CloseFilePnt(global)
        CALL RFLU_TEC_CloseFilePatchStats(global)
      END IF ! global%plagUsed
#endif
#endif
    CASE ( POST_OUTPUT_FORMAT_ENSIGHT ) 
      pRegion => levels(1)%regions(0)
      CALL RFLU_ENS_BuildDataInfo(pRegion,1) 
      CALL RFLU_ENS_WriteFileCase(global,1) 
      
      CALL RFLU_ENS_OpenFileGeometry(global,1)
      CALL RFLU_ENS_OpenFileFlowWrapper(global,1)        
              
      CALL RFLU_ENS_StorePartNumber(global)   
      CALL RFLU_ENS_WriteGridWrapper(pRegionSerial)    
      CALL RFLU_ENS_WriteFlowWrapper(pRegionSerial)        
      
      CALL RFLU_ENS_CloseFileGeometry(global)
      CALL RFLU_ENS_CloseFileFlowWrapper(global)      
      CALL RFLU_ENS_DestroyDataInfo(global)      
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%postOutputFormat

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  IF ( RFLU_DecideComputePlottingVars(pRegionSerial) .EQV. .TRUE. ) THEN
    CALL RFLU_DestroyPlottingVarMaps(pRegionSerial)
    CALL RFLU_DestroyPlottingVars(pRegionSerial)
  END IF ! RFLU_DecideComputePlottingVars
  
  IF ( RFLU_DecideBuildGeometry(pRegionSerial%global) .EQV. .TRUE. ) THEN   
    CALL RFLU_DestroyGeometry(pRegionSerial)
  END IF ! RFLU_DecideBuildGeometry  

  CALL RFLU_DestroyFaceList(pRegionSerial)                        
  CALL RFLU_DestroyBVertexLists(pRegionSerial)
  CALL RFLU_DestroyCellMapping(pRegionSerial)            
  CALL RFLU_DeallocMemSolWrapper(pRegionSerial)
  CALL RFLU_DeallocateMemoryGSpeeds(pRegionSerial)
  
  IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN  
    CALL RFLU_DestroyPatchCoeffs(pRegionSerial)
  END IF ! global%patchCoeffFlag     
    
  CALL RFLU_DestroyGrid(pRegionSerial)   

  IF ( (global%postPlotType == PLOT_GRID_FLOW) .AND. & 
       (global%postInterpType /= INTERP_TYPE_NONE) ) THEN 
    CALL RFLU_DeallocMemVertWrapper(pRegionSerial)
  END IF ! global%postPlotType
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_MergePostProcessRegions

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_MergePostProcessRegions.F90,v $
! Revision 1.29  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.28  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.27  2007/03/27 00:46:48  haselbac
! Completed merging of particle solution
!
! Revision 1.26  2007/03/19 21:42:09  haselbac
! Adapted to changes related to plotting variables
!
! Revision 1.25  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.24  2006/01/24 21:16:20  mparmar
! Now pass pRegion instead of global into RFLU_ENS_BuildDataInfo
!
! Revision 1.23  2006/01/06 22:21:40  haselbac
! Adapted to name changes
!
! Revision 1.22  2005/12/10 16:56:03  haselbac
! Added RFLU_DecideComputePlottingVars and made accomp changes
!
! Revision 1.21  2005/11/10 02:48:41  haselbac
! Removed superfluous region pointers, now print grid info after merging
!
! Revision 1.20  2005/11/04 15:13:48  haselbac
! Added priniting of flow info after merging
!
! Revision 1.19  2005/10/27 19:22:30  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.18  2005/10/05 20:13:14  haselbac
! Added support for ENSIGHT
!
! Revision 1.17  2005/10/05 16:20:11  haselbac
! Bug fix: Added missing USE RFLU_ModStencilCells
!
! Revision 1.16  2005/10/05 14:27:44  haselbac
! Adapted to changes in stencil modules, added use of vertex list module
!
! Revision 1.15  2005/09/21 19:41:51  haselbac
! Bug fix: Incorrect region pointer for RFLU_SetPatchPlotFlags
!
! Revision 1.14  2005/08/18 18:48:16  haselbac
! Added postVortCoreFlag to IF statements
!
! Revision 1.13  2005/08/10 00:36:13  haselbac
! Adapted to changes in RFLU_ModPlottingVars
!
! Revision 1.12  2005/08/09 01:10:49  haselbac
! Added IFs for patch coeffs, writing serial files, plotting patches
!
! Revision 1.11  2005/07/25 12:23:30  haselbac
! Added vorticity plotting variables
!
! Revision 1.10  2005/06/11 20:35:03  haselbac
! Bug fix: Added IFs for postPlotType, now read dims for serial region at t = 0 always
!
! Revision 1.9  2005/05/03 20:39:34  haselbac
! Only compute plotting vars if postDiscFlag is TRUE
!
! Revision 1.8  2005/05/01 14:22:50  haselbac
! Added postprocessing of plotting vars
!
! Revision 1.7  2005/04/29 12:53:15  haselbac
! Added calls to compute errors at probe locations
!
! Revision 1.6  2005/04/15 15:08:32  haselbac
! Added merging of patch coeffs, modified call to RFLU_SetVarsWrapper
!
! Revision 1.5  2005/01/30 22:04:55  haselbac
! Bug fix for no interpolation of data
!
! Revision 1.4  2005/01/20 14:53:10  haselbac
! Adapted to changed RNMB routines
!
! Revision 1.3  2005/01/17 19:51:44  haselbac
! Bug fix: Now pass proper pointer to RFLU_AllocateMemoryGSpeeds
!
! Revision 1.2  2005/01/03 16:05:49  haselbac
! Adapted to changes in RFLU_ModStencils
!
! Revision 1.1  2004/12/29 20:58:54  haselbac
! Initial revision
!
! ******************************************************************************







