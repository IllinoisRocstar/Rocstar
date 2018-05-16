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
! Purpose: Wrapper for postprocessing results without merging the regions.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PostProcessRegionsCommon1.F90,v 1.12 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PostProcessRegionsCommon1(pRegion,postInfoFileExists)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModParameters

  USE RFLU_ModAllocateMemory 
  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModDeallocateMemory
  USE RFLU_ModDimensions
  USE RFLU_ModExtractFlowData  
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModInterpolation
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModPlottingVars
  USE RFLU_ModProbes
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsVert
  USE RFLU_ModVertexLists
  USE RFLU_ModWeights

#ifdef PLAG
  USE PLAG_ModSurfStats
#endif

  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, &
                           RFLU_AllocMemVertWrapper, & 
                           RFLU_ClosePostInfo, &
                           RFLU_ComputeExactFlowError, &
                           RFLU_ComputeExactFlowProbeError, &
                           RFLU_CreateGrid, & 
                           RFLU_DeallocMemSolWrapper, &
                           RFLU_DeallocMemVertWrapper, &
                           RFLU_DecideBuildGeometry, &
                           RFLU_DecideBuildStencilsWeights, & 
                           RFLU_DestroyGrid, &
                           RFLU_InterpolateWrapper, &
                           RFLU_OpenPostInfo, &
                           RFLU_PrintFlowInfoWrapper, &
                           RFLU_PrintGridInfo, &
                           RFLU_PrintHeader, &
                           RFLU_ReadPostInfo, & 
                           RFLU_SetPatchPlotFlags, & 
                           RFLU_SetVarInfoWrapper, &
                           RFLU_SetVarsWrapper
    
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(IN) :: postInfoFileExists
  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iPatch,iReg
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PostProcessRegionsCommon1.F90,v $ $Revision: 1.12 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PostProcessRegionsCommon1',&
  'RFLU_PostProcessRegionsCommon1.F90')

! ******************************************************************************
! Read dimensions file, create grid, and read bc file
! ******************************************************************************

  CALL RFLU_ReadDimensionsWrapper(pRegion)
  CALL RFLU_CreateGrid(pRegion) 

  IF ( pRegion%grid%nPatches > 0 ) THEN      
    CALL RFLU_ReadBCInputFileWrapper(pRegion)    
  END IF ! pRegion%grid%nPatches

! ******************************************************************************
! Read grid, build data structures
! ******************************************************************************

  CALL RFLU_ReadGridWrapper(pRegion)  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    CALL RFLU_PrintGridInfo(pRegion)
  END IF ! global%verbLevel

  CALL RFLU_CreateCellMapping(pRegion)
  CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
  CALL RFLU_BuildGlob2LocCellMapping(pRegion)

  CALL RFLU_CreateBVertexLists(pRegion)
  CALL RFLU_BuildBVertexLists(pRegion)              

  CALL RFLU_CreateFaceList(pRegion)
  CALL RFLU_BuildFaceList(pRegion) 
  CALL RFLU_RenumberBFaceLists(pRegion)      

  IF ( RFLU_DecideBuildGeometry(global) .EQV. .TRUE. ) THEN
    CALL RFLU_CreateGeometry(pRegion)       
    CALL RFLU_BuildGeometry(pRegion)      
  END IF ! RFLU_DecideBuildGeometry

! ******************************************************************************
! Read flow solution, define properties of mixture
! ******************************************************************************

  CALL RFLU_AllocMemSolWrapper(pRegion)
  CALL RFLU_SetVarInfoWrapper(pRegion)
  CALL RFLU_AllocateMemoryGSpeeds(pRegion) 

  IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN           
    CALL RFLU_CreatePatchCoeffs(pRegion)
  END IF ! global%patchCoeffFlag

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN           
    CALL PLAG_CreateSurfStats(pRegion)
  END IF ! global%plagUsed
#endif

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN 
    CALL RFLU_ReadFlowWrapper(pRegion)
    CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      CALL RFLU_PrintFlowInfoWrapper(pRegion)
    END IF ! global%verbLevel

    IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN      
      CALL RFLU_ReadPatchCoeffsWrapper(pRegion)
    END IF ! global%patchCoeffFlag

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN           
      CALL PLAG_ReadSurfStatsWrapper(pRegion)
    END IF ! global%plagUsed
#endif     
  END IF ! global%postPlotType

! ******************************************************************************
!  Compute plotting variables and errors. NOTE need to be done together because
!  need gradients and weights for both and would be wasteful to do twice. 
!  NOTE always use order 1 for cell gradients.
! ******************************************************************************

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
  
! ==============================================================================
!   Build stencils and compute weights  
! ==============================================================================
  
    IF ( RFLU_DecideBuildStencilsWeights(global) .EQV. .TRUE. ) THEN
      CALL RFLU_CreateVert2CellList(pRegion)
      CALL RFLU_BuildVert2CellList(pRegion)

      IF ( pRegion%mixtInput%stencilDimensCells == 1 ) THEN 
        CALL RFLU_CreateCell2FaceList(pRegion)
        CALL RFLU_BuildCell2FaceList(pRegion)
      END IF ! pRegion%mixtInput%stencilDimensCells

      CALL RFLU_SetInfoC2CStencilWrapper(pRegion,1)
      CALL RFLU_CreateC2CStencilWrapper(pRegion)
      CALL RFLU_BuildC2CStencilWrapper(pRegion,constrInput=CONSTR_NONE)

      IF ( pRegion%mixtInput%stencilDimensCells == 1 ) THEN 
        CALL RFLU_DestroyCell2FaceList(pRegion)
      END IF ! pRegion%mixtInput%stencilDimensCells

      CALL RFLU_DestroyVert2CellList(pRegion) 

      CALL RFLU_CreateWtsC2CWrapper(pRegion,1)
      CALL RFLU_ComputeWtsC2CWrapper(pRegion,1)
    END IF ! RFLU_DecideBuildStencilsWeights  

! ==============================================================================
!   Compute plotting variables
! ==============================================================================
    
    IF ( RFLU_DecideComputePlottingVars(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_CountPlottingVars(pRegion)
      CALL RFLU_CreatePlottingVarMaps(pRegion)
      CALL RFLU_BuildPlottingVarMaps(pRegion)
      CALL RFLU_PrintPlottingVarsInfo(pRegion)
      CALL RFLU_CreatePlottingVars(pRegion)         
      CALL RFLU_ComputePlottingVarsWrapper(pRegion)       
    END IF ! RFLU_DecideComputePlottingVars

! ==============================================================================
!   Compute errors
! ==============================================================================
    
    IF ( global%initFlowFlag == INITFLOW_FROMHARDCODE ) THEN   
      CALL RFLU_ComputeExactFlowError(pRegion)  
    END IF ! global%initFlowFlag

! ==============================================================================
!   Destroy stencils and weights  
! ==============================================================================
  
    IF ( RFLU_DecideBuildStencilsWeights(global) .EQV. .TRUE. ) THEN
      CALL RFLU_DestroyWtsC2CWrapper(pRegion)   
      CALL RFLU_DestroyC2CStencilWrapper(pRegion)         
    END IF ! RFLU_DecideBuildStencilsWeights      
  END IF ! global%postPlotType

! ******************************************************************************
! Compute errors at probe locations
! ******************************************************************************

  IF ( (global%nProbes > 0) .AND. & 
       (global%postCompErrFlag .EQV. .TRUE.) .AND. & 
       (global%flowType == FLOW_UNSTEADY) .AND. &          
       (global%initFlowFlag == INITFLOW_FROMHARDCODE) ) THEN
    CALL RFLU_CreateCell2FaceList(pRegion)
    CALL RFLU_BuildCell2FaceList(pRegion)
    CALL RFLU_FindProbeCells(pRegion) 
    CALL RFLU_DestroyCell2FaceList(pRegion)
    CALL RFLU_OpenProbeFiles(pRegion)
    CALL RFLU_ComputeExactFlowProbeError(pRegion)
    CALL RFLU_CloseProbeFiles(pRegion)
  END IF ! global%nProbes

! ******************************************************************************
! Extract data from solution
! ******************************************************************************

  IF ( global%postExtractFlag .EQV. .TRUE. ) THEN 
    CALL RFLU_ExtractFlowData(pRegion)
  END IF ! global%postExtractFlag

! ******************************************************************************
! Interpolate data from cell centers to vertices
! ******************************************************************************

  IF ( (global%postPlotType == PLOT_GRID_FLOW) .AND. & 
       (global%postInterpType /= INTERP_TYPE_NONE) ) THEN 
    CALL RFLU_CreateVert2CellList(pRegion)
    CALL RFLU_BuildVert2CellList(pRegion) 

    IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN   
      CALL RFLU_SetInfoStencilVert2Cell(pRegion,global%postInterpOrder)                                                                                                      
      CALL RFLU_CreateStencilVert2Cell(pRegion)
      CALL RFLU_BuildStencilVert2Cell(pRegion) 
    END IF ! global%postInterpType

    CALL RFLU_AllocMemVertWrapper(pRegion)           
    CALL RFLU_InterpolateWrapper(pRegion)

    IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN
      CALL RFLU_DestroyStencilVert2Cell(pRegion)
    END IF ! global%postInterpType  

    CALL RFLU_DestroyVert2CellList(pRegion)   
  END IF ! global%postPlotType

! ******************************************************************************
! Read information on special cells and faces
! ******************************************************************************

  IF ( postInfoFileExists .EQV. .TRUE. ) THEN 
    CALL RFLU_ReadPostInfo(pRegion,INFOFILE_READMODE_DATA)    
  END IF ! postInfoFileExists       
        
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PostProcessRegionsCommon1

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PostProcessRegionsCommon1.F90,v $
! Revision 1.12  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2007/03/19 21:45:19  haselbac
! Adapted to changes related to plotting variables
!
! Revision 1.9  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.8  2006/04/07 14:57:43  haselbac
! Adapted to new stencilDimens param
!
! Revision 1.7  2006/01/06 22:20:51  haselbac
! Rewrote so can compute errors of gradients
!
! Revision 1.6  2005/12/10 16:57:30  haselbac
! Added use of RFLU_DecideComputePlottingVars
!
! Revision 1.5  2005/12/10 13:57:17  haselbac
! Bug fix: Only need to build geometry when necessary
!
! Revision 1.4  2005/11/10 16:51:29  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.3  2005/10/27 19:22:30  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.2  2005/10/05 20:53:00  haselbac
! Fixed missing module bugs
!
! Revision 1.1  2005/10/05 20:23:34  haselbac
! Initial revision
!
! ******************************************************************************







