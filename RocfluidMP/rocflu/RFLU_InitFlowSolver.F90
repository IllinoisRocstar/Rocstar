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
! Purpose: Main initialization routine of Rocflu-MP.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowSolver.F90,v 1.123 2008/12/06 08:44:30 mtcampbe Exp $
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

#ifdef GENX
SUBROUTINE RFLU_InitFlowSolver(globalGenx,initialTime,communicator,genxHandle, &
                               surfWinNameInput,volWinNameInput,handleObtain)
#else
SUBROUTINE RFLU_InitFlowSolver(casename,verbLevel,global,levels)
#endif

  USE ModDataTypes
#ifdef GENX
  USE ModRocstar, ONLY: t_globalGenx
#endif  
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModMixture, ONLY: t_mixt_input  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModMPI
  USE ModParameters  
     
  USE RFLU_ModBFaceGradAccessList
  USE RFLU_ModBoundLists   
  USE RFLU_ModBoundXvUtils
  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
  USE RFLU_ModCellMapping
  USE RFLU_ModCommLists
  USE RFLU_ModDimensionality  
  USE RFLU_ModDimensions
  USE RFLU_ModEdgeList   
  USE RFLU_ModFaceList
  USE RFLU_ModForcesMoments
  USE RFLU_ModGeometry
  USE RFLU_ModGridSpeedUtils
  USE RFLU_ModInCellTest
  USE RFLU_ModInterpolation
  USE RFLU_ModMPI
  USE RFLU_ModOLES
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModPatchUtils
  USE RFLU_ModProbes
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteBcDataFile
  USE RFLU_ModReadWriteFlow   
  USE RFLU_ModReadWriteGrid 
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModRegionMapping
  USE RFLU_ModRenumberings
  USE RFLU_ModStencilsBFaces
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsFaces
  USE RFLU_ModStencilsUtils
  USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_HaveSyPePatches, & 
                                      RFLU_SYPE_ReadTransforms
  USE RFLU_ModVertexLists
  USE RFLU_ModWeights
    
#ifdef GENX
  USE RFLU_ModRocstarAdmin
  USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_GetDataInterf, & 
                            RFLU_GENX_GetDimensionsDerived,&
                            RFLU_GENX_GetGlobalData
#endif  
    
#ifdef PETSC
  USE RFLU_ModPETScAdmin
!  USE RFLU_ModPETScPoisson
  USE RFLU_ModPETScNewtonKrylov
#endif
  
#ifdef PLAG
  USE PLAG_ModSurfStats, ONLY: PLAG_CreateSurfStats, & 
                               PLAG_ReadSurfStatsWrapper
#endif    
    
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformWrapper
    
  USE ModInterfaces, ONLY: RFLU_AllocateMemoryTbc, &
                           RFLU_AllocateMemoryWrapper, &
                           RFLU_BuildDataStruct, & 
                           RFLU_ComputeIntegrals1245OLES, &
                           RFLU_ComputeIntegralValues, &
                           RFLU_CreateGrid, &
                           RFLU_DecideNeedBGradFace, &
                           RFLU_GetUserInput, & 
                           RFLU_InitGlobal, &
                           RFLU_MoveGridWrapper, &
                           RFLU_OpenConverFile, &
                           RFLU_OpenStatsFileOLES, &
                           RFLU_OpenTotalMassFile, &
                           RFLU_PrintFlowInfo, &
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintHeader, &
                           RFLU_RandomInit, &
                           RFLU_ReadIntegrals1245OLES, &
                           RFLU_ReadRestartInfo, &
                           RFLU_ReadTbcInputFile, &
                           RFLU_SetModuleType, & 
                           RFLU_SetMoveGridOptions, & 
                           RFLU_SetRestartTimeFlag, &
                           RFLU_SetVarsContWrapper, &
                           RFLU_SetVarsDiscWrapper, & 
                           RFLU_SetVarInfoWrapper, & 
                           RFLU_SetVarsWrapper, & 
                           RFLU_WriteIntegrals1245OLES, &
                           RFLU_WriteVersionString, &
                           WriteTotalMass
                                                
#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_CheckCouplingInput, &
                           RFLU_GetBoundaryValues, &
                           RFLU_PutBoundaryValues
  USE RFLU_ModRocstarTools, ONLY: RFLU_GENX_InitBFLAG
#endif
#ifdef STATS
  USE ModStatsRoutines, ONLY: StatMapping, InitStatistics
#endif
  
  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ******************************************************************************
! Arguments
! ******************************************************************************

#ifdef GENX
  CHARACTER(*), INTENT(IN) :: surfWinNameInput,volWinNameInput
  INTEGER, INTENT(IN) :: communicator,genxHandle,handleObtain
  DOUBLE PRECISION, INTENT(IN) :: initialTime
  TYPE(t_globalGenx), POINTER :: globalGenx
#else
  CHARACTER(CHRLEN), INTENT(IN) :: casename
  INTEGER, INTENT(IN) :: verbLevel
  TYPE(t_global), POINTER :: global
#endif
  TYPE(t_level), POINTER :: levels(:)  
  
! ******************************************************************************
! Locals
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString,msg
  LOGICAL :: fileExists,moveGrid
  INTEGER :: errorFlag,iPatch,iReg
  REAL(RFREAL) :: currentTime
  TYPE(t_grid) :: grid
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_level), POINTER :: pLevel
  TYPE(t_mixt_input), POINTER :: pMixtInput  
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

#ifdef GENX
  INTEGER :: verbLevel
  CHARACTER(CHRLEN) :: casename
  TYPE(t_global), POINTER :: global  
#endif 

! ******************************************************************************
! Start, initialize some variables
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowSolver.F90,v $ $Revision: 1.123 $'
 
  moveGrid = .FALSE.
 
! ******************************************************************************
! Set global pointer and initialize global type, register function
! ****************************************************************************** 
 
#ifdef GENX
  global => globalGenx%global 
     
  CALL RFLU_InitGlobal(" ",VERBOSE_LOW,communicator,global)
  CALL RFLU_GENX_ReadCtrlFile(global)
  CALL RFLU_GENX_StoreCommunicator(global,communicator)

! ******************************************************************************
! Initialize MPI - ROCSTAR
! ****************************************************************************** 

  CALL MPI_Comm_size(global%mpiComm,global%nProcAlloc,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error
    
  CALL MPI_Comm_rank(global%mpiComm,global%myProcid,errorFlag)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error

  CALL RFLU_GENX_StoreNamesHandles(global,surfWinNameInput,volWinNameInput, & 
                                   handleObtain)

  CALL COM_set_verbose(global%verbLevelCOM)

  global%timeStamp   = initialTime
  global%currentTime = initialTime
#else
  CALL RFLU_InitGlobal(casename,verbLevel,MPI_COMM_WORLD,global)

! ******************************************************************************
! Initialize MPI Native
! ****************************************************************************** 

  CALL MPI_Init(errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error
  
  CALL MPI_Comm_size(global%mpiComm,global%nProcAlloc,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error
    
  CALL MPI_Comm_rank(global%mpiComm,global%myProcid,errorFlag)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error

#endif

  CALL RFLU_SetModuleType(global,MODULE_TYPE_SOLVER)

  CALL RegisterFunction(global,'RFLU_InitFlowSolver',&
  'RFLU_InitFlowSolver.F90')

#ifndef GENX
#ifdef ROCPROF
! ******************************************************************************
! Initialize Rocprof. NOTE done here so everything can be profiled and also bcos
! RFLU_SetVars is called here.
! ****************************************************************************** 

  CALL Rocprof_Init("Rocflu",global%myProcid)
#endif
#endif

! ******************************************************************************
! Print header and check for stop file 
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN 
    CALL RFLU_WriteVersionString(global)
  
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      CALL RFLU_PrintHeader(global)
    END IF ! global%verbLevel
  END IF ! global

  IF ( global%myProcid == MASTERPROC ) THEN 
    INQUIRE(FILE="STOP",EXIST=fileExists)
    IF ( fileExists .EQV. .TRUE. ) THEN
      CALL ErrorStop(global,ERR_STOPFILE_FOUND,__LINE__)
    END IF ! fileExists
  END IF ! global%myProcid  

! ******************************************************************************
! Read processor mapping file, prepare data structure, set level pointer
! ******************************************************************************
 
  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_ALL,global%myProcId)
  
  IF ( global%nRegions == 1 ) THEN 
    CALL RFLU_ImposeRegionMappingSerial(global)
  END IF ! global

  IF ( global%nProcs /= global%nProcAlloc ) THEN 
    CALL ErrorStop(global,ERR_PROC_MISMATCH,__LINE__)
  END IF ! global%nProcs
  
  CALL RFLU_BuildDataStruct(global,levels)
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)
    
#ifdef GENX    
  globalGenx%levels => levels
  pLevel => globalGenx%levels(1) ! single-level grids for now      
#else
  pLevel => levels(1) ! single-level grids for now   
#endif

#ifdef GENX
! ******************************************************************************
! Create windows and new attributes. NOTE must only be done once for every 
! process.
! ******************************************************************************
  
  pRegion => pLevel%regions(1)
  CALL RFLU_GENX_CreateWindows(pRegion,communicator)
  CALL RFLU_GENX_CreateAttrWrapper(pRegion)
  
! ******************************************************************************
! Load modules
! ******************************************************************************

  CALL Rocmop_load_module(TRIM(global%winName)//'-MOP')  
#endif           
           
! ******************************************************************************
! Read user input
! ******************************************************************************

  CALL RFLU_GetUserInput(pLevel%regions)

! ******************************************************************************
! Initialize random number generator
! ******************************************************************************

  CALL RFLU_RandomInit(pLevel%regions)

! ******************************************************************************
! Read restart info. NOTE must not read restart info file for GENx runs because 
! will overwrite global%currentTime.
! ******************************************************************************
  
#ifndef GENX
  CALL RFLU_ReadRestartInfo(global)
#endif
  CALL RFLU_SetRestartTimeFlag(global)
  
! ******************************************************************************
! Read dimensions file
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_ReadDimensionsWrapper(pRegion)     
  END DO ! iReg

! ******************************************************************************
! Determine whether have moving grid
! ******************************************************************************

  moveGridLoop: DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( pRegion%mixtInput%movegrid .EQV. .TRUE. ) THEN
      moveGrid = .TRUE.       
      EXIT moveGridLoop
    END IF ! pRegion
  END DO moveGridLoop

! ******************************************************************************
! Allocate memory for grid and borders
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    CALL RFLU_CreateGrid(pRegion)
    CALL RFLU_COMM_CreateBorders(pRegion,CREATE_BORDERS_MODE_DIM_KNOWN)
  END DO ! iReg

! ******************************************************************************
! Read boundary condition file, must be done after creation of grid, because 
! grid and patches are created there.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( pRegion%grid%nPatches > 0 ) THEN      
      CALL RFLU_ReadBcInputFileWrapper(pRegion)
      CALL RFLU_AllocateMemoryTbc(pRegion)
      CALL RFLU_ReadTbcInputFile(pRegion)
    END IF ! pRegion%grid%nPatches    
  END DO ! iReg

#ifdef GENX
! ******************************************************************************
! Create surface and and register volume and surface grid. NOTE must be done 
! after having read comm lists because of pconn.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    CALL RFLU_GENX_GetDimensionsDerived(pRegion)
    CALL RFLU_GENX_CreateGridSurf(pRegion)    
    CALL RFLU_GENX_CreatePConn(pRegion)
    CALL RFLU_GENX_RegisterGridVol(pRegion)
    CALL RFLU_GENX_RegisterGridSurf(pRegion)
  END DO ! iReg  
#endif
  
! ******************************************************************************
! Read grid file
! ******************************************************************************
  
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_ReadGridWrapper(pRegion)
  END DO ! iReg
 
! ******************************************************************************
! Print grid information
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_PrintGridInfo(pRegion)
    END DO ! iReg
  END IF ! global%verbLevel

!*******************************************************************************
! Register Time Zooming Parameter
!*******************************************************************************
   
#ifdef GENX
  CALL RFLU_GENX_GetGlobalData(global)
#endif

! ******************************************************************************
! Build data stucture, part 1
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)      
  END DO ! iReg 

! ******************************************************************************
! Create and read comm lists, set proc ids for borders. Must be done before 
! face list is constructed because reorientation of actual-virtual faces
! requires knowledge of communication lists.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    pGrid   => pRegion%grid   
    
    IF ( pGrid%nBorders > 0 ) THEN 
      CALL RFLU_COMM_CreateCommLists(pRegion)
      CALL RFLU_COMM_ReadCommLists(pRegion)
      CALL RFLU_COMM_GetProcLocRegIds(pRegion)
    END IF ! pRegion%grid%nBorders
  END DO ! iReg  

! ******************************************************************************
! Build data stucture, part 2. Building the face list must be done after the 
! communication lists have been updated to take into account renumbering of 
! cells. Build cell-to-face list anyway, although only needed when have probes 
! or when running with higher-order scheme.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    CALL RFLU_CreateBVertexLists(pRegion)
    CALL RFLU_BuildBVertexLists(pRegion)        
  END DO ! iReg

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    CALL RFLU_CreateFaceList(pRegion)
    CALL RFLU_BuildFaceList(pRegion)
    CALL RFLU_RenumberBFaceLists(pRegion)

! TEMPORARY
!    IF ( pRegion%grid%nCellsTot > pRegion%grid%nCells ) THEN 
!      CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)
!      CALL RFLU_RNMB_CreatePC2SCMap(pRegion)
!      CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
!      CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)
!      CALL RFLU_ReorientFaces(pRegion) 
!! TEMPORARY
!!      CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)
!!      CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)
!!      CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
!! END TEMPORARY 
!    END IF ! pGrid%nCellsTot
! END TEMPORARY
  END DO ! iReg 

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_CreateCell2FaceList(pRegion)
    CALL RFLU_BuildCell2FaceList(pRegion)
  END DO ! iReg

#ifdef PLAG
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    IF ( (pRegion%global%plagUsed .EQV. .TRUE.) .AND. & 
         (pRegion%grid%nFacesAV > 0) ) THEN
      CALL RFLU_CreateAVFace2BorderList(pRegion)
      CALL RFLU_BuildAVFace2BorderList(pRegion)
      
      CALL RFLU_CreateAVFace2PatchList(pRegion)
      CALL RFLU_BuildAVFace2PatchList(pRegion)      
    END IF ! pRegion%global%plagUsed
  END DO ! iReg
#endif

! ******************************************************************************
! Build data stucture, part 3
! ******************************************************************************

  IF ( moveGrid .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_CreateEdgeList(pRegion)       
      CALL RFLU_BuildEdgeList(pRegion)      
      CALL RFLU_CreateEdge2CellList(pRegion)              
      CALL RFLU_BuildEdge2CellList(pRegion)
      CALL RFLU_DestroyEdge2CellList(pRegion)      
    END DO ! iReg 
  END IF ! moveGrid

! ******************************************************************************
! Check topology
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => pLevel%regions(iReg)              
    CALL RFLU_123D_CheckTopology(pRegion)
  END DO ! iReg  

! ******************************************************************************
! Build boundary-face gradient access list for viscous flows
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => pLevel%regions(iReg) 
    
    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN          
      CALL RFLU_CreateBFaceGradAccessList(pRegion)
      CALL RFLU_BuildBFaceGradAccessList(pRegion)
    END IF ! pRegion%mixtInput
  END DO ! iReg  

! ******************************************************************************
! Allocate memory. NOTE must be done after having built face lists and after 
! reading boundary-condition file because of particle module.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_AllocateMemoryWrapper(pRegion)
    
    CALL RFLU_CreatePatchCoeffs(pRegion)
   
    IF ( global%forceFlag .EQV. .TRUE. ) THEN 
      CALL RFLU_CreateForcesMoments(pRegion)
      CALL RFLU_CreateGlobalThrustFlags(pRegion)
    END IF ! global%forceFlag
  END DO ! iReg
    
! ------------------------------------------------------------------------------
! Set Global thrust flags from local thrust flags 
! ------------------------------------------------------------------------------

  CALL RFLU_SetGlobalThrustFlags(pLevel%regions)

! ------------------------------------------------------------------------------
! Allocate memory for boundary variables if required, i.e., if NSCBC==1 
! ------------------------------------------------------------------------------

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_CreateVarsCv(pRegion)
      CALL RFLU_BXV_CreateVarsDv(pRegion)
      CALL RFLU_BXV_CreateVarsTStep(pRegion)
    END IF ! RFLU_NSCBC_DecideHaveNSCBC(pRegion)    

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN 
      CALL PLAG_CreateSurfStats(pRegion)
    END IF ! global%plagUsed
#endif    
  END DO ! iReg

! ******************************************************************************
! Set variable info. NOTE must be done after allocation of memory for solution.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_SetVarInfoWrapper(pRegion)
  END DO ! iReg  

#ifdef STATS
! ******************************************************************************
! Statistics mapping (must be done before RFLU_GENX_RegisterData)
! ******************************************************************************

  IF (global%flowType == FLOW_UNSTEADY .AND. global%doStat==ACTIVE) THEN
    CALL StatMapping( global )
  ENDIF
#endif

#ifdef GENX
! ******************************************************************************
! Register data with Roccom for interaction in GENX. NOTE must be done after 
! memory allocation.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    CALL RFLU_GENX_RegisterDataWrapper(pRegion)
    CALL RFLU_GENX_CreateWindowsDone(pRegion)
  END DO ! iReg  
#endif  

! ******************************************************************************
! Compute geometry
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)   
    CALL RFLU_CreateGeometry(pRegion)    
    CALL RFLU_BuildGeometry(pRegion)
    CALL RFLU_123D_CheckGeometryWrapper(pRegion)    
  END DO ! iReg 

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)   
    CALL RFLU_ComputePatchNormalsLocal(pRegion)       
  END DO ! iReg 

  CALL RFLU_ComputePatchNormalsGlobal(pLevel%regions)
  
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)   
  CALL RFLU_CheckPatchBcConsistency(pRegion)
  END DO ! iReg 

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    pMixtInput => pRegion%mixtInput

    IF ( pMixtInput%fluidModel == FLUID_MODEL_INCOMP ) THEN
      CALL RFLU_CreateFaceDist(pRegion)
      CALL RFLU_ComputeFaceDist(pRegion)
    END IF ! pMixtInput%fluidModel
  END DO ! iReg

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    
    IF ( RFLU_SYPE_HaveSyPePatches(pRegion) .EQV. .TRUE. ) THEN   
      CALL RFLU_SYPE_ReadTransforms(pRegion)       
    END IF ! RFLU_SYPE_HaveSyPePatches
  END DO ! iReg

! ******************************************************************************
! Build stencils. NOTE This must be done after the geometry is computed 
! because the stencils are now tested for singularity, which requires the
! computation of the weights and hence the geometry.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => pLevel%regions(iReg)
    pMixtInput => pRegion%mixtInput 

    CALL RFLU_CreateVert2CellList(pRegion)
    CALL RFLU_BuildVert2CellList(pRegion)
    
    IF ( pMixtInput%spaceOrder > 1 ) THEN 
      CALL RFLU_SetInfoC2CStencilWrapper(pRegion,pMixtInput%spaceOrder-1) 
      CALL RFLU_CreateC2CStencilWrapper(pRegion)
      CALL RFLU_BuildC2CStencilWrapper(pRegion)
      CALL RFLU_BuildListCC2CStencil(pRegion) 
    END IF ! pMixtInput%spaceOrder      
    
    IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
      CALL RFLU_SetInfoF2CStencilWrapper(pRegion,pMixtInput%spaceOrder-1)    
      CALL RFLU_CreateF2CStencilWrapper(pRegion)      
      CALL RFLU_BuildF2CStencilWrapper(pRegion)   
      CALL RFLU_BuildListCF2CStencil(pRegion) 
    END IF ! pMixtInput%flowModel
      
    DO iPatch = 1,pRegion%grid%nPatches 
      pPatch => pRegion%patches(iPatch)

      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        CALL RFLU_SetInfoBF2CStencilWrapper(pRegion,pPatch, &
                                            pPatch%spaceOrder)    
        CALL RFLU_CreateBF2CStencilWrapper(pRegion,pPatch)      
        CALL RFLU_BuildBF2CStencilWrapper(pRegion,pPatch) 
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
    
    CALL RFLU_DestroyVert2CellList(pRegion)
  END DO ! iReg

! ******************************************************************************
! Compute weights. NOTE must be done after stencils are built.
! ******************************************************************************

! ==============================================================================
! Weights for cell and face gradients
! ==============================================================================

  DO iReg = 1,global%nRegionsLocal 
    pRegion => pLevel%regions(iReg)
    pMixtInput => pRegion%mixtInput   

    IF ( pMixtInput%spaceOrder > 1 ) THEN 
      CALL RFLU_CreateWtsC2CWrapper(pRegion,pMixtInput%spaceOrder-1)
      CALL RFLU_ComputeWtsC2CWrapper(pRegion,pMixtInput%spaceOrder-1)
    END IF ! pMixtInput%spaceOrder      

    IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
      CALL RFLU_CreateWtsF2CWrapper(pRegion,pMixtInput%spaceOrder-1)    
      CALL RFLU_ComputeWtsF2CWrapper(pRegion,pMixtInput%spaceOrder-1)
    END IF ! pMixtInput%flowModel
      
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        CALL RFLU_CreateWtsBF2CWrapper(pRegion,pPatch,pPatch%spaceOrder)    
        CALL RFLU_ComputeWtsBF2CWrapper(pRegion,pPatch,pPatch%spaceOrder)    
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
  END DO ! iReg   

! ==============================================================================
! Weights for optimal LES approach. NOTE that at present, this will work only
! for non-moving grids.
! ==============================================================================

  DO iReg = 1,global%nRegionsLocal   
    pRegion => pLevel%regions(iReg)
    pMixtInput => pRegion%mixtInput 

    IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN            
      CALL RFLU_FindPrototypeFacesOLES(pRegion)
      CALL RFLU_BuildStencilsOLES(pRegion)
      CALL RFLU_ComputeGeometricTermsOLES(pRegion)
      CALL RFLU_BuildSymmetryMapsOLES(pRegion)
      
!      CALL RFLU_ComputeIntegrals1245OLES(pRegion)
!      CALL RFLU_EnforceSymmetryOLES(pRegion)
!      CALL RFLU_WriteIntegrals1245OLES(pRegion)
      
      CALL RFLU_ReadIntegrals1245OLES(pRegion)
      CALL RFLU_EnforceSymmetryOLES(pRegion) 
!      CALL RFLU_WriteIntegrals1245OLES(pRegion)           
    END IF ! pMixtInput 
  END DO ! iReg 

#ifdef PLAG
! ******************************************************************************
! Compute in-cell test tolerance
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    IF ( pRegion%global%plagUsed .EQV. .TRUE. ) THEN 
      CALL RFLU_ICT_ComputeTolerance(pRegion)
    END IF ! pRegion%global%plagUsed
  END DO ! iReg
#endif

! ******************************************************************************
! Read solution data files
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)     
    CALL RFLU_ReadFlowWrapper(pRegion)
    
#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      CALL PLAG_ReadSurfStatsWrapper(pRegion)
    END IF ! plagUsed
#endif          
  END DO ! iReg

#ifdef GENX
! ******************************************************************************
! Get interface data
! ******************************************************************************
  
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)     
    CALL RFLU_GENX_GetDataInterf(pRegion)
  END DO ! iReg  
#endif

! ******************************************************************************
! Read boundary condition data (only relevant if distribution specified)
! Also Read boundary variables if required, i.e., if NSCBC==1 
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    CALL RFLU_BXV_ReadVarsWrapper(pRegion)

    IF ( RFLU_DecideReadWriteBcDataFile(pRegion) .EQV. .TRUE. ) THEN 
      CALL RFLU_ReadBcDataFile(pRegion)
    END IF ! RFLU_DecideReadWriteBcDataFile
  END DO ! iReg
  
! ******************************************************************************
! Read grid speeds. NOTE must be done after face list is constructed, because 
! need to know number of faces to be able to read grid speeds. 
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( RFLU_DecideNeedGridSpeeds(pRegion) .EQV. .TRUE. ) THEN 
      CALL RFLU_ReadGridSpeedsWrapper(pRegion)    
    END IF ! RFLU_DecideNeedGridSpeeds
  END DO ! iReg

! ******************************************************************************
! Initialize grid speed scaling factor. NOTE needs to be done for steady flows
! to work properly because grid speed is used together with scaling routines.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    CALL RFLU_InitGridSpeedScaleFactor(pRegion)
  END DO ! iReg

#ifdef STATS
! ******************************************************************************
! Statistics initialization
! ******************************************************************************

  CALL InitStatistics(pLevel%regions)
#endif

! ******************************************************************************
! Create buffers and tags for communication
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    pGrid   => pRegion%grid   
    
    IF ( pGrid%nBorders > 0 ) THEN 
      CALL RFLU_MPI_CreateBuffersWrapper(pRegion)
      CALL RFLU_MPI_SetTagsWrapper(pRegion)
      
#ifdef PLAG
      IF ( pRegion%global%plagUsed .EQV. .TRUE. ) THEN 
        CALL RFLU_MPI_CreateBufferIPclSend(pRegion)
      END IF ! pRegion%global%plagUsed
#endif      
    END IF ! pRegion%grid%nBorders
  END DO ! iReg  

! *****************************************************************************
! Initialize virtual-cell data via communication (done for GENX runs after 
! remeshing because remeshing tool does not initialize ghosts)
! *****************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_MPI_ISendWrapper(pRegion)
  END DO ! iReg 

  CALL RFLU_MPI_CopyWrapper(pLevel%regions)

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_MPI_RecvWrapper(pRegion)
  END DO ! iReg 

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_MPI_ClearRequestWrapper(pRegion)
  END DO ! iReg 

! ******************************************************************************
! Initialize dependent variables 
! ******************************************************************************
  
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)   

! TEMPORARY : move it to proper place ....
    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_SetDependentVars(pRegion)
    END IF !
! END TEMPORARY
  END DO ! iReg

#ifdef GENX
! ******************************************************************************
! Initialize Rocman
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)        
    CALL RFLU_GENX_InitRocman(pRegion,genxHandle)     
  END DO ! iReg

! ******************************************************************************
! Check that have coupled boundaries defined, must be done after field flags
! are created
! ******************************************************************************

! TEMPORARY - Breaks intermittently, without this you need to make sure that
!             your input files are correct; there is no checking...
!  CALL RFLU_CheckCouplingInput(pLevel%regions)
! END TEMPORARY

! ******************************************************************************
! Put initial data into GENX buffers. NOTE no longer get data to fill incoming
! buffers to make sure restarts are always the same regardless of the time. 
! Incoming buffers are filled at beginning of time-stepping routine
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    IF ( global%currentTime /= 0.0_RFREAL ) THEN
      CALL RFLU_GENX_InitBFLAG(pRegion)
      CALL RFLU_GetBoundaryValues(pLevel%regions(iReg))
    END IF
    
    CALL RFLU_PutBoundaryValues(pLevel%regions(iReg))            
  END DO ! iReg
#endif

! ******************************************************************************
! Set options for grid motion, must be done after boundary normals are built 
! and data fields for communication are created
! ******************************************************************************

! TEMPORARY - Disabled until merged with changes for periodic and symmetry
!             boundaries
!  IF ( moveGrid .EQV. .TRUE. ) THEN 
!    DO iReg = 1,global%nRegionsLocal
!      pRegion => pLevel%regions(iReg)
!      CALL RFLU_SetMoveGridOptions(pRegion) 
!    END DO ! iReg
!  END IF ! moveGrid
! END TEMPORARY

! ******************************************************************************
! Find probe cells and print information on probe locations
! ******************************************************************************

  IF ( global%nProbes > 0 ) THEN 
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_FindProbeCells(pRegion)
    END DO ! iReg
    
    CALL RFLU_PrintProbeInfo(global)
  END IF ! global 

! ******************************************************************************
! Open files for convergence data, probe data, and OLES output
! ******************************************************************************

  CALL RFLU_OpenConverFile(global)
  
#ifndef GENX  
  IF ( moveGrid .EQV. .TRUE. ) THEN 
    CALL RFLU_OpenTotalMassFile(global)
  END IF ! moveGrid
#endif

  IF ( global%nProbes > 0 ) THEN 
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_OpenProbeFiles(pRegion)
    END DO ! iReg
  END IF ! global      

! BEGIN TEMPORARY - only works for single regions at the moment
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( pRegion%mixtInput%spaceDiscr == DISCR_OPT_LES ) THEN 
      CALL RFLU_OpenStatsFileOLES(global)
    END IF ! pRegion%mixtInput
  END DO ! iReg  
! END TEMPORARY

#ifndef GENX
! ******************************************************************************
! Compute integral values and write to file
! ******************************************************************************
      
  global%massIn  = 0.0_RFREAL
  global%massOut = 0.0_RFREAL    
      
  IF ( moveGrid .EQV. .TRUE. ) THEN
    CALL RFLU_ComputeIntegralValues(pLevel%regions) 
    CALL WriteTotalMass(pLevel%regions)
  END IF ! moveGrid
#endif  

#ifdef PETSC
! ******************************************************************************
! Initialize PETSc-related arrays and quantities
! ******************************************************************************
      
!  pRegion => pLevel%regions(1)
!  CALL RFLU_PETSC_Init(global)
!  CALL RFLU_PETSC_CreatePoisson(pRegion)
!  CALL RFLU_PETSC_BuildPoisson(pRegion)
!  CALL RFLU_PETSC_SetSolverPoisson(pRegion)
!  CALL RFLU_PETSC_SolvePressurePoisson(pRegion)
!  CALL RFLU_PETSC_Finalize(global)

  IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
    pRegion => pLevel%regions(1)
    CALL RFLU_PETSC_Init(global)
    CALL RFLU_PETSC_CreateVectors(pRegion)
    CALL RFLU_PETSC_CreateJacobian(pRegion)
  END IF ! global%solverType
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initialization done.'
    WRITE(STDOUT,'(A)') SOLVER_NAME 
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowSolver

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowSolver.F90,v $
! Revision 1.123  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.122  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.121  2007/04/20 16:07:49  mtcampbe
! Updating for burnout support function RFLU_GENX_InitBFLAG
!
! Revision 1.120  2007/04/14 14:33:35  mtcampbe
! Updated for TZ
!
! Revision 1.119  2007/02/27 13:12:29  haselbac
! Adapted to changes in RFLU_ModDimensionality
!
! Revision 1.118  2006/12/21 12:22:12  haselbac
! Added call to check patch-bc consistency
!
! Revision 1.117  2006/10/20 21:32:23  mparmar
! Added calls to create and compute global thrustFlags
!
! Revision 1.116  2006/08/19 15:46:27  mparmar
! Changed because of NSCBC implementation
!
! Revision 1.115  2006/08/18 14:04:50  haselbac
! Added calls to build AVFace2Patch list, read transforms, cosmetics
!
! Revision 1.114  2006/06/23 21:37:47  mtcampbe
! Added a GetBoundary call for genx restarts - makes surface data correct.
!
! Revision 1.113  2006/06/16 19:39:30  mtcampbe
! Fixed MPI_Init sequencing for Rocstar/Native initialization
!
! Revision 1.112  2006/06/06 21:36:04  mtcampbe
! Remeshing/Rocstar mods DummyCell population, time 0 inits
!
! Revision 1.111  2006/04/07 16:04:03  haselbac
! Adapted to changes in bf2c wts computation
!
! Revision 1.110  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.109  2006/04/07 14:53:27  haselbac
! Adapted to changes in bface stencil routines
!
! Revision 1.108  2006/03/25 22:03:26  haselbac
! Fix comment
!
! Revision 1.107  2006/03/25 22:02:51  haselbac
! Changes bcos of sype patches: Comment out Rxx2Sxx maps and face 
! reorientation, compute patch normals
!
! Revision 1.106  2006/03/09 14:10:03  haselbac
! Now call wrapper routines for stencils
!
! Revision 1.105  2006/02/06 23:55:54  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.104  2006/01/10 05:04:58  wasistho
! moved InitStatistics to before InitRocman
!
! Revision 1.103  2006/01/06 22:15:36  haselbac
! Adapted to name changes
!
! Revision 1.102  2006/01/03 06:32:03  wasistho
! moved Genx new-attr and registr. routines to wrappers
!
! Revision 1.101  2005/12/24 21:34:49  haselbac
! Added computation of ICT tolerance
!
! Revision 1.100  2005/12/24 02:25:31  wasistho
! moved statistics mapping befor Genx registration
!
! Revision 1.99  2005/12/24 02:15:22  wasistho
! activated statistics mapping
!
! Revision 1.98  2005/12/01 17:14:36  fnajjar
! Moved call to RFLU_RandomInit after reading user input
!
! Revision 1.97  2005/11/10 16:51:29  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.96  2005/11/04 14:07:24  haselbac
! Renamed dim check routine, added dim geom check routine
!
! Revision 1.95  2005/10/28 19:18:15  haselbac
! Added check for nProcs not being equal to nProcAlloc
!
! Revision 1.94  2005/10/27 19:20:19  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.93  2005/10/25 19:39:23  haselbac
! Added IF on forceFlag
!
! Revision 1.92  2005/10/05 14:19:23  haselbac
! Adapted to changes in stencil mods, added call to create and comp bface wts
!
! Revision 1.91  2005/09/14 15:59:33  haselbac
! Minor clean-up
!
! Revision 1.90  2005/09/13 20:39:39  mtcampbe
! Moved profiling call into this file
!
! Revision 1.89  2005/08/03 18:42:39  hdewey2
! Enclosed PETSc init calls inside IF
!
! Revision 1.88  2005/08/02 18:57:58  hdewey2
! Temporarily commented out USE RFLU_ModPETScPoisson
!
! Revision 1.87  2005/08/02 18:29:17  hdewey2
! Added init of PETSc data and quantities
!
! Revision 1.86  2005/07/01 15:15:17  haselbac
! Added setting of Roccom verbosity level
!
! Revision 1.85  2005/06/09 20:30:14  haselbac
! Disabled checking of move grid options
!
! Revision 1.84  2005/05/18 22:13:27  fnajjar
! ACH: Added creation of iPclSend buffers, now use nFacesAV
!
! Revision 1.83  2005/05/12 18:02:53  haselbac
! Removed call to COM_set_verbose
!
! Revision 1.82  2005/04/29 23:02:54  haselbac
! Added building of avf2b list
!
! Revision 1.81  2005/04/29 13:00:56  haselbac
! Fixed bug in name of routine for opening probe files
!
! Revision 1.80  2005/04/29 12:49:31  haselbac
! Added USE RFLU_ModProbes, removed interfaces for probe routines
!
! Revision 1.79  2005/04/15 16:31:18  haselbac
! Removed calls to XyzEdge2RegionDegrList routines
!
! Revision 1.78  2005/04/15 15:07:18  haselbac
! Converted to MPI, integrated MPI code with GENx
!
! Revision 1.77  2005/03/09 15:08:25  haselbac
! Added dimensionality check
!
! Revision 1.76  2005/01/18 15:18:18  haselbac
! Commented out COMM calls for now
!
! Revision 1.75  2005/01/14 21:35:29  haselbac
! Added calls to create and read comm lists
!
! Revision 1.74  2005/01/13 21:40:29  haselbac
! Bug fix in setting pRegion for PETSc testing
!
! Revision 1.73  2005/01/07 19:25:56  fnajjar
! Added call to PLAG_ReadSurfStatsWrapper
!
! Revision 1.72  2005/01/03 15:58:34  haselbac
! Adapted to changes in RFLU_ModStencils
!
! Revision 1.71  2004/12/28 20:28:13  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.70  2004/12/21 23:34:32  fnajjar
! Added definition of pMixtInput pointer in incompressible part
!
! Revision 1.69  2004/12/21 15:05:11  fnajjar
! Included calls for PLAG surface statistics
!
! Revision 1.68  2004/12/19 15:49:50  haselbac
! Added incompressible stuff
!
! Revision 1.67  2004/12/04 03:34:04  haselbac
! Adapted to changes in RFLU_ModCellMapping
!
! Revision 1.66  2004/11/29 17:17:20  wasistho
! use ModInterfacesStatistics
!
! Revision 1.65  2004/11/14 19:47:15  haselbac
! Now call RFLU_SetVarsWrapper instead of UpdateDependentVarsMP
!
! Revision 1.64  2004/11/03 17:05:08  haselbac
! Removed HACK_PERIODIC ifdef and call to RFLU_GENX_CreateAttrGridVol
!
! Revision 1.63  2004/11/02 02:32:51  haselbac
! Added call to RFLU_SetVarInfoWrapper
!
! Revision 1.62  2004/10/21 15:54:13  haselbac
! Added ifdef for static linking
!
! Revision 1.61  2004/10/20 15:01:13  haselbac
! Bug fix: Changes made to allow compilation within GENx without CHARM=1
!
! Revision 1.60  2004/10/19 19:29:19  haselbac
! Substantial changes because of new GENX logic
!
! Revision 1.59  2004/07/06 15:14:52  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.58  2004/06/25 20:08:25  haselbac
! Added call to RFLU_SetRestartTimeFlag
!
! Revision 1.57  2004/06/22 15:55:05  haselbac
! Bug fix for running serial jobs with Charm code
!
! Revision 1.56  2004/06/16 20:01:09  haselbac
! Added allocation of patch and force and moment coeffs, cosmetics
!
! Revision 1.55  2004/06/07 23:09:57  wasistho
! moved statistics mapping from initStatistics to before initGenxInterfaces
!
! Revision 1.54  2004/04/14 02:09:24  haselbac
! Added initialization of grid-speed factor
!
! Revision 1.53  2004/02/26 21:02:09  haselbac
! Added PLAG support, changed alloc logic, added updateDependentVarsMP
!
! Revision 1.52  2004/01/29 22:59:22  haselbac
! Added call to reading of bc data file, clean-up
!
! Revision 1.51  2004/01/22 16:04:33  haselbac
! Changed declaration to eliminate warning on ALC
!
! Revision 1.50  2003/12/07 04:59:22  jiao
! When GENX is defined, added a local variable "levels" to be
! used when calling RFLU_BuildDataStruct and RFLU_AssignRegionMapping
! to work around a complaint of PGI compiler.
!
! Revision 1.49  2003/12/04 03:30:04  haselbac
! Added calls for stencils and weights, cleaned up
!
! Revision 1.48  2003/11/25 21:04:41  haselbac
! Added call to RFLU_ReadBcInputFileWrapper, commented out communication
!
! Revision 1.47  2003/11/03 03:51:31  haselbac
! Added call to build boundary-face gradient access list
!
! Revision 1.46  2003/09/12 21:35:51  haselbac
! Fixed bug: Needed fem files for 1 proc cases
!
! Revision 1.45  2003/08/29 22:48:33  haselbac
! mpif.h should not be included directly here...
!
! Revision 1.44  2003/08/28 20:27:19  olawlor
! Minor tweaks for non-genx compilation--
!    - We need mpif.h
!    - Read input files from current directory, not Rocflu/inDir/.
!
! Revision 1.43  2003/08/27 15:38:00  haselbac
! Removed single and double quotes from last comments
!
! Revision 1.42  2003/08/26 22:48:56  olawlor
! Changes to startup sequence for latest Charm++/FEM framework:
!   - Call FEM_Init with an MPI communicator to set up FEM.
!     This replaces the old FEM_Attach call, and corresponding
!     -fem parameter in genx.C.
!
!   - Get partition number and number of partitions from MPI.
!     The FEM partition numbers are not needed any more.
!
!   - Read in the FEM input files yourself, using FEM_Mesh_read.
!     This allows the FEM input files to reside in the Rocflu
!     directory.
!
! Revision 1.41  2003/08/13 20:29:01  haselbac
! Fixed bug with writing probe data within GENx
!
! Revision 1.40  2003/07/22 15:39:50  haselbac
! Added Nullify routines
!
! Revision 1.39  2003/06/20 22:35:43  haselbac
! Added call to RFLU_ReadRestartInfo
!
! Revision 1.38  2003/06/04 20:05:53  jferry
! re-worked implementation of TBCs in unstructured code
!
! Revision 1.37  2003/05/16 21:51:29  mtcampbe
! ACH: Temporarily disabled checking of coupling input
!
! Revision 1.36  2003/05/13 23:49:46  haselbac
! Changed format for writing out number of procs
!
! Revision 1.35  2003/05/01 14:11:19  haselbac
! Added call to RFLU_CheckCouplingInput
!
! Revision 1.34  2003/04/24 15:43:13  haselbac
! Adapted interface to RFLU_PutBoundaryValues
!
! Revision 1.33  2003/04/12 21:38:08  haselbac
! Added setting of FEMRocfluGrid, verb and check level now read
!
! Revision 1.32  2003/04/07 14:26:50  haselbac
! Added cell-to-face list and probe info calls
!
! Revision 1.31  2003/03/31 16:17:21  haselbac
! Added writing of version string
!
! Revision 1.30  2003/03/25 19:16:58  haselbac
! Changed calling sequence for reorientation of AV faces
!
! Revision 1.29  2003/03/18 21:34:28  haselbac
! Modified call to RFLU_AllocateMemoryWrapper
!
! Revision 1.28  2003/03/15 18:38:12  haselbac
! Some adaptations for || gm, added serial mapping imposition
!
! Revision 1.27  2003/02/17 19:31:12  jferry
! Implemented portable random number generator ModRandom
!
! Revision 1.26  2003/01/31 13:53:57  haselbac
! Removed call to RFLU_GetBoundaryValues if not restarting
!
! Revision 1.25  2003/01/28 14:40:33  haselbac
! Extensive reorganization to get GENx restart to work with gm
!
! Revision 1.24  2002/11/27 20:25:14  haselbac
! Moved RFLU_SetMoveGridOptions to after having read bc file
!
! Revision 1.23  2002/11/26 15:27:55  haselbac
! Added RFLU_SetMoveGridOptions
!
! Revision 1.22  2002/11/26 15:18:39  haselbac
! Added ifdef to opening of total mass file
!
! Revision 1.21  2002/11/08 21:32:01  haselbac
! Added opening of total-mass file
!
! Revision 1.20  2002/10/27 19:13:03  haselbac
! New calls for grid motion, cosmetic redesign
!
! Revision 1.19  2002/10/19 16:15:19  haselbac
! Cosmetic changes to output
!
! Revision 1.18  2002/10/17 22:35:03  jiao
! ACH: Deactivate RFLU_PrintFlowInfo with GENX, gives small discrepancy for restart
!
! Revision 1.17  2002/10/17 14:14:00  haselbac
! Removed RFLU_GetBValues: Moved to RFLU_FlowSolver (discussion with Jim J.)
!
! Revision 1.16  2002/10/16 21:17:37  haselbac
! Added writing of header when running
!
! Revision 1.15  2002/10/12 14:59:13  haselbac
! Changed order of BC and GENX, added call to RFLU_GetBValues
!
! Revision 1.14  2002/10/05 19:22:53  haselbac
! GENX integration, added interface for RFLU_CreateFields
!
! Revision 1.13  2002/09/09 16:27:12  haselbac
! CVS appears to have got confused...
!
! Revision 1.12  2002/09/09 15:51:56  haselbac
! RFLU_HACK_PeriodicCellMapping.F90
!
! Revision 1.11  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.10  2002/07/25 14:20:51  haselbac
! Added call to OLES routines, completed CHARM code segment
!
! Revision 1.9  2002/06/27 15:28:22  haselbac
! Added CHARM stuff, deleted MPI stuff
!
! Revision 1.8  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.7  2002/06/14 20:21:26  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal
!
! Revision 1.6  2002/06/10 21:31:14  haselbac
! Added call to RFLU_PrintHeader, changed flag to CHECK_UNIFLOW
!
! Revision 1.5  2002/06/05 18:59:53  haselbac
! Added RFLU_PrintGridInfo and changed version number
!
! Revision 1.4  2002/05/04 17:10:20  haselbac
! Added memory allocation, mixture properties, and checking
!
! Revision 1.3  2002/04/11 19:03:45  haselbac
! Added calls and cosmetic changes
!
! Revision 1.2  2002/03/26 19:23:39  haselbac
! Some cleaning and added reading of input files
!
! Revision 1.1  2002/03/14 19:12:48  haselbac
! Initial revision
!
! ******************************************************************************







