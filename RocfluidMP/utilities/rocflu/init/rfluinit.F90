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
! Purpose: Driver routine for rfluinit. 
!
! Description: None.
!
! Input: 
!   caseString  String with casename
!   verbLevel   Verbosity level
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rfluinit.F90,v 1.21 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rfluinit(caseString,verbLevel)

  USE ModError
  USE ModDataTypes
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModMPI
#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
#endif
  
  USE RFLU_ModAllocateMemory
  USE RFLU_ModBoundLists
  USE RFLU_ModBoundXvUtils 
  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
  USE RFLU_ModCellMapping
  USE RFLU_ModDeallocateMemory
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModGridSpeedUtils
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteBcDataFile
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModRegionMapping
  USE RFLU_ModRenumberings

#ifdef GENX
  USE RFLU_ModRocstarAdmin
  USE RFLU_ModRocstarIO
#endif

  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, & 
                           RFLU_BuildDataStruct, &   
                           RFLU_CreateGrid, & 
                           RFLU_DeallocMemSolWrapper, & 
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, &
                           RFLU_InitBcDataHardCode, & 
                           RFLU_InitFlowHardCodeLimWrapper, &  
                           RFLU_InitFlowHardCodeWrapper, & 
                           RFLU_InitFlowScratchWrapper, &
                           RFLU_InitFlowSerialWrapper, & 
                           RFLU_PrintFlowInfoWrapper, & 
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintHeader, & 
                           RFLU_PrintWarnInfo, &
                           RFLU_RandomInit, &
                           RFLU_ReadRestartInfo, & 
                           RFLU_SetDependentVars, &
                           RFLU_SetModuleType, &  
                           RFLU_SetRestartTimeFlag, &
                           RFLU_SetVarInfoWrapper, & 
                           RFLU_WriteVersionString, &
                           ScaleRotateVector
                           
#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_SetGasVars
#endif                  

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_InitPatchData, &
                                     PLAG_RFLU_AllocMemSol, & 
                                     PLAG_RFLU_AllocMemSolTile, & 
                                     PLAG_RFLU_DeallocMemSol, & 
                                     PLAG_RFLU_DeallocMemSolTile, &
                                     PLAG_RFLU_InitSolutionRandom, &
                                     PLAG_RFLU_InitSolutionScratch, &
                                     PLAG_RFLU_InitSolFromSerial, &
                                     PLAG_RFLU_InitSolSerialWrapper, & 
                                     PLAG_RFLU_WriteSolutionASCII, &
                                     PLAG_RFLU_WriteSolutionBinary

  USE PLAG_ModDataStruct
  USE PLAG_ModDimensions,      ONLY: PLAG_RFLU_WriteDimensions, & 
                                     PLAG_SetDimensions, &
                                     PLAG_SetMaxDimensions
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString
  INTEGER, INTENT(IN) :: verbLevel

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: casename
#ifdef GENX
  CHARACTER(CHRLEN) :: surfWinNameInput,volWinNameInput
#endif
  INTEGER :: errorFlag,iLev,iReg,iRegLow,iRegUpp
#ifdef GENX
  INTEGER :: handleObtain
#endif
  TYPE(t_region), POINTER :: pRegion,pRegionSerial
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)
#ifdef PLAG
  TYPE(t_plag), POINTER :: pPlag,pPlagSerial
#endif

! ******************************************************************************
! Initialize global data
! ******************************************************************************  
  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag 

  casename = caseString(1:LEN(caseString))
  
  CALL RFLU_InitGlobal(casename,verbLevel,MPI_COMM_WORLD,global)

  CALL RegisterFunction(global,'rfluinit', & 
                        'rfluinit.F90')

  CALL RFLU_SetModuleType(global,MODULE_TYPE_INIT)

#ifdef GENX
! ******************************************************************************
! Read GENX control file, store communicator and hardcode window name
! ****************************************************************************** 

  CALL RFLU_GENX_ReadCtrlFile(global)
  CALL RFLU_GENX_StoreCommunicator(global,MPI_COMM_WORLD)
  CALL RFLU_GENX_HardCodeWindowName(global)
#endif

! ******************************************************************************
! Print header and write version string
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL RFLU_WriteVersionString(global)     
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      CALL RFLU_PrintHeader(global)
    END IF ! global%verbLevel
  END IF ! global%myProcid

! ******************************************************************************
! Read mapping file, impose serial mapping, and build basic data structure
! ******************************************************************************

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_SetRegionMappingSerial(global)  
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)

  CALL RFLU_BuildDataStruct(global,levels) 
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)  

#ifdef GENX
! ******************************************************************************
! Initialize Roccom and load Rocin and Rocout
! ******************************************************************************

  CALL COM_init
! TEMPORARY
!  CALL COM_set_verbose(10)
! END TEMPORARY

  surfWinNameInput = 'RocfluInputSurf'
  volWinNameInput  = 'RocfluInputVol TEMP'

  global%winNameIn  = TRIM(global%winName)//'-IN'
  global%winNameOut = TRIM(global%winName)//'-OUT'

  CALL SimIN_load_module(TRIM(global%winNameIn))
  CALL SimOUT_load_module(TRIM(global%winNameOut))
  
  handleObtain = COM_get_function_handle(TRIM(global%winNameIn)// &
                                         '.obtain_dataitem')
  CALL RFLU_GENX_StoreNamesHandles(global,surfWinNameInput,volWinNameInput, & 
                                   handleObtain)
  
! ******************************************************************************
! Create windows and dataitems 
! ******************************************************************************
  
  pRegionSerial => levels(1)%regions(0)

  CALL RFLU_GENX_CreateWindows(pRegionSerial) 
  CALL RFLU_GENX_CreateAttrGridSurf(pRegionSerial)
  CALL RFLU_GENX_CreateAttrFlow(pRegionSerial)    
  CALL RFLU_GENX_CreateAttrGSpeeds(pRegionSerial)
  CALL RFLU_GENX_CreateAttrInterf(pRegionSerial)
#endif

! ******************************************************************************
! Initialize random number generator. NOTE needed in order to write sensible 
! data when writing Rocpart solution files. 
! ******************************************************************************

  CALL RFLU_RandomInit(levels(1)%regions)

! ******************************************************************************
! Read input file and restart info. NOTE need restart info for GENX runs to 
! determine whether have a restart. 
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions,.TRUE.) 
  CALL RFLU_ReadRestartInfo(global)
  CALL RFLU_SetRestartTimeFlag(global)

! ******************************************************************************
! Initialize solutions
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN 
    iRegLow = 0
    iRegUpp = 0
  ELSE 
    iRegLow = 1
    iRegUpp = global%nRegions
  END IF ! global%nRegions

#ifdef PLAG
! ==============================================================================
! Temporarily disable particles so can use wrapper routines to initialize only
! Eulerian solution fields
! ==============================================================================

  global%plagUsedSave = global%plagUsed

  global%plagUsed = .FALSE.
#endif

! ==============================================================================
! Initialize Eulerian solution fields
! ==============================================================================

  SELECT CASE ( global%initFlowFlag ) 
    
! ------------------------------------------------------------------------------
!   Initialize from scratch
! ------------------------------------------------------------------------------

    CASE ( INITFLOW_FROMSCRATCH )
      DO iReg = iRegLow,iRegUpp
        pRegion => levels(1)%regions(iReg)

        CALL RFLU_ReadDimensions(pRegion)                          
        CALL RFLU_CreateGrid(pRegion)

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches                   

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_ReadGridWrapper(pRegion)

          CALL RFLU_CreateCellMapping(pRegion)
          CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
          CALL RFLU_BuildGlob2LocCellMapping(pRegion)         

          CALL RFLU_CreateBVertexLists(pRegion)
          CALL RFLU_BuildBVertexLists(pRegion)

          CALL RFLU_CreateFaceList(pRegion)
          CALL RFLU_BuildFaceList(pRegion)
          CALL RFLU_RenumberBFaceLists(pRegion)

          CALL RFLU_CreateGeometry(pRegion)
          CALL RFLU_BuildGeometry(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

        CALL RFLU_AllocMemSolWrapper(pRegion)  
        CALL RFLU_SetVarInfoWrapper(pRegion) 

        CALL RFLU_InitFlowScratchWrapper(pRegion)
        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)

! Initializing boundary array would need solution in domain to be initialized first
        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_CreateVarsCv(pRegion)
          CALL RFLU_BXV_CreateVarsDv(pRegion)
          CALL RFLU_BXV_InitVars(pRegion)
          CALL RFLU_BXV_WriteVarsWrapper(pRegion)
          CALL RFLU_BXV_DestroyVarsCv(pRegion)
          CALL RFLU_BXV_DestroyVarsDv(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

        IF ( global%verbLevel > VERBOSE_NONE ) THEN   
          CALL RFLU_PrintFlowInfoWrapper(pRegion)    
        END IF ! global%verbLevel       

#ifdef GENX
        CALL RFLU_GENX_SetConnSize(pRegion)
        CALL RFLU_GENX_RegisterDataFlow(pRegion)               
#endif
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_DestroyFaceList(pRegion)
          CALL RFLU_DestroyBVertexLists(pRegion)    
          CALL RFLU_DestroyCellMapping(pRegion)
          CALL RFLU_DestroyGeometry(pRegion)  
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

        CALL RFLU_DeallocMemSolWrapper(pRegion)
        CALL RFLU_DestroyGrid(pRegion)
      END DO ! iReg  

! ------------------------------------------------------------------------------
!   Initialize parallel run by reading solution from serial file.
! ------------------------------------------------------------------------------
  
    CASE ( INITFLOW_FROMFILE ) 
      IF ( global%nRegions > 1 ) THEN 
        pRegionSerial => levels(1)%regions(0)

! ----- Read serial solution ---------------------------------------------------

        CALL RFLU_ReadDimensions(pRegionSerial)               
        CALL RFLU_CreateGrid(pRegionSerial)

        IF ( pRegionSerial%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)    
        END IF ! pRegionSerial%grid%nPatches         

        CALL RFLU_AllocMemSolWrapper(pRegionSerial)  
        CALL RFLU_SetVarInfoWrapper(pRegionSerial)     

        CALL RFLU_ReadFlowWrapper(pRegionSerial)

! ----- Loop over regions and initialize ---------------------------------------

        DO iReg = 1,global%nRegions
          pRegion => levels(1)%regions(iReg)

          CALL RFLU_ReadDimensionsWrapper(pRegion)               
          CALL RFLU_CreateGrid(pRegion)

          IF ( pRegion%grid%nPatches > 0 ) THEN        
            CALL RFLU_ReadBCInputFileWrapper(pRegion)    
          END IF ! pRegion%grid%nPatches         

          CALL RFLU_AllocMemSolWrapper(pRegion)  
          CALL RFLU_SetVarInfoWrapper(pRegion)

          CALL RFLU_RNMB_CreatePC2SCMap(pRegion)    
          CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
          CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)

          CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)

          CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
          CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)        

          CALL RFLU_InitFlowSerialWrapper(pRegion,pRegionSerial)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN   
            CALL RFLU_PrintFlowInfoWrapper(pRegion)    
          END IF ! global%verbLevel  

#ifdef GENX
          CALL RFLU_GENX_RegisterDataFlow(pRegion)
          CALL RFLU_SetGasVars(pRegion,1,pRegion%grid%nCellsTot)               
          CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
          CALL RFLU_WriteFlowWrapper(pRegion)

          CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)               
          CALL RFLU_DeallocMemSolWrapper(pRegion) 
          CALL RFLU_DestroyGrid(pRegion) 
        END DO ! iReg

! ----- Deallocate memory ------------------------------------------------------

        CALL RFLU_DeallocMemSolWrapper(pRegionSerial) 
        CALL RFLU_DestroyGrid(pRegionSerial)    
      END IF ! global%nRegions

! ------------------------------------------------------------------------------
!   Initialize from hardcode
! ------------------------------------------------------------------------------

    CASE ( INITFLOW_FROMHARDCODE ) 
      DO iReg = iRegLow,iRegUpp
        pRegion => levels(1)%regions(iReg)

        CALL RFLU_ReadDimensions(pRegion)               
        CALL RFLU_CreateGrid(pRegion)

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches       

#ifdef GENX
        CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_SURF)
        CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_VOL)
        CALL RFLU_GENX_GetDimensionsDerived(pRegion)
        CALL RFLU_GENX_CreateGridSurf(pRegion) 
        CALL RFLU_GENX_RegisterGridSurf(pRegion) 
        CALL RFLU_GENX_RegisterGridVol(pRegion)    
#endif

        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

#ifdef GENX
        CALL RFLU_GENX_DestroyGridSurf(pRegion) 
#endif

        CALL RFLU_CreateCellMapping(pRegion)
        CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
        CALL RFLU_BuildGlob2LocCellMapping(pRegion)         

        CALL RFLU_CreateBVertexLists(pRegion)
        CALL RFLU_BuildBVertexLists(pRegion)

        CALL RFLU_CreateFaceList(pRegion)
        CALL RFLU_BuildFaceList(pRegion)
        CALL RFLU_RenumberBFaceLists(pRegion)

        CALL RFLU_CreateGeometry(pRegion)    
        CALL RFLU_BuildGeometry(pRegion)   

        CALL RFLU_AllocMemSolWrapper(pRegion)  
        CALL RFLU_SetVarInfoWrapper(pRegion) 
        
        CALL RFLU_InitFlowHardCodeWrapper(pRegion)

        IF ( RFLU_DecideReadWriteBcDataFile(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_InitBcDataHardCode(pRegion)      
        END IF ! RFLU_DecideReadWriteBcDataFile      

        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_CreateVarsCv(pRegion)
          CALL RFLU_BXV_CreateVarsDv(pRegion)
          CALL RFLU_BXV_InitVars(pRegion)
          CALL RFLU_BXV_WriteVarsWrapper(pRegion)
          CALL RFLU_BXV_DestroyVarsCv(pRegion)
          CALL RFLU_BXV_DestroyVarsDv(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

        IF ( global%verbLevel > VERBOSE_NONE ) THEN   
          CALL RFLU_PrintFlowInfoWrapper(pRegion)    
        END IF ! global%verbLevel       

#ifdef GENX
        CALL RFLU_GENX_RegisterDataFlow(pRegion)              
        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
        CALL RFLU_WriteFlowWrapper(pRegion)       

        IF ( RFLU_DecideReadWriteBcDataFile(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_WriteBcDataFile(pRegion)
        END IF ! RFLU_DecideReadWriteBcDataFile

        CALL RFLU_DestroyGeometry(pRegion)  
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)    
        CALL RFLU_DestroyCellMapping(pRegion)
        CALL RFLU_DeallocMemSolWrapper(pRegion)
        CALL RFLU_DestroyGrid(pRegion)                                  
      END DO ! iReg   

! ------------------------------------------------------------------------------
!   Initialize parallel run from combo using serial solution
! ------------------------------------------------------------------------------
  
    CASE ( INITFLOW_FROMCOMBO_SERIAL ) 
      IF ( global%nRegions > 1 ) THEN 
        pRegionSerial => levels(1)%regions(0)

! ----- Read serial solution ---------------------------------------------------

        CALL RFLU_ReadDimensions(pRegionSerial)               
        CALL RFLU_CreateGrid(pRegionSerial)

        IF ( pRegionSerial%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)    
        END IF ! pRegionSerial%grid%nPatches         

        CALL RFLU_AllocMemSolWrapper(pRegionSerial)  
        CALL RFLU_SetVarInfoWrapper(pRegionSerial)     

        CALL RFLU_ReadFlowWrapper(pRegionSerial)

! ----- Loop over regions and initialize ---------------------------------------

        DO iReg = 1,global%nRegions
          pRegion => levels(1)%regions(iReg)

          CALL RFLU_ReadDimensionsWrapper(pRegion)               
          CALL RFLU_CreateGrid(pRegion)

          IF ( pRegion%grid%nPatches > 0 ) THEN        
            CALL RFLU_ReadBCInputFileWrapper(pRegion)    
          END IF ! pRegion%grid%nPatches         

          CALL RFLU_ReadGridWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_LOW ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          CALL RFLU_AllocMemSolWrapper(pRegion)  
          CALL RFLU_SetVarInfoWrapper(pRegion)

          CALL RFLU_RNMB_CreatePC2SCMap(pRegion)    
          CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
          CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)

          CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)

          CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
          CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)        

          CALL RFLU_InitFlowSerialWrapper(pRegion,pRegionSerial)
          
          CALL RFLU_CreateCellMapping(pRegion)
          CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
          CALL RFLU_BuildGlob2LocCellMapping(pRegion)         

          CALL RFLU_CreateBVertexLists(pRegion)
          CALL RFLU_BuildBVertexLists(pRegion)

          CALL RFLU_CreateFaceList(pRegion)
          CALL RFLU_BuildFaceList(pRegion)
          CALL RFLU_RenumberBFaceLists(pRegion)

          CALL RFLU_CreateGeometry(pRegion)    
          CALL RFLU_BuildGeometry(pRegion)
                  
          CALL RFLU_InitFlowHardCodeLimWrapper(pRegion)

          CALL RFLU_DestroyGeometry(pRegion)  
          CALL RFLU_DestroyFaceList(pRegion)
          CALL RFLU_DestroyBVertexLists(pRegion)    
          CALL RFLU_DestroyCellMapping(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN   
            CALL RFLU_PrintFlowInfoWrapper(pRegion)    
          END IF ! global%verbLevel  

#ifdef GENX
          CALL RFLU_GENX_RegisterDataFlow(pRegion)               
          CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
          CALL RFLU_WriteFlowWrapper(pRegion)

          CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)               
          CALL RFLU_DeallocMemSolWrapper(pRegion) 
          CALL RFLU_DestroyGrid(pRegion) 
        END DO ! iReg

! ----- Deallocate memory ------------------------------------------------------

        CALL RFLU_DeallocMemSolWrapper(pRegionSerial) 
        CALL RFLU_DestroyGrid(pRegionSerial)    
      END IF ! global%nRegions

! ------------------------------------------------------------------------------
!   Initialize parallel run from combo using parallel solution
! ------------------------------------------------------------------------------
  
    CASE ( INITFLOW_FROMCOMBO_PARALLEL ) 

! --- Loop over regions and initialize -----------------------------------------

      DO iReg = 1,global%nRegions
        pRegion => levels(1)%regions(iReg)

        CALL RFLU_ReadDimensionsWrapper(pRegion)               
        CALL RFLU_CreateGrid(pRegion)

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches         

        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        CALL RFLU_AllocMemSolWrapper(pRegion)  
        CALL RFLU_SetVarInfoWrapper(pRegion)

        CALL RFLU_ReadFlowWrapper(pRegion)

        CALL RFLU_CreateCellMapping(pRegion)
        CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
        CALL RFLU_BuildGlob2LocCellMapping(pRegion)         

        CALL RFLU_CreateBVertexLists(pRegion)
        CALL RFLU_BuildBVertexLists(pRegion)

        CALL RFLU_CreateFaceList(pRegion)
        CALL RFLU_BuildFaceList(pRegion)
        CALL RFLU_RenumberBFaceLists(pRegion)

        CALL RFLU_CreateGeometry(pRegion)    
        CALL RFLU_BuildGeometry(pRegion)

        CALL RFLU_InitFlowHardCodeLimWrapper(pRegion)

        CALL RFLU_DestroyGeometry(pRegion)  
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)    
        CALL RFLU_DestroyCellMapping(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN   
          CALL RFLU_PrintFlowInfoWrapper(pRegion)    
        END IF ! global%verbLevel  

#ifdef GENX
        CALL RFLU_GENX_RegisterDataFlow(pRegion)               
        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
        CALL RFLU_WriteFlowWrapper(pRegion)
      
        CALL RFLU_DeallocMemSolWrapper(pRegion) 
        CALL RFLU_DestroyGrid(pRegion) 
      END DO ! iReg    

! ------------------------------------------------------------------------------
!   Default
! ------------------------------------------------------------------------------
  
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%initFlowFlag  

#ifdef PLAG
! ==============================================================================
! Re-enable particles
! ==============================================================================

  global%plagUsed = global%plagUsedSave

! ==============================================================================
! Initialize particle solution 
! ==============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    pRegionSerial => levels(1)%regions(0)
    pPlagSerial   => pRegionSerial%plag

    CALL RFLU_ReadDimensions(pRegionSerial)               
    CALL RFLU_CreateGrid(pRegionSerial)       

    IF ( pRegionSerial%grid%nPatches > 0 ) THEN        
      CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)    
    END IF ! pRegionSerial%grid%nPatches
    
! ------------------------------------------------------------------------------
!   Initialize particle data
! ------------------------------------------------------------------------------

    CALL PLAG_SetDimensions(pRegionSerial,pRegionSerial%plagInput%nPclsIni)
    CALL PLAG_SetMaxDimensions(pRegionSerial)
    CALL PLAG_RFLU_AllocMemSol(pRegionSerial,pPlagSerial)
    CALL PLAG_RFLU_AllocMemSolTile(pRegionSerial)

    SELECT CASE ( global%initPlagFlag )
      CASE ( PLAG_INIT_FROMSCRATCH ) 
        CALL PLAG_RFLU_InitSolutionScratch(pRegionSerial)
      CASE ( PLAG_INIT_FROMRANDOMSTATE )
        CALL PLAG_RFLU_InitSolutionRandom(pRegionSerial)
      CASE ( PLAG_INIT_FROMFILE ) 
! TO DO 
! Read particle dimension file
! Read solution from file - needed for restart on diff number of procs
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
! END TO DO 
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! global%initPlagFlag

    CALL RFLU_ReadGridWrapper(pRegionSerial)    
    
    CALL RFLU_CreateCellMapping(pRegionSerial)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegionSerial)
    CALL RFLU_BuildGlob2LocCellMapping(pRegionSerial)         

    CALL RFLU_CreateBVertexLists(pRegionSerial)
    CALL RFLU_BuildBVertexLists(pRegionSerial)

    CALL RFLU_CreateFaceList(pRegionSerial)
    CALL RFLU_BuildFaceList(pRegionSerial)
    CALL RFLU_RenumberBFaceLists(pRegionSerial)

    CALL RFLU_CreateGeometry(pRegionSerial)    
    CALL RFLU_BuildGeometry(pRegionSerial)    
    
    CALL PLAG_RFLU_InitSolSerialWrapper(pRegionSerial)
    CALL PLAG_InitPatchData(pRegionSerial)

    CALL PLAG_RFLU_WriteDimensions(pRegionSerial)
 
    IF ( global%nRegions == 1 ) THEN 
! TO DO  Call wrapper routine eventually
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL PLAG_RFLU_WriteSolutionASCII(pRegionSerial)
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL PLAG_RFLU_WriteSolutionBinary(pRegionSerial)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
! END TO DO Call wrapper routine eventually

      CALL RFLU_DestroyGeometry(pRegionSerial)
      CALL RFLU_DestroyFaceList(pRegionSerial) 
      CALL RFLU_DestroyBVertexLists(pRegionSerial)
      CALL RFLU_DestroyCellMapping(pRegionSerial) 
    END IF ! global%nRegions

! ------------------------------------------------------------------------------
!   Initialize parallel regions and write solution files
! ------------------------------------------------------------------------------

    IF ( global%nRegions > 1 ) THEN
      CALL PLAG_DSTR_CreatePclListCSR(pRegionSerial)
      CALL PLAG_DSTR_BuildCell2PclList(pRegionSerial)

      DO iReg = 1,global%nRegions 
        pRegion => levels(1)%regions(iReg)
        pPlag   => pRegion%plag

        CALL RFLU_ReadDimensions(pRegion)               
        CALL RFLU_CreateGrid(pRegion)       

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches

        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        CALL PLAG_SetDimensions(pRegion,pRegion%plagInput%nPclsIni)
        CALL PLAG_SetMaxDimensions(pRegion)
        CALL PLAG_RFLU_AllocMemSol(pRegion,pPlag)
        CALL PLAG_RFLU_AllocMemSolTile(pRegion)

        CALL PLAG_RFLU_InitSolFromSerial(pRegion,pRegionSerial)
        CALL PLAG_InitPatchData(pRegion)
  
        CALL PLAG_SetMaxDimensions(pRegion)
        CALL PLAG_RFLU_WriteDimensions(pRegion)

! TO DO Call wrapper routine eventually
        IF ( global%solutFormat == FORMAT_ASCII ) THEN
          CALL PLAG_RFLU_WriteSolutionASCII(pRegion)
        ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
          CALL PLAG_RFLU_WriteSolutionBinary(pRegion)
        ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! global%solutFormat
! END TO DO Call wrapper routine eventually

        CALL PLAG_RFLU_DeallocMemSolTile(pRegion)
        CALL PLAG_RFLU_DeallocMemSol(pRegion,pPlag)

        CALL RFLU_DestroyGrid(pRegion) 
      END DO ! iReg

      CALL PLAG_DSTR_DestroyPclListCSR(pRegionSerial)
      CALL PLAG_DSTR_DestroyCell2PclList(pRegionSerial)

! TO DO      
!      CALL PLAG_RFLU_GetPclsSizeSumLocal(regions)
!      CALL PLAG_RFLU_CheckPclsSize(regions)
! END TO DO 
    END IF ! global%nRegions

! ------------------------------------------------------------------------------
!   Deallocate memory for serial region
! ------------------------------------------------------------------------------

    CALL PLAG_RFLU_DeallocMemSolTile(pRegionSerial)
    CALL PLAG_RFLU_DeallocMemSol(pRegionSerial,pPlagSerial)
       
    CALL RFLU_DestroyGrid(pRegionSerial)     
  END IF ! global%plagUsed
#endif  

! ******************************************************************************
! Write grid speed files. NOTE separated from above because need number of 
! faces, which is not always known above.
! ******************************************************************************

  DO iReg = iRegLow,iRegUpp
    pRegion => levels(1)%regions(iReg)

    IF ( RFLU_DecideNeedGridSpeeds(pRegion) .EQV. .TRUE. ) THEN 
      CALL RFLU_ReadDimensions(pRegion)                          
      CALL RFLU_CreateGrid(pRegion)

      IF ( pRegion%grid%nPatches > 0 ) THEN        
        CALL RFLU_ReadBCInputFileWrapper(pRegion)    
      END IF ! pRegion%grid%nPatches                   

#ifdef GENX
      CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_SURF)
      CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_VOL)
      CALL RFLU_GENX_GetDimensionsDerived(pRegion)
      CALL RFLU_GENX_CreateGridSurf(pRegion) 
      CALL RFLU_GENX_RegisterGridSurf(pRegion) 
      CALL RFLU_GENX_RegisterGridVol(pRegion)    
#endif

      CALL RFLU_ReadGridWrapper(pRegion)

#ifdef GENX
      CALL RFLU_GENX_DestroyGridSurf(pRegion) 
#endif

      CALL RFLU_CreateCellMapping(pRegion)
      CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
      CALL RFLU_BuildGlob2LocCellMapping(pRegion)         

      CALL RFLU_CreateBVertexLists(pRegion)
      CALL RFLU_BuildBVertexLists(pRegion)

      CALL RFLU_CreateFaceList(pRegion)
      CALL RFLU_BuildFaceList(pRegion)
      CALL RFLU_RenumberBFaceLists(pRegion)

      CALL RFLU_AllocateMemoryGSpeeds(pRegion)
#ifdef GENX
      CALL RFLU_GENX_RegisterDataGSpeeds(pRegion) 
#endif
      CALL RFLU_WriteGridSpeedsWrapper(pRegion)

      CALL RFLU_DeallocateMemoryGSpeeds(pRegion) 

      CALL RFLU_DestroyFaceList(pRegion)
      CALL RFLU_DestroyBVertexLists(pRegion)    
      CALL RFLU_DestroyCellMapping(pRegion) 
      CALL RFLU_DestroyGrid(pRegion)  

#ifdef GENX
      CALL COM_delete_window(TRIM(global%volWinNameInput))
      CALL COM_delete_window(TRIM(global%surfWinNameInput))
#endif
    END IF ! RFLU_DecideNeedGridSpeeds 
  END DO ! iReg  

! ******************************************************************************
! Print info about warnings
! ******************************************************************************
 
  CALL RFLU_PrintWarnInfo(global)
                              
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)
  WRITE (*,*) "rfluinit finished sucessfully"

END SUBROUTINE rfluinit

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rfluinit.F90,v $
! Revision 1.21  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.20  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.19  2007/03/27 01:31:45  haselbac
! Remove superfluous USE PLAG_ModParameters statement (bad check-in)
!
! Revision 1.18  2007/03/27 00:46:14  haselbac
! Adapted to changes in RFLU_SetDimensions call
!
! Revision 1.17  2007/03/27 00:23:23  haselbac
! PLAG init completely revamped to speed up 1d cases substantially
!
! Revision 1.16  2007/03/20 17:35:16  fnajjar
! Modified USE call to streamline with new module PLAG_ModDimensions
!
! Revision 1.15  2007/03/15 22:00:58  haselbac
! Adapted to changes in PLAG init for serial runs
!
! Revision 1.14  2006/08/19 15:41:13  mparmar
! Added calls to create, init, write, and destroy patch arrays
!
! Revision 1.13  2006/05/05 18:23:47  haselbac
! Changed PLAG init so do not need serial region anymore
!
! Revision 1.12  2006/02/06 23:55:55  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.11  2005/11/10 02:44:45  haselbac
! Added support for variable properties
!
! Revision 1.10  2005/09/23 19:00:45  haselbac
! Bug fix: When init PLAG, did not know about bc
!
! Revision 1.9  2005/09/13 21:37:30  haselbac
! Added new init option
!
! Revision 1.8  2005/05/18 22:23:59  fnajjar
! Added capability of init particles
!
! Revision 1.7  2005/05/05 18:38:39  haselbac
! Removed MPI calls after bug in Rocin/out fixed
!
! Revision 1.6  2005/05/04 03:37:50  haselbac
! Commented out COM_set_verbose call
!
! Revision 1.5  2005/05/04 03:35:58  haselbac
! Added init and finalize MPI when running within GENX
!
! Revision 1.4  2005/05/03 03:10:06  haselbac
! Converted to C++ reading of command-line
!
! Revision 1.3  2005/04/22 15:20:24  haselbac
! Fixed bug in combo init: grid and geom was missing; added grid info calls
!
! Revision 1.2  2005/04/18 20:33:27  haselbac
! Removed USE RFLU_ModCommLists
!
! Revision 1.1  2005/04/15 15:08:15  haselbac
! Initial revision
!
! ******************************************************************************







