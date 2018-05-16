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
! Purpose: Driver routine for rflupick.
!
! Description: None.
!
! Input: 
!   caseString	String with casename
!   stampString	String with iteration or time stamp
!   verbLevel	Verbosity level
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rflupick.F90,v 1.14 2008/12/06 08:45:04 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rflupick(caseString,stampString,verbLevel)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModMixture, ONLY: t_mixt_input
  USE ModBndPatch, ONLY: t_patch
  USE ModParameters
  USE ModMPI
  
  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModRegionMapping
  USE RFLU_ModStencilsBFaces
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsFaces
  USE RFLU_ModStencilsUtils
  USE RFLU_ModStencilsVert
  USE RFLU_ModVertexLists
  
  USE ModInterfaces, ONLY: RFLU_BuildDataStruct, &
                           RFLU_ClosePostInfo, &
                           RFLU_CreateGrid, &
                           RFLU_DecideNeedBGradFace, &
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, &
                           RFLU_OpenPostInfo, &
                           RFLU_PickRegionsCoord, &
                           RFLU_PickRegionsManual, &
                           RFLU_PickSpecialCells, &
                           RFLU_PickSpecialFaces, &
                           RFLU_PrintGridInfo, &
                           RFLU_PrintHeader, &
                           RFLU_PrintWarnInfo, &
                           RFLU_ReadInputFile, & 
                           RFLU_WritePostInfo, & 
                           RFLU_WriteVersionString
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString,stampString
  INTEGER, INTENT(IN) :: verbLevel
  
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: fileExists
  CHARACTER(CHRLEN) :: casename,choice,nRegions,RCSIdentString,stamp
  INTEGER :: errorFlag,iPatch,iReg
  TYPE(t_global), POINTER :: global
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_level), DIMENSION(:), POINTER :: levels  
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: rflupick.F90,v $ $Revision: 1.14 $'

! ******************************************************************************
! Initialize global data
! ******************************************************************************
  
  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag   

  casename = caseString(1:LEN(caseString))
  stamp    = stampString(1:LEN(stampString))
  
  CALL RFLU_InitGlobal(casename,verbLevel,CRAZY_VALUE_INT,global)

  CALL RegisterFunction(global,'rflupick', &
                        'rflupick.F90')

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
! Read mapping file and impose serial mapping
! ******************************************************************************

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_SetRegionMappingSerial(global)
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)

! ******************************************************************************
! Prepare data structure
! ******************************************************************************
  
  CALL RFLU_BuildDataStruct(global,levels)
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)  
  
! ****************************************************************************** 
! Read input file
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions)  

  IF ( global%flowType == FLOW_STEADY ) THEN
    READ(stamp,*) global%currentIter
  ELSE
    READ(stamp,*) global%currentTime
  END IF ! global%flowType  

! ****************************************************************************** 
! Get additional input
! ******************************************************************************

  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking regions based on bounding box?'
  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'n - No'
  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'y - Yes'    
  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter choice:'
  READ(STDIN,'(A)') choice
  
  SELECT CASE ( TRIM(choice) ) 
    CASE ( 'y' ) 
      global%pickCoordFlag = .TRUE.
      
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter coordinates of bounding box:'
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter x-coordinate range (low,high):'
      READ(STDIN,*) global%pickXCoordLow,global%pickXCoordUpp 
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter y-coordinate range (low,high):'
      READ(STDIN,*) global%pickYCoordLow,global%pickYCoordUpp 
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter z-coordinate range (low,high):'
      READ(STDIN,*) global%pickZCoordLow,global%pickZCoordUpp 
    CASE ( 'n' )
      global%pickCoordFlag = .FALSE.
    CASE DEFAULT
      global%warnCounter = global%warnCounter + 1
      
      WRITE(STDOUT,'(A,3X,A)')  SOLVER_NAME,'*** WARNING *** Invalid input.'
      WRITE(STDOUT,'(A,17X,A)') SOLVER_NAME,'Continuing assuming choice of no.'
      
      global%pickCoordFlag = .FALSE.
  END SELECT ! TRIM(choice)      

! ******************************************************************************
! Pick regions based on bounding box
! ******************************************************************************
    
  IF ( global%nRegionsLocal > 1 ) THEN 
    IF ( global%pickCoordFlag .EQV. .TRUE. ) THEN  
      DO iReg = 1,global%nRegionsLocal
        pRegion => levels(1)%regions(iReg)

        CALL RFLU_ReadDimensions(pRegion)         
        CALL RFLU_CreateGrid(pRegion)
        
        IF ( pRegion%grid%nPatches > 0 ) THEN      
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches        
        
        CALL RFLU_ReadGridWrapper(pRegion)
        
        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel

        CALL RFLU_PickRegionsCoord(pRegion)

        CALL RFLU_DestroyGrid(pRegion)     
      END DO ! iReg  
    END IF ! global%pickCoordFlag    
  END IF ! global%nRegionsLocal    
      
! ******************************************************************************
! Pick regions based on manual input
! ******************************************************************************

  IF ( global%nRegionsLocal > 1 ) THEN 
    CALL RFLU_PickRegionsManual(levels(1)%regions)
  END IF ! global%nRegions

! ******************************************************************************
! Open post-processing information file
! ******************************************************************************

  CALL RFLU_OpenPostInfo(global,FILE_STATUS_UNKNOWN,fileExists) 

! ******************************************************************************
! Pick special cells and write to file
! ******************************************************************************

  IF ( global%postSpecFlag .EQV. .TRUE. ) THEN 
    DO iReg = 1,global%nRegionsLocal
      IF ( global%nRegionsLocal /= 1 ) THEN 
        pRegion => levels(1)%regions(iReg)
      ELSE 
        pRegion => levels(1)%regions(0)
      END IF ! global%nRegionsLocal

      pMixtInput => pRegion%mixtInput
    
! ==============================================================================
!     If region active, proceed to allow user to pick cells
! ==============================================================================    
    
      IF ( pRegion%postActiveFlag .EQV. .TRUE. ) THEN
        CALL RFLU_ReadDimensions(pRegion) ! Must be done again
        CALL RFLU_CreateGrid(pRegion)
        
        IF ( pRegion%grid%nPatches > 0 ) THEN      
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches          
        
        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!       Build data structure
! ------------------------------------------------------------------------------

        CALL RFLU_CreateCellMapping(pRegion)
        CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
        CALL RFLU_BuildGlob2LocCellMapping(pRegion)
        
        CALL RFLU_CreateBVertexLists(pRegion)
        CALL RFLU_BuildBVertexLists(pRegion)
        
        CALL RFLU_CreateFaceList(pRegion)       
        CALL RFLU_BuildFaceList(pRegion)
        CALL RFLU_RenumberBFaceLists(pRegion) 

! ------------------------------------------------------------------------------
!       Build stencils. NOTE: Need geometry.
! ------------------------------------------------------------------------------

        CALL RFLU_CreateGeometry(pRegion)
        CALL RFLU_BuildGeometry(pRegion)

        CALL RFLU_CreateVert2CellList(pRegion)
        CALL RFLU_BuildVert2CellList(pRegion)

        CALL RFLU_CreateCell2FaceList(pRegion)
        CALL RFLU_BuildCell2FaceList(pRegion)

        IF ( pMixtInput%spaceOrder > 1 ) THEN 
          CALL RFLU_SetInfoC2CStencilWrapper(pRegion,pMixtInput%spaceOrder-1)        
          CALL RFLU_CreateC2CStencilWrapper(pRegion)
          CALL RFLU_BuildC2CStencilWrapper(pRegion)
        END IF ! pMixtInput%spaceOrder      

        IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
          CALL RFLU_SetInfoF2CStencilWrapper(pRegion,pMixtInput%spaceOrder-1)        
          CALL RFLU_CreateF2CStencilWrapper(pRegion)
          CALL RFLU_BuildF2CStencilWrapper(pRegion)
        END IF ! pMixtInput%flowModel 
          
        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)
          
          IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
            CALL RFLU_SetInfoBF2CStencilWrapper(pRegion,pPatch,pPatch%spaceOrder)
            CALL RFLU_CreateBF2CStencilWrapper(pRegion,pPatch)
            CALL RFLU_BuildBF2CStencilWrapper(pRegion,pPatch)          
          END IF ! RFLU_DecideNeedBGradFace

        END DO ! iPatch

        CALL RFLU_SetInfoStencilVert2Cell(pRegion,global%postInterpOrder)
        CALL RFLU_CreateStencilVert2Cell(pRegion)
        CALL RFLU_BuildStencilVert2Cell(pRegion)

        CALL RFLU_DestroyCell2FaceList(pRegion)
        CALL RFLU_DestroyVert2CellList(pRegion)

        CALL RFLU_DestroyGeometry(pRegion)    

! ------------------------------------------------------------------------------
!       Pick special cells and write to file
! ------------------------------------------------------------------------------

        CALL RFLU_PickSpecialCells(pRegion) 
        CALL RFLU_PickSpecialFaces(pRegion)          
        CALL RFLU_WritePostInfo(pRegion)

! ------------------------------------------------------------------------------
!       Deallocate memory
! ------------------------------------------------------------------------------
      
        CALL RFLU_DestroyStencilVert2Cell(pRegion)

        IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
          CALL RFLU_DestroyF2CStencilWrapper(pRegion)
        END IF ! pMixtInput%flowModel
          
        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)

          IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
            CALL RFLU_DestroyBF2CStencilWrapper(pRegion,pPatch)
          END IF ! RFLU_DecideNeedBGradFace
        END DO ! iPatch

        IF ( pMixtInput%spaceOrder > 1 ) THEN 
          CALL RFLU_DestroyC2CStencilWrapper(pRegion)
        END IF ! pMixtInput%spaceOrder      

        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)
        CALL RFLU_DestroyCellMapping(pRegion)              
      ELSE 
        CALL RFLU_WritePostInfo(pRegion)
      END IF ! pRegion%postActiveFlag    
    END DO ! iReg
  ELSE 
    DO iReg = 1,global%nRegionsLocal
      IF ( global%nRegionsLocal /= 1 ) THEN 
        pRegion => levels(1)%regions(iReg)
      ELSE 
        pRegion => levels(1)%regions(0)
      END IF ! global%nRegionsLocal  
  
      CALL RFLU_WritePostInfo(pRegion)
    END DO ! iReg      
  END IF ! global%postSpecFlag 

! ******************************************************************************
! Close post-processing information file
! ******************************************************************************

  CALL RFLU_ClosePostInfo(global)
  
! *****************************************************************************
! Print info about warnings
! *****************************************************************************

  CALL RFLU_PrintWarnInfo(global)
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE rflupick

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rflupick.F90,v $
! Revision 1.14  2008/12/06 08:45:04  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:18:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2006/08/19 15:41:15  mparmar
! Used pPatch%spaceOrder, RFLU_NSCBC_DecideNeedBGradFace
!
! Revision 1.11  2006/04/07 14:57:17  haselbac
! Adapted to changes in bface stencil routines
!
! Revision 1.10  2006/03/09 14:11:06  haselbac
! Now call wrapper routine for F2C stencils
!
! Revision 1.9  2006/02/06 23:55:55  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.8  2006/01/06 22:18:36  haselbac
! Adapted to name changes
!
! Revision 1.7  2005/12/10 23:30:23  haselbac
! Added user input for bbox
!
! Revision 1.6  2005/10/27 19:21:49  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.5  2005/10/09 15:13:29  haselbac
! Bug fix: Added bc reading, needed for bface stencils
!
! Revision 1.4  2005/10/05 14:26:13  haselbac
! Adapted to changes in stencil modules, added use of vertex list module
!
! Revision 1.3  2005/07/19 19:17:48  haselbac
! Bug fix: Added calls to build c2f list for picking stencils
!
! Revision 1.2  2005/05/03 03:11:54  haselbac
! Converted to C++ reading of command-line
!
! Revision 1.1  2005/04/18 14:57:56  haselbac
! Initial revision
!
! ******************************************************************************







