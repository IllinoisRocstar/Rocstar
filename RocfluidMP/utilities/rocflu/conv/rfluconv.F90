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
! Purpose: Driver routine for rfluconv.
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
! $Id: rfluconv.F90,v 1.6 2008/12/06 08:44:55 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rfluconv(caseString,stampString,verbLevel)

  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModParameters
           
  USE RFLU_ModBoundLists             
  USE RFLU_ModCellMapping
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList 
  USE RFLU_ModGeometry 
  USE RFLU_ModReadBcInputFile 
  USE RFLU_ModReadWriteFlow 
  USE RFLU_ModReadWriteGrid 
  USE RFLU_ModRegionMapping        
  USE RFLU_ModSTL
  USE RFLU_ModTETMESH
               
  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, & 
                           RFLU_BuildDataStruct, &
                           RFLU_CreateGrid, &
                           RFLU_DeallocMemSolWrapper, &
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, &
                           RFLU_PrintFlowInfo, &
                           RFLU_PrintGridInfo, &
                           RFLU_PrintHeader, &
                           RFLU_PrintWarnInfo, &
                           RFLU_SetVarInfoWrapper, &   
                           RFLU_WriteVersionString
                              
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString,stampString
  INTEGER, INTENT(IN) :: verbLevel

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: casename,nRegions,RCSIdentString,stamp
  INTEGER, PARAMETER :: CONV_RFLU2RFLU_ASC2BIN_GRID     = 10, & 
                        CONV_RFLU2RFLU_ASC2BIN_GRIDFLOW = 11, & 
                        CONV_RFLU2RFLU_BIN2ASC_GRID     = 20, & 
                        CONV_RFLU2RFLU_BIN2ASC_GRIDFLOW = 21, & 
                        CONV_RFLU2TETM_GRID             = 40, & 
			CONV_RFLU2STL_GRID              = 50                        
  INTEGER :: convOption,errorFlag,iReg
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: rfluconv.F90,v $ $Revision: 1.6 $'

! ******************************************************************************
! Start, initialize global data
! ******************************************************************************

  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag 

  casename = caseString(1:LEN(caseString))
  stamp    = stampString(1:LEN(stampString))
  
  CALL RFLU_InitGlobal(casename,verbLevel,CRAZY_VALUE_INT,global)

  CALL RegisterFunction(global,'main', &
                        'rfluconv.F90')

! ******************************************************************************
! Print header
! ******************************************************************************

  CALL RFLU_WriteVersionString(global)
  CALL RFLU_PrintHeader(global)

! ******************************************************************************
! Get conversion option
! ******************************************************************************

  WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Options:'
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
                              'Convert Rocflu files from ASCII to binary format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BIN_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BIN_GRIDFLOW, &
                              ': Grid and flow file'   
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
                              'Convert Rocflu files from binary to ASCII format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BIN2ASC_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BIN2ASC_GRIDFLOW, &
                              ': Grid and flow file' 
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &                                  
                              'Convert Rocflu grid file to surface grid:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2TETM_GRID, &
                              ': Tetmesh/YAMS surface grid (msh2 format)' 
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2STL_GRID, &
                              ': STL surface grid' 			                                                                 
  WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Enter selection:'

  READ(STDIN,*) convOption

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

  CALL RFLU_GetUserInput(levels(1)%regions,.TRUE.) 

  IF ( global%flowType == FLOW_STEADY ) THEN
    READ(stamp,*) global%currentIter
  ELSE
    READ(stamp,*) global%currentTime
  END IF ! global%flowType  

! ******************************************************************************
! Read dimensions and allocate memory
! ****************************************************************************** 

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    
    CALL RFLU_ReadDimensions(pRegion)      
    CALL RFLU_CreateGrid(pRegion)       
 
    IF ( pRegion%grid%nPatches > 0 ) THEN      
      CALL RFLU_ReadBcInputFileWrapper(pRegion)    
    END IF ! pRegion%grid%nPatches     

    CALL RFLU_AllocMemSolWrapper(pRegion)
    CALL RFLU_SetVarInfoWrapper(pRegion)    
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_ReadDimensions(pRegion)         
      CALL RFLU_CreateGrid(pRegion)

      IF ( pRegion%grid%nPatches > 0 ) THEN      
        CALL RFLU_ReadBCInputFileWrapper(pRegion)    
      END IF ! pRegion%grid%nPatches
      
      CALL RFLU_AllocMemSolWrapper(pRegion)
      CALL RFLU_SetVarInfoWrapper(pRegion)             
    END DO ! iReg  
  END IF ! global%nRegions

! ******************************************************************************
! Convert
! ******************************************************************************

  SELECT CASE ( convOption )

! ==============================================================================
!   Rocflu files from ASCII to binary - NOTE that the gridFormat and solutFormat 
!   variables are overridden.
! ==============================================================================  
  
! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------
  
    CASE (CONV_RFLU2RFLU_ASC2BIN_GRID)        
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          
      
        pRegion%global%gridFormat = FORMAT_ASCII           
        CALL RFLU_ReadGridWrapper(pRegion)  
      
        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_BINARY
        CALL RFLU_WriteGridWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 
         
          pRegion%global%gridFormat = FORMAT_ASCII              
          CALL RFLU_ReadGridWrapper(pRegion)  
      
          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_BINARY
          CALL RFLU_WriteGridWrapper(pRegion)                
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_ASC2BIN_GRIDFLOW)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%gridFormat  = FORMAT_ASCII
        pRegion%global%solutFormat = FORMAT_ASCII   
        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY      
        pRegion%global%solutFormat = FORMAT_BINARY      

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%gridFormat  = FORMAT_ASCII
          pRegion%global%solutFormat = FORMAT_ASCII   
          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY      
          pRegion%global%solutFormat = FORMAT_BINARY      

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 
        END DO ! iReg       
      END IF ! global%nRegions    
      
! ==============================================================================
!   Rocflu files from binary to ASCII - NOTE that the gridFormat and solutFormat 
!   variables are overridden.
! ==============================================================================        

! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------
      
    CASE (CONV_RFLU2RFLU_BIN2ASC_GRID)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          
      
        pRegion%global%gridFormat = FORMAT_BINARY          
        CALL RFLU_ReadGridWrapper(pRegion)  
      
        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_ASCII
        CALL RFLU_WriteGridWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 
         
          pRegion%global%gridFormat = FORMAT_BINARY             
          CALL RFLU_ReadGridWrapper(pRegion)  
      
          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_ASCII
          CALL RFLU_WriteGridWrapper(pRegion)                
        END DO ! iReg       
      END IF ! global%nRegions 
      
! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BIN2ASC_GRIDFLOW) 
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%gridFormat  = FORMAT_BINARY
        pRegion%global%solutFormat = FORMAT_BINARY  
        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_ASCII      
        pRegion%global%solutFormat = FORMAT_ASCII      

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%gridFormat  = FORMAT_BINARY
          pRegion%global%solutFormat = FORMAT_BINARY  
          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_ASCII     
          pRegion%global%solutFormat = FORMAT_ASCII     

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 
        END DO ! iReg       
      END IF ! global%nRegions    

! ==============================================================================
!   Convert Rocflu grid to surface grid
! ==============================================================================

! ------------------------------------------------------------------------------
!   Tetmesh/YAMS .msh2 format
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2TETM_GRID) 
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

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

        CALL RFLU_ConvROCFLU2TETMESH(pRegion)       
        CALL RFLU_WriteSurfGridTETMESH(pRegion)                

        CALL RFLU_DestroyFaceList(pRegion)                    
        CALL RFLU_DestroyBVertexLists(pRegion)	
        CALL RFLU_DestroyCellMapping(pRegion)                            
      ELSE ! multiple regions, must not happen
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)       
      END IF ! global%nRegions  

! ------------------------------------------------------------------------------
!   STL format
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2STL_GRID) 
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

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

        CALL RFLU_CreateGeometry(pRegion)
	CALL RFLU_BuildGeometry(pRegion)			   

        CALL RFLU_STL_WriteSurfGridASCII(pRegion)                

        CALL RFLU_DestroyGeometry(pRegion)
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)	                    
        CALL RFLU_DestroyCellMapping(pRegion)                            
      ELSE ! multiple regions, must not happen
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)       
      END IF ! global%nRegions  

! ==============================================================================
!   Default
! ==============================================================================        

    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! convOption

! ******************************************************************************
! Deallocate memory for solution and grid
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    
    CALL RFLU_DeallocMemSolWrapper(pRegion) 
    CALL RFLU_DestroyGrid(pRegion)
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_DeallocMemSolWrapper(pRegion) 
      CALL RFLU_DestroyGrid(pRegion)
    END DO ! iReg  
  END IF ! global%nRegions

! *****************************************************************************
! Print info about warnings
! *****************************************************************************

  CALL RFLU_PrintWarnInfo(global)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE rfluconv

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rfluconv.F90,v $
! Revision 1.6  2008/12/06 08:44:55  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/02/06 23:55:54  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.3  2005/07/07 03:51:09  haselbac
! Added STL conversion, some clean-up
!
! Revision 1.2  2005/05/03 03:09:09  haselbac
! Converted to C++ reading of command-line, fixed bug in formats
!
! Revision 1.1  2005/04/18 14:57:56  haselbac
! Initial revision
!
! ******************************************************************************







