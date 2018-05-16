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
! Purpose: Driver to extract special cells for postprocessor.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: 
!   1. This utility is needed because rflupost_charm cannot accept input
!      from the command-line. Hence extract the information here and write to 
!      text file.
!
!*******************************************************************************
!
! $Id: main.F90,v 1.5 2008/12/06 08:44:54 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!*******************************************************************************

PROGRAM rflucells

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModParameters
  
  USE RFLU_ModCellMapping
  USE RFLU_ModFaceList
  USE RFLU_ModFEM  
  
  USE ModInterfaces, ONLY: RFLU_AllocateMemoryWrapper, & 
                           RFLU_AssignRegionMapping,RFLU_BuildDataStruct, &
                           RFLU_CloseSpecialCellsFile, &
                           RFLU_ComputeDimensions, &
                           RFLU_DeallocateMemoryWrapper,RFLU_GetSpecialCells, &                       
                           RFLU_GetUserInput,RFLU_ImposeSerialMapping, & 
                           RFLU_InitGlobal,RFLU_OpenSpecialCellsFile, &
                           RFLU_PrintGridInfo,RFLU_PrintHeader, &
                           RFLU_ReadDimensions,RFLU_ReadGridWrapper, &
                           RFLU_ReadInputFile,RFLU_ReadMappingFile, &
                           RFLU_WriteSpecialCells,RFLU_WriteVersionString                                                     
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: specialCellFileExists
  CHARACTER(CHRLEN) :: casename,nRegions,RCSIdentString,stampStr,verbosity
  INTEGER :: errorFlag,iReg,verbLevel
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: main.F90,v $ $Revision: 1.5 $'

! ******************************************************************************
! Get command-line options
! ******************************************************************************

  CALL GETARG(1,casename)
  CALL GETARG(2,stampStr)  
  CALL GETARG(3,verbosity)

  IF ( LEN_TRIM(casename)  == 0 .OR. & 
       LEN_TRIM(stampStr)  == 0 .OR. & 
       LEN_TRIM(verbosity) == 0         ) THEN
    WRITE(STDOUT,'(A     )') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Usage: rflucells <casename> <stamp> <verb>'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'       stamp = iteration (steady flow)'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'               time (unsteady flow)'    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'       verb  = 0 - no verbose output'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'               1 - some verbose output'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'               2 - a lot of verbose output'        
    WRITE(STDOUT,'(A     )') SOLVER_NAME
    STOP
  END IF ! LEN_TRIM

  READ(verbosity,*) verbLevel

! ******************************************************************************
! Initialize global data
! ******************************************************************************
  
  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag   
  
  CALL RFLU_InitGlobal(casename,verbLevel,global)

  CALL RegisterFunction(global,'driver', &
                        'main.F90')

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

  CALL RFLU_ReadMappingFile(global,MAPFILE_READMODE_PEEK)
  CALL RFLU_ImposeSerialMapping(global)

! ******************************************************************************
! Prepare data structure
! ******************************************************************************
  
  CALL RFLU_BuildDataStruct(global,levels)
  CALL RFLU_AssignRegionMapping(global,levels)
  
! ******************************************************************************  
! Read input file
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions)  

  IF ( global%flowType == FLOW_STEADY ) THEN
    READ(stampStr,*) global%currentIter
  ELSE
    READ(stampStr,*) global%currentTime
  END IF ! global%flowType  

! ******************************************************************************  
! Read dimensions file, allocate memory, and read bc file
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    
    CALL RFLU_ReadDimensions(pRegion)         
    CALL RFLU_AllocateMemoryWrapper(pRegion,ALLOC_MODE_ALL)      
    
    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_CreateFaceList(pRegion)
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)
      
      CALL RFLU_ReadDimensions(pRegion) 
      CALL RFLU_AllocateMemoryWrapper(pRegion,ALLOC_MODE_ALL)
           
      CALL RFLU_CreateCellMapping(pRegion)
      CALL RFLU_CreateFaceList(pRegion)                                       
    END DO ! iReg    
  END IF ! global%nRegions
  
! ******************************************************************************  
! Read grid
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    CALL RFLU_ReadGridWrapper(pRegion)  
    
    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      CALL RFLU_PrintGridInfo(pRegion)
    END IF ! global%verbLevel    
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)
      CALL RFLU_ReadGridWrapper(pRegion)  

      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        CALL RFLU_PrintGridInfo(pRegion)
      END IF ! global%verbLevel     
    END DO ! iReg            
  END IF ! global%nRegions

! ******************************************************************************
! Build data structure 
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    CALL RFLU_BuildCellMapping(pRegion)     
    CALL RFLU_BuildFaceList(pRegion)  
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)
      CALL RFLU_BuildCellMapping(pRegion)     
      CALL RFLU_BuildFaceList(pRegion)      
    END DO ! iReg            
  END IF ! global%nRegions 

! ******************************************************************************
! Get special cells and write to file
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    CALL RFLU_GetSpecialCells(pRegion)    
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)
      CALL RFLU_GetSpecialCells(pRegion)            
    END DO ! iReg            
  END IF ! global%nRegions   

  CALL RFLU_OpenSpecialCellsFile(global,FILE_STATUS_UNKNOWN, &
                                 specialCellFileExists)
  
  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    CALL RFLU_WriteSpecialCells(pRegion)     
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)
      CALL RFLU_WriteSpecialCells(pRegion)            
    END DO ! iReg            
  END IF ! global%nRegions   
  
  CALL RFLU_CloseSpecialCellsFile(global)

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN 
    pRegion => levels(1)%regions(0)  
    CALL RFLU_DeallocateMemoryWrapper(pRegion)
  ELSE 
    DO iReg = 1,global%nRegions
      pRegion => levels(1)%regions(iReg)          
      CALL RFLU_DeallocateMemoryWrapper(pRegion)   
    END DO ! iReg
  END IF ! global%nRegions
  
! ******************************************************************************
! End
! ******************************************************************************

  WRITE(STDOUT,'(A)')      SOLVER_NAME
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finished.'
  WRITE(STDOUT,'(A)')      SOLVER_NAME

  CALL DeregisterFunction(global)

END PROGRAM rflucells

!*******************************************************************************
!
! RCS Revision history:
!
! $Log: main.F90,v $
! Revision 1.5  2008/12/06 08:44:54  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2003/04/07 14:28:34  haselbac
! Fixed bug: Deleted superfluous call
!
! Revision 1.2  2003/04/01 19:39:20  haselbac
! Changed call to open special cells file
!
! Revision 1.1  2003/04/01 17:02:46  haselbac
! Initial revision
!
!*******************************************************************************







