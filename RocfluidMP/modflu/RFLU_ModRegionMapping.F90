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
! Purpose: Collection of routines relating to mapping of regions to processes. 
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
! $Id: RFLU_ModRegionMapping.F90,v 1.7 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRegionMapping

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_level
  USE ModMPI 

  USE ModBuildFileNames, ONLY: BuildFileNamePlain
    
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_ApplyRegionMapping, &             
            RFLU_BuildRegionMappingSimple, & 
            RFLU_CheckRegionMapping, & 
            RFLU_CloseRegionMappingFile, &
            RFLU_CreateRegionMapping, & 
            RFLU_DestroyRegionMapping, &
            RFLU_GetProcLocRegIds, & 
            RFLU_ImposeRegionMappingSerial, & 
            RFLU_OpenRegionMappingFile, &
            RFLU_ReadRegionMappingFile, & 
            RFLU_SetRegionMappingSerial, & 
            RFLU_WriteRegionMappingFile

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModRegionMapping.F90,v $ $Revision: 1.7 $'
    
  INTEGER, PARAMETER, PUBLIC :: MAPTYPE_REG      = 1, & 
                                MAPTYPE_PROC2REG = 2        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS









! ******************************************************************************
!
! Purpose: Apply region mapping.
!
! Description: None.
!
! Input: 
!   global      Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ApplyRegionMapping(global,levels)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_level), POINTER :: levels(:)
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: errorFlag,iLev,iReg
  
! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ApplyRegionMapping',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Applying region mapping...'
  END IF ! global%verbLevel

! ******************************************************************************
! Copy global region number into region type
! ******************************************************************************
  
  iLev = 1
  
  DO iReg = 0,global%nRegionsLocal
    levels(iLev)%regions(iReg)%iRegionGlobal = global%regMap(iReg)
  END DO ! iReg  
  
  DO iLev = 2,global%nLevels    
    DO iReg = 1,global%nRegionsLocal
      levels(iLev)%regions(iReg)%iRegionGlobal = global%regMap(iReg)
    END DO ! iReg
  END DO ! iLev
    
! ******************************************************************************
! End
! ******************************************************************************
  
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Applying region mapping done.'    
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)  
    
END SUBROUTINE RFLU_ApplyRegionMapping







! ******************************************************************************
!
! Purpose: Build simple region mapping.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: 
!   1. This mapping simply maps an equal number of regions to each process 
!      without looking at the dimensions of each region.
!   2. If the number of regions is not divisible by the number of processes, 
!      the last process gets the remainder of regions. 
!
! ******************************************************************************

SUBROUTINE RFLU_BuildRegionMappingSimple(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================    
! Arguments
! ==============================================================================    
   
  TYPE(t_global), POINTER :: global 
   
! ==============================================================================    
! Locals
! ==============================================================================    
  
  INTEGER :: i,iBeg,iEnd,iOffset,iProc,nRegsPerProc

! ******************************************************************************
! Start
! ******************************************************************************  
    
  CALL RegisterFunction(global,'RFLU_BuildRegionMappingSimple',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_MED) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building region mapping...'
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Method: Simple'
  END IF ! global%verbLevel  
    
! ******************************************************************************
! Build simple region mapping
! ******************************************************************************  

  nRegsPerProc = global%nRegions/global%nProcs ! NOTE integer division

  IF ( global%nRegions > 1 ) THEN 
    iOffset = 0
  ELSE 
    iOffset = 1
  END IF ! global%nRegions

  DO iProc = 1,global%nProcs  
    IF ( iProc > 1 ) THEN 
      iBeg = global%proc2RegMapInfo(2,iProc-1) + 1
    ELSE 
      iBeg = 1 
    END IF ! iProc 
    
    IF ( iProc /= global%nProcs ) THEN 
      iEnd = iBeg + nRegsPerProc - 1      
    ELSE ! Add remainder to last process
      iEnd = global%nRegions
    END IF ! iProc

    global%proc2RegMapInfo(1,iProc) = iBeg
    global%proc2RegMapInfo(2,iProc) = iEnd
  
    DO i = iBeg,iEnd
      global%proc2RegMap(i) = i - iOffset
    END DO ! i
  END DO ! iProc
    
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building region mapping done.'
  END IF ! global%verbLevel 

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_BuildRegionMappingSimple








! ******************************************************************************
!
! Purpose: Check region mapping.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_CheckRegionMapping(global)

  USE ModSortSearch, ONLY: QuickSortInteger 

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================    
! Arguments
! ==============================================================================    
   
  TYPE(t_global), POINTER :: global 
   
! ==============================================================================    
! Locals
! ==============================================================================    
  
  INTEGER :: checkSum,errorFlag,i
  INTEGER, DIMENSION(:), ALLOCATABLE :: proc2RegMapTemp

! ******************************************************************************
! Start
! ******************************************************************************  
    
  CALL RegisterFunction(global,'RFLU_CheckRegionMapping',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking region mapping...'
  END IF ! global%verbLevel  
    
! ******************************************************************************
! Check region mapping
! ******************************************************************************  

  IF ( global%nRegions > 1 ) THEN 
    ALLOCATE(proc2RegMapTemp(global%nRegions),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'proc2RegMapTemp')
    END IF ! global%error

    DO i = 1,global%nRegions
      proc2RegMapTemp(i) = global%proc2RegMap(i)
    END DO ! i

    CALL QuickSortInteger(proc2RegMapTemp,global%nRegions)

    checkSum = 0
    
    DO i = 1,global%nRegions
      checkSum = checkSum + proc2RegMapTemp(i)
    END DO ! i
  
    IF ( checkSum /= global%nRegions*(global%nRegions+1)/2 ) THEN 
      CALL ErrorStop(global,ERR_PROC2REG_MAPPING,__LINE__)  
    END IF ! checkSum
    
    DEALLOCATE(proc2RegMapTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'proc2RegMapTemp')
    END IF ! global%error    
  END IF ! global%nRegions
    
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking region mapping done.'
  END IF ! global%verbLevel 

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_CheckRegionMapping








! ******************************************************************************
!
! Purpose: Close region mapping file.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_CloseRegionMappingFile(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================    
! Arguments
! ==============================================================================    
   
  TYPE(t_global), POINTER :: global 
   
! ==============================================================================    
! Locals
! ==============================================================================    
  
  CHARACTER(CHRLEN) :: iFileName
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************  
    
  CALL RegisterFunction(global,'RFLU_CloseRegionMappingFile',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing region mapping file...'
  END IF ! global%verbLevel  
    
! ******************************************************************************
! Close file
! ******************************************************************************  

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.map',iFileName)  
    
  CLOSE(IF_REGMAP,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName) 
  END IF ! global%error    
    
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing region mapping file done.'
  END IF ! global%verbLevel 

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_CloseRegionMappingFile








! ******************************************************************************
!
! Purpose: Create region mapping.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_CreateRegionMapping(global,mapType)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================    
! Arguments
! ==============================================================================    
   
  INTEGER, INTENT(IN) :: mapType 
  TYPE(t_global), POINTER :: global 
   
! ==============================================================================    
! Locals
! ==============================================================================    
  
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************  
    
  CALL RegisterFunction(global,'RFLU_CreateRegionMapping',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating region mapping...'
  END IF ! global%verbLevel  
    
! ******************************************************************************
! Create region mapping
! ******************************************************************************  

  IF ( mapType == MAPTYPE_REG ) THEN 
    ALLOCATE(global%regMap(0:global%nRegionsLocal),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'global%regMap')
    END IF ! global%error
  ELSE IF ( mapType == MAPTYPE_PROC2REG ) THEN 
    ALLOCATE(global%proc2RegMap(global%nRegions),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'global%proc2RegMap')
    END IF ! global%error

    ALLOCATE(global%proc2RegMapInfo(2,global%nProcs),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'global%proc2RegMapInfo')
    END IF ! global%error    
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! mapType
    
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating region mapping done.'
  END IF ! global%verbLevel 

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_CreateRegionMapping








! ******************************************************************************
!
! Purpose: Destroy region mapping.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_DestroyRegionMapping(global,mapType)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================    
! Arguments
! ==============================================================================    
   
  INTEGER, INTENT(IN) :: mapType 
  TYPE(t_global), POINTER :: global 
   
! ==============================================================================    
! Locals
! ==============================================================================    
  
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************  
    
  CALL RegisterFunction(global,'RFLU_DestroyRegionMapping',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying region mapping...'
  END IF ! global%verbLevel  
    
! ******************************************************************************
! Create region mapping
! ******************************************************************************  

  IF ( mapType == MAPTYPE_REG ) THEN 
    DEALLOCATE(global%regMap,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'global%regMap')
    END IF ! global%error
  ELSE IF ( mapType == MAPTYPE_PROC2REG ) THEN 
    DEALLOCATE(global%proc2RegMap,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'global%proc2RegMap')
    END IF ! global%error

    DEALLOCATE(global%proc2RegMapInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'global%proc2RegMapInfo')
    END IF ! global%error    
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! mapType
    
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying region mapping done.'
  END IF ! global%verbLevel 

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_DestroyRegionMapping







! ******************************************************************************
!
! Purpose: Get process and local region ids given a set of global region ids.
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
!   regIds              Set of region ids
!   nRegIds             Number of region ids
!
! Output: 
!   procIds             Set of process ids
!   locRegIds           Set of local region ids
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_GetProcLocRegIds(global,regIds,nRegIds,procIds,locRegIds)
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  INTEGER, INTENT(IN) :: nRegIds
  INTEGER, INTENT(IN) :: regIds(nRegIds)
  INTEGER, INTENT(OUT) :: locRegIds(nRegIds),procIds(nRegIds)
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: sectionString,targetString
  INTEGER :: errorFlag,iFile,iProc,iReg,iRegId,iRegId2,nProcs,nRegionsGlobal, & 
             nRegionsLocal

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_GetProcLocRegIds',&
  'RFLU_ModRegionMapping.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  iFile = IF_REGMAP

! ******************************************************************************
! Initialize output arrays
! ******************************************************************************

  DO iRegId = 1,nRegIds
    procIds(iRegId)   = CRAZY_VALUE_INT
    locRegIds(iRegId) = CRAZY_VALUE_INT
  END DO ! iRegIds

! ******************************************************************************
! Loop over region ids
! ******************************************************************************

  DO iRegId = 1,nRegIds
    REWIND(iFile)
  
! ==============================================================================
!   Read header
! ==============================================================================
  
    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU region mapping file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
    END IF ! TRIM

! ==============================================================================
!   Read number of regions and processors 
! ==============================================================================

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Number of regions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
    END IF ! TRIM

    READ(iFile,*) nRegionsGlobal

    IF ( nRegionsGlobal /= global%nRegions ) THEN 
      CALL ErrorStop(global,ERR_NREGIONS_MISMATCH,__LINE__)
    END IF ! nRegionsGlobal

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Number of processes' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
    END IF ! TRIM    

    READ(iFile,*) nProcs

    IF ( nProcs /= global%nProcs ) THEN 
      CALL ErrorStop(global,ERR_NPROCS_MISMATCH,__LINE__)
    END IF ! nProcs
    
! ==============================================================================
!   Loop over processes and look for region id
! ==============================================================================

    procLoop: DO iProc = 0,global%nProcs-1
      WRITE(targetString,'(A,1X,I6.6)') '# Process',iProc+1       

      READ(iFile,'(A)') sectionString        

      IF ( TRIM(sectionString) /= TRIM(targetString) ) THEN 
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
      END IF ! TRIM        

      READ(iFile,*) nRegionsLocal

      DO iReg = 1,nRegionsLocal
        READ(iFile,*) iRegId2

        IF ( iRegId2 == regIds(iRegId) ) THEN 
          procIds(iRegId)   = iProc+1
          locRegIds(iRegId) = iReg
          
          EXIT procLoop          
        END IF ! iRegId
      END DO ! iReg
    END DO procLoop
  END DO ! iRegId   

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetProcLocRegIds








! ******************************************************************************
!
! Purpose: Impose serial region mapping.
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
!
! Output: None.
!
! Notes: 
!   1. This routine must only be called if RFLU_ReadRegionMappingFile
!      was called with readMode equal to MAPFILE_READMODE_PEEK. There is no 
!      check for this here, so it is up to the user to make sure this is so.
!
! ******************************************************************************

SUBROUTINE RFLU_ImposeRegionMappingSerial(global)
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,iReg

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ImposeRegionMappingSerial',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Imposing serial region mapping...' 
  END IF ! global%myProcid

! ******************************************************************************
! Set region mapping
! ******************************************************************************

  DO iReg = 0,global%nRegionsLocal
    global%regMap(iReg) = iReg        
  END DO ! iReg

  IF ( global%nRegionsLocal == 1 ) THEN ! for single-processor jobs
    global%regMap(1) = 0 
  END IF ! global%nRegionsLocal

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Imposing serial region mapping done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ImposeRegionMappingSerial








! ******************************************************************************
!
! Purpose: Open region mapping file.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_OpenRegionMappingFile(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================    
! Arguments
! ==============================================================================    
   
  TYPE(t_global), POINTER :: global 
   
! ==============================================================================    
! Locals
! ==============================================================================    
  
  CHARACTER(CHRLEN) :: iFileName 
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************  
    
  CALL RegisterFunction(global,'RFLU_OpenRegionMappingFile',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening region mapping file...'
  END IF ! global%verbLevel  
    
! ******************************************************************************
! Open file
! ******************************************************************************  
    
  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.map',iFileName)  

  OPEN(IF_REGMAP,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName) 
  END IF ! global%error    
    
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening region mapping file done.'
  END IF ! global%verbLevel 

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_OpenRegionMappingFile









! ******************************************************************************
!
! Purpose: Read region mapping file for parallel processing
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
!   readMode            Reading mode of file, can take only two values:
!                         MAPFILE_READMODE_ALL:  Read whole file
!                         MAPFILE_READMODE_PEEK: Read only number of regions
!   myProcId            Processor id
!
! Output: None.
!
! Notes: 
!   1. The complete file is only read if the readMode argument is set to 
!      MAPFILE_READMODE_ALL. If the readMode argument is set to 
!      MAPFILE_READMODE_PEEK, then only the number of regions is read, and 
!      the rest of the file is ignored. 
!
! ******************************************************************************

SUBROUTINE RFLU_ReadRegionMappingFile(global,readMode,myProcId)
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  INTEGER, INTENT(IN) :: myProcId,readMode 
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: fileExists,foundEntry
  CHARACTER(CHRLEN) :: dummyString,iFileName,sectionString,targetString
  INTEGER :: errorFlag,iFile,iProc,iReg,nRegions

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ReadRegionMappingFile',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC ) THEN 
    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading region mapping file...'
      
      IF ( global%verbLevel >= VERBOSE_MED ) THEN 
        IF ( readMode == MAPFILE_READMODE_ALL ) THEN 
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Read mode: All'
        ELSE IF ( readMode == MAPFILE_READMODE_PEEK ) THEN 
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Read mode: Peek'    
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! readMode
      END IF ! global%verbLevel  
    END IF ! global%verbLevel
  END IF ! global%myProcid

! ******************************************************************************
! Test for existence of mapping file
! ******************************************************************************

  iFile = IF_REGMAP
  
  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.map',iFileName)

  INQUIRE(FILE=iFileName,EXIST=fileExists)

! ==============================================================================
! File exists
! ==============================================================================
  
  IF ( fileExists .EQV. .TRUE. ) THEN     
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mapping file found.' 
    END IF ! global%verbLevel  
  
! ------------------------------------------------------------------------------
!   Open file and read header
! ------------------------------------------------------------------------------  
  
    CALL RFLU_OpenRegionMappingFile(global)

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU region mapping file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
    END IF ! TRIM

! ------------------------------------------------------------------------------
!   Read number of regions and processors 
! ------------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Number of regions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
    END IF ! TRIM

    READ(iFile,*) global%nRegions

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Number of processes' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
    END IF ! TRIM    

    READ(iFile,*) global%nProcs

! ------------------------------------------------------------------------------
!   Read rest of file
! ------------------------------------------------------------------------------

    DO iProc = 0,global%nProcs-1
      foundEntry = .FALSE. 
    
      WRITE(targetString,'(A,1X,I6.6)') '# Process',iProc+1       
      READ(iFile,'(A)') sectionString 
             
      IF ( TRIM(sectionString) /= TRIM(targetString) ) THEN 
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
      END IF ! TRIM        

      READ(iFile,*) nRegions

      IF ( iProc == myProcId ) THEN ! Found entry for my processor
        foundEntry = .TRUE.
      
        global%nRegionsLocal = nRegions                        
      END IF ! iProc

      IF ( (foundEntry .EQV. .TRUE.) .AND. & 
           (readMode == MAPFILE_READMODE_ALL) ) THEN
        global%regMap(0) = 0         

        DO iReg = 1,global%nRegionsLocal
          READ(iFile,*) global%regMap(iReg)
        END DO ! iReg
      ELSE       
        DO iReg = 1,nRegions
          READ(iFile,'(A)') dummyString ! Read irrelevant lines
        END DO ! iReg                    
      END IF ! readMode  
    END DO ! iProc

    READ(iFile,'(A)') sectionString        
    IF ( TRIM(sectionString) /= '# End' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,TRIM(sectionString)) 
    END IF ! TRIM 

    CALL RFLU_CloseRegionMappingFile(global)       
    
! ==============================================================================
! File does not exist, initialize for single-processor run
! ==============================================================================    
    
  ELSE IF ( (fileExists .EQV. .FALSE.) .AND. (global%nProcAlloc == 1) ) THEN
    global%warnCounter = global%warnCounter + 1  
   
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
                                    'Mapping file not found.'
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Initializing for', & 
                                    'single-process run.'
    END IF ! global%verbLevel  
 
    global%nProcs        = 1
    global%nRegions      = 1
    global%nRegionsLocal = 1
    
! ==============================================================================
! File should exist but does not 
! ==============================================================================    
        
  ELSE 
    CALL ErrorStop(global,ERR_FILE_EXIST,__LINE__,TRIM(iFileName))
  END IF ! fileExists

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading region mapping file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadRegionMappingFile









! ******************************************************************************
!
! Purpose: Set serial region mapping.
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
!
! Output: None.
!
! Notes: 
!   1. This routine must only be called if RFLU_ReadRegionMappingFile
!      was called with readMode equal to MAPFILE_READMODE_PEEK. There is no 
!      check for this here, so it is up to the user to make sure this is so.
!
! ******************************************************************************

SUBROUTINE RFLU_SetRegionMappingSerial(global)
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_SetRegionMappingSerial',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting serial region mapping...' 
  END IF ! global%myProcid

! ******************************************************************************
! Set serial region mapping
! ******************************************************************************

  global%nRegionsLocal = global%nRegions

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting serial region mapping done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetRegionMappingSerial









! ******************************************************************************
!
! Purpose: Write region mapping file.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_WriteRegionMappingFile(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================    
! Arguments
! ==============================================================================    
   
  TYPE(t_global), POINTER :: global 
   
! ==============================================================================    
! Locals
! ==============================================================================    
  
  INTEGER :: i,iBeg,iEnd,iProc,nRegions

! ******************************************************************************
! Start
! ******************************************************************************  
    
  CALL RegisterFunction(global,'RFLU_WriteRegionMappingFile',&
  'RFLU_ModRegionMapping.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing region mapping file...'
  END IF ! global%verbLevel  
    
! ******************************************************************************
! Write mapping file
! ******************************************************************************  
    
  WRITE(IF_REGMAP,'(A)') '# ROCFLU region mapping file'
  
  WRITE(IF_REGMAP,'(A)') '# Number of regions'
  WRITE(IF_REGMAP,'(I7)') global%nRegions
    
  WRITE(IF_REGMAP,'(A)') '# Number of processes'
  WRITE(IF_REGMAP,'(I7)') global%nProcs
    
  DO iProc = 1,global%nProcs       
    WRITE(IF_REGMAP,'(A,1X,I6.6)') '# Process',iProc

    iBeg = global%proc2RegMapInfo(1,iProc)
    iEnd = global%proc2RegMapInfo(2,iProc)

    WRITE(IF_REGMAP,'(I7)') iEnd - iBeg + 1
    
    DO i = iBeg,iEnd
      WRITE(IF_REGMAP,'(I7)') global%proc2RegMap(i)
    END DO ! i 
  END DO ! iProc  
    
  WRITE(IF_REGMAP,'(A)') '# End'    
    
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing region mapping file done.'
  END IF ! global%verbLevel 

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_WriteRegionMappingFile







! ******************************************************************************
! End
! ******************************************************************************  

END MODULE RFLU_ModRegionMapping

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRegionMapping.F90,v $
! Revision 1.7  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/03/31 23:53:13  haselbac
! Changed global region index for region 0 in parallel jobs
!
! Revision 1.4  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.3  2005/04/15 15:07:01  haselbac
! Modified routine to get local proc ids
!
! Revision 1.2  2005/01/14 21:25:53  haselbac
! Added routine to get pid given regid, changed reading of map file
!
! Revision 1.1  2004/10/19 19:26:57  haselbac
! Initial revision
!
! ******************************************************************************


















