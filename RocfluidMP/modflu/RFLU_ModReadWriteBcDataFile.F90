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
! Purpose: Collection of routines to read boundary condition data file.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadWriteBcDataFile.F90,v 1.6 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadWriteBcDataFile

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_DecideReadWriteBcDataFile, & 
            RFLU_ReadBcDataFile, & 
            RFLU_WriteBcDataFile
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModReadWriteBcDataFile.F90,v $ $Revision: 1.6 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Decide whether need to read and/or write boundary-condition data 
!   file.
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

  LOGICAL FUNCTION RFLU_DecideReadWriteBcDataFile(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: RCSIdentString
    INTEGER :: iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DecideReadWriteBcDataFile',&
  'RFLU_ModReadWriteBcDataFile.F90')
    
! ******************************************************************************
!   Set pointers and initialize
! ******************************************************************************    
    
    pGrid => pRegion%grid  
    
    RFLU_DecideReadWriteBcDataFile = .FALSE.  
    
! ******************************************************************************
!   Loop over patches
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
  
! ==============================================================================
!     If not coupled and have distribution, write data
! ==============================================================================    
  
      IF ( (pPatch%bcCoupled == BC_NOT_COUPLED) .AND. & 
           (pPatch%mixt%distrib == BCDAT_DISTRIB) ) THEN
! TEMPORARY
!        RFLU_DecideReadWriteBcDataFile = .TRUE.
        RFLU_DecideReadWriteBcDataFile = .FALSE.
! END TEMPORARY

        EXIT
      END IF ! pPatch%bcCoupled
    END DO ! iPatch 

! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)

  END FUNCTION RFLU_DecideReadWriteBcDataFile








! ******************************************************************************
!
! Purpose: Read boundary-condition data file.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. At present restricted to mixture. Best way to extend to MP not clear.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcDataFile(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,iFileName,sectionString
    INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch,iPatchGlobal, & 
               loopCounter,nBFaces,nData
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcDataFile',&
  'RFLU_ModReadWriteBcDataFile.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading boundary-condition data...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Open file
! ******************************************************************************

    iFile = IF_DISTR

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.bcd', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   Read header
! ******************************************************************************

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU boundary-condition data file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read data
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Patch data
! ==============================================================================
        
        CASE ( '# Patch data' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch data...'
          END IF ! global%verbLevel       

          READ(iFile,*) iPatch,iPatchGlobal,nBFaces,nData  
            
! ------------------------------------------------------------------------------
!         Check that input correct
! ------------------------------------------------------------------------------            
                    
          IF ( iPatch > pGrid%nPatches ) THEN
            WRITE(errorString,'(A,1X,I3)') 'Patch index invalid:',iPatch 
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! iPatch                                                  

          pPatch => pRegion%patches(iPatch)

          IF ( nBFaces /= pPatch%nBFaces ) THEN 
            WRITE(errorString,'(A,1X,I6)') 'Number of faces invalid:',nBFaces
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! nBFaces          

          IF ( iPatchGlobal /= pPatch%iPatchGlobal ) THEN 
            WRITE(errorString,'(A,1X,I3)') 'Global patch index invalid:', & 
                                           iPatchGlobal
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! iPatchGlobal           

          IF ( nData /= pPatch%mixt%nData ) THEN 
            WRITE(errorString,'(A,1X,I3)') & 
              'Number of pieces of data invalid:',nData
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! iPatchGlobal    
    
! ------------------------------------------------------------------------------
!         Read data
! ------------------------------------------------------------------------------    
    
          DO ifl = 1,pPatch%nBFaces
            DO iData = 1,pPatch%mixt%nData
              READ(iFile,'(1X,I6,1X,I2,1X,E23.16)') & 
                dummyInteger,dummyInteger,pPatch%mixt%vals(iData,ifl)
            END DO ! iData  
          END DO ! ifl

! ==============================================================================
!       End marker
! ==============================================================================
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

          EXIT

! ==============================================================================
!       Invalid section string
! ==============================================================================

        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel           

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)      

      END SELECT ! TRIM

! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter

    END DO ! <empty> 

! ******************************************************************************
!   Close file
! ******************************************************************************
    
    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading boundary-condition data done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcDataFile








! ******************************************************************************
!
! Purpose: Write boundary-condition data file.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. At present restricted to mixture. Best way to extend to MP not clear.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteBcDataFile(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteBcDataFile',&
  'RFLU_ModReadWriteBcDataFile.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing boundary-condition data...'
    END IF ! global%verbLevel
  
! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
  
! ******************************************************************************
!   Open file
! ******************************************************************************

    iFile = IF_DISTR

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.bcd', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='UNKNOWN', & 
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   Write header
! ******************************************************************************
 
    sectionString = '# ROCFLU boundary-condition data file'   
    WRITE(iFile,'(A)') sectionString

! ******************************************************************************
!   Loop over patches
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
  
! ==============================================================================
!     If not coupled and have distribution, write data
! ==============================================================================    
  
      IF ( (pPatch%bcCoupled == BC_NOT_COUPLED) .AND. & 
           (pPatch%mixt%distrib == BCDAT_DISTRIB) ) THEN
        WRITE(iFile,'(A)') '# Patch data' 
        WRITE(iFile,'(2(1X,I3),1X,I6,1X,I2)') iPatch,pPatch%iPatchGlobal, & 
                                              pPatch%nBFaces, & 
                                              pPatch%mixt%nData

        DO ifl = 1,pPatch%nBFaces
          DO iData = 1,pPatch%mixt%nData
            WRITE(iFile,'(1X,I6,1X,I2,1X,E23.16)') ifl,iData, &
              pPatch%mixt%vals(iData,ifl)
          END DO ! iData  
        END DO ! ifl
      END IF ! pPatch%bcCoupled
    END DO ! iPatch

! ******************************************************************************
!   Write footer
! ******************************************************************************
 
    sectionString = '# End'   
    WRITE(iFile,'(A)') sectionString

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
        'Writing boundary-condition data done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteBcDataFile









! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModReadWriteBcDataFile


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadWriteBcDataFile.F90,v $
! Revision 1.6  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:39:14  mparmar
! Renamed patch variables
!
! Revision 1.3  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.2  2005/12/23 13:47:01  haselbac
! Temporarily disable bc data files to avoid hardcoding for runs with tbc
!
! Revision 1.1  2004/07/06 15:14:30  haselbac
! Initial revision
!
! ******************************************************************************









