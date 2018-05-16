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
! Purpose: Suite of routines to read and write species solution files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_ModReadWriteVars.F90,v 1.4 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE SPEC_RFLU_ModReadWriteVars

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE ModBuildFileNames, ONLY: BuildFileNameSteady, &
                               BuildFileNameUnsteady

  USE SPEC_ModParameters

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: SPEC_RFLU_ReadCvASCII, & 
            SPEC_RFLU_ReadCvBinary, & 
            SPEC_RFLU_ReadEEvASCII, & 
            SPEC_RFLU_ReadEEvBinary, & 
            SPEC_RFLU_WriteCvASCII, & 
            SPEC_RFLU_WriteCvBinary, & 
            SPEC_RFLU_WriteEEvASCII, & 
            SPEC_RFLU_WriteEEvBinary            
                        

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString = &
    '$RCSfile: SPEC_RFLU_ModReadWriteVars.F90,v $ $Revision: 1.4 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Read flow file for species in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!   2. For GENX runs, read file from time zero if restarting. This is for 
!      convenience and will have to be changed once grid adaptation is used.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_ReadCvASCII(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: fileExists
  CHARACTER(CHRLEN) :: errorString,iFileName,iFileNameOld,sectionString, &
                       timeString1,timeString2
  INTEGER :: errorFlag,i,iFile,iVars,j,loopCounter,nCellsTot,nCellsExpected, & 
             nVars,nVarsExpected,precActual,precExpected,rangeActual, & 
             rangeExpected
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadCvASCII',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII species cv file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.spec.cva', & 
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.spca', & 
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileNameOld)                               
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.spec.cva', & 
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.spca', & 
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileNameOld)                               
  END IF ! global%flowType

  iFile = IF_SOLUT

  INQUIRE(FILE=iFileName,EXIST=fileExists)

  IF ( fileExists .EQV. .TRUE. ) THEN 
    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error
  ELSE 
    OPEN(iFile,FILE=iFileNameOld,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileNameOld)
    END IF ! global%error    
  END IF ! fileExists

! ==============================================================================
! Set state vector state: Solution always stored in conservative form
! ==============================================================================

  pRegion%spec%cvState = CV_MIXT_STATE_CONS

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# ROCFLU species file' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
! -----------------------------------------------------------------------------
! Precision and range
! -----------------------------------------------------------------------------
    
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
  precExpected  = PRECISION(1.0_RFREAL)
  rangeExpected = RANGE(1.0_RFREAL)
  
  READ(iFile,'(2(I8))') precActual,rangeActual
  IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN 
    CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
  END IF ! precActual
  
! -----------------------------------------------------------------------------
! Physical time
! -----------------------------------------------------------------------------
  
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile,'(E23.16)') currentTime 
 
  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    IF ( global%currentTime < 0.0_RFREAL ) THEN
      global%currentTime = currentTime
    ELSE
      WRITE(timeString1,'(1PE11.5)') global%currentTime
      WRITE(timeString2,'(1PE11.5)') currentTime          
      IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
        CALL ErrorStop(global,ERR_TIME_SOLUTION,__LINE__,TRIM(iFileName))
      END IF ! global%currentTime 
    END IF ! global%currentTime
  END IF ! global%flowType  
  
! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  
  nVarsExpected  = pRegion%specInput%nSpecies
  nCellsExpected = pGrid%nCellsTot    
  
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM
    
  READ(iFile,'(2(I8))') nCellsTot,nVars
  IF ( nCellsTot /= nCellsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, & 
                                              'but expected:',nCellsExpected
    CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,errorString)
  END IF ! nCellsExpected     
  
  IF ( nVars /= nVarsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, & 
                                              'but expected:',nVarsExpected  
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVarsExpected   
  
! ==============================================================================
! Rest of file
! ==============================================================================

  iVars       = 0
  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1
  
    READ(iFile,'(A)') sectionString

    SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!     Species concentration
! ------------------------------------------------------------------------------

      CASE ( '# Density' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species concentration...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%spec%cv
      
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') (pCv(iVars,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     End marker
! ------------------------------------------------------------------------------ 
      
      CASE ( '# End' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel           
      
        EXIT
      
! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------ 
      
      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel           
      
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)        
             
    END SELECT ! TRIM
  
! ------------------------------------------------------------------------------
!   Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------  
  
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  
  END DO ! <empty>

! ==============================================================================
! Check and information about number of variables read
! ==============================================================================

  IF ( iVars /= nVars ) THEN 
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! iVar

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)   
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
    
! ******************************************************************************
! End
! ******************************************************************************

   IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII species cv file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
 
END SUBROUTINE SPEC_RFLU_ReadCvASCII






! ******************************************************************************
!
! Purpose: Read flow file for species in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!   2. For GENX runs, read file from time zero if restarting. This is for 
!      convenience and will have to be changed once grid adaptation is used.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_ReadCvBinary(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: fileExists
  CHARACTER(CHRLEN) :: errorString,iFileName,iFileNameOld,sectionString, &
                       timeString1,timeString2
  INTEGER :: errorFlag,i,iFile,iVars,j,loopCounter,nCellsTot,nCellsExpected, & 
             nVars,nVarsExpected,precActual,precExpected,rangeActual, & 
             rangeExpected
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadCvBinary',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary species cv file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.spec.cv', & 
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.spc', & 
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileNameOld)                               
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.spec.cv', & 
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)  
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.spc', & 
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileNameOld)                               
  END IF ! global%flowType

  iFile = IF_SOLUT

  INQUIRE(FILE=iFileName,EXIST=fileExists)

  IF ( fileExists .EQV. .TRUE. ) THEN 
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error
  ELSE 
    OPEN(iFile,FILE=iFileNameOld,FORM="UNFORMATTED",STATUS="OLD", &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileNameOld)
    END IF ! global%error    
  END IF ! fileExists

! ==============================================================================
! Set state vector state: Solution always stored in conservative form
! ==============================================================================

  pRegion%spec%cvState = CV_MIXT_STATE_CONS
  
! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# ROCFLU species file' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
! -----------------------------------------------------------------------------
! Precision and range
! -----------------------------------------------------------------------------
    
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
  precExpected  = PRECISION(1.0_RFREAL)
  rangeExpected = RANGE(1.0_RFREAL)
  
  READ(iFile) precActual,rangeActual
  IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN 
    CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
  END IF ! precActual
  
! -----------------------------------------------------------------------------
! Physical time
! -----------------------------------------------------------------------------
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile) currentTime 
 
  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    IF ( global%currentTime < 0.0_RFREAL ) THEN
      global%currentTime = currentTime
    ELSE
      WRITE(timeString1,'(1PE11.5)') global%currentTime
      WRITE(timeString2,'(1PE11.5)') currentTime          
      IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
        CALL ErrorStop(global,ERR_TIME_SOLUTION,__LINE__,TRIM(iFileName))
      END IF ! global%currentTime 
    END IF ! global%currentTime
  END IF ! global%flowType  
  
! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  
  nVarsExpected  = pRegion%specInput%nSpecies
  nCellsExpected = pGrid%nCellsTot    
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM
    
  READ(iFile) nCellsTot,nVars
  IF ( nCellsTot /= nCellsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, & 
                                              'but expected:',nCellsExpected
    CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,errorString)
  END IF ! nCellsExpected     
  
  IF ( nVars /= nVarsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, & 
                                              'but expected:',nVarsExpected  
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVarsExpected   
  
! ==============================================================================
! Rest of file
! ==============================================================================

  iVars       = 0
  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1
  
    READ(iFile) sectionString

    SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!     Species concentration
! ------------------------------------------------------------------------------

      CASE ( '# Density' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species concentration...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%spec%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(iVars,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     End marker
! ------------------------------------------------------------------------------ 
      
      CASE ( '# End' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel           
      
        EXIT
      
! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------ 
      
      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel           
      
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)        
             
    END SELECT ! TRIM
  
! ------------------------------------------------------------------------------
!   Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------  
  
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  
  END DO ! <empty>

! ==============================================================================
! Check and information about number of variables read
! ==============================================================================

  IF ( iVars /= nVars ) THEN 
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! iVar

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)   
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
    
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary species cv file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
   
END SUBROUTINE SPEC_RFLU_ReadCvBinary









! ******************************************************************************
!
! Purpose: Read eev file for species in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_ReadEEvASCII(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,iFileName,sectionString,timeString1, &
                       timeString2
  INTEGER :: errorFlag,iFile,iSpecEEvTemp,iSpecEEvXVel,iSpecEEvYVel, &
             iSpecEEvZVel,iVars,j,loopCounter,nCellsTot,nCellsExpected, & 
             nVars,nVarsExpected,precActual,precExpected,rangeActual, & 
             rangeExpected
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadEEvASCII',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII species eev file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.spec.eeva', & 
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.spec.eeva', & 
                             pRegion%iRegionGlobal,global%currentIter,iFileName)  
  END IF ! global%flowType

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
  global%error = errorFlag       
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# ROCFLU species eev file' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
! -----------------------------------------------------------------------------
! Precision and range
! -----------------------------------------------------------------------------
    
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
  precExpected  = PRECISION(1.0_RFREAL)
  rangeExpected = RANGE(1.0_RFREAL)
  
  READ(iFile,'(2(I8))') precActual,rangeActual
  IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN 
    CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
  END IF ! precActual
  
! -----------------------------------------------------------------------------
! Physical time
! -----------------------------------------------------------------------------
  
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile,'(E23.16)') currentTime 
 
  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    IF ( global%currentTime < 0.0_RFREAL ) THEN
      global%currentTime = currentTime
    ELSE
      WRITE(timeString1,'(1PE11.5)') global%currentTime
      WRITE(timeString2,'(1PE11.5)') currentTime          
      IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
        CALL ErrorStop(global,ERR_TIME_SOLUTION,__LINE__,TRIM(iFileName))
      END IF ! global%currentTime 
    END IF ! global%currentTime
  END IF ! global%flowType  
  
! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  
  nVarsExpected  = pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
  nCellsExpected = pGrid%nCellsTot    
  
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM
    
  READ(iFile,'(2(I8))') nCellsTot,nVars
  IF ( nCellsTot /= nCellsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, & 
                                              'but expected:',nCellsExpected
    CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,errorString)
  END IF ! nCellsExpected     
  
  IF ( nVars /= nVarsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, & 
                                              'but expected:',nVarsExpected  
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVarsExpected   
  
! ==============================================================================
! Rest of file
! ==============================================================================

  iSpecEEvXVel = 0
  iSpecEEvYVel = 0
  iSpecEEvZVel = 0    
  iSpecEEvTemp = 0
  loopCounter  = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1
  
    READ(iFile,'(A)') sectionString

    SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!     Species x-velocity
! ------------------------------------------------------------------------------

      CASE ( '# x-velocity' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species x-velocity...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvXVel = iSpecEEvXVel + 1
        READ(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_XVEL,iSpecEEvXVel,j), &
                                  j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Species y-velocity
! ------------------------------------------------------------------------------

      CASE ( '# y-velocity' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species y-velocity...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvYVel = iSpecEEvYVel + 1
        READ(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_YVEL,iSpecEEvYVel,j), &
                                  j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Species x-velocity
! ------------------------------------------------------------------------------

      CASE ( '# z-velocity' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species z-velocity...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvZVel = iSpecEEvZVel + 1
        READ(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_ZVEL,iSpecEEvZVel,j), &
                                  j=1,pGrid%nCellsTot)
                                  
! ------------------------------------------------------------------------------
!     Species temperature
! ------------------------------------------------------------------------------

      CASE ( '# Temperature' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species temperature...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvTemp = iSpecEEvTemp + 1
        READ(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_TEMP,iSpecEEvTemp,j), &
                                   j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     End marker
! ------------------------------------------------------------------------------ 
      
      CASE ( '# End' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel           
      
        EXIT
      
! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------ 
      
      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel           
      
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)        
                    
    END SELECT ! TRIM
  
! ------------------------------------------------------------------------------
!   Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------  
  
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  END DO ! <empty>

! ==============================================================================
! Check and information about number of variables read
! ==============================================================================

  nVars = iSpecEEvXVel + iSpecEEvYVel + iSpecEEvZVel + iSpecEEvTemp

  IF ( nVars /= nVarsExpected ) THEN 
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVars

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)   
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
    
! ******************************************************************************
! End
! ******************************************************************************
 
   IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII species eev file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
  
END SUBROUTINE SPEC_RFLU_ReadEEvASCII







! ******************************************************************************
!
! Purpose: Read eev file for species in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_ReadEEvBinary(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,iFileName,sectionString,timeString1, & 
                       timeString2
  INTEGER :: errorFlag,iFile,iSpecEEvTemp,iSpecEEvXVel,iSpecEEvYVel, &
             iSpecEEvZVel,iVars,j,loopCounter,nCellsTot,nCellsExpected, & 
             nVars,nVarsExpected,precActual,precExpected,rangeActual, & 
             rangeExpected
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadEEvBinary',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary species eev file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.spec.eev', & 
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.spec.eev', & 
                             pRegion%iRegionGlobal,global%currentIter,iFileName)  
  END IF ! global%flowType

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
  global%error = errorFlag       
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# ROCFLU species eev file' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
! -----------------------------------------------------------------------------
! Precision and range
! -----------------------------------------------------------------------------
    
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
  precExpected  = PRECISION(1.0_RFREAL)
  rangeExpected = RANGE(1.0_RFREAL)
  
  READ(iFile) precActual,rangeActual
  IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN 
    CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
  END IF ! precActual
  
! -----------------------------------------------------------------------------
! Physical time
! -----------------------------------------------------------------------------
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile) currentTime 
 
  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    IF ( global%currentTime < 0.0_RFREAL ) THEN
      global%currentTime = currentTime
    ELSE
      WRITE(timeString1,'(1PE11.5)') global%currentTime
      WRITE(timeString2,'(1PE11.5)') currentTime          
      IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
        CALL ErrorStop(global,ERR_TIME_SOLUTION,__LINE__,TRIM(iFileName))
      END IF ! global%currentTime 
    END IF ! global%currentTime
  END IF ! global%flowType  
  
! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  
  nVarsExpected  = pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
  nCellsExpected = pGrid%nCellsTot    
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM
    
  READ(iFile) nCellsTot,nVars
  IF ( nCellsTot /= nCellsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, & 
                                              'but expected:',nCellsExpected
    CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,errorString)
  END IF ! nCellsExpected     
  
  IF ( nVars /= nVarsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, & 
                                              'but expected:',nVarsExpected  
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVarsExpected   
  
! ==============================================================================
! Rest of file
! ==============================================================================

  iSpecEEvXVel = 0
  iSpecEEvYVel = 0
  iSpecEEvZVel = 0    
  iSpecEEvTemp = 0
  loopCounter  = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1
  
    READ(iFile) sectionString

    SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!     Species x-velocity
! ------------------------------------------------------------------------------

      CASE ( '# x-velocity' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species x-velocity...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvXVel = iSpecEEvXVel + 1
        READ(iFile) (pEEv(EEV_SPEC_XVEL,iSpecEEvXVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Species y-velocity
! ------------------------------------------------------------------------------

      CASE ( '# y-velocity' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species y-velocity...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvYVel = iSpecEEvYVel + 1
        READ(iFile) (pEEv(EEV_SPEC_YVEL,iSpecEEvYVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Species x-velocity
! ------------------------------------------------------------------------------

      CASE ( '# z-velocity' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species z-velocity...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvZVel = iSpecEEvZVel + 1
        READ(iFile) (pEEv(EEV_SPEC_ZVEL,iSpecEEvZVel,j),j=1,pGrid%nCellsTot)
                                  
! ------------------------------------------------------------------------------
!     Species temperature
! ------------------------------------------------------------------------------

      CASE ( '# Temperature' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Species temperature...'
        END IF ! global%verbLevel    
      
        pEEv => pRegion%spec%eev
      
        iSpecEEvTemp = iSpecEEvTemp + 1
        READ(iFile) (pEEv(EEV_SPEC_TEMP,iSpecEEvTemp,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     End marker
! ------------------------------------------------------------------------------ 
      
      CASE ( '# End' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel           
      
        EXIT
      
! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------ 
      
      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel           
      
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)        
                    
    END SELECT ! TRIM
  
! ------------------------------------------------------------------------------
!   Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------  
  
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  END DO ! <empty>

! ==============================================================================
! Check and information about number of variables read
! ==============================================================================

  nVars = iSpecEEvXVel + iSpecEEvYVel + iSpecEEvZVel + iSpecEEvTemp

  IF ( nVars /= nVarsExpected ) THEN 
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVars

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)   
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
    
! ******************************************************************************
! End
! ******************************************************************************
 
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary species eev file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
   
END SUBROUTINE SPEC_RFLU_ReadEEvBinary








! ******************************************************************************
!
! Purpose: Write flow file for species in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Write physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_WriteCvASCII(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: errorFlag,iFile,iVar,j,nVars
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_WriteCvASCII',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII species cv file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.spec.cva', & 
                               pRegion%iRegionGlobal,global%currentTime, & 
                               iFileName)  
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.spec.cva', & 
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)
  END IF ! global%flowType

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  sectionString = '# ROCFLU species file'
  WRITE(iFile,'(A)') sectionString
  
  sectionString = '# Precision and range'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)
  
  sectionString = '# Physical time'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(E23.16)') global%currentTime 

! ==============================================================================
! Dimensions
! ==============================================================================
  
  nVars = pRegion%specInput%nSpecies
  
  pGrid => pRegion%grid  
  
  sectionString = '# Dimensions'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(2(I8))') pGrid%nCellsTot,nVars 
  
! ==============================================================================
! Species concentration
! ==============================================================================

  pCv => pRegion%spec%cv

  DO iVar = 1,nVars
    sectionString = '# Density'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pCv(iVar,j),j=1,pGrid%nCellsTot)
  END DO ! iVar
 
! ==============================================================================
! End marker
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF ! global%verbLevel

  sectionString = '# End'
  WRITE(iFile,'(A)') sectionString  

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
    
! ******************************************************************************
! End
! ******************************************************************************
 
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII species cv file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
   
END SUBROUTINE SPEC_RFLU_WriteCvASCII







! ******************************************************************************
!
! Purpose: Write flow file for species in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Write physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_WriteCvBinary(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: errorFlag,iFile,iVar,j,nVars
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_WriteCvBinary',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary species cv  file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.spec.cv', & 
                               pRegion%iRegionGlobal,global%currentTime, & 
                               iFileName)                             
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.spec.cv', & 
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)                                
  END IF ! global%flowType

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  sectionString = '# ROCFLU species file'
  WRITE(iFile) sectionString
  
  sectionString = '# Precision and range'
  WRITE(iFile) sectionString
  WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)
  
  sectionString = '# Physical time'
  WRITE(iFile) sectionString
  WRITE(iFile) global%currentTime 

! ==============================================================================
! Dimensions
! ==============================================================================
  
  nVars = pRegion%specInput%nSpecies
  
  pGrid => pRegion%grid  
  
  sectionString = '# Dimensions'
  WRITE(iFile) sectionString
  WRITE(iFile) pGrid%nCellsTot,nVars 
  
! ==============================================================================
! Species concentration
! ==============================================================================

  pCv => pRegion%spec%cv

  DO iVar = 1,nVars
    sectionString = '# Density'
    WRITE(iFile) sectionString
    WRITE(iFile) (pCv(iVar,j),j=1,pGrid%nCellsTot)
  END DO ! iVar
 
! ==============================================================================
! End marker
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF ! global%verbLevel

  sectionString = '# End'
  WRITE(iFile) sectionString  

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
   
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary species cv file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
 
END SUBROUTINE SPEC_RFLU_WriteCvBinary






! ******************************************************************************
!
! Purpose: Write eev file for species in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Write physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_WriteEEvASCII(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: errorFlag,iFile,iSpecEE,j,nVars
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_WriteEEvASCII',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII species eev file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.spec.eeva', & 
                               pRegion%iRegionGlobal,global%currentTime, & 
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.spec.eeva', & 
                             pRegion%iRegionGlobal,global%currentIter,iFileName)  
  ENDIF ! global%flowType

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)
  global%error = errorFlag          
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  sectionString = '# ROCFLU species eev file'
  WRITE(iFile,'(A)') sectionString
  
  sectionString = '# Precision and range'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)
  
  sectionString = '# Physical time'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(E23.16)') global%currentTime 

! ==============================================================================
! Dimensions
! ==============================================================================
  
  nVars = pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
  
  pGrid => pRegion%grid  
  
  sectionString = '# Dimensions'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(2(I8))') pGrid%nCellsTot,nVars 
  
! ==============================================================================
! Variables
! ==============================================================================

  pEEv => pRegion%spec%eev

  DO iSpecEE = 1,pRegion%specInput%nSpeciesEE
    sectionString = '# x-velocity'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_XVEL,iSpecEE,j), &
                               j=1,pGrid%nCellsTot)
                            
    sectionString = '# y-velocity'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_YVEL,iSpecEE,j), &
                               j=1,pGrid%nCellsTot)
                               
    sectionString = '# z-velocity'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_ZVEL,iSpecEE,j), &
                               j=1,pGrid%nCellsTot)  
                               
    sectionString = '# Temperature'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pEEv(EEV_SPEC_TEMP,iSpecEE,j), &
                               j=1,pGrid%nCellsTot)                                                                                            
  END DO ! iSpecEE
 
! ==============================================================================
! End marker
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF ! global%verbLevel

  sectionString = '# End'
  WRITE(iFile,'(A)') sectionString  

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
       
! ******************************************************************************
! End
! ******************************************************************************
 
   IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII species eev file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_WriteEEvASCII





! ******************************************************************************
!
! Purpose: Write eev file for species in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Write physical time for both steady and unsteady flows so that could 
!      use steady solution as input for unsteady run and vice versa.
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_WriteEEvBinary(pRegion)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: errorFlag,iFile,iSpecEE,j,nVars
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_WriteEEvBinary',&
  'SPEC_RFLU_ModReadWriteVars.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary species eev file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.spec.eev', & 
                               pRegion%iRegionGlobal,global%currentTime, & 
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.spec.eev', & 
                             pRegion%iRegionGlobal,global%currentIter,iFileName)  
  ENDIF ! global%flowType

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)
  global%error = errorFlag          
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  sectionString = '# ROCFLU species eev file'
  WRITE(iFile) sectionString
  
  sectionString = '# Precision and range'
  WRITE(iFile) sectionString
  WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)
  
  sectionString = '# Physical time'
  WRITE(iFile) sectionString
  WRITE(iFile) global%currentTime 

! ==============================================================================
! Dimensions
! ==============================================================================
  
  nVars = pRegion%specInput%nSpeciesEE*EEV_SPEC_NVAR
  
  pGrid => pRegion%grid  
  
  sectionString = '# Dimensions'
  WRITE(iFile) sectionString
  WRITE(iFile) pGrid%nCellsTot,nVars 
  
! ==============================================================================
! Variables
! ==============================================================================

  pEEv => pRegion%spec%eev

  DO iSpecEE = 1,pRegion%specInput%nSpeciesEE
    sectionString = '# x-velocity'
    WRITE(iFile) sectionString
    WRITE(iFile) (pEEv(EEV_SPEC_XVEL,iSpecEE,j),j=1,pGrid%nCellsTot)
                            
    sectionString = '# y-velocity'
    WRITE(iFile) sectionString
    WRITE(iFile) (pEEv(EEV_SPEC_YVEL,iSpecEE,j),j=1,pGrid%nCellsTot)
                               
    sectionString = '# z-velocity'
    WRITE(iFile) sectionString
    WRITE(iFile) (pEEv(EEV_SPEC_ZVEL,iSpecEE,j),j=1,pGrid%nCellsTot)  
                               
    sectionString = '# Temperature'
    WRITE(iFile) sectionString
    WRITE(iFile) (pEEv(EEV_SPEC_TEMP,iSpecEE,j),j=1,pGrid%nCellsTot)                                                                                            
  END DO ! iSpecEE
 
! ==============================================================================
! End marker
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF ! global%verbLevel

  sectionString = '# End'
  WRITE(iFile) sectionString  

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
       
! ******************************************************************************
! End
! ******************************************************************************
 
   IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII species eev file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_WriteEEvBinary




! ******************************************************************************
! End
! ******************************************************************************

END MODULE SPEC_RFLU_ModReadWriteVars


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ModReadWriteVars.F90,v $
! Revision 1.4  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/05 12:44:52  haselbac
! Removed superfluous close parentheses - found by ifort compiler
!
! Revision 1.1  2005/11/27 01:47:26  haselbac
! Initial revision
!
! ******************************************************************************














