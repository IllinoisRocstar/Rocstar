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
! Purpose: Suite of routines to read and write mixture solution files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadWriteFlow.F90,v 1.13 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadWriteFlow

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

  USE ModInterfaces, ONLY: RFLU_GetCvLoc

  IMPLICIT NONE
  INCLUDE 'comf90.h'

  PRIVATE
  PUBLIC :: RFLU_ReadFlowWrapper, &
            RFLU_WriteFlowWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString = &
    '$RCSfile: RFLU_ModReadWriteFlow.F90,v $ $Revision: 1.13 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Read flow file for mixture in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowASCII(pRegion)

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

    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: errorString,iFileName,iFileNameOld,sectionString, &
                         timeString1,timeString2
    INTEGER :: cvMixtDens,cvMixtEner,cvMixtPres,cvMixtXMom,cvMixtXVel, & 
               cvMixtYMom,cvMixtYVel,cvMixtZMom,cvMixtZVel,errorFlag,i, &
               iFile,iVars,j,loopCounter,nCellsTot,nCellsExpected,nVars, &
               nVarsExpected,precActual,precExpected,rangeActual, &
               rangeExpected
    REAL(RFREAL) :: currentTime
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowASCII',&
  'RFLU_ModReadWriteFlow.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.mixt.cva', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName) 
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.floa', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileNameOld)
                                
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.mixt.cva', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)      
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.floa', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileNameOld)                             

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
       
    INQUIRE(FILE=iFileName,EXIST=fileExists)
    
    IF ( fileExists .EQV. .TRUE. ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
      END IF ! global%error
    ELSE 
      OPEN(iFile,FILE=iFileNameOld,FORM="FORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileNameOld)
      END IF ! global%error    
    END IF ! fileExists

! ==============================================================================
!   Set state vector state depending on fluid model
! ==============================================================================

    SELECT CASE ( pRegion%mixtInput%fluidModel ) 
      CASE ( FLUID_MODEL_COMP )
        pRegion%mixt%cvState = CV_MIXT_STATE_CONS
      CASE ( FLUID_MODEL_INCOMP ) 
        pRegion%mixt%cvState = CV_MIXT_STATE_PRIM
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%fluidModel

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU flow file' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
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
!   Initial residual and physical time
! -----------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Initial residual' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile,'(E23.16)') global%resInit

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
          WRITE(STDOUT,'(A,4(1X,A))') SOLVER_NAME, & 
                                      '*** WARNING *** Time mismatch:', & 
                                      TRIM(timeString1),'vs.',TRIM(timeString2)
                                      
          global%warnCounter = global%warnCounter + 1
        END IF ! global%currentTime
      END IF ! global%currentTime
    END IF ! global%flowType  

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    nVarsExpected  = pRegion%mixtInput%nCv
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
!   Rest of file
! ==============================================================================

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       Mixture density - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture density' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture density...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtDens = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_DENS)          
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtDens,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtXMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_XMOM)
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtXMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)          
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtXVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture y-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtYMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_YMOM)            
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtYMom,j),j=1,pGrid%nCellsTot)
          
! ------------------------------------------------------------------------------
!       Mixture y-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)          
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtYVel,j),j=1,pGrid%nCellsTot)          

! ------------------------------------------------------------------------------
!       Mixture z-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtZMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ZMOM)            
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtZMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture z-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)          
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtZVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture total internal energy - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture total internal energy' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture total internal '// &
                                     'energy...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtEner = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ENER)            
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtEner,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture pressure - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture pressure' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture pressure...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)          
          
          READ(iFile,'(5(E23.16))') (pCv(cvMixtPres,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( '# End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------

        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)

      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
    END IF ! iVar

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowASCII







! ******************************************************************************
!
! Purpose: Read flow file for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowBinary(pRegion)

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

    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: errorString,iFileName,iFileNameOld,sectionString, &
                         timeString1,timeString2
    INTEGER :: cvMixtDens,cvMixtEner,cvMixtPres,cvMixtXMom,cvMixtXVel, & 
               cvMixtYMom,cvMixtYVel,cvMixtZMom,cvMixtZVel,errorFlag,i, &
               iFile,iVars,j,loopCounter,nCellsTot,nCellsExpected,nVars, &
               nVarsExpected,precActual,precExpected,rangeActual, &
               rangeExpected
    REAL(RFREAL) :: currentTime
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowBinary',&
  'RFLU_ModReadWriteFlow.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.mixt.cv', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.flo', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileNameOld)                                 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.mixt.cv', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.flo', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileNameOld)                               

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
       
    INQUIRE(FILE=iFileName,EXIST=fileExists)
    
    IF ( fileExists .EQV. .TRUE. ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
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
!   Set state vector state depending on fluid model
! ==============================================================================

    SELECT CASE ( pRegion%mixtInput%fluidModel ) 
      CASE ( FLUID_MODEL_INCOMP ) 
        pRegion%mixt%cvState = CV_MIXT_STATE_PRIM
      CASE ( FLUID_MODEL_COMP )
        pRegion%mixt%cvState = CV_MIXT_STATE_CONS
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%fluidModel

! ==============================================================================
! Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU flow file' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
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
!   Initial residual and physical time
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Initial residual' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile) global%resInit

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
          WRITE(STDOUT,'(A,4(1X,A))') SOLVER_NAME, & 
                                      '*** WARNING *** Time mismatch:', & 
                                      TRIM(timeString1),'vs.',TRIM(timeString2)
                                      
          global%warnCounter = global%warnCounter + 1
        END IF ! global%currentTime
      END IF ! global%currentTime
    END IF ! global%flowType  

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    nVarsExpected  = pRegion%mixtInput%nCv
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
!   Rest of file
! ==============================================================================

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       Mixture density - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture density' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture density...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtDens = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_DENS)          
          
          READ(iFile) (pCv(cvMixtDens,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtXMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_XMOM)
          
          READ(iFile) (pCv(cvMixtXMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)          
          
          READ(iFile) (pCv(cvMixtXVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture y-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtYMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_YMOM)            
          
          READ(iFile) (pCv(cvMixtYMom,j),j=1,pGrid%nCellsTot)
          
! ------------------------------------------------------------------------------
!       Mixture y-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)          
          
          READ(iFile) (pCv(cvMixtYVel,j),j=1,pGrid%nCellsTot)          

! ------------------------------------------------------------------------------
!       Mixture z-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtZMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ZMOM)            
          
          READ(iFile) (pCv(cvMixtZMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture z-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)          
          
          READ(iFile) (pCv(cvMixtZVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture total internal energy - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture total internal energy' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture total internal '// &
                                     'energy...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtEner = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ENER)            
          
          READ(iFile) (pCv(cvMixtEner,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture pressure - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture pressure' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture pressure...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1
          
          cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)          
          
          READ(iFile) (pCv(cvMixtPres,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( '# End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------

        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)

      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter

    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
    END IF ! iVar

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowBinary








! ******************************************************************************
!
! Purpose: Wrapper for reading of flow files in ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowWrapper(pRegion)

#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_DecideReadFile, & 
                              RFLU_GENX_GetDataFlow
#endif

#ifdef PLAG
    USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_ReadSolutionASCII, &
                                       PLAG_RFLU_ReadSolutionBinary
#endif

#ifdef SPEC
    USE SPEC_RFLU_ModReadWriteVars, ONLY: SPEC_RFLU_ReadCvASCII, &
                                          SPEC_RFLU_ReadCvBinary, &     
                                          SPEC_RFLU_ReadEEvASCII, & 
                                          SPEC_RFLU_ReadEEvBinary
#endif

#ifdef TURB
    USE ModInterfacesTurbulence, ONLY: TURB_RFLU_ReadSolutionASCII, &
                                       TURB_RFLU_ReadSolutionBinary, &
                                       TURB_InitSolution
#ifdef GENX
    USE TURB_RFLU_ModRocstarAdmin, ONLY: TURB_RFLU_GenxGetData
#endif
#endif

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

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowWrapper',&
  'RFLU_ModReadWriteFlow.F90')

! ******************************************************************************
!   Read mixture solution files
! ******************************************************************************

#ifdef GENX
    IF ( RFLU_GENX_DecideReadFile(global) .EQV. .FALSE. ) THEN 
#endif
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL RFLU_ReadFlowASCII(pRegion)
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL RFLU_ReadFlowBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
#ifdef GENX
    ELSE 
      CALL RFLU_GENX_GetDataFlow(pRegion) 
    END IF ! RFLU_GENX_DecideReadFile         
#endif

! ******************************************************************************
!   Read physical module solution files
! ******************************************************************************

#ifdef PLAG
! ==============================================================================
!   Particles
! ==============================================================================

    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL PLAG_RFLU_ReadSolutionASCII(pRegion)
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL PLAG_RFLU_ReadSolutionBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
    END IF ! plagUsed
#endif

#ifdef SPEC
! ==============================================================================
!   Species
! ==============================================================================

    IF ( global%specUsed .EQV. .TRUE. ) THEN
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL SPEC_RFLU_ReadCvASCII(pRegion)
                
        IF ( global%moduleType == MODULE_TYPE_POSTPROC ) THEN
          IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN  
            CALL SPEC_RFLU_ReadEEvASCII(pRegion)
          END IF ! pRegion%specInput%nSpeciesEE
        END IF ! global%moduleType 
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL SPEC_RFLU_ReadCvBinary(pRegion)
        
        IF ( global%moduleType == MODULE_TYPE_POSTPROC ) THEN
          IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN  
            CALL SPEC_RFLU_ReadEEvBinary(pRegion)
          END IF ! pRegion%specInput%nSpeciesEE
        END IF ! global%moduleType        
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
    END IF ! global%specUsed
#endif

#ifdef TURB
! ==============================================================================
!   Turbulence
! ==============================================================================

    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST .AND. &
         pRegion%mixtInput%turbModel /= TURB_MODEL_NONE ) THEN
      IF ( (global%flowType == FLOW_UNSTEADY .AND. &
            global%currentTime > 0._RFREAL)     .OR. &
           (global%flowType == FLOW_STEADY   .AND. &
            global%currentIter > 0) ) THEN
#ifndef GENX
        IF ( global%solutFormat == FORMAT_ASCII ) THEN
          CALL TURB_RFLU_ReadSolutionASCII(pRegion)
        ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
          CALL TURB_RFLU_ReadSolutionBinary(pRegion)
        ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! global%solutFormat
#else
        CALL TURB_RFLU_GenxGetData(pRegion)
#endif
      END IF ! global%flowType
      CALL TURB_InitSolution(pRegion)  ! temporary
    END IF ! pRegion%mixtInput%flowModel
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowWrapper








! ******************************************************************************
!
! Purpose: Write flow file for mixture in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteFlowASCII(pRegion)

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
    INTEGER :: cvMixtDens,cvMixtEner,cvMixtPres,cvMixtXMom,cvMixtXVel, & 
               cvMixtYMom,cvMixtYVel,cvMixtZMom,cvMixtZVel,errorFlag,i,iFile,j
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteFlowASCII',&
  'RFLU_ModReadWriteFlow.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.mixt.cva', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.mixt.cva', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT      
    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU flow file'
    WRITE(iFile,'(A)') sectionString

    sectionString = '# Precision and range'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Initial residual'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(E23.16)') global%resInit

    sectionString = '# Physical time'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(E23.16)') global%currentTime
    
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(2(I8))') pGrid%nCellsTot,pRegion%mixtInput%nCv

! ==============================================================================
!   Data
! ==============================================================================

    pCv => pRegion%mixt%cv

    SELECT CASE ( pRegion%mixtInput%fluidModel )

! ------------------------------------------------------------------------------
!     Compressible fluid 
! ------------------------------------------------------------------------------

      CASE ( FLUID_MODEL_COMP )
      
! ----- Density ----------------------------------------------------------------                   

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture density...'
        END IF ! global%verbLevel

        cvMixtDens = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_DENS)

        sectionString = '# Mixture density'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtDens,j),j=1,pGrid%nCellsTot)

! ----- X-momentum -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-momentum...'
        END IF ! global%verbLevel

        cvMixtXMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_XMOM)

        sectionString = '# Mixture x-momentum'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtXMom,j),j=1,pGrid%nCellsTot)

! ----- Y-momentum -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-momentum...'
        END IF ! global%verbLevel

        cvMixtYMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_YMOM)

        sectionString = '# Mixture y-momentum'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtYMom,j),j=1,pGrid%nCellsTot)

! ----- Z-momentum -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-momentum...'
        END IF ! global%verbLevel

        cvMixtZMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ZMOM)

        sectionString = '# Mixture z-momentum'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtZMom,j),j=1,pGrid%nCellsTot)

! ----- Total internal energy --------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                   'Mixture total internal energy...'
        END IF ! global%verbLevel

        cvMixtEner = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ENER)

        sectionString = '# Mixture total internal energy'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtEner,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Incompressible fluid 
! ------------------------------------------------------------------------------

      CASE ( FLUID_MODEL_INCOMP )

! ----- X-velocity -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-velocity...'
        END IF ! global%verbLevel

        cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)

        sectionString = '# Mixture x-velocity'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtXVel,j),j=1,pGrid%nCellsTot)

! ----- Y-velocity -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-velocity...'
        END IF ! global%verbLevel

        cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)

        sectionString = '# Mixture y-velocity'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtYVel,j),j=1,pGrid%nCellsTot)

! ----- Z-velocity -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-velocity...'
        END IF ! global%verbLevel

        cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)

        sectionString = '# Mixture z-velocity'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtZVel,j),j=1,pGrid%nCellsTot)

! ----- Pressure ---------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                   'Mixture pressure...'
        END IF ! global%verbLevel

        cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)

        sectionString = '# Mixture pressure'
        WRITE(iFile,'(A)') sectionString
        WRITE(iFile,'(5(E23.16))') (pCv(cvMixtPres,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Default
! ------------------------------------------------------------------------------

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%fluidModel            
                
! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') sectionString

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteFlowASCII







! ******************************************************************************
!
! Purpose: Write flow file for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteFlowBinary(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: cvMixtDens,cvMixtEner,cvMixtPres,cvMixtXMom,cvMixtXVel, & 
               cvMixtYMom,cvMixtYVel,cvMixtZMom,cvMixtZVel,errorFlag,i,iFile,j
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteFlowBinary',&
  'RFLU_ModReadWriteFlow.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.mixt.cv', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.mixt.cv', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU flow file'
    WRITE(iFile) sectionString

    sectionString = '# Precision and range'
    WRITE(iFile) sectionString
    WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Initial residual'
    WRITE(iFile) sectionString
    WRITE(iFile) global%resInit

    sectionString = '# Physical time'
    WRITE(iFile) sectionString
    WRITE(iFile) global%currentTime
    
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    sectionString = '# Dimensions'
    WRITE(iFile) sectionString
    WRITE(iFile) pGrid%nCellsTot,pRegion%mixtInput%nCv

! ==============================================================================
!   Data
! ==============================================================================

    pCv => pRegion%mixt%cv

    SELECT CASE ( pRegion%mixtInput%fluidModel )

! ------------------------------------------------------------------------------
!     Compressible fluid 
! ------------------------------------------------------------------------------

      CASE ( FLUID_MODEL_COMP )
      
! ----- Density ----------------------------------------------------------------                   

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture density...'
        END IF ! global%verbLevel

        cvMixtDens = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_DENS)

        sectionString = '# Mixture density'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtDens,j),j=1,pGrid%nCellsTot)

! ----- X-momentum -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-momentum...'
        END IF ! global%verbLevel

        cvMixtXMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_XMOM)

        sectionString = '# Mixture x-momentum'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtXMom,j),j=1,pGrid%nCellsTot)

! ----- Y-momentum -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-momentum...'
        END IF ! global%verbLevel

        cvMixtYMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_YMOM)

        sectionString = '# Mixture y-momentum'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtYMom,j),j=1,pGrid%nCellsTot)

! ----- Z-momentum -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-momentum...'
        END IF ! global%verbLevel

        cvMixtZMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ZMOM)

        sectionString = '# Mixture z-momentum'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtZMom,j),j=1,pGrid%nCellsTot)

! ----- Total internal energy --------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                   'Mixture total internal energy...'
        END IF ! global%verbLevel

        cvMixtEner = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ENER)

        sectionString = '# Mixture total internal energy'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtEner,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Incompressible fluid 
! ------------------------------------------------------------------------------

      CASE ( FLUID_MODEL_INCOMP )

! ----- X-velocity -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-velocity...'
        END IF ! global%verbLevel

        cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)

        sectionString = '# Mixture x-velocity'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtXVel,j),j=1,pGrid%nCellsTot)

! ----- Y-velocity -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-velocity...'
        END IF ! global%verbLevel

        cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)

        sectionString = '# Mixture y-velocity'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtYVel,j),j=1,pGrid%nCellsTot)

! ----- Z-velocity -------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-velocity...'
        END IF ! global%verbLevel

        cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)

        sectionString = '# Mixture z-velocity'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtZVel,j),j=1,pGrid%nCellsTot)

! ----- Pressure ---------------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                   'Mixture pressure...'
        END IF ! global%verbLevel

        cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)

        sectionString = '# Mixture pressure'
        WRITE(iFile) sectionString
        WRITE(iFile) (pCv(cvMixtPres,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!     Default
! ------------------------------------------------------------------------------

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%fluidModel  

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteFlowBinary








! ******************************************************************************
!
! Purpose: Wrapper for writing of flow files in ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteFlowWrapper(pRegion)

#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_DecideWriteFile, & 
                              RFLU_GENX_PutDataFlow
#endif

#ifdef PLAG
    USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_WriteSolutionASCII, &
                                       PLAG_RFLU_WriteSolutionBinary
#endif

#ifdef SPEC
    USE SPEC_RFLU_ModReadWriteVars, ONLY: SPEC_RFLU_WriteCvASCII, &
                                          SPEC_RFLU_WriteCvBinary, & 
                                          SPEC_RFLU_WriteEEvASCII, &
                                          SPEC_RFLU_WriteEEvBinary
#endif

#ifdef TURB
    USE ModInterfacesTurbulence, ONLY: TURB_RFLU_WriteSolutionASCII, &
                                       TURB_RFLU_WriteSolutionBinary
#endif

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

    TYPE(t_global), POINTER :: global
    INTEGER :: sz, ng 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteFlowWrapper',&
  'RFLU_ModReadWriteFlow.F90')

! ******************************************************************************
!   Write mixture solution files
! ******************************************************************************

#ifdef GENX
    IF ( RFLU_GENX_DecideWriteFile(global) .EQV. .FALSE. ) THEN 
#endif
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL RFLU_WriteFlowASCII(pRegion)
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL RFLU_WriteFlowBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
#ifdef GENX
    ELSE 
      CALL COM_get_size(TRIM(global%volWinName)//'.rhof', 101, sz, ng)
      CALL RFLU_GENX_PutDataFlow(pRegion) 
    END IF ! RFLU_GENX_DecideReadFile         
#endif

! ******************************************************************************
!   Write physical module solution files
! ******************************************************************************

#ifdef PLAG
! ==============================================================================
!   Particles
! ==============================================================================

    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL PLAG_RFLU_WriteSolutionASCII(pRegion)
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL PLAG_RFLU_WriteSolutionBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
    END IF ! plagUsed
#endif

#ifdef SPEC
! ==============================================================================
!   Species
! ==============================================================================

    IF ( global%specUsed .EQV. .TRUE. ) THEN
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL SPEC_RFLU_WriteCvASCII(pRegion)
                       
        IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
          CALL SPEC_RFLU_WriteEEvASCII(pRegion)
        END IF ! pRegion%specInput%nSpeciesEE
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL SPEC_RFLU_WriteCvBinary(pRegion)

        IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
          CALL SPEC_RFLU_WriteEEvBinary(pRegion)
        END IF ! pRegion%specInput%nSpeciesEE
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
    END IF ! global%specUsed
#endif

#ifdef TURB
! ==============================================================================
!   Turbulence
! ==============================================================================

    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST .AND. &
         pRegion%mixtInput%turbModel /= TURB_MODEL_NONE ) THEN
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL TURB_RFLU_WriteSolutionASCII(pRegion)
      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
        CALL TURB_RFLU_WriteSolutionBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
    END IF ! NS-turb
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteFlowWrapper







! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModReadWriteFlow


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadWriteFlow.F90,v $
! Revision 1.13  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.12  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.11  2006/12/15 13:25:27  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.10  2006/03/26 20:22:07  haselbac
! Removed error traps for GL model
!
! Revision 1.9  2006/01/12 09:40:50  wasistho
! timeStamp to currentTime in turb readFlow
!
! Revision 1.8  2006/01/10 05:04:23  wasistho
! Get turbulence data from Genx
!
! Revision 1.7  2005/12/29 19:54:51  wasistho
! modified Rocturb part in ReadFlowWrapper
!
! Revision 1.6  2005/11/27 01:51:39  haselbac
! Added calls to EEv routines, changed extensions in backw-compatible way
!
! Revision 1.5  2004/11/06 03:18:26  haselbac
! Substantial additions to allow reading/writing of data for other fluid models
!
! Revision 1.4  2004/10/19 19:28:24  haselbac
! Adapted to changes in GENX logic
!
! Revision 1.3  2004/08/23 23:08:43  fnajjar
! Activated binary IO routine calls
!
! Revision 1.2  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.1  2004/07/06 15:14:31  haselbac
! Initial revision
!
! ******************************************************************************












