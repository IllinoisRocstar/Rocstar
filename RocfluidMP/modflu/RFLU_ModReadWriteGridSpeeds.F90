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
! Purpose: Suite of routines to read and write grid speed files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadWriteGridSpeeds.F90,v 1.3 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadWriteGridSpeeds

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
  PUBLIC :: RFLU_ReadGridSpeedsWrapper, & 
            RFLU_WriteGridSpeedsWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModReadWriteGridSpeeds.F90,v $ $Revision: 1.3 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS








! ******************************************************************************
!
! Purpose: Read grid speeds in ASCII ROCFLU format.
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

  SUBROUTINE RFLU_ReadGridSpeedsASCII(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, & 
                                 BuildFileNameUnsteady     

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

    CHARACTER(CHRLEN) :: iFileName,sectionString,timeString1,timeString2
    INTEGER :: errorFlag,iFile,ifg,ifl,iPatch,loopCounter,nFaces,p,r
    REAL(RFREAL) :: currentTime
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global  
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridSpeedsASCII',&
  'RFLU_ModReadWriteGridSpeeds.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII grid speeds file...'
    END IF ! global%verbLevel

    iFile = IF_GRID
  
    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.gspa', & 
                                 pRegion%iRegionGlobal,global%currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN    
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal       
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime
      END IF ! global%verbLevel                                      
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.gspa', & 
                              pRegion%iRegionGlobal,iFileName)    

        IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Read header stuff
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU grid speeds file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
! -----------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString  
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile,'(2(I8))') p,r
    IF ( p < PRECISION(1.0_RFREAL) .OR. r < RANGE(1.0_RFREAL) ) THEN 
      CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
    END IF ! p

! -----------------------------------------------------------------------------
!   Initial residual and physical time
! -----------------------------------------------------------------------------
  
    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM    
  
    READ(iFile,'(E23.16)') currentTime 

    IF ( global%flowType == FLOW_UNSTEADY .AND. & 
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN
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
!   Dimensions
! ==============================================================================

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM
  
    pGrid => pRegion%grid 

    READ(iFile,'(I8)') nFaces

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================

    IF ( nFaces /= pGrid%nFaces ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nFaces

! ==============================================================================
!   Rest of file
! ==============================================================================

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1
    
      READ(iFile,'(A)') sectionString
    
      SELECT CASE ( TRIM(sectionString) ) 
   
! ------------------------------------------------------------------------------
!       Grid speeds
! ------------------------------------------------------------------------------

        CASE ( '# Grid speeds' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH) THEN   
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Grid speeds...'
          END IF ! global%verbLevel

          IF ( pGrid%nFaces > 0 ) THEN 
            READ(iFile,'(5(E23.16))') (pGrid%gs(ifg),ifg=1,pGrid%nFaces)
          END IF ! pGrid%nFaces

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)  

            IF ( pPatch%nBFaces > 0 ) THEN 
              READ(iFile,'(5(E23.16))') (pPatch%gs(ifl),ifl=1,pPatch%nBFaces)
            END IF ! pPatch%nFaces
          END DO ! iPatch     
              
! ------------------------------------------------------------------------------
!     End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel       

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(3X,A)') sectionString
          END IF ! global%verbLevel        

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
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII grid speeds file done.'
    END IF ! global%verbLevel           

    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_ReadGridSpeedsASCII







! ******************************************************************************
!
! Purpose: Read grid speeds in binary ROCFLU format.
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

  SUBROUTINE RFLU_ReadGridSpeedsBinary(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, & 
                                 BuildFileNameUnsteady     

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

    CHARACTER(CHRLEN) :: iFileName,sectionString,timeString1,timeString2 
    INTEGER :: errorFlag,iFile,ifg,ifl,iPatch,loopCounter,nFaces,p,r
    REAL(RFREAL) :: currentTime
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global  
  
! ******************************************************************************
! Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridSpeedsBinary',&
  'RFLU_ModReadWriteGridSpeeds.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary grid speeds file...'
    END IF ! global%verbLevel

    iFile = IF_GRID
  
    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN 
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.gsp', & 
                                 pRegion%iRegionGlobal,global%currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN    
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal       
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime
      END IF ! global%verbLevel                                      
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.gsp', & 
                              pRegion%iRegionGlobal,iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Read header stuff
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU grid speeds file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile) p,r
    IF ( p < PRECISION(1.0_RFREAL) .OR. r < RANGE(1.0_RFREAL) ) THEN 
      CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
    END IF ! p

! -----------------------------------------------------------------------------
!   Initial residual and physical time
! -----------------------------------------------------------------------------
  
    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM    

    READ(iFile) currentTime 

    IF ( global%flowType == FLOW_UNSTEADY .AND. & 
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN
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
!   Dimensions
! ==============================================================================

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    pGrid => pRegion%grid 

    READ(iFile) nFaces

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================

    IF ( nFaces /= pGrid%nFaces ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nVertTot

! ==============================================================================
!   Rest of file
! ==============================================================================

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

      SELECT CASE ( TRIM(sectionString) ) 
   
! ------------------------------------------------------------------------------
!       Grid speeds
! ------------------------------------------------------------------------------

        CASE ( '# Grid speeds' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH) THEN   
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Grid speeds...'
          END IF ! global%verbLevel

          IF ( pGrid%nFaces > 0 ) THEN 
            READ(iFile) (pGrid%gs(ifg),ifg=1,pGrid%nFaces)
          END IF ! pGrid%nFaces

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)  

            IF ( pPatch%nBFaces > 0 ) THEN 
              READ(iFile) (pPatch%gs(ifl),ifl=1,pPatch%nBFaces)    
            END IF ! pPatch%nBFaces
          END DO ! iPatch       
      
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel       

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(3X,A)') sectionString
          END IF ! global%verbLevel        

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
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary grid speeds file done.'
    END IF ! global%verbLevel           

    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_ReadGridSpeedsBinary







! ******************************************************************************
!
! Purpose: Wrapper for reading of grid speeds files.
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

  SUBROUTINE RFLU_ReadGridSpeedsWrapper(pRegion)
  
#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_GetDataGSpeedsSurf, & 
                              RFLU_GENX_GetDataGSpeedsVol
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

    CALL RegisterFunction(global,'RFLU_ReadGridSpeedsWrapper',&
  'RFLU_ModReadWriteGridSpeeds.F90')

! ******************************************************************************
!   Read grid speed files
! ******************************************************************************

#ifndef GENX
    IF ( global%gridFormat == FORMAT_ASCII ) THEN 
      CALL RFLU_ReadGridSpeedsASCII(pRegion)
    ELSE IF ( global%gridFormat == FORMAT_BINARY ) THEN
      CALL RFLU_ReadGridSpeedsBinary(pRegion)
    ELSE  
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%gridFormat
#else
    CALL RFLU_GENX_GetDataGSpeedsSurf(pRegion)
    CALL RFLU_GENX_GetDataGSpeedsVol(pRegion)    
#endif
    
! ******************************************************************************
!   End
! ******************************************************************************
 
    CALL DeregisterFunction(global) 
 
  END SUBROUTINE RFLU_ReadGridSpeedsWrapper






! ******************************************************************************
!
! Purpose: Write grid speeds in ASCII ROCFLU format.
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

  SUBROUTINE RFLU_WriteGridSpeedsASCII(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, & 
                                 BuildFileNameUnsteady     

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
    INTEGER :: errorFlag,iFile,ifg,ifl,iPatch
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteGridSpeedsASCII',&
  'RFLU_ModReadWriteGridSpeeds.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII grid speeds file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN 
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.gspa', & 
                                 pRegion%iRegionGlobal,global%currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime   
      END IF ! global%verbLevel
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_OUTDIR,'.gspa', & 
                              pRegion%iRegionGlobal,iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal                                 
      END IF ! global%verbLevel
    END IF ! global

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
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    sectionString = '# ROCFLU grid speeds file'  
    WRITE(iFile,'(A)') TRIM(sectionString)  

    sectionString = '# Precision and range'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Physical time'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(E23.16)') global%currentTime 
  
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid  

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pGrid%nFaces

! ==============================================================================
!   Grid speeds
! ==============================================================================

    IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_MED) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Grid speeds...'
      END IF ! global%verbLevel

      sectionString = '# Grid speeds'
      WRITE(iFile,'(A)') TRIM(sectionString) 

      IF ( pGrid%nFaces > 0 ) THEN 
        WRITE(iFile,'(5(E23.16))') (pGrid%gs(ifg),ifg=1,pGrid%nFaces)
      END IF ! pGrid%nFaces

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)  

        IF ( pPatch%nBFaces > 0 ) THEN 
          WRITE(iFile,'(5(E23.16))') (pPatch%gs(ifl),ifl=1,pPatch%nBFaces)    
        END IF ! pPatch%nBFaces
      END DO ! iPatch    
    END IF ! pRegion%mixtInput%moveGrid

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString) 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
   
! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing ASCII grid speeds file done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_WriteGridSpeedsASCII






! ******************************************************************************
!
! Purpose: Write grid speeds in binary ROCFLU format.
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

  SUBROUTINE RFLU_WriteGridSpeedsBinary(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, & 
                                 BuildFileNameUnsteady     

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
    INTEGER :: errorFlag,iFile,ifg,ifl,iPatch
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteGridSpeedsBinary',&
  'RFLU_ModReadWriteGridSpeeds.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary grid speeds file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN 
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.gsp', & 
                                 pRegion%iRegionGlobal,global%currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime   
      END IF ! global%verbLevel
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_OUTDIR,'.gsp', & 
                              pRegion%iRegionGlobal,iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal                                 
      END IF ! global%verbLevel
    END IF ! global

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
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    sectionString = '# ROCFLU grid speeds file'  
    WRITE(iFile) sectionString  

    sectionString = '# Precision and range'
    WRITE(iFile) sectionString 
    WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Physical time'
    WRITE(iFile) sectionString 
    WRITE(iFile) global%currentTime 
  
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid  

    sectionString = '# Dimensions'
    WRITE(iFile) sectionString  
    WRITE(iFile) pGrid%nFaces
  
! ==============================================================================
!   Grid speeds
! ==============================================================================

    IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_MED) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Grid speeds...'
      END IF ! global%verbLevel

      sectionString = '# Grid speeds'
      WRITE(iFile) sectionString 

      IF ( pGrid%nFaces > 0 ) THEN   
        WRITE(iFile) (pGrid%gs(ifg),ifg=1,pGrid%nFaces)
      END IF ! pGrid%nFaces

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)  

        IF ( pPatch%nBFaces > 0 ) THEN 
          WRITE(iFile) (pPatch%gs(ifl),ifl=1,pPatch%nBFaces)    
        END IF ! pPatch%nBFaces
      END DO ! iPatch    
    END IF ! pRegion%mixtInput%moveGrid
  
! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
   
! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing binary grid speeds file done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_WriteGridSpeedsBinary






! ******************************************************************************
!
! Purpose: Wrapper for writing of grid speeds files in ROCFLU format.
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

  SUBROUTINE RFLU_WriteGridSpeedsWrapper(pRegion)

#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_DecideWriteFile, & 
                              RFLU_GENX_PutDataGSpeedsSurf, & 
                              RFLU_GENX_PutDataGSpeedsVol
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

    CALL RegisterFunction(global,'RFLU_WriteGridSpeedsWrapper',&
  'RFLU_ModReadWriteGridSpeeds.F90')

! ******************************************************************************
!   Read grid speed file
! ******************************************************************************

#ifdef GENX
    IF ( RFLU_GENX_DecideWriteFile(global) .EQV. .FALSE. ) THEN 
#endif
      IF ( global%gridFormat == FORMAT_ASCII ) THEN 
        CALL RFLU_WriteGridSpeedsASCII(pRegion)
      ELSE IF ( global%gridFormat == FORMAT_BINARY ) THEN
        CALL RFLU_WriteGridSpeedsBinary(pRegion)
      ELSE  
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%gridFormat
#ifdef GENX
    ELSE 
      CALL RFLU_GENX_PutDataGSpeedsSurf(pRegion)
      CALL RFLU_GENX_PutDataGSpeedsVol(pRegion) 
    END IF ! RFLU_GENX_DecideReadFile         
#endif
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
 
  END SUBROUTINE RFLU_WriteGridSpeedsWrapper





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModReadWriteGridSpeeds


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadWriteGridSpeeds.F90,v $
! Revision 1.3  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/10/19 19:27:10  haselbac
! Initial revision
!
! ******************************************************************************












