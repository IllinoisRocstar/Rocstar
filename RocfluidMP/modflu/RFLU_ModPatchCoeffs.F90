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
! Purpose: Collection of routines for patch coefficients.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModPatchCoeffs.F90,v 1.6 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPatchCoeffs

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI

  USE ModBuildFileNames, ONLY: BuildFileNameSteady, &
                               BuildFileNameUnsteady

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPatchCoeffs.F90,v $ $Revision: 1.6 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: RFLU_CreatePatchCoeffs, & 
            RFLU_DestroyPatchCoeffs, & 
            RFLU_ReadPatchCoeffsWrapper, &
            RFLU_WritePatchCoeffsWrapper

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_ClosePatchCoeffs, & 
             RFLU_InitPatchCoeffs, &
             RFLU_NullifyPatchCoeffs, &
             RFLU_OpenPatchCoeffsASCII, &
             RFLU_OpenPatchCoeffsBinary, &
             RFLU_ReadPatchCoeffsASCII, & 
             RFLU_ReadPatchCoeffsBinary, & 
             RFLU_WritePatchCoeffsASCII, & 
             RFLU_WritePatchCoeffsBinary

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Close patch-coefficients file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ClosePatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ClosePatchCoeffs',&
  'RFLU_ModPatchCoeffs.F90')
             
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Closing patch-coefficients file...'
    END IF ! global%verbLevel             
             
! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(IF_PATCH_COEF,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Closing patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ClosePatchCoeffs  
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Create patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreatePatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_CreatePatchCoeffs',&
  'RFLU_ModPatchCoeffs.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      ALLOCATE(pPatch%cp(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cp')
      END IF ! global%error

      ALLOCATE(pPatch%cf(XCOORD:ZCOORD,pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cf')
      END IF ! global%error
      
      ALLOCATE(pPatch%ch(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%ch')
      END IF ! global%error                              
      
      ALLOCATE(pPatch%cmass(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cmass')
      END IF ! global%error
      
      ALLOCATE(pPatch%cmom(XCOORD:ZCOORD,pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cmom')
      END IF ! global%error
    END DO ! iPatch

! ******************************************************************************
!   Initialize memory
! ******************************************************************************

    CALL RFLU_InitPatchCoeffs(pRegion)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreatePatchCoeffs
 








! *******************************************************************************
!
! Purpose: Destroy patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyPatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_DestroyPatchCoeffs',&
  'RFLU_ModPatchCoeffs.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      DEALLOCATE(pPatch%cp,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cp')
      END IF ! global%error

      DEALLOCATE(pPatch%cf,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cf')
      END IF ! global%error
      
      DEALLOCATE(pPatch%ch,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%ch')
      END IF ! global%error              
      
      DEALLOCATE(pPatch%cmass,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cmass')
      END IF ! global%error
      
      DEALLOCATE(pPatch%cmom,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cmom')
      END IF ! global%error
    END DO ! iPatch

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyPatchCoeffs(pRegion)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyPatchCoeffs





! *******************************************************************************
!
! Purpose: Initialize patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_InitPatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,ifl,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_InitPatchCoeffs',&
  'RFLU_ModPatchCoeffs.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      DO ifl = 1,pPatch%nBFaces
        pPatch%cp(ifl)        = 0.0_RFREAL
        pPatch%cf(XCOORD,ifl) = 0.0_RFREAL
        pPatch%cf(YCOORD,ifl) = 0.0_RFREAL
        pPatch%cf(ZCOORD,ifl) = 0.0_RFREAL
        pPatch%ch(ifl)        = 0.0_RFREAL                
        pPatch%cmass(ifl)     = 0.0_RFREAL                
        pPatch%cmom(XCOORD,ifl) = 0.0_RFREAL                
        pPatch%cmom(YCOORD,ifl) = 0.0_RFREAL                
        pPatch%cmom(ZCOORD,ifl) = 0.0_RFREAL                
      END DO ! ifl                      
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InitPatchCoeffs







! *******************************************************************************
!
! Purpose: Nullify patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_NullifyPatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_NullifyPatchCoeffs',&
  'RFLU_ModPatchCoeffs.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      NULLIFY(pPatch%cp)
      NULLIFY(pPatch%cf)
      NULLIFY(pPatch%ch)           
      NULLIFY(pPatch%cmass)           
      NULLIFY(pPatch%cmom)           
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyPatchCoeffs







! *******************************************************************************
!
! Purpose: Open patch-coefficient file in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: 
!   fileExists          Flag indicating existence of file
!
! Notes:
!   1. File may not exist when it is opened in postprocessor before solver 
!      was run, so need to deal with this gracefully.
!
! ******************************************************************************

  SUBROUTINE RFLU_OpenPatchCoeffsASCII(pRegion,fileStatus,fileExists)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: fileStatus
    LOGICAL, INTENT(OUT), OPTIONAL :: fileExists
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global 
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_OpenPatchCoeffsASCII',&
  'RFLU_ModPatchCoeffs.F90')
         
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Opening ASCII patch-coefficients file...'
    END IF ! global%verbLevel     
     
    iFile = IF_PATCH_COEF
     
! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.pcoa', & 
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
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.pcoa', & 
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)  

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN                                      
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter             
      END IF ! global%verbLevel
    ENDIF ! global%flowType

! ******************************************************************************
!   Open file
! ******************************************************************************

    IF ( fileStatus == FILE_STATUS_OLD ) THEN 
      INQUIRE(FILE=iFileName,EXIST=fileExists)
    
      IF ( fileExists .EQV. .TRUE. ) THEN     
        OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", &
             IOSTAT=errorFlag)
        global%error = errorFlag          
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
        END IF ! global%error
      END IF ! fileExists
    ELSE IF ( fileStatus == FILE_STATUS_UNKNOWN ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
           IOSTAT=errorFlag)
      global%error = errorFlag          
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
      END IF ! global%error      
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
    END IF ! fileStatus       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Opening ASCII patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_OpenPatchCoeffsASCII







! *******************************************************************************
!
! Purpose: Open patch-coefficient file in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: 
!   fileExists          Flag indicating existence of file
!
! Notes: 
!   1. File may not exist when it is opened in postprocessor before solver 
!      was run, so need to deal with this gracefully.
!
! ******************************************************************************

  SUBROUTINE RFLU_OpenPatchCoeffsBinary(pRegion,fileStatus,fileExists)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: fileStatus
    LOGICAL, INTENT(OUT), OPTIONAL :: fileExists
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global 
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_OpenPatchCoeffsBinary',&
  'RFLU_ModPatchCoeffs.F90')
         
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Opening binary patch-coefficients file...'
    END IF ! global%verbLevel     
     
    iFile = IF_PATCH_COEF     
     
! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.pco', & 
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
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.pco', & 
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)  

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN                                      
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter             
      END IF ! global%verbLevel
    ENDIF ! global%flowType

! ******************************************************************************
!   Open file
! ******************************************************************************

    IF ( fileStatus == FILE_STATUS_OLD ) THEN 
      INQUIRE(FILE=iFileName,EXIST=fileExists)
    
      IF ( fileExists .EQV. .TRUE. ) THEN     
        OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
             IOSTAT=errorFlag)
        global%error = errorFlag          
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
        END IF ! global%error
      END IF ! fileExists
    ELSE IF ( fileStatus == FILE_STATUS_UNKNOWN ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
           IOSTAT=errorFlag)
      global%error = errorFlag          
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
      END IF ! global%error      
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
    END IF ! fileStatus         

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Opening binary patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_OpenPatchCoeffsBinary







! *******************************************************************************
!
! Purpose: Read patch coefficients in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadPatchCoeffsASCII(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch,loopCounter
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ReadPatchCoeffsASCII',&
  'RFLU_ModPatchCoeffs.F90')
        
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading ASCII patch-coefficients file...'
    END IF ! global%verbLevel        
   
    pGrid => pRegion%grid
    
    iFile = IF_PATCH_COEF   
        
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU patch-coefficients file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read rest of file
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Pressure coefficient
! ==============================================================================

        CASE ( '# Pressure coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile,'(5(E23.16))') (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch
        
! ==============================================================================
!       Skin-friction coefficient
! ==============================================================================

        CASE ( '# Skin-friction coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile,'(5(E23.16))') (pPatch%cf(XCOORD,ifl), & 
                                       ifl=1,pPatch%nBFaces)       
            READ(iFile,'(5(E23.16))') (pPatch%cf(YCOORD,ifl), &
                                       ifl=1,pPatch%nBFaces)       
            READ(iFile,'(5(E23.16))') (pPatch%cf(ZCOORD,ifl), &
                                       ifl=1,pPatch%nBFaces) 
          END DO ! iPatch

! ==============================================================================
!       Heat-transfer coefficient
! ==============================================================================

        CASE ( '# Heat-transfer coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile,'(5(E23.16))') (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch

! ==============================================================================
!       End marker
! ==============================================================================
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

          EXIT
      
! ==============================================================================
!       Invalid section string
! ==============================================================================
      
        CASE DEFAULT
          IF ( global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel           

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)                               
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading ASCII patch-coefficients file done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadPatchCoeffsASCII






! *******************************************************************************
!
! Purpose: Read patch coefficients in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadPatchCoeffsBinary(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch,loopCounter
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ReadPatchCoeffsBinary',&
  'RFLU_ModPatchCoeffs.F90')
   
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reading binary patch-coefficients file...'
    END IF ! global%verbLevel
   
    pGrid => pRegion%grid   
     
    iFile = IF_PATCH_COEF
         
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU patch-coefficients file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read rest of file
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Pressure coefficient
! ==============================================================================

        CASE ( '# Pressure coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile) (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch
        
! ==============================================================================
!       Skin-friction coefficient
! ==============================================================================

        CASE ( '# Skin-friction coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile) (pPatch%cf(XCOORD,ifl),ifl=1,pPatch%nBFaces)       
            READ(iFile) (pPatch%cf(YCOORD,ifl),ifl=1,pPatch%nBFaces)       
            READ(iFile) (pPatch%cf(ZCOORD,ifl),ifl=1,pPatch%nBFaces) 
          END DO ! iPatch

! ==============================================================================
!       Heat-transfer coefficient
! ==============================================================================

        CASE ( '# Heat-transfer coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile) (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch

! ==============================================================================
!       End marker
! ==============================================================================
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

          EXIT
      
! ==============================================================================
!       Invalid section string
! ==============================================================================
      
        CASE DEFAULT
          IF ( global%verbLevel >= VERBOSE_MED) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel           

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)                               
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading binary patch-coefficients file done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadPatchCoeffsBinary





! *******************************************************************************
!
! Purpose: Wrapper for reading patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadPatchCoeffsWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    LOGICAL :: fileExists
    TYPE(t_global), POINTER :: global  
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ReadPatchCoeffsWrapper',&
  'RFLU_ModPatchCoeffs.F90')

! ******************************************************************************
!   Read file (if it exists)
! ******************************************************************************
    
    IF ( global%solutFormat == FORMAT_ASCII ) THEN
      CALL RFLU_OpenPatchCoeffsASCII(pRegion,FILE_STATUS_OLD,fileExists)
      
      IF ( fileExists .EQV. .TRUE. ) THEN
        CALL RFLU_ReadPatchCoeffsASCII(pRegion)
        CALL RFLU_ClosePatchCoeffs(pRegion)
      END IF ! fileExists
    ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
      CALL RFLU_OpenPatchCoeffsBinary(pRegion,FILE_STATUS_OLD,fileExists)
      
      IF ( fileExists .EQV. .TRUE. ) THEN 
        CALL RFLU_ReadPatchCoeffsBinary(pRegion)
        CALL RFLU_ClosePatchCoeffs(pRegion)        
      END IF ! fileExists
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%solutFormat   

! ******************************************************************************
!   Write warning if file is missing
! ******************************************************************************
    
    IF ( fileExists .EQV. .FALSE. ) THEN 
      global%warnCounter = global%warnCounter + 1
    
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,2(1X,A))') SOLVER_NAME,'*** WARNING ***', &
                                    'Patch coefficient file missing, not read.'        
      END IF ! global%myProcid
    END IF ! fileExists 
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadPatchCoeffsWrapper






! *******************************************************************************
!
! Purpose: Write patch coefficients in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WritePatchCoeffsASCII(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_WritePatchCoeffsASCII',&
  'RFLU_ModPatchCoeffs.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing ASCII patch-coefficients file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PATCH_COEF
            
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU patch-coefficients file'
    WRITE(iFile,'(A)') sectionString
              
! ******************************************************************************
!   Write patch coefficients to file
! ******************************************************************************
    
! ==============================================================================
!   Pressure coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Pressure coefficient'
    WRITE(iFile,'(A)') sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile,'(5(E23.16))') (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch
    
! ==============================================================================
!   Skin-friction coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Skin-friction coefficient'
    WRITE(iFile,'(A)') sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile,'(5(E23.16))') (pPatch%cf(XCOORD,ifl),ifl=1,pPatch%nBFaces)
      WRITE(iFile,'(5(E23.16))') (pPatch%cf(YCOORD,ifl),ifl=1,pPatch%nBFaces) 
      WRITE(iFile,'(5(E23.16))') (pPatch%cf(ZCOORD,ifl),ifl=1,pPatch%nBFaces)                    
    END DO ! iPatch    
    
! ==============================================================================
!   Heat-transfer coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Heat-transfer coefficient'
    WRITE(iFile,'(A)') sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile,'(5(E23.16))') (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch    

! ******************************************************************************
!   End marker
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') sectionString  

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Writing ASCII patch-coefficients file done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WritePatchCoeffsASCII






! *******************************************************************************
!
! Purpose: Write patch coefficients in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WritePatchCoeffsBinary(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_WritePatchCoeffsBinary',&
  'RFLU_ModPatchCoeffs.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing binary patch-coefficients file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PATCH_COEF
             
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU patch-coefficients file'
    WRITE(iFile) sectionString
              
! ******************************************************************************
!   Write patch coefficients to file
! ******************************************************************************
    
! ==============================================================================
!   Pressure coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Pressure coefficient'
    WRITE(iFile) sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile) (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch
    
! ==============================================================================
!   Skin-friction coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Skin-friction coefficient'
    WRITE(iFile) sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile) (pPatch%cf(XCOORD,ifl),ifl=1,pPatch%nBFaces)
      WRITE(iFile) (pPatch%cf(YCOORD,ifl),ifl=1,pPatch%nBFaces) 
      WRITE(iFile) (pPatch%cf(ZCOORD,ifl),ifl=1,pPatch%nBFaces)                    
    END DO ! iPatch    
    
! ==============================================================================
!   Heat-transfer coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Heat-transfer coefficient'
    WRITE(iFile) sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile) (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch    

! ******************************************************************************
!   End marker
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString  

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Writing binary patch-coefficients file done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WritePatchCoeffsBinary






! *******************************************************************************
!
! Purpose: Wrapper for writing patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WritePatchCoeffsWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    TYPE(t_global), POINTER :: global  
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_WritePatchCoeffsWrapper',&
  'RFLU_ModPatchCoeffs.F90')
    
    IF ( global%solutFormat == FORMAT_ASCII ) THEN
      CALL RFLU_OpenPatchCoeffsASCII(pRegion,FILE_STATUS_UNKNOWN)
      CALL RFLU_WritePatchCoeffsASCII(pRegion)
      CALL RFLU_ClosePatchCoeffs(pRegion)
    ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
      CALL RFLU_OpenPatchCoeffsBinary(pRegion,FILE_STATUS_UNKNOWN)
      CALL RFLU_WritePatchCoeffsBinary(pRegion)
      CALL RFLU_ClosePatchCoeffs(pRegion)      
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%solutFormat   
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WritePatchCoeffsWrapper






END MODULE RFLU_ModPatchCoeffs

! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModPatchCoeffs.F90,v $
!   Revision 1.6  2008/12/06 08:44:23  mtcampbe
!   Updated license.
!
!   Revision 1.5  2008/11/19 22:17:34  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.4  2006/10/20 21:18:56  mparmar
!   Added allocation/deallocation of cmass and cmom
!
!   Revision 1.3  2006/04/07 15:19:20  haselbac
!   Removed tabs
!
!   Revision 1.2  2004/12/07 19:39:14  haselbac
!   Bug fix: Incorrect WRITE statements when writing binary data
!
!   Revision 1.1  2004/06/16 20:00:59  haselbac
!   Initial revision
!
! ******************************************************************************



















