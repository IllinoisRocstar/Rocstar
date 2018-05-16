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
! Purpose: Suite of routines to compute, read, and write dimension files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModDimensions.F90,v 1.15 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModDimensions

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModBndPatch, ONLY: t_patch
  USE ModBorder, ONLY: t_border
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI
    
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_ReadDimensions, & 
            RFLU_ReadDimensionsWrapper, &
            RFLU_SetMaxDimension, & 
            RFLU_SetMaxDimensions, & 
            RFLU_WriteDimensions, &
            RFLU_WriteDimensionsBorders, &  
            RFLU_WriteDimensionsWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Public
! ==============================================================================

  INTEGER, PARAMETER, PUBLIC :: WRITE_DIMENS_MODE_FORCE = 0, &
                                WRITE_DIMENS_MODE_MAYBE = 1

! ==============================================================================
! Private
! ==============================================================================
   
  CHARACTER(CHRLEN), PRIVATE :: RCSIdentString = &
    '$RCSfile: RFLU_ModDimensions.F90,v $ $Revision: 1.15 $' 

  INTEGER, PARAMETER, PRIVATE :: ratioMax2Tot = 4
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Read dimensions.
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

  SUBROUTINE RFLU_ReadDimensions(pRegion)

    USE RFLU_ModGrid

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
    INTEGER :: errorFlag,dummy,iBorder,iPatch,iFile,loopCounter
    REAL(RFREAL) :: currentTime
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadDimensions',&
  'RFLU_ModDimensions.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading dimensions...'
    END IF ! global%verbLevel

    iFile = IF_DIMS
  
    IF ( (global%flowType == FLOW_UNSTEADY) .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN
         
! TEMPORARY - Read files from time zero because 1) for GENx do not write 
! dimension files and 2) because for standalone computations in which rflupost
! is to merge the files, the dimensions of the serial region at the given 
! time-stamp at which postprocessing occurs is not known because that file 
! does not exist. When need to merge files, need to write new routines which 
! will compute dimensions of serial region from parallel regions. This is not 
! as straightforward as may appear at first sight because vertices cannot be 
! summed... Do not actually need to build vertex lists, but would need to 
! estimate number of vertices. 
      IF ( global%timeStamp > 0.0_RFREAL ) THEN
        global%warnCounter = global%warnCounter + 1    

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel >= VERBOSE_NONE ) THEN        
          WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
                                        'Hard-code - read file from time zero.'
        END IF ! global%myProcid
      END IF ! global%timeStamp

      currentTime = 0.0_RFREAL 
! END TEMPORARY      
      
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.dim', & 
                                 pRegion%iRegionGlobal,currentTime, & 
                                 iFileName)
                                 
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            currentTime   
      END IF ! global%verbLevel
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.dim', & 
                              pRegion%iRegionGlobal,iFileName) 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", & 
         IOSTAT=errorFlag) 
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(iFileName))
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    READ(iFile,'(A)') sectionString  
    IF ( TRIM(sectionString) /= '# ROCFLU dimensions file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM  
    
! ==============================================================================
!   Rest of file
! ==============================================================================
  
    pGrid => pRegion%grid  

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!       Vertices
! ------------------------------------------------------------------------------ 

        CASE ( '# Vertices' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Vertices...'
          END IF ! global%verbLevel          

          READ(iFile,'(3(I8))') pGrid%nVert,pGrid%nVertTot,pGrid%nVertMax        

! ------------------------------------------------------------------------------
!       Cells
! ------------------------------------------------------------------------------ 
    
        CASE ( '# Cells' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
          END IF ! global%verbLevel          

          READ(iFile,'(3(I8))') pGrid%nCells,pGrid%nCellsTot,pGrid%nCellsMax         
        
! ------------------------------------------------------------------------------
!       Tetrahedra
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Tetrahedra' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
          END IF ! global%verbLevel          

          READ(iFile,'(3(I8))') pGrid%nTets,pGrid%nTetsTot,pGrid%nTetsMax         
    
! ------------------------------------------------------------------------------
!       Hexahedra
! ------------------------------------------------------------------------------ 

        CASE ( '# Hexahedra' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...'
          END IF ! global%verbLevel          

          READ(iFile,'(3(I8))') pGrid%nHexs,pGrid%nHexsTot,pGrid%nHexsMax        
    
! ------------------------------------------------------------------------------
!       Prisms
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Prisms' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Prisms...'
          END IF ! global%verbLevel          

          READ(iFile,'(3(I8))') pGrid%nPris,pGrid%nPrisTot,pGrid%nPrisMax      
    
! ------------------------------------------------------------------------------
!       Pyramids
! ------------------------------------------------------------------------------ 

        CASE ( '# Pyramids' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pyramids...'
          END IF ! global%verbLevel          

          READ(iFile,'(3(I8))') pGrid%nPyrs,pGrid%nPyrsTot,pGrid%nPyrsMax       
    
! ------------------------------------------------------------------------------
!       Patches (old style, retained for backward compatibility). NOTE 
!       initialize new variables here so as to make sure that all new data 
!       comes out of this routine with sensible values.
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Patches' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patches...'
          END IF ! global%verbLevel          

          READ(iFile,'(2(I8))') pGrid%nPatches,global%nPatches        

          IF ( pGrid%nPatches > PATCH_DIMENS_NPATCHMAX ) THEN 
            CALL ErrorStop(global,ERR_PATCH_DIMENS,__LINE__)
          END IF ! pGrid%nPatches

          DO iPatch = 1,pGrid%nPatches 
            READ(iFile,'(5(I8))',IOSTAT=errorFlag) & 
              pGrid%patchDimens(PATCH_DIMENS_IPGLOBAL   ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBTRIS     ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT  ,iPatch), &
              pGrid%patchDimens(PATCH_DIMENS_NBQUADS    ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT ,iPatch)
              
            pGrid%patchDimens(PATCH_DIMENS_NBTRISMAX,iPatch) = & 
              pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT,iPatch)
            pGrid%patchDimens(PATCH_DIMENS_NBQUADSMAX,iPatch) = &
              pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT,iPatch)
            pGrid%patchDimens(PATCH_DIMENS_NBCELLSVIRT,iPatch) = 0           
          END DO ! iPatch
          
! ------------------------------------------------------------------------------
!       Patches (new style)
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Patches (v2)' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patches...'
          END IF ! global%verbLevel          

          READ(iFile,'(2(I8))') pGrid%nPatches,global%nPatches        

          IF ( pGrid%nPatches > PATCH_DIMENS_NPATCHMAX ) THEN 
            CALL ErrorStop(global,ERR_PATCH_DIMENS,__LINE__)
          END IF ! pGrid%nPatches

          DO iPatch = 1,pGrid%nPatches 
            READ(iFile,'(8(I8))',IOSTAT=errorFlag) & 
              pGrid%patchDimens(PATCH_DIMENS_IPGLOBAL   ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBTRIS     ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT  ,iPatch), &
              pGrid%patchDimens(PATCH_DIMENS_NBTRISMAX  ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBQUADS    ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT ,iPatch), &
              pGrid%patchDimens(PATCH_DIMENS_NBQUADSMAX ,iPatch), & 
              pGrid%patchDimens(PATCH_DIMENS_NBCELLSVIRT,iPatch)
          END DO ! iPatch          
    
! ------------------------------------------------------------------------------
!       Borders
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Borders' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Borders...'
          END IF ! global%verbLevel          

          READ(iFile,'(I8)') pGrid%nBorders        

          DO iBorder = 1,pGrid%nBorders   
            READ(iFile,'(7(I8))') & 
              pGrid%borderInfo(BORDER_INFO_IRGLOB,iBorder), & 
              pGrid%borderInfo(BORDER_INFO_IBORD ,iBorder), & 
              pGrid%borderInfo(BORDER_INFO_NCSEND,iBorder), & 
              pGrid%borderInfo(BORDER_INFO_NCRECV,iBorder), & 
              pGrid%borderInfo(BORDER_INFO_NVSEND,iBorder), & 
              pGrid%borderInfo(BORDER_INFO_NVRECV,iBorder), & 
              pGrid%borderInfo(BORDER_INFO_NVSHAR,iBorder)
          END DO ! iBorder   
    
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 

        CASE DEFAULT
          IF ( global%verbLevel >= VERBOSE_HIGH ) THEN  
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
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading dimensions done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_ReadDimensions







! ******************************************************************************
!
! Purpose: Wrapper function for reading dimensions.
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

  SUBROUTINE RFLU_ReadDimensionsWrapper(pRegion)

#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_GetDimensions
#endif

#ifdef PLAG
    USE PLAG_ModDimensions, ONLY: PLAG_RFLU_ReadDimensions
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

    CALL RegisterFunction(global,'RFLU_ReadDimensionsWrapper',&
  'RFLU_ModDimensions.F90')

! ******************************************************************************
!   Call routines
! ******************************************************************************

! TEMPORARY
!#ifndef GENX
    CALL RFLU_ReadDimensions(pRegion)
!#else
!    CALL RFLU_GENX_GetDimensions(pRegion)
!#endif
! END TEMPORARY

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      CALL PLAG_RFLU_ReadDimensions(pRegion)
    END IF ! plagUsed
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadDimensionsWrapper







! ******************************************************************************
!
! Purpose: Set maximum dimension.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   nXyzTot     Some total dimension
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  INTEGER FUNCTION RFLU_SetMaxDimension(global,nXyzTot)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nXyzTot
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    IF ( (global%moduleType == MODULE_TYPE_PART) .AND. & 
         (global%syPePatchesFlag .EQV. .TRUE.) ) THEN 
      RFLU_SetMaxDimension = INT(4*ratioMax2Tot*nXyzTot)    
    ELSE 
      RFLU_SetMaxDimension = ratioMax2Tot*nXyzTot        
    END IF ! global%moduleType

! ******************************************************************************
!   End
! ******************************************************************************
  
  END FUNCTION RFLU_SetMaxDimension








! ******************************************************************************
!
! Purpose: Set maximum dimensions.
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

  SUBROUTINE RFLU_SetMaxDimensions(pRegion)
  
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

    INTEGER :: iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetMaxDimensions',&
  'RFLU_ModDimensions.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting maximum dimensions...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************
  
    pGrid => pRegion%grid  

! ******************************************************************************
!   Print information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Ratio:',ratioMax2Tot
    END IF ! global%verbLevel

! ******************************************************************************
!   Set maximum dimensions
! ******************************************************************************
      
    pGrid%nVertMax = ratioMax2Tot*pGrid%nVertTot    
   
    pGrid%nTetsMax = ratioMax2Tot*pGrid%nTetsTot
    pGrid%nHexsMax = ratioMax2Tot*pGrid%nHexsTot
    pGrid%nPrisMax = ratioMax2Tot*pGrid%nPrisTot
    pGrid%nPyrsMax = ratioMax2Tot*pGrid%nPyrsTot

    pGrid%nCellsMax = pGrid%nTetsMax & 
                    + pGrid%nHexsMax & 
                    + pGrid%nPrisMax & 
                    + pGrid%nPyrsMax
                    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      pPatch%nBTrisMax  = ratioMax2Tot*pPatch%nBTrisTot
      pPatch%nBQuadsMax = ratioMax2Tot*pPatch%nBQuadsTot    
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting maximum dimensions done.'
    END IF ! global%verbLevel
    
    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_SetMaxDimensions





! ******************************************************************************
!
! Purpose: Write dimensions.
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

  SUBROUTINE RFLU_WriteDimensions(pRegion)

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
    INTEGER :: errorFlag,iBorder,iPatch,iFile
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteDimensions',&
  'RFLU_ModDimensions.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing dimensions...'
    END IF ! global%verbLevel

    iFile = IF_DIMS

    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN          
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.dim', & 
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
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.dim', & 
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
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    sectionString = '# ROCFLU dimensions file'  
    WRITE(iFile,'(A)') TRIM(sectionString)  
      
! ==============================================================================
!   Write dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    sectionString = '# Vertices'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(3(I8))') pGrid%nVert,pGrid%nVertTot,pGrid%nVertMax

    sectionString = '# Cells'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(3(I8))') pGrid%nCells,pGrid%nCellsTot,pGrid%nCellsMax

    sectionString = '# Tetrahedra'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(3(I8))') pGrid%nTets,pGrid%nTetsTot,pGrid%nTetsMax

    sectionString = '# Hexahedra'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(3(I8))') pGrid%nHexs,pGrid%nHexsTot,pGrid%nHexsMax    

    sectionString = '# Prisms'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(3(I8))') pGrid%nPris,pGrid%nPrisTot,pGrid%nPrisMax  

    sectionString = '# Pyramids'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(3(I8))') pGrid%nPyrs,pGrid%nPyrsTot,pGrid%nPyrsMax  

    sectionString = '# Patches (v2)'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(2(I8))') pGrid%nPatches,global%nPatches
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  
 
      WRITE(iFile,'(8(I8))') pPatch%iPatchGlobal, &
                             pPatch%nBTris,pPatch%nBTrisTot, &
                             pPatch%nBTrisMax,pPatch%nBQuads, & 
                             pPatch%nBQuadsTot,pPatch%nBQuadsMax, & 
                             pPatch%nBCellsVirt
    END DO ! iPatch


    sectionString = '# Borders'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pGrid%nBorders
    
    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)  

      WRITE(iFile,'(7(I8))') pBorder%iRegionGlobal,pBorder%iBorder, & 
                             pBorder%nCellsSend,pBorder%nCellsRecv, &
                             pBorder%nVertSend,pBorder%nVertRecv, &
                             pBorder%nVertShared
    END DO ! iBorder    

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN  
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
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing dimensions done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_WriteDimensions
  








! ******************************************************************************
!
! Purpose: Write border dimensions.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine must only be called from the preprocessor to write the 
!      border information. 
!   2. This routine reads the dimensions file until the borders section is 
!      found, then the borders information is written. 
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteDimensionsBorders(pRegion)

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

    CHARACTER(CHRLEN) :: dummyString,iFileName,sectionString
    INTEGER :: errorFlag,iBorder,iPatch,iFile,loopCounter
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteDimensionsBorders',&
  'RFLU_ModDimensions.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing border dimensions...'
    END IF ! global%verbLevel

    iFile = IF_DIMS

    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN          
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.dim', & 
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
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.dim', & 
                              pRegion%iRegionGlobal,iFileName) 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal                                 
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error
      
! ==============================================================================
!   Read dimensions file line by line until hit border section, then write it.
! ==============================================================================
  
    pGrid => pRegion%grid  

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') dummyString

      SELECT CASE ( TRIM(dummyString) )        

! ------------------------------------------------------------------------------
!       Borders section. NOTE once found borders section, write information, 
!       end string, and exit infinite loop.
! ------------------------------------------------------------------------------ 

        CASE ( '# Borders' )     
          WRITE(iFile,'(I8)') pGrid%nBorders        

          DO iBorder = 1,pGrid%nBorders   
            pBorder => pGrid%borders(iBorder)
            
            WRITE(iFile,'(7(I8))') pBorder%iRegionGlobal,pBorder%iBorder, & 
                                   pBorder%nCellsSend,pBorder%nCellsRecv, &
                                   pBorder%nVertSend,pBorder%nVertRecv, &
                                   pBorder%nVertShared
          END DO ! iBorder  
          
          sectionString = '# End'
          WRITE(iFile,'(A)') TRIM(sectionString) 
          
          EXIT 
    
! ------------------------------------------------------------------------------
!       End marker. NOTE this must not happen in present context, so trap as 
!       an error.
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_HIGH ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

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
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing border dimensions done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_WriteDimensionsBorders








! ******************************************************************************
!
! Purpose: Wrapper function for writing dimensions.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   writeMode   Writing mode
!
! Output: None.
!
! Notes: 
!   1. Need to call this subroutine even if grid is not moving because 
!      Rocpart dimension file needs to be written. 
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteDimensionsWrapper(pRegion,writeMode)

#ifdef PLAG
    USE PLAG_ModDimensions, ONLY: PLAG_RFLU_WriteDimensions
#endif

#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_DecideWriteFile
#endif

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: writeMode
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteDimensionsWrapper',&
  'RFLU_ModDimensions.F90')

! ******************************************************************************
!   Call routines
! ******************************************************************************

#ifdef GENX
    IF ( RFLU_GENX_DecideWriteFile(global) .EQV. .TRUE. ) THEN     
#endif
      IF ( writeMode == WRITE_DIMENS_MODE_FORCE ) THEN 
        CALL RFLU_WriteDimensions(pRegion)
      ELSE 
        IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN
          CALL RFLU_WriteDimensions(pRegion)
        END IF ! pRegion%mixtInput%moveGrid 
      END IF ! writeMode
#ifdef GENX
    END IF ! RFLU_GENX_DecideWriteFile
#endif

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      CALL PLAG_RFLU_WriteDimensions(pRegion)
    END IF ! plagUsed
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteDimensionsWrapper








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModDimensions


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModDimensions.F90,v $
! Revision 1.15  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2007/03/20 17:34:22  fnajjar
! Modified USE calls to streamline with new module PLAG_ModDimensions
!
! Revision 1.12  2006/12/15 13:21:18  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.11  2006/10/08 18:21:31  haselbac
! Changes to make code work with PathScale compilers on hpc cluster at UF
!
! Revision 1.10  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.9  2006/03/25 23:52:15  haselbac
! Restored backward compatibility thru new keyword for patches
!
! Revision 1.8  2006/03/25 21:51:47  haselbac
! Substantial changes because of sype patches
!
! Revision 1.7  2005/06/11 20:52:27  haselbac
! Fixed comment and note
!
! Revision 1.6  2005/06/11 20:33:55  haselbac
! Removed ifdef GENX because now always read from 0 if grid moves
!
! Revision 1.5  2005/04/15 15:06:49  haselbac
! Added vertex dimensions to border section
!
! Revision 1.4  2004/12/29 21:06:30  haselbac
! Modified reading of border dims to mirror patch dims
!
! Revision 1.3  2004/12/04 03:29:22  haselbac
! Added reading/writing of border dims, routine to write border dims
!
! Revision 1.2  2004/10/19 19:27:49  haselbac
! Added/removed procs, adapted to changes in other procs
!
! Revision 1.1  2004/07/06 15:14:26  haselbac
! Initial revision
!
! ******************************************************************************












