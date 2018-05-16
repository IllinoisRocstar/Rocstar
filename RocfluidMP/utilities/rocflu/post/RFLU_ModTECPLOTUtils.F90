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

!
! Purpose: Collection of utility routines for writing TECPLOT files.
!
! Description: None.
!
! Input: None.
! 
! Output: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModTECPLOTUtils.F90,v 1.14 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModTECPLOTUtils

  USE ModDataTypes
  USE ModParameters
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

  PRIVATE

! ==============================================================================
! Private data
! ==============================================================================
 
  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModTECPLOTUtils.F90,v $ $Revision: 1.14 $'

! ==============================================================================
! Public data
! ==============================================================================

  CHARACTER(1), PARAMETER, PUBLIC :: NULLCHAR = CHAR(0)
  CHARACTER(2*CHRLEN), PUBLIC :: headerTEC
  INTEGER, PARAMETER, PUBLIC :: ZONE_TYPE_TRI  = 2, &
                                ZONE_TYPE_QUAD = 3, & 
                                ZONE_TYPE_TET  = 4, & 
                                ZONE_TYPE_HEX  = 5 
  INTEGER, PARAMETER, PUBLIC :: ZONE_FORM_BLOCK = 1
  INTEGER, PARAMETER, PUBLIC :: VAR_POS_CELL = 0, & 
                                VAR_POS_FACE = 0, & 
                                VAR_POS_VERT = 1
  INTEGER, PARAMETER, PUBLIC :: FILE_TYPE_FIELD       = 1, & 
                                FILE_TYPE_PATCH       = 2, & 
                                FILE_TYPE_PATCH_STATS = 3
  INTEGER, PARAMETER, PUBLIC :: FILE_CNTR_TEC_MAX = 10
  INTEGER, PUBLIC :: debugFlagTEC,doubleFlagTEC,fileCntrTEC,nVarsTEC, &
                     nVarsTotTEC,nVarsCellTEC,nVarsFaceTEC,nVarsVertTEC
  INTEGER, PUBLIC :: fileType2CntrTEC(FILE_CNTR_TEC_MAX)
  INTEGER, DIMENSION(8), PARAMETER, PUBLIC :: tet2vTEC = (/1,2,3,3,4,4,4,4/)
  INTEGER, DIMENSION(8), PARAMETER, PUBLIC :: hex2vTEC = (/1,2,3,4,5,6,7,8/)
  INTEGER, DIMENSION(8), PARAMETER, PUBLIC :: pri2vTEC = (/1,2,3,1,4,5,6,4/)
  INTEGER, DIMENSION(8), PARAMETER, PUBLIC :: pyr2vTEC = (/1,2,3,4,5,5,5,5/)
  INTEGER, ALLOCATABLE, DIMENSION(:), PUBLIC :: posTEC
  REAL(RFREAL), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: varCellTEC,varFaceTEC, &
                                                       varVertTEC

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_TEC_CloseFile, & 
            RFLU_TEC_OpenFile, & 
            RFLU_TEC_WriteZoneCellsSpecial, &
            RFLU_TEC_WriteZoneFacesSpecial, & 
            RFLU_TEC_WriteZoneInterf, &
            RFLU_TEC_WriteZoneSurf,& 
            RFLU_TEC_WriteZoneSurfMixed, &
            RFLU_TEC_WriteZoneVol


! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS 
 
 
 
 
 
 
! ******************************************************************************
!
! Purpose: Close TECPLOT file
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
! 
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!
! ******************************************************************************  

SUBROUTINE RFLU_TEC_CloseFile(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECEND100

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TEC_CloseFile', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ******************************************************************************
! Close TECPLOT file
! ******************************************************************************

  errorFlag = TECEND100()
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_CloseFile
  
  





! ******************************************************************************
!
! Purpose: Open TECPLOT file
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   title               Title of TECPLOT dataset
!   fileName            Name of TECPLOT file
! 
! Output: None.
!
! Notes:
!   1. Isolated into separate routine so that can write several zones to same
!      TECPLOT file.
!
! ******************************************************************************  

SUBROUTINE RFLU_TEC_OpenFile(global,title,fileName)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(CHRLEN) :: fileName,title
  TYPE(t_global), POINTER :: global
 
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: scratchDir
  INTEGER :: errorFlag

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECINI100

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TEC_OpenFile', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ==============================================================================
! Set some TECPLOT variables
! ==============================================================================

  doubleFlagTEC = 1 ! Flag for single/double precision data
  debugFlagTEC  = 0 ! Flag for debugging
     
  scratchDir = '.'   
     
! ==============================================================================
! Open TECPLOT file
! ==============================================================================

  errorFlag = TECINI100(TRIM(title)//NULLCHAR, &
                        TRIM(headerTEC)//NULLCHAR, &
                        TRIM(fileName)//NULLCHAR, &
                        TRIM(scratchDir)//NULLCHAR, & 
                        debugFlagTEC, &
                        doubleFlagTEC)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error
   
! ******************************************************************************
! End
! ******************************************************************************  
    
  CALL DeregisterFunction(global)  
  
END SUBROUTINE RFLU_TEC_OpenFile







! ******************************************************************************
!
! Purpose: Write zone for special cells.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. Special cells are written as degenerate hexahedra.
!
! ******************************************************************************  

SUBROUTINE RFLU_TEC_WriteZoneCellsSpecial(pRegion)

  USE ModSortSearch
  
  USE RFLU_ModGrid

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
 
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: zonePrefix,zoneTitle 
  INTEGER :: cs2vCard,errorFlag,icg,icl,ics,ict,iLoc,iVar,iVarCell,iVarVert, &
             ivCntr,ivl,nVert,zoneType
  INTEGER :: cs2v(8),cs2vCopy(8),cs2vFlag(8),cs2vTEC(8) 
  INTEGER, DIMENSION(:), ALLOCATABLE :: nullArray
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: varCellTmp  
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: varVertTmp
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECDAT100,TECNOD100,TECZNE100
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteZoneCellsSpecial', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Create null array
! ******************************************************************************

  ALLOCATE(nullArray(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

  DO iVar = 1,nVarsTEC
    nullArray(iVar) = 0
  END DO ! iVar 
 
! ******************************************************************************
! Loop over special cells
! ******************************************************************************
  
  DO ics = 1,pGrid%nCellsSpecial
    icg = pGrid%cellsSpecial(ics)
    
! ==============================================================================
!   Determine cell type and local cell index
! ==============================================================================
        
    ict = pGrid%cellGlob2Loc(1,icg)
    icl = pGrid%cellGlob2Loc(2,icg)
        
    IF ( global%verbLevel > VERBOSE_NONE ) THEN     
      WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Writing special cell:',ics
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Global cell index:',icg
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Local cell index: ',icl  
    END IF ! global%verbLevel         
        
! ==============================================================================
!   Set local cell connectivity. Note the local cell connectivity for TECPLOT 
!   contains duplicated vertices...
! ==============================================================================
    
    SELECT CASE ( ict ) 
      CASE ( CELL_TYPE_TET )
        IF ( global%verbLevel > VERBOSE_NONE ) THEN     
          WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Local cell type: Tetrahedron'        
        END IF ! global%verbLevel     
             
        nVert = 4
             
        DO ivl = 1,8
          cs2v(ivl)    = pGrid%tet2v(tet2vTEC(ivl),icl)
          cs2vTEC(ivl) = tet2vTEC(ivl)        
        END DO ! ivl            
      CASE ( CELL_TYPE_HEX )      
        IF ( global%verbLevel > VERBOSE_NONE ) THEN     
          WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Local cell type: Hexahedron'        
        END IF ! global%verbLevel
                        
        nVert = 8
             
        DO ivl = 1,8
          cs2v(ivl)    = pGrid%hex2v(hex2vTEC(ivl),icl)
          cs2vTEC(ivl) = hex2vTEC(ivl)                   
        END DO ! ivl                   
      CASE ( CELL_TYPE_PRI ) 
        IF ( global%verbLevel > VERBOSE_NONE ) THEN     
          WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Local cell type: Prism'        
        END IF ! global%verbLevel      
      
        nVert = 6
             
        DO ivl = 1,8
          cs2v(ivl)    = pGrid%pri2v(pri2vTEC(ivl),icl)
          cs2vTEC(ivl) = pri2vTEC(ivl)                   
        END DO ! ivl            
      CASE ( CELL_TYPE_PYR )
        IF ( global%verbLevel > VERBOSE_NONE ) THEN     
          WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Local cell type: Pyramid'        
        END IF ! global%verbLevel      
       
        nVert = 5
             
        DO ivl = 1,8
          cs2v(ivl)    = pGrid%pyr2v(pyr2vTEC(ivl),icl)
          cs2vTEC(ivl) = pyr2vTEC(ivl)                   
        END DO ! ivl  
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ict
    
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
       WRITE(STDOUT,'(A,7X,A,8(1X,I9))') SOLVER_NAME,'Vertices:',cs2v(1:8) 
    END IF ! global%verbLevel
    
! ==============================================================================
!   Allocate temporary memory
! ==============================================================================
    
    ALLOCATE(varVertTmp(nVert,nVarsVertTEC),STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varVertTmp')
    END IF ! global%error  
    
    IF ( nVarsCellTEC > 0 ) THEN 
      ALLOCATE(varCellTmp(nVarsCellTEC),STAT=errorFlag)
      global%error = errorFlag      
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varCellTmp')
      END IF ! global%error      
    END IF ! nVarsCellTEC
    
! ==============================================================================
!   Build additional data structures: Sort cs2v array and delete duplicated 
!   vertices. This is required because solution is not written in duplicated 
!   manner. Also check that vertex list without duplicated vertices has correct 
!   length (same as nVertPerCell defined above).
! ==============================================================================
    
    cs2vCopy(1:8) = cs2v(1:8)
    
    CALL QuickSortInteger(cs2vCopy,8)
    CALL SimplifySortedIntegers(cs2vCopy,8,cs2vCard)
    
    IF ( cs2vCard /= nVert ) THEN 
      CALL ErrorStop(global,ERR_SVERT_LIST_INVALID,__LINE__)
    END IF ! cs2vCard
        
! ==============================================================================
!   Extract vertex variables
! ==============================================================================
             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN             
      WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Cell vertex coordinates:' 
      WRITE(STDOUT,'(A,9X,A,1X,A,4X,A,2X,6X,A,2(12X,A))') SOLVER_NAME,'#', &
            'Local','Global','x-coordinate','y-coordinate','z-coordinate'       
    END IF ! global%verbLevel             
             
    cs2vFlag(1:8) = 0
    ivCntr = 0
    
    DO ivl = 1,8                    
      CALL BinarySearchInteger(cs2vCopy(1:nVert),nVert,cs2v(ivl),iLoc)
                                  
      IF ( cs2vFlag(iLoc) == 0 ) THEN
        ivCntr = ivCntr + 1   
                    
        cs2vFlag(iLoc) = ivCntr        
        
        iVarVert = 0
        
        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
            iVarVert = iVarVert + 1
            varVertTmp(ivCntr,iVarVert) = pRegion%varVertTEC(cs2v(ivl),iVarVert)
          END IF ! posTEC        
        END DO ! iVar        
        
        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,9X,I1,3X,I1,3X,I9,1X,3(1X,E23.16))') SOLVER_NAME, &
                ivCntr,ivl,cs2v(ivl),pRegion%varVertTEC(cs2v(ivl),1:3)
        END IF ! global%verbLevel
      END IF ! cs2vFlag
    END DO ! ivl
    
! ==============================================================================
!   Extract cell variables
! ==============================================================================
    
    IF ( nVarsCellTEC > 0 ) THEN 
      iVarCell = 0

      DO iVar = 1,nVarsTEC
        IF ( posTEC(iVar) == VAR_POS_CELL ) THEN 
          iVarCell = iVarCell + 1
          varCellTmp(iVarCell) = pRegion%varCellTEC(icg,iVarCell)
        END IF ! posTEC        
      END DO ! iVar 
    END IF ! nVarsCellTEC    
    
! ==============================================================================
!   Write to TECPLOT file
! ==============================================================================
    
! ------------------------------------------------------------------------------
!   Zone information
! ------------------------------------------------------------------------------
    
    WRITE(zoneTitle,'(A,I8.8,A,I5.5)') 'C_',icg,'_',pRegion%iRegionGlobal      
    errorFlag = TECZNE100(TRIM(zoneTitle)//NULLCHAR,ZONE_TYPE_HEX,nVert,1, & 
                          0,0,0,0,ZONE_FORM_BLOCK,0,0,posTEC,nullArray,0)    
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error                          

! ------------------------------------------------------------------------------
!   Variables
! ------------------------------------------------------------------------------

    iVarCell = 0
    iVarVert = 0

    DO iVar = 1,nVarsTEC
      IF ( posTEC(iVar) == VAR_POS_CELL ) THEN 
        iVarCell = iVarCell + 1    
        errorFlag = TECDAT100(1,varCellTmp(iVarCell),doubleFlagTEC)     
      ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
        iVarVert = iVarVert + 1    
        errorFlag = TECDAT100(nVert,varVertTmp(1,iVarVert),doubleFlagTEC)      
      ELSE ! defensive coding
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! posTEC

      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
      END IF ! global%error       
    END DO ! iVar
        
! ------------------------------------------------------------------------------
!   Connectivity
! ------------------------------------------------------------------------------

    errorFlag = TECNOD100(cs2vTEC)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error         
        
! ==============================================================================
!   Deallocate temporary memory
! ==============================================================================
        
    DEALLOCATE(varVertTmp,STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varVertTmp')
    END IF ! global%error  
    
    IF ( nVarsCellTEC > 0 ) THEN 
      DEALLOCATE(varCellTmp,STAT=errorFlag)
      global%error = errorFlag      
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varCellTmp')
      END IF ! global%error      
    END IF ! nVarsCellTEC  
  END DO ! ics

! ******************************************************************************
! Destroy null array 
! ******************************************************************************

  DEALLOCATE(nullArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteZoneCellsSpecial









! ******************************************************************************
!
! Purpose: Write zone for special faces.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. Special faces are written as degenerate quadrilaterals.
!   2. Cannot write special faces when have cell-centered data because do not
!      have any data for interior faces. 
!
! ******************************************************************************  

SUBROUTINE RFLU_TEC_WriteZoneFacesSpecial(pRegion)

  USE ModSortSearch
  
  USE RFLU_ModGrid

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
 
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: zonePrefix,zoneTitle 
  INTEGER :: fs2vCard,errorFlag,ifl,ifs,iLoc,iPatch,iVar,iVarCell,iVarVert, &
             ivCntr,ivl,nVert,zoneType
  INTEGER :: fs2v(4),fs2vCopy(4),fs2vFlag(4),fs2vTEC(4) 
  INTEGER, DIMENSION(:), ALLOCATABLE :: nullArray 
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: varVertTmp
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECDAT100,TECNOD100,TECZNE100
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteZoneFacesSpecial', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Create null array
! ******************************************************************************

  ALLOCATE(nullArray(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

  DO iVar = 1,nVarsTEC
    nullArray(iVar) = 0
  END DO ! iVar 
 
! ******************************************************************************
! Loop over special cells
! ******************************************************************************
  
  DO ifs = 1,pGrid%nFacesSpecial
    iPatch = pGrid%facesSpecial(1,ifs)
    ifl    = pGrid%facesSpecial(2,ifs)
            
    IF ( global%verbLevel > VERBOSE_NONE ) THEN     
      WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Writing special face:',ifs
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Patch index:',iPatch
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Face index: ',ifl  
    END IF ! global%verbLevel         
        
! ==============================================================================
!   Set local face connectivity. Note the local face connectivity for TECPLOT 
!   contains duplicated vertices for triangles
! ==============================================================================
    
    IF ( iPatch == 0 ) THEN 
      DO ivl = 1,4
        fs2v(ivl)    = pGrid%f2v(ivl,ifl)
        fs2vTEC(ivl) = ivl
      END DO ! ivl
    ELSE IF ( iPatch > 0 .AND. iPatch <= pGrid%nPatches ) THEN 
      pPatch => pRegion%patches(iPatch)
      
      DO ivl = 1,4            
        IF ( pPatch%bf2v(ivl,ifl) /= VERT_NONE ) THEN 
          fs2v(ivl) = pPatch%bv(pPatch%bf2v(ivl,ifl))
        ELSE 
          fs2v(ivl) = VERT_NONE 
        END IF ! pPatch%bf2v
        
        fs2vTEC(ivl) = ivl              
      END DO ! ivl      
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! iPatch
    
    IF ( fs2v(4) == VERT_NONE ) THEN ! Triangle
      nVert = 3
    
      fs2v(4)    = fs2v(1)
      fs2vTEC(4) = 1      
    ELSE ! Quadrilateral
      nVert = 4      
    END IF ! fs2v
        
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
       WRITE(STDOUT,'(A,7X,A,4(1X,I9))') SOLVER_NAME,'Vertices:',fs2v(1:4) 
    END IF ! global%verbLevel
    
! ==============================================================================
!   Allocate temporary memory
! ==============================================================================
    
    ALLOCATE(varVertTmp(nVert,nVarsVertTEC),STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varVertTmp')
    END IF ! global%error  
        
! ==============================================================================
!   Build additional data structures: Sort cs2v array and delete duplicated 
!   vertices. This is required because solution is not written in duplicated 
!   manner. Also check that vertex list without duplicated vertices has correct 
!   length (same as nVertPerCell defined above).
! ==============================================================================
    
    fs2vCopy(1:4) = fs2v(1:4)
    
    CALL QuickSortInteger(fs2vCopy,4)
    CALL SimplifySortedIntegers(fs2vCopy,4,fs2vCard)
    
    IF ( fs2vCard /= nVert ) THEN 
      CALL ErrorStop(global,ERR_SVERT_LIST_INVALID,__LINE__)
    END IF ! fs2vCard
        
! ==============================================================================
!   Extract vertex variables
! ==============================================================================
             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN             
      WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Cell vertex coordinates:' 
      WRITE(STDOUT,'(A,9X,A,1X,A,4X,A,2X,6X,A,2(12X,A))') SOLVER_NAME,'#', &
            'Local','Global','x-coordinate','y-coordinate','z-coordinate'       
    END IF ! global%verbLevel             
             
    fs2vFlag(1:4) = 0
    ivCntr = 0
    
    DO ivl = 1,4                    
      CALL BinarySearchInteger(fs2vCopy(1:nVert),nVert,fs2v(ivl),iLoc)
                                  
      IF ( fs2vFlag(iLoc) == 0 ) THEN
        ivCntr = ivCntr + 1   
                    
        fs2vFlag(iLoc) = ivCntr        
        
        iVarVert = 0
        
        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
            iVarVert = iVarVert + 1
            varVertTmp(ivCntr,iVarVert) = pRegion%varVertTEC(fs2v(ivl),iVarVert)
          END IF ! posTEC        
        END DO ! iVar        
        
        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,9X,I1,3X,I1,3X,I9,1X,3(1X,E23.16))') SOLVER_NAME, &
                ivCntr,ivl,fs2v(ivl),pRegion%varVertTEC(fs2v(ivl),1:3)
        END IF ! global%verbLevel
      END IF ! cs2vFlag
    END DO ! ivl
         
! ==============================================================================
!   Write to TECPLOT file
! ==============================================================================
    
! ------------------------------------------------------------------------------
!   Zone information
! ------------------------------------------------------------------------------
    
    WRITE(zoneTitle,'(A,I8.8,A,I3.3,A,I5.5)') 'F_',ifl,'_PAT_',iPatch,'_', & 
                                              pRegion%iRegionGlobal      
    errorFlag = TECZNE100(TRIM(zoneTitle)//NULLCHAR,ZONE_TYPE_QUAD,nVert,1, & 
                          0,0,0,0,ZONE_FORM_BLOCK,0,0,posTEC,nullArray,0)    
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error                          

! ------------------------------------------------------------------------------
!   Variables
! ------------------------------------------------------------------------------

    iVarVert = 0

    DO iVar = 1,nVarsTEC
      IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
        iVarVert = iVarVert + 1    
        errorFlag = TECDAT100(nVert,varVertTmp(1,iVarVert),doubleFlagTEC)      
      ELSE ! defensive coding
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! posTEC

      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
      END IF ! global%error       
    END DO ! iVar
        
! ------------------------------------------------------------------------------
!   Connectivity
! ------------------------------------------------------------------------------

    errorFlag = TECNOD100(fs2vTEC)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error         
        
! ==============================================================================
!   Deallocate temporary memory
! ==============================================================================
        
    DEALLOCATE(varVertTmp,STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varVertTmp')
    END IF ! global%error  
  END DO ! ifs

! ******************************************************************************
! Destroy null array 
! ******************************************************************************

  DEALLOCATE(nullArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteZoneFacesSpecial









! ******************************************************************************
!
! Purpose: Write zone for surface data.
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

SUBROUTINE RFLU_TEC_WriteZoneInterf(pRegion)

  USE RFLU_ModCellFaceEdgeInfo

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
 
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: zonePrefix,zoneTitle 
  INTEGER :: c1,c1k,c2,c2k,errorFlag,ifg,ifl,iQuad,iTri,iVar, &
             iVarFace,iVarVert,ivl,nFaces,nFacesEst,nQuads,nTris,nVert, &
             nVertEst,v1,v2,v3,v4,zoneType
  INTEGER, DIMENSION(:), ALLOCATABLE :: bvFlag,nullArray
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: iBFace,f2v,f2vTmp 
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: varFaceTmp,varVertTmp
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECDAT100,TECNOD100,TECZNE100
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteZoneInterf', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Create null array
! ******************************************************************************

  ALLOCATE(nullArray(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

  DO iVar = 1,nVarsTEC
    nullArray(iVar) = 0
  END DO ! iVar 
 
! ******************************************************************************
! Allocate temporary memory. NOTE include pretty large safety factor, because 
! the difference between the total and actual number of vertices can lead to 
! underestimation of interface vertices if have small number of virtual cells,
! (as for first-order scheme with no grid motion), and hence small number of 
! virtual vertices.
! ******************************************************************************

  nVertEst = INT(10.0_RFREAL*(pGrid%nVertTot - pGrid%nVert))

  ALLOCATE(bvFlag(pGrid%nVertTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bvFlag')
  END IF ! global%error   
            
  ALLOCATE(varVertTmp(nVertEst,nVarsVertTEC),STAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varVertTmp')
  END IF ! global%error  
    
! ******************************************************************************
! Allocate memory for interpartition face list
! ******************************************************************************
  
  IF ( pGrid%nFacesTot /= pGrid%nFaces ) THEN 
    nFacesEst = pGrid%nFacesTot - pGrid%nFaces
  ELSE ! May have just one face
    nFacesEst = 1  
  END IF ! pGrid%nFacesTot

  ALLOCATE(f2v(4,nFacesEst),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2v')
  END IF ! global%error  

! ******************************************************************************
! Build temporary interpartition face list
! ******************************************************************************
 
  nFaces = 0
  nTris  = 0
  nQuads = 0    

  DO ifg = 1,pGrid%nFaces ! NOTE upper limit
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

    c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
    c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

    IF ( RFLU_GetFaceKind(global,c1k,c2k) == FACE_KIND_AV ) THEN 
      nFaces = nFaces + 1       

      IF ( nFaces > nFacesEst ) THEN 
        CALL ErrorStop(global,ERR_NFACES_ESTIMATE,__LINE__)
      END IF ! nFaces

      f2v(1:4,nFaces) = pGrid%f2v(1:4,ifg)

      IF ( f2v(4,nFaces) == VERT_NONE ) THEN 
        nTris  = nTris  + 1
      ELSE 
        nQuads = nQuads + 1
      END IF ! f2v
    END IF ! RFLU_GetFaceKind
  END DO ! ifg

! ******************************************************************************
! Build final interpartition face list (note sorting)
! ******************************************************************************
 
  ALLOCATE(f2vTmp(4,nFaces),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2vTmp')
  END IF ! global%error  

  f2vTmp(1:4,1:nFaces) = f2v(1:4,1:nFaces)

  DEALLOCATE(f2v,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2v')
  END IF ! global%error    

  ALLOCATE(f2v(4,nFaces),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2v')
  END IF ! global%error   

  iTri  = 0
  iQuad = nTris   

  DO ifg = 1,nFaces
    IF ( f2vTmp(4,ifg) == VERT_NONE ) THEN 
      iTri = iTri + 1
      f2v(1:4,iTri) = f2vTmp(1:4,ifg)      
    ELSE
      iQuad = iQuad + 1
      f2v(1:4,iQuad) = f2vTmp(1:4,ifg) 
    END IF ! f2vTmp      
  END DO ! ifg

  DEALLOCATE(f2vTmp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2vTmp')
  END IF ! global%error     

! ******************************************************************************
! Allocate memory for face variables
! ******************************************************************************

  IF ( nVarsFaceTEC > 0 ) THEN 
    ALLOCATE(varFaceTmp(nFaces,nVarsFaceTEC),STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varFaceTmp')
    END IF ! global%error  
  END IF ! nVarsFaceTEC

! ******************************************************************************
! Triangular interpartition faces
! ******************************************************************************
 
  IF ( nTris > 0 ) THEN 
    ALLOCATE(f2vTmp(3,nTris),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2vTmp')
    END IF ! global%error   
    
    DO ivl = 1,SIZE(bvFlag,1)
      bvFlag(ivl) = 0 
    END DO ! ivl

    nVert = 0

    DO ifl = 1,nTris
      v1 = f2v(1,ifl)
      v2 = f2v(2,ifl)
      v3 = f2v(3,ifl)

      IF ( bvFlag(v1) == 0 ) THEN
        nVert = nVert + 1
        bvFlag(v1) = nVert

        iVarVert = 0

        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
            iVarVert = iVarVert + 1 
            
            IF ( nVert > nVertEst ) THEN 
              CALL ErrorStop(global,ERR_NBVERT_ESTIMATE,__LINE__)
            END IF ! nVert
                                    
            varVertTmp(nVert,iVarVert) = pRegion%varVertTEC(v1,iVarVert)
          END IF ! posTEC        
        END DO ! iVar
      END IF ! bvFlag

      IF ( bvFlag(v2) == 0 ) THEN
        nVert = nVert + 1
        bvFlag(v2) = nVert

        iVarVert = 0

        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
            iVarVert = iVarVert + 1 
                                
            IF ( nVert > nVertEst ) THEN 
              CALL ErrorStop(global,ERR_NBVERT_ESTIMATE,__LINE__)
            END IF ! nVert                                
                                    
            varVertTmp(nVert,iVarVert) = pRegion%varVertTEC(v2,iVarVert)
          END IF ! posTEC        
        END DO ! iVar
      END IF ! bvFlag

      IF ( bvFlag(v3) == 0 ) THEN
        nVert = nVert + 1
        bvFlag(v3) = nVert

        iVarVert = 0

        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
            iVarVert = iVarVert + 1 
            
            IF ( nVert > nVertEst ) THEN 
              CALL ErrorStop(global,ERR_NBVERT_ESTIMATE,__LINE__)
            END IF ! nVert            
            
            varVertTmp(nVert,iVarVert) = pRegion%varVertTEC(v3,iVarVert)
          END IF ! posTEC        
        END DO ! iVar
      END IF ! bvFlag               

      f2vTmp(1,ifl) = bvFlag(v1)
      f2vTmp(2,ifl) = bvFlag(v2)
      f2vTmp(3,ifl) = bvFlag(v3)
    END DO ! ifl

! ******************************************************************************
!   Write data
! ******************************************************************************

! ==============================================================================
!   Zone information
! ==============================================================================
 
    WRITE(zoneTitle,'(A,I5.5)') 'INT_TRI_',pRegion%iRegionGlobal
    errorFlag = TECZNE100(TRIM(zoneTitle)//NULLCHAR,ZONE_TYPE_TRI,nVert, & 
                          nTris,0,0,0,0,ZONE_FORM_BLOCK,0,0,posTEC,nullArray,0)    
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error                           

! ==============================================================================
!   Variables
! ==============================================================================
 
    iVarFace = 0
    iVarVert = 0    

    DO iVar = 1,nVarsTEC
      IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
        iVarVert = iVarVert + 1
        errorFlag = TECDAT100(nVert,varVertTmp(1,iVarVert),doubleFlagTEC)
      ELSE IF ( posTEC(iVar) == VAR_POS_FACE ) THEN 
        iVarFace = iVarFace + 1
        errorFlag = TECDAT100(nTris,varFaceTmp(1,iVarFace),doubleFlagTEC)        
      ELSE 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! posTEC
    
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
      END IF ! global%error       
    END DO ! iVar
   
! ==============================================================================
!   Connectivity
! ==============================================================================

    errorFlag = TECNOD100(f2vTmp)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error 

    DEALLOCATE(f2vTmp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2vTmp')
    END IF ! global%error    
  END IF ! nTris

! ******************************************************************************
! Quadrilateral interpartition faces
! ******************************************************************************
 
  IF ( nQuads > 0 ) THEN
    ALLOCATE(f2vTmp(4,nQuads),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2vTmp')
    END IF ! global%error   
    
    DO ivl = 1,SIZE(bvFlag,1)
      bvFlag(ivl) = 0 
    END DO ! ivl
    
    nVert = 0
 
    DO ifl = 1,nQuads
      v1 = f2v(1,nTris+ifl)
      v2 = f2v(2,nTris+ifl)
      v3 = f2v(3,nTris+ifl)
      v4 = f2v(4,nTris+ifl)

      IF ( bvFlag(v1) == 0 ) THEN
        nVert = nVert + 1
        bvFlag(v1) = nVert

        iVarVert = 0

        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
            iVarVert = iVarVert + 1 
            
            IF ( nVert > nVertEst ) THEN 
              CALL ErrorStop(global,ERR_NBVERT_ESTIMATE,__LINE__)
            END IF ! nVert            
            
            varVertTmp(nVert,iVarVert) = pRegion%varVertTEC(v1,iVarVert)
          END IF ! posTEC        
        END DO ! iVar          
      END IF ! bvFlag

      IF ( bvFlag(v2) == 0 ) THEN
        nVert = nVert + 1
        bvFlag(v2) = nVert

        iVarVert = 0

        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
            iVarVert = iVarVert + 1 
            
            IF ( nVert > nVertEst ) THEN 
              CALL ErrorStop(global,ERR_NBVERT_ESTIMATE,__LINE__)
            END IF ! nVert            
            
            varVertTmp(nVert,iVarVert) = pRegion%varVertTEC(v2,iVarVert)
          END IF ! posTEC        
        END DO ! iVar 
      END IF ! bvFlag

      IF ( bvFlag(v3) == 0 ) THEN
        nVert = nVert + 1
        bvFlag(v3) = nVert

        iVarVert = 0

        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
            iVarVert = iVarVert + 1 
            varVertTmp(nVert,iVarVert) = pRegion%varVertTEC(v3,iVarVert)
          END IF ! posTEC        
        END DO ! iVar     
      END IF ! bvFlag

      IF ( bvFlag(v4) == 0 ) THEN
        nVert = nVert + 1
        bvFlag(v4) = nVert

        iVarVert = 0

        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_VERT ) THEN
            iVarVert = iVarVert + 1 
            
            IF ( nVert > nVertEst ) THEN 
              CALL ErrorStop(global,ERR_NBVERT_ESTIMATE,__LINE__)
            END IF ! nVert
                        
            varVertTmp(nVert,iVarVert) = pRegion%varVertTEC(v4,iVarVert)
          END IF ! posTEC        
        END DO ! iVar 
      END IF ! bvFlag

      f2vTmp(1,ifl) = bvFlag(v1)
      f2vTmp(2,ifl) = bvFlag(v2)
      f2vTmp(3,ifl) = bvFlag(v3)
      f2vTmp(4,ifl) = bvFlag(v4)
    END DO ! ifl

! ******************************************************************************
!   Write data
! ******************************************************************************

! ==============================================================================
!   Zone information
! ==============================================================================
 
    WRITE(zoneTitle,'(A,I5.5)') 'INT_QUAD_',pRegion%iRegionGlobal
    errorFlag = TECZNE100(TRIM(zoneTitle)//NULLCHAR,ZONE_TYPE_QUAD,nVert, & 
                          nQuads,0,0,0,0,ZONE_FORM_BLOCK,0,0,posTEC,nullArray,0)    
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error                           

! ==============================================================================
!   Variables
! ==============================================================================
 
    iVarFace = 0
    iVarVert = 0    

    DO iVar = 1,nVarsTEC
      IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
        iVarVert = iVarVert + 1
        errorFlag = TECDAT100(nVert,varVertTmp(1,iVarVert),doubleFlagTEC)
      ELSE IF ( posTEC(iVar) == VAR_POS_FACE ) THEN 
        iVarFace = iVarFace + 1
        errorFlag = TECDAT100(nQuads,varFaceTmp(1,iVarFace),doubleFlagTEC)        
      ELSE 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! posTEC
    
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
      END IF ! global%error       
    END DO ! iVar

! ==============================================================================
!   Connectivity
! ==============================================================================
 
    errorFlag = TECNOD100(f2vTmp)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error 

    DEALLOCATE(f2vTmp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2vTmp')
    END IF ! global%error   
  END IF ! nQuads  

  DEALLOCATE(f2v,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2v')
  END IF ! global%error  

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************
 
  DEALLOCATE(bvFlag,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bvFlag')
  END IF ! global%error   

  DEALLOCATE(varVertTmp,STAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varVertTmp')
  END IF ! global%error  

  IF ( nVarsFaceTEC > 0 ) THEN 
    DEALLOCATE(varFaceTmp,STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varFaceTmp')
    END IF ! global%error  
  END IF ! nVarsFaceTEC
  
! ******************************************************************************
! Destroy null array 
! ******************************************************************************
 
  DEALLOCATE(nullArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteZoneInterf









! ******************************************************************************
!
! Purpose: Write zone for surface data.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!   faceType            Face type
!   faceKind            Face kind
! 
! Output: None.
!
! Notes: 
!   1. This routine combines writing of surface data irrespective of whether
!      volume data is written.
!
! ******************************************************************************  

SUBROUTINE RFLU_TEC_WriteZoneSurf(pRegion,pPatch,faceType,faceKind)

  USE ModSortSearch

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
 
  INTEGER, INTENT(IN) :: faceKind,faceType
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: zonePrefix,zoneTitle 
  INTEGER :: errorFlag,ifg,ifgBeg,ifgEnd,ifgOff,ifg2,ifl,iLoc,iVar,iVarFace, &
             iVarVert,ivg,ivlf,ivlp,nBFaces,nBVert,zoneType
  INTEGER, DIMENSION(:), ALLOCATABLE :: bvFlag,nullArray
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: bf2vTmp 
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: varFaceTmp,varVertTmp
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECDAT100,TECNOD100,TECZNE100
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteZoneSurf', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid
  
! ******************************************************************************
! Create null array
! ******************************************************************************

  ALLOCATE(nullArray(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

  DO iVar = 1,nVarsTEC
    nullArray(iVar) = 0
  END DO ! iVar 
 
! ******************************************************************************
! Set zone type, number of faces, and zone prefix. Face offset must be defined 
! so that boundary face data can be accessed properly. The difference arises 
! because the connectivity is stored for each face type, whereas the data is 
! stored for the combined face list (actual triangles, actual quads, virtual
! triangles, virtual quads).
! ******************************************************************************
  
  SELECT CASE ( faceKind ) 
    CASE ( FACE_KIND_AB ) 
      SELECT CASE ( faceType ) 
        CASE ( FACE_TYPE_TRI )
          ifgBeg     = 1
          ifgEnd     = pPatch%nBTris
          ifgOff     = 0 
          zoneType   = ZONE_TYPE_TRI 
          zonePrefix = 'TRI-A_'           
        CASE ( FACE_TYPE_QUAD ) 
          ifgBeg     = 1
          ifgEnd     = pPatch%nBQuads
          ifgOff     = pPatch%nBTris             
          zoneType   = ZONE_TYPE_QUAD            
          zonePrefix = 'QUAD-A_'
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
      END SELECT ! faceType      
    CASE ( FACE_KIND_VB ) 
      SELECT CASE ( faceType ) 
        CASE ( FACE_TYPE_TRI ) 
          ifgBeg     = pPatch%nBTris+1
          ifgEnd     = pPatch%nBTrisTot
          ifgOff     = pPatch%nBTris+pPatch%nBQuads 
          zoneType   = ZONE_TYPE_TRI 
          zonePrefix = 'TRI-V_'           
        CASE ( FACE_TYPE_QUAD ) 
          ifgBeg     = pPatch%nBQuads+1
          ifgEnd     = pPatch%nBQuadsTot           
          ifgOff     = pPatch%nBTrisTot+pPatch%nBQuads
          zoneType   = ZONE_TYPE_QUAD            
          zonePrefix = 'QUAD-V_'
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
      END SELECT ! cellType
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)       
  END SELECT ! cellKind   
  
  nBFaces = ifgEnd - ifgBeg + 1

! ******************************************************************************
! Allocate temporary memory 
! ******************************************************************************

  ALLOCATE(bvFlag(pPatch%nBVertTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bvFlag')
  END IF ! global%error   
 
  IF ( nVarsFaceTEC > 0 ) THEN 
    ALLOCATE(varFaceTmp(nBFaces,nVarsFaceTEC),STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varFaceTmp')
    END IF ! global%error  
  END IF ! nVarsFaceTEC 
 
! ******************************************************************************
! Construct local connectivity array. NOTE this loop computes nBVert, which will
! differ from pPatch%nBVert if have triangular and quadrilateral faces on same 
! patch or actual and virtual faces on same patch.
! ******************************************************************************

! ==============================================================================
! Triangles
! ==============================================================================

  IF ( faceType == FACE_TYPE_TRI ) THEN 
    ALLOCATE(bf2vTmp(3,nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bf2vTmp')
    END IF ! global%error   
    
    DO ivlp = 1,pPatch%nBVertTot
      bvFlag(ivlp) = 0 
    END DO ! ivlp
    
    ifl    = 0
    nBVert = 0

! ------------------------------------------------------------------------------
!   Loop over faces
! ------------------------------------------------------------------------------

    DO ifg = ifgBeg,ifgEnd
      ifl  = ifl + 1
      ifg2 = ifg - ifgBeg + ifgOff + 1
      
      DO ivlf = 1,3
        ivg = pPatch%bTri2v(ivlf,ifg)
        
        CALL BinarySearchInteger(pPatch%bv,pPatch%nBVertTot,ivg,ivlp)
        
        IF ( ivlp == ELEMENT_NOT_FOUND ) THEN 
          CALL ErrorStop(global,ERR_BINARY_SEARCH,__LINE__)
        END IF ! ivlp            

        IF ( bvFlag(ivlp) == 0 ) THEN
          nBVert = nBVert + 1
          bvFlag(ivlp) = nBVert               
        END IF ! bvFlag
      
! ----- Define connectivity ----------------------------------------------------      
      
        bf2vTmp(ivlf,ifl) = bvFlag(ivlp)

! ----- Copy face data ---------------------------------------------------------
     
        iVarFace = 0
    
        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_FACE ) THEN 
            iVarFace = iVarFace + 1
            varFaceTmp(ifl,iVarFace) = pPatch%varFaceTEC(ifg2,iVarFace)
          END IF ! posTEC
        END DO ! iVar                            
      END DO ! ivlf
    END DO ! ifg

! ==============================================================================
! Quadrilaterals
! ==============================================================================

  ELSE 
    ALLOCATE(bf2vTmp(4,nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bf2vTmp')
    END IF ! global%error   

    DO ivlp = 1,pPatch%nBVertTot
      bvFlag(ivlp) = 0 
    END DO ! ivlp    
    
    ifl    = 0 
    nBVert = 0

! ------------------------------------------------------------------------------
!   Loop over faces
! ------------------------------------------------------------------------------

    DO ifg = ifgBeg,ifgEnd
      ifl  = ifl + 1
      ifg2 = ifg - ifgBeg + ifgOff + 1
      
      DO ivlf = 1,4
        ivg = pPatch%bQuad2v(ivlf,ifg)
        
        CALL BinarySearchInteger(pPatch%bv,pPatch%nBVertTot,ivg,ivlp)
        
        IF ( ivlp == ELEMENT_NOT_FOUND ) THEN 
          CALL ErrorStop(global,ERR_BINARY_SEARCH,__LINE__)
        END IF ! ivlp             

        IF ( bvFlag(ivlp) == 0 ) THEN
          nBVert = nBVert + 1
          bvFlag(ivlp) = nBVert               
        END IF ! bvFlag
      
! ----- Define connectivity ----------------------------------------------------      
      
        bf2vTmp(ivlf,ifl) = bvFlag(ivlp)

! ----- Copy face data ---------------------------------------------------------

        iVarFace = 0
    
        DO iVar = 1,nVarsTEC
          IF ( posTEC(iVar) == VAR_POS_FACE ) THEN 
            iVarFace = iVarFace + 1
            varFaceTmp(ifl,iVarFace) = pPatch%varFaceTEC(ifg2,iVarFace)
          END IF ! posTEC
        END DO ! iVar
      END DO ! ivlf
    END DO ! ifg
  END IF ! faceType

! ******************************************************************************
! Allocate memory for and copy vertex data
! ******************************************************************************

  ALLOCATE(varVertTmp(nBVert,nVarsVertTEC),STAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varVertTmp')
  END IF ! global%error  
    
  DO ivlp = 1,pPatch%nBVertTot
    IF ( bvFlag(ivlp) /= 0 ) THEN
      iVarVert = 0
    
      DO iVar = 1,nVarsTEC
        IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
          iVarVert = iVarVert + 1
          varVertTmp(bvFlag(ivlp),iVarVert) = pPatch%varVertTEC(ivlp,iVarVert)
        END IF ! posTEC
      END DO ! iVar
    END IF ! bvFlag
  END DO ! ivlp

! ******************************************************************************
! Write to TECPLOT file
! ******************************************************************************

! ==============================================================================
! Zone information
! ==============================================================================

  WRITE(zoneTitle,'(A,I3.3,A,I5.5)') 'PAT_',pPatch%iPatchGlobal, & 
                                      '_'//TRIM(zonePrefix), &
                                      pRegion%iRegionGlobal

  errorFlag = TECZNE100(TRIM(zoneTitle)//NULLCHAR,zoneType,nBVert,nBFaces, &
                        0,0,0,0,ZONE_FORM_BLOCK,0,0,posTEC,nullArray,0)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error              

! ==============================================================================
! Variables
! ==============================================================================

  iVarFace = 0
  iVarVert = 0

  DO iVar = 1,nVarsTEC
    IF ( posTEC(iVar) == VAR_POS_FACE ) THEN 
      iVarFace = iVarFace + 1    
      errorFlag = TECDAT100(nBFaces,varFaceTmp(1,iVarFace),doubleFlagTEC)     
    ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
      iVarVert = iVarVert + 1    
      errorFlag = TECDAT100(nBVert,varVertTmp(1,iVarVert),doubleFlagTEC)      
    ELSE ! defensive coding
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! posTEC
  
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error       
  END DO ! iVar

! ==============================================================================
! Connectivity
! ==============================================================================

  errorFlag = TECNOD100(bf2vTmp)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error 

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(bvFlag,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bvFlag')
  END IF ! global%error   

  IF ( ALLOCATED(varFaceTmp) .EQV. .TRUE. ) THEN 
    DEALLOCATE(varFaceTmp,STAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varFaceTmp')
    END IF ! global%error  
  END IF ! ALLOCATED

  DEALLOCATE(varVertTmp,STAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varVertTmp')
  END IF ! global%error  

  DEALLOCATE(bf2vTmp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bf2vTmp')
  END IF ! global%error   

! ******************************************************************************
! Destroy null array 
! ******************************************************************************

  DEALLOCATE(nullArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteZoneSurf









! ******************************************************************************
!
! Purpose: Write zone for surface data in mixed format.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
! 
! Output: None.
!
! Notes: 
!   1. The term mixed format refers to the list of faces being written out in
!      a single block, without distinguishing between triangular and 
!      quadrilateral faces. 
!   2. Triangular faces are written as degenerate quadrilateral faces.
!   3. Only actual faces are written because face data is typically not defined 
!      for virtual faces.  
!
! ******************************************************************************  

SUBROUTINE RFLU_TEC_WriteZoneSurfMixed(pRegion,pPatch)

  USE ModSortSearch

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
 
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: zoneTitle 
  INTEGER :: errorFlag,ifl,iVar,iVarFace,iVarVert,ivl
  INTEGER, DIMENSION(:), ALLOCATABLE :: nullArray
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: bf2vTmp  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECDAT100,TECNOD100,TECZNE100
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteZoneSurfMixed', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid
  
! ******************************************************************************
! Create null array
! ******************************************************************************

  ALLOCATE(nullArray(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

  DO iVar = 1,nVarsTEC
    nullArray(iVar) = 0
  END DO ! iVar 
 
! ******************************************************************************
! Allocate temporary memory 
! ******************************************************************************
 
  ALLOCATE(bf2vTmp(4,pPatch%nBFaces),STAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bf2vTmp')
  END IF ! global%error  

! ******************************************************************************
! Convert connectivity array
! ******************************************************************************

  DO ifl = 1,pPatch%nBFaces
    bf2vTmp(1,ifl) = pPatch%bf2v(1,ifl)
    bf2vTmp(2,ifl) = pPatch%bf2v(2,ifl)
    bf2vTmp(3,ifl) = pPatch%bf2v(3,ifl)
  
    IF ( pPatch%bf2v(4,ifl) /= VERT_NONE ) THEN
      bf2vTmp(4,ifl) = pPatch%bf2v(4,ifl)
    ELSE 
      bf2vTmp(4,ifl) = bf2vTmp(1,ifl)
    END IF ! v4l            
  END DO ! ifl  

! ******************************************************************************
! Write to TECPLOT file
! ******************************************************************************

! ==============================================================================
! Zone information
! ==============================================================================

  WRITE(zoneTitle,'(A,I3.3,A,I5.5)') 'PAT_',pPatch%iPatchGlobal, & 
                                      '_MIX-A_',pRegion%iRegionGlobal
                                                                                                              
  errorFlag = TECZNE100(TRIM(zoneTitle)//NULLCHAR,ZONE_TYPE_QUAD, &
                        pPatch%nBVert,pPatch%nBFaces,0,0,0,0,ZONE_FORM_BLOCK, &
                        0,0,posTEC,nullArray,0)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error              

! ==============================================================================
! Variables
! ==============================================================================

  iVarFace = 0
  iVarVert = 0

  DO iVar = 1,nVarsTEC 
    IF ( posTEC(iVar) == VAR_POS_FACE ) THEN 
      iVarFace = iVarFace + 1
      errorFlag = TECDAT100(pPatch%nBFaces,pPatch%varFaceTEC(1,iVarFace), & 
                            doubleFlagTEC)      
    ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
      iVarVert = iVarVert + 1    
      errorFlag = TECDAT100(pPatch%nBVert,pPatch%varVertTEC(1,iVarVert), & 
                            doubleFlagTEC)      
    ELSE ! defensive coding
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! posTEC
  
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error       
  END DO ! iVar

! ==============================================================================
! Connectivity
! ==============================================================================

  errorFlag = TECNOD100(bf2vTmp)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error 

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(bf2vTmp,STAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bf2vTmp')
  END IF ! global%error  
    
! ******************************************************************************
! Destroy null array 
! ******************************************************************************

  DEALLOCATE(nullArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteZoneSurfMixed






! ******************************************************************************
!
! Purpose: Write zone for volume data.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   cellType            Cell type
!   cellKind            Cell kind
! 
! Output: None.
!
! Notes:
!   1. Prisms and pyramids are treated as degenerate hexahedra.
!
! ******************************************************************************  

SUBROUTINE RFLU_TEC_WriteZoneVol(pRegion,cellType,cellKind)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
 
  INTEGER, INTENT(IN) :: cellKind,cellType
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: zonePrefix,zoneTitle 
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,icl,iVar,iVarCell,iVarVert,nCells, &
             nVert,zoneType
  INTEGER, DIMENSION(:), ALLOCATABLE :: nullArray
  INTEGER, DIMENSION(:), POINTER :: pXyz2CellGlob
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: c2v 
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: varCellTmp
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
! Externals: TECPLOT functions
! ==============================================================================

  INTEGER, EXTERNAL :: TECDAT100,TECNOD100,TECZNE100
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TEC_WriteZoneVol', & 
                        'RFLU_ModTECPLOTUtils.F90')

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Create null array
! ******************************************************************************

  ALLOCATE(nullArray(nVarsTEC),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

  DO iVar = 1,nVarsTEC
    nullArray(iVar) = 0
  END DO ! iVar 
 
! ******************************************************************************
! Set zone type, number of cells, and zone title
! ******************************************************************************
  
  SELECT CASE ( cellKind ) 
    CASE ( CELL_KIND_ACTUAL ) 
      SELECT CASE ( cellType ) 
        CASE ( CELL_TYPE_TET )
          icgBeg     = 1
          icgEnd     = pGrid%nTets
          nVert      = pGrid%nVertTot
          zoneType   = ZONE_TYPE_TET 
          zonePrefix = 'TET-A_'           
        CASE ( CELL_TYPE_HEX ) 
          icgBeg     = 1
          icgEnd     = pGrid%nHexs       
          nVert      = pGrid%nVertTot          
          zoneType   = ZONE_TYPE_HEX            
          zonePrefix = 'HEX-A_'
        CASE ( CELL_TYPE_PRI ) 
          icgBeg     = 1
          icgEnd     = pGrid%nPris                
          nVert      = pGrid%nVertTot           
          zoneType   = ZONE_TYPE_HEX          
          zonePrefix = 'PRI-A_'
        CASE ( CELL_TYPE_PYR ) 
          icgBeg     = 1
          icgEnd     = pGrid%nPyrs                        
          nVert      = pGrid%nVertTot          
          zoneType   = ZONE_TYPE_HEX            
          zonePrefix = 'PYR-A_'
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
      END SELECT ! cellType      
    CASE ( CELL_KIND_VIRTUAL ) 
      SELECT CASE ( cellType ) 
        CASE ( CELL_TYPE_TET ) 
          icgBeg     = pGrid%nTets+1
          icgEnd     = pGrid%nTetsTot                                 
          nVert      = pGrid%nVertTot
          zoneType   = ZONE_TYPE_TET 
          zonePrefix = 'TET-V_'           
        CASE ( CELL_TYPE_HEX ) 
          icgBeg     = pGrid%nHexs+1
          icgEnd     = pGrid%nHexsTot                                 
          nVert      = pGrid%nVertTot          
          zoneType   = ZONE_TYPE_HEX            
          zonePrefix = 'HEX-V_'
        CASE ( CELL_TYPE_PRI ) 
          icgBeg     = pGrid%nPris+1
          icgEnd     = pGrid%nPrisTot                                 
          nVert      = pGrid%nVertTot           
          zoneType   = ZONE_TYPE_HEX          
          zonePrefix = 'PRI-V_'
        CASE ( CELL_TYPE_PYR ) 
          icgBeg     = pGrid%nPyrs+1
          icgEnd     = pGrid%nPyrsTot                                 
          nVert      = pGrid%nVertTot          
          zoneType   = ZONE_TYPE_HEX            
          zonePrefix = 'PYR-V_'
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
      END SELECT ! cellType 
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
  END SELECT ! cellKind   
  
  nCells = icgEnd - icgBeg + 1
  
! ******************************************************************************
! Convert cell array. This is required whenever have mixed grid because nCells
! would then not match pGrid%nCellsTot, which is the dimension of the array 
! pRegion%varCellTEC.
! ******************************************************************************
  
  IF ( nVarsCellTEC > 0 ) THEN 
    ALLOCATE(varCellTmp(nCells,nVarsCellTEC),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varCellTmp')
    END IF ! global%error    
    
    SELECT CASE ( cellType )     
      CASE ( CELL_TYPE_TET ) 
        pXyz2CellGlob => pGrid%tet2CellGlob
      CASE ( CELL_TYPE_HEX ) 
        pXyz2CellGlob => pGrid%hex2CellGlob      
      CASE ( CELL_TYPE_PRI ) 
        pXyz2CellGlob => pGrid%pri2CellGlob       
      CASE ( CELL_TYPE_PYR ) 
        pXyz2CellGlob => pGrid%pyr2CellGlob       
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! cellType
    
    DO icg = icgBeg,icgEnd
      icl = icg - icgBeg + 1 

      DO iVarCell = 1,nVarsCellTEC
        varCellTmp(icl,iVarCell) = & 
          pRegion%varCellTEC(pXyz2CellGlob(icg),iVarCell)
      END DO ! iVarCell
    END DO ! icl               
  END IF ! nVarsCellTEC
  
! ******************************************************************************
! Convert connectivity for prisms and pyramids
! ******************************************************************************

  IF ( cellType == CELL_TYPE_PRI .OR. cellType == CELL_TYPE_PYR ) THEN 
    ALLOCATE(c2v(8,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'c2v')
    END IF ! global%error                   

    IF ( cellType == CELL_TYPE_PRI ) THEN ! Prisms
      DO icg = icgBeg,icgEnd
        icl = icg - icgBeg + 1
        c2v(1:3,icl) = pGrid%pri2v(1:3,icg)
        c2v(4  ,icl) = pGrid%pri2v(1  ,icg)
        c2v(5:7,icl) = pGrid%pri2v(4:6,icg)
        c2v(8  ,icl) = pGrid%pri2v(4  ,icg)
      END DO ! icg 
    ELSE ! Pyramids
      DO icg = icgBeg,icgEnd
        icl = icg - icgBeg + 1      
        c2v(1:5,icl) = pGrid%pyr2v(1:5,icg)
        c2v(6  ,icl) = pGrid%pyr2v(5  ,icg)
        c2v(7  ,icl) = pGrid%pyr2v(5  ,icg)
        c2v(8  ,icl) = pGrid%pyr2v(5  ,icg)
      END DO ! icg      
    END IF ! cellType   
  END IF ! cellType

! ******************************************************************************
! Write to TECPLOT file
! ******************************************************************************
  
! ==============================================================================
! Zone information
! ==============================================================================

  WRITE(zoneTitle,'(A,I5.5)') TRIM(zonePrefix),pRegion%iRegionGlobal

  errorFlag = TECZNE100(TRIM(zoneTitle)//NULLCHAR,zoneType,nVert,nCells, &
                        0,0,0,0,ZONE_FORM_BLOCK,0,0,posTEC,nullArray,0)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error              

! ==============================================================================
! Variables
! ==============================================================================

  iVarCell = 0
  iVarVert = 0

  DO iVar = 1,nVarsTEC
    IF ( posTEC(iVar) == VAR_POS_CELL ) THEN 
      iVarCell = iVarCell + 1
      errorFlag = TECDAT100(nCells,varCellTmp(1,iVarCell),doubleFlagTEC)                                  
    ELSE IF ( posTEC(iVar) == VAR_POS_VERT ) THEN 
      iVarVert = iVarVert + 1    
      errorFlag = TECDAT100(pGrid%nVertTot,pRegion%varVertTEC(1,iVarVert), & 
                            doubleFlagTEC)      
    ELSE ! defensive coding
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! posTEC
  
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
    END IF ! global%error       
  END DO ! iVar

! ==============================================================================
! Connectivity
! ==============================================================================

  SELECT CASE ( cellType ) 
    CASE ( CELL_TYPE_TET ) 
      errorFlag = TECNOD100(pGrid%tet2v(1:4,icgBeg:icgEnd))
    CASE ( CELL_TYPE_HEX ) 
      errorFlag = TECNOD100(pGrid%hex2v(1:8,icgBeg:icgEnd))
    CASE ( CELL_TYPE_PRI ) 
      errorFlag = TECNOD100(c2v)
    CASE ( CELL_TYPE_PYR ) 
      errorFlag = TECNOD100(c2v)
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END SELECT ! cellType

  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_TECPLOT_OUTPUT,__LINE__)
  END IF ! global%error

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  IF ( nVarsCellTEC > 0 ) THEN 
    DEALLOCATE(varCellTmp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varCellTmp')
    END IF ! global%error                  
  END IF ! nVarsCellTEC
  
  IF ( cellType == CELL_TYPE_PRI .OR. cellType == CELL_TYPE_PYR ) THEN 
    DEALLOCATE(c2v,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'c2v')
    END IF ! global%error                      
  END IF ! cellType

! ******************************************************************************
! Destroy null array 
! ******************************************************************************

  DEALLOCATE(nullArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nullArray')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TEC_WriteZoneVol






END MODULE RFLU_ModTECPLOTUtils 

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModTECPLOTUtils.F90,v $
! Revision 1.14  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.11  2005/09/23 19:02:21  haselbac
! Added FILE_TYPE_PATCH_STATS
!
! Revision 1.10  2005/03/16 20:31:44  haselbac
! Changed naming of special face zones to be more consistent
!
! Revision 1.9  2005/01/06 04:43:17  haselbac
! Two bug fixes to allow cell data to be written correctly for mixed grids
!
! Revision 1.8  2004/12/27 23:34:42  haselbac
! Added writing of field cell and face data
!
! Revision 1.7  2004/12/21 15:11:04  fnajjar
! Increased Tecplot header size to accommodate for PLAG surface statistics
!
! Revision 1.6  2004/10/19 19:30:18  haselbac
! Modified RFLU_TEC_WriteZoneSurf, improved estimate of interface vertices
!
! Revision 1.5  2004/09/27 01:42:54  haselbac
! Added procedure to write zone with special faces
!
! Revision 1.4  2004/08/09 16:22:28  haselbac
! Bug fix: Increased safety margin for alloc
!
! Revision 1.3  2004/07/06 15:15:09  haselbac
! Added xyz2vTEC mapping arrays
!
! Revision 1.2  2004/07/02 03:05:39  haselbac
! Bug fix: Changed pNull to nullArray
!
! Revision 1.1  2004/06/16 20:01:40  haselbac
! Initial revision
!
! ******************************************************************************














