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
! Purpose: Collection of routines to read and convert GAMBIT grids.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModGAMBIT.F90,v 1.4 2008/12/06 08:45:03 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGAMBIT

  USE ModParameters
  USE ModDataTypes  
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModGrid

  USE ModBuildFileNames, ONLY: BuildFileNamePlain

  USE RFLU_ModDimensions, ONLY: RFLU_SetMaxDimension

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  PRIVATE

! ==============================================================================
! Data
! ==============================================================================

  TYPE t_patchGAMBIT
    INTEGER :: bcType,nFaces
    INTEGER, DIMENSION(:), POINTER :: bf2ct,bf2cgi,bf2fli
  END TYPE t_patchGAMBIT

  TYPE t_gridGAMBIT
    CHARACTER(CHRLEN) :: title
    INTEGER :: MTYP,NDFCD,NDFVL,NDP,NELGP,NFLAGS,NGP,NGRPS,NTYPE
    INTEGER :: nMappings,nPatches   
    INTEGER, DIMENSION(:), POINTER :: ct
    INTEGER, DIMENSION(:,:), POINTER :: c2v,patch2bc 
    TYPE(t_patchGAMBIT), DIMENSION(:), POINTER :: patches  
  END TYPE t_gridGAMBIT

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGAMBIT.F90,v $ $Revision: 1.4 $'
            
! ------------------------------------------------------------------------------
! GAMBIT element types
! ------------------------------------------------------------------------------    

  INTEGER, PARAMETER :: GAMBIT_NTYPE_EDGE = 1, & 
                        GAMBIT_NTYPE_QUAD = 2, & 
                        GAMBIT_NTYPE_TRI  = 3, & 
                        GAMBIT_NTYPE_HEX  = 4, & 
                        GAMBIT_NTYPE_PRI  = 5, &
                        GAMBIT_NTYPE_TET  = 6, & 
                        GAMBIT_NTYPE_PYR  = 7
    
! ------------------------------------------------------------------------------
! Mapping of GAMBIT face local to cell to corresponding ROCFLU face
! ------------------------------------------------------------------------------    
                               
  INTEGER, DIMENSION(4), PARAMETER :: f2fTetGAMBIT = (/4,1,2,3/)
  INTEGER, DIMENSION(6), PARAMETER :: f2fHexGAMBIT = (/2,3,4,5,1,6/)
  INTEGER, DIMENSION(5), PARAMETER :: f2fPriGAMBIT = (/2,3,4,1,5/) 
  INTEGER, DIMENSION(5), PARAMETER :: f2fPyrGAMBIT = (/1,2,3,4,5/) 

! ------------------------------------------------------------------------------
! Mapping of faces to vertices for GAMBIT
! ------------------------------------------------------------------------------    
                                                   
  INTEGER, DIMENSION(4,4), PARAMETER :: f2vTetGAMBIT = &
    RESHAPE((/2,1,3,VERT_NONE,1,2,4,VERT_NONE,2,3,4,VERT_NONE,3,1,4, & 
              VERT_NONE/), (/4,4/))
  INTEGER, DIMENSION(4,6), PARAMETER :: f2vHexGAMBIT = &
    RESHAPE((/1,2,6,5,2,4,8,6,4,3,7,8,3,1,5,7,2,1,3,4,5,6,7,8/), (/4,6/))
  INTEGER, DIMENSION(4,5), PARAMETER :: f2vPriGAMBIT = &
    RESHAPE((/1,2,5,4,2,3,6,5,3,1,4,6,1,3,2,VERT_NONE,4,5,6, & 
              VERT_NONE/), (/4,5/))
  INTEGER, DIMENSION(4,5), PARAMETER :: f2vPyrGAMBIT = &
    RESHAPE((/1,3,4,2,1,2,5,VERT_NONE,2,4,5,VERT_NONE,4,3,5,VERT_NONE,3,1,5, & 
              VERT_NONE/), (/4,5/))                        
                        
  TYPE(t_gridGAMBIT) :: gridGAMBIT

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_ConvGAMBIT2ROCFLU, & 
            RFLU_ReadGridGAMBITNeutral

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS








! ******************************************************************************
!
! Purpose: Check connectivity.
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

  SUBROUTINE RFLU_CheckGridGAMBIT(pRegion)
     
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

    INTEGER :: errorFlag,icg,ivgMax,ivgMin
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CheckGridGAMBIT', &
                          'RFLU_ModGAMBIT.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Checking connectivity arrays...'    
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Volume grid. NOTE only do volume grid because GAMBIT format specifies 
!   boundary grid in terms of volume grid faces, so checking volume grid only 
!   is ok.
! ******************************************************************************

    DO icg = 1,pGrid%nCellsTot
      SELECT CASE ( gridGAMBIT%ct(icg) ) 
        CASE ( GAMBIT_NTYPE_HEX ) 
          ivgMin = MINVAL(gridGAMBIT%c2v(1:8,icg))
          ivgMax = MINVAL(gridGAMBIT%c2v(1:8,icg))
        CASE ( GAMBIT_NTYPE_TET ) 
          ivgMin = MINVAL(gridGAMBIT%c2v(1:4,icg))
          ivgMax = MINVAL(gridGAMBIT%c2v(1:4,icg))        
        CASE ( GAMBIT_NTYPE_PRI ) 
          ivgMin = MINVAL(gridGAMBIT%c2v(1:6,icg))
          ivgMax = MINVAL(gridGAMBIT%c2v(1:6,icg))      
        CASE ( GAMBIT_NTYPE_PYR ) 
          ivgMin = MINVAL(gridGAMBIT%c2v(1:5,icg))
          ivgMax = MINVAL(gridGAMBIT%c2v(1:5,icg))        
        CASE DEFAULT 
      END SELECT ! gridGAMBIT%ct
            
      IF ( ivgMin < 1 .OR. ivgMax > pGrid%nVertTot ) THEN 
        global%error = ERR_VERTEX_NUMBER      
      END IF ! ivgMin

      IF ( global%error /= ERR_NONE ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN          
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel           
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END DO ! icg

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Checking connectivity arrays done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_CheckGridGAMBIT







! ******************************************************************************
!
! Purpose: Convert grid format from GAMBIT to ROCFLU.
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

  SUBROUTINE RFLU_ConvGAMBIT2ROCFLU(pRegion)
     
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

    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iBegMax,iBegMin,iBeg1,iBeg2,icg,icl,ict,iEndMax, & 
               iEndMin,iEnd1,iEnd2,ifg,iFile,ifl,ifl2,iMap,iMap2,iPatch, &
               iPatch2,ivg,ivl,j
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_patchGAMBIT), POINTER :: pPatchGAMBIT    
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvGAMBIT2ROCFLU', &
                          'RFLU_ModGAMBIT.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from GAMBIT to ROCFLU format...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer and initialize variables
! ******************************************************************************

    pGrid => pRegion%grid

    pGrid%nEdges    = 0
    pGrid%nEdgesTot = 0

    pGrid%nFaces    = 0
    pGrid%nFacesTot = 0

! ******************************************************************************
!   Allocate memory for cell connectivity arrays
! ******************************************************************************

    IF ( pGrid%nTetsMax > 0 ) THEN 
      ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%tet2v)
    END IF ! pGrid%nTetsTot

    IF ( pGrid%nHexsMax > 0 ) THEN 
      ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%hex2v)
    END IF ! pGrid%nHexsTot
    
    IF ( pGrid%nPrisMax > 0 ) THEN 
      ALLOCATE(pGrid%pri2v(6,pGrid%nPrisMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pri2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%pri2v)
    END IF ! pGrid%nPrisTot    
    
    IF ( pGrid%nPyrsMax > 0 ) THEN 
      ALLOCATE(pGrid%pyr2v(5,pGrid%nPyrsMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pyr2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%pyr2v)
    END IF ! pGrid%nPyrsTot    

! ==============================================================================
!   Allocate memory for cell mapping, needed to convert GAMBITs boundary data
!   structure, because for given boundary face, the data structure gives global
!   cell number, not local one.
! ==============================================================================

    ALLOCATE(pGrid%cellGlob2Loc(2,pGrid%nCellsTot),STAT=errorFlag)      
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%cellGlob2Loc')
    END IF ! global%error

! ******************************************************************************
!   Copy connectivity information from GAMBIT format to Rocflu format. NOTE 
!   renumbering of vertices.
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting connectivity...'  
    END IF ! global%verbLevel
        
    pGrid%nTets = 0
    pGrid%nHexs = 0
    pGrid%nPris = 0
    pGrid%nPyrs = 0    

    DO icg = 1,pGrid%nCellsTot
      SELECT CASE ( gridGAMBIT%ct(icg) ) 
        CASE ( GAMBIT_NTYPE_TET ) 
          pGrid%nTets = pGrid%nTets + 1
          
          pGrid%tet2v(1,pGrid%nTets) = gridGAMBIT%c2v(1,icg)
          pGrid%tet2v(2,pGrid%nTets) = gridGAMBIT%c2v(2,icg)
          pGrid%tet2v(3,pGrid%nTets) = gridGAMBIT%c2v(4,icg)
          pGrid%tet2v(4,pGrid%nTets) = gridGAMBIT%c2v(3,icg)
          
          pGrid%cellGlob2Loc(1,icg) = CELL_TYPE_TET
          pGrid%cellGlob2Loc(2,icg) = pGrid%nTets                                                           
        CASE ( GAMBIT_NTYPE_HEX ) 
          pGrid%nHexs = pGrid%nHexs + 1    
          
          pGrid%hex2v(1,pGrid%nHexs) = gridGAMBIT%c2v(1,icg)
          pGrid%hex2v(2,pGrid%nHexs) = gridGAMBIT%c2v(2,icg)
          pGrid%hex2v(3,pGrid%nHexs) = gridGAMBIT%c2v(4,icg)
          pGrid%hex2v(4,pGrid%nHexs) = gridGAMBIT%c2v(3,icg)
          pGrid%hex2v(5,pGrid%nHexs) = gridGAMBIT%c2v(5,icg)
          pGrid%hex2v(6,pGrid%nHexs) = gridGAMBIT%c2v(6,icg)
          pGrid%hex2v(7,pGrid%nHexs) = gridGAMBIT%c2v(8,icg)
          pGrid%hex2v(8,pGrid%nHexs) = gridGAMBIT%c2v(7,icg)                    
                       
          pGrid%cellGlob2Loc(1,icg) = CELL_TYPE_HEX
          pGrid%cellGlob2Loc(2,icg) = pGrid%nHexs                                 
        CASE ( GAMBIT_NTYPE_PRI ) 
          pGrid%nPris = pGrid%nPris + 1    
          
          pGrid%pri2v(1,pGrid%nPris) = gridGAMBIT%c2v(1,icg)
          pGrid%pri2v(2,pGrid%nPris) = gridGAMBIT%c2v(2,icg)  
          pGrid%pri2v(3,pGrid%nPris) = gridGAMBIT%c2v(3,icg)  
          pGrid%pri2v(4,pGrid%nPris) = gridGAMBIT%c2v(4,icg)  
          pGrid%pri2v(5,pGrid%nPris) = gridGAMBIT%c2v(5,icg)  
          pGrid%pri2v(6,pGrid%nPris) = gridGAMBIT%c2v(6,icg)
          
          pGrid%cellGlob2Loc(1,icg) = CELL_TYPE_PRI
          pGrid%cellGlob2Loc(2,icg) = pGrid%nPris                                                                            
        CASE ( GAMBIT_NTYPE_PYR ) 
          pGrid%nPyrs = pGrid%nPyrs + 1      
          
          pGrid%pyr2v(1,pGrid%nPyrs) = gridGAMBIT%c2v(1,icg)
          pGrid%pyr2v(2,pGrid%nPyrs) = gridGAMBIT%c2v(2,icg)
          pGrid%pyr2v(3,pGrid%nPyrs) = gridGAMBIT%c2v(4,icg)
          pGrid%pyr2v(4,pGrid%nPyrs) = gridGAMBIT%c2v(3,icg)
          pGrid%pyr2v(5,pGrid%nPyrs) = gridGAMBIT%c2v(5,icg)
          
          pGrid%cellGlob2Loc(1,icg) = CELL_TYPE_PYR
          pGrid%cellGlob2Loc(2,icg) = pGrid%nPyrs                                          
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! gridGAMBIT%ct
    END DO ! icg    

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting connectivity done.'  
    END IF ! global%verbLevel

! ******************************************************************************
!   Convert patch data structure
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting patch data structure...'  
    END IF ! global%verbLevel

! ==============================================================================
!   Read patch mapping file
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Reading patch mapping file...'  
    END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!   Open file
! ------------------------------------------------------------------------------

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.ggi',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

! ------------------------------------------------------------------------------
!   Read file
! ------------------------------------------------------------------------------
  
    READ(iFile,*) pGrid%nPatches
    READ(iFile,*) gridGAMBIT%nMappings

    ALLOCATE(gridGAMBIT%patch2bc(3,gridGAMBIT%nMappings),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridGAMBIT%patch2bc')
    END IF ! global%error   

    DO iMap = 1,gridGAMBIT%nMappings
      READ(iFile,*) (gridGAMBIT%patch2bc(j,iMap),j=1,3)
    END DO ! iMap

! ------------------------------------------------------------------------------
!   Close file
! ------------------------------------------------------------------------------

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Reading patch mapping file done.'  
    END IF ! global%verbLevel    

! ==============================================================================
!   Check for consistent input - somewhat complicated...
! ==============================================================================

    IF ( global%checkLevel > CHECK_NONE ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Checking patch mapping entries...'
      END IF ! global%verbLevel

      DO iMap = 1,gridGAMBIT%nMappings
        IF ( gridGAMBIT%patch2bc(2,iMap) < gridGAMBIT%patch2bc(1,iMap) ) THEN
          IF ( global%verbLevel > VERBOSE_NONE ) THEN  
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.' 
          END IF ! global%verbLevel   
          CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
        END IF ! gridGAMBIT
      END DO ! iMap   

      IF ( MINVAL(gridGAMBIT%patch2bc(3,:)) /= 1 .OR. & 
           MAXVAL(gridGAMBIT%patch2bc(3,:)) /= pGrid%nPatches ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN              
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel               
        CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
      END IF ! gridGAMBIT

      DO iMap = 1,gridGAMBIT%nMappings
        DO iMap2 = 1,gridGAMBIT%nMappings

          IF ( iMap /= iMap2 ) THEN 
            iBeg1 = gridGAMBIT%patch2bc(1,iMap)
            iEnd1 = gridGAMBIT%patch2bc(2,iMap)

            iBeg2 = gridGAMBIT%patch2bc(1,iMap2)
            iEnd2 = gridGAMBIT%patch2bc(2,iMap2)        

            IF ( iBeg1 < iBeg2 ) THEN 
              iBegMin = iBeg1
              iEndMin = iEnd1
              iBegMax = iBeg2
              iEndMax = iEnd2
            ELSE IF ( iBeg1 > iBeg2 ) THEN 
              iBegMin = iBeg2
              iEndMin = iEnd2
              iBegMax = iBeg1
              iEndMax = iEnd1         
            ELSE ! iBeg1 and iBeg2 have the same value
              IF ( global%verbLevel > VERBOSE_NONE ) THEN
                WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel    
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iBeg1

            IF ( iEndMin >= iBegMax ) THEN
              IF ( global%verbLevel > VERBOSE_NONE ) THEN          
                WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel    
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iEndMin
          END IF ! iMap

        END DO ! iMap2
      END DO ! iMap

      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                 'Checking patch mapping entries done.'
      END IF ! global%verbLevel  
    END IF ! global%checkLevel

! ==============================================================================
!   Allocate patch memory and initialize patch structure
! ==============================================================================

    ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%error       

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      pPatch%nBTris  = 0
      pPatch%nBQuads = 0 
      pPatch%nBVert  = 0  
        
      pPatch%iPatchGlobal = iPatch
      pPatch%iBorder      = PATCH_IBORDER_DEFAULT      
      pPatch%renumFlag    = .FALSE.
    END DO ! iPatch

    global%nPatches = pGrid%nPatches

! ==============================================================================
!   Determine number of faces on each patch and set number of boundary faces
! ==============================================================================
    
    DO iPatch = 1,gridGAMBIT%nPatches
      pPatchGAMBIT => gridGAMBIT%patches(iPatch)
    
      DO iMap = 1,gridGAMBIT%nMappings
        IF ( iPatch >= gridGAMBIT%patch2bc(1,iMap) .AND. & 
             iPatch <= gridGAMBIT%patch2bc(2,iMap) ) THEN 
          iPatch2 = gridGAMBIT%patch2bc(3,iMap)
        END IF ! iPatch
      END DO ! iMap
    
      pPatch => pRegion%patches(iPatch2)
    
      DO ifg = 1,pPatchGAMBIT%nFaces
        ict = pPatchGAMBIT%bf2ct(ifg)
        ifl = pPatchGAMBIT%bf2fli(ifg) 
        
        SELECT CASE ( ict )         
          CASE ( GAMBIT_NTYPE_TET ) 
            pPatch%nBTris = pPatch%nBTris + 1
          CASE ( GAMBIT_NTYPE_HEX ) 
            pPatch%nBQuads = pPatch%nBQuads + 1
          CASE ( GAMBIT_NTYPE_PRI ) 
            IF ( f2vPriGAMBIT(4,ifl) == VERT_NONE ) THEN 
              pPatch%nBTris = pPatch%nBTris + 1           
            ELSE 
              pPatch%nBQuads = pPatch%nBQuads + 1
            END IF ! f2vPriGAMBIT    
          CASE ( GAMBIT_NTYPE_PYR ) 
            IF ( f2vPyrGAMBIT(4,ifl) == VERT_NONE ) THEN 
              pPatch%nBTris = pPatch%nBTris + 1       
            ELSE 
              pPatch%nBQuads = pPatch%nBQuads + 1
            END IF ! f2vPyrGAMBIT        
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! ict
      END DO ! ifg
    END DO ! iPatch
          
! ==============================================================================
!   Set total boundary patch quantities and number of boundary faces
! ==============================================================================

    pGrid%nBFaces = 0 
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      pPatch%nBFaces = pPatch%nBTris + pPatch%nBQuads
      pGrid%nBFaces  = pGrid%nBFaces + pPatch%nBFaces

      pPatch%nBFacesTot = pPatch%nBFaces
      pPatch%nBQuadsTot = pPatch%nBQuads  
      pPatch%nBTrisTot  = pPatch%nBTris    
      pPatch%nBVertTot  = pPatch%nBVert  
      
      pPatch%nBTrisMax  = RFLU_SetMaxDimension(global,pPatch%nBTrisTot)
      pPatch%nBQuadsMax = RFLU_SetMaxDimension(global,pPatch%nBQuadsTot) 
      pPatch%nBFacesMax = RFLU_SetMaxDimension(global,pPatch%nBFacesTot)             
      pPatch%nBVertMax  = RFLU_SetMaxDimension(global,pPatch%nBVertTot)
      
      pPatch%nBCellsVirt = 0                    
    END DO ! iPatch  
    
    pGrid%nBFacesTot = pGrid%nBFaces    
    
! ==============================================================================
!   Allocate memory for boundary-face connectivity
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( pPatch%nBTrisMax > 0 ) THEN 
        ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
        END IF ! global%error 
      ELSE 
        NULLIFY(pPatch%bTri2v)
      END IF ! pPatch%nBTrisMax
      
      IF ( pPatch%nBQuadsMax > 0 ) THEN      
        ALLOCATE(pPatch%bQuad2v(4,pPatch%nBQuadsMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bQuad2v')
        END IF ! global%error                       
      ELSE 
        NULLIFY(pPatch%bQuad2v)
      END IF ! pPatch%nBQuadsMax
    END DO ! iPatch    
    
! ==============================================================================
!   Build boundary face connectivity
! ==============================================================================
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      pPatch%nBTris  = 0
      pPatch%nBQuads = 0              
    END DO ! iPatch    
    
! ------------------------------------------------------------------------------
!   Loop over patches
! ------------------------------------------------------------------------------    
    
    DO iPatch = 1,gridGAMBIT%nPatches
      pPatchGAMBIT => gridGAMBIT%patches(iPatch)

! --- Get mapped new patch number ----------------------------------------------
    
      DO iMap = 1,gridGAMBIT%nMappings
        IF ( iPatch >= gridGAMBIT%patch2bc(1,iMap) .AND. & 
             iPatch <= gridGAMBIT%patch2bc(2,iMap) ) THEN 
          iPatch2 = gridGAMBIT%patch2bc(3,iMap)
        END IF ! iPatch
      END DO ! iMap
            
      pPatch => pRegion%patches(iPatch2)

! --- Loop over faces on patch -------------------------------------------------
    
      DO ifg = 1,pPatchGAMBIT%nFaces
        ict = pPatchGAMBIT%bf2ct(ifg)
        icg = pPatchGAMBIT%bf2cgi(ifg)
        ifl = pPatchGAMBIT%bf2fli(ifg)              
        
! ----- Check that cell types agree (defensive coding)
             
        SELECT CASE ( pGrid%cellGlob2Loc(1,icg) ) 
          CASE ( CELL_TYPE_TET ) 
            IF ( ict /= GAMBIT_NTYPE_TET ) THEN 
              CALL ErrorStop(global,ERR_NTYPE_INVALID,__LINE__)
            END IF ! ict
          CASE ( CELL_TYPE_HEX ) 
            IF ( ict /= GAMBIT_NTYPE_HEX ) THEN 
              CALL ErrorStop(global,ERR_NTYPE_INVALID,__LINE__)
            END IF ! ict          
          CASE ( CELL_TYPE_PRI ) 
            IF ( ict /= GAMBIT_NTYPE_PRI ) THEN 
              CALL ErrorStop(global,ERR_NTYPE_INVALID,__LINE__)
            END IF ! ict          
          CASE ( CELL_TYPE_PYR ) 
            IF ( ict /= GAMBIT_NTYPE_PYR ) THEN 
              CALL ErrorStop(global,ERR_NTYPE_INVALID,__LINE__)
            END IF ! ict          
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pGrid%cellGlob2Loc    
        
! ----- Get local cell index         
        
        icl = pGrid%cellGlob2Loc(2,icg)
             
! ----- Store boundary-face connectivity             
                        
        SELECT CASE ( ict )         
          CASE ( GAMBIT_NTYPE_TET ) 
            pPatch%nBTris = pPatch%nBTris + 1
            
            DO ivl = 1,3                                   
              ifl2 = f2fTetGAMBIT(ifl)  
              ivg = pGrid%tet2v(f2vTet(ivl,ifl2),icl)
                            
              pPatch%bTri2v(ivl,pPatch%nBTris) = ivg
            END DO ! ivl
          CASE ( GAMBIT_NTYPE_HEX ) 
            pPatch%nBQuads = pPatch%nBQuads + 1
            
            DO ivl = 1,4
              ifl2 = f2fHexGAMBIT(ifl)
              ivg = pGrid%hex2v(f2vHex(ivl,ifl2),icl)
              
              pPatch%bQuad2v(ivl,pPatch%nBQuads) = ivg
            END DO ! ivl            
          CASE ( GAMBIT_NTYPE_PRI ) 
            IF ( f2vPriGAMBIT(4,ifl) == VERT_NONE ) THEN 
              pPatch%nBTris = pPatch%nBTris + 1

              DO ivl = 1,3  
                ifl2 = f2fPriGAMBIT(ifl)           
                ivg = pGrid%pri2v(f2vPri(ivl,ifl2),icl)

                pPatch%bTri2v(ivl,pPatch%nBTris) = ivg
              END DO ! ivl              
            ELSE 
              pPatch%nBQuads = pPatch%nBQuads + 1

              DO ivl = 1,4
                ifl2 = f2fPriGAMBIT(ifl)
                ivg = pGrid%pri2v(f2vPri(ivl,ifl2),icl)

                pPatch%bQuad2v(ivl,pPatch%nBQuads) = ivg
              END DO ! ivl             
            END IF ! f2vPriGAMBIT            
          CASE ( GAMBIT_NTYPE_PYR ) 
            IF ( f2vPyrGAMBIT(4,ifl) == VERT_NONE ) THEN 
              pPatch%nBTris = pPatch%nBTris + 1

              DO ivl = 1,3  
                ifl2 = f2fPyrGAMBIT(ifl)           
                ivg = pGrid%pyr2v(f2vPyr(ivl,ifl2),icl)

                pPatch%bTri2v(ivl,pPatch%nBTris) = ivg
              END DO ! ivl              
            ELSE 
              pPatch%nBQuads = pPatch%nBQuads + 1

              DO ivl = 1,4
                ifl2 = f2fPyrGAMBIT(ifl)
                ivg = pGrid%pyr2v(f2vPyr(ivl,ifl2),icl)

                pPatch%bQuad2v(ivl,pPatch%nBQuads) = ivg
              END DO ! ivl             
            END IF ! f2vPyrGAMBIT         
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! ict        
      END DO ! ifg                  
    END DO ! iPatch       
    
    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
        'Converting patch data structure done.'  
    END IF ! global%verbLevel    
    
! ******************************************************************************
!   Dellocate temporary memory 
! ******************************************************************************
    
    DEALLOCATE(pGrid%cellGlob2Loc,STAT=errorFlag)      
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%cellGlob2Loc')
    END IF ! global%error    
    
! ******************************************************************************
!   Allocate memory for other Rocflu data structures 
! ******************************************************************************
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      ALLOCATE(pPatch%bf2c(pPatch%nBFacesMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2c')
      END IF ! global%error

      ALLOCATE(pPatch%bf2v(4,pPatch%nBFacesMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2v')
      END IF ! global%error  

      DO ifl = 1,pPatch%nBFacesMax
        pPatch%bf2v(1,ifl) = VERT_NONE
        pPatch%bf2v(2,ifl) = VERT_NONE 
        pPatch%bf2v(3,ifl) = VERT_NONE 
        pPatch%bf2v(4,ifl) = VERT_NONE                                    
      END DO ! ifl
    END DO ! iPatch    
        
! ******************************************************************************
!   Deallocate GAMBIT memory
! ******************************************************************************

    DEALLOCATE(gridGAMBIT%ct,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridGAMBIT%ct')
    END IF ! global%error
    
    DEALLOCATE(gridGAMBIT%c2v,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridGAMBIT%c2v')
    END IF ! global%error
  
    DEALLOCATE(gridGAMBIT%patches,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridGAMBIT%patches')
    END IF ! global%error        
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Converting from GAMBIT to ROCFLU format done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_ConvGAMBIT2ROCFLU









! *******************************************************************************
!
! Purpose: Print GAMBIT grid information.
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
! *******************************************************************************

  SUBROUTINE RFLU_PrintGridGAMBITInfo(pRegion)

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
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patchGAMBIT), POINTER :: pPatchGAMBIT    

! ******************************************************************************
!   Start, set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Write information
! ******************************************************************************

    WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Grid Statistics:'
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Vertices:       ', & 
                                   pGrid%nVertTot
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Cells:          ', & 
                                   pGrid%nCellsTot                   
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Patches:        ', & 
                                   gridGAMBIT%nPatches   

    IF ( gridGAMBIT%nPatches > 0 ) THEN 
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Patch statistics:'    
    
      DO iPatch = 1,gridGAMBIT%nPatches
        pPatchGAMBIT => gridGAMBIT%patches(iPatch)
        
        WRITE(STDOUT,'(A,7X,A,2X,I4)') SOLVER_NAME,'Patch:',iPatch       
        WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Faces:',pPatchGAMBIT%nFaces
      END DO ! iPatch    
    END IF ! gridGAMBIT%nPatches

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_PrintGridGAMBITInfo










! *******************************************************************************
!
! Purpose: Read grid file from GAMBIT in neutral format.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes:
!   1. CENTAUR cell and node pointers are not read in - read into the 
!      dummy integer idum.
!
! *******************************************************************************

  SUBROUTINE RFLU_ReadGridGAMBITNeutral(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNamePlain
  
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

    CHARACTER(CHRLEN) :: dummyString,iFileName,sectionString,versionString  
    INTEGER :: ITYPE,NENTRY
    INTEGER :: dummyInteger,errorFlag,icg,icl,iFile,ifl,iPatch,ivg,ivl, & 
               loopCounter
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    TYPE(t_patchGAMBIT), POINTER :: pPatchGAMBIT    

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridGAMBITNeutral', &
                          'RFLU_ModGAMBIT.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading GAMBIT neutral grid file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.neu',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

! ******************************************************************************
!   Set grid pointer and initialize variables
! ******************************************************************************

    pGrid => pRegion%grid
    
    pGrid%nVert = 0
    pGrid%nTets = 0        
    pGrid%nHexs = 0        
    pGrid%nPris = 0        
    pGrid%nPyrs = 0     
                  
    pGrid%nPatches = 0

    iPatch = 0

! ******************************************************************************
!   Read header
! ******************************************************************************

    READ(iFile,'(2(A20))') dummyString,versionString
    READ(iFile,*) dummyString    
    READ(iFile,*) gridGAMBIT%title
        
    READ(iFile,*) dummyString 
    READ(iFile,*) dummyString
    READ(iFile,*) dummyString
    
    READ(iFile,*) pGrid%nVert,pGrid%nCells,gridGAMBIT%NGRPS, &
                  gridGAMBIT%nPatches,gridGAMBIT%NDFCD,gridGAMBIT%NDFVL

    READ(iFile,*) dummyString
      
    IF ( TRIM(dummyString) /= "ENDOFSECTION" ) THEN 
      CALL ErrorStop(global,ERR_STRING_INVALID,__LINE__)
    END IF ! TRIM(sectionString)     
  
    pGrid%nVertTot  = pGrid%nVert
    pGrid%nCellsTot = pGrid%nCells

    pGrid%nVertMax  = RFLU_SetMaxDimension(global,pGrid%nVertTot)
    pGrid%nCellsMax = RFLU_SetMaxDimension(global,pGrid%nCellsTot)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pGrid%xyz(XCOORD:ZCOORD,pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grid%xyz')
    END IF ! global%error

    ALLOCATE(gridGAMBIT%ct(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridGAMBIT%ct')
    END IF ! global%error
    
    ALLOCATE(gridGAMBIT%c2v(8,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridGAMBIT%c2v')
    END IF ! global%error
  
    DO icg = 1,pGrid%nCellsTot        
      gridGAMBIT%c2v(1,icg) = C2V_INIT
      gridGAMBIT%c2v(2,icg) = C2V_INIT
      gridGAMBIT%c2v(3,icg) = C2V_INIT
      gridGAMBIT%c2v(4,icg) = C2V_INIT
      gridGAMBIT%c2v(5,icg) = C2V_INIT
      gridGAMBIT%c2v(6,icg) = C2V_INIT
      gridGAMBIT%c2v(7,icg) = C2V_INIT
      gridGAMBIT%c2v(8,icg) = C2V_INIT            
    END DO ! icg
  
    ALLOCATE(gridGAMBIT%patches(gridGAMBIT%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridGAMBIT%patches')
    END IF ! global%error    
  
! ******************************************************************************
!   Read sections
! ******************************************************************************
  
    loopCounter = 0
  
    DO ! set up infinite loop
      loopCounter = loopCounter + 1

! ==============================================================================
!     Read section string and take appropriate action
! ==============================================================================
  
      READ(iFile,'(A32)',IOSTAT=errorFlag,END=100) sectionString

! ==============================================================================
!     Coordinates
! ==============================================================================

      IF ( ADJUSTL(TRIM(sectionString)) == & 
           "NODAL COORDINATES"//" "//ADJUSTL(TRIM(versionString)) ) THEN           
        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading coordinate section...'
        END IF ! global%verbLevel           
             
        DO ivl = 1,pGrid%nVertTot
          READ(iFile,*) ivg
          BACKSPACE(iFile)
            
          READ(iFile,*) dummyInteger,pGrid%xyz(XCOORD:ZCOORD,ivg)
        END DO ! ivl

! ==============================================================================
!     Element connectivity
! ==============================================================================
      
      ELSE IF ( ADJUSTL(TRIM(sectionString)) == & 
                "ELEMENTS/CELLS"//" "//ADJUSTL(TRIM(versionString)) ) THEN 
        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading element connectivity...'
        END IF ! global%verbLevel  

        DO icl = 1,pGrid%nCellsTot
          READ(iFile,*) icg
          BACKSPACE(iFile)
            
          READ(iFile,*) dummyString,gridGAMBIT%NTYPE,gridGAMBIT%NDP
          BACKSPACE(iFile)
            
          gridGAMBIT%ct(icg) = gridGAMBIT%NTYPE

! ------------------------------------------------------------------------------
!         Read element connectivity
! ------------------------------------------------------------------------------
          
          SELECT CASE ( gridGAMBIT%NTYPE ) 
            
! --------- Edge ---------------------------------------------------------------   
            
            CASE ( GAMBIT_NTYPE_EDGE ) 
              CALL ErrorStop(global,ERR_NTYPE_INVALID,__LINE__)               

! --------- Quadrilateral ------------------------------------------------------        

            CASE ( GAMBIT_NTYPE_QUAD )
              CALL ErrorStop(global,ERR_NTYPE_INVALID,__LINE__)               

! --------- Triangle -----------------------------------------------------------      

            CASE ( GAMBIT_NTYPE_TRI )
              CALL ErrorStop(global,ERR_NTYPE_INVALID,__LINE__)

! --------- Hexahedron ---------------------------------------------------------           

            CASE ( GAMBIT_NTYPE_HEX ) 
              IF ( gridGAMBIT%NDP == 8 ) THEN
                READ(iFile,*) dummyString,dummyString,dummyString, & 
                              gridGAMBIT%c2v(1:7,icg)
                        
                IF ( gridGAMBIT%NDP == 8 ) THEN
                  pGrid%nHexs = pGrid%nHexs + 1
                 
                  READ(iFile,*) gridGAMBIT%c2v(8,icg)              
                END IF ! iFile                                
              ELSE 
                CALL ErrorStop(global,ERR_NDP_INVALID,__LINE__)
              END IF ! gridGAMBIT%NDP

! --------- Prism --------------------------------------------------------------            

            CASE ( GAMBIT_NTYPE_PRI ) 
              IF ( gridGAMBIT%NDP == 6 ) THEN
                pGrid%nPris = pGrid%nPris + 1              
              
                READ(iFile,*) dummyString,dummyString,dummyString, & 
                              gridGAMBIT%c2v(1:6,icg)                
              ELSE 
                CALL ErrorStop(global,ERR_NDP_INVALID,__LINE__)
              END IF ! gridGAMBIT%NDP              

! --------- Tetrahedron --------------------------------------------------------             

            CASE ( GAMBIT_NTYPE_TET ) 
              IF ( gridGAMBIT%NDP == 4 ) THEN
                pGrid%nTets = pGrid%nTets + 1              
              
                READ(iFile,*) dummyString,dummyString,dummyString, & 
                              gridGAMBIT%c2v(1:4,icg)                
              ELSE                  
                CALL ErrorStop(global,ERR_NDP_INVALID,__LINE__)
              END IF ! gridGAMBIT%NDP                   

! --------- Pyramid ------------------------------------------------------------             

            CASE ( GAMBIT_NTYPE_PYR ) 
              IF ( gridGAMBIT%NDP == 5 ) THEN
                pGrid%nPyrs = pGrid%nPyrs + 1              
              
                READ(iFile,*) dummyString,dummyString,dummyString, & 
                              gridGAMBIT%c2v(1:5,icg)                
              ELSE                  
                CALL ErrorStop(global,ERR_NDP_INVALID,__LINE__)
              END IF ! gridGAMBIT%NDP
              
! --------- Unknown ------------------------------------------------------------              
                            
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)   
          END SELECT ! gridGAMBIT%NTYPE                                 
        END DO ! icl

! ==============================================================================
!     Element group
! ==============================================================================

      ELSE IF ( ADJUSTL(TRIM(sectionString)) == & 
                "ELEMENT GROUP"//" "//ADJUSTL(TRIM(versionString)) ) THEN        
        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
            'Reading element group information...'
        END IF ! global%verbLevel 

        READ(iFile,*) dummyString,gridGAMBIT%NGP,dummyString,gridGAMBIT%NELGP, & 
                      dummyString,gridGAMBIT%MTYP,dummyString,gridGAMBIT%NFLAGS 

        READ(iFile,*) dummyString          

        READ(iFile,'(10(I8))') (dummyInteger,ifl=1,gridGAMBIT%NFLAGS)
        READ(iFile,'(10(I8))') (dummyInteger,icl=1,gridGAMBIT%NELGP)

! ==============================================================================
!     Boundary information
! ==============================================================================

      ELSE IF ( ADJUSTL(TRIM(sectionString)) == & 
                "BOUNDARY CONDITIONS"//" "//ADJUSTL(TRIM(versionString)) ) THEN 
        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
            'Reading boundary condition information...'
        END IF ! global%verbLevel 

        iPatch = iPatch + 1

        pPatchGAMBIT => gridGAMBIT%patches(iPatch)   

        READ(iFile,'(A32,2(I10))') dummyString,ITYPE,NENTRY

! ------------------------------------------------------------------------------
!       Read face or vertex data. NOTE in either case, values are not read. 
! ------------------------------------------------------------------------------

        SELECT CASE ( ITYPE )

! ------- Face data ------------------------------------------------------------  

          CASE ( 1 ) ! Face data              
            pPatchGAMBIT%nFaces = NENTRY

            ALLOCATE(pPatchGAMBIT%bf2cgi(pPatchGAMBIT%nFaces),STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN 
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                             'pPatchGAMBIT%bf2cgi')
            END IF ! global%error

            ALLOCATE(pPatchGAMBIT%bf2ct(pPatchGAMBIT%nFaces),STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN 
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                             'pPatchGAMBIT%bf2ct')
            END IF ! global%error

            ALLOCATE(pPatchGAMBIT%bf2fli(pPatchGAMBIT%nFaces),STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN 
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                             'pPatchGAMBIT%bf2fli')
            END IF ! global%error              

            DO ifl = 1,pPatchGAMBIT%nFaces
              READ(iFile,*) pPatchGAMBIT%bf2cgi(ifl), & 
                            pPatchGAMBIT%bf2ct(ifl), & 
                            pPatchGAMBIT%bf2fli(ifl)                              
            END DO ! ifl

! ------- Vertex data ----------------------------------------------------------

          CASE ( 0 ) ! Vertex data
            DO ivl = 1,NENTRY
              READ(iFile,*) dummyInteger
            END DO ! ivl

! ------- Unknown --------------------------------------------------------------             

          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)          
        END SELECT ! ITYPE 
          
! ==============================================================================
!     Unknown section string
! ==============================================================================
      
      ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
      END IF ! TRIM(sectionString)

! ==============================================================================
!     Read end of section string 
! ==============================================================================

      READ(iFile,*) dummyString
      
      IF ( TRIM(dummyString) /= "ENDOFSECTION" ) THEN 
        CALL ErrorStop(global,ERR_STRING_INVALID,__LINE__)
      END IF ! TRIM(sectionString)

! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>

! ******************************************************************************
!   EOF condition. NOTE assume that this is ok, because GAMBIT grid files do not
!   have a dedicated EOF marker, and hence need to use READ statement to detect
!   EOF condition.
! ******************************************************************************

    global%warnCounter = global%warnCounter + 1

100 WRITE(STDOUT,*) SOLVER_NAME,'*** WARNING *** Encountered EOF.'

! ******************************************************************************
!   Set sizes
! ******************************************************************************

    pGrid%nTetsTot = pGrid%nTets       
    pGrid%nHexsTot = pGrid%nHexs         
    pGrid%nPrisTot = pGrid%nPris          
    pGrid%nPyrsTot = pGrid%nPyrs       

    pGrid%nTetsMax = RFLU_SetMaxDimension(global,pGrid%nTetsTot)
    pGrid%nHexsMax = RFLU_SetMaxDimension(global,pGrid%nHexsTot)
    pGrid%nPrisMax = RFLU_SetMaxDimension(global,pGrid%nPrisTot)    
    pGrid%nPyrsMax = RFLU_SetMaxDimension(global,pGrid%nPyrsTot)    
    
! ******************************************************************************
!   Check validity of connectivity arrays
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN
      CALL RFLU_CheckGridGAMBIT(pRegion)
    END IF ! global%checkLevel

! ******************************************************************************
!   Print grid statistics
! ******************************************************************************
 
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      CALL RFLU_PrintGridGAMBITInfo(pRegion) 
    END IF ! global%verbLevel   

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

    IF ( global%verbLevel > VERBOSE_NONE ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading GAMBIT neutral grid file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
   
  END SUBROUTINE RFLU_ReadGridGAMBITNeutral






! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModGAMBIT

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGAMBIT.F90,v $
! Revision 1.4  2008/12/06 08:45:03  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/25 22:04:29  haselbac
! Changes because of sype patches
!
! Revision 1.1  2005/04/15 15:09:09  haselbac
! Initial revision
!
! Revision 1.3  2005/01/20 14:54:56  haselbac
! Added setting of nBFaces and nBFacesTot
!
! Revision 1.2  2004/11/03 17:12:15  haselbac
! Removed setting of vertex and cell flags
!
! Revision 1.1  2004/11/03 15:06:57  haselbac
! Initial revision
!
! ******************************************************************************









