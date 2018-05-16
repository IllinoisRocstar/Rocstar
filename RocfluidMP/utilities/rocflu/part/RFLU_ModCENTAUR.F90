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
! Purpose: Collection of routines to read and convert CENTAUR grids.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModCENTAUR.F90,v 1.4 2008/12/06 08:45:03 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModCENTAUR

  USE ModParameters
  USE ModDataTypes  
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModDimensions, ONLY: RFLU_SetMaxDimension

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  PRIVATE

! ==============================================================================
! Data
! ==============================================================================

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModCENTAUR.F90,v $ $Revision: 1.4 $'        

  TYPE t_gridCENTAUR 
    INTEGER :: nBQuads,nBTris
    INTEGER, DIMENSION(:,:), POINTER :: bInfo,bTri2v,bQuad2v  
    CHARACTER(CHRLEN) :: title
    CHARACTER(CHRLEN), DIMENSION(:), POINTER :: bName    
  END TYPE t_gridCENTAUR

  TYPE(t_gridCENTAUR) :: gridCENTAUR

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_ConvCENTAUR2ROCFLU, & 
            RFLU_ReadGridCENTAURASCII, & 
            RFLU_ReadGridCENTAURBinary

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

  SUBROUTINE RFLU_CheckGridCENTAUR(pRegion)
     
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

    INTEGER :: errorFlag,cvmax,cvmin
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CheckGridCENTAUR', &
                          'RFLU_ModCENTAUR.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Checking connectivity arrays...'    
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Volume grid
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                               'Volume grid...'    
    END IF ! global%verbLevel

! ==============================================================================
!   Tetrahedra
! ==============================================================================

    IF ( pGrid%nTetsTot > 0 ) THEN 
      cvmin = MINVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))
      cvmax = MAXVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))

      IF ( pGrid%nTetsTot == pGrid%nCellsTot ) THEN       
        IF ( cvmin /= 1 .OR. cvmax /= pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! cvmin
      ELSE  
        IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! vmin
      END IF ! cvmin

      IF ( global%error /= ERR_NONE ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN          
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel           
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END IF ! pGrid

! ==============================================================================
!   Hexahedra
! ==============================================================================

    IF ( pGrid%nHexsTot > 0 ) THEN 
      cvmin = MINVAL(pGrid%hex2v(1:8,1:pGrid%nHexsTot))
      cvmax = MAXVAL(pGrid%hex2v(1:8,1:pGrid%nHexsTot))

      IF ( pGrid%nHexsTot == pGrid%nCellsTot ) THEN       
        IF ( cvmin /= 1 .OR. cvmax /= pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! cvmin
      ELSE  
        IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! vmin
      END IF ! cvmin

      IF ( global%error /= ERR_NONE ) THEN 
        IF ( global%verbLevel > VERBOSE_NONE ) THEN     
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel     
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END IF ! pGrid    

! ==============================================================================
!   Prisms
! ==============================================================================

    IF ( pGrid%nPrisTot > 0 ) THEN 
      cvmin = MINVAL(pGrid%pri2v(1:6,1:pGrid%nPrisTot))
      cvmax = MAXVAL(pGrid%pri2v(1:6,1:pGrid%nPrisTot))

      IF ( pGrid%nPrisTot == pGrid%nCellsTot ) THEN       
        IF ( cvmin /= 1 .OR. cvmax /= pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! cvmin
      ELSE  
        IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! vmin
      END IF ! cvmin

      IF ( global%error /= ERR_NONE ) THEN   
        IF ( global%verbLevel > VERBOSE_NONE ) THEN         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel              
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END IF ! pGrid     

! ==============================================================================
!   Pyramids
! ==============================================================================

    IF ( pGrid%nPyrsTot > 0 ) THEN 
      cvmin = MINVAL(pGrid%pyr2v(1:5,1:pGrid%nPyrsTot))
      cvmax = MAXVAL(pGrid%pyr2v(1:5,1:pGrid%nPyrsTot))

      IF ( pGrid%nPyrsTot == pGrid%nCellsTot ) THEN       
        IF ( cvmin /= 1 .OR. cvmax /= pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! cvmin
      ELSE  
        IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! vmin
      END IF ! cvmin

      IF ( global%error /= ERR_NONE ) THEN 
        IF ( global%verbLevel > VERBOSE_NONE ) THEN         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel   
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END IF ! pGrid

! ******************************************************************************
!   Surface grid
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN     
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Surface grid...'    
    END IF ! global%verbLevel

    IF ( gridCENTAUR%nBTris > 0 ) THEN 
      cvmin = MINVAL(gridCENTAUR%bTri2v(:,:))
      cvmax = MAXVAL(gridCENTAUR%bTri2v(:,:))

      IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN 
        global%error = ERR_VERTEX_NUMBER
      END IF ! cvmin

      IF ( global%error /= ERR_NONE ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN  
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel 
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END IF ! gridCENTAUR

    IF ( gridCENTAUR%nBQuads > 0 ) THEN 
      cvmin = MINVAL(gridCENTAUR%bQuad2v(:,:))
      cvmax = MAXVAL(gridCENTAUR%bQuad2v(:,:))

      IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN 
        global%error = ERR_VERTEX_NUMBER
      END IF ! cvmin

      IF ( global%error /= ERR_NONE ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN  
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel 
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error      
    END IF ! gridCENTAUR

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Checking connectivity arrays done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_CheckGridCENTAUR





! ******************************************************************************
!
! Purpose: Convert grid format from CENTAUR to ROCFLU.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. bcName from CENTAUR file is discarded. It is NOT transferred to the 
!      Rocflu data structure, instead, the name is read from the bc file.
!   2. IMPORTANT: If empty patches are detected (the case for multi-zone grids
!      when hybconvert is instructed to ignore the interzone boundary), these
!      are deleted automatically. This means that the number of patches changes
!      and hence the mapping of boundary conditions to boundary patches changes
!      also. Therefore if empty patches were eliminated, the user will have to 
!      edit the boundary-condition file accordingly.
!
! ******************************************************************************

  SUBROUTINE RFLU_ConvCENTAUR2ROCFLU(pRegion)
     
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

    INTEGER :: errorFlag,ifl,iPatch,iPatch2,iqbeg,iqend,itbeg,itend,nBQuads, &
               nBTris,nPatchesOld
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvCENTAUR2ROCFLU', &
                          'RFLU_ModCENTAUR.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from CENTAUR to ROCFLU format...'
    END IF ! global%verbLevel

! ==============================================================================
!   Set grid pointer and initialize variables
! ==============================================================================

    pGrid => pRegion%grid

    pGrid%nEdges    = 0
    pGrid%nEdgesTot = 0

    pGrid%nFaces    = 0
    pGrid%nFacesTot = 0

! ==============================================================================
!   Check for and eliminate empty patches. NOTE: boundary names are not copied 
!   and the user must adjust the boundary-condition file if empty patches were 
!   indeed eliminated...
! ==============================================================================

    nPatchesOld = pGrid%nPatches

    iPatch = 0

    IF ( pGrid%nPatches > 0 ) THEN 
      iPatch = iPatch + 1

      DO iPatch2 = 1,nPatchesOld
        IF ( iPatch2 /= 1 ) THEN 
          itbeg = gridCENTAUR%bInfo(2,iPatch2-1) + 1
          iqbeg = gridCENTAUR%bInfo(3,iPatch2-1) + 1
        ELSE 
          itbeg = 1
          iqbeg = 1
        END IF ! iPatch

        itend = gridCENTAUR%bInfo(2,iPatch2)
        iqend = gridCENTAUR%bInfo(3,iPatch2)

        nBTris  = itend - itbeg + 1
        nBQuads = iqend - iqbeg + 1             

        IF ( (nBTris + nBQuads) /= 0 ) THEN
          gridCENTAUR%bInfo(1,iPatch) = gridCENTAUR%bInfo(1,iPatch2)
          gridCENTAUR%bInfo(2,iPatch) = gridCENTAUR%bInfo(2,iPatch2)
          gridCENTAUR%bInfo(3,iPatch) = gridCENTAUR%bInfo(3,iPatch2)

          iPatch = iPatch + 1      
        ELSE 
          WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A)') SOLVER_NAME, & 
            '*** WARNING *** Patch',iPatch2,'is empty and will be deleted!' 

          pGrid%nPatches = pGrid%nPatches - 1     
        END IF ! nBTris    
      END DO ! iPatch2
    END IF ! pGrid%nPatches  
    
! ==============================================================================
!   Convert patch data structure
! ==============================================================================

    global%nPatches = pGrid%nPatches

    IF ( pGrid%nPatches > 0 ) THEN 
      ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
      END IF ! global%error

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        pPatch%bcType = gridCENTAUR%bInfo(1,iPatch)

        pPatch%iPatchGlobal = iPatch
        pPatch%iBorder      = PATCH_IBORDER_DEFAULT
        pPatch%renumFlag    = .FALSE.
                   
        IF ( iPatch /= 1 ) THEN 
          itbeg = gridCENTAUR%bInfo(2,iPatch-1) + 1
          iqbeg = gridCENTAUR%bInfo(3,iPatch-1) + 1
        ELSE 
          itbeg = 1
          iqbeg = 1
        END IF ! iPatch

        itend = gridCENTAUR%bInfo(2,iPatch)
        iqend = gridCENTAUR%bInfo(3,iPatch)

        pPatch%nBTrisTot  = itend - itbeg + 1
        pPatch%nBQuadsTot = iqend - iqbeg + 1
        pPatch%nBVertTot  = 0

        pPatch%nBTris  = pPatch%nBTrisTot  
        pPatch%nBQuads = pPatch%nBQuadsTot 
        pPatch%nBVert  = pPatch%nBVertTot
                
        pPatch%nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot 
        pPatch%nBFaces    = pPatch%nBFacesTot
          
        pPatch%nBTrisMax  = RFLU_SetMaxDimension(global,pPatch%nBTrisTot)
        pPatch%nBQuadsMax = RFLU_SetMaxDimension(global,pPatch%nBQuadsTot)        
        pPatch%nBVertMax  = RFLU_SetMaxDimension(global,pPatch%nBVertTot)
   
        pPatch%nBFacesMax = RFLU_SetMaxDimension(global,pPatch%nBFacesTot)

        pPatch%nBCellsVirt = 0  

        IF ( pPatch%nBTrisMax > 0 ) THEN 
          ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
          END IF ! global%error
        ELSE 
          NULLIFY(pPatch%bTri2v)
        END IF ! pPatch%nBTrisMax

        IF ( pPatch%nBTrisTot > 0 ) THEN 
          pPatch%bTri2v(1:3,1:pPatch%nBTrisTot) = & 
            gridCENTAUR%bTri2v(1:3,itbeg:itend)                   
        END IF ! pPatch

        IF ( pPatch%nBQuadsMax > 0 ) THEN 
          ALLOCATE(pPatch%bQuad2v(4,pPatch%nBQuadsMax),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bQuad2v')
          END IF ! global%error
        ELSE
          NULLIFY(pPatch%bQuad2v)
        END IF ! pPatch%nBQuadsMax
        
        IF ( pPatch%nBQuadsTot > 0 ) THEN 
          pPatch%bQuad2v(1:4,1:pPatch%nBQuadsTot) = & 
            gridCENTAUR%bQuad2v(1:4,iqbeg:iqend)  
        END IF ! pPatch 
      END DO ! iPatch
    END IF ! pGrid%nPatches

! ==============================================================================
!   Compute number of faces on each patch and set patch quantities
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      pPatch%nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot 

      pPatch%nBFaces = pPatch%nBFacesTot
      pPatch%nBQuads = pPatch%nBQuadsTot  
      pPatch%nBTris  = pPatch%nBTrisTot   
      pPatch%nBVert  = pPatch%nBVertTot    

      pPatch%nBFacesMax = RFLU_SetMaxDimension(global,pPatch%nBFacesTot)
    END DO ! iPatch   

! ==============================================================================
!   Check that number of triangular and quadrilateral faces correct and set 
!   number of boundary faces
! ==============================================================================

    nBTris  = 0
    nBQuads = 0

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      nBTris  = nBTris  + pPatch%nBTris
      nBQuads = nBQuads + pPatch%nBQuads
    END DO ! iPatch 

    IF ( nBTris  /= gridCENTAUR%nBTris .OR. & 
         nBQuads /= gridCENTAUR%nBQuads ) THEN 
      CALL ErrorStop(global,ERR_NBFACES_WRONG,__LINE__)
    END IF ! nBTris

    pGrid%nBFaces    = nBTris + nBQuads
    pGrid%nBFacesTot = pGrid%nBFaces    

! ******************************************************************************
!   Deallocate CENTAUR memory
! ******************************************************************************

    IF ( pGrid%nPatches > 0 ) THEN 
      IF ( gridCENTAUR%nBTris > 0 ) THEN 
        DEALLOCATE(gridCENTAUR%bTri2v,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridCENTAUR%bTri2v')
        END IF ! global%error  
      END IF ! gridCENTAUR

      IF ( gridCENTAUR%nBQuads > 0 ) THEN 
        DEALLOCATE(gridCENTAUR%bQuad2v,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridCENTAUR%bQuad2v')
        END IF ! global%error 
      END IF ! gridCENTAUR

      DEALLOCATE(gridCENTAUR%bInfo,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridCENTAUR%bInfo')
      END IF ! global%error

      DEALLOCATE(gridCENTAUR%bName,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridCENTAUR%bName')
      END IF ! global%error
    END IF ! pGrid%nPatches
  
! *****************************************************************************
!   Allocate memory for boundary face lists bf2c and bf2v
! *****************************************************************************

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
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Converting from CENTAUR to ROCFLU format done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_ConvCENTAUR2ROCFLU







! *******************************************************************************
!
! Purpose: Print CENTAUR grid information.
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

  SUBROUTINE RFLU_PrintGridCENTAURInfo(pRegion)

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

    TYPE(t_grid), POINTER :: pGrid

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
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Tetrahedra:     ', & 
                                   pGrid%nTetsTot
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Hexahedra:      ', & 
                                   pGrid%nHexsTot
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Prisms:         ', & 
                                   pGrid%nPrisTot
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Pyramids:       ', & 
                                   pGrid%nPyrsTot
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Patches:        ', & 
                                   pGrid%nPatches   

    IF ( pGrid%nPatches > 0 ) THEN 
      WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Patch faces:    ', & 
            gridCENTAUR%bInfo(2,pGrid%nPatches) &
          + gridCENTAUR%bInfo(3,pGrid%nPatches)
      WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Triangles:      ', & 
            gridCENTAUR%bInfo(2,pGrid%nPatches)
      WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Quadrilaterals: ', & 
            gridCENTAUR%bInfo(3,pGrid%nPatches)
    END IF ! pGrid%nPatches

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_PrintGridCENTAURInfo








! *******************************************************************************
!
! Purpose: Read grid file from CENTAUR in ASCII format.
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

  SUBROUTINE RFLU_ReadGridCENTAURASCII(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNamePlain

    USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_SetSyPePatchesFlag

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
    INTEGER :: dummyInteger,errorFlag,i,iFile,j
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file and read title
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridCENTAURASCII', &
                          'RFLU_ModCENTAUR.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII CENTAUR grid file...'
    END IF ! global%verbLevel

    iFile  = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.hyb.asc',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

    READ(iFile,'(A80)') gridCENTAUR%title
  
! ==============================================================================
!   Coordinates
! ==============================================================================

    pGrid => pRegion%grid

    READ(iFile,'(I8)') pGrid%nVertTot

    pGrid%nVertMax = RFLU_SetMaxDimension(global,pGrid%nVertTot)

    ALLOCATE(pGrid%xyz(3,pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grid%xyz')
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
    END IF ! global%verbLevel

    DO i = 1,3 
      READ(iFile,'(5E16.9)') (pGrid%xyz(i,j),j=1,pGrid%nVertTot)  
    END DO ! i

    READ(iFile,'(10I8)') (dummyInteger,i=1,pGrid%nVertTot)

! ==============================================================================
!   Cell connectivity
! ==============================================================================

    READ(iFile,'(I8)') pGrid%nTetsTot

    pGrid%nTetsMax = RFLU_SetMaxDimension(global,pGrid%nTetsTot)

    IF ( pGrid%nTetsMax > 0 ) THEN 
      ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%tet2v)
    END IF ! pGrid%nTetsMax

    IF ( pGrid%nTetsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel 

      DO i = 1,4         
        READ(iFile,'(10I8)') (pGrid%tet2v(i,j),j=1,pGrid%nTetsTot)
      END DO ! i

      READ(iFile,'(10I8)') (dummyInteger,i=1,pGrid%nTetsTot) 
    END IF ! pGrid%nTetsTot


    READ(iFile,'(I8)') pGrid%nHexsTot

    pGrid%nHexsMax = RFLU_SetMaxDimension(global,pGrid%nHexsTot)

    IF ( pGrid%nHexsMax > 0 ) THEN 
      ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%hex2v)
    END IF ! pGrid%nHexsMax

    IF ( pGrid%nHexsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...' 
      END IF ! global%verbLevel

      DO i = 1,8
        READ(iFile,'(10I8)') (pGrid%hex2v(i,j),j=1,pGrid%nHexsTot)
      END DO ! i

      READ(iFile,'(10I8)') (dummyInteger,i=1,pGrid%nHexsTot)
    END IF ! pGrid%nHexsTot


    READ(iFile,'(I8)') pGrid%nPrisTot

    pGrid%nPrisMax = RFLU_SetMaxDimension(global,pGrid%nPrisTot)

    IF ( pGrid%nPrisMax > 0 ) THEN 
      ALLOCATE(pGrid%pri2v(6,pGrid%nPrisMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pri2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%pri2v)
    END IF ! pGrid%nPrisMax
    
    IF ( pGrid%nPrisTot > 0 ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Prisms...'    
      END IF ! global%verbLevel

      DO i = 1,6
        READ(iFile,'(10I8)') (pGrid%pri2v(i,j),j=1,pGrid%nPrisTot)
      END DO ! i

      READ(iFile,'(10I8)') (dummyInteger,i=1,pGrid%nPrisTot)    
    END IF ! pGrid%nPrisTot


    READ(iFile,'(I8)') pGrid%nPyrsTot

    pGrid%nPyrsMax = RFLU_SetMaxDimension(global,pGrid%nPyrsTot)

    IF ( pGrid%nPyrsMax > 0 ) THEN 
      ALLOCATE(pGrid%pyr2v(5,pGrid%nPyrsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pyr2v')
      END IF ! global%error
    ELSE
      NULLIFY(pGrid%pyr2v)
    END IF ! pGrid%nPyrsMax

    IF ( pGrid%nPyrsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pyramids...'    
      END IF ! global%verbLevel

      DO i = 1,5
        READ(iFile,'(10I8)') (pGrid%pyr2v(i,j),j=1,pGrid%nPyrsTot)
      END DO ! i

      READ(iFile,'(10I8)') (dummyInteger,i=1,pGrid%nPyrsTot)
    END IF ! pGrid%nPyrsTot

! ******************************************************************************
!   Set grid size variables
! ******************************************************************************

    pGrid%nCellsTot = pGrid%nTetsTot + pGrid%nHexsTot + pGrid%nPrisTot & 
                    + pGrid%nPyrsTot
    pGrid%nCellsMax = RFLU_SetMaxDimension(global,pGrid%nCellsTot)

    pGrid%nVert  = pGrid%nVertTot
    pGrid%nCells = pGrid%nCellsTot
    pGrid%nTets  = pGrid%nTetsTot
    pGrid%nHexs  = pGrid%nHexsTot
    pGrid%nPris  = pGrid%nPrisTot
    pGrid%nPyrs  = pGrid%nPyrsTot

! ==============================================================================
!   Boundary types
! ==============================================================================
  
    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
    END IF ! global%verbLevel

    READ(ifile,'(I8)') pGrid%nPatches  

    IF ( pGrid%nPatches > 0 ) THEN 
      ALLOCATE(gridCENTAUR%bInfo(3,pGrid%nPatches),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bInfo')
      END IF ! global%error

      DO i = 1,3
        READ(ifile,'(10I8)') (gridCENTAUR%bInfo(i,j),j=1,pGrid%nPatches)
      END DO ! i

      ALLOCATE(gridCENTAUR%bName(pGrid%nPatches),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bName')
      END IF ! global%error

      READ(ifile,'(A80)') (gridCENTAUR%bName(i),i=1,pGrid%nPatches) 
    END IF ! pGrid%nPatches

! ==============================================================================
!   Boundary face connectivity
! ==============================================================================

    IF ( pGrid%nPatches > 0 ) THEN 
      READ(iFile,'(I8)') gridCENTAUR%nBTris

      IF ( gridCENTAUR%nBTris > 0 ) THEN 
        ALLOCATE(gridCENTAUR%bTri2v(3,gridCENTAUR%nBTris),STAT=errorFlag)
        global%error = errorFlag 
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bTri2v')
        END IF ! global%error

        IF ( global%verbLevel > VERBOSE_NONE ) THEN    
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary triangles...'    
        END IF ! global%verbLevel

        DO i = 1,3
          READ(iFile,'(10I8)') (gridCENTAUR%bTri2v(i,j), & 
                                j=1,gridCENTAUR%nBTris)
        END DO ! i
      END IF ! gridCENTAUR

      READ(iFile,'(I8)') gridCENTAUR%nBQuads

      IF ( gridCENTAUR%nBQuads > 0 ) THEN 
        ALLOCATE(gridCENTAUR%bQuad2v(4,gridCENTAUR%nBQuads),STAT=errorFlag)
        global%error = errorFlag 
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bQuad2v')
        END IF ! global%error

        IF ( global%verbLevel > VERBOSE_NONE ) THEN     
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary quadrilaterals...'    
        END IF ! global%verbLevel

        DO i = 1,4
          READ(iFile,'(10I8)') (gridCENTAUR%bQuad2v(i,j), & 
                                j=1,gridCENTAUR%nBQuads)
        END DO ! i
      END IF ! gridCENTAUR
    ELSE 
      gridCENTAUR%nBTris  = 0
      gridCENTAUR%nBQuads = 0  
    END IF ! pGrid%nPatches

! ******************************************************************************
!   Check validity of connectivity arrays
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN
      CALL RFLU_CheckGridCENTAUR(pRegion)
    END IF ! global%checkLevel

! ******************************************************************************
!   Print grid statistics
! ******************************************************************************
 
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      CALL RFLU_PrintGridCENTAURInfo(pRegion) 
    END IF ! global%verbLevel  

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile, IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reading ASCII CENTAUR grid file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
   
  END SUBROUTINE RFLU_ReadGridCENTAURASCII







! *******************************************************************************
!
! Purpose: Read grid file from CENTAUR in binary format.
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

  SUBROUTINE RFLU_ReadGridCENTAURBinary(pRegion)

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

    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: dummyInteger,errorFlag,i,iFile,j
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file and read title
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridCENTAURBinary', &
                          'RFLU_ModCENTAUR.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary CENTAUR grid file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.hyb.bin',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

    READ(iFile) gridCENTAUR%title
  
! ==============================================================================
!   Coordinates
! ==============================================================================

    pGrid => pRegion%grid

    READ(iFile) pGrid%nVertTot

    pGrid%nVertMax = RFLU_SetMaxDimension(global,pGrid%nVertTot)

    ALLOCATE(pGrid%xyz(3,pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grid%xyz')
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
    END IF ! global%verbLevel

    READ(iFile) ((pGrid%xyz(i,j),j=1,pGrid%nVertTot),i=1,3)  
    READ(iFile) (dummyInteger,i=1,pGrid%nVertTot)

! ==============================================================================
!   Cell connectivity
! ==============================================================================

    READ(iFile) pGrid%nTetsTot

    pGrid%nTetsMax = RFLU_SetMaxDimension(global,pGrid%nTetsTot)

    IF ( pGrid%nTetsMax > 0 ) THEN 
      ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
      END IF ! global%error
    ELSE
      NULLIFY(pGrid%tet2v)
    END IF ! pGrid%nTetsMax

    IF ( pGrid%nTetsTot > 0 ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel          

      READ(iFile) ((pGrid%tet2v(i,j),j=1,pGrid%nTetsTot),i=1,4)

      READ(iFile) (dummyInteger,i=1,pGrid%nTetsTot) 
    END IF ! pGrid%nTetsTot


    READ(iFile) pGrid%nHexsTot

    pGrid%nHexsMax = RFLU_SetMaxDimension(global,pGrid%nHexsTot)

    IF ( pGrid%nHexsMax > 0 ) THEN 
      ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%hex2v)
    END IF ! pGrid%nHexsMax

    IF ( pGrid%nHexsTot > 0 ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...' 
      END IF ! global%verbLevel

      READ(iFile) ((pGrid%hex2v(i,j),j=1,pGrid%nHexsTot),i=1,8)

      READ(iFile) (dummyInteger,i=1,pGrid%nHexsTot)   
    END IF ! pGrid%nHexsTot


    READ(iFile) pGrid%nPrisTot

    pGrid%nPrisMax = RFLU_SetMaxDimension(global,pGrid%nPrisTot)

    IF ( pGrid%nPrisMax > 0 ) THEN 
      ALLOCATE(pGrid%pri2v(6,pGrid%nPrisMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pri2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%pri2v)
    END IF ! pGrid%nPrisMax

    IF ( pGrid%nPrisTot > 0 ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Prisms...'    
      END IF ! global%verbLevel

      READ(iFile) ((pGrid%pri2v(i,j),j=1,pGrid%nPrisTot),i=1,6)

      READ(iFile) (dummyInteger,i=1,pGrid%nPrisTot)    
    END IF ! pGrid%nPrisTot


    READ(iFile) pGrid%nPyrsTot

    pGrid%nPyrsMax = RFLU_SetMaxDimension(global,pGrid%nPyrsTot)

    IF ( pGrid%nPyrsMax > 0 ) THEN 
      ALLOCATE(pGrid%pyr2v(5,pGrid%nPyrsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pyr2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%pyr2v)
    END IF ! pGrid%nPyrsMax

    IF ( pGrid%nPyrsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pyramids...'    
      END IF ! global%verbLevel

      READ(iFile) ((pGrid%pyr2v(i,j),j=1,pGrid%nPyrsTot),i=1,5)

      READ(iFile) (dummyInteger,i=1,pGrid%nPyrsTot)
    END IF ! pGrid%nPyrsTot

! ******************************************************************************
!   Set grid size variables
! ******************************************************************************

    pGrid%nCellsTot = pGrid%nTetsTot + pGrid%nHexsTot + pGrid%nPrisTot & 
                    + pGrid%nPyrsTot

    pGrid%nVert  = pGrid%nVertTot
    pGrid%nCells = pGrid%nCellsTot
    pGrid%nTets  = pGrid%nTetsTot
    pGrid%nHexs  = pGrid%nHexsTot
    pGrid%nPris  = pGrid%nPrisTot
    pGrid%nPyrs  = pGrid%nPyrsTot

! ==============================================================================
!   Boundary types
! ==============================================================================
  
    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
    END IF ! global%verbLevel

    READ(ifile) pGrid%nPatches  

    IF ( pGrid%nPatches > 0 ) THEN 
      ALLOCATE(gridCENTAUR%bInfo(3,pGrid%nPatches),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bInfo')
      END IF ! global%error

      READ(ifile) ((gridCENTAUR%bInfo(i,j),j=1,pGrid%nPatches),i=1,3)

      ALLOCATE(gridCENTAUR%bName(pGrid%nPatches),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bName')
      END IF ! global%error

      READ(ifile) (gridCENTAUR%bName(i),i=1,pGrid%nPatches) 
    END IF ! pGrid%nPatches

! ==============================================================================
!   Boundary face connectivity
! ==============================================================================

    IF ( pGrid%nPatches > 0 ) THEN 
      READ(iFile) gridCENTAUR%nBTris

      IF ( gridCENTAUR%nBTris > 0 ) THEN 
        ALLOCATE(gridCENTAUR%bTri2v(3,gridCENTAUR%nBTris),STAT=errorFlag)
        global%error = errorFlag 
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bTri2v')
        END IF ! global%error

        IF ( global%verbLevel > VERBOSE_NONE ) THEN    
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary triangles...'    
        END IF ! global%verbLevel

        READ(iFile) ((gridCENTAUR%bTri2v(i,j),j=1,gridCENTAUR%nBTris),i=1,3)    
      END IF ! gridCENTAUR

      READ(iFile) gridCENTAUR%nBQuads

      IF ( gridCENTAUR%nBQuads > 0 ) THEN 
        ALLOCATE(gridCENTAUR%bQuad2v(4,gridCENTAUR%nBQuads),STAT=errorFlag)
        global%error = errorFlag 
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCENTAUR%bQuad2v')
        END IF ! global%error

        IF ( global%verbLevel > VERBOSE_NONE ) THEN     
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary quadrilaterals...'    
        END IF ! global%verbLevel

        READ(iFile) ((gridCENTAUR%bQuad2v(i,j),j=1,gridCENTAUR%nBQuads),i=1,4)    
      END IF ! gridCENTAUR
    ELSE 
      gridCENTAUR%nBTris  = 0
      gridCENTAUR%nBQuads = 0  
    END IF ! pGrid%nPatches

! ******************************************************************************
!   Check validity of connectivity arrays
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN
      CALL RFLU_CheckGridCENTAUR(pRegion)
    END IF ! global%checkLevel

! ******************************************************************************
!   Print grid statistics
! ******************************************************************************
 
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      CALL RFLU_PrintGridCENTAURInfo(pRegion) 
    END IF ! global%verbLevel   

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile, IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading binary CENTAUR grid file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
   
  END SUBROUTINE RFLU_ReadGridCENTAURBinary






! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModCENTAUR

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModCENTAUR.F90,v $
! Revision 1.4  2008/12/06 08:45:03  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:13  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/25 22:04:28  haselbac
! Changes because of sype patches
!
! Revision 1.1  2005/04/15 15:09:06  haselbac
! Initial revision
!
! Revision 1.4  2005/01/20 14:54:56  haselbac
! Added setting of nBFaces and nBFacesTot
!
! Revision 1.3  2004/11/03 17:09:24  haselbac
! Removed setting of vertex and cell flags
!
! Revision 1.2  2004/10/19 19:31:02  haselbac
! Removed renumbering of bface lists, cosmetics
!
! Revision 1.1  2004/07/06 15:15:46  haselbac
! Initial revision
!
! ******************************************************************************










