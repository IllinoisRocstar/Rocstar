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
! Purpose: Collection of routines to read and convert TETMESH grids.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModTETMESH.F90,v 1.8 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModTETMESH

  USE ModParameters
  USE ModDataTypes  
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

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

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModTETMESH.F90,v $ $Revision: 1.8 $'        

  TYPE t_gridTETMESH
    INTEGER :: nBFaces,nBQuads,nBTris,nBVert
    INTEGER, DIMENSION(:), POINTER :: bf2p,bv
    INTEGER, DIMENSION(:,:), POINTER :: bf2v    
  END TYPE t_gridTETMESH

  TYPE(t_gridTETMESH) :: gridTETMESH

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_ConvROCFLU2TETMESH, & 
            RFLU_ConvTETMESH2ROCFLU, & 
            RFLU_ReadGridTETMESH, & 
            RFLU_WriteSurfGridTETMESH

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  







! ******************************************************************************
!
! Purpose: Convert grid format from ROCFLU to TETMESH.
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

  SUBROUTINE RFLU_ConvROCFLU2TETMESH(pRegion)
  
    USE RFLU_ModBoundLists

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

    INTEGER :: errorFlag,ifl,iPatch,iv
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid,pGridTETMESH
    TYPE(t_patch), POINTER :: pPatch,pPatchTETMESH
    TYPE(t_region), TARGET :: regionTETMESH
    TYPE(t_region), POINTER :: pRegionTETMESH
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvROCFLU2TETMESH', &
                          'RFLU_ModTETMESH.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from ROCFLU to TETMESH format...'
    END IF ! global%verbLevel

! ==============================================================================
!   Set grid pointer and initialize variables
! ==============================================================================

    pGrid => pRegion%grid
  
! ******************************************************************************
!   Convert patch data structure. NOTE treat combined boundary within the 
!   Rocflu data structure so that can use the boundary-vertex list routines 
!   (which require pRegion as argument and corresponding data), then copy over 
!   into gridTETMESH data structure.
! ******************************************************************************

    pRegionTETMESH => regionTETMESH
    pGridTETMESH   => pRegionTETMESH%grid

    pRegionTETMESH%global => global
    
! ==============================================================================
!   Fill relevant data of pRegionTETMESH so can use boundary-vertex list 
!   routines
! ==============================================================================  
  
! ------------------------------------------------------------------------------
!   Patch data structure
! ------------------------------------------------------------------------------  

    pGridTETMESH%nPatches = 1

    ALLOCATE(pRegionTETMESH%patches(pGridTETMESH%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegionTETMESH%patches')
    END IF ! global%error    

! ------------------------------------------------------------------------------
!   Patch dimensions
! ------------------------------------------------------------------------------  

    gridTETMESH%nBTris  = 0 
    gridTETMESH%nBQuads = 0

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      gridTETMESH%nBTris  = gridTETMESH%nBTris  + pPatch%nBTris
      gridTETMESH%nBQuads = gridTETMESH%nBQuads + pPatch%nBQuads    
    END DO ! iPatch  

    gridTETMESH%nBFaces = gridTETMESH%nBTris + gridTETMESH%nBQuads        

    pPatchTETMESH => regionTETMESH%patches(pGridTETMESH%nPatches)

    pPatchTETMESH%nBTris     = gridTETMESH%nBTris
    pPatchTETMESH%nBQuads    = gridTETMESH%nBQuads  
    pPatchTETMESH%nBFaces    = gridTETMESH%nBFaces
    pPatchTETMESH%nBVert     = 0     
    pPatchTETMESH%nBTrisTot  = pPatchTETMESH%nBTris
    pPatchTETMESH%nBQuadsTot = pPatchTETMESH%nBQuads  
    pPatchTETMESH%nBFacesTot = pPatchTETMESH%nBFaces  
    pPatchTETMESH%nBVertTot  = 0  
  
! ------------------------------------------------------------------------------
!   Allocate memory for face connectivity and fill in
! ------------------------------------------------------------------------------    

    IF ( pPatchTETMESH%nBTris > 0 ) THEN 
      ALLOCATE(pPatchTETMESH%bTri2v(3,gridTETMESH%nBTris),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatchTETMESH%bTri2v')
      END IF ! global%error    
    END IF ! pPatchTETMESH%nBTris

    IF ( pPatchTETMESH%nBQuads > 0 ) THEN 
      ALLOCATE(pPatchTETMESH%bQuad2v(4,gridTETMESH%nBQuads),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatchTETMESH%bQuad2v')
      END IF ! global%error    
    END IF ! pPatchTETMESH%nBQuads   

    ALLOCATE(pPatchTETMESH%bf2v(4,gridTETMESH%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatchTETMESH%bf2v')
    END IF ! global%error   

    ALLOCATE(gridTETMESH%bf2p(gridTETMESH%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridTETMESH%bf2p')
    END IF ! global%error    
    
! ------------------------------------------------------------------------------
!   Denumber face list
! ------------------------------------------------------------------------------    
    
    CALL RFLU_DenumberBFaceLists(pRegion)  
    
! ------------------------------------------------------------------------------
!   Fill in face connectivity
! ------------------------------------------------------------------------------    

    gridTETMESH%nBTris  = 0 
    gridTETMESH%nBQuads = 0  
    gridTETMESH%nBFaces = 0  

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      IF ( pPatch%nBTris > 0 ) THEN 
        DO ifl = 1,pPatch%nBTris 
          gridTETMESH%nBTris = gridTETMESH%nBTris + 1
          pPatchTETMESH%bTri2v(1,gridTETMESH%nBTris) = pPatch%bTri2v(1,ifl)
          pPatchTETMESH%bTri2v(2,gridTETMESH%nBTris) = pPatch%bTri2v(2,ifl)
          pPatchTETMESH%bTri2v(3,gridTETMESH%nBTris) = pPatch%bTri2v(3,ifl)                
        END DO ! ifl
      END IF ! pPatch%nBTris

      IF ( pPatch%nBQuads > 0 ) THEN 
        DO ifl = 1,pPatch%nBQuads 
          gridTETMESH%nBQuads = gridTETMESH%nBQuads + 1
          pPatchTETMESH%bQuad2v(1,gridTETMESH%nBQuads) = pPatch%bQuad2v(1,ifl)
          pPatchTETMESH%bQuad2v(2,gridTETMESH%nBQuads) = pPatch%bQuad2v(2,ifl)
          pPatchTETMESH%bQuad2v(3,gridTETMESH%nBQuads) = pPatch%bQuad2v(3,ifl)
          pPatchTETMESH%bQuad2v(4,gridTETMESH%nBQuads) = pPatch%bQuad2v(4,ifl)                        
        END DO ! ifl    
      END IF ! pPatch%nBQuads

      DO ifl = 1,pPatch%nBFaces 
        gridTETMESH%nBFaces = gridTETMESH%nBFaces + 1
        pPatchTETMESH%bf2v(1,gridTETMESH%nBFaces) = pPatch%bf2v(1,ifl)
        pPatchTETMESH%bf2v(2,gridTETMESH%nBFaces) = pPatch%bf2v(2,ifl)
        pPatchTETMESH%bf2v(3,gridTETMESH%nBFaces) = pPatch%bf2v(3,ifl)
        pPatchTETMESH%bf2v(4,gridTETMESH%nBFaces) = pPatch%bf2v(4,ifl) 

        gridTETMESH%bf2p(gridTETMESH%nBFaces) = iPatch                       
      END DO ! ifl    
    END DO ! iPatch    
        
! ==============================================================================
!   Build patch vertex lists and renumber boundary face lists 
! ==============================================================================  

    CALL RFLU_CreateBVertexLists(pRegionTETMESH)
    CALL RFLU_BuildBVertexLists(pRegionTETMESH)
    CALL RFLU_RenumberBFaceLists(pRegionTETMESH)

    gridTETMESH%nBVert = pRegionTETMESH%patches(1)%nBVert    
      
! ==============================================================================
!   Copy lists into gridTETMESH data structure
! ==============================================================================  
    
    ALLOCATE(gridTETMESH%bf2v(4,gridTETMESH%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridTETMESH%bf2v')
    END IF ! global%error

    ALLOCATE(gridTETMESH%bv(gridTETMESH%nBVert),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridTETMESH%bv')
    END IF ! global%error  

    DO ifl = 1,gridTETMESH%nBFaces
      gridTETMESH%bf2v(1,ifl) = pPatchTETMESH%bf2v(1,ifl)
      gridTETMESH%bf2v(2,ifl) = pPatchTETMESH%bf2v(2,ifl)
      gridTETMESH%bf2v(3,ifl) = pPatchTETMESH%bf2v(3,ifl)
      gridTETMESH%bf2v(4,ifl) = pPatchTETMESH%bf2v(4,ifl)            
    END DO ! ifl

    DO iv = 1,gridTETMESH%nBVert
      gridTETMESH%bv(iv) = pPatchTETMESH%bv(iv)
    END DO ! iv

    CALL RFLU_DestroyBVertexLists(pRegionTETMESH)
 
! ==============================================================================
!   Deallocate temporary memory
! ==============================================================================

    IF ( pPatchTETMESH%nBTris > 0 ) THEN 
      DEALLOCATE(pPatchTETMESH%bTri2v,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatchTETMESH%bTri2v')
      END IF ! global%error    
    END IF ! pPatchTETMESH%nBTris

    IF ( pPatchTETMESH%nBQuads > 0 ) THEN 
      DEALLOCATE(pPatchTETMESH%bQuad2v,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatchTETMESH%bQuad2v')
      END IF ! global%error    
    END IF ! pPatchTETMESH%nBQuads   

    DEALLOCATE(pPatchTETMESH%bf2v,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatchTETMESH%bf2v')
    END IF ! global%error  

    DEALLOCATE(pRegionTETMESH%patches,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegionTETMESH%patches')
    END IF ! global%error    
  
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from ROCFLU to TETMESH format done.'
    END IF ! global%verbLevel     

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ConvROCFLU2TETMESH








! ******************************************************************************
!
! Purpose: Convert grid format from TETMESH to ROCFLU.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. The .tmi file is written by rfluconv and contains a single integer 
!      which is equal to the number of boundary patches. Knowing this quantity
!      makes the conversion easier because it is not necessary to determine 
!      the number of boundary patches from the .faces file.
!
! ******************************************************************************

  SUBROUTINE RFLU_ConvTETMESH2ROCFLU(pRegion)
  
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
    INTEGER :: dummyInteger,errorFlag,i,ic,ifc,iFile,ifl,iPatch
    INTEGER, DIMENSION(3) :: dummyIntegerArray(3)
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************  

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvTETMESH2ROCFLU', &
                          'RFLU_ModTETMESH.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Converting from TETMESH to ROCFLU format...'
    END IF ! global%verbLevel

! ==============================================================================
!   Set grid pointer and initialize variables
! ==============================================================================

    pGrid => pRegion%grid

    pGrid%nEdges    = 0
    pGrid%nEdgesTot = 0

    pGrid%nFaces    = 0
    pGrid%nFacesTot = 0

! ******************************************************************************  
!   Read information file (generated by rfluconv)
! ******************************************************************************  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading TETMESH information file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.tmi',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

    READ(iFile,*) pGrid%nPatches

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Reading TETMESH information file done.'
    END IF ! global%verbLevel

! ******************************************************************************  
!   Allocate boundary patch data structure
! ******************************************************************************  

    ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%error

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch) 

      pPatch%nBFaces = 0 
      pPatch%nBTris  = 0
      pPatch%nBQuads = 0
      pPatch%nBVert  = 0

      pPatch%iPatchGlobal = iPatch
      pPatch%iBorder      = PATCH_IBORDER_DEFAULT      
      pPatch%renumFlag    = .FALSE.    
    END DO ! iPatch

! ******************************************************************************  
!   Read boundary face file (generated by rfluconv). NOTE do not read second
!   (reserved integer) on first line from file because YAMS does not write it.
! ******************************************************************************  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading .faces file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.faces',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

    READ(iFile,*) gridTETMESH%nBFaces

    ALLOCATE(gridTETMESH%bf2v(3,gridTETMESH%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridTETMESH%bf2v')
    END IF ! global%error

    ALLOCATE(gridTETMESH%bf2p(gridTETMESH%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridTETMESH%bf2v')
    END IF ! global%error

    gridTETMESH%nBTris  = gridTETMESH%nBFaces 
    gridTETMESH%nBQuads = 0

    DO ifc = 1,gridTETMESH%nBFaces    
      READ(IF_GRID,*) dummyInteger,gridTETMESH%bf2v(1:3,ifc), &
                      gridTETMESH%bf2p(ifc),dummyIntegerArray(1:3) 
    END DO ! ifc  

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading .faces file done.'
    END IF ! global%verbLevel

! ******************************************************************************  
!   Convert boundary face data structure
! ******************************************************************************  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting boundary faces...'
    END IF ! global%verbLevel

    DO ifc = 1,gridTETMESH%nBFaces
      pPatch => pRegion%patches(gridTETMESH%bf2p(ifc))  

      pPatch%nBTris  = pPatch%nBTris  + 1
      pPatch%nBFaces = pPatch%nBFaces + 1
    END DO ! ifc

! ------------------------------------------------------------------------------
!   Set total boundary patch quantities and number of boundary faces
! ------------------------------------------------------------------------------

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

! ------------------------------------------------------------------------------
!   Build face list
! ------------------------------------------------------------------------------

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
      END IF ! global%error   

      pPatch%nBTris = 0             
    END DO ! iPatch

    DO ifc = 1,gridTETMESH%nBFaces
      pPatch => pRegion%patches(gridTETMESH%bf2p(ifc))  

      pPatch%nBTris = pPatch%nBTris + 1

      pPatch%bTri2v(1,pPatch%nBTris) = gridTETMESH%bf2v(1,ifc)
      pPatch%bTri2v(2,pPatch%nBTris) = gridTETMESH%bf2v(2,ifc)
      pPatch%bTri2v(3,pPatch%nBTris) = gridTETMESH%bf2v(3,ifc)        
    END DO ! ifc

    DEALLOCATE(gridTETMESH%bf2v,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridTETMESH%bf2v')
    END IF ! global%error

    DEALLOCATE(gridTETMESH%bf2p,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridTETMESH%bf2v')
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting boundary faces done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Convert connectivity
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting connectivity...'
    END IF ! global%verbLevel

    DO ic = 1,pGrid%nTetsTot
      dummyInteger      = pGrid%tet2v(2,ic)
      pGrid%tet2v(2,ic) = pGrid%tet2v(3,ic)
      pGrid%tet2v(3,ic) = dummyInteger
    END DO ! ic

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting connectivity done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set/allocate remaining quantities
! ******************************************************************************
 
    global%nPatches = pGrid%nPatches
 
! ==============================================================================
!   Initialize other cell types
! ==============================================================================

    pGrid%nHexsTot = 0
    pGrid%nPrisTot = 0
    pGrid%nPyrsTot = 0

    pGrid%nHexs = pGrid%nHexsTot
    pGrid%nPris = pGrid%nPrisTot
    pGrid%nPyrs = pGrid%nPyrsTot

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
                               'Converting from TETMESH to ROCFLU format done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ConvTETMESH2ROCFLU

  
  
  
  
  
! ******************************************************************************
!
! Purpose: Read grid file from TETMESH in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Part of the file is read directly into the grid data structure, while
!      the remainder needs to be processed in RFLU_ConvTETMESH2ROCFLU.F90.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadGridTETMESH(pRegion)

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
  
    CHARACTER(CHRLEN) :: dummyString,iFileName  
    INTEGER :: cvmax,cvmin,dummyInteger,errorFlag,i,ic,iFile,iv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridTETMESH', &
                          'RFLU_ModTETMESH.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading TETMESH grid file...'
    END IF ! global%verbLevel

! ==============================================================================
!   Open grid file
! ==============================================================================

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.noboite',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

! ==============================================================================
!   Header
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid

    READ(iFile,*) pGrid%nTetsTot,pGrid%nVertTot

    pGrid%nTetsMax = RFLU_SetMaxDimension(global,pGrid%nTetsTot)
    pGrid%nVertMax = RFLU_SetMaxDimension(global,pGrid%nVertTot)

    pGrid%nCellsTot = pGrid%nTetsTot ! TETMESH only generates tetrahedral grids
    pGrid%nCells    = pGrid%nCellsTot  
    pGrid%nTets     = pGrid%nTetsTot
    pGrid%nVert     = pGrid%nVertTot

! ==============================================================================
!   Connectivity
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Connectivity...'
    END IF ! global%verbLevel

    ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
    END IF ! global%error  

    READ(iFile,*) ((pGrid%tet2v(i,ic),i=1,4),ic=1,pGrid%nTetsTot)

! ==============================================================================
!   Coordinates
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
    END IF ! global%verbLevel

    ALLOCATE(pGrid%xyz(3,pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%xyz')
    END IF ! global%error  

    READ(iFile,*) ((pGrid%xyz(i,iv),i=1,3),iv=1,pGrid%nVertTot)

! ==============================================================================
!   Close grid file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error 

! ==============================================================================
!   Checking - only valid for tetrahedral grids, cf. RFLU_ReadGridCENTAUR.F90
! ==============================================================================

    IF ( global%checkLevel > CHECK_NONE ) THEN  
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                 'Checking cell connectivity array entries...'
      END IF ! global%verbLevel

      cvmin = MINVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))
      cvmax = MAXVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))    

      IF ( cvmin /= 1 .OR. cvmax /= pGrid%nVertTot ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN  
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel
        CALL ErrorStop(global,ERR_VERTEX_NUMBER,__LINE__)  
      END IF ! cvmin 
      
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
              'Checking cell connectivity array entries done.'
      END IF ! global%verbLevel      
    END IF ! global%checkLevel

! ******************************************************************************
!   Print grid statistics
! ******************************************************************************
 
    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)')        SOLVER_NAME,'Grid Statistics:'
      WRITE(STDOUT,'(A,5X,A,11X,I9)') SOLVER_NAME,'Vertices:  ',pGrid%nVertTot                                     
      WRITE(STDOUT,'(A,5X,A,11X,I9)') SOLVER_NAME,'Tetrahedra:',pGrid%nTetsTot  
    END IF ! global%verbLevel 

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading TETMESH grid file done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)        

  END SUBROUTINE RFLU_ReadGridTETMESH  

  
  
  
  
  
! ******************************************************************************
!
! Purpose: Write surface grid in ASCII TETMESH/YAMS format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine has the capability to write quadrilateral faces to the
!      TETMESH file, but the converter in rfluprep cannot yet handle 
!      quadrilateral boundary faces. The problem is that TETMESH will convert
!      the quadrilateral faces to triangular faces internally, but it does not
!      appear to write out a new .faces file, so it is impossible for rfluprep
!      to know in which way the faces were split. This knowledge is needed to
!      correctly reconstruct the boundary-face list so that it is consistent
!      with the tetrahedra.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteSurfGridTETMESH(pRegion)

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
    INTEGER :: dummyInteger,errorFlag,ifc,ifcType,iFile,iv
    INTEGER, DIMENSION(4) :: dummyIntegerArray
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteSurfGridTETMESH', &
                          'RFLU_ModTETMESH.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing TETMESH surface grid...'
    END IF ! global%verbLevel

! ==============================================================================
!   Set pointers and variables
! ==============================================================================

    pGrid => pRegion%grid

    dummyInteger = 1
    dummyIntegerArray(:) = dummyInteger

! ******************************************************************************
!   Writing faces file
! ******************************************************************************

! ==============================================================================
!   Open file
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing .faces file...'
    END IF ! global%verbLevel  

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.faces',iFileName)  

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Write file
! ==============================================================================

    WRITE(IF_GRID,*) gridTETMESH%nBFaces,dummyInteger

    DO ifc = 1,gridTETMESH%nBFaces 
      IF ( gridTETMESH%bf2v(4,ifc) == VERT_NONE ) THEN 
        ifcType = 3
      ELSE
        ifcType = 4
      END IF ! gridTETMESH%bf2v

      WRITE(IF_GRID,*) ifcType,gridTETMESH%bf2v(1:ifcType,ifc), &
                       gridTETMESH%bf2p(ifc),dummyIntegerArray(1:ifcType)    
    END DO ! ifc

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing .faces file done.'
    END IF ! global%verbLevel     
     
! ******************************************************************************
!   Writing points file
! ******************************************************************************

! ==============================================================================
!   Open file
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing .points file...'
    END IF ! global%verbLevel  

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.points',iFileName)  

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Write file
! ==============================================================================

    WRITE(IF_GRID,*) gridTETMESH%nBVert

    DO iv = 1,gridTETMESH%nBVert
      WRITE(IF_GRID,*) pGrid%xyz(XCOORD:ZCOORD,iv),dummyInteger
    END DO ! iv

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error   

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing .points file done.'
    END IF ! global%verbLevel  

! ******************************************************************************
!   Writing information file (for ROCFLU use only, in rfluprep)
! ******************************************************************************

! ==============================================================================
!   Open file
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing .tmi file...'
    END IF ! global%verbLevel  

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.tmi',iFileName)  

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Write file
! ==============================================================================

    WRITE(IF_GRID,*) pGrid%nPatches
  
! ==============================================================================
!   Close file
! ==============================================================================
  
    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing .tmi file done.'
    END IF ! global%verbLevel 
  
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing TETMESH surface grid done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteSurfGridTETMESH  
  
  
  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModTETMESH

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModTETMESH.F90,v $
! Revision 1.8  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/18 21:02:37  haselbac
! Fixed checks; bug arose bcos of max-dimensioned arrays for sype
!
! Revision 1.5  2006/03/25 21:58:45  haselbac
! Changes because of sype patches
!
! Revision 1.4  2005/01/20 14:52:29  haselbac
! Added setting of nBFaces and nBFacesTot
!
! Revision 1.3  2004/11/03 17:04:32  haselbac
! Removed code related to vertex and tet flags
!
! Revision 1.2  2004/10/19 19:28:33  haselbac
! Rewrote to use new renumbering module
!
! Revision 1.1  2004/07/06 15:14:33  haselbac
! Initial revision
!
! ******************************************************************************
  
  










