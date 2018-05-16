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
! Purpose: Collection of routines to read and convert COBALT grids.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModCOBALT.F90,v 1.5 2008/12/06 08:45:03 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModCOBALT

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
    RCSIdentString = '$RCSfile: RFLU_ModCOBALT.F90,v $ $Revision: 1.5 $'        

  TYPE t_gridCOBALT
    INTEGER :: nFaces,nFacesPerCellMax,nMappings,nPatches,nQuads,nTris, &
               nVertPerFaceMax
    INTEGER, DIMENSION(:), POINTER :: nvpf
    INTEGER, DIMENSION(:,:), POINTER :: f2c,f2v,patch2bc   
  END TYPE t_gridCOBALT

  TYPE(t_gridCOBALT) :: gridCOBALT

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_ConvCOBALT2ROCFLU, & 
            RFLU_ReadGridCOBALT

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS








! ******************************************************************************
!
! Purpose: Convert grid format from COBALT to ROCFLU.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine is very long, but easily understood due to structure and 
!      extensive comments. Breaking it into pieces might be more confusing.
!
! ******************************************************************************

  SUBROUTINE RFLU_ConvCOBALT2ROCFLU(pRegion)

    USE ModSortSearch 

    USE RFLU_ModCellMapping
    USE RFLU_ModDimensions, ONLY: RFLU_SetMaxDimensions
    USE RFLU_ModGeometry, ONLY: RFLU_ComputeApproxCentroids, & 
                                RFLU_CreateApproxCentroids, &
                                RFLU_DestroyApproxCentroids
    USE RFLU_ModGrid

    USE ModInterfaces, ONLY: FaceCentroidTria, & 
                             FaceCentroidQuad, & 
                             FaceVectorTria, &
                             FaceVectorQuad

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
    INTEGER :: cntr,c1,c2,errorFlag,iBegMax,iBegMin,iBeg1,iBeg2,icg,icl, &
               iEndMax,iEndMin,iEnd1,iEnd2,ifg,iFile,ifl,iloc,iMap,iMap2, &
               iPatch,ivg,ivl,j,v1,v2,v3,v4
    INTEGER :: flag(4),key(4),vert(4),vertTemp(4)
    INTEGER, DIMENSION(:), ALLOCATABLE :: cellMapp,cellType
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: cellDegr
    REAL(RFREAL) :: cofgAppX,cofgAppY,cofgAppZ,dotProd,fCenX,fCenY,fCenZ, &
                    fVecX,fVecY,fVecZ
    REAL(RFREAL) :: xyzNodes(XCOORD:ZCOORD,4)
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvCOBALT2ROCFLU', &
                          'RFLU_ModCOBALT.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Converting from COBALT to ROCFLU format...'
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
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(cellDegr(FACE_TYPE_TRI:FACE_TYPE_QUAD,gridCOBALT%nFaces), & 
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cellDegr')
    END IF ! global%error

    DO icg = 1,pGrid%nCellsTot
      cellDegr(FACE_TYPE_TRI,icg)  = 0
      cellDegr(FACE_TYPE_QUAD,icg) = 0
    END DO ! icg

    ALLOCATE(cellType(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cellType')
    END IF ! global%error

    ALLOCATE(cellMapp(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cellMapp')
    END IF ! global%error

! ******************************************************************************
!   Determine number of cells for each cell type
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining cell types...'
    END IF ! global%verbLevel

! ==============================================================================
!   Determine number of triangular and quadrilateral faces per cell
! ==============================================================================

    DO ifg = 1,gridCOBALT%nFaces
      c1 = gridCOBALT%f2c(1,ifg)
      c2 = gridCOBALT%f2c(2,ifg)

      IF ( gridCOBALT%nvpf(ifg) == 3 ) THEN 
        IF ( c1 > 0 ) THEN 
          cellDegr(FACE_TYPE_TRI,c1) = cellDegr(FACE_TYPE_TRI,c1) + 1
        END IF ! c1

        IF ( c2 > 0 ) THEN 
          cellDegr(FACE_TYPE_TRI,c2) = cellDegr(FACE_TYPE_TRI,c2) + 1
        END IF ! c2
      ELSE IF ( gridCOBALT%nvpf(ifg) == 4 ) THEN 
        IF ( c1 > 0 ) THEN 
          cellDegr(FACE_TYPE_QUAD,c1) = cellDegr(FACE_TYPE_QUAD,c1) + 1
        END IF ! c1

        IF ( c2 > 0 ) THEN 
          cellDegr(FACE_TYPE_QUAD,c2) = cellDegr(FACE_TYPE_QUAD,c2) + 1
        END IF ! c2    
      ELSE ! Should be impossible
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! gridCOBALT%nvpf
    END DO ! ifg

! ==============================================================================
!   Determine cell types using number of triangular and quadrilateral faces
! ==============================================================================

    pGrid%nTetsTot = 0
    pGrid%nHexsTot = 0
    pGrid%nPrisTot = 0
    pGrid%nPyrsTot = 0

    DO icg = 1,pGrid%nCellsTot
      IF ( cellDegr(FACE_TYPE_TRI ,icg) == 4 .AND. & 
           cellDegr(FACE_TYPE_QUAD,icg) == 0 ) THEN ! Tetrahedron
        pGrid%nTetsTot = pGrid%nTetsTot + 1 
        cellType(icg) = CELL_TYPE_TET
        cellMapp(icg) = pGrid%nTetsTot
      ELSE IF ( cellDegr(FACE_TYPE_TRI ,icg) == 0 .AND. & 
                cellDegr(FACE_TYPE_QUAD,icg) == 6 ) THEN ! Hexahedron 
        pGrid%nHexsTot = pGrid%nHexsTot + 1 
        cellType(icg) = CELL_TYPE_HEX 
        cellMapp(icg) = pGrid%nHexsTot               
      ELSE IF ( cellDegr(FACE_TYPE_TRI ,icg) == 2 .AND. & 
                cellDegr(FACE_TYPE_QUAD,icg) == 3 ) THEN ! Prism
        pGrid%nPrisTot = pGrid%nPrisTot + 1 
        cellType(icg) = CELL_TYPE_PRI 
        cellMapp(icg) = pGrid%nPrisTot            
      ELSE IF ( cellDegr(FACE_TYPE_TRI ,icg) == 4 .AND. & 
                cellDegr(FACE_TYPE_QUAD,icg) == 1 ) THEN ! Pyramid
        pGrid%nPyrsTot = pGrid%nPyrsTot + 1 
        cellType(icg) = CELL_TYPE_PYR
        cellMapp(icg) = pGrid%nPyrsTot       
      ELSE 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! cellDegr 
    END DO ! icg

    pGrid%nTets = pGrid%nTetsTot
    pGrid%nHexs = pGrid%nHexsTot
    pGrid%nPris = pGrid%nPrisTot 
    pGrid%nPyrs = pGrid%nPyrsTot 

    pGrid%nTetsMax = RFLU_SetMaxDimension(global,pGrid%nTetsTot)
    pGrid%nHexsMax = RFLU_SetMaxDimension(global,pGrid%nHexsTot)    
    pGrid%nPrisMax = RFLU_SetMaxDimension(global,pGrid%nPrisTot)
    pGrid%nPyrsMax = RFLU_SetMaxDimension(global,pGrid%nPyrsTot)        

! ===============================================================================
!   Write statistics
! ===============================================================================

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Cell statistics:'
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Tetrahedra:',pGrid%nTetsTot
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Hexahedra: ',pGrid%nHexsTot
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Prisms:    ',pGrid%nPrisTot 
      WRITE(STDOUT,'(A,7X,A,1X,I9)') SOLVER_NAME,'Pyramids:  ',pGrid%nPyrsTot
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining cell types done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set maximum dimensions. NOTE this is necessary because dimensioning of 
!   arrays is in terms of maximum dimensions.
! ******************************************************************************
 
    CALL RFLU_SetMaxDimensions(pRegion)

! ******************************************************************************
!   Renumber cells. NOTE this is necessary because Rocflu assumes that cells are 
!   numbered in the order tetrahedra, hexahedra, prisms, and pyramids, but 
!   Cobalt violates this assumption.
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Renumbering cells...'
    END IF ! global%verbLevel

! ==============================================================================
!   Build cell mapping
! ==============================================================================

    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_BuildLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)    

! ==============================================================================
!   Renumber cells, make face list consistent
! ==============================================================================

    DO ifg = 1,gridCOBALT%nFaces
      DO j = 1,2
        c1 = gridCOBALT%f2c(j,ifg) 

        IF ( c1 > 0 ) THEN 
          icl = cellMapp(c1)

          SELECT CASE ( cellType(c1) ) 
            CASE ( CELL_TYPE_TET ) 
              icg = pGrid%tet2CellGlob(icl)
            CASE ( CELL_TYPE_HEX ) 
              icg = pGrid%hex2CellGlob(icl)
            CASE ( CELL_TYPE_PRI ) 
              icg = pGrid%pri2CellGlob(icl)
            CASE ( CELL_TYPE_PYR ) 
              icg = pGrid%pyr2CellGlob(icl)
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)        
          END SELECT ! cellType

          gridCOBALT%f2c(j,ifg) = icg
        END IF ! c1 
      END DO ! j
    END DO ! ifg

! ==============================================================================
!   Deallocate part of temporary memory
! ==============================================================================
  
    DEALLOCATE(cellMapp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cellMapp')
    END IF ! global%error  

    DEALLOCATE(cellType,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cellType')
    END IF ! global%error  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Renumbering cells done.'
    END IF ! global%verbLevel    

! ******************************************************************************
!   Building cell connectivity. NOTE in the following, the cell connectivity is 
!   built first without checking for correctness for convenience. The 
!   connectivity will be checked later in two passes.  
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building cell connectivity...'
    END IF ! global%verbLevel

    IF ( pGrid%nTetsMax > 0 ) THEN 
      ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
      END IF ! global%error

      DO icl = 1,pGrid%nTetsMax
        pGrid%tet2v(1,icl) = C2V_INIT 
        pGrid%tet2v(2,icl) = C2V_INIT 
        pGrid%tet2v(3,icl) = C2V_INIT 
        pGrid%tet2v(4,icl) = C2V_INIT                    
      END DO ! icl    
    END IF ! pGrid%nTetsMax

    IF ( pGrid%nHexsMax > 0 ) THEN 
      ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2v')
      END IF ! global%error

      DO icl = 1,pGrid%nHexsMax
        pGrid%hex2v(1,icl) = C2V_INIT 
        pGrid%hex2v(2,icl) = C2V_INIT 
        pGrid%hex2v(3,icl) = C2V_INIT 
        pGrid%hex2v(4,icl) = C2V_INIT
        pGrid%hex2v(5,icl) = C2V_INIT 
        pGrid%hex2v(6,icl) = C2V_INIT 
        pGrid%hex2v(7,icl) = C2V_INIT 
        pGrid%hex2v(8,icl) = C2V_INIT                          
      END DO ! icl     
    END IF ! pGrid%nHexsMax

    IF ( pGrid%nPrisMax > 0 ) THEN 
      ALLOCATE(pGrid%pri2v(6,pGrid%nPrisMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pri2v')
      END IF ! global%error

      DO icl = 1,pGrid%nPrisMax
        pGrid%pri2v(1,icl) = C2V_INIT 
        pGrid%pri2v(2,icl) = C2V_INIT 
        pGrid%pri2v(3,icl) = C2V_INIT 
        pGrid%pri2v(4,icl) = C2V_INIT
        pGrid%pri2v(5,icl) = C2V_INIT 
        pGrid%pri2v(6,icl) = C2V_INIT 
      END DO ! icl     
    END IF ! pGrid%nPrisMax

    IF ( pGrid%nPyrsMax > 0 ) THEN 
      ALLOCATE(pGrid%pyr2v(5,pGrid%nPyrsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pyr2v')
      END IF ! global%error

      DO icl = 1,pGrid%nPyrsMax
        pGrid%pyr2v(1,icl) = C2V_INIT 
        pGrid%pyr2v(2,icl) = C2V_INIT 
        pGrid%pyr2v(3,icl) = C2V_INIT 
        pGrid%pyr2v(4,icl) = C2V_INIT
        pGrid%pyr2v(5,icl) = C2V_INIT 
      END DO ! icl    
    END IF ! pGrid%nPyrsMax

! ==============================================================================
!   Loop over faces and visit interior cells for each face. NOTE do not use 
!   quadrilateral faces to define prism because it is easier to use the 
!   triangular faces (no overlap of vertices) and then use the quadrilateral 
!   faces to get the correct relative orientation. NOTE for pyramids, if first 
!   face to be added is a triangular face, store in positions 1 to 3 of pyr2v 
!   array despite these being the wrong locations. When quadrilateral face is 
!   added to existing triangular face, make sure that these 3 vertices are 
!   taken into account.
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot ! Reinitialize
      cellDegr(FACE_TYPE_TRI,icg)  = 0
      cellDegr(FACE_TYPE_QUAD,icg) = 0
    END DO ! icg

    DO ifg = 1,gridCOBALT%nFaces
      DO j = 1,2
        c1 = gridCOBALT%f2c(j,ifg)
      
        IF ( c1 > 0 ) THEN ! interior cell
      
! ------------------------------------------------------------------------------
!         Triangular face
! ------------------------------------------------------------------------------

          IF ( gridCOBALT%nvpf(ifg) == 3 ) THEN         
            SELECT CASE ( pGrid%cellGlob2Loc(1,c1) )
                   
! ----------- Tetrahedron ------------------------------------------------------
           
              CASE ( CELL_TYPE_TET ) 
                IF ( cellDegr(FACE_TYPE_TRI,c1) < 2 ) THEN 
                  cellDegr(FACE_TYPE_TRI,c1) = cellDegr(FACE_TYPE_TRI,c1) + 1

                  icl = pGrid%cellGlob2Loc(2,c1)

                  IF ( cellDegr(FACE_TYPE_TRI,c1) == 1 ) THEN 
                    pGrid%tet2v(1:3,icl) = gridCOBALT%f2v(1:3,ifg)
                  ELSE 
                    vert(1:3) = pGrid%tet2v(1:3,icl)
                    CALL QuickSortInteger(vert(1:3),3)

                    DO ivl = 1,3
                      CALL BinarySearchInteger(vert,3,gridCOBALT%f2v(ivl,ifg), &
                                               iloc)

                      IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                        pGrid%tet2v(4,icl) = gridCOBALT%f2v(ivl,ifg)
                        EXIT
                      END IF ! iloc                      
                    END DO ! ivl

                    IF ( pGrid%tet2v(4,icl) == C2V_INIT ) THEN
                      CALL ErrorStop(global,ERR_C2VLIST_INVALID,__LINE__)          
                    END IF ! pGrid%tet2v
                  END IF ! cellDegr
                END IF ! cellDegr 
              
! ----------- Prism ------------------------------------------------------------
                             
              CASE ( CELL_TYPE_PRI ) 
                IF ( cellDegr(FACE_TYPE_TRI,c1) < 2 ) THEN 
                  cellDegr(FACE_TYPE_TRI,c1) = cellDegr(FACE_TYPE_TRI,c1) + 1

                  icl = pGrid%cellGlob2Loc(2,c1)

                  IF ( cellDegr(FACE_TYPE_TRI,c1) == 1 ) THEN 
                    pGrid%pri2v(1:3,icl) = gridCOBALT%f2v(1:3,ifg)
                  ELSE 
                    pGrid%pri2v(4:6,icl) = gridCOBALT%f2v(1:3,ifg)             
                  END IF ! cellDegr
                END IF ! cellDegr 
              
! ----------- Pyramid ----------------------------------------------------------

              CASE ( CELL_TYPE_PYR ) 
                IF ( cellDegr(FACE_TYPE_TRI,c1) < 1 ) THEN 
                  cellDegr(FACE_TYPE_TRI,c1) = cellDegr(FACE_TYPE_TRI,c1) + 1

                  icl = pGrid%cellGlob2Loc(2,c1)

                  IF ( cellDegr(FACE_TYPE_QUAD,c1) == 0 ) THEN 
                    pGrid%pyr2v(1,icl) = gridCOBALT%f2v(1,ifg)
                    pGrid%pyr2v(2,icl) = gridCOBALT%f2v(2,ifg)
                    pGrid%pyr2v(3,icl) = gridCOBALT%f2v(3,ifg)
                  ELSE 
                    vert(1:4) = pGrid%pyr2v(1:4,icl)   
                    CALL QuickSortInteger(vert(1:4),4)

                    DO ivl = 1,3
                      CALL BinarySearchInteger(vert(1:4),4, &
                                               gridCOBALT%f2v(ivl,ifg),iloc)

                      IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                        pGrid%pyr2v(5,icl) = gridCOBALT%f2v(ivl,ifg)
                      END IF ! iloc                      
                    END DO ! ivl  

                    IF ( pGrid%pyr2v(5,icl) == C2V_INIT ) THEN
                      CALL ErrorStop(global,ERR_C2VLIST_INVALID,__LINE__)      
                    END IF ! pGrid%pyr2v 
                  END IF ! cellDegr                
                END IF ! cellDegr                                            
            END SELECT ! pGrid%cellGlob2Loc
          
! ------------------------------------------------------------------------------
!         Quadrilateral face
! ------------------------------------------------------------------------------
           
          ELSE IF ( gridCOBALT%nvpf(ifg) == 4 ) THEN 
            SELECT CASE ( pGrid%cellGlob2Loc(1,c1) )

! ----------- Hexahedron -------------------------------------------------------

              CASE ( CELL_TYPE_HEX ) 
                IF ( cellDegr(FACE_TYPE_QUAD,c1) < 2 ) THEN 
                  icl = pGrid%cellGlob2Loc(2,c1)

                  IF ( cellDegr(FACE_TYPE_QUAD,c1) == 0 ) THEN ! Face 1
                    cellDegr(FACE_TYPE_QUAD,c1) = cellDegr(FACE_TYPE_QUAD,c1) & 
                                                + 1
                    pGrid%hex2v(1:4,icl) = gridCOBALT%f2v(1:4,ifg)
                  ELSE 
                    vert(1:4) = pGrid%hex2v(1:4,icl)   
                    CALL QuickSortInteger(vert(1:4),4) 

                    cntr = 0

                    DO ivl = 1,4
                      ivg = gridCOBALT%f2v(ivl,ifg) 

                      CALL BinarySearchInteger(vert(1:4),4,ivg,iloc)

                      IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                        cntr = cntr + 1
                      END IF ! iloc
                    END DO ! ivl

                    IF ( cntr == 4 ) THEN ! Face 6
                      cellDegr(FACE_TYPE_QUAD,c1) = & 
                        cellDegr(FACE_TYPE_QUAD,c1) + 1                  
                      pGrid%hex2v(5:8,icl) = gridCOBALT%f2v(1:4,ifg)             
                    END IF ! cntr
                  END IF ! cellDegr
                END IF ! cellDegr              

! ----------- Pyramid ----------------------------------------------------------

              CASE ( CELL_TYPE_PYR ) 
                IF ( cellDegr(FACE_TYPE_QUAD,c1) < 1 ) THEN 
                  cellDegr(FACE_TYPE_QUAD,c1) = cellDegr(FACE_TYPE_QUAD,c1) + 1

                  icl = pGrid%cellGlob2Loc(2,c1)               

                  IF ( cellDegr(FACE_TYPE_TRI,c1) == 0 ) THEN 
                    pGrid%pyr2v(1:4,icl) = gridCOBALT%f2v(1:4,ifg)
                  ELSE 
                    vertTemp(1:3) = pGrid%pyr2v(1:3,icl)

                    pGrid%pyr2v(1:4,icl) = gridCOBALT%f2v(1:4,ifg)      
                    vert(1:4) = pGrid%pyr2v(1:4,icl)   
                    CALL QuickSortInteger(vert(1:4),4)                            

                    DO ivl = 1,3
                      CALL BinarySearchInteger(vert(1:4),4, &
                                               vertTemp(ivl),iloc)

                      IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                        pGrid%pyr2v(5,icl) = vertTemp(ivl)
                      END IF ! iloc                      
                    END DO ! ivl  

                    IF ( pGrid%pyr2v(5,icl) == C2V_INIT ) THEN
                      CALL ErrorStop(global,ERR_C2VLIST_INVALID,__LINE__)
                    END IF ! pGrid%pyr2v
                  END IF ! cellDegr                
                END IF ! cellDegr 
            END SELECT ! pGrid%cellGlob2Loc

          ELSE ! Should be impossible by now, but leave for defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
          END IF ! gridCOBALT%nvpf
        END IF ! c1
      END DO ! j
    END DO ! ifg

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building cell connectivity done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Checking connectivity, pass 1 (checking orientation). For cells other than 
!   tetrahedra, this is only done for those faces which were added above, and 
!   for pyramids, it is only done for the quadrilateral face. It is not 
!   necessary for the triangular faces of the pyramids to be checked because 
!   reversing the order would mean that the order of the quadrilateral face is 
!   wrong.
! ******************************************************************************

! ==============================================================================
!   Compute approximate centroids. NOTE cell mapping (built above) required for
!   computation of approximate centroids.
! ==============================================================================
   
    CALL RFLU_CreateApproxCentroids(pRegion)
    CALL RFLU_ComputeApproxCentroids(pRegion)
  
! ==============================================================================
!   Check connectivity by checking orientation of faces. This is done by 
!   computing the dot product between the face normal and the relative position 
!   vector between the cell centroid and the face centroid. This must be 
!   positive if the face normal is outward pointing, which means that the face 
!   vertices are oriented properly.
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Checking cell connectivity...'

      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Pass 1...'
      END IF ! global%verbLevel
    END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!   Tetrahedra
! ------------------------------------------------------------------------------
    
    IF ( pGrid%nTetsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel  

      DO icl = 1,pGrid%nTetsTot
        icg = pGrid%tet2CellGlob(icl)

        cofgAppX = pGrid%cofgApp(XCOORD,icg)
        cofgAppY = pGrid%cofgApp(YCOORD,icg)
        cofgAppZ = pGrid%cofgApp(ZCOORD,icg)        

        cntr = 0

        DO ifl = 1,4
          v1 = pGrid%tet2v(f2vTet(1,ifl),icl)
          v2 = pGrid%tet2v(f2vTet(2,ifl),icl) 
          v3 = pGrid%tet2v(f2vTet(3,ifl),icl)        

          xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
          xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
          xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)

          CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fVecX,fVecY,fVecZ)
          CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fCenX,fCenY,fCenZ)

          dotProd = fVecX*(fCenX - cofgAppX) & 
                  + fVecY*(fCenY - cofgAppY) &
                  + fVecZ*(fCenZ - cofgAppZ)

          IF ( dotProd < 0.0_RFREAL ) THEN ! Only one face can be wrong
            cntr = cntr + 1

            IF ( cntr > 1 ) THEN 
              CALL ErrorStop(global,ERR_FACE_ORIENT,__LINE__)
            END IF ! cntr

            pGrid%tet2v(f2vTet(1,ifl),icl) = v1
            pGrid%tet2v(f2vTet(2,ifl),icl) = v3
            pGrid%tet2v(f2vTet(3,ifl),icl) = v2               
          END IF ! dotProd
        END DO ! ifl
      END DO ! icl
    END IF ! pGrid%nTetsTot
  
! ------------------------------------------------------------------------------
!   Hexahedra   
! ------------------------------------------------------------------------------
  
    IF ( pGrid%nHexsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Hexahedra...'
      END IF ! global%verbLevel  

      DO icl = 1,pGrid%nHexsTot
        icg = pGrid%hex2CellGlob(icl)

        cofgAppX = pGrid%cofgApp(XCOORD,icg)
        cofgAppY = pGrid%cofgApp(YCOORD,icg)
        cofgAppZ = pGrid%cofgApp(ZCOORD,icg)       

        DO ifl = 1,6,5 ! Note loop over opposing faces
          v1 = pGrid%hex2v(f2vHex(1,ifl),icl)
          v2 = pGrid%hex2v(f2vHex(2,ifl),icl) 
          v3 = pGrid%hex2v(f2vHex(3,ifl),icl)        
          v4 = pGrid%hex2v(f2vHex(4,ifl),icl)    

          xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
          xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
          xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
          xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

          CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fVecX,fVecY,fVecZ)
          CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fCenX,fCenY,fCenZ)

          dotProd = fVecX*(fCenX - cofgAppX) & 
                  + fVecY*(fCenY - cofgAppY) &
                  + fVecZ*(fCenZ - cofgAppZ)

          IF ( dotProd < 0.0_RFREAL ) THEN
            pGrid%hex2v(f2vHex(1,ifl),icl) = v1
            pGrid%hex2v(f2vHex(2,ifl),icl) = v4
            pGrid%hex2v(f2vHex(3,ifl),icl) = v3
            pGrid%hex2v(f2vHex(4,ifl),icl) = v2                         
          END IF ! dotProd        
        END DO ! ifl
      END DO ! icl
    END IF ! pGrid%nHexsTot 
  
! ------------------------------------------------------------------------------
!   Prisms  
! ------------------------------------------------------------------------------
  
    IF ( pGrid%nPrisTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Prisms...'
      END IF ! global%verbLevel  

      DO icl = 1,pGrid%nPrisTot
        icg = pGrid%pri2CellGlob(icl)

        cofgAppX = pGrid%cofgApp(XCOORD,icg)
        cofgAppY = pGrid%cofgApp(YCOORD,icg)
        cofgAppZ = pGrid%cofgApp(ZCOORD,icg)        

        DO ifl = 1,5,4 ! NOTE loop only over triangular faces
          v1 = pGrid%pri2v(f2vPri(1,ifl),icl)
          v2 = pGrid%pri2v(f2vPri(2,ifl),icl) 
          v3 = pGrid%pri2v(f2vPri(3,ifl),icl)        

          xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
          xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
          xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)

          CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fVecX,fVecY,fVecZ)
          CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fCenX,fCenY,fCenZ)

          dotProd = fVecX*(fCenX - cofgAppX) & 
                  + fVecY*(fCenY - cofgAppY) &
                  + fVecZ*(fCenZ - cofgAppZ)

          IF ( dotProd < 0.0_RFREAL ) THEN
            pGrid%pri2v(f2vPri(1,ifl),icl) = v1
            pGrid%pri2v(f2vPri(2,ifl),icl) = v3
            pGrid%pri2v(f2vPri(3,ifl),icl) = v2               
          END IF ! dotProd
        END DO ! ifl
      END DO ! icl
    END IF ! pGrid%nPrisTot   
   
! ------------------------------------------------------------------------------
!   Pyramids   
! ------------------------------------------------------------------------------
  
    IF ( pGrid%nPyrsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Pyramids...'
      END IF ! global%verbLevel  

      DO icl = 1,pGrid%nPyrsTot
        icg = pGrid%pyr2CellGlob(icl)

        cofgAppX = pGrid%cofgApp(XCOORD,icg)
        cofgAppY = pGrid%cofgApp(YCOORD,icg)
        cofgAppZ = pGrid%cofgApp(ZCOORD,icg)       

        DO ifl = 1,1 ! Note loop only over first face
          v1 = pGrid%pyr2v(f2vPyr(1,ifl),icl)
          v2 = pGrid%pyr2v(f2vPyr(2,ifl),icl) 
          v3 = pGrid%pyr2v(f2vPyr(3,ifl),icl)        
          v4 = pGrid%pyr2v(f2vPyr(4,ifl),icl)    

          xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
          xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
          xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
          xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

          CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fVecX,fVecY,fVecZ)
          CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fCenX,fCenY,fCenZ)

          dotProd = fVecX*(fCenX - cofgAppX) & 
                  + fVecY*(fCenY - cofgAppY) &
                  + fVecZ*(fCenZ - cofgAppZ)

          IF ( dotProd < 0.0_RFREAL ) THEN
            pGrid%pyr2v(f2vPyr(1,ifl),icl) = v1
            pGrid%pyr2v(f2vPyr(2,ifl),icl) = v4
            pGrid%pyr2v(f2vPyr(3,ifl),icl) = v3       
            pGrid%pyr2v(f2vPyr(4,ifl),icl) = v2   
          END IF ! dotProd
        END DO ! ifl
      END DO ! icl
    END IF ! pGrid%nPyrsTot    
  
! ******************************************************************************
!   Checking orientation, pass 2. At this point, tetrahedra and pyramids are 
!   oriented properly. However, hexahedra and prisms are not necessarily 
!   oriented properly, because although the two individual faces are oriented 
!   properly, they may not be aligned correctly relative to each other. To find  
!   correct alignment, use third face.
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Pass 2...'
    END IF ! global%verbLevel  
  
! ==============================================================================
!   Loop over faces
! ==============================================================================
  
    DO ifg = 1,gridCOBALT%nFaces
      DO j = 1,2
        c1 = gridCOBALT%f2c(j,ifg)

        IF ( c1 > 0 ) THEN 
          IF ( gridCOBALT%nvpf(ifg) == 4 ) THEN         
            SELECT CASE ( pGrid%cellGlob2Loc(1,c1) )  
          
! ------------------------------------------------------------------------------
!             Hexahedra
! ------------------------------------------------------------------------------

              CASE ( CELL_TYPE_HEX )               
                IF ( cellDegr(FACE_TYPE_QUAD,c1) < 3 ) THEN             
                  icl = pGrid%cellGlob2Loc(2,c1)
                      
! --------------- Find faces which can be used to determine orientation --------

                  vertTemp(1:4) = pGrid%hex2v(1:4,icl)
                  CALL QuickSortInteger(vertTemp(1:4),4)

                  cntr = 0

                  DO ivl = 1,4
                    ivg = gridCOBALT%f2v(ivl,ifg) 

                    CALL BinarySearchInteger(vertTemp(1:4),4,ivg,iloc)

                    IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                      cntr = cntr + 1
                    END IF ! iloc
                  END DO ! ivl              
                                  
! --------------- Need faces which share 2 vertices with face 1 ----------------
            
                  IF ( cntr == 2 ) THEN ! Have proper candidate face
                    cellDegr(FACE_TYPE_QUAD,c1) = cellDegr(FACE_TYPE_QUAD,c1) & 
                                                + 1

                    v1 = gridCOBALT%f2v(1,ifg)
                    v2 = gridCOBALT%f2v(2,ifg)              
                    v3 = gridCOBALT%f2v(3,ifg)              
                    v4 = gridCOBALT%f2v(4,ifg)                          

                    xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
                    xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
                    xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
                    xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

                    CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                                        fVecX,fVecY,fVecZ)
                    CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                                          fCenX,fCenY,fCenZ)

! ----------------- Check orientation, reverse if necessary --------------------

                    dotProd = fVecX*(fCenX - pGrid%cofgApp(XCOORD,c1)) & 
                            + fVecY*(fCenY - pGrid%cofgApp(YCOORD,c1)) &
                            + fVecZ*(fCenZ - pGrid%cofgApp(ZCOORD,c1))

                    IF ( dotProd < 0.0_RFREAL ) THEN
                      vert(1) = v1
                      vert(2) = v4
                      vert(3) = v3
                      vert(4) = v2 
                    ELSE 
                      vert(1) = v1
                      vert(2) = v2
                      vert(3) = v3
                      vert(4) = v4  
                    END IF ! dotProd

! ----------------- Determine local face number --------------------------------

                    vertTemp(1:4) = pGrid%hex2v(1:4,icl)
                    CALL QuickSortInteger(vertTemp(1:4),4)

                    DO ivl = 1,4
                      ivg = pGrid%hex2v(ivl,icl)

                      CALL BinarySearchInteger(vertTemp(1:4),4,ivg,iloc)

                      IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                        key(iloc) = ivl
                      ELSE
                        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                      END IF ! iloc
                    END DO ! ivl

                    cntr = 0
                    flag(1:4) = 0

                    DO ivl = 1,4
                      ivg = vert(ivl)

                      CALL BinarySearchInteger(vertTemp(1:4),4,ivg,iloc)

                      IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                        flag(key(iloc)) = 1
                      ELSE 
                        cntr = cntr + 1
                      END IF ! iloc
                    END DO ! ivl

! ----------------- Cycle vertices of face 6 to match ordering on face 1 -------

                    IF ( flag(1) == 1 .AND. flag(2) == 1 ) THEN ! Face 2 
                      CALL CycleList(vert(1:4),4,1,pGrid%hex2v(1,icl))
                      CALL CycleList(pGrid%hex2v(5:8,icl),4,1,vert(4))
                    ELSE IF ( flag(2) == 1 .AND. flag(3) == 1 ) THEN ! Face 3
                      CALL CycleList(vert(1:4),4,1,pGrid%hex2v(2,icl)) 
                      CALL CycleList(pGrid%hex2v(5:8,icl),4,2,vert(4))
                    ELSE IF ( flag(3) == 1 .AND. flag(4) == 1 ) THEN ! Face 4
                      CALL CycleList(vert(1:4),4,1,pGrid%hex2v(3,icl)) 
                      CALL CycleList(pGrid%hex2v(5:8,icl),4,3,vert(4)) 
                    ELSE IF ( flag(4) == 1 .AND. flag(1) == 1 ) THEN ! Face 5
                      CALL CycleList(vert(1:4),4,1,pGrid%hex2v(4,icl)) 
                      CALL CycleList(pGrid%hex2v(5:8,icl),4,4,vert(4))
                    ELSE 
                      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
                    END IF ! flag
                  ELSE IF ( cntr /= 4 .AND. cntr /= 0 ) THEN 
                    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                  END IF ! cntr
                END IF ! cellDegr
              
! ------------------------------------------------------------------------------
!             Prism
! ------------------------------------------------------------------------------

              CASE ( CELL_TYPE_PRI ) 
                IF ( cellDegr(FACE_TYPE_QUAD,c1) < 1 ) THEN             
                  icl = pGrid%cellGlob2Loc(2,c1)
                      
! --------------- Find faces which can be used to determine orientation --------
                          
                  vertTemp(1:4) = pGrid%pri2v(1:4,icl)
                  CALL QuickSortInteger(vertTemp(1:4),4)

                  cntr = 0

                  DO ivl = 1,4
                    ivg = gridCOBALT%f2v(ivl,ifg) 

                    CALL BinarySearchInteger(vertTemp(1:4),4,ivg,iloc)

                    IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                      cntr = cntr + 1
                    END IF ! iloc
                  END DO ! ivl              
                                  
! --------------- Need faces which share 2 vertices with face 1 ----------------

                  IF ( cntr == 2 ) THEN ! Have proper candidate face
                    cellDegr(FACE_TYPE_QUAD,c1) = cellDegr(FACE_TYPE_QUAD,c1) & 
                                                + 1

                    v1 = gridCOBALT%f2v(1,ifg)
                    v2 = gridCOBALT%f2v(2,ifg)              
                    v3 = gridCOBALT%f2v(3,ifg)              
                    v4 = gridCOBALT%f2v(4,ifg)                          

                    xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
                    xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
                    xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
                    xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

                    CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                                        fVecX,fVecY,fVecZ)
                    CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                                          fCenX,fCenY,fCenZ)

! ----------------- Check orientation, reverse if necessary --------------------

                    dotProd = fVecX*(fCenX - pGrid%cofgApp(XCOORD,c1)) & 
                            + fVecY*(fCenY - pGrid%cofgApp(YCOORD,c1)) &
                            + fVecZ*(fCenZ - pGrid%cofgApp(ZCOORD,c1))

                    IF ( dotProd < 0.0_RFREAL ) THEN
                      vert(1) = v1
                      vert(2) = v4
                      vert(3) = v3
                      vert(4) = v2 
                    ELSE 
                      vert(1) = v1
                      vert(2) = v2
                      vert(3) = v3
                      vert(4) = v4  
                    END IF ! dotProd

! ----------------- Determine local face number --------------------------------

                    vertTemp(1:3) = pGrid%pri2v(1:3,icl)
                    CALL QuickSortInteger(vertTemp(1:3),3)

                    DO ivl = 1,3
                      ivg = pGrid%pri2v(ivl,icl)

                      CALL BinarySearchInteger(vertTemp(1:3),3,ivg,iloc)

                      IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                        key(iloc) = ivl
                      ELSE
                        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                      END IF ! iloc
                    END DO ! ivl

                    cntr = 0
                    flag(1:3) = 0

                    DO ivl = 1,4
                      ivg = vert(ivl)

                      CALL BinarySearchInteger(vertTemp(1:3),3,ivg,iloc)

                      IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                        flag(key(iloc)) = 1
                      ELSE 
                        cntr = cntr + 1
                      END IF ! iloc
                    END DO ! ivl

! ----------------- Cycle vertices of face 5 to match ordering on face 1 -------

                    IF ( flag(1) == 1 .AND. flag(2) == 1 ) THEN ! Face 2 
                      CALL CycleList(vert(1:4),4,1,pGrid%pri2v(1,icl))
                      CALL CycleList(pGrid%pri2v(4:6,icl),3,1,vert(4))
                    ELSE IF ( flag(2) == 1 .AND. flag(3) == 1 ) THEN ! Face 3
                      CALL CycleList(vert(1:4),4,1,pGrid%pri2v(2,icl)) 
                      CALL CycleList(pGrid%pri2v(4:6,icl),3,2,vert(4)) 
                    ELSE IF ( flag(3) == 1 .AND. flag(1) == 1 ) THEN ! Face 4
                      CALL CycleList(vert(1:4),4,1,pGrid%pri2v(3,icl)) 
                      CALL CycleList(pGrid%pri2v(4:6,icl),3,3,vert(4))
                    ELSE 
                      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
                    END IF ! flag
                  ELSE IF ( cntr /= 3 .AND. cntr /= 0 ) THEN 
                    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                  END IF ! cntr
                END IF ! cellDegr

            END SELECT ! pGrid%cellGlob2Loc                        
          END IF ! gridCOBALT%nvpf 
        END IF ! c1

      END DO ! j     
    END DO ! ifg

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Checking cell connectivity done.'
    END IF ! global%verbLevel  
  
! ==============================================================================
!   Deallocate remaining temporary memory
! ==============================================================================
  
    DEALLOCATE(cellDegr,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cellDegr')
    END IF ! global%error  
  
! ******************************************************************************
!   Building boundary data structure. Cobalt has boundary flags in the face 
!   lists as negative cell indices. The boundary flags do not have to be 
!   contiguous. Use mapping file to map the boundary flags to Rocflu patches. 
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building boundary face lists...'  
    END IF ! global%verbLevel   
    
! ==============================================================================
!   Read additional info on grid, required for mapping of patches
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Reading information file...'  
    END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!   Open file
! ------------------------------------------------------------------------------

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.cgi',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

! ------------------------------------------------------------------------------
!   Read file
! ------------------------------------------------------------------------------
  
    READ(iFile,*) pGrid%nPatches
    READ(iFile,*) gridCOBALT%nMappings

    ALLOCATE(gridCOBALT%patch2bc(3,gridCOBALT%nMappings),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCOBALT%patch2bc')
    END IF ! global%error   

    DO iMap = 1,gridCOBALT%nMappings
      READ(iFile,*) (gridCOBALT%patch2bc(j,iMap),j=1,3)
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
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Reading information file done.'  
    END IF ! global%verbLevel

! ==============================================================================
!   Check for consistent input - somewhat complicated...
! ==============================================================================

    IF ( global%checkLevel > CHECK_NONE ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Checking patch mapping entries...'
      END IF ! global%verbLevel

      DO iMap = 1,gridCOBALT%nMappings
        IF ( gridCOBALT%patch2bc(2,iMap) < gridCOBALT%patch2bc(1,iMap) ) THEN
          IF ( global%verbLevel > VERBOSE_NONE ) THEN  
            WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Check failed.' 
          END IF ! global%verbLevel   
          CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
        END IF ! gridCOBALT
      END DO ! iMap   

      IF ( MINVAL(gridCOBALT%patch2bc(3,:)) /= 1 .OR. & 
           MAXVAL(gridCOBALT%patch2bc(3,:)) /= pGrid%nPatches ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN              
          WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel               
        CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
      END IF ! gridCOBALT

      DO iMap = 1,gridCOBALT%nMappings
        DO iMap2 = 1,gridCOBALT%nMappings

          IF ( iMap /= iMap2 ) THEN 
            iBeg1 = gridCOBALT%patch2bc(1,iMap)
            iEnd1 = gridCOBALT%patch2bc(2,iMap)

            iBeg2 = gridCOBALT%patch2bc(1,iMap2)
            iEnd2 = gridCOBALT%patch2bc(2,iMap2)        

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
                WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel    
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iBeg1

            IF ( iEndMin >= iBegMax ) THEN
              IF ( global%verbLevel > VERBOSE_NONE ) THEN          
                WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel    
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iebmax
          END IF ! iMap

        END DO ! iMap2
      END DO ! iMap

      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                 'Checking patch mapping entries done.'
      END IF ! global%verbLevel  
    END IF ! global%checkLevel

! ******************************************************************************
!   Convert boundary information to ROCFLU format
! ******************************************************************************

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
!   Determine number of triangular and quadrilateral faces on each patch
! ==============================================================================
  
    DO ifg = 1,gridCOBALT%nFaces
      DO j = 1,2
        c1 = gridCOBALT%f2c(j,ifg)

        IF ( c1 < 0 ) THEN ! boundary cell
          NULLIFY(pPatch)

          DO iMap = 1,gridCOBALT%nMappings
            IF ( ABS(c1) >= gridCOBALT%patch2bc(1,iMap) .AND. & 
                 ABS(c1) <= gridCOBALT%patch2bc(2,iMap) ) THEN 
              iPatch = gridCOBALT%patch2bc(3,iMap)
              pPatch => pRegion%patches(iPatch)
              EXIT
            END IF ! ABS
          END DO ! iMap

          IF ( ASSOCIATED(pPatch) .EQV. .TRUE. ) THEN 
            IF ( gridCOBALT%nvpf(ifg) == 3 ) THEN 
              pPatch%nBTris = pPatch%nBTris + 1
            ELSE IF ( gridCOBALT%nvpf(ifg) == 4 ) THEN 
              pPatch%nBQuads = pPatch%nBQuads + 1
            ELSE ! Should be impossible, keep for defensive programming
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! gridCOBALT%nvpf
          ELSE 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! ASSOCIATED      
        END IF ! c1            
      END DO ! j    
    END DO ! ifg
    
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
!   Allocate memory for boundary face lists
! ==============================================================================
  
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%nBTrisMax > 0 ) THEN 
        ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
        global%error = errorFlag 
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
        END IF ! global%error 

        pPatch%nBTris = 0 ! Reused as counter below
      END IF ! pPatch%nBTrisMax

      IF ( pPatch%nBQuadsMax > 0 ) THEN 
        ALLOCATE(pPatch%bQuad2v(4,pPatch%nBQuadsMax),STAT=errorFlag)
        global%error = errorFlag 
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bQuad2v')
        END IF ! global%error 

        pPatch%nBQuads = 0 ! Reused as counter below      
      END IF ! pPatch%nBQuadsMax    
    END DO ! iPatch    

! ==============================================================================
!   Extract boundary faces and check orientation (best done here because have
!   access to other cell index)
! ==============================================================================
  
    DO ifg = 1,gridCOBALT%nFaces
      DO j = 1,2
        c1 = gridCOBALT%f2c(j,ifg)

        IF ( c1 < 0 ) THEN ! boundary cell      
          DO iMap = 1,gridCOBALT%nMappings
            IF ( ABS(c1) >= gridCOBALT%patch2bc(1,iMap) .AND. & 
                 ABS(c1) <= gridCOBALT%patch2bc(2,iMap) ) THEN 
              iPatch = gridCOBALT%patch2bc(3,iMap)
              pPatch => pRegion%patches(iPatch)
              EXIT
            END IF ! ABS
          END DO ! iMap

          IF ( ASSOCIATED(pPatch) .EQV. .TRUE. ) THEN
            c2 = gridCOBALT%f2c(3-j,ifg)

            IF ( gridCOBALT%nvpf(ifg) == 3 ) THEN 
              pPatch%nBTris = pPatch%nBTris + 1                                   

              v1 = gridCOBALT%f2v(1,ifg)
              v2 = gridCOBALT%f2v(2,ifg)
              v3 = gridCOBALT%f2v(3,ifg)                       

              xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
              xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
              xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)

              CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3), &
                                  fVecX,fVecY,fVecZ)
              CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3), &
                                    fCenX,fCenY,fCenZ)

              dotProd = fVecX*(fCenX - pGrid%cofgApp(XCOORD,c2)) & 
                      + fVecY*(fCenY - pGrid%cofgApp(YCOORD,c2)) &
                      + fVecZ*(fCenZ - pGrid%cofgApp(ZCOORD,c2))

              IF ( dotProd < 0.0_RFREAL ) THEN
                pPatch%bTri2v(1,pPatch%nBTris) = v1
                pPatch%bTri2v(2,pPatch%nBTris) = v3
                pPatch%bTri2v(3,pPatch%nBTris) = v2 
              ELSE              
                pPatch%bTri2v(1,pPatch%nBTris) = v1
                pPatch%bTri2v(2,pPatch%nBTris) = v2
                pPatch%bTri2v(3,pPatch%nBTris) = v3
              END IF ! dotProd                                                
            ELSE IF ( gridCOBALT%nvpf(ifg) == 4 ) THEN 
              pPatch%nBQuads = pPatch%nBQuads + 1

              v1 = gridCOBALT%f2v(1,ifg)
              v2 = gridCOBALT%f2v(2,ifg)
              v3 = gridCOBALT%f2v(3,ifg)                        
              v4 = gridCOBALT%f2v(4,ifg)

              xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
              xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
              xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
              xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

              CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                                  fVecX,fVecY,fVecZ)
              CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                                    fCenX,fCenY,fCenZ)

              dotProd = fVecX*(fCenX - pGrid%cofgApp(XCOORD,c2)) & 
                      + fVecY*(fCenY - pGrid%cofgApp(YCOORD,c2)) &
                      + fVecZ*(fCenZ - pGrid%cofgApp(ZCOORD,c2))

              IF ( dotProd < 0.0_RFREAL ) THEN
                pPatch%bQuad2v(1,pPatch%nBQuads) = v1
                pPatch%bQuad2v(2,pPatch%nBQuads) = v4
                pPatch%bQuad2v(3,pPatch%nBQuads) = v3                       
                pPatch%bQuad2v(4,pPatch%nBQuads) = v2
              ELSE 
                pPatch%bQuad2v(1,pPatch%nBQuads) = v1
                pPatch%bQuad2v(2,pPatch%nBQuads) = v2
                pPatch%bQuad2v(3,pPatch%nBQuads) = v3                       
                pPatch%bQuad2v(4,pPatch%nBQuads) = v4            
              END IF ! dotProd                      
            ELSE ! Should be impossible, keep for defensive programming
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! gridCOBALT%nvpf
          ELSE 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! ASSOCIATED      
        END IF ! c1

      END DO ! j    
    END DO ! ifg  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building boundary face lists done.'  
    END IF ! global%verbLevel 

! ==============================================================================
!   Deallocate temporary memory
! ==============================================================================

    DEALLOCATE(gridCOBALT%patch2bc,STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridCOBALT%patch2bc')
    END IF ! global%error       
    
! ==============================================================================
!   Destroy approximate centroids and cell mapping.
! ==============================================================================

    CALL RFLU_DestroyApproxCentroids(pRegion)
    CALL RFLU_DestroyCellMapping(pRegion)      
  
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
                               'Converting from COBALT to ROCFLU format done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ConvCOBALT2ROCFLU










! ******************************************************************************
!
! Purpose: Read grid file from COBALT in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. COBALT grid file is usually called <something>.inp, which conflicts 
!      with the Rocflu input file, so the COBALT grid file needs to be renamed
!      before being run through rfluprep. The same applies to the COBALT 
!      boundary condition file, which is called <something>.bc. 
!   2. Part of the file is read directly into the grid data structure, while
!      the remainder needs to be processed in RFLU_ConvCOBALT2ROCFLU.F90.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadGridCOBALT(pRegion)
 
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
    INTEGER :: errorFlag,ic,ifc,iFile,iv,j,nDim,nVertPerFace,nZones
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
     
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridCOBALT', &
                          'RFLU_ModCOBALT.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading COBALT grid file...'
    END IF ! global%verbLevel

! ==============================================================================
! Open grid file
! ==============================================================================

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.cgr',iFileName)

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

    READ(iFile,*) nDim,nZones,gridCOBALT%nPatches

    IF ( nDim /= 3 ) THEN 
      CALL ErrorStop(global,ERR_NDIMENS_INVALID,__LINE__)
    END IF ! nDim

    IF ( nZones /= 1 ) THEN 
      CALL ErrorStop(global,ERR_NZONES_INVALID,__LINE__)
    END IF ! nZones

    READ(iFile,*) pGrid%nVertTot,gridCOBALT%nFaces,pGrid%nCellsTot, & 
                  gridCOBALT%nVertPerFaceMax,gridCOBALT%nFacesPerCellMax

    pGrid%nVert  = pGrid%nVertTot
    pGrid%nCells = pGrid%nCellsTot  
    
    pGrid%nVertMax  = RFLU_SetMaxDimension(global,pGrid%nVertTot)
    pGrid%nCellsMax = RFLU_SetMaxDimension(global,pGrid%nCellsTot)

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

    DO iv = 1,pGrid%nVertTot
      READ(iFile,*) (pGrid%xyz(j,iv),j=1,3) 
    END DO ! iv

! ==============================================================================
!   Connectivity
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Connectivity...'
    END IF ! global%verbLevel

    ALLOCATE(gridCOBALT%f2v(gridCOBALT%nVertPerFaceMax,gridCOBALT%nFaces), & 
             STAT=errorFlag)
    global%error = errorFlag           
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCOBALT%f2v')
    END IF ! global%error

    DO ifc = 1,gridCOBALT%nFaces
      DO j = 1,gridCOBALT%nVertPerFaceMax
        gridCOBALT%f2v(j,ifc) = 0
      END DO ! j
    END DO ! ifc

    ALLOCATE(gridCOBALT%f2c(2,gridCOBALT%nFaces),STAT=errorFlag)
    global%error = errorFlag           
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCOBALT%f2c')
    END IF ! errorFlag

    ALLOCATE(gridCOBALT%nvpf(gridCOBALT%nFaces),STAT=errorFlag)
    global%error = errorFlag           
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridCOBALT%nvpf')
    END IF ! errorFlag

    gridCOBALT%nTris  = 0
    gridCOBALT%nQuads = 0

    DO ifc = 1,gridCOBALT%nFaces
      READ(iFile,*) nVertPerFace
      BACKSPACE(iFile)

      IF ( nVertPerFace == 3 ) THEN
        gridCOBALT%nTris = gridCOBALT%nTris + 1
      ELSE IF ( nVertPerFace == 4 ) THEN
        gridCOBALT%nQuads = gridCOBALT%nQuads + 1
      ELSE  
        CALL ErrorStop(global,ERR_FACETYPE_INVALID,__LINE__)
      END IF ! nVertPerFace

      READ(iFile,*) gridCOBALT%nvpf(ifc), &
                   (gridCOBALT%f2v(j,ifc),j=1,nVertPerFace), & 
                   (gridCOBALT%f2c(j,ifc),j=1,2)  
    END DO ! ifc

! ==============================================================================
!   Close grid file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error 

! ******************************************************************************
!   Print grid statistics
! ******************************************************************************
 
    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)')        SOLVER_NAME,'Grid Statistics:'
      WRITE(STDOUT,'(A,5X,A,11X,I9)') SOLVER_NAME,'Vertices:           ', &
                                      pGrid%nVertTot                                     
      WRITE(STDOUT,'(A,5X,A,11X,I9)') SOLVER_NAME,'Cells:              ', &
                                      pGrid%nCellsTot
      WRITE(STDOUT,'(A,5X,A,11X,I9)') SOLVER_NAME,'Faces:              ', &
                                      gridCOBALT%nFaces
      WRITE(STDOUT,'(A,7X,A,9X,I9)')  SOLVER_NAME,'Triangular faces:   ', &
                                      gridCOBALT%nTris
      WRITE(STDOUT,'(A,7X,A,9X,I9)')  SOLVER_NAME,'Quadrilateral faces:', &
                                      gridCOBALT%nQuads
      WRITE(STDOUT,'(A,5X,A,11X,I9)') SOLVER_NAME,'Patches:            ', &
                                      gridCOBALT%nPatches      
    END IF ! global%verbLevel 

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading COBALT grid file done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)        

  END SUBROUTINE RFLU_ReadGridCOBALT









! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModCOBALT

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModCOBALT.F90,v $
! Revision 1.5  2008/12/06 08:45:03  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/28 19:48:58  haselbac
! Bug fix: Cell connectivity arrays were only dimensioned to Tot, not Max
!
! Revision 1.2  2006/03/25 22:04:29  haselbac
! Changes because of sype patches
!
! Revision 1.1  2005/04/15 15:09:08  haselbac
! Initial revision
!
! Revision 1.7  2005/01/20 14:54:56  haselbac
! Added setting of nBFaces and nBFacesTot
!
! Revision 1.6  2004/12/10 15:22:14  haselbac
! Bug fix: Added missing RFLU_SetMaxDimensions call
!
! Revision 1.5  2004/12/04 03:37:35  haselbac
! Adapted to changes in RFLU_ModCellMapping
!
! Revision 1.4  2004/11/03 17:09:24  haselbac
! Removed setting of vertex and cell flags
!
! Revision 1.3  2004/11/03 15:06:08  haselbac
! Bug fix: Corrected loop limits
!
! Revision 1.2  2004/10/19 19:31:08  haselbac
! Removed renumbering of bface lists
!
! Revision 1.1  2004/07/06 15:15:47  haselbac
! Initial revision
!
! ******************************************************************************








