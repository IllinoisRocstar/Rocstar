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
! Purpose: Suite of routines to build edge list.
!
! Description: None.
!
! Notes:
!   1. There is no routine called RFLU_CreateEdgeList - as might be expected - 
!      because the sizes of the various arrays is not known in advance and 
!      splitting the routine becomes cumbersome.
!
! ******************************************************************************
!
! $Id: RFLU_ModEdgeList.F90,v 1.13 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModEdgeList

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_NullifyEdgeList, & 
            RFLU_NullifyEdge2CellList, &          
            RFLU_CreateEdgeList, & 
            RFLU_CreateEdge2CellList, &  
            RFLU_BuildEdgeList, & 
            RFLU_BuildEdge2CellList, & 
            RFLU_DestroyEdgeList, &  
            RFLU_DestroyEdge2CellList
  
  SAVE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
    
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModEdgeList.F90,v $ $Revision: 1.13 $' 
  INTEGER, PRIVATE :: ekCntr(EDGE_KIND_AA:EDGE_KIND_VV), & 
                      ekOffs(EDGE_KIND_AA:EDGE_KIND_VV), & 
                      ekStrt(EDGE_KIND_AA:EDGE_KIND_VV)
  INTEGER, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: degr,strt                       
                                   
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  


! ******************************************************************************
!
! Purpose: Nullify edge list.
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
    
  SUBROUTINE RFLU_NullifyEdgeList(pRegion) 

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
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyEdgeList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Nullifying edge list...'         
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid  
  
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pGrid%e2v)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Nullifying edge list done.'
    END IF ! global%verbLevel
      
    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyEdgeList  





! ******************************************************************************
!
! Purpose: Nullify edge-to-cell list.
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
    
  SUBROUTINE RFLU_NullifyEdge2CellList(pRegion) 

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
    TYPE(t_global), POINTER :: global

! ******************************************************************************  
!   Start
! ******************************************************************************  

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyEdge2CellList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Nullifying edge-to-cell list...'         
    END IF ! global%verbLevel

! ******************************************************************************  
!   Set grid pointer
! ******************************************************************************  

    pGrid => pRegion%grid                               
  
! ******************************************************************************  
!   Nullify memory
! ******************************************************************************  

    NULLIFY(pGrid%e2cStrt)
    NULLIFY(pGrid%e2cDegr)
        
! ******************************************************************************  
!   End
! ******************************************************************************  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Nullifying edge-to-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyEdge2CellList   








! ******************************************************************************
!
! Purpose: Create edge list.
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

  SUBROUTINE RFLU_CreateEdgeList(pRegion) 

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

    INTEGER :: errorFlag,ieg,iPatch,nBFaces,nFacesEst
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateEdgeList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating edge list...'         
    END IF ! global%verbLevel

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyEdgeList(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid  
  
! ******************************************************************************
!   Estimate number of edges for allocation of hash table. The estimation 
!   formula shown below assumes that boundary effects are negligible, so for 
!   very small grids, where boundary edges dominate the total number of edges, 
!   need a kludge. Kludge is also needed if running in parallel and boundary
!   edges become a significant fraction of total number of edges, and when 
!   running code with periodic hack, where the missing boundaries mean that 
!   number of edges is underestimated.
! ******************************************************************************

    nBFaces = 0

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      nBFaces = nBFaces + pPatch%nBTrisTot + pPatch%nBQuadsTot       
    END DO ! iPatch

    nFacesEst  = nBFaces + 2*pGrid%nTetsTot + 3*pGrid%nHexsTot & 
               + 5*pGrid%nPrisTot/2

    pGrid%nEdgesEst = nBFaces + pGrid%nTetsTot + 2*pGrid%nHexsTot + & 
                      3*pGrid%nPrisTot/2 + pGrid%nVertTot - pGrid%nPyrsTot

    IF ( nBFaces/REAL(nFacesEst,KIND=RFREAL) > 0.8_RFREAL .OR. & 
         nBFaces/REAL(nFacesEst,KIND=RFREAL) < 0.2_RFREAL ) THEN 
      pGrid%nEdgesEst = 2*pGrid%nEdgesEst

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Corrected estimate of '// & 
                                 'number of edges.' 
      END IF ! global%verbLevel       
    END IF ! nBFaces

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A,3X,I9)') SOLVER_NAME,'Estimated number of '// & 
                                     'edges: ',pGrid%nEdgesEst        
    END IF ! global%verbLevel  

    ALLOCATE(pGrid%e2v(2,pGrid%nEdgesEst),STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%e2vTemp')
    END IF ! global%error  

    DO ieg = 1,pGrid%nEdgesEst ! Explicit loop because of ASCI White
      pGrid%e2v(1,ieg) = 0 
      pGrid%e2v(2,ieg) = 0               
    END DO ! ieg                        
      
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating edge list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateEdgeList  
  
  
  


! ******************************************************************************
!
! Purpose: Create edge-to-cell list.
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
  
  SUBROUTINE RFLU_CreateEdge2CellList(pRegion) 

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
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************  
!   Start
! ******************************************************************************  

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateEdge2CellList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating edge-to-cell list...'         
    END IF ! global%verbLevel

! ******************************************************************************  
!   Set grid pointer
! ******************************************************************************  

    pGrid => pRegion%grid                               
  
! ******************************************************************************  
!   Nullify memory
! ******************************************************************************  

    CALL RFLU_NullifyEdge2CellList(pRegion)

! ******************************************************************************  
!   Allocate memory
! ******************************************************************************  

    ALLOCATE(pGrid%e2cStrt(pGrid%nEdgesTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%e2cStrt')
    END IF ! global%error
      
    pGrid%e2cStrt(1) = 1 ! Initial value important

    ALLOCATE(pGrid%e2cDegr(pGrid%nEdgesTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%e2cDegr')
    END IF ! global%error     
        
! ******************************************************************************  
!   End
! ******************************************************************************  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating edge-to-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateEdge2CellList      
  
  





! ******************************************************************************
!
! Purpose: Build edge list.
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
      
  SUBROUTINE RFLU_BuildEdgeList(pRegion) 

    USE ModSortSearch

    USE RFLU_ModCellFaceEdgeInfo
    USE RFLU_ModGrid
    USE RFLU_ModHashTable

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

    INTEGER :: eCntr,edgeType,ekSum,errorFlag,ieg,iegb,iege,iek,iel,icl, & 
               ivg,key,nEdges,vSize,v1,v2
    INTEGER :: v(2)

    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildEdgeList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building edge list...' 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building hash table...'         
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer and initialize nEdgesTot
! ******************************************************************************

    pGrid => pRegion%grid

    pGrid%nEdgesTot = 0

! ******************************************************************************
!   Create hash table
! ******************************************************************************

    CALL RFLU_CreateHashTable(global,pGrid%nEdgesEst)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN                            
      WRITE(STDOUT,'(A,5X,A,1X,I9)') SOLVER_NAME,'Hash table size:    '// &
                                     '       ',hashTableSize
    END IF ! global%verbLevel

! ******************************************************************************
!   Loop over cell types, construct hash table of faces
! ******************************************************************************

    vSize  = 2 ! Edges have only two vertices
    nEdges = 0

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN                              
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Looping over cell types...'
    END IF ! global%verbLevel  

! ==============================================================================
!   Tetrahedra
! ==============================================================================

    IF ( pGrid%nTetsTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel
    END IF ! pGrid%nTetsTot

    DO icl = 1,pGrid%nTetsTot            
      DO iel = 1,6    
        v(1) = pGrid%tet2v(ce2vTet(1,iel),icl)
        v(2) = pGrid%tet2v(ce2vTet(2,iel),icl)         

        CALL QuickSortInteger(v(1:vSize),vSize)
        CALL RFLU_HashBuildKey(v(1:2),2,key)              
        CALL RFLU_HashEdge(global,key,pGrid,v(1:2),edgeType) 

        IF ( edgeType == EDGE_TYPE_NEW ) THEN 
          nEdges = nEdges + 1
        END IF ! edgeType
      END DO ! iel
    END DO ! icl

! ==============================================================================
!   Hexahedra
! ==============================================================================

    IF ( pGrid%nHexsTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Hexahedra...'
      END IF ! global%verbLevel
    END IF ! pGrid%nHexsTot

    DO icl = 1,pGrid%nHexsTot
      DO iel = 1,12
        v(1) = pGrid%hex2v(ce2vHex(1,iel),icl)
        v(2) = pGrid%hex2v(ce2vHex(2,iel),icl)   

        CALL QuickSortInteger(v(1:vSize),vSize)
        CALL RFLU_HashBuildKey(v(1:2),2,key)    
        CALL RFLU_HashEdge(global,key,pGrid,v(1:2),edgeType)  

        IF ( edgeType == EDGE_TYPE_NEW ) THEN 
          nEdges = nEdges + 1
        END IF ! edgeType
      END DO ! ifl
    END DO ! icl

! ==============================================================================
!   Prisms
! ==============================================================================

    IF ( pGrid%nPrisTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Prisms...'
      END IF ! global%verbLevel
    END IF ! pGrid%nPrisTot

    DO icl = 1,pGrid%nPrisTot
      DO iel = 1,9
        v(1) = pGrid%pri2v(ce2vPri(1,iel),icl) 
        v(2) = pGrid%pri2v(ce2vPri(2,iel),icl)     

        CALL QuickSortInteger(v(1:vSize),vSize)
        CALL RFLU_HashBuildKey(v(1:2),2,key)    
        CALL RFLU_HashEdge(global,key,pGrid,v(1:2),edgeType)  

        IF ( edgeType == EDGE_TYPE_NEW ) THEN 
          nEdges = nEdges + 1
        END IF ! edgeType    
      END DO ! ifl
    END DO ! icl  

! ==============================================================================
!   Pyramids
! ==============================================================================

    IF ( pGrid%nPyrsTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Pyramids...'
      END IF ! global%verbLevel
    END IF ! pGrid%nPyrsTot

    DO icl = 1,pGrid%nPyrsTot
      DO iel = 1,8
        v(1) = pGrid%pyr2v(ce2vPyr(1,iel),icl) 
        v(2) = pGrid%pyr2v(ce2vPyr(2,iel),icl)    

        CALL QuickSortInteger(v(1:vSize),vSize)
        CALL RFLU_HashBuildKey(v(1:2),2,key)    
        CALL RFLU_HashEdge(global,key,pGrid,v(1:2),edgeType)  

        IF ( edgeType == EDGE_TYPE_NEW ) THEN 
          nEdges = nEdges + 1
        END IF ! edgeType    
      END DO ! ifl
    END DO ! icl    

! ******************************************************************************
!   Set number of edges, print some info, and destroy hash table
! ******************************************************************************

    pGrid%nEdgesTot = nEdges

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN
      WRITE(STDOUT,'(A,5X,A,6X,I9)') SOLVER_NAME,'Hash table collisions:', & 
                                     hashTableCollisions  
    END IF ! global%verbLevel

    CALL RFLU_DestroyHashTable(global)

! ******************************************************************************
!   Determine edge kinds and set counter offsets
! ******************************************************************************

    ekCntr(EDGE_KIND_AA) = 0
    ekCntr(EDGE_KIND_AV) = 0
    ekCntr(EDGE_KIND_VV) = 0    

    DO ieg = 1,pGrid%nEdgesTot
      v1 = pGrid%e2v(1,ieg)
      v2 = pGrid%e2v(2,ieg)

      SELECT CASE (RFLU_GetEdgeKind(global,pGrid,v1,v2))         
        CASE (EDGE_KIND_AA) 
          ekCntr(EDGE_KIND_AA) = ekCntr(EDGE_KIND_AA) + 1
        CASE (EDGE_KIND_AV)
          ekCntr(EDGE_KIND_AV) = ekCntr(EDGE_KIND_AV) + 1
        CASE (EDGE_KIND_VV) 
          ekCntr(EDGE_KIND_VV) = ekCntr(EDGE_KIND_VV) + 1          
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)     
      END SELECT ! RFLU_GetEdgeKind      
    END DO ! ie        

    ekOffs(EDGE_KIND_AA) = 0
    ekOffs(EDGE_KIND_AV) = ekCntr(EDGE_KIND_AA)
    ekOffs(EDGE_KIND_VV) = ekCntr(EDGE_KIND_AV) + ekOffs(EDGE_KIND_AV)

! ==============================================================================
!   Set nEdges. NOTE that it contains only those edges which link two actual
!   vertices. This is done for parallel grid motion, in which each region only
!   smooths by looping over the edges linking actual vertices. 
! ==============================================================================

!    pGrid%nEdges = ekCntr(EDGE_KIND_AA) + ekCntr(EDGE_KIND_AV)
    pGrid%nEdges = ekCntr(EDGE_KIND_AA)

! ==============================================================================
!   Print info 
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      !IF ( global%verbLevel >= VERBOSE_LOW ) THEN   
        WRITE(STDOUT,'(A,5X,A)')       SOLVER_NAME,'Edge-type statistics:' 
        WRITE(STDOUT,'(A,7X,A,4X,I9)') SOLVER_NAME,'Total edges:          ', & 
                                       pGrid%nEdgesTot        
        WRITE(STDOUT,'(A,5X,A)')       SOLVER_NAME,'Edge-kind statistics:' 
        WRITE(STDOUT,'(A,7X,A,4X,I9)') SOLVER_NAME,'Actual-actual edges:  ', & 
                                       ekCntr(EDGE_KIND_AA)
        WRITE(STDOUT,'(A,7X,A,4X,I9)') SOLVER_NAME,'Actual-virtual edges: ', & 
                                       ekCntr(EDGE_KIND_AV)
        WRITE(STDOUT,'(A,7X,A,4X,I9)') SOLVER_NAME,'Virtual-virtual edges:', & 
                                       ekCntr(EDGE_KIND_VV)
      !END IF ! global%verbLevel
    END IF ! global%myProcid     

! ******************************************************************************
!   Build actual edgelist
! ******************************************************************************

! ==============================================================================
!   Allocate and build helper arrays for sorting. The helper arrays are needed
!   so that within each edge kind, have edges sorted according to increasing
!   origination vertices and increasing destination vertices.
! ==============================================================================

    ALLOCATE(strt(EDGE_KIND_AA:EDGE_KIND_VV,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'strt')
    END IF ! global%error

    ALLOCATE(degr(EDGE_KIND_AA:EDGE_KIND_VV,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'degr')
    END IF ! global%error

    DO ivg = 1,pGrid%nVertTot ! Explicit loop because of ASCI White problems
      degr(EDGE_KIND_AA,ivg) = 0
      degr(EDGE_KIND_AV,ivg) = 0
      degr(EDGE_KIND_VV,ivg) = 0        
    END DO ! ivg

    DO ieg = 1,pGrid%nEdgesTot
      v1 = pGrid%e2v(1,ieg)
      v2 = pGrid%e2v(2,ieg)

      iek = RFLU_GetEdgeKind(global,pGrid,v1,v2)

      degr(iek,v1) = degr(iek,v1) + 1
    END DO ! ieg      

    strt(EDGE_KIND_AA:EDGE_KIND_VV,1) = 1

    DO ivg = 2,pGrid%nVertTot
      strt(EDGE_KIND_AA,ivg) = strt(EDGE_KIND_AA,ivg-1) & 
                             + degr(EDGE_KIND_AA,ivg-1)
      strt(EDGE_KIND_AV,ivg) = strt(EDGE_KIND_AV,ivg-1) & 
                             + degr(EDGE_KIND_AV,ivg-1)
      strt(EDGE_KIND_VV,ivg) = strt(EDGE_KIND_VV,ivg-1) & 
                             + degr(EDGE_KIND_VV,ivg-1)
    END DO ! ivg

! ==============================================================================
!   Transfer provisional edge-array into temporary array. Sort so edge kinds 
!   actual-actual and actual-virtual are listed before virtual-virtual ones. 
!   Within each edge-kind segment, sort according to increasing origination
!   vertices. 
! ==============================================================================

    ekCntr(EDGE_KIND_AA) = 0 ! Reinitialize
    ekCntr(EDGE_KIND_AV) = 0
    ekCntr(EDGE_KIND_VV) = 0    

    DO ivg = 1,pGrid%nVertTot ! Explicit loop because of ASCI White 
      degr(EDGE_KIND_AA,ivg) = 0
      degr(EDGE_KIND_AV,ivg) = 0
      degr(EDGE_KIND_VV,ivg) = 0        
    END DO ! ivg

    ALLOCATE(pGrid%e2vTemp(2,pGrid%nEdgesTot),STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%e2vTemp')
    END IF ! global%error  

    DO ieg = 1,pGrid%nEdgesTot
      v1 = pGrid%e2v(1,ieg)
      v2 = pGrid%e2v(2,ieg)

      iek = RFLU_GetEdgeKind(global,pGrid,v1,v2)

      eCntr = strt(iek,v1) + degr(iek,v1) + ekOffs(iek)        
      degr(iek,v1) = degr(iek,v1) + 1        
      ekCntr(iek)  = ekCntr(iek)  + 1 ! used below for consistency check     

      pGrid%e2vTemp(1,eCntr) = pGrid%e2v(1,ieg)
      pGrid%e2vTemp(2,eCntr) = pGrid%e2v(2,ieg) 
    END DO ! ieg

    ekSum = ekCntr(EDGE_KIND_AA) + ekCntr(EDGE_KIND_AV) + ekCntr(EDGE_KIND_VV)

    IF ( ekSum /= pGrid%nEdgesTot ) THEN ! Consistency check
      CALL ErrorStop(global,ERR_NEDGES_WRONG,__LINE__) 
    END IF ! ekSum

! ==============================================================================
!   Finally sort within each edge-kind segment and each origination vertex
!   segment so that destination vertices are listed in increasing order. Note
!   because edges within each segment are already sorted according to 
!   origination vertices, need to sort here only according to destination
!   vertices for each range of origination vertices.
! ==============================================================================

    DO ivg = 1,pGrid%nVertTot
      
! ------------------------------------------------------------------------------
!     Actual-actual edges
! ------------------------------------------------------------------------------
      
      iek  = EDGE_KIND_AA      
      iegb = ekOffs(iek) + strt(iek,ivg)                      
      iege = ekOffs(iek) + strt(iek,ivg) + degr(iek,ivg) - 1          

      IF ( iege > iegb ) THEN 
        CALL QuickSortInteger(pGrid%e2vTemp(2,iegb:iege),iege-iegb+1)
      END IF ! iege

! ------------------------------------------------------------------------------
!     Actual-virtual edges 
! ------------------------------------------------------------------------------
      
      iek  = EDGE_KIND_AV
      iegb = ekOffs(iek) + strt(iek,ivg)                      
      iege = ekOffs(iek) + strt(iek,ivg) + degr(iek,ivg) - 1   

      IF ( iege > iegb ) THEN 
        CALL QuickSortInteger(pGrid%e2vTemp(2,iegb:iege),iege-iegb+1)
      END IF ! iege
        
! ------------------------------------------------------------------------------
!      Virtual-virtual edges        
! ------------------------------------------------------------------------------
        
      iek  = EDGE_KIND_VV
      iegb = ekOffs(iek) + strt(iek,ivg)                      
      iege = ekOffs(iek) + strt(iek,ivg) + degr(iek,ivg) - 1   

      IF ( iege > iegb ) THEN 
        CALL QuickSortInteger(pGrid%e2vTemp(2,iegb:iege),iege-iegb+1)
      END IF ! iege                            
    END DO ! ivg

! ==============================================================================
!   Deallocate helper arrays
! ==============================================================================

    DEALLOCATE(strt,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'strt')
    END IF ! global%error

    DEALLOCATE(degr,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'degr')
    END IF ! global%error

! ==============================================================================
!   Allocate actual edge array and copy from temporary array
! ==============================================================================

    DEALLOCATE(pGrid%e2v,STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%e2v')
    END IF ! global%error  

    ALLOCATE(pGrid%e2v(2,pGrid%nEdgesTot),STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%e2v')
    END IF ! global%error  

    DO ieg = 1,pGrid%nEdgesTot ! Explicit loop because of ASCI White
      pGrid%e2v(1,ieg) = pGrid%e2vTemp(1,ieg)
      pGrid%e2v(2,ieg) = pGrid%e2vTemp(2,ieg)
    END DO ! ieg

! ==============================================================================
!   Deallocate temporary array
! ==============================================================================

    DEALLOCATE(pGrid%e2vTemp,STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%e2v')
    END IF ! global%error  
      
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building edge list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildEdgeList
  
  
  






! ******************************************************************************
!
! Purpose: Build edge-to-cell list.
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
  
  SUBROUTINE RFLU_BuildEdge2CellList(pRegion)
      
    USE ModSortSearch

    USE RFLU_ModCellFaceEdgeInfo
    USE RFLU_ModGrid      
      
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
  
    INTEGER :: errorFlag,e2cSize,icl,ieg,iegb,iege,iek,iel,iloc,iPass, & 
               ivg,vSize,v1,v2
    INTEGER :: v(2)
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************  
!   Start
! ******************************************************************************  
      
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildEdge2CellList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building edge-to-cell list...'         
    END IF ! global%verbLevel      

    pGrid => pRegion%grid
      
! ******************************************************************************  
!   Build access lists for edge kinds, allows easy access of e2v list based on
!   edge kinds
! ******************************************************************************  
      
    ekCntr(EDGE_KIND_AA) = 0
    ekCntr(EDGE_KIND_AV) = 0
    ekCntr(EDGE_KIND_VV) = 0    

    DO ieg = 1,pGrid%nEdgesTot
      v1 = pGrid%e2v(1,ieg)
      v2 = pGrid%e2v(2,ieg)

      SELECT CASE (RFLU_GetEdgeKind(global,pGrid,v1,v2))         
        CASE (EDGE_KIND_AA) 
          ekCntr(EDGE_KIND_AA) = ekCntr(EDGE_KIND_AA) + 1
        CASE (EDGE_KIND_AV)
          ekCntr(EDGE_KIND_AV) = ekCntr(EDGE_KIND_AV) + 1
        CASE (EDGE_KIND_VV) 
          ekCntr(EDGE_KIND_VV) = ekCntr(EDGE_KIND_VV) + 1          
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)     
      END SELECT ! RFLU_GetEdgeKind      
    END DO ! ie        

    ekOffs(EDGE_KIND_AA) = 0
    ekOffs(EDGE_KIND_AV) = ekCntr(EDGE_KIND_AA)
    ekOffs(EDGE_KIND_VV) = ekCntr(EDGE_KIND_AV) + ekOffs(EDGE_KIND_AV)

    ekStrt(EDGE_KIND_AA) = 1
    ekStrt(EDGE_KIND_AV) = ekStrt(EDGE_KIND_AA) + ekCntr(EDGE_KIND_AA)
    ekStrt(EDGE_KIND_VV) = ekStrt(EDGE_KIND_AV) + ekCntr(EDGE_KIND_AV)     
      
! ******************************************************************************  
!   Allocate and build helper arrays for sorting. The helper arrays are needed
!   so that within each edge kind, have edges sorted according to increasing
!   origination vertices and increasing destination vertices.
! ******************************************************************************  

    ALLOCATE(strt(EDGE_KIND_AA:EDGE_KIND_VV,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'strt')
    END IF ! global%error

    ALLOCATE(degr(EDGE_KIND_AA:EDGE_KIND_VV,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'degr')
    END IF ! global%error

    DO ivg = 1,pGrid%nVertTot ! Explicit loop because of ASCI White
      degr(EDGE_KIND_AA,ivg) = 0
      degr(EDGE_KIND_AV,ivg) = 0
      degr(EDGE_KIND_VV,ivg) = 0        
    END DO ! ivg

    DO ieg = 1,pGrid%nEdgesTot
      v1 = pGrid%e2v(1,ieg)
      v2 = pGrid%e2v(2,ieg)

      iek = RFLU_GetEdgeKind(global,pGrid,v1,v2)

      degr(iek,v1) = degr(iek,v1) + 1
    END DO ! ieg      

    strt(EDGE_KIND_AA,1) = 1
    strt(EDGE_KIND_AV,1) = 1
    strt(EDGE_KIND_VV,1) = 1

    DO ivg = 2,pGrid%nVertTot
      strt(EDGE_KIND_AA,ivg) = strt(EDGE_KIND_AA,ivg-1) & 
                             + degr(EDGE_KIND_AA,ivg-1)
      strt(EDGE_KIND_AV,ivg) = strt(EDGE_KIND_AV,ivg-1) & 
                             + degr(EDGE_KIND_AV,ivg-1)
      strt(EDGE_KIND_VV,ivg) = strt(EDGE_KIND_VV,ivg-1) & 
                             + degr(EDGE_KIND_VV,ivg-1)
    END DO ! ivg

! ******************************************************************************  
!   Loop over passes. In first pass, determine how many cells meet at a given
!   edge, and use this information to allocate the properly dimensioned e2c
!   array. In the second pass, fill the e2c array.
! ******************************************************************************  

    DO iPass = 1,2
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH) THEN   
        WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME,'Pass:',iPass
      END IF ! global%verbLevel

! ==============================================================================
!     Loop over cell types
! ==============================================================================

      vSize = 2

      DO ieg = 1,pGrid%nEdgesTot ! Explicit loop because of ASCI White
        pGrid%e2cDegr(ieg) = 0
      END DO ! ieg

! ------------------------------------------------------------------------------
!     Tetrahedra
! ------------------------------------------------------------------------------
    
      IF ( pGrid%nTetsTot /= 0 ) THEN
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel >= VERBOSE_HIGH) THEN   
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Tetrahedra...'
        END IF ! global%verbLevel
      END IF ! pGrid%nTetsTot

      DO icl = 1,pGrid%nTetsTot            
        DO iel = 1,6    
          v(1) = pGrid%tet2v(ce2vTet(1,iel),icl)
          v(2) = pGrid%tet2v(ce2vTet(2,iel),icl)         

          CALL QuickSortInteger(v(1:vSize),vSize)

          iek  = RFLU_GetEdgeKind(global,pGrid,v(1),v(2))          
          iegb = ekOffs(iek) + strt(iek,v(1))                      
          iege = ekOffs(iek) + strt(iek,v(1)) + degr(iek,v(1)) - 1    

          CALL BinarySearchInteger(pGrid%e2v(2,iegb:iege),iege-iegb+1,v(2), &
                                   iloc)

          IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
            ieg = iloc + iegb - 1

            IF ( pGrid%e2v(1,ieg) /= v(1) .OR. pGrid%e2v(2,ieg) /= v(2) ) THEN
              CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)
            END IF ! pGrid%e2v

            pGrid%e2cDegr(ieg) = pGrid%e2cDegr(ieg) + 1

            IF ( iPass == 2 ) THEN 
              iloc = pGrid%e2cStrt(ieg) + pGrid%e2cDegr(ieg) - 1                            
              pGrid%e2c(iloc) = pGrid%tet2CellGlob(icl)
            END IF ! iPass   
          ELSE 
            CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)                
          END IF ! iloc
        END DO ! iel
      END DO ! icl      
    
! ------------------------------------------------------------------------------
!     Hexahedra 
! ------------------------------------------------------------------------------
    
      IF ( pGrid%nHexsTot /= 0 ) THEN
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel >= VERBOSE_HIGH) THEN   
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Hexahedra...'
        END IF ! global%verbLevel
      END IF ! pGrid%nHexsTot

      DO icl = 1,pGrid%nHexsTot            
        DO iel = 1,12    
          v(1) = pGrid%hex2v(ce2vHex(1,iel),icl)
          v(2) = pGrid%hex2v(ce2vHex(2,iel),icl)         

          CALL QuickSortInteger(v(1:vSize),vSize)

          iek  = RFLU_GetEdgeKind(global,pGrid,v(1),v(2))         
          iegb = ekOffs(iek) + strt(iek,v(1))                      
          iege = ekOffs(iek) + strt(iek,v(1)) + degr(iek,v(1)) - 1    

          CALL BinarySearchInteger(pGrid%e2v(2,iegb:iege),iege-iegb+1,v(2), &
                                   iloc)

          IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
            ieg = iloc + iegb - 1

            IF ( pGrid%e2v(1,ieg) /= v(1) .OR. pGrid%e2v(2,ieg) /= v(2) ) THEN
              CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)
            END IF ! pGrid%e2v

            pGrid%e2cDegr(ieg) = pGrid%e2cDegr(ieg) + 1

            IF ( iPass == 2 ) THEN               
              iloc = pGrid%e2cStrt(ieg) + pGrid%e2cDegr(ieg) - 1                            
              pGrid%e2c(iloc) = pGrid%hex2CellGlob(icl)
            END IF ! iPass                               
          ELSE 
            CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)                
          END IF ! iloc
        END DO ! iel
      END DO ! icl          

! ------------------------------------------------------------------------------
!     Prisms
! ------------------------------------------------------------------------------
    
      IF ( pGrid%nPrisTot /= 0 ) THEN
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel >= VERBOSE_HIGH) THEN   
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Prisms...'
        END IF ! global%verbLevel
      END IF ! pGrid%nPrisTot

      DO icl = 1,pGrid%nPrisTot            
        DO iel = 1,9    
          v(1) = pGrid%pri2v(ce2vPri(1,iel),icl)
          v(2) = pGrid%pri2v(ce2vPri(2,iel),icl)         

          CALL QuickSortInteger(v(1:vSize),vSize)

          iek  = RFLU_GetEdgeKind(global,pGrid,v(1),v(2))         
          iegb = ekOffs(iek) + strt(iek,v(1))                      
          iege = ekOffs(iek) + strt(iek,v(1)) + degr(iek,v(1)) - 1    

          CALL BinarySearchInteger(pGrid%e2v(2,iegb:iege),iege-iegb+1,v(2), &
                                   iloc)

          IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
            ieg = iloc + iegb - 1

            IF ( pGrid%e2v(1,ieg) /= v(1) .OR. pGrid%e2v(2,ieg) /= v(2) ) THEN
              CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)
            END IF ! pGrid%e2v

            pGrid%e2cDegr(ieg) = pGrid%e2cDegr(ieg) + 1              

            IF ( iPass == 2 ) THEN               
              iloc = pGrid%e2cStrt(ieg) + pGrid%e2cDegr(ieg) - 1                            
              pGrid%e2c(iloc) = pGrid%pri2CellGlob(icl)
            END IF ! iPass            
          ELSE 
            CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)                
          END IF ! iloc
        END DO ! iel
      END DO ! icl     
    
! ------------------------------------------------------------------------------
!     Pyramids 
! ------------------------------------------------------------------------------
    
      IF ( pGrid%nPyrsTot /= 0 ) THEN
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel >= VERBOSE_HIGH) THEN   
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Pyramids...'
        END IF ! global%verbLevel
      END IF ! pGrid%nPyrTot

      DO icl = 1,pGrid%nPyrsTot            
        DO iel = 1,8
          v(1) = pGrid%pyr2v(ce2vPyr(1,iel),icl)
          v(2) = pGrid%pyr2v(ce2vPyr(2,iel),icl)         

          CALL QuickSortInteger(v(1:vSize),vSize)

          iek  = RFLU_GetEdgeKind(global,pGrid,v(1),v(2))         
          iegb = ekOffs(iek) + strt(iek,v(1))                      
          iege = ekOffs(iek) + strt(iek,v(1)) + degr(iek,v(1)) - 1    

          CALL BinarySearchInteger(pGrid%e2v(2,iegb:iege),iege-iegb+1,v(2), & 
                                   iloc)

          IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
            ieg = iloc + iegb - 1

            IF ( pGrid%e2v(1,ieg) /= v(1) .OR. pGrid%e2v(2,ieg) /= v(2) ) THEN
              CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)
            END IF ! pGrid%e2v

            pGrid%e2cDegr(ieg) = pGrid%e2cDegr(ieg) + 1              

            IF ( iPass == 2 ) THEN               
              iloc = pGrid%e2cStrt(ieg) + pGrid%e2cDegr(ieg) - 1                            
              pGrid%e2c(iloc) = pGrid%pyr2CellGlob(icl)
            END IF ! iPass
          ELSE 
            CALL ErrorStop(global,ERR_EDGELIST_INVALID,__LINE__)                
          END IF ! iloc
        END DO ! iel
      END DO ! icl         
    
! ==============================================================================
!     Determine total size of e2c array, build e2cStrt array
! ==============================================================================
      
      IF ( iPass == 1 ) THEN
        e2cSize = 0

        DO ieg = 1,pGrid%nEdgesTot
          e2cSize = e2cSize + pGrid%e2cDegr(ieg)

          IF ( ieg > 1 ) THEN 
            pGrid%e2cStrt(ieg) = pGrid%e2cStrt(ieg-1) + pGrid%e2cDegr(ieg-1)
          END IF ! ieg
        END DO ! ieg

        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel >= VERBOSE_HIGH) THEN 
          WRITE(STDOUT,'(A,5X,A,1X,I9)') SOLVER_NAME,'Size:',e2cSize
        END IF ! global%myProcid

        ALLOCATE(pGrid%e2c(e2cSize),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%e2c')
        END IF ! global%error

        pGrid%e2c(1:e2cSize) = 0                   
      END IF ! iPass    
    END DO ! iPass

! ******************************************************************************  
!   Deallocate helper arrays
! ******************************************************************************  
    
    DEALLOCATE(strt,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'strt')
    END IF ! global%error

    DEALLOCATE(degr,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'degr')
    END IF ! global%error    
    
! ******************************************************************************  
!   End 
! ******************************************************************************  
    
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building edge-to-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)    

  END SUBROUTINE RFLU_BuildEdge2CellList

  
  
     
  
  


! ******************************************************************************
!
! Purpose: Destroy edge list.
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
  
  SUBROUTINE RFLU_DestroyEdgeList(pRegion) 

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
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_global), POINTER :: global

! ******************************************************************************  
!   Start
! ******************************************************************************  

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyEdgeList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying edge list...'         
    END IF ! global%verbLevel

! ******************************************************************************  
!   Set grid pointer
! ******************************************************************************  

    pGrid => pRegion%grid  
  
! ******************************************************************************  
!   Deallocate memory
! ******************************************************************************  

    DEALLOCATE(pGrid%e2v,STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%e2v')
    END IF ! global%error  
        
! ******************************************************************************  
!   Nullify memory
! ******************************************************************************  

    CALL RFLU_NullifyEdgeList(pRegion)

! ******************************************************************************  
!   End
! ******************************************************************************  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying edge list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyEdgeList     
  





! ******************************************************************************
!
! Purpose: Destroy edge-to-cell list.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Test for association status of e2c because it is not created with other 
!      arrays in RFLU_CreateEdge2CellList, but simply allocated in 
!      RFLU_BuildEdge2CellList. This is done because the size is not known at 
!      time of creation.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_DestroyEdge2CellList(pRegion) 

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
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_global), POINTER :: global

! ==============================================================================
!     Start
! ==============================================================================

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyEdge2CellList',&
  'RFLU_ModEdgeList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying edge-to-cell list...'         
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid  
  
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(pGrid%e2cDegr,STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%e2cDegr')
    END IF ! global%error

    DEALLOCATE(pGrid%e2cStrt,STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%e2cStrt')
    END IF ! global%error       

    IF ( ASSOCIATED(pGrid%e2c) .EQV. .TRUE. ) THEN 
      DEALLOCATE(pGrid%e2c,STAT=errorFlag) 
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%e2c')
      END IF ! global%error                 
    END IF ! ASSOCIATED
       
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyEdge2CellList(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying edge-to-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyEdge2CellList    





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModEdgeList


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModEdgeList.F90,v $
! Revision 1.13  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.12  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.11  2005/04/15 15:06:51  haselbac
! Removed Charm/FEM stuff and RFLU_XyzEdge2RegionDegrList routines
!
! Revision 1.10  2004/10/19 19:27:52  haselbac
! Substantial clean-up
!
! Revision 1.9  2004/07/06 15:14:36  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!                                                         
! Revision 1.8  2004/01/22 16:03:59  haselbac                                            
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC 
! and titan  
!
! Revision 1.7  2004/01/11 02:17:35  jiao                                                
! Eliminated some redundant trailing spaces that made some lines too long.               
! This changed was necessary to compile with NAG F90 compiler.                           
!
! Revision 1.6  2003/08/18 15:30:36  haselbac                                            
! Changed initialization to solve possible problem on Frost                              
!
! Revision 1.5  2003/06/04 22:08:30  haselbac                                            
! Added Nullify routines, some cosmetics                                                 
!
! Revision 1.4  2003/04/02 17:27:35  haselbac                                            
! Changed limits of modified estimate of no of edges                                     
!
! Revision 1.3  2003/03/15 18:05:15  haselbac                                            
! New routines (|| gm), deleted params, bug fix for pyrs                                 
!
! Revision 1.2  2003/01/28 16:27:58  haselbac                                            
! Added creation and destruction, removed renumbering (bcos of                           
! RFLU_InitFlowSolver changes)                                                           
!
! Revision 1.1  2002/10/27 19:08:52  haselbac                                            
! Initial revision                                                                       
!
! ******************************************************************************














