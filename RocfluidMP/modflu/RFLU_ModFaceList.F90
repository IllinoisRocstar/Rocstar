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
! Purpose: Suite of routines to build face list.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModFaceList.F90,v 1.43 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModFaceList

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModBorder, ONLY: t_border
  USE ModMPI

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_BuildAVFace2BorderList, &
            RFLU_BuildAVFace2PatchList, &  
            RFLU_BuildCell2FaceList, &
            RFLU_BuildFaceList, & 
            RFLU_CreateAVFace2BorderList, &
            RFLU_CreateAVFace2PatchList, & 
            RFLU_CreateCell2FaceList, &   
            RFLU_CreateFaceList, & 
            RFLU_DestroyAVFace2BorderList, &
            RFLU_DestroyAVFace2PatchList, &            
            RFLU_DestroyCell2FaceList, & 
            RFLU_DestroyFaceList, &
            RFLU_GetOpposingFaces, & 
            RFLU_InsertIntoCell2FaceList, &           
            RFLU_ReorientFaces

  SAVE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  INTEGER, PARAMETER, PUBLIC :: DESTROY_FACE_INT = 1, & 
                                DESTROY_FACE_EXT = 2, & 
                                DESTROY_FACE_ALL = 4     
    
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModFaceList.F90,v $ $Revision: 1.43 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  






! ******************************************************************************
!
! Purpose: Build actual-virtual-face-to-border list.
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
    
  SUBROUTINE RFLU_BuildAVFace2BorderList(pRegion)

    USE ModSortSearch

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
      
    INTEGER :: c1,c2,errorFlag,iBorder,icg,icgMax,icl,ifg,ifl,iLoc, &
               nCellsRecvTot
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: vc,vc2b
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildAVFace2BorderList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building av-face-to-border list...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Build temporary lists
! ******************************************************************************

! ==============================================================================
!   Determine maximum cell index of actual-virtual faces so can limit searching  
! ==============================================================================

    icgMax = 0

    DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)
      
      IF ( c1 > pGrid%nCells ) THEN ! c1 is virtual cell
        icgMax = MAX(c1,icgMax)
      ELSE IF ( c2 > pGrid%nCells ) THEN ! c2 is virtual cell
        icgMax = MAX(c2,icgMax)
      ELSE ! Error
! TEMPORARY      
        WRITE(*,*) 'ERROR 1!'
        STOP
! END TEMPORARY        
      END IF ! c1
    END DO ! ifg

! ==============================================================================
!   Determine number of border cells whose global cell indices are smaller than 
!   maximum cell index determined above and allocate temporary memory 
! ==============================================================================

    nCellsRecvTot = 0

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
      
      DO icl = 1,pBorder%nCellsRecv
        IF ( pBorder%icgRecv(icl) <= icgMax ) THEN 
          nCellsRecvTot = nCellsRecvTot + 1
        END IF ! pBorder%icgRecv
      END DO ! icl
    END DO ! iBorder
    
    ALLOCATE(vc(2,nCellsRecvTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vc')
    END IF ! global%error

    ALLOCATE(vc2b(2,nCellsRecvTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vc2b')
    END IF ! global%error    
    
! ==============================================================================
!   Build temporary lists and sort to allow searching 
! ==============================================================================

    nCellsRecvTot = 0 ! Reinitialize
    
    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
      
      DO icl = 1,pBorder%nCellsRecv
        IF ( pBorder%icgRecv(icl) <= icgMax ) THEN 
          nCellsRecvTot = nCellsRecvTot + 1
          
          vc(1,nCellsRecvTot) = pBorder%icgRecv(icl)
          vc(2,nCellsRecvTot) = nCellsRecvTot
          
          vc2b(1,nCellsRecvTot) = iBorder
          vc2b(2,nCellsRecvTot) = icl
        END IF ! pBorder%icgRecv
      END DO ! icl
    END DO ! iBorder    
    
    CALL QuickSortIntegerInteger(vc(1:1,1:nCellsRecvTot), &
                                 vc(2:2,1:nCellsRecvTot), &
                                 nCellsRecvTot)

! ******************************************************************************
!   Build 
! ******************************************************************************
               
    DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)
      
      IF ( c1 > pGrid%nCells ) THEN ! c1 is virtual cell
        icg = c1
      ELSE IF ( c2 > pGrid%nCells ) THEN ! c2 is virtual cell
        icg = c2
      ELSE ! Error
! TEMPORARY      
        WRITE(*,*) 'ERROR 2!'
        STOP
! END TEMPORARY        
      END IF ! c1

      CALL BinarySearchInteger(vc(1:1,1:nCellsRecvTot),nCellsRecvTot,icg,iLoc)
      
      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        icl = vc(2,iLoc)

        ifl = ifg - (pGrid%nFaces-pGrid%nFacesAV)

        pGrid%avf2b(1,ifl) = vc2b(1,icl)
        pGrid%avf2b(2,ifl) = vc2b(2,icl)

      ELSE 
! TEMPORARY
        WRITE(*,*) 'ERROR 3!'
        STOP
! END TEMPORARY      
      END IF ! iLoc
    END DO ! ifg                                 

! ==============================================================================
!   Deallocate temporary memory 
! ==============================================================================
 
    DEALLOCATE(vc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vc')
    END IF ! global%error

    DEALLOCATE(vc2b,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vc2b')
    END IF ! global%error          
        
! ******************************************************************************
!   End 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building av-face-to-border list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildAVFace2BorderList








! ******************************************************************************
!
! Purpose: Build actual-virtual-face-to-patch list.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Acess avf2b lists, so these must be built first.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_BuildAVFace2PatchList(pRegion)

    USE ModSortSearch

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
      
    INTEGER :: errorFlag,iBorder,icg,icl,ifg,ifl,iLoc,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildAVFace2PatchList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building av-face-to-patch list...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Build temporary lists of sorted virtual cells adjacent to patches
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( pPatch%nBCellsVirt > 0 ) THEN 
        ALLOCATE(pPatch%bvcSorted(pPatch%nBCellsVirt),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvcSorted')
        END IF ! global%error

        DO icl = 1,pPatch%nBCellsVirt
          pPatch%bvcSorted(icl) = pPatch%bvc(icl)
        END DO ! icl

        CALL QuickSortInteger(pPatch%bvcSorted,pPatch%nBCellsVirt)
      ELSE 
        NULLIFY(pPatch%bvcSorted)
      END IF ! pPatch%nBCellsVirt
    END DO ! iPatch

! ******************************************************************************
!   Build lists of patch indices for each av face
! ******************************************************************************

    DO ifl = 1,pGrid%nFacesAV
      iBorder =  pGrid%avf2b(1,ifl)
      pBorder => pGrid%borders(iBorder)
      
      icl = pGrid%avf2b(2,ifl)
      
      IF ( icl <= pBorder%nCellsRecv ) THEN 
        icg = pBorder%icgRecv(icl)
      ELSE 
! TEMPORARY
        WRITE(*,*) 'ERROR! Exceeding dims of pBorder%icgRecv'
        STOP
! END TEMPORARY        
      END IF ! icl
      
      pGrid%avf2p(ifl) = CRAZY_VALUE_INT
      
      patchLoop: DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        IF ( pPatch%nBCellsVirt > 0 ) THEN 
          CALL BinarySearchInteger(pPatch%bvcSorted,pPatch%nBCellsVirt, & 
                                   icg,iLoc)

          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
            pGrid%avf2p(ifl) = iPatch

            EXIT patchLoop
          END IF ! iLoc
        END IF ! pPatch%nBCellsVirt                               
      END DO patchLoop
    END DO ! ifl
      
! ******************************************************************************
!   Destroy temporary lists 
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( pPatch%nBCellsVirt > 0 ) THEN 
        DEALLOCATE(pPatch%bvcSorted,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvcSorted')
        END IF ! global%error
      END IF ! pPatch%nBCellsVirt
    END DO ! iPatch
          
! ******************************************************************************
!   End 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building av-face-to-patch list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildAVFace2PatchList








! ******************************************************************************
!
! Purpose: Build cell-to-face list.
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
    
  SUBROUTINE RFLU_BuildCell2FaceList(pRegion)

    USE RFLU_ModCellFaceEdgeInfo

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
      
    INTEGER :: errorFlag,icl,icg,ick,ifg,ifl,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildCell2FaceList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-face list...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Initialize
! ******************************************************************************

    IF ( pGrid%nTetsTot > 0 ) THEN 
      DO icl = 1,pGrid%nTetsTot ! Explicit loops because of ASCI White
        DO ifl = 1,4
          pGrid%tet2f(1,ifl,icl) = C2F_INIT 
          pGrid%tet2f(2,ifl,icl) = C2F_INIT 
        END DO ! ifl
      END DO ! icl
    END IF ! pGrid%nTetsTot

    IF ( pGrid%nHexsTot > 0 ) THEN 
      DO icl = 1,pGrid%nHexsTot ! Explicit loops because of ASCI White
        DO ifl = 1,6
          pGrid%hex2f(1,ifl,icl) = C2F_INIT 
          pGrid%hex2f(2,ifl,icl) = C2F_INIT 
        END DO ! ifl 
      END DO ! icl
    END IF ! pGrid%nHexsTot

    IF ( pGrid%nPrisTot > 0 ) THEN 
      DO icl = 1,pGrid%nPrisTot ! Explicit loops because of ASCI White
        DO ifl = 1,5
          pGrid%pri2f(1,ifl,icl) = C2F_INIT 
          pGrid%pri2f(2,ifl,icl) = C2F_INIT 
        END DO ! ifl 
      END DO ! icl
    END IF ! pGrid%nPrisTot

    IF ( pGrid%nPyrsTot > 0 ) THEN 
      DO icl = 1,pGrid%nPyrsTot ! Explicit loops because of ASCI White
        DO ifl = 1,5
          pGrid%pyr2f(1,ifl,icl) = C2F_INIT 
          pGrid%pyr2f(2,ifl,icl) = C2F_INIT 
        END DO ! ifl 
      END DO ! icl                                           
    END IF ! pGrid%nPyrsTot          

! ******************************************************************************
!   Build cell-to-face list
! ******************************************************************************

! ==============================================================================
!   Interior faces
! ==============================================================================

    iPatch = 0 ! Indicates interior face

    DO ifg = 1,pGrid%nFacesTot
      DO icl = 1,2      
        icg = pGrid%f2c(icl,ifg)
        ick = RFLU_GetGlobalCellKind(global,pGrid,icg)        

        IF ( ick /= CELL_KIND_EXT .AND. ick /= CELL_KIND_BND ) THEN
          CALL RFLU_InsertIntoCell2FaceList(pRegion,iPatch,icg,ifg)
        END IF ! ick
      END DO ! icl     
    END DO ! ifg

! ==============================================================================
!   Boundary faces. NOTE skip periodic and symmetry patches because otherwise
!   get duplicated entries and hence error message from algorithm in routine
!   RFLU_InsertIntoCell2FaceList
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%bcType /= BC_PERIODIC .AND. & 
           pPatch%bcType /= BC_SYMMETRY ) THEN
        DO ifl = 1,pPatch%nBFacesTot            
          icg = pPatch%bf2c(ifl)
          ick = RFLU_GetGlobalCellKind(global,pGrid,icg)  

          CALL RFLU_InsertIntoCell2FaceList(pRegion,iPatch,icg,ifl)          
        END DO ! ifl
      END IF ! pPatch%bcType
    END DO ! iPatch

! ******************************************************************************
!   Check that all cells have all faces defined
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN 
      IF ( pGrid%nTetsTot > 0 ) THEN 
        DO icl = 1,pGrid%nTetsTot ! Explicit loops because of ASCI White
          DO ifl = 1,4
            IF ( pGrid%tet2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%tet2f(2,ifl,icl) == C2F_INIT ) THEN 
              CALL ErrorStop(global,ERR_C2FLIST_INVALID,__LINE__)   
            END IF ! pGrid%tet2f
          END DO ! ifl
        END DO ! icl
      END IF ! pGrid%nTetsTot

      IF ( pGrid%nHexsTot > 0 ) THEN       
        DO icl = 1,pGrid%nHexsTot ! Explicit loops because of ASCI White
          DO ifl = 1,6
            IF ( pGrid%hex2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%hex2f(2,ifl,icl) == C2F_INIT ) THEN 
              CALL ErrorStop(global,ERR_C2FLIST_INVALID,__LINE__)   
            END IF ! pGrid%hex2f
          END DO ! ifl
        END DO ! icl
      END IF ! pGrid%nHexsTot

      IF ( pGrid%nPrisTot > 0 ) THEN 
        DO icl = 1,pGrid%nPrisTot ! Explicit loops because of ASCI White
          DO ifl = 1,5
            IF ( pGrid%pri2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%pri2f(2,ifl,icl) == C2F_INIT ) THEN 
              CALL ErrorStop(global,ERR_C2FLIST_INVALID,__LINE__)   
            END IF ! pGrid%pri2f
          END DO ! ifl
        END DO ! icl
      END IF ! pGrid%nPrisTot

      IF ( pGrid%nPyrsTot > 0 ) THEN 
        DO icl = 1,pGrid%nPyrsTot ! Explicit loops because of ASCI White
          DO ifl = 1,5
            IF ( pGrid%pyr2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%pyr2f(2,ifl,icl) == C2F_INIT ) THEN 
              CALL ErrorStop(global,ERR_C2FLIST_INVALID,__LINE__)   
            END IF ! pGrid%pyr2f
          END DO ! ifl
        END DO ! icl
      END IF ! pGrid%nPyrsTot                                
    END IF ! global%checkLevel
        
#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Data structure output for checking
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Face-to-cell list'

    IF ( pGrid%nTetsTot > 0 ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Tetrahedra:'
      DO icl = 1,pGrid%nTetsTot
        DO ifl = 1,4
          WRITE(STDOUT,'(A,4(1X,I6))') SOLVER_NAME,icl,ifl, &
                                       pGrid%tet2f(1:2,ifl,icl)
        END DO ! ifc
      END DO ! icl
    END IF ! pGrid%nTetsTot

    IF ( pGrid%nHexsTot > 0 ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Hexahedra:'
      DO icl = 1,pGrid%nHexsTot
        DO ifl = 1,6
          WRITE(STDOUT,'(A,4(1X,I6))') SOLVER_NAME,icl,ifl, &
                                       pGrid%hex2f(1:2,ifl,icl)
        END DO ! ifl
      END DO ! icl
    END IF ! pGrid%nHexsTot      

    IF ( pGrid%nPrisTot > 0 ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Prisms:'
      DO icl = 1,pGrid%nPrisTot
        DO ifl = 1,5
          WRITE(STDOUT,'(A,4(1X,I6))') SOLVER_NAME,icl,ifl, &
                                       pGrid%pri2f(1:2,ifl,icl)
        END DO ! ifc
      END DO ! ic
    END IF ! pGrid%nPrisTot        

    IF ( pGrid%nPyrsTot > 0 ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Pyramids:'
      DO icl = 1,pGrid%nPyrsTot
        DO ifl = 1,5
          WRITE(STDOUT,'(A,4(1X,I6))') SOLVER_NAME,icl,ifl, &
                                       pGrid%pyr2f(1:2,ifl,icl)
        END DO ! ifl
      END DO ! icl
    END IF ! pGrid%nPyrsTot      

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME
#endif

! ******************************************************************************
!   End 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-face list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildCell2FaceList
    






  
! ******************************************************************************
!
! Purpose: Build face list.
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
    
  SUBROUTINE RFLU_BuildFaceList(pRegion)

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

    INTEGER :: c1,c1k,c1t,c2,c2k,c2t,errorFlag,fCntr,fType,faceType,fkSum, & 
               fvSize,icg,icl,ict,ifc,ifg,ifk,ifl,iPatch,iq,it,key, & 
               nBFacesQuad,nBFacesTri,nFacesInt,nFacesQuad,nFacesTri, & 
               nHexsLow,nHexsUpp,nPrisLow,nPrisUpp,nPyrsLow,nPyrsUpp, & 
               nTetsLow,nTetsUpp,v1g,v2g,v3g,v4g
    INTEGER :: fv(4) 
    INTEGER :: fkCntr(FACE_KIND_AA:FACE_KIND_AB,FACE_TYPE_TRI:FACE_TYPE_QUAD), & 
               fkOffs(FACE_KIND_AA:FACE_KIND_AB)
    INTEGER, DIMENSION(:), ALLOCATABLE :: fKind
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global

! TEMPORARY
    INTEGER :: nBFaces
! END TEMPORARY

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildFaceList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building face list...' 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal        
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building hash table...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer and initialize 
! ******************************************************************************

    pGrid => pRegion%grid
            
    pGrid%nFacesTot = 0
 
    nBFacesQuad = 0
    nBFacesTri  = 0
    nFacesTri   = 0
    nFacesQuad  = 0

! ******************************************************************************
!   Create hash table
! ******************************************************************************
 
    CALL RFLU_CreateHashTable(global,pGrid%nFacesEst)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN                              
      WRITE(STDOUT,'(A,5X,A,12X,I9)') SOLVER_NAME,'Hash table size:', &
                                      hashTableSize
    END IF ! global%verbLevel   

! ******************************************************************************
!   Loop over cell types, construct hash table of faces
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN                              
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Looping over cell types...'
    END IF ! global%verbLevel  

! ==============================================================================
!   Tetrahedra
! ==============================================================================

    IF ( pGrid%nTetsTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel
    END IF ! pGrid%nTetsTot

    fvSize = 3 ! Tetrahedra only have 3 vertices per face

    DO icl = 1,pGrid%nTetsTot      
      icg = pGrid%tet2CellGlob(icl)      

      DO ifl = 1,4   
        fv(1) = pGrid%tet2v(f2vTet(1,ifl),icl)
        fv(2) = pGrid%tet2v(f2vTet(2,ifl),icl)
        fv(3) = pGrid%tet2v(f2vTet(3,ifl),icl)          

        CALL QuickSortInteger(fv(1:fvSize),fvSize)         
        CALL RFLU_HashBuildKey(fv(1:3),3,key)             
        CALL RFLU_HashFace(global,key,pGrid,icg,ifl,fv(1:3),fvSize,faceType)  

        IF ( faceType == FACE_TYPE_NEW ) THEN 
          nFacesTri = nFacesTri + 1
        END IF ! faceType          
      END DO ! ifl
    END DO ! icl

! ==============================================================================
!   Hexahedra
! ==============================================================================

    IF ( pGrid%nHexsTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Hexahedra...'
      END IF ! global%verbLevel
    END IF ! pGrid%nHexsTot

    fvSize = 4 ! Hexahedra only have 4 vertices per face

    DO icl = 1,pGrid%nHexsTot
      icg = pGrid%hex2CellGlob(icl)  

      DO ifl = 1,6
        fv(1) = pGrid%hex2v(f2vHex(1,ifl),icl)
        fv(2) = pGrid%hex2v(f2vHex(2,ifl),icl)
        fv(3) = pGrid%hex2v(f2vHex(3,ifl),icl)          
        fv(4) = pGrid%hex2v(f2vHex(4,ifl),icl)     

        CALL QuickSortInteger(fv(1:fvSize),fvSize)
        CALL RFLU_HashBuildKey(fv(1:3),3,key)                    
        CALL RFLU_HashFace(global,key,pGrid,icg,ifl,fv(1:3),fvSize,faceType)

        IF ( faceType == FACE_TYPE_NEW ) THEN 
          nFacesQuad = nFacesQuad + 1
        END IF ! faceType                
      END DO ! ifl
    END DO ! icl

! ==============================================================================
!   Prisms
! ==============================================================================

    IF ( pGrid%nPrisTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Prisms...'
      END IF ! global%verbLevel
    END IF ! pGrid%nPrisTot

    DO icl = 1,pGrid%nPrisTot
      icg = pGrid%pri2CellGlob(icl)

      DO ifl = 1,5
        fvSize = 3
        fv(1)  = pGrid%pri2v(f2vPri(1,ifl),icl) 
        fv(2)  = pGrid%pri2v(f2vPri(2,ifl),icl) 
        fv(3)  = pGrid%pri2v(f2vPri(3,ifl),icl)       

        IF ( f2vPri(4,ifl) /= VERT_NONE ) THEN
          fvSize = 4
          fv(4)  = pGrid%pri2v(f2vPri(4,ifl),icl)
        END IF ! f2v%pri      

        CALL QuickSortInteger(fv(1:fvSize),fvSize)
        CALL RFLU_HashBuildKey(fv(1:3),3,key)    
        CALL RFLU_HashFace(global,key,pGrid,icg,ifl,fv(1:3),fvSize,faceType)

        IF ( faceType == FACE_TYPE_NEW ) THEN
          IF ( fvSize == 4 ) THEN   
            nFacesQuad = nFacesQuad + 1
          ELSE 
            nFacesTri  = nFacesTri  + 1
          END IF ! fvSize
        END IF ! faceType
      END DO ! ifl
    END DO ! icl  

! ==============================================================================
!   Pyramids
! ==============================================================================

    IF ( pGrid%nPyrsTot /= 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN   
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Pyramids...'
      END IF ! global%verbLevel
    END IF ! pGrid%nPyrsTot

    DO icl = 1,pGrid%nPyrsTot
      icg = pGrid%pyr2CellGlob(icl)

      DO ifl = 1,5
        fvSize = 3
        fv(1)  = pGrid%pyr2v(f2vPyr(1,ifl),icl) 
        fv(2)  = pGrid%pyr2v(f2vPyr(2,ifl),icl) 
        fv(3)  = pGrid%pyr2v(f2vPyr(3,ifl),icl)    

        IF ( f2vPyr(4,ifl) /= VERT_NONE ) THEN
          fvSize = 4
          fv(4)  = pGrid%pyr2v(f2vPyr(4,ifl),icl)
        END IF ! f2v%pyr      

        CALL QuickSortInteger(fv(1:fvSize),fvSize)
        CALL RFLU_HashBuildKey(fv(1:3),3,key)    
        CALL RFLU_HashFace(global,key,pGrid,icg,ifl,fv(1:3),fvSize,faceType)

        IF ( faceType == FACE_TYPE_NEW ) THEN
          IF ( fvSize == 3 ) THEN   
            nFacesTri  = nFacesTri  + 1
          ELSE 
            nFacesQuad = nFacesQuad + 1
          END IF ! fvSize
        END IF ! faceType 
      END DO ! ifl
    END DO ! icl    

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,5X,A,5X,I9)') SOLVER_NAME,'Hash table collisions: ', & 
                                     hashTableCollisions  
    END IF ! global%myProcid

! ******************************************************************************
!   Number of internal faces (not on boundaries). NOTE nFacesTot is computed
!   in RFLU_HashFace. At this point it is the total number of faces. Further
!   below, it will be redefined to include only those faces which are not on
!   boundaries.
! ******************************************************************************

    nFacesInt = pGrid%nFacesTot - pGrid%nBFacesTot    
                                         
! ******************************************************************************
!   Build actual facelist from hash table
! ****************************************************************************** 

! ==============================================================================
!   Actual facelist - boundary faces. NOTE combined face list is ordered such
!   that actual triangular faces appear first, then quadrilateral actual
!   faces, virtual triangular faces, and quadrilateral virtual faces last.
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building boundary face lists...'
    END IF ! global%verbLevel  

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      nBFacesTri  = nBFacesTri  + pPatch%nBTrisTot
      nBFacesQuad = nBFacesQuad + pPatch%nBQuadsTot      

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,5X,A,I4,1X,A,1X,I4,A)') SOLVER_NAME,'Patch: ',iPatch, & 
                                    '(Global patch:',pPatch%iPatchGlobal,')'

        WRITE(STDOUT,'(A,7X,A,8X,2(1X,I9))') SOLVER_NAME, & 
                                             'Triangular faces:', & 
                                             pPatch%nBTris, &
                                             pPatch%nBTrisTot-pPatch%nBTris 
        WRITE(STDOUT,'(A,7X,A,5X,2(1X,I9))') SOLVER_NAME, & 
                                             'Quadrilateral faces:', & 
                                             pPatch%nBQuads, &
                                             pPatch%nBQuadsTot-pPatch%nBQuads 
        WRITE(STDOUT,'(A,7X,A,3X,2(1X,I9))') SOLVER_NAME, & 
                                             'Total number of faces:', &
                                             pPatch%nBFaces, &
                                             pPatch%nBFacesTot-pPatch%nBFaces
      END IF ! global%verbLevel     

! ------------------------------------------------------------------------------
!     Triangular faces - NOTE sorting of actual and virtual faces...
! ------------------------------------------------------------------------------

      fvSize = 3

      DO it = 1,pPatch%nBTrisTot
        fv(1) = pPatch%bTri2v(1,it)
        fv(2) = pPatch%bTri2v(2,it)
        fv(3) = pPatch%bTri2v(3,it)

        CALL QuickSortInteger(fv(1:fvSize),fvSize)
        CALL RFLU_HashBuildKey(fv(1:3),3,key)              
        CALL RFLU_UnHashBFace(global,key,pGrid,fv(1:3),fvSize,pPatch%bcType, & 
                              icg,ifg)

        fCntr = it

        IF ( it > pPatch%nBTris ) THEN 
          fCntr = it + pPatch%nBQuads
        END IF ! it

        pPatch%bf2c(  fCntr) = icg
        pPatch%bf2v(1,fCntr) = pPatch%bTri2v(1,it)
        pPatch%bf2v(2,fCntr) = pPatch%bTri2v(2,it)
        pPatch%bf2v(3,fCntr) = pPatch%bTri2v(3,it)                            
      END DO ! it

! ------------------------------------------------------------------------------
!     Quadrilateral faces - NOTE sorting of actual and virtual faces...
! ------------------------------------------------------------------------------

      fvSize = 4   

      DO iq = 1,pPatch%nBQuadsTot
        fv(1) = pPatch%bQuad2v(1,iq)
        fv(2) = pPatch%bQuad2v(2,iq)
        fv(3) = pPatch%bQuad2v(3,iq)
        fv(4) = pPatch%bQuad2v(4,iq)                 

        CALL QuickSortInteger(fv(1:fvSize),fvSize)
        CALL RFLU_HashBuildKey(fv(1:3),3,key)
        CALL RFLU_UnHashBFace(global,key,pGrid,fv(1:3),fvSize,pPatch%bcType, & 
                              icg,ifg)

        IF ( iq <= pPatch%nBQuads ) THEN 
          fCntr = iq + pPatch%nBTris
        ELSE 
          fCntr = iq + pPatch%nBTrisTot
        END IF ! iq

        pPatch%bf2c(  fCntr) = icg
        pPatch%bf2v(1,fCntr) = pPatch%bQuad2v(1,iq)
        pPatch%bf2v(2,fCntr) = pPatch%bQuad2v(2,iq)
        pPatch%bf2v(3,fCntr) = pPatch%bQuad2v(3,iq)
        pPatch%bf2v(4,fCntr) = pPatch%bQuad2v(4,iq)
      END DO ! iq
    END DO ! iPatch

! ==============================================================================
!   Destroy hash table
! ==============================================================================

    CALL RFLU_DestroyHashTable(global)

! ==============================================================================
!   Deallocate f2v array, not needed, regenerated below automatically
! ==============================================================================

    DEALLOCATE(pGrid%f2v,STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2v')
    END IF ! global%error 
           
! ******************************************************************************
!   Determine face kinds. NOTE do not modify fkCntr array, it is used below. 
!   NOTE also: These statistics must be computed AFTER the boundary face lists
!   are built, because before cannot distinguish between faces of kinds 
!   FACE_KIND_VB and FACE_KIND_VX. This is because f2c array is initialized 
!   with CELL_TYPE_EXT, and boundary-face lists, which are given explicitly
!   in grid file, have not yet been used to mark those faces. 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining face statistics...'
    END IF ! global%verbLevel  
    
    fkCntr(FACE_KIND_AA,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0 
    fkCntr(FACE_KIND_AV,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0 
    fkCntr(FACE_KIND_VV,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0 
    fkCntr(FACE_KIND_VB,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0   
    fkCntr(FACE_KIND_VX,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0 
    fkCntr(FACE_KIND_AB,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0

    DO ifc = 1,pGrid%nFacesTot
      c1 = pGrid%f2c(1,ifc)
      c2 = pGrid%f2c(2,ifc)

      ifl = pGrid%f2c(3,ifc)
      c1t = RFLU_GetGlobalCellType(global,pGrid,c1)

      SELECT CASE ( c1t ) 
        CASE ( CELL_TYPE_TET ) ! Tetrahedral cell 
          fType = FACE_TYPE_TRI           
        CASE ( CELL_TYPE_HEX ) ! Hexahedral cell
          fType = FACE_TYPE_QUAD            
        CASE ( CELL_TYPE_PRI ) ! Prismatic cell
          IF ( f2vPri(4,ifl) /= VERT_NONE ) THEN 
            fType = FACE_TYPE_QUAD
          ELSE 
            fType = FACE_TYPE_TRI 
          END IF ! f2vPri            
        CASE ( CELL_TYPE_PYR ) ! Pyramidal cell
          IF ( f2vPyr(4,ifl) == VERT_NONE ) THEN 
            fType = FACE_TYPE_TRI 
          ELSE 
            fType = FACE_TYPE_QUAD 
          END IF ! f2vPyr            
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! cellType

      c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
      c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)                

      SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) )       
        CASE ( FACE_KIND_AA )
          fkCntr(FACE_KIND_AA,fType) = fkCntr(FACE_KIND_AA,fType) + 1     
        CASE ( FACE_KIND_AV )
          fkCntr(FACE_KIND_AV,fType) = fkCntr(FACE_KIND_AV,fType) + 1
        CASE ( FACE_KIND_VV )
          fkCntr(FACE_KIND_VV,fType) = fkCntr(FACE_KIND_VV,fType) + 1
        CASE ( FACE_KIND_VB )
          fkCntr(FACE_KIND_VB,fType) = fkCntr(FACE_KIND_VB,fType) + 1
        CASE ( FACE_KIND_AB )
          fkCntr(FACE_KIND_AB,fType) = fkCntr(FACE_KIND_AB,fType) + 1
        CASE ( FACE_KIND_VX )
          fkCntr(FACE_KIND_VX,fType) = fkCntr(FACE_KIND_VX,fType) + 1            
        CASE DEFAULT        
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! fKind
    END DO ! ifc
            
! ==============================================================================
!   Write info 
! ==============================================================================
            
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,5X,A)')       SOLVER_NAME,'Face-type statistics:' 
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Total faces:           ', & 
                                     pGrid%nFacesTot
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     nFacesTri
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     nFacesQuad
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Total boundary faces:  ', & 
                                     pGrid%nBFacesTot 
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     nBFacesTri
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     nBFacesQuad         
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Non-boundary faces:    ', & 
                                     nFacesInt 
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     nFacesTri - nBFacesTri
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     nFacesQuad - nBFacesQuad
      WRITE(STDOUT,'(A,5X,A)')       SOLVER_NAME,'Face-kind statistics:' 
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Actual-actual faces:   ', & 
                                     fkCntr(FACE_KIND_AA,FACE_TYPE_TRI ) +  & 
                                     fkCntr(FACE_KIND_AA,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     fkCntr(FACE_KIND_AA,FACE_TYPE_TRI)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     fkCntr(FACE_KIND_AA,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Actual-virtual faces:  ', & 
                                     fkCntr(FACE_KIND_AV,FACE_TYPE_TRI ) +  & 
                                     fkCntr(FACE_KIND_AV,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     fkCntr(FACE_KIND_AV,FACE_TYPE_TRI)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     fkCntr(FACE_KIND_AV,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Virtual-virtual faces: ', & 
                                     fkCntr(FACE_KIND_VV,FACE_TYPE_TRI ) +  & 
                                     fkCntr(FACE_KIND_VV,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     fkCntr(FACE_KIND_VV,FACE_TYPE_TRI)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     fkCntr(FACE_KIND_VV,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Virtual-boundary faces:', & 
                                     fkCntr(FACE_KIND_VB,FACE_TYPE_TRI ) +  & 
                                     fkCntr(FACE_KIND_VB,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     fkCntr(FACE_KIND_VB,FACE_TYPE_TRI)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     fkCntr(FACE_KIND_VB,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Actual-boundary faces: ', & 
                                     fkCntr(FACE_KIND_AB,FACE_TYPE_TRI ) +  & 
                                     fkCntr(FACE_KIND_AB,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     fkCntr(FACE_KIND_AB,FACE_TYPE_TRI)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     fkCntr(FACE_KIND_AB,FACE_TYPE_QUAD) 
      WRITE(STDOUT,'(A,7X,A,5X,I9)') SOLVER_NAME,'Virtual-external faces:', & 
                                     fkCntr(FACE_KIND_VX,FACE_TYPE_TRI ) +  & 
                                     fkCntr(FACE_KIND_VX,FACE_TYPE_QUAD)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Triangular faces:      ', & 
                                     fkCntr(FACE_KIND_VX,FACE_TYPE_TRI)
      WRITE(STDOUT,'(A,9X,A,3X,I9)') SOLVER_NAME,'Quadrilateral faces:   ', & 
                                     fkCntr(FACE_KIND_VX,FACE_TYPE_QUAD)
    END IF ! global%myProcid             

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining face statistics done.'
    END IF ! global%verbLevel 
          
! ******************************************************************************
!   Actual facelist - interior faces
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building non-boundary face '// & 
                               'lists...'
    END IF ! global%verbLevel

! ==============================================================================
!   Determine face kind offsets for sorting of faces 
! ==============================================================================

    fkOffs(FACE_KIND_AA) = 0
    fkOffs(FACE_KIND_AV) = fkCntr(FACE_KIND_AA,FACE_TYPE_TRI ) & 
                         + fkCntr(FACE_KIND_AA,FACE_TYPE_QUAD)
    fkOffs(FACE_KIND_VV) = fkOffs(FACE_KIND_AV) &
                         + fkCntr(FACE_KIND_AV,FACE_TYPE_TRI ) & 
                         + fkCntr(FACE_KIND_AV,FACE_TYPE_QUAD)
    fkOffs(FACE_KIND_VX) = fkOffs(FACE_KIND_VV) & 
                         + fkCntr(FACE_KIND_VV,FACE_TYPE_TRI ) & 
                         + fkCntr(FACE_KIND_VV,FACE_TYPE_QUAD)
               
! ==============================================================================
!   Determine number of flux faces (nFaces) and total internal faces. NOTE need
!   to set nBFaces and nBFacesTot here (in addition to RFLU_CreateGrid) because 
!   when adding virtual cells for symmetry or periodic patches, face-list is 
!   regenerated because of new cells, but RFLU_CreateGrid is not called again.
! ==============================================================================
               
    pGrid%nFaces = fkCntr(FACE_KIND_AA,FACE_TYPE_TRI ) & 
                 + fkCntr(FACE_KIND_AA,FACE_TYPE_QUAD) &       
                 + fkCntr(FACE_KIND_AV,FACE_TYPE_TRI ) &                 
                 + fkCntr(FACE_KIND_AV,FACE_TYPE_QUAD)   
                           
    pGrid%nFacesAV = fkCntr(FACE_KIND_AV,FACE_TYPE_TRI ) &                 
                   + fkCntr(FACE_KIND_AV,FACE_TYPE_QUAD) 
    pGrid%nFacesVV = fkCntr(FACE_KIND_VV,FACE_TYPE_TRI ) &                 
                   + fkCntr(FACE_KIND_VV,FACE_TYPE_QUAD) 

    nFacesInt = fkCntr(FACE_KIND_AA,FACE_TYPE_TRI ) & 
              + fkCntr(FACE_KIND_AA,FACE_TYPE_QUAD) &       
              + fkCntr(FACE_KIND_AV,FACE_TYPE_TRI ) &                 
              + fkCntr(FACE_KIND_AV,FACE_TYPE_QUAD) & 
              + fkCntr(FACE_KIND_VV,FACE_TYPE_TRI ) &                 
              + fkCntr(FACE_KIND_VV,FACE_TYPE_QUAD) & 
              + fkCntr(FACE_KIND_VX,FACE_TYPE_TRI ) &                 
              + fkCntr(FACE_KIND_VX,FACE_TYPE_QUAD)                 

    pGrid%nBFaces = fkCntr(FACE_KIND_AB,FACE_TYPE_TRI ) & 
                  + fkCntr(FACE_KIND_AB,FACE_TYPE_QUAD)
                  
    pGrid%nBFacesTot = fkCntr(FACE_KIND_AB,FACE_TYPE_TRI ) & 
                     + fkCntr(FACE_KIND_AB,FACE_TYPE_QUAD) &
                     + fkCntr(FACE_KIND_VB,FACE_TYPE_TRI ) & 
                     + fkCntr(FACE_KIND_VB,FACE_TYPE_QUAD)                  
                 
! ==============================================================================
!   Allocate temporary face arrays and initialize arrays and counters. NOTE 
!   in the following, write counter information into fkCntr(ifk,FACE_TYPE_TRI)
!   for convenience even if face is not a triangle
! ==============================================================================

    ALLOCATE(pGrid%f2cTemp(2,nFacesInt),STAT=errorFlag) ! NOTE nFacesInt
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cTemp')
    END IF ! global%error

    DO ifc = 1,nFacesInt ! Explicit loop because of Frost
      pGrid%f2cTemp(1,ifc) = 0 ! Initial value immaterial
      pGrid%f2cTemp(2,ifc) = 0 ! Initial value immaterial
    END DO ! ifc

    ALLOCATE(pGrid%f2vTemp(4,nFacesInt),STAT=errorFlag) ! NOTE nFacesInt 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2vTemp')
    END IF ! global%error       

    DO ifc = 1,nFacesInt ! Explicit loop because of Frost
      pGrid%f2vTemp(1,ifc) = 0 ! Initial value immaterial            
      pGrid%f2vTemp(2,ifc) = 0 ! Initial value immaterial    
      pGrid%f2vTemp(3,ifc) = 0 ! Initial value immaterial    
      pGrid%f2vTemp(4,ifc) = 0 ! Initial value immaterial    
    END DO ! ifc

    fCntr = 0  

    fkCntr(FACE_KIND_AA,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0       
    fkCntr(FACE_KIND_AV,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0  
    fkCntr(FACE_KIND_VV,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0
    fkCntr(FACE_KIND_VX,FACE_TYPE_TRI:FACE_TYPE_QUAD) = 0                    

! ==============================================================================
!   Loop over all faces
! ==============================================================================

    DO ifg = 1,pGrid%nFacesTot
      
! ------------------------------------------------------------------------------
!     Determine face kind
! ------------------------------------------------------------------------------
      
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

      c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
      c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)         
      ifk = RFLU_GetFaceKind(global,c1k,c2k)
      
! ------------------------------------------------------------------------------
!     If not on boundary, add to face list 
! ------------------------------------------------------------------------------
      
      IF ( (ifk /= FACE_KIND_AB) .AND. (ifk /= FACE_KIND_VB) ) THEN 
        ict = pGrid%cellGlob2Loc(1,c1) ! cell type
        icl = pGrid%cellGlob2Loc(2,c1) ! local cell index
        ifl = pGrid%f2c(3,ifg)         ! local face index

        SELECT CASE ( ict ) 
          CASE ( CELL_TYPE_TET ) ! Tetrahedral cell 
            v1g = pGrid%tet2v(f2vTet(1,ifl),icl)
            v2g = pGrid%tet2v(f2vTet(2,ifl),icl)
            v3g = pGrid%tet2v(f2vTet(3,ifl),icl)
            v4g = VERT_NONE          
          CASE ( CELL_TYPE_HEX ) ! Hexahedral cell
            v1g = pGrid%hex2v(f2vHex(1,ifl),icl)
            v2g = pGrid%hex2v(f2vHex(2,ifl),icl)
            v3g = pGrid%hex2v(f2vHex(3,ifl),icl)
            v4g = pGrid%hex2v(f2vHex(4,ifl),icl)              
          CASE ( CELL_TYPE_PRI ) ! Prismatic cell
            v1g = pGrid%pri2v(f2vPri(1,ifl),icl)
            v2g = pGrid%pri2v(f2vPri(2,ifl),icl)
            v3g = pGrid%pri2v(f2vPri(3,ifl),icl)
            v4g = VERT_NONE

            IF ( f2vPri(4,ifl) /= VERT_NONE ) THEN 
              v4g = pGrid%pri2v(f2vPri(4,ifl),icl)  
            END IF ! f2vPri            
          CASE ( CELL_TYPE_PYR ) ! Pyramidal cell
            v1g = pGrid%pyr2v(f2vPyr(1,ifl),icl)
            v2g = pGrid%pyr2v(f2vPyr(2,ifl),icl)
            v3g = pGrid%pyr2v(f2vPyr(3,ifl),icl)
            v4g = VERT_NONE

            IF ( f2vPyr(4,ifl) /= VERT_NONE ) THEN 
              v4g = pGrid%pyr2v(f2vPyr(4,ifl),icl)    
            END IF ! f2vPyr            
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! cellType

! ------------------------------------------------------------------------------
!       Insert into face list 
! ------------------------------------------------------------------------------

        fkCntr(ifk,FACE_TYPE_TRI) = fkCntr(ifk,FACE_TYPE_TRI) + 1
        fCntr = fkCntr(ifk,FACE_TYPE_TRI) + fkOffs(ifk)       

        pGrid%f2cTemp(1,fCntr) = pGrid%f2c(1,ifg)
        pGrid%f2cTemp(2,fCntr) = pGrid%f2c(2,ifg)      

        pGrid%f2vTemp(1,fCntr) = v1g
        pGrid%f2vTemp(2,fCntr) = v2g
        pGrid%f2vTemp(3,fCntr) = v3g
        pGrid%f2vTemp(4,fCntr) = v4g
      END IF ! ifk
    END DO ! ifg

! ------------------------------------------------------------------------------
!   Check for consistency 
! ------------------------------------------------------------------------------

    fkSum = fkCntr(FACE_KIND_AA,FACE_TYPE_TRI) & 
          + fkCntr(FACE_KIND_AV,FACE_TYPE_TRI) &
          + fkCntr(FACE_KIND_VV,FACE_TYPE_TRI) & 
          + fkCntr(FACE_KIND_VX,FACE_TYPE_TRI)

    IF ( fkSum /= nFacesInt ) THEN ! Incorrect number of internal faces
      CALL ErrorStop(global,ERR_NFACES_WRONG,__LINE__) 
    END IF ! fkSum

    pGrid%nFacesTot = nFacesInt

! ==============================================================================
!   Deallocate arrays
! ==============================================================================

    DEALLOCATE(pGrid%f2c,STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2c')
    END IF ! global%error     

! ==============================================================================
!   Reallocate arrays, copy from temporary arrays, and deallocate those
! ==============================================================================

    ALLOCATE(pGrid%f2c(2,pGrid%nFacesTot),STAT=errorFlag) ! NOTE nFacesInt
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2c')
    END IF ! global%error

    DO ifc = 1,pGrid%nFacesTot ! Explicit loop because of Frost
      pGrid%f2c(1,ifc) = pGrid%f2cTemp(1,ifc)
      pGrid%f2c(2,ifc) = pGrid%f2cTemp(2,ifc)
    END DO ! ifc

    DEALLOCATE(pGrid%f2cTemp,STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2cTemp')
    END IF ! global%error     

    ALLOCATE(pGrid%f2v(4,pGrid%nFacesTot),STAT=errorFlag) ! NOTE nFacesInt 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2v')
    END IF ! global%error       

    DO ifc = 1,pGrid%nFacesTot ! Explicit loop because of Frost
      pGrid%f2v(1,ifc) = pGrid%f2vTemp(1,ifc)   
      pGrid%f2v(2,ifc) = pGrid%f2vTemp(2,ifc)   
      pGrid%f2v(3,ifc) = pGrid%f2vTemp(3,ifc)   
      pGrid%f2v(4,ifc) = pGrid%f2vTemp(4,ifc)   
    END DO ! ifc

    DEALLOCATE(pGrid%f2vTemp,STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2vTemp')
    END IF ! global%error             

#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Data structure output for checking
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Face neighbours and vertices:' 
                                      
    DO ifc = 1,pGrid%nFacesTot
      WRITE(STDOUT,'(A,7(1X,I6))') SOLVER_NAME,ifc,pGrid%f2c(1:2,ifc), & 
                                   pGrid%f2v(1:4,ifc)
    END DO ! ifc
    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME

    IF ( pGrid%nPatches > 0 ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Patch:',iPatch
        
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Face neighbour and vertices:'
      
        DO ifc = 1,pPatch%nBFacesTot
          WRITE(STDOUT,'(A,6(1X,I6))') SOLVER_NAME,ifc,pPatch%bf2c(ifc), & 
                                       pPatch%bf2v(1:4,ifc)
        END DO ! ifc
      END DO ! iPatch

      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME               
    END IF ! pGrid%nPatches      
#endif
            
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building face list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildFaceList
    
    
   
   
   
   
   
   
! ******************************************************************************
!
! Purpose: Create actual-virtual-face-to-border list.
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
  
  SUBROUTINE RFLU_CreateAVFace2BorderList(pRegion)

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

    CALL RegisterFunction(global,'RFLU_CreateAVFace2BorderList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating av-face-to-border list...'   
    END IF ! global%verbLevel  
     
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyAVFace2BorderList(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
          
    ALLOCATE(pGrid%avf2b(2,pGrid%nFacesAV),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%avf2b')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating av-face-to-border list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_CreateAVFace2BorderList
   
   
   
   
   
 
 
! ******************************************************************************
!
! Purpose: Create actual-virtual-face-to-patch list.
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
  
  SUBROUTINE RFLU_CreateAVFace2PatchList(pRegion)

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

    CALL RegisterFunction(global,'RFLU_CreateAVFace2PatchList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating av-face-to-patch list...'   
    END IF ! global%verbLevel  
     
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyAVFace2PatchList(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
          
    ALLOCATE(pGrid%avf2p(pGrid%nFacesAV),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%avf2p')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating av-face-to-patch list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_CreateAVFace2PatchList   
   
   
   
    
    
    
    
! ******************************************************************************
!
! Purpose: Create cell-to-face lists.
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
  
  SUBROUTINE RFLU_CreateCell2FaceList(pRegion)

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

    CALL RegisterFunction(global,'RFLU_CreateCell2FaceList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating cell-to-face list...'   
    END IF ! global%verbLevel  
     
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyCell2FaceList(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
          
    IF ( pGrid%nTetsTot > 0 ) THEN 
      ALLOCATE(pGrid%tet2f(2,4,pGrid%nTetsTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2f')
      END IF ! global%error
    END IF ! pGrid%nTetsTot    

    IF ( pGrid%nHexsTot > 0 ) THEN 
      ALLOCATE(pGrid%hex2f(2,6,pGrid%nHexsTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2f')
      END IF ! global%error
    END IF ! pGrid%nHexsTot          

    IF ( pGrid%nPrisTot > 0 ) THEN 
      ALLOCATE(pGrid%pri2f(2,5,pGrid%nPrisTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pri2f')
      END IF ! global%error
    END IF ! pGrid%nPrisTot          

    IF ( pGrid%nPyrsTot > 0 ) THEN 
      ALLOCATE(pGrid%pyr2f(2,5,pGrid%nPyrsTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pyr2f')
      END IF ! global%error
    END IF ! pGrid%nPyrsTot          
          
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating cell-to-face list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_CreateCell2FaceList










! ******************************************************************************
!
! Purpose: Create face list.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. The bf2c and bf2v arrays are not created here, because it is more 
!      convenient to do so when creating the grid.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_CreateFaceList(pRegion)    

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
  
    INTEGER :: errorFlag,ifc,iPatch,nBFaces
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch  

! ******************************************************************************
!   Start
! ******************************************************************************
    
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateFaceList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating face list...'   
    END IF ! global%verbLevel  
     
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyFaceList(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Estimate number of faces for allocation of hash table. The estimation      
!   formula shown below assumes that boundary effects are negligible, so for   
!   very small grids, where boundary faces dominate the total number of faces  , 
!   need a kludge. Kludge is also needed if running in parallel and boundary   
!   faces become a significant fraction of total number of faces, and when     
!   running code with periodic hack, where the missing boundaries mean that    
!   number of faces is underestimated.                                         
! ******************************************************************************
     
    nBFaces = 0

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      nBFaces = nBFaces + pPatch%nBTrisTot + pPatch%nBQuadsTot     
    END DO ! iPatch     
          
! ==============================================================================
!   Estimate total number of faces
! ==============================================================================

    pGrid%nFacesEst = nBFaces + 2*pGrid%nTetsTot + 3*pGrid%nHexsTot & 
                    + 5*pGrid%nPrisTot/2

    IF ( nBFaces/REAL(pGrid%nFacesEst,KIND=RFREAL) > 0.8_RFREAL .OR. & 
         nBFaces/REAL(pGrid%nFacesEst,KIND=RFREAL) < 0.3_RFREAL ) THEN 
      pGrid%nFacesEst = 2*pGrid%nFacesEst

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Corrected estimate', & 
                                      'of number of faces.' 
      END IF ! global%verbLevel       
    END IF ! nBFaces

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A,3X,I9)') SOLVER_NAME,'Estimated number of '// & 
                                     'faces: ',pGrid%nFacesEst       
    END IF ! global%verbLevel      
             
! ******************************************************************************
!   Allocate memory for interior faces
! ******************************************************************************
             
    ALLOCATE(pGrid%f2c(4,pGrid%nFacesEst),STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2c')
    END IF ! global%error

    DO ifc = 1,pGrid%nFacesEst 
      pGrid%f2c(1,ifc) = CELL_TYPE_EXT ! Initial value NOT immaterial 
      pGrid%f2c(2,ifc) = CELL_TYPE_EXT ! Initial value NOT immaterial 
      pGrid%f2c(3,ifc) = 0             ! Initial value     immaterial
      pGrid%f2c(4,ifc) = 0             ! Initial value     immaterial
    END DO ! ifc

    ALLOCATE(pGrid%f2v(3,pGrid%nFacesEst),STAT=errorFlag) 
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2v')
    END IF ! global%error  

    DO ifc = 1,pGrid%nFacesEst
      pGrid%f2v(1,ifc) = VERT_NONE ! Initial value NOT immaterial 
      pGrid%f2v(2,ifc) = VERT_NONE ! Initial value NOT immaterial
      pGrid%f2v(3,ifc) = VERT_NONE ! Initial value NOT immaterial
    END DO ! ifc

! ******************************************************************************
!   Allocate memory for boundary faces 
! ******************************************************************************
             
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      ALLOCATE(pPatch%bf2c(pPatch%nBFacesTot),STAT=errorFlag) 
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2c')
      END IF ! global%error                    

      ALLOCATE(pPatch%bf2v(4,pPatch%nBFacesTot),STAT=errorFlag) 
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2v')
      END IF ! global%error                    

      DO ifc = 1,pPatch%nBFacesTot                 
        pPatch%bf2c(ifc)   = 0         ! Initial value immaterial, not used
        pPatch%bf2v(1,ifc) = VERT_NONE ! Initial value NOT immaterial 
        pPatch%bf2v(2,ifc) = VERT_NONE ! Initial value NOT immaterial 
        pPatch%bf2v(3,ifc) = VERT_NONE ! Initial value NOT immaterial 
        pPatch%bf2v(4,ifc) = VERT_NONE ! Initial value NOT immaterial
      END DO ! ifc
    END DO ! iPatch             
                   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating face list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_CreateFaceList
  
  
   



! ******************************************************************************
!
! Purpose: Destroy actual-virtual-face-to-border list.
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

  SUBROUTINE RFLU_DestroyAVFace2BorderList(pRegion)
  
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

    CALL RegisterFunction(global,'RFLU_DestroyAVFace2BorderList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying av-face-to-border list...'   
    END IF ! global%verbLevel  
     
    pGrid => pRegion%grid
          
! ******************************************************************************
!   Deallocate memory
! ****************************************************************************** 
                   
    DEALLOCATE(pGrid%avf2b,STAT=errorFlag)    
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%avf2b')
    END IF ! global%error         

! ******************************************************************************
!   Nullify memory
! ******************************************************************************
          
    CALL RFLU_NullifyAVFace2BorderList(pRegion)
          
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying av-face-to-border list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_DestroyAVFace2BorderList
  
  





! ******************************************************************************
!
! Purpose: Destroy actual-virtual-face-to-patch list.
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

  SUBROUTINE RFLU_DestroyAVFace2PatchList(pRegion)
  
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

    CALL RegisterFunction(global,'RFLU_DestroyAVFace2PatchList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying av-face-to-patch list...'   
    END IF ! global%verbLevel  
     
    pGrid => pRegion%grid
          
! ******************************************************************************
!   Deallocate memory
! ****************************************************************************** 
                   
    DEALLOCATE(pGrid%avf2p,STAT=errorFlag)    
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%avf2p')
    END IF ! global%error         

! ******************************************************************************
!   Nullify memory
! ******************************************************************************
          
    CALL RFLU_NullifyAVFace2PatchList(pRegion)
          
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying av-face-to-patch list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_DestroyAVFace2PatchList





! ******************************************************************************
!
! Purpose: Destroy cell-to-face lists.
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

  SUBROUTINE RFLU_DestroyCell2FaceList(pRegion)
  
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

    CALL RegisterFunction(global,'RFLU_DestroyCell2FaceList',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying cell-to-face list...'   
    END IF ! global%verbLevel  
     
    pGrid => pRegion%grid
          
! ******************************************************************************
!   Deallocate memory
! ****************************************************************************** 
          
    IF ( pGrid%nTetsTot > 0 ) THEN 
      DEALLOCATE(pGrid%tet2f,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%tet2f')
      END IF ! global%error
    END IF ! pGrid%nTetsTot    

    IF ( pGrid%nHexsTot > 0 ) THEN 
      DEALLOCATE(pGrid%hex2f,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%hex2f')
      END IF ! global%error
    END IF ! pGrid%nHexsTot          

    IF ( pGrid%nPrisTot > 0 ) THEN 
      DEALLOCATE(pGrid%pri2f,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pri2f')
      END IF ! global%error
    END IF ! pGrid%nPrisTot          

    IF ( pGrid%nPyrsTot > 0 ) THEN 
      DEALLOCATE(pGrid%pyr2f,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pyr2f')
      END IF ! global%error
    END IF ! pGrid%nPyrsTot          

! ******************************************************************************
!   Nullify memory
! ******************************************************************************
          
    CALL RFLU_NullifyCell2FaceList(pRegion)
          
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying cell-to-face list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_DestroyCell2FaceList










  
! ******************************************************************************
!
! Purpose: Destroy face list.
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
  
  SUBROUTINE RFLU_DestroyFaceList(pRegion)

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
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global  

    CALL RegisterFunction(global,'RFLU_DestroyFaceList',&
  'RFLU_ModFaceList.F90')        

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying face lists...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate and nullify memory for interior face lists
! ******************************************************************************

    DEALLOCATE(pGrid%f2c,STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2c')
    END IF ! global%error  

    DEALLOCATE(pGrid%f2v,STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2v')
    END IF ! global%error   
  
! ******************************************************************************
!   Deallocate and nullify memory for boundary face lists. NOTE need to set 
!   renumbering flag to false because when have symmetry or periodic patches
!   need to recreate boundary face lists after adding virtual cells on patches
!   and renumber them. If not reset here, lists do not get renumbered because
!   will get skipped on account of flag already indicating renumbering having 
!   been done.
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DEALLOCATE(pPatch%bf2c,STAT=errorFlag) 
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2c')
      END IF ! global%error                    

      DEALLOCATE(pPatch%bf2v,STAT=errorFlag) 
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2v')
      END IF ! global%error  
      
      pPatch%renumFlag = .FALSE.                                    
    END DO ! iPatch     
    
! ******************************************************************************
!   Nullify memory
! ******************************************************************************
  
    CALL RFLU_NullifyFaceList(pRegion)
    
! ******************************************************************************
!   End
! ******************************************************************************
      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying face lists done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)    

  END SUBROUTINE RFLU_DestroyFaceList    








! ******************************************************************************
!
! Purpose: Get faces opposing given face. 
!
! Description: Compare vertices of given face with vertices of the faces of the 
!   two cells adjacent to given face. 
!
! Input:
!   pRegion     Pointer to region
!   iPatch      Patch index
!   ifg         Face index
!
! Output: 
!   nFacesOpp   Number of opposing faces
!   oppFaceInfo Information about opposing faces
!
! Notes:
!   1. There are in general two faces opposing a given face because two cells
!      share a face. For faces of a boundary cell, there will be at most one 
!      opposing face. 
!   2. The term opposing is used in the sense that two faces of a cell are 
!      opposing if they do not share any vertices. Hence tetrahedra or pyramids
!      do not have any opposing faces. 
!   3. The array oppFaceInfo contains two entries for each opposing face. The
!      first entry is the patch number of the opposing face (0 if the opposing
!      face is an internal face), and the second entry is the index of the face 
!      in the respective face list (internal face list or patch face list). If 
!      a given face has no or only one opposing face, the corresponding entries
!      are set to OPP_FACE_NONE.
!   4. Looping down when determining whether global and local face are identical
!      because missing vertex (which could have entry 0 depending on parameter 
!      VERT_NONE) would mean that could only check for future unnecessary 
!      comparisons for local face 2, so by starting from end can save one
!      comparison.
!
! ******************************************************************************

  SUBROUTINE RFLU_GetOpposingFaces(pRegion,iPatch,ifg,nFacesOpp,faceOppInfo)

    USE ModSortSearch

    USE RFLU_ModGrid

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: ifg,iPatch
    INTEGER, INTENT(OUT) :: nFacesOpp
    INTEGER, INTENT(OUT) :: faceOppInfo(2,2)
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: ic,icg,icl,ict,ifgOpp,ifl,iflOpp,iPatchOpp,ivl,nCells,term  
    INTEGER :: c(2),f2vSort(4),f2vSort2(4)
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_GetOpposingFaces',&
  'RFLU_ModFaceList.F90')

    pGrid => pRegion%grid      

! ******************************************************************************
!   Set patch pointer if face is located on a patch
! ******************************************************************************

    IF ( iPatch > 0 ) THEN 
      pPatch => pRegion%patches(iPatch)        
    ELSE IF ( iPatch == 0 ) THEN 
      NULLIFY(pPatch)
    ELSE  
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! iPatch

! ******************************************************************************
!   Get cells sharing face and the vertices of face
! ******************************************************************************

    IF ( iPatch == 0 ) THEN ! Interior face
      nCells = 2

      c(1) = pGrid%f2c(1,ifg)
      c(2) = pGrid%f2c(2,ifg)

      DO ivl = 1,4
        f2vSort(ivl) = pGrid%f2v(ivl,ifg) 
      END DO ! ivl
    ELSE IF ( iPatch > 0 ) THEN ! Boundary face
      nCells = 1

      c(1) = pPatch%bf2c(ifg)

      DO ivl = 1,4
        IF ( pPatch%bf2v(ivl,ifg) /= VERT_NONE ) THEN 
          f2vSort(ivl) = pPatch%bv(pPatch%bf2v(ivl,ifg))
        ELSE 
          f2vSort(ivl) = VERT_NONE
        END IF ! pPatch%bf2v
      END DO ! ivl        
    END IF ! iPatch

    CALL QuickSortInteger(f2vSort,4)

! ******************************************************************************
!   Find local face of cells which matches vertices of face in question
! ******************************************************************************

    nFacesOpp = 0

    DO ic = 1,nCells
      icg = c(ic)

      ict = pGrid%cellGlob2Loc(1,icg) ! cell type
      icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

      iPatchOpp = OPP_FACE_NONE
      ifgOpp    = OPP_FACE_NONE

! ==============================================================================
!     Select cell type 
! ==============================================================================

      SELECT CASE ( ict )         

! ------------------------------------------------------------------------------
!       Tetrahedron and pyramid (do not have opposing faces)
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_TET, CELL_TYPE_PYR )

! ------------------------------------------------------------------------------
!       Hexahedron
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_HEX )
          hexFaceLoop: DO ifl = 1,6
            DO ivl = 1,4 ! Get vertices of local face
              f2vSort2(ivl) = pGrid%hex2v(f2vHex(ivl,ifl),icl) 
            END DO ! ivl
            
            CALL QuickSortInteger(f2vSort2,4)

            term = 0

            DO ivl = 4,1,-1 ! See note about looping down
              term = term + ABS(f2vSort(ivl) - f2vSort2(ivl)) 
              
              IF ( term /= 0 ) THEN 
                CYCLE hexFaceLoop
              END IF ! term
            END DO ! ivl

            IF ( term == 0 ) THEN ! Found opposing face
              iflOpp = f2fOppHex(ifl)
              
              nFacesOpp = nFacesOpp + 1
              
              faceOppInfo(1,nFacesOpp) = pGrid%hex2f(1,iflOpp,icl)
              faceOppInfo(2,nFacesOpp) = pGrid%hex2f(2,iflOpp,icl)              
              
              EXIT hexFaceLoop
            END IF ! term
          END DO hexFaceLoop

! ------------------------------------------------------------------------------
!       Prism
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_PRI ) 
          priFaceLoop: DO ifl = 1,5 ! Get vertices of local face
            DO ivl = 1,4 
              IF ( f2vPri(ivl,ifl) /= VERT_NONE ) THEN 
                f2vSort2(ivl) = pGrid%pri2v(f2vPri(ivl,ifl),icl) 
              ELSE 
                f2vSort2(ivl) = VERT_NONE
              END IF ! f2vPri                   
            END DO ! ivl

            CALL QuickSortInteger(f2vSort2,4)

            term = 0

            DO ivl = 4,1,-1 ! See note about looping down
              term = term + ABS(f2vSort(ivl) - f2vSort2(ivl)) 
              
              IF ( term /= 0 ) THEN 
                CYCLE priFaceLoop
              END IF ! term
            END DO ! ivl

            IF ( term == 0 ) THEN ! Found opposing face
              iflOpp = f2fOppPri(ifl)

              nFacesOpp = nFacesOpp + 1
              
              faceOppInfo(1,nFacesOpp) = pGrid%pri2f(1,iflOpp,icl)
              faceOppInfo(2,nFacesOpp) = pGrid%pri2f(2,iflOpp,icl)              

              EXIT priFaceLoop
            END IF ! term
          END DO priFaceLoop          

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! ict
    END DO ! ic

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GetOpposingFaces



    
    
    
    
! ******************************************************************************
!
! Purpose: Insert face into cell-to-face list.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iPatch      Patch index
!   icg         Cell index
!   ifg         Face index
!
! Output: None.
!
! Notes: 
!   1. Looping down when determining whether global and local face are identical
!      because missing vertex (which could have entry 0 depending on parameter 
!      VERT_NONE) would mean that could only check for future unnecessary 
!      comparisons for local face 2, so by starting from end can save one
!      comparison.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_InsertIntoCell2FaceList(pRegion,iPatch,icg,ifg)
    
    USE ModSortSearch, ONLY: QuickSortInteger
    
    USE RFLU_ModGrid
    
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: icg,ifg,iPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icl,ict,ifl,ivl,nFaces,term
    INTEGER, DIMENSION(4) :: f2vSort,f2vSort2
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_patch), POINTER :: pPatch   
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_InsertIntoCell2FaceList',&
  'RFLU_ModFaceList.F90')
    
    pGrid => pRegion%grid
    
! ******************************************************************************
!   Get cell type and local cell index
! ******************************************************************************
    
    ict = pGrid%cellGlob2Loc(1,icg) ! cell type
    icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

! ******************************************************************************
!   Get vertices of face
! ******************************************************************************

    IF ( iPatch == 0 ) THEN
      NULLIFY(pPatch)
     
      DO ivl = 1,4
        f2vSort(ivl) = pGrid%f2v(ivl,ifg)
      END DO ! ivl
    ELSE IF ( iPatch > 0 ) THEN 
      pPatch => pRegion%patches(iPatch)
    
      DO ivl = 1,4
        IF ( pPatch%bf2v(ivl,ifg) /= VERT_NONE ) THEN 
          f2vSort(ivl) = pPatch%bv(pPatch%bf2v(ivl,ifg))
        ELSE 
          f2vSort(ivl) = VERT_NONE
        END IF ! pPatch%bf2v
      END DO ! ivl  
    ELSE     
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! iPatch

    CALL QuickSortInteger(f2vSort,4) 

! ******************************************************************************
!   Insert face into cell-to-face list
! ******************************************************************************

    SELECT CASE ( ict )

! ==============================================================================
!     Tetrahedra
! ==============================================================================

      CASE ( CELL_TYPE_TET )
        nFaces = 4
                    
        tetFaceLoop: DO ifl = 1,nFaces                                
          DO ivl = 1,4
            IF ( f2vTet(ivl,ifl) /= VERT_NONE ) THEN 
              f2vSort2(ivl) = pGrid%tet2v(f2vTet(ivl,ifl),icl)          
            ELSE 
              f2vSort2(ivl) = VERT_NONE
            END IF ! f2vTet
          END DO ! ivl

          CALL QuickSortInteger(f2vSort2,4)        

          term = 0

          tetFaceVertLoop: DO ivl = 4,1,-1 ! See note about looping down
            term = term + ABS(f2vSort(ivl) - f2vSort2(ivl))

            IF ( term /= 0 ) THEN 
              EXIT tetFaceVertLoop
            END IF ! term 
          END DO tetFaceVertLoop

          IF ( term == 0 ) THEN
            IF ( pGrid%tet2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%tet2f(2,ifl,icl) == C2F_INIT ) THEN 
              pGrid%tet2f(1,ifl,icl) = iPatch
              pGrid%tet2f(2,ifl,icl) = ifg           

              EXIT tetFaceLoop
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! pGrid%tet2f
          ELSE 
            IF ( ifl == nFaces ) THEN 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! ifl
          END IF ! term          
        END DO tetFaceLoop                 

! ==============================================================================
!     Hexahedra
! ==============================================================================

      CASE ( CELL_TYPE_HEX ) 
        nFaces = 6

        hexFaceLoop: DO ifl = 1,nFaces                                
          DO ivl = 1,4
            f2vSort2(ivl) = pGrid%hex2v(f2vHex(ivl,ifl),icl)          
          END DO ! ivl

          CALL QuickSortInteger(f2vSort2,4)        

          term = 0

          hexFaceVertLoop: DO ivl = 4,1,-1 ! See note about looping down
            term = term + ABS(f2vSort(ivl) - f2vSort2(ivl))

            IF ( term /= 0 ) THEN 
              EXIT hexFaceVertLoop
            END IF ! term 
          END DO hexFaceVertLoop

          IF ( term == 0 ) THEN
            IF ( pGrid%hex2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%hex2f(2,ifl,icl) == C2F_INIT ) THEN 
              pGrid%hex2f(1,ifl,icl) = iPatch
              pGrid%hex2f(2,ifl,icl) = ifg           

              EXIT hexFaceLoop
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! pGrid%hex2f
          ELSE 
            IF ( ifl == nFaces ) THEN 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! ifl            
          END IF ! term          
        END DO hexFaceLoop            

! ==============================================================================
!     Prisms
! ==============================================================================

      CASE ( CELL_TYPE_PRI )
        nFaces = 5
                       
        priFaceLoop: DO ifl = 1,nFaces                                
          DO ivl = 1,4
            IF ( f2vPri(ivl,ifl) /= VERT_NONE ) THEN 
              f2vSort2(ivl) = pGrid%pri2v(f2vPri(ivl,ifl),icl)          
            ELSE 
              f2vSort2(ivl) = VERT_NONE
            END IF ! f2vPri
          END DO ! ivl

          CALL QuickSortInteger(f2vSort2,4)        

          term = 0

          priFaceVertLoop: DO ivl = 4,1,-1 ! See note about looping down
            term = term + ABS(f2vSort(ivl) - f2vSort2(ivl))

            IF ( term /= 0 ) THEN 
              EXIT priFaceVertLoop
            END IF ! term 
          END DO priFaceVertLoop

          IF ( term == 0 ) THEN
            IF ( pGrid%pri2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%pri2f(2,ifl,icl) == C2F_INIT ) THEN 
              pGrid%pri2f(1,ifl,icl) = iPatch
              pGrid%pri2f(2,ifl,icl) = ifg           

              EXIT priFaceLoop
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! pGrid%pri2f
          ELSE 
            IF ( ifl == nFaces ) THEN 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! ifl                        
          END IF ! term          
        END DO priFaceLoop              

! ==============================================================================
!     Pyramids
! ==============================================================================

      CASE ( CELL_TYPE_PYR ) 
        nFaces = 5
      
        pyrFaceLoop: DO ifl = 1,nFaces                                
          DO ivl = 1,4
            IF ( f2vPyr(ivl,ifl) /= VERT_NONE ) THEN
              f2vSort2(ivl) = pGrid%pyr2v(f2vPyr(ivl,ifl),icl)          
            ELSE 
              f2vSort2(ivl) = VERT_NONE
            END IF ! f2vPyr
          END DO ! ivl

          CALL QuickSortInteger(f2vSort2,4)        

          term = 0

          pyrFaceVertLoop: DO ivl = 4,1,-1 ! See note about looping down
            term = term + ABS(f2vSort(ivl) - f2vSort2(ivl))

            IF ( term /= 0 ) THEN 
              EXIT pyrFaceVertLoop
            END IF ! term 
          END DO pyrFaceVertLoop

          IF ( term == 0 ) THEN
            IF ( pGrid%pyr2f(1,ifl,icl) == C2F_INIT .OR. & 
                 pGrid%pyr2f(2,ifl,icl) == C2F_INIT ) THEN 
              pGrid%pyr2f(1,ifl,icl) = iPatch
              pGrid%pyr2f(2,ifl,icl) = ifg           

              EXIT pyrFaceLoop
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! pGrid%pyr2f
          ELSE 
            IF ( ifl == nFaces ) THEN 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! ifl             
          END IF ! term          
        END DO pyrFaceLoop  
        
! ==============================================================================
!     Default
! ==============================================================================
        
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                      
    END SELECT ! ict     

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)
    
  END SUBROUTINE RFLU_InsertIntoCell2FaceList     
    
    
    
    
    


! ******************************************************************************
!
! Purpose: Nullify actual-virtual-face-to-border lists.
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

  SUBROUTINE RFLU_NullifyAVFace2BorderList(pRegion)

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

    CALL RegisterFunction(global,'RFLU_NullifyAVFace2BorderList',&
  'RFLU_ModFaceList.F90')

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Nullify memory
! ******************************************************************************
          
    NULLIFY(pGrid%avf2b)
          
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_NullifyAVFace2BorderList 

    
    
    
    

! ******************************************************************************
!
! Purpose: Nullify actual-virtual-face-to-patch lists.
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

  SUBROUTINE RFLU_NullifyAVFace2PatchList(pRegion)

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

    CALL RegisterFunction(global,'RFLU_NullifyAVFace2PatchList',&
  'RFLU_ModFaceList.F90')
  
    pGrid => pRegion%grid
          
! ******************************************************************************
!   Nullify memory
! ******************************************************************************
          
    NULLIFY(pGrid%avf2p)
          
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_NullifyAVFace2PatchList     
    
    




! ******************************************************************************
!
! Purpose: Nullify cell-to-face lists.
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

  SUBROUTINE RFLU_NullifyCell2FaceList(pRegion)

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

    CALL RegisterFunction(global,'RFLU_NullifyCell2FaceList',&
  'RFLU_ModFaceList.F90')

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Nullify memory
! ******************************************************************************
          
    NULLIFY(pGrid%tet2f)
    NULLIFY(pGrid%hex2f)
    NULLIFY(pGrid%pri2f)
    NULLIFY(pGrid%pyr2f)
          
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_NullifyCell2FaceList    
    
    




    
    
    
! ******************************************************************************
!
! Purpose: Nullify face list.
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

  SUBROUTINE RFLU_NullifyFaceList(pRegion)    

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
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch  
  
! ******************************************************************************
!   Start
! ******************************************************************************
    
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyFaceList',&
  'RFLU_ModFaceList.F90')

! ******************************************************************************
!   Set grid pointer 
! ******************************************************************************

    pGrid => pRegion%grid
          
! ******************************************************************************
!   Nullify memory 
! ******************************************************************************

    NULLIFY(pGrid%f2c)
    NULLIFY(pGrid%f2v)

! TEMPORARY
!    DO iPatch = 1,pGrid%nPatches
!      pPatch => pRegion%patches(iPatch)
!
!      NULLIFY(pPatch%bf2c) 
!      NULLIFY(pPatch%bf2v) 
!    END DO ! iPatch
! END TEMPORARY

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_NullifyFaceList







    
! ******************************************************************************
!
! Purpose: Reorient faces.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Need to make sure that normals of interpartition faces are identical for 
!      the two domains to which such faces belong. This is because, for 
!      example, the Roe dissipative flux depends on ABS(qh-ah), and this will 
!      give different results for the two domains if qh has different signs in 
!      the two domains, which it will if the normal vector does not point in 
!      the same direction. Here make normal vector point in the same direction
!      as on serial grid.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ReorientFaces(pRegion)

    USE RFLU_ModCellFaceEdgeInfo

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

    INTEGER :: cntr,c1,c1s,c2,c2s,ifg,v1g,v2g,v3g,v4g
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReorientFaces',&
  'RFLU_ModFaceList.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reorienting faces...'  
    END IF ! global%verbLevel  
  
! ******************************************************************************
!   Set grid pointer and initialize counter
! ******************************************************************************

    pGrid => pRegion%grid

    cntr = 0
  
! ******************************************************************************
!   Loop over faces and reorient depending on indices in serial grid
! ******************************************************************************

    DO ifg = 1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

      c1s = pGrid%pc2sc(c1)
      c2s = pGrid%pc2sc(c2)

      IF ( c1s > c2s ) THEN
        pGrid%f2c(1,ifg) = c2
        pGrid%f2c(2,ifg) = c1
                                                                                                                      
        v1g = pGrid%f2v(1,ifg)
        v2g = pGrid%f2v(2,ifg)
        v3g = pGrid%f2v(3,ifg)
        v4g = pGrid%f2v(4,ifg)
                                                                                                                      
        pGrid%f2v(1,ifg) = v3g
        pGrid%f2v(2,ifg) = v2g
        pGrid%f2v(3,ifg) = v1g
        pGrid%f2v(4,ifg) = v4g
                                                                                                                    
        cntr = cntr + 1
      END IF ! c1s
    END DO ! ifg

! ******************************************************************************
!   Write info
! ******************************************************************************
      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, & 
                                     'Number of reoriented faces:',cntr
    END IF ! global%verbLevel      
      
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reorienting faces done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReorientFaces  


   
    
    

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModFaceList


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModFaceList.F90,v $
! Revision 1.43  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.42  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.41  2006/08/18 14:01:48  haselbac
! Added routines for AVFace2Patch list, removed Nullify routines from IF
!
! Revision 1.40  2006/03/25 21:52:12  haselbac
! Substantial changes because of sype patches
!
! Revision 1.39  2006/03/20 13:54:04  haselbac
! Added output of region index
!
! Revision 1.38  2005/05/18 22:11:10  fnajjar
! ACH: Fixed bug in building of av-face list, added setting of nFacesAV
!
! Revision 1.37  2005/04/29 23:00:11  haselbac
! Added routines for AVFace2Border list
!
! Revision 1.36  2005/04/15 15:06:52  haselbac
! Removed Charm/FEM stuff and changed routine to reorient faces
!
! Revision 1.35  2004/12/29 21:07:33  haselbac
! Removed setting of pGrid%nBFaces, now done when reading dims
!
! Revision 1.34  2004/12/04 03:29:48  haselbac
! Cosmetics only
!
! Revision 1.33  2004/11/09 00:28:37  haselbac
! Cosmetics only
!
! Revision 1.32  2004/11/03 15:04:53  haselbac
! Changed correction of estimation of boundary faces, cosmetics
!
! Revision 1.31  2004/10/19 19:40:37  haselbac
! Adapted to changes in boundary connectivity, GEN3
!
! Revision 1.30  2004/10/07 14:55:55  haselbac
! Bug fix: Prism limit for pyramid loop
!
! Revision 1.29  2004/10/06 15:35:19  haselbac
! Bug fixes: VERT_NONE was used to access face and cell conn
!
! Revision 1.28  2004/09/27 02:43:22  haselbac
! Bug fix in RFLU_InsertIntoCell2FaceList
!
! Revision 1.27  2004/09/27 01:39:03  haselbac
! Added proper headers and routine for getting opposing faces
!
! Revision 1.26  2004/07/06 15:14:38  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!                                                
! Revision 1.25  2004/06/16 20:01:01  haselbac                                  
! Now store nBFaces directly in grid data type                                  
!
! Revision 1.24  2004/01/22 16:03:59  haselbac                                  
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC
! and titan
!
! Revision 1.23  2003/12/04 03:28:47  haselbac                                  
! Cleaned up                                                                    
!
! Revision 1.22  2003/11/03 03:49:51  haselbac                                  
! Removed building and deallocation of bf2bg list                               
!
! Revision 1.21  2003/08/20 02:09:58  haselbac                                  
! Changed verbosity conditions to reduce solver output in GENx runs             
!
! Revision 1.20  2003/08/19 22:48:43  haselbac                                  
! Added explicit initialization (to avoid Frost problems)                       
!
! Revision 1.19  2003/06/04 22:08:30  haselbac                                  
! Added Nullify routines, some cosmetics                                        
!
! Revision 1.18  2003/04/16 19:14:24  mtcampbe                                  
! ACH: Changed implied to explicit loops bcos of Frost problems                 
!
! Revision 1.17  2003/04/07 14:24:45  haselbac                                  
! Added cell-to-face lists                                                      
!
! Revision 1.16  2003/04/02 17:27:22  haselbac                                  
! Changed limits of modified estimate of no of faces                            
!
! Revision 1.15  2003/03/25 19:14:46  haselbac                                  
! Added new subroutine for reorienting AV faces                                 
!
! Revision 1.14  2003/03/15 18:07:29  haselbac                                  
! Added creation, removed face splitting, numerous other fixes                  
!
! Revision 1.13  2003/01/28 16:29:53  haselbac                                  
! Added creation, removed renumbering (bcos of RFLU_InitFlowSolver changes),    
! and more                                                                      
!
! Revision 1.12  2002/10/27 19:04:50  haselbac                                  
! Bug fix, added CHECK_DATASTRUCT output, and use explicit copies (ASCI White)  
!
! Revision 1.11  2002/10/08 15:49:21  haselbac                                  
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem            
!
! Revision 1.10  2002/09/10 20:26:24  haselbac                                  
! Corrected bug in RFLU_BuildExtraBFaceList                                     
!
! Revision 1.9  2002/09/09 15:05:09  haselbac                                   
! global now under regions, compute nBFaces, added construction of bf2bg access
! array
!
! Revision 1.8  2002/07/27 18:10:39  haselbac                                   
! More conservative estimate of nFacesEst when no boundaries present, got 
! overflow
!
! Revision 1.7  2002/07/25 15:00:19  haselbac                                   
! Only write out for MASTERPROC, more output for CHECK_DATASTRUCT               
!
! Revision 1.6  2002/06/27 15:49:47  haselbac                                   
! Modifications for parallelization, flagging of boundary faces, added
! CHECK_DATASTRUCT flag
!
! Revision 1.5  2002/06/17 13:39:45  haselbac                                   
! Prefixed SOLVER_NAME to all screen output                                     
!
! Revision 1.4  2002/06/14 20:13:18  haselbac                                   
! Added destroy flag, routines for parallelization (for Charm)                  
!
! Revision 1.3  2002/06/05 18:51:58  haselbac                                   
! Cosmetic change only, added empty line after destroying face list             
!
! Revision 1.2  2002/05/04 16:39:51  haselbac                                   
! Cosmetic changes                                                              
!
! Revision 1.1  2002/04/11 18:48:48  haselbac                                   
! Initial revision                                                              
!
! ******************************************************************************

























