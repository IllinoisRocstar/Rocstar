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
! Purpose: Collection of routines for symmetry and periodic patch operations.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModSymmetryPeriodic.F90,v 1.8 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModSymmetryPeriodic

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModBorder, ONLY: t_border
  USE ModError
  USE ModMPI

  USE ModSortSearch

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModSymmetryPeriodic.F90,v $ $Revision: 1.8 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: RFLU_SYPE_AddVirtualCells, &
            RFLU_SYPE_BuildP2VCList, &
            RFLU_SYPE_BuildP2VCListSerial, &     
            RFLU_SYPE_BuildTransforms, &  
            RFLU_SYPE_BuildVertexMaps, &
            RFLU_SYPE_CreateVertexMaps, &
            RFLU_SYPE_DestroyP2VCList, &                 
            RFLU_SYPE_GetActualSerialCell, & 
            RFLU_SYPE_GetRelatedVertex, &  
            RFLU_SYPE_HaveSyPePatches, & 
            RFLU_SYPE_ReadTransforms, & 
            RFLU_SYPE_SetSyPePatchesFlag, & 
            RFLU_SYPE_WriteTransforms

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_SYPE_AddVirtualCellsInv2


! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  

! ******************************************************************************
!
! Purpose: Add virtual cells arising from sype patches.
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

SUBROUTINE RFLU_SYPE_AddVirtualCells(pRegion)                                   

  USE RFLU_ModHashTable

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

  INTEGER :: errorFlag,i,icg,icg2,icl,ict,iLayer,iLoc,iPatch,iReg,j,key, &
             nCellsVirt,nCellsVirtMax,nLayers
  INTEGER, DIMENSION(:), ALLOCATABLE :: vc
  TYPE(t_global), POINTER :: global   
  TYPE(t_grid), POINTER :: pGrid      
  TYPE(t_patch), POINTER :: pPatch,pPatchRelated    

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_AddVirtualCells',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Adding s/p virtual cells...'
  END IF ! global%verbLevel

  IF ( global%verbLevel > VERBOSE_LOW ) THEN        
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal          
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers 
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Allocate temporary for virtual cells
! ******************************************************************************

  nCellsVirtMax = pGrid%nCellsMax - pGrid%nCells

  ALLOCATE(vc(nCellsVirtMax),STAT=errorFlag)
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vc')
  END IF ! global%error 

! ******************************************************************************
! Add virtual cells
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%bcType ) 
      CASE ( BC_SYMMETRY )       
        IF ( global%verbLevel > VERBOSE_LOW ) THEN        
          WRITE(STDOUT,'(A,5X,A,2(1X,I3))') SOLVER_NAME,'Patch:',iPatch, & 
                                            pPatch%iPatchGlobal          
        END IF ! global%verbLevel      
      
        SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
          CASE ( 1 ) 
            CALL RFLU_SYPE_AddVirtualCellsInv1(pRegion,pPatch,vc, &
                                               nCellsVirtMax,nCellsVirt)
          CASE ( 2 ) 
            CALL RFLU_SYPE_AddVirtualCellsInv2(pRegion,pPatch,vc, & 
                                               nCellsVirtMax,nCellsVirt)
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%spaceOrder                                     
                                               
        CALL RFLU_SYPE_AugmentCellLists(pRegion,pPatch,vc,nCellsVirt)
      CASE ( BC_PERIODIC ) 
        pPatchRelated => pRegion%patches(pPatch%iPatchRelated)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN        
          WRITE(STDOUT,'(A,5X,A,2(1X,I3))') SOLVER_NAME,'Patch:',iPatch, & 
                                            pPatch%iPatchGlobal          
        END IF ! global%verbLevel

        SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
          CASE ( 1 ) 
            CALL RFLU_SYPE_AddVirtualCellsInv1(pRegion,pPatchRelated,vc, &
                                               nCellsVirtMax,nCellsVirt)
          CASE ( 2 ) 
            CALL RFLU_SYPE_AddVirtualCellsInv2(pRegion,pPatchRelated,vc, & 
                                               nCellsVirtMax,nCellsVirt)
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%spaceOrder    
                                                         
        CALL RFLU_SYPE_AugmentCellLists(pRegion,pPatch,vc,nCellsVirt)                                           
    END SELECT ! pPatch%bcType
  END DO ! iPatch             

! ******************************************************************************
! Deallocate temporary memory for virtual cells
! ******************************************************************************

  DEALLOCATE(vc,STAT=errorFlag)
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'virtCells')
  END IF ! global%error 

! ******************************************************************************
! Write information about numbers of virtual cells
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Virtual cell statistics:'
    WRITE(STDOUT,'(A,5X,A,1X,I6)') SOLVER_NAME,'Tetrahedra:', &
                                   pGrid%nTetsTot-pGrid%nTets
    WRITE(STDOUT,'(A,5X,A,1X,I6)') SOLVER_NAME,'Hexahedra: ', &
                                   pGrid%nHexsTot-pGrid%nHexs
    WRITE(STDOUT,'(A,5X,A,1X,I6)') SOLVER_NAME,'Prisms:    ', &
                                   pGrid%nPrisTot-pGrid%nPris
    WRITE(STDOUT,'(A,5X,A,1X,I6)') SOLVER_NAME,'Pyramids:  ', &
                                   pGrid%nPyrsTot-pGrid%nPyrs 
  END IF ! global%verbLevel 

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Adding s/p virtual cells done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SYPE_AddVirtualCells
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Add virtual cells for inviscid fluxes based on first-order scheme.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pRegionSerial       Pointer to serial region        
!   vc                  List of virtual cells (empty)
!   nCellsVirtMax       Maximum allowable number of virtual cells  
!
! Output: 
!   vc                  List of virtual cells 
!   nCellsVirt          Number of virtual cells
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_SYPE_AddVirtualCellsInv1(pRegion,pPatch,vc,nCellsVirtMax, & 
                                         nCellsVirt)                                   

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: nCellsVirtMax
  INTEGER, INTENT(OUT) :: nCellsVirt
  INTEGER, INTENT(INOUT) :: vc(nCellsVirtMax)
  TYPE(t_patch), POINTER :: pPatch    
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifl
  TYPE(t_global), POINTER :: global     

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_AddVirtualCellsInv1',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Adding s/p virtual cells for inviscid first-order stencil...' 
  END IF ! global%myProcid

! ******************************************************************************
! Set pointers and initialize variables
! ******************************************************************************

  nCellsVirt = 0
                                  
! ******************************************************************************
! Loop over faces on this patch and add adjacent cells
! ******************************************************************************

  DO ifl = 1,pPatch%nBFaces
    IF ( nCellsVirt < nCellsVirtMax ) THEN     
      nCellsVirt = nCellsVirt + 1
      vc(nCellsVirt) = pPatch%bf2c(ifl)      
    ELSE 
      CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vc') 
    END IF ! nCellsVirt
  END DO ! ifl                         

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Adding s/p virtual cells for first-order stencil done.' 
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SYPE_AddVirtualCellsInv1  
  
  
  
  
  
    
  

! ******************************************************************************
!
! Purpose: Add virtual cells for inviscid fluxes based on higher-order scheme.
!
! Description: Building list of virtual cells proceeds in several steps:
!   1. Build list of vertices from list of actual-virtual faces
!   2. Build list of source cells adjacent to actual-virtual faces which are 
!      in same region as that for which virtual cells are constructed. One 
!      layer of source cells is sufficient. If more layers of source cells 
!      are specified here, this should not change the number of virtual cells.
!   3. Build stencils for source cells and add stencil members to list of 
!      virtual cells if not in same region
!   4. Loop over existing virtual cells and add layers.
!
! Input:
!   pRegion             Pointer to region
!   pRegionSerial       Pointer to serial region        
!   vc                  List of virtual cells (empty)
!   nCellsVirtMax       Maximum allowable number of virtual cells  
!
! Output: 
!   vc                  List of virtual cells 
!   nCellsVirt          Number of virtual cells
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_SYPE_AddVirtualCellsInv2(pRegion,pPatch,vc,nCellsVirtMax, &
                                         nCellsVirt)                                   

  USE ModSortSearch

  USE RFLU_ModStencilsCells, ONLY: RFLU_BuildC2CStencilWrapper, & 
                                   RFLU_CreateC2CStencilWrapper, & 
                                   RFLU_DestroyC2CStencilWrapper, & 
                                   RFLU_SetInfoC2CStencilWrapper
  USE RFLU_ModTopologyUtils, ONLY: RFLU_BuildFaceVertList, & 
                                   RFLU_BuildVertCellNghbList

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: nCellsVirtMax
  INTEGER, INTENT(OUT) :: nCellsVirt
  INTEGER, INTENT(INOUT) :: vc(nCellsVirtMax)
  TYPE(t_patch), POINTER :: pPatch    
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,i,icg,icg2,iflBeg,iflEnd,iLayer,iLoc,j,nLayers, &
             nVert,scDim,vcNewDim,vcNewDimMax,vcOldDim
  INTEGER, DIMENSION(:), ALLOCATABLE :: avv,sc,vcNew,vcOld
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global     

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_AddVirtualCellsInv2',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Adding s/p virtual cells for inviscid higher-order stencil...' 
  END IF ! global%myProcid

! ******************************************************************************
! Set pointers and initialize variables
! ******************************************************************************

  pGrid => pRegion%grid    

  nCellsVirt = 0

! ******************************************************************************
! Based on list of patch vertices, build list of source cells adjacent to 
! actual-virtual faces for construction of virtual cells. One layer of source 
! cells is sufficient. If more layers of source cells are specified here, this
! should not change the number of virtual cells.
! ******************************************************************************

  ALLOCATE(sc(nCellsVirtMax),STAT=errorFlag)
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'sc')
  END IF ! global%error             

  nLayers = 1    

  CALL RFLU_BuildVertCellNghbList(global,pGrid,pPatch%bv,pPatch%nBVert, & 
                                  nLayers,pRegion%iRegionGlobal,sc, &
                                  nCellsVirtMax,scDim)
                                  
! ******************************************************************************
! Loop over list of source cells, build stencil for each cell, and add cells
! in stencil as virtual cells if not in same region and not already added
! ******************************************************************************

  CALL RFLU_SetInfoC2CStencilWrapper(pRegion,pRegion%mixtInput%spaceOrder-1)        
  CALL RFLU_CreateC2CStencilWrapper(pRegion)   
  
  DO i = 1,scDim
    icg = sc(i)

    CALL RFLU_BuildC2CStencilWrapper(pRegion,icg,CONSTR_NONE)

    DO j = 1,pGrid%c2cs(icg)%nCellMembs
      icg2 = pGrid%c2cs(icg)%cellMembs(j)

      IF ( nCellsVirt > 0 ) THEN 
        CALL BinarySearchInteger(vc(1:nCellsVirt),nCellsVirt,icg2,iLoc)
      ELSE 
        iLoc = ELEMENT_NOT_FOUND
      END IF ! nCellsVirt

      IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
        IF ( nCellsVirt < nCellsVirtMax ) THEN
          nCellsVirt = nCellsVirt + 1       
          vc(nCellsVirt) = icg2
        ELSE
          CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vc')           
        END IF ! nCellsVirt
        
        IF ( nCellsVirt > 1 ) THEN    
          CALL QuickSortInteger(vc(1:nCellsVirt),nCellsVirt)              
        END IF ! nCellsVirt 
      END IF ! iLoc                    
    END DO ! j      
  END DO ! i      

  DEALLOCATE(sc,STAT=errorFlag)
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'sc')
  END IF ! global%error             

! ******************************************************************************
! Loop over current list of virtual cells and add layers of cells. For the 
! current scheme, two additional layers are required.
! ******************************************************************************

! ==============================================================================
! Allocate temporary memory 
! ==============================================================================

  vcOldDim    = nCellsVirt
  vcNewDimMax = nCellsVirtMax

  ALLOCATE(vcOld(vcOldDim),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vcOld')
  END IF ! global%error

  DO i = 1,vcOldDim
    vcOld(i) = vc(i)
  END DO ! i

  ALLOCATE(vcNew(vcNewDimMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vcNew')
  END IF ! global%error

! ==============================================================================
! Loop over layers
! ==============================================================================

  nLayers = 1

  DO iLayer = 1,nLayers
    vcNewDim = 0            

    DO i = 1,vcOldDim
      icg = vcOld(i)

      CALL RFLU_BuildC2CStencilWrapper(pRegion,icg,CONSTR_NONE)

      DO j = 1,pGrid%c2cs(icg)%nCellMembs
        icg2 = pGrid%c2cs(icg)%cellMembs(j)

        CALL BinarySearchInteger(vc(1:nCellsVirt),nCellsVirt,icg2,iLoc)

        IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
          IF ( nCellsVirt < nCellsVirtMax ) THEN  
            nCellsVirt = nCellsVirt + 1              
            vc(nCellsVirt) = icg2
          ELSE 
            CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vc')          
          END IF ! nCellsVirt

          IF ( vcNewDim < vcNewDimMax ) THEN 
            vcNewDim = vcNewDim + 1              
            vcNew(vcNewDim) = icg2
          ELSE
            CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vcNew')         
          END IF ! vcNewDim 

          CALL QuickSortInteger(vc(1:nCellsVirt),nCellsVirt) 
        END IF ! iLoc
      END DO ! j                      
    END DO ! i

    DEALLOCATE(vcOld,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vcOld')
    END IF ! global%error

    IF ( iLayer /= nLayers ) THEN 
      vcOldDim = vcNewDim

      ALLOCATE(vcOld(vcOldDim),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vcOld')
      END IF ! global%error    

      DO i = 1,vcOldDim
        vcOld(i) = vcNew(i)
      END DO ! i
    END IF ! iLayer      
  END DO ! iLayer   

! ==============================================================================
! Deallocate temporary memory 
! ==============================================================================

  DEALLOCATE(vcNew,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vcNew')
  END IF ! global%error

  CALL RFLU_DestroyC2CStencilWrapper(pRegion)                           

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Adding s/p virtual cells for inviscid higher-order stencil done.' 
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SYPE_AddVirtualCellsInv2







! ******************************************************************************
!
! Purpose: Add virtual cells to cell and patch data structures.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pRegionSerial       Pointer to serial region        
!   vc                  List of virtual cells 
!   nCellsVirt          Number of virtual cells  
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_SYPE_AugmentCellLists(pRegion,pPatch,vc,nCellsVirt)                              

  USE ModTools, ONLY: SwapIntegers, & 
                      SwapRFREALs

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: nCellsVirt
  INTEGER, INTENT(IN) :: vc(nCellsVirt)
  TYPE(t_patch), POINTER :: pPatch    
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,iBorder,icg,icg2,icl,ict,icv,ifl,iflBeg,iflEnd,ifl2, & 
             iLayer,iLoc,iLoc2,iPatchRelated,iPatch2,ivg,ivgm,ivl,ivlm,ivl2, & 
             j,nLayers,nVert,nVertVirtEst,scDim,vcNewDim,vcNewDimMax,vcOldDim
  INTEGER, DIMENSION(:), ALLOCATABLE :: ivgRecvSorted,ivgSendSorted, & 
                                        ivgSharedSorted
  TYPE(t_border), POINTER :: pBorder,pBorderRelated
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global  
  TYPE(t_patch), POINTER :: pPatchRelated,pPatch2,pPatch3         

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_AugmentCellLists',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Adding s/p virtual cells to data structure...' 
  END IF ! global%myProcid

! ******************************************************************************
! Set pointers and initialize variables
! ******************************************************************************

  pGrid => pRegion%grid    

  nVertVirtEst = 8*nCellsVirt

  pBorder => pGrid%borders(pPatch%iBorder)
    
  iPatchRelated  =  pPatch%iPatchRelated    
  pPatchRelated  => pRegion%patches(iPatchRelated)        
    
  IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
    pBorderRelated => pGrid%borders(pPatch%iBorder)
  ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN 
    pBorderRelated => pGrid%borders(pPatchRelated%iBorder)  
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! pPatch%bcType
  
  pBorder%iBorder       = pPatchRelated%iBorder    
  pBorder%iRegionGlobal = pRegion%iRegionGlobal
  
  pBorderRelated%nCellsSend = nCellsVirt
  pBorderRelated%nVertSend  = 0
          
  pBorder%nCellsRecv  = nCellsVirt
  pBorder%nVertRecv   = 0
  pBorder%nVertShared = 0  
  
  ALLOCATE(pBorderRelated%icgSend(pBorderRelated%nCellsSend),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorderRelated%icgSend')
  END IF ! global%error
  
  ALLOCATE(pBorderRelated%ivgSend(nVertVirtEst),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorderRelated%ivgSend')
  END IF ! global%error
    
  ALLOCATE(pBorder%icgRecv(pBorder%nCellsRecv),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%icgRecv')
  END IF ! global%error  

  ALLOCATE(pBorder%ivgRecv(nVertVirtEst),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgRecv')
  END IF ! global%error

  ALLOCATE(pBorder%ivgShared(nVertVirtEst),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgShared')
  END IF ! global%error

  ALLOCATE(ivgSendSorted(nVertVirtEst),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ivgSendSorted')
  END IF ! global%error
  
  ALLOCATE(ivgRecvSorted(nVertVirtEst),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ivgRecvSorted')
  END IF ! global%error  
  
  ALLOCATE(ivgSharedSorted(nVertVirtEst),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ivgSharedSorted')
  END IF ! global%error  

! ******************************************************************************
! Loop over list of virtual cells, add to connectivity lists according to cell 
! type
! ******************************************************************************

  DO icv = 1,nCellsVirt
    icg = vc(icv)

    ict = pGrid%cellGlob2Loc(1,icg)
    icl = pGrid%cellGlob2Loc(2,icg)

! ==============================================================================
!   Specify connectivity and set local-to-global mapping
! ==============================================================================

    SELECT CASE ( ict )           

! ------------------------------------------------------------------------------
!     Tetrahedra 
! ------------------------------------------------------------------------------

      CASE ( CELL_TYPE_TET )
        pGrid%nCellsTot = pGrid%nCellsTot + 1  
        pGrid%nTetsTot  = pGrid%nTetsTot  + 1

        IF ( pGrid%nTetsTot <= pGrid%nTetsMax ) THEN           
          CALL RFLU_SYPE_AugmentVertexLists(pRegion,pGrid,pPatch,pBorder, &
                                            pBorderRelated,pGrid%tet2v,icl, &
                                            pGrid%nTetsTot,ivgSendSorted, &
                                            ivgRecvSorted,ivgSharedSorted)                                    
                                            
          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_TET
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nTetsTot                

          pGrid%tet2CellGlob(pGrid%nTetsTot) = pGrid%nCellsTot  
          
          IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
            CALL SwapIntegers(pGrid%tet2v(1,pGrid%nTetsTot), & 
                              pGrid%tet2v(2,pGrid%nTetsTot))
          END IF ! pPatch%bcType
        ELSE 
          CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'pGrid%tet2v') 
        END IF ! pGrid%nTetsTot       

! ------------------------------------------------------------------------------
!     Hexahedra 
! ------------------------------------------------------------------------------

      CASE ( CELL_TYPE_HEX ) 
        pGrid%nCellsTot = pGrid%nCellsTot + 1  
        pGrid%nHexsTot  = pGrid%nHexsTot  + 1

        IF ( pGrid%nHexsTot <= pGrid%nHexsMax ) THEN   
          CALL RFLU_SYPE_AugmentVertexLists(pRegion,pGrid,pPatch,pBorder, &
                                            pBorderRelated,pGrid%hex2v,icl, & 
                                            pGrid%nHexsTot,ivgSendSorted, & 
                                            ivgRecvSorted,ivgSharedSorted)                                    
                                        
          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_HEX
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nHexsTot                

          pGrid%hex2CellGlob(pGrid%nHexsTot) = pGrid%nCellsTot  
          
          IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
            CALL SwapIntegers(pGrid%hex2v(2,pGrid%nHexsTot), & 
                              pGrid%hex2v(4,pGrid%nHexsTot))
            CALL SwapIntegers(pGrid%hex2v(6,pGrid%nHexsTot), & 
                              pGrid%hex2v(8,pGrid%nHexsTot)) 
          END IF ! pPatch%bcType
        ELSE 
          CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'pGrid%hex2v') 
        END IF ! pGrid%nHexsTot              

! ------------------------------------------------------------------------------
!     Prisms
! ------------------------------------------------------------------------------

      CASE ( CELL_TYPE_PRI ) 
        pGrid%nCellsTot = pGrid%nCellsTot + 1  
        pGrid%nPrisTot  = pGrid%nPrisTot  + 1

        IF ( pGrid%nPrisTot <= pGrid%nPrisMax ) THEN   
          CALL RFLU_SYPE_AugmentVertexLists(pRegion,pGrid,pPatch,pBorder, &
                                            pBorderRelated,pGrid%pri2v,icl, & 
                                            pGrid%nPrisTot,ivgSendSorted, & 
                                            ivgRecvSorted,ivgSharedSorted)                                    
                    
          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_PRI
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nPrisTot                

          pGrid%pri2CellGlob(pGrid%nPrisTot) = pGrid%nCellsTot  
          
          IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
            CALL SwapIntegers(pGrid%pri2v(2,pGrid%nPrisTot), & 
                              pGrid%pri2v(3,pGrid%nPrisTot))
            CALL SwapIntegers(pGrid%pri2v(5,pGrid%nPrisTot), & 
                              pGrid%pri2v(6,pGrid%nPrisTot))                              
          END IF ! pPatch%bcType
        ELSE 
          CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'pGrid%pri2v') 
        END IF ! pGrid%nPrisTot                         

! ------------------------------------------------------------------------------
!     Pyramids
! ------------------------------------------------------------------------------

      CASE ( CELL_TYPE_PYR ) 
        pGrid%nCellsTot = pGrid%nCellsTot + 1  
        pGrid%nPyrsTot  = pGrid%nPyrsTot  + 1

        IF ( pGrid%nPyrsTot <= pGrid%nPyrsMax ) THEN   
          CALL RFLU_SYPE_AugmentVertexLists(pRegion,pGrid,pPatch,pBorder, &
                                            pBorderRelated,pGrid%pyr2v,icl, &
                                            pGrid%nPyrsTot,ivgSendSorted, &
                                            ivgRecvSorted,ivgSharedSorted)                                    
                    
          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_PYR
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nPyrsTot                

          pGrid%pyr2CellGlob(pGrid%nPyrsTot) = pGrid%nCellsTot  
          
          IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
            CALL SwapIntegers(pGrid%pyr2v(2,pGrid%nPyrsTot), & 
                              pGrid%pyr2v(4,pGrid%nPyrsTot))                             
          END IF ! pPatch%bcType
        ELSE 
          CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'pGrid%pyr2v') 
        END IF ! pGrid%nPyrsTot 

! ------------------------------------------------------------------------------
!     Default
! ------------------------------------------------------------------------------

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)          
    END SELECT ! pGridSerial%cellGlob2Loc          

! ==============================================================================
!   Build cell communication lists
! ==============================================================================

    pBorderRelated%icgSend(icv) = icg
    
    pBorder%icgRecv(icv) = pGrid%nCellsTot 
  END DO ! icl                           

! ******************************************************************************
! Add faces of virtual cells to patch face lists
! ******************************************************************************

  DO iPatch2 = 1,pGrid%nPatches
    pPatch2 => pRegion%patches(iPatch2)
    
    ALLOCATE(pPatch2%bf2cSorted(pPatch2%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2cSorted')
    END IF ! global%error
    
    ALLOCATE(pPatch2%bf2cSortedKeys(pPatch2%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2cSortedKeys')
    END IF ! global%error
    
    ALLOCATE(pPatch2%bf2vSorted(4,pPatch2%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2vSorted')
    END IF ! global%error            
    
    DO ifl = 1,pPatch2%nBFaces
      pPatch2%bf2cSorted(ifl)     = pPatch2%bf2c(ifl)
      pPatch2%bf2cSortedKeys(ifl) = ifl      
    END DO ! ifl
    
    CALL QuickSortIntegerInteger(pPatch2%bf2cSorted,pPatch2%bf2cSortedKeys, &
                                 pPatch2%nBFaces)
                                 
    DO ifl = 1,pPatch2%nBFaces
      ifl2 = pPatch2%bf2cSortedKeys(ifl)
        
      pPatch2%bf2vSorted(1,ifl) = pPatch2%bf2v(1,ifl2)
      pPatch2%bf2vSorted(2,ifl) = pPatch2%bf2v(2,ifl2)
      pPatch2%bf2vSorted(3,ifl) = pPatch2%bf2v(3,ifl2)
      pPatch2%bf2vSorted(4,ifl) = pPatch2%bf2v(4,ifl2)                  
    END DO ! ifl                             
  END DO ! iPatch2

! ******************************************************************************
! Loop over cells to be sent 
! ******************************************************************************

  DO icl = 1,pBorderRelated%nCellsSend
    icg = pBorderRelated%icgSend(icl)

! ==============================================================================
!   Loop over patches         
! ==============================================================================

    DO iPatch2 = 1,pGrid%nPatches
      pPatch2 => pRegion%patches(iPatch2)

! ------------------------------------------------------------------------------
!     If current patch is not the same as patch associated with symmetry or 
!     periodic patch, find out whether cell to be sent is adjacent to current 
!     patch
! ------------------------------------------------------------------------------

      IF ( (pPatch2%iPatchGlobal /= pPatch%iPatchGlobal) .AND. & 
           (pPatch2%nbMap(pPatch%iPatchGlobal) .EQV. .TRUE. ) ) THEN  
        CALL BinarySearchInteger(pPatch2%bf2cSorted(1:pPatch2%nBFaces), & 
                                 pPatch2%nBFaces,icg,iLoc)
                
! ----- Cell found, so adjacent to current patch -------------------------------        
        
        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        
! ------- Triangle face       
        
          IF ( pPatch2%bf2vSorted(4,iLoc) == VERT_NONE ) THEN                            
            IF ( pPatch2%nBTrisTot < pPatch2%nBTrisMax ) THEN
              pPatch2%nBTrisTot  = pPatch2%nBTrisTot  + 1
              pPatch2%nBFacesTot = pPatch2%nBFacesTot + 1
            ELSE 
              CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'tri2v')
            END IF ! pPatch2%nBTrisTot
               
            DO ivl2 = 1,3
              ivl = pPatch2%bf2vSorted(ivl2,iLoc)
              ivg = pPatch2%bv(ivl)
            
              CALL BinarySearchInteger(ivgSendSorted(1:pBorderRelated%nVertSend), & 
                                       pBorderRelated%nVertSend,ivg,iLoc2)

              IF ( iLoc2 /= ELEMENT_NOT_FOUND ) THEN                            
                pPatch2%bTri2v(ivl2,pPatch2%nBTrisTot) = ivgRecvSorted(iLoc2)                   
              ELSE
                IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
                  CALL BinarySearchInteger(ivgSharedSorted(1:pBorder%nVertShared), & 
                                           pBorder%nVertShared,ivg,iLoc2)
                ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN 
                  iPatchRelated = pPatch%iPatchRelated
                    
                  pPatch3 => pRegion%patches(iPatchRelated)
                  
                  CALL BinarySearchInteger(pPatch3%bv(1:pPatch3%nBVert), & 
                                           pPatch3%nBVert,ivg,iLoc2)                                    
                ELSE 
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                END IF ! pPatch%bcType
                                         
                IF ( iLoc2 /= ELEMENT_NOT_FOUND ) THEN
                  IF ( pPatch%bcType == BC_SYMMETRY ) THEN  
                    pPatch2%bTri2v(ivl2,pPatch2%nBTrisTot) = ivg
                  ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN 
                    ivlm = pPatch3%bv2bv(iLoc2)
                    ivgm = pPatch%bv(ivlm) 

                    pPatch2%bTri2v(ivl2,pPatch2%nBTrisTot) = ivgm
                  END IF ! pPatch%bcType                  
                ELSE
                  CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)
                END IF ! iLoc2             
              END IF ! iLoc2
            END DO ! ivl2   
                        
            IF ( pPatch%bcType == BC_SYMMETRY ) THEN
              CALL SwapIntegers(pPatch2%bTri2v(2,pPatch2%nBTrisTot), & 
                                pPatch2%bTri2v(3,pPatch2%nBTrisTot))
            END IF ! pPatch%bcType
            
! ------- Quadrilateral face           
            
          ELSE 
            IF ( pPatch2%nBQuadsTot < pPatch2%nBQuadsMax ) THEN
              pPatch2%nBQuadsTot = pPatch2%nBQuadsTot + 1
              pPatch2%nBFacesTot = pPatch2%nBFacesTot + 1
            ELSE 
              CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'quad2v')
            END IF ! pPatch2%nBQuadsTot
               
            DO ivl2 = 1,4
              ivl = pPatch2%bf2vSorted(ivl2,iLoc)
              ivg = pPatch2%bv(ivl)
            
              CALL BinarySearchInteger(ivgSendSorted(1:pBorderRelated%nVertSend), & 
                                       pBorderRelated%nVertSend,ivg,iLoc2)
            
              IF ( iLoc2 /= ELEMENT_NOT_FOUND ) THEN                            
                pPatch2%bQuad2v(ivl2,pPatch2%nBQuadsTot) = ivgRecvSorted(iLoc2)                    
              ELSE
                IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
                  CALL BinarySearchInteger(ivgSharedSorted(1:pBorder%nVertShared), & 
                                           pBorder%nVertShared,ivg,iLoc2)
                ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN 
                  iPatchRelated = pPatch%iPatchRelated
                    
                  pPatch3 => pRegion%patches(iPatchRelated)
                  
                  CALL BinarySearchInteger(pPatch3%bv(1:pPatch3%nBVert), & 
                                           pPatch3%nBVert,ivg,iLoc2)                                    
                ELSE 
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                END IF ! pPatch%bcType
                                         
                IF ( iLoc2 /= ELEMENT_NOT_FOUND ) THEN
                  IF ( pPatch%bcType == BC_SYMMETRY ) THEN  
                    pPatch2%bQuad2v(ivl2,pPatch2%nBQuadsTot) = ivg
                  ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN 
                    ivlm = pPatch3%bv2bv(iLoc2)
                    ivgm = pPatch%bv(ivlm) 

                    pPatch2%bQuad2v(ivl2,pPatch2%nBQuadsTot) = ivgm
                  END IF ! pPatch%bcType                 
                ELSE
                  CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)
                END IF ! iLoc2
              END IF ! iLoc2
            END DO ! ivl2   
                        
            IF ( pPatch%bcType == BC_SYMMETRY ) THEN
              CALL SwapIntegers(pPatch2%bQuad2v(2,pPatch2%nBQuadsTot), & 
                                pPatch2%bQuad2v(4,pPatch2%nBQuadsTot))
            END IF ! pPatch%bcType
          END IF ! pPatch2%bf2v
        END IF ! iLoc
      END IF ! pPatch2%iPatchGlobal
    END DO ! iPatch
  END DO ! icl

  DO iPatch2 = 1,pGrid%nPatches
    pPatch2 => pRegion%patches(iPatch2)
    
    DEALLOCATE(pPatch2%bf2cSorted,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2cSorted')
    END IF ! global%error
  END DO ! iPatch2

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(ivgSendSorted,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ivgSendSorted')
  END IF ! global%error
  
  DEALLOCATE(ivgRecvSorted,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ivgRecvSorted')
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Adding s/p virtual cells to data structure done.' 
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SYPE_AugmentCellLists
 


  
  
  
! ******************************************************************************
!
! Purpose: Loop over vertices of given cell and update lists relating to 
!   vertices.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pGrid               Pointer to grid 
!   pPatch              Pointer to patch
!   pBorder             Pointer to border
!   pBorderRelated      Pointer to related border
!   x2v                 Cell-to-vertex connectivity
!   icl                 Local cell index
!   nXTot               Total number of cells of given type
!   ivgSendSorted       List of vertices to be sent for given border
!   ivgRecvSorted       List of vertices to be received for given border
!   ivgSharedSorted     List of vertices to be shared for given border
!
! Output: 
!   x2v                 Cell-to-vertex connectivity
!   nXTot               Total number of cells of given type
!   ivgSendSorted       List of vertices to be sent for given border
!   ivgRecvSorted       List of vertices to be received for given border
!   ivgSharedSorted     List of vertices to be shared for given border
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SYPE_AugmentVertexLists(pRegion,pGrid,pPatch,pBorder, &
                                        pBorderRelated,x2v,icl,nXTot, &
                                        ivgSendSorted,ivgRecvSorted, &
                                        ivgSharedSorted)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icl
  INTEGER, INTENT(INOUT) :: nXTot
  INTEGER, DIMENSION(:) :: ivgSendSorted,ivgSharedSorted,ivgRecvSorted
  INTEGER, DIMENSION(:,:), POINTER :: x2v
  REAL(RFREAL), DIMENSION(XCOORD:XYZMAG) :: xyz,xyzt
  TYPE(t_border), POINTER :: pBorder,pBorderRelated
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch 
  TYPE(t_region), POINTER :: pRegion   

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iLoc,iLoc2,iPatchRelated,ivg,ivgm,ivl
  REAL(RFREAL) :: nxpc,nypc,nzpc,xpc,xvg,ypc,yvg,zpc,zvg
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatchRelated

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_AugmentVertexLists',&
  'RFLU_ModSymmetryPeriodic.F90')
 
! ******************************************************************************
! Get patch geometry and set patch index and pointer
! ******************************************************************************

  xpc = pPatch%pc(XCOORD)
  ypc = pPatch%pc(YCOORD)
  zpc = pPatch%pc(ZCOORD)

  nxpc = pPatch%pn(XCOORD)
  nypc = pPatch%pn(YCOORD)
  nzpc = pPatch%pn(ZCOORD) 

  IF ( pPatch%bcType == BC_PERIODIC ) THEN 
    iPatchRelated = pPatch%iPatchRelated
    
    pPatchRelated => pRegion%patches(iPatchRelated)
  ELSE
    iPatchRelated = CRAZY_VALUE_INT
    
    NULLIFY(pPatchRelated)
  END IF ! pPatch%bcType

! ******************************************************************************
! Loop over vertices of cell
! ******************************************************************************

  DO ivl = 1,SIZE(x2v,1)
    ivg = x2v(ivl,icl)

! ==============================================================================
!   Check whether vertex is on patch for symmetry patches or on related patch
!   for periodic patches
! ==============================================================================

    IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
      CALL BinarySearchInteger(pPatch%bv(1:pPatch%nBVert),pPatch%nBVert,ivg, & 
                               iLoc)
    ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN 
      CALL BinarySearchInteger(pPatchRelated%bv(1:pPatchRelated%nBVert), & 
                               pPatchRelated%nBVert,ivg,iLoc)
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! pPatch%bcType

! ==============================================================================
!   Vertex not on patch
! ==============================================================================

    IF ( iLoc == ELEMENT_NOT_FOUND ) THEN                           

! ------------------------------------------------------------------------------
!     Check whether already in list of vertices to be sent
! ------------------------------------------------------------------------------

      IF ( pBorderRelated%nVertSend > 1 ) THEN        
        CALL BinarySearchInteger(ivgSendSorted(1:pBorderRelated%nVertSend), & 
                                 pBorderRelated%nVertSend,ivg,iLoc2)
      ELSE 
        iLoc2 = ELEMENT_NOT_FOUND            
      END IF ! pBorderRelated%nVertSend

! ------------------------------------------------------------------------------
!     Vertex not in list of vertices to be sent, so add. Compute new coordinates 
!     and add vertex to lists
! ------------------------------------------------------------------------------

      IF ( iLoc2 == ELEMENT_NOT_FOUND ) THEN      
        IF ( pGrid%nVertTot < pGrid%nVertMax ) THEN 
          pGrid%nVertTot = pGrid%nVertTot + 1
          x2v(ivl,nXTot) = pGrid%nVertTot      
        ELSE 
          CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'x2v')        
        END IF ! pGrid%nVertTot

        IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
          xvg = pGrid%xyz(XCOORD,ivg)
          yvg = pGrid%xyz(YCOORD,ivg)
          zvg = pGrid%xyz(ZCOORD,ivg)

          CALL ReflectPosition(nxpc,nypc,nzpc,xpc,ypc,zpc,xvg,yvg,zvg)

          pGrid%xyz(XCOORD,pGrid%nVertTot) = xvg
          pGrid%xyz(YCOORD,pGrid%nVertTot) = yvg
          pGrid%xyz(ZCOORD,pGrid%nVertTot) = zvg              
        ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN         
          xyz(XCOORD) = pGrid%xyz(XCOORD,ivg)
          xyz(YCOORD) = pGrid%xyz(YCOORD,ivg)
          xyz(ZCOORD) = pGrid%xyz(ZCOORD,ivg)
          xyz(XYZMAG) = 1.0_RFREAL
                    
          xyzt = MATMUL(pPatchRelated%tm,xyz)
        
          pGrid%xyz(XCOORD,pGrid%nVertTot) = xyzt(XCOORD)
          pGrid%xyz(YCOORD,pGrid%nVertTot) = xyzt(YCOORD)
          pGrid%xyz(ZCOORD,pGrid%nVertTot) = xyzt(ZCOORD)   
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)        
        END IF ! pPatch%bcType

        pBorderRelated%nVertSend = pBorderRelated%nVertSend + 1
        
        pBorderRelated%ivgSend(pBorderRelated%nVertSend) = ivg
        ivgSendSorted(pBorderRelated%nVertSend)          = ivg
        
        pBorder%nVertRecv = pBorder%nVertRecv + 1
        
        pBorder%ivgRecv(pBorder%nVertRecv) = pGrid%nVertTot              
        ivgRecvSorted(pBorder%nVertRecv)   = pGrid%nVertTot

        IF ( pBorder%nVertRecv > 1 ) THEN 
          CALL QuickSortIntegerInteger(ivgSendSorted(1:pBorder%nVertRecv), & 
                                       ivgRecvSorted(1:pBorder%nVertRecv), &
                                       pBorder%nVertRecv)                                       
        END IF ! pBorder%nVertRecv

! ------------------------------------------------------------------------------
!     Vertex already in list of vertices to be sent
! ------------------------------------------------------------------------------

      ELSE 
        x2v(ivl,nXTot) = ivgRecvSorted(iLoc2)        
      END IF ! iLoc2

! ==============================================================================
!   Vertex on patch 
! ==============================================================================

    ELSE     
      IF ( pPatch%bcType == BC_SYMMETRY ) THEN 
        x2v(ivl,nXTot) = ivg
      ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN                        
        CALL RFLU_SYPE_ImposeVertexMap(pRegion,pPatchRelated,ivg,ivgm)    
      
        x2v(ivl,nXTot) = ivgm
      ELSE 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! pPatch%bcType
      
! ------------------------------------------------------------------------------
!     Check whether already in list of vertices to be shared
! ------------------------------------------------------------------------------

      IF ( pBorder%nVertShared > 1 ) THEN 
        IF ( pPatch%bcType == BC_SYMMETRY ) THEN
          CALL BinarySearchInteger(ivgSharedSorted(1:pBorder%nVertShared), & 
                                   pBorder%nVertShared,ivg,iLoc2)
        ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN
          CALL BinarySearchInteger(ivgSharedSorted(1:pBorder%nVertShared), & 
                                   pBorder%nVertShared,ivgm,iLoc2)        
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! pPatch%bcType                  
      ELSE 
        iLoc2 = ELEMENT_NOT_FOUND 
      END IF ! pBorder%nVertShared

! ------------------------------------------------------------------------------
!     Vertex not in list of vertices to be shared, so add
! ------------------------------------------------------------------------------

      IF ( iLoc2 == ELEMENT_NOT_FOUND ) THEN 
        pBorder%nVertShared = pBorder%nVertShared + 1
 
        IF ( pPatch%bcType == BC_SYMMETRY ) THEN
          pBorder%ivgShared(pBorder%nVertShared) = ivg 
          ivgSharedSorted(pBorder%nVertShared)   = ivg 
        ELSE IF ( pPatch%bcType == BC_PERIODIC ) THEN
          pBorder%ivgShared(pBorder%nVertShared) = ivgm 
          ivgSharedSorted(pBorder%nVertShared)   = ivgm
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! pPatch%bcType     
        
        IF ( pBorder%nVertShared > 1 ) THEN 
          CALL QuickSortInteger(ivgSharedSorted(1:pBorder%nVertShared), & 
                                pBorder%nVertShared)                
        END IF ! pBorder%nVertShared
      END IF ! iLoc2       
    END IF ! iLoc
  END DO ! ivl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SYPE_AugmentVertexLists 
  
  
  
  
  
 
 
! ******************************************************************************
!
! Purpose: Build list of virtual cells adjacent to periodic and symmetry patches
!   on a parallel region.
!
! Description: None.
!
! Input:
!   pRegion           Pointer to parallel region
!   pRegionSerial     Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_SYPE_BuildP2VCList(pRegion,pRegionSerial)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ==============================================================================
!   Locals
! ==============================================================================
      
    INTEGER :: errorFlag,iBorder,icg,icgs,icl,iLoc,iPatch,iPatchSerial
    INTEGER, DIMENSION(:), ALLOCATABLE :: iPatchGlob2iPatchLoc
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
    TYPE(t_patch), POINTER :: pPatch,pPatchSerial
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SYPE_BuildP2VCList',&
  'RFLU_ModSymmetryPeriodic.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building patch-to-virtual-cell list...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal                               
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Check and set grid pointer
! ******************************************************************************

    IF ( pRegion%iRegionGlobal == 0 ) THEN ! Defensive coding
      CALL ErrorStop(global,ERR_REGION_ID_INVALID,__LINE__)
    END IF ! pRegion%iRegionGlobal

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid    
                          
! ******************************************************************************
!   Allocate temporary memory for and build mapping of global to local patch 
!   indices
! ******************************************************************************

    ALLOCATE(iPatchGlob2iPatchLoc(pGridSerial%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'iPatchGlob2iPatchLoc')
    END IF ! global%error
    
    DO iPatchSerial = 1,pGridSerial%nPatches
      iPatchGlob2iPatchLoc(iPatchSerial) = CRAZY_VALUE_INT
    END DO ! iPatchSerial    
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      iPatchGlob2iPatchLoc(pPatch%iPatchGlobal) = iPatch                      
    END DO ! iPatch
    
! ******************************************************************************
!   Initialize number of virtual cells for each patch. NOTE should not be 
!   necessary, but done for safety
! ******************************************************************************
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      pPatch%nBCellsVirt = 0      
    END DO ! iPatch    
    
! ******************************************************************************
!   Determine number of virtual cells for periodic and symmetry patches
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders 
      pBorder => pGrid%borders(iBorder)

      DO icl = 1,pBorder%nCellsRecv      
        icg  = pBorder%icgRecv(icl)                
        icgs = pGrid%pc2sc(icg)

        patchSerialLoop: DO iPatchSerial = 1,pGridSerial%nPatches
          pPatchSerial => pRegionSerial%patches(iPatchSerial)

          IF ( pPatchSerial%bcType == BC_PERIODIC .OR. & 
               pPatchSerial%bcType == BC_SYMMETRY ) THEN
            CALL BinarySearchInteger(pPatchSerial%bvc, &
                                     pPatchSerial%nBCellsVirt,icgs,iLoc) 

            IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
              iPatch =  iPatchGlob2iPatchLoc(iPatchSerial)              
              pPatch => pRegion%patches(iPatch)

              pPatch%nBCellsVirt = pPatch%nBCellsVirt + 1
              
              EXIT patchSerialLoop      
            END IF ! iLoc                                  
          END IF ! pPatchSerial%bcType
        END DO patchSerialLoop        
      END DO ! icl
    END DO ! iBorder

! ******************************************************************************
!   Allocate memory and set counters back to zero
! ******************************************************************************
        
    CALL RFLU_SYPE_CreateP2VCList(pRegion)
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
            
      pPatch%nBCellsVirt = 0 
    END DO ! iPatch    
        
! ******************************************************************************
!   Find and store virtual cells for periodic and symmetry patches
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)    
                  
      DO icl = 1,pBorder%nCellsRecv
        icg  = pBorder%icgRecv(icl)
        icgs = pGrid%pc2sc(icg)
              
        patchSerialLoop2: DO iPatchSerial = 1,pGridSerial%nPatches
          pPatchSerial => pRegionSerial%patches(iPatchSerial)
          
          IF ( pPatchSerial%bcType == BC_PERIODIC .OR. & 
               pPatchSerial%bcType == BC_SYMMETRY ) THEN
            CALL BinarySearchInteger(pPatchSerial%bvc, &
                                     pPatchSerial%nBCellsVirt,icgs,iLoc) 
                                     
            IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
              iPatch =  iPatchGlob2iPatchLoc(iPatchSerial)              
              pPatch => pRegion%patches(iPatch)
              
              pPatch%nBCellsVirt = pPatch%nBCellsVirt + 1              
              pPatch%bvc(pPatch%nBCellsVirt) = icg

              EXIT patchSerialLoop2      
            END IF ! iLoc                                  
          END IF ! pPatchSerial%bcType
        END DO patchSerialLoop2
      END DO ! icl
    END DO ! iBorder        
               
! ******************************************************************************
!   Deallocate temporary memory 
! ******************************************************************************

    DEALLOCATE(iPatchGlob2iPatchLoc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'iPatchGlob2iPatchLoc')
    END IF ! global%error        
        
! ******************************************************************************
!   End 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building patch-to-virtual-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_SYPE_BuildP2VCList  
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Build list of virtual cells adjacent to periodic and symmetry patches
!   on serial region.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine can only be used on the serial region.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_SYPE_BuildP2VCListSerial(pRegion)

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
      
    INTEGER :: errorFlag,iBorder,icl,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SYPE_BuildP2VCListSerial',&
  'RFLU_ModSymmetryPeriodic.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building patch-to-virtual-cell list...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal                               
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Check and set grid pointer
! ******************************************************************************

    IF ( pRegion%iRegionGlobal /= 0 ) THEN ! Defensive coding
      CALL ErrorStop(global,ERR_REGION_ID_INVALID,__LINE__)
    END IF ! pRegion%iRegionGlobal

    pGrid => pRegion%grid
                          
! ******************************************************************************
!   Build list
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%bcType == BC_PERIODIC .OR. & 
           pPatch%bcType == BC_SYMMETRY ) THEN
        iBorder =  pPatch%iBorder
        pBorder => pGrid%borders(iBorder)   
                         
        pPatch%nBCellsVirt = pBorder%nCellsRecv                    
              
        ALLOCATE(pPatch%bvc(pPatch%nBCellsVirt),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvc')
        END IF ! global%error

        DO icl = 1,pPatch%nBCellsVirt
          pPatch%bvc(icl) = pBorder%icgRecv(icl)
        END DO ! icl
      ELSE 
        pPatch%nBCellsVirt = 0
        
        NULLIFY(pPatch%bvc)
      END IF ! pPatch%bcType
    END DO ! iPatch
        
! ******************************************************************************
!   End 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building patch-to-virtual-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_SYPE_BuildP2VCListSerial  
 
 



! ******************************************************************************
!
! Purpose: Build transforms between periodic patches.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SYPE_BuildTransforms(pRegion)

  USE ModTools

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

  INTEGER :: iPatch,iPatchRelated
  REAL(RFREAL) :: ct,dotp,eqTol,ex,ey,ez,nx,nxr,ny,nyr,nz,nzr,st,theta
  TYPE(t_border), POINTER :: pBorder
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch,pPatchRelated

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_BuildTransforms',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building periodicity transforms...' 
  END IF ! global%myProcid

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
 
  pGrid => pRegion%grid
 
  eqTol = 1.0E-6_RFREAL
 
! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
       
! ==============================================================================       
!   Periodic patch, find related patch       
! ==============================================================================       

    IF ( pPatch%bcType == BC_PERIODIC ) THEN 
      iPatchRelated = pPatch%iPatchRelated
      
      pPatchRelated => pRegion%patches(iPatchRelated)
      
! ------------------------------------------------------------------------------
!     Check whether the two patches have the same flatness flag
! ------------------------------------------------------------------------------      
      
      IF ( pPatch%flatFlag .EQV. pPatchRelated%flatFlag ) THEN
        nx = pPatch%pn(XCOORD)
        ny = pPatch%pn(YCOORD)
        nz = pPatch%pn(ZCOORD)
        
        nxr = pPatchRelated%pn(XCOORD)
        nyr = pPatchRelated%pn(YCOORD)
        nzr = pPatchRelated%pn(ZCOORD)        
       
        dotp = nx*nxr + ny*nyr + nz*nzr
       
! ----- Linear periodicity -----------------------------------------------------       
       
        IF ( (FloatEqual(dotp,-1.0_RFREAL,eqTol) .EQV. .TRUE.) ) THEN
          pPatch%tm(XCOORD,XCOORD) = 1.0_RFREAL
          pPatch%tm(YCOORD,XCOORD) = 0.0_RFREAL
          pPatch%tm(ZCOORD,XCOORD) = 0.0_RFREAL    
          pPatch%tm(XYZMAG,XCOORD) = 0.0_RFREAL 
          
          pPatch%tm(XCOORD,YCOORD) = 0.0_RFREAL
          pPatch%tm(YCOORD,YCOORD) = 1.0_RFREAL
          pPatch%tm(ZCOORD,YCOORD) = 0.0_RFREAL    
          pPatch%tm(XYZMAG,YCOORD) = 0.0_RFREAL              
              
          pPatch%tm(XCOORD,ZCOORD) = 0.0_RFREAL
          pPatch%tm(YCOORD,ZCOORD) = 0.0_RFREAL
          pPatch%tm(ZCOORD,ZCOORD) = 1.0_RFREAL    
          pPatch%tm(XYZMAG,ZCOORD) = 0.0_RFREAL              
              
          pPatch%tm(XCOORD,XYZMAG) = pPatchRelated%pc(XCOORD) - pPatch%pc(XCOORD)
          pPatch%tm(YCOORD,XYZMAG) = pPatchRelated%pc(YCOORD) - pPatch%pc(YCOORD)
          pPatch%tm(ZCOORD,XYZMAG) = pPatchRelated%pc(ZCOORD) - pPatch%pc(ZCOORD)    
          pPatch%tm(XYZMAG,XYZMAG) = 1.0_RFREAL              
              
! ----- Rotational periodicity -------------------------------------------------              
              
        ELSE 
          theta = pPatch%angleRelated

          ct = COS(theta)
          st = SIN(theta)

          SELECT CASE ( pPatch%axisRelated ) 
            CASE ( 1 ) 
              ex = 1.0_RFREAL
              ey = 0.0_RFREAL
              ez = 0.0_RFREAL
            CASE ( 2 ) 
              ex = 0.0_RFREAL
              ey = 1.0_RFREAL
              ez = 0.0_RFREAL
            CASE ( 3 ) 
              ex = 0.0_RFREAL
              ey = 0.0_RFREAL
              ez = 1.0_RFREAL
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%axisRelated

          pPatch%tm(XCOORD,XCOORD) = ct + (1.0_RFREAL-ct)*ex*ex
          pPatch%tm(YCOORD,XCOORD) =      (1.0_RFREAL-ct)*ex*ey - st*ez
          pPatch%tm(ZCOORD,XCOORD) =      (1.0_RFREAL-ct)*ex*ez + st*ey
          pPatch%tm(XYZMAG,XCOORD) = 0.0_RFREAL 

          pPatch%tm(XCOORD,YCOORD) =      (1.0_RFREAL-ct)*ey*ex + st*ez
          pPatch%tm(YCOORD,YCOORD) = ct + (1.0_RFREAL-ct)*ey*ey
          pPatch%tm(ZCOORD,YCOORD) =      (1.0_RFREAL-ct)*ey*ez - st*ex
          pPatch%tm(XYZMAG,YCOORD) = 0.0_RFREAL 

          pPatch%tm(XCOORD,ZCOORD) =      (1.0_RFREAL-ct)*ez*ex - st*ey
          pPatch%tm(YCOORD,ZCOORD) =      (1.0_RFREAL-ct)*ez*ey + st*ex
          pPatch%tm(ZCOORD,ZCOORD) = ct + (1.0_RFREAL-ct)*ez*ez  
          pPatch%tm(XYZMAG,ZCOORD) = 0.0_RFREAL 
          
          pPatch%tm(XCOORD,XYZMAG) = 0.0_RFREAL 
          pPatch%tm(YCOORD,XYZMAG) = 0.0_RFREAL 
          pPatch%tm(ZCOORD,XYZMAG) = 0.0_RFREAL  
          pPatch%tm(XYZMAG,XYZMAG) = 1.0_RFREAL                        
        END IF ! FloatEqual  
        
! ------------------------------------------------------------------------------
!     Invalid combination - related patches must have same flatness flag
! ------------------------------------------------------------------------------        
        
      ELSE 
        CALL ErrorStop(global,ERR_FLATFLAG_INCONSISTENT,__LINE__)
      END IF ! pPatch%flatFlag
    END IF ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building periodicity transforms done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SYPE_BuildTransforms
 









! ******************************************************************************
!
! Purpose: Build vertex maps which link local vertex lists on related periodic 
!   patches.
!
! Description: Coordinates on current patch are transformed to related patch, 
!   then use Octree to accelerate search for matching vertex pairs.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SYPE_BuildVertexMaps(pRegion)

  USE ModTools
  
  USE RFLU_ModOctree

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

  INTEGER :: cvSize,errorFlag,icv,iPatch,iPatchRelated,ivg,ivg2,ivl,ivlMin,ivl2
  INTEGER, DIMENSION(:), ALLOCATABLE :: cv
  REAL(RFREAL) :: delFrac,dist,distMax,distMin,distMinSave,eqTol,maxFrac, &
                  xDel,xMax,xMin,yDel,yMax,yMin,zDel,zMax,zMin
  REAL(RFREAL) :: xyz(XCOORD:XYZMAG),xyzt(XCOORD:XYZMAG)
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: xyzr
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch,pPatchRelated

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_BuildVertexMaps',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building vertex maps...' 
  END IF ! global%myProcid

! ******************************************************************************
! Set pointers and constants
! ******************************************************************************
 
  pGrid => pRegion%grid
 
  cvSize = 10
  
  delFrac = 0.01_RFREAL  
  maxFrac = 0.1_RFREAL

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************
  
  ALLOCATE(cv(cvSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cv')
  END IF ! global%error 
 
! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

! ==============================================================================       
!   Set tolerance for vertex matching
! ==============================================================================       
       
       
! ==============================================================================       
!   Periodic patch, find related patch       
! ==============================================================================       

    IF ( pPatch%bcType == BC_PERIODIC ) THEN     
      iPatchRelated =  pPatch%iPatchRelated            
      pPatchRelated => pRegion%patches(iPatchRelated)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN   
        WRITE(STDOUT,'(A,2X,2(1X,A,1X,I2),A)') SOLVER_NAME,'Matching patches', &
                                               iPatch,'and',iPatchRelated,'...'
      END IF ! global%myProcid      

! ------------------------------------------------------------------------------
!     Index of related patch greater than that of current patch, so match 
!     vertices
! ------------------------------------------------------------------------------
                  
      IF ( iPatchRelated > iPatch ) THEN 

! ----- Initialize and write info  ---------------------------------------------

        eqTol = 0.1_RFREAL*SQRT(MINVAL(pPatch%fn(XYZMAG,1:pPatch%nBFaces)))
      
        distMax     = -HUGE(1.0_RFREAL)
        distMinSave =  HUGE(1.0_RFREAL)            
      
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN           
          WRITE(STDOUT,'(A,5X,A,24X,E13.6)') SOLVER_NAME, &
                                             'Distance tolerance:',eqTol      
        END IF ! global%myProcid      
                  
! ----- Allocate temporary memory ----------------------------------------------      
      
        ALLOCATE(xyzr(XCOORD:ZCOORD,pPatchRelated%nBVert),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'xyzr')
        END IF ! global%error
        
! ----- Set coordinates of vertices on related patch ---------------------------        
        
        DO ivl = 1,pPatchRelated%nBVert
          ivg = pPatchRelated%bv(ivl)
          
          xyzr(XCOORD,ivl) = pGrid%xyz(XCOORD,ivg)
          xyzr(YCOORD,ivl) = pGrid%xyz(YCOORD,ivg)
          xyzr(ZCOORD,ivl) = pGrid%xyz(ZCOORD,ivg) 
        END DO ! ivl
        
! ----- Build Octree based on vertices on related patch ------------------------        
        
        xMin = MINVAL(xyzr(XCOORD,1:pPatchRelated%nBVert))
        xMax = MAXVAL(xyzr(XCOORD,1:pPatchRelated%nBVert))  
        yMin = MINVAL(xyzr(YCOORD,1:pPatchRelated%nBVert))
        yMax = MAXVAL(xyzr(YCOORD,1:pPatchRelated%nBVert))                 
        zMin = MINVAL(xyzr(ZCOORD,1:pPatchRelated%nBVert))
        zMax = MAXVAL(xyzr(ZCOORD,1:pPatchRelated%nBVert))           
               
        xDel = xMax - xMin 
        yDel = yMax - yMin
        zDel = zMax - zMin
        
        xDel = MAX(xDel,MIN(maxFrac*yDel,maxFrac*zDel))
        yDel = MAX(yDel,MIN(maxFrac*xDel,maxFrac*zDel))
        zDel = MAX(zDel,MIN(maxFrac*xDel,maxFrac*yDel))                

        xMin = xMin - delFrac*xDel
        xMax = xMax + delFrac*xDel 
        yMin = yMin - delFrac*yDel
        yMax = yMax + delFrac*yDel 
        zMin = zMin - delFrac*zDel
        zMax = zMax + delFrac*zDel  
        
        CALL RFLU_CreateOctree(global,pPatchRelated%nBVert)
        CALL RFLU_BuildOctree(xyzr(XCOORD,1:pPatchRelated%nBVert), & 
                              xyzr(YCOORD,1:pPatchRelated%nBVert), &
                              xyzr(ZCOORD,1:pPatchRelated%nBVert), &
                              xMin,xMax,yMin,yMax,zMin,zMax)
                              
! ----- Loop over vertices on current patch ------------------------------------                              
                              
        DO ivl = 1,pPatch%nBVert
          ivg = pPatch%bv(ivl)
                    
! ------- Transform coordinates of current vertex          
          
          xyz(XCOORD) = pGrid%xyz(XCOORD,ivg)
          xyz(YCOORD) = pGrid%xyz(YCOORD,ivg)
          xyz(ZCOORD) = pGrid%xyz(ZCOORD,ivg)
          xyz(XYZMAG) = 1.0_RFREAL
                              
          xyzt = MATMUL(pPatch%tm,xyz)
          
! ------- Find closest <cvSize> vertices on related patch           
          
          CALL RFLU_QueryOctree(xyzt(XCOORD),xyzt(YCOORD),xyzt(ZCOORD), & 
                                cvSize,cv)
                                
! ------- Find closest vertex                                
                                
          distMin = HUGE(1.0_RFREAL)                      
                                
          cvLoop: DO icv = 1,cvSize
            ivl2 = cv(icv)
            
            ivg2 = pPatchRelated%bv(ivl2)
    
            dist = (xyzt(XCOORD)-pGrid%xyz(XCOORD,ivg2))**2 & 
                 + (xyzt(YCOORD)-pGrid%xyz(YCOORD,ivg2))**2 &
                 + (xyzt(ZCOORD)-pGrid%xyz(ZCOORD,ivg2))**2
       
            IF ( dist > distMax ) THEN 
              distMax = dist
            END IF ! dist 
              
            IF ( dist <= distMin ) THEN 
              distMin = dist
              ivlMin  = ivl2
              
              IF ( distMin <= eqTol ) THEN 
                EXIT cvLoop
              END IF ! distMin           
            END IF ! dist 
          END DO cvLoop  
                          
! ------- Store closest vertex in mapping array                           
                                                    
          IF ( distMin <= eqTol ) THEN 
            pPatch%bv2bv(ivl) = ivlMin            
          ELSE 
            CALL ErrorStop(global,ERR_VERTEX_MATCH_FAILED,__LINE__)
          END IF ! distMin
          
          distMinSave = MIN(distMin,distMinSave)                          
        END DO ! ivl                      
                     
! ----- Destroy Octree and deallocate temporary memory -------------------------                     
                     
        CALL RFLU_DestroyOctree(global)             
                              
        DEALLOCATE(xyzr,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'xyzr')
        END IF ! global%error        

! ----- Write some info --------------------------------------------------------

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN   
          WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                'Maximum distance between matched vertices:',distMax
          WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                'Minimum distance between matched vertices:',distMinSave                
        END IF ! global%myProcid   

! ------------------------------------------------------------------------------
!     Index of related patch less than that of current patch, so copy previously
!     built vertex matching array
! ------------------------------------------------------------------------------
         
      ELSE            
        DO ivl = 1,pPatch%nBVert
          ivl2 = pPatchRelated%bv2bv(ivl)

          pPatch%bv2bv(ivl2) = ivl
        END DO ! ivl
      END IF ! iPatchRelated
         
    END IF ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************
  
  DEALLOCATE(cv,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cv')
  END IF ! global%error 

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building vertex maps done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SYPE_BuildVertexMaps






! ******************************************************************************
!
! Purpose: Create list of virtual cells adjacent to periodic and symmetry 
!   patches.
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
    
  SUBROUTINE RFLU_SYPE_CreateP2VCList(pRegion)

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
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SYPE_CreateP2VCList',&
  'RFLU_ModSymmetryPeriodic.F90')
    
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating patch-to-virtual-cell list...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal                                         
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
                          
! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( pPatch%nBCellsVirt > 0 ) THEN
        ALLOCATE(pPatch%bvc(pPatch%nBCellsVirt),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvc')
        END IF ! global%error
      ELSE 
        NULLIFY(pPatch%bvc)        
      END IF ! pPatch%nBCellsVirt
    END DO ! iPatch
        
! ******************************************************************************
!   End 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating patch-to-virtual-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_SYPE_CreateP2VCList   
  
  
   
  


! ******************************************************************************
!
! Purpose: Create vertex maps.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SYPE_CreateVertexMaps(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_CreateVertexMaps',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating vertex maps...' 
  END IF ! global%myProcid
 
  pGrid => pRegion%grid
 
! ******************************************************************************
! Loop over patches and allocate memory
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
       
    IF ( pPatch%bcType == BC_PERIODIC ) THEN 
      ALLOCATE(pPatch%bv2bv(pPatch%nBVert),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bv2bv')
      END IF ! global%error
    ELSE 
      NULLIFY(pPatch%bv2bv)
    END IF ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating vertex maps done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SYPE_CreateVertexMaps
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Destroy list of virtual cells adjacent to periodic and symmetry 
!   patches.
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
    
  SUBROUTINE RFLU_SYPE_DestroyP2VCList(pRegion)

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
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SYPE_DestroyP2VCList',&
  'RFLU_ModSymmetryPeriodic.F90')
    
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying patch-to-virtual-cell list...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
                          
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( pPatch%nBCellsVirt > 0 ) THEN 
        pPatch%nBCellsVirt = 0                                                            
                                                                          
        DEALLOCATE(pPatch%bvc,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bvc')
        END IF ! global%error
      ELSE       
        NULLIFY(pPatch%bvc)  
      END IF ! pPatch%nBCellsVirt
    END DO ! iPatch
        
! ******************************************************************************
!   End 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying patch-to-virtual-cell list done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_SYPE_DestroyP2VCList   
    
  
  
  
   
 
 
 
! ******************************************************************************
!
! Purpose: For given virtual cell index on serial region, find corresponding 
!   actual cell index on serial region.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   icgs        Cell index (must be virtual cell index) 
!
! Output: 
!   icgs2       Cell index (actual cell index)  
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_SYPE_GetActualSerialCell(pRegionSerial,icgs,icgs2)          

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icgs
  INTEGER, INTENT(OUT) :: icgs2
  TYPE(t_region), POINTER :: pRegionSerial

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: foundFlag
  INTEGER :: iBorderSerial,iBorderSerial2,iLoc
  TYPE(t_border), POINTER :: pBorderSerial,pBorderSerial2
  TYPE(t_global), POINTER :: global   
  TYPE(t_grid), POINTER :: pGridSerial   

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegionSerial%global

  CALL RegisterFunction(global,'RFLU_SYPE_GetActualSerialCell',&
  'RFLU_ModSymmetryPeriodic.F90')

! ******************************************************************************
! Set pointers 
! ******************************************************************************

  pGridSerial => pRegionSerial%grid

! ******************************************************************************
! Check region index and cell kind
! ******************************************************************************

  IF ( pRegionSerial%iRegionGlobal /= 0 ) THEN 
    CALL ErrorStop(global,ERR_REGION_ID_INVALID,__LINE__)
  END IF ! pRegionSerial%iRegionGlobal

  IF ( (icgs <= pGridSerial%nCells   ) .OR. & 
       (icgs >  pGridSerial%nCellsTot) ) THEN 
    CALL ErrorStop(global,ERR_CELL_KIND_INVALID,__LINE__)
  END IF ! icgs

! ******************************************************************************
! Loop over serial borders and search for cell 
! ******************************************************************************

  borderLoop2: DO iBorderSerial = 1,pGridSerial%nBorders
    pBorderSerial => pGridSerial%borders(iBorderSerial)

    CALL BinarySearchInteger(pBorderSerial%icgRecv(1:pBorderSerial%nCellsRecv), &
                             pBorderSerial%nCellsRecv,icgs,iLoc)

    IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
      foundFlag = .TRUE.

      iBorderSerial2 = pBorderSerial%iBorder

      IF ( iBorderSerial2 > pGridSerial%nBorders ) THEN 
        CALL ErrorStop(global,ERR_BORDER_INDEX_INVALID,__LINE__)            
      END IF ! iBorderSerial2

      pBorderSerial2 => pGridSerial%borders(iBorderSerial2)

      icgs2 = pBorderSerial2%icgSend(iLoc)

      IF ( icgs2 > pGridSerial%nCells ) THEN 
        CALL ErrorStop(global,ERR_CELL_KIND_INVALID,__LINE__)      
      END IF ! icgs2

      EXIT borderLoop2
    END IF ! iLoc                                    
  END DO borderLoop2

  IF ( foundFlag .EQV. .FALSE. ) THEN 
    CALL ErrorStop(global,ERR_CELL_NOT_FOUND,__LINE__)          
  END IF ! foundFlag
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SYPE_GetActualSerialCell 
 
 
  
 
 




! ******************************************************************************
!
! Purpose: Given serial vertex index, find index of related vertex in another
!   region. 
!
! Description: None.
!
! Input:
!   pRegionSerial       Pointer to serial region data
!   pRegion             Pointer to region data
!   ivgs                Serial vertex index
!
! Output: 
!   ivg                 Related vertex index
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_SYPE_GetRelatedVertex(pRegionSerial,pRegion,ivgs,ivg)          

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: ivgs
  INTEGER, INTENT(OUT) :: ivg
  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: exitBorderLoopFlag,foundFlagRecv,foundFlagShared
  INTEGER :: errorFlag,iBorderSerial,iBorderSerial2,iLoc,ivgs2,ivl,ivls2
  INTEGER, DIMENSION(:), ALLOCATABLE :: ivgSharedSorted,sortKey 
  TYPE(t_border), POINTER :: pBorderSerial,pBorderSerial2
  TYPE(t_global), POINTER :: global   
  TYPE(t_grid), POINTER :: pGrid,pGridSerial   

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegionSerial%global

  CALL RegisterFunction(global,'RFLU_SYPE_GetRelatedVertex',&
  'RFLU_ModSymmetryPeriodic.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGridSerial => pRegionSerial%grid
  pGrid       => pRegion%grid

  foundFlagRecv   = .FALSE.
  foundFlagShared = .FALSE.

! ******************************************************************************
! Check region index
! ******************************************************************************

  IF ( pRegionSerial%iRegionGlobal /= 0 ) THEN 
    CALL ErrorStop(global,ERR_REGION_ID_INVALID,__LINE__)
  END IF ! pRegionSerial%iRegionGlobal
  
! ******************************************************************************
! Search for vertex in vertices to be shared
! ******************************************************************************

  exitBorderLoopFlag = .FALSE.
 
  borderLoopShared: DO iBorderSerial = 1,pGridSerial%nBorders
    pBorderSerial => pGridSerial%borders(iBorderSerial)

    ALLOCATE(ivgSharedSorted(pBorderSerial%nVertShared),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ivgSharedSorted')
    END IF ! global%error

    ALLOCATE(sortKey(pBorderSerial%nVertShared),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'sortKey')
    END IF ! global%error                  

    DO ivl = 1,pBorderSerial%nVertShared
      ivgSharedSorted(ivl) = pBorderSerial%ivgShared(ivl)
      sortKey(ivl)         = ivl                  
    END DO ! ivl

    CALL QuickSortIntegerInteger(ivgSharedSorted,sortKey, & 
                                 pBorderSerial%nVertShared)

    CALL BinarySearchInteger(ivgSharedSorted,pBorderSerial%nVertShared,ivgs, &
                             iLoc)

    IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
      iBorderSerial2 = pBorderSerial%iBorder

      pBorderSerial2 => pGridSerial%borders(iBorderSerial2)

      ivls2 = sortKey(iLoc)
      ivgs2 = pBorderSerial2%ivgShared(ivls2)

      CALL BinarySearchInteger(pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                               pGrid%nVertTot,ivgs2,iLoc)

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN
        ivg = pGrid%sv2pv(2,iLoc)

        IF ( ivg > pGrid%nVert ) THEN ! Not actual vertex
          CALL ErrorStop(global,ERR_VERTEX_KIND_INVALID,__LINE__)               
        END IF ! ivg
                         
        exitBorderLoopFlag = .TRUE.                            
        foundFlagShared    = .TRUE.                                                               
      ELSE
        CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)             
      END IF ! iLoc
    END IF ! iLoc                                           

    DEALLOCATE(ivgSharedSorted,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ivgSharedSorted')
    END IF ! global%error 

    DEALLOCATE(sortKey,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'sortKey')
    END IF ! global%error                                

    IF ( exitBorderLoopFlag .EQV. .TRUE. ) THEN 
      EXIT borderLoopShared
    END IF ! exitBorderLoopFlag                                                      
  END DO borderLoopShared
  
! ******************************************************************************
! Search for vertex in vertices to be received
! ******************************************************************************
    
  IF ( foundFlagShared .EQV. .FALSE. ) THEN 
    exitBorderLoopFlag = .FALSE.

    borderLoopRecv: DO iBorderSerial = 1,pGridSerial%nBorders
      pBorderSerial => pGridSerial%borders(iBorderSerial)
                
      CALL BinarySearchInteger(pBorderSerial%ivgRecv(1:pBorderSerial%nVertRecv), &
                               pBorderSerial%nVertRecv,ivgs,iLoc)

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        iBorderSerial2 = pBorderSerial%iBorder

        pBorderSerial2 => pGridSerial%borders(iBorderSerial2)

        ivgs2 = pBorderSerial2%ivgSend(iLoc)

        CALL BinarySearchInteger(pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                                 pGrid%nVertTot,ivgs2,iLoc)

        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN
          ivg = pGrid%sv2pv(2,iLoc)

          IF ( ivg > pGrid%nVert ) THEN ! Not actual vertex
            CALL ErrorStop(global,ERR_VERTEX_KIND_INVALID,__LINE__)               
          END IF ! ivg

          exitBorderLoopFlag = .TRUE.                            
          foundFlagRecv      = .TRUE.                                                               
        ELSE
          CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)              
        END IF ! iLoc
      END IF ! iLoc                                           
                   
      IF ( exitBorderLoopFlag .EQV. .TRUE. ) THEN 
        EXIT borderLoopRecv
      END IF ! exitBorderLoopFlag                                                      
    END DO borderLoopRecv
  END IF ! foundFlagShared
  
  
  IF ( (foundFlagShared .EQV. .FALSE.) .AND. & 
       (foundFlagRecv   .EQV. .FALSE.) ) THEN 
    CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)
  END IF ! foundFlagShared  

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SYPE_GetRelatedVertex

 
 
 
 
 
 
 
 
 
! ******************************************************************************
!
! Purpose: Determine whether have symmetry and/or periodic patches.
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

LOGICAL FUNCTION RFLU_SYPE_HaveSyPePatches(pRegion)                                   

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

  INTEGER :: errorFlag,iPatch,nSyPePatches
  TYPE(t_global), POINTER :: global   
  TYPE(t_grid), POINTER :: pGrid      
  TYPE(t_patch), POINTER :: pPatch    

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_HaveSyPePatches',&
  'RFLU_ModSymmetryPeriodic.F90')

! ******************************************************************************
! Set pointers 
! ******************************************************************************

  pGrid => pRegion%grid

  nSyPePatches = 0

! ******************************************************************************
! Determine whether have symmetry or periodic patches
! ******************************************************************************

  iPatchLoop: DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    IF ( pPatch%bcType == BC_SYMMETRY .OR. pPatch%bcType == BC_PERIODIC ) THEN 
      nSyPePatches = nSyPePatches + 1
      
      EXIT iPatchLoop
    END IF ! pPatch%bcType
  END DO iPatchLoop

  IF ( nSyPePatches > 0 ) THEN 
    RFLU_SYPE_HaveSyPePatches = .TRUE.
  ELSE
    RFLU_SYPE_HaveSyPePatches = .FALSE.
  END IF ! nPatches

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_SYPE_HaveSyPePatches
 
 







! ******************************************************************************
!
! Purpose: Impose vertex maps.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch
!   ivg         Vertex index
!
! Output: 
!   ivgm        Vertex index
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SYPE_ImposeVertexMap(pRegion,pPatch,ivg,ivgm)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: ivg
  INTEGER, INTENT(OUT) :: ivgm
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion 

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,iLoc,iPatchRelated,ivlm
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatchRelated  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_ImposeVertexMap',&
  'RFLU_ModSymmetryPeriodic.F90')
 
! ******************************************************************************
! Get related patch
! ******************************************************************************
 
  iPatchRelated =  pPatch%iPatchRelated
  pPatchRelated => pRegion%patches(iPatchRelated)

! ******************************************************************************
! Search for vertex and impose mapping
! ******************************************************************************

  CALL BinarySearchInteger(pPatch%bv(1:pPatch%nBVert),pPatch%nBVert,ivg,iLoc)
                                   
  IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
    ivlm = pPatch%bv2bv(iLoc)
    ivgm = pPatchRelated%bv(ivlm)
  ELSE 
    CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)       
  END IF ! iLoc
    
! ******************************************************************************
! End
! ******************************************************************************
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SYPE_ImposeVertexMap
  





! ******************************************************************************
!
! Purpose: Read and impose transforms between periodic patches.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SYPE_ReadTransforms(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNamePlain  
 
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

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: errorFlag,iFile,iPatch,iPatchGlobal
  REAL(RFREAL) :: tm(XCOORD:XYZMAG,XCOORD:XYZMAG)
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_ReadTransforms',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading periodicity transforms...'
      
      IF ( global%verbLevel > VERBOSE_LOW ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel
    END IF ! global%verbLevel 
  END IF ! global%myProcid

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
 
  pGrid => pRegion%grid

! ******************************************************************************
! Open file
! ******************************************************************************

  iFile = IF_PTMATRIX

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.ptm',iFileName)

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
  global%error = errorFlag        
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Read header
! ******************************************************************************

  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= & 
       '# ROCFLU periodic-transformation matrix file' ) THEN
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM
     
! ******************************************************************************
! Read remainder of file
! ******************************************************************************

  emptyLoop: DO 
    READ(iFile,'(A)') sectionString

    SELECT CASE ( TRIM(sectionString) ) 
    
! ==============================================================================
!     Patch section
! ==============================================================================
    
      CASE ( '# Patch' ) 
        READ(iFile,'(I3)') iPatchGlobal
      
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME,'Global patch:', &
                                         iPatchGlobal
        END IF ! global%myProcid      
      
        READ(iFile,'(4(E23.16))') tm(XCOORD,XCOORD),tm(XCOORD,YCOORD), & 
                                  tm(XCOORD,ZCOORD),tm(XCOORD,XYZMAG)
        READ(iFile,'(4(E23.16))') tm(YCOORD,XCOORD),tm(YCOORD,YCOORD), & 
                                  tm(YCOORD,ZCOORD),tm(YCOORD,XYZMAG)
        READ(iFile,'(4(E23.16))') tm(ZCOORD,XCOORD),tm(ZCOORD,YCOORD), & 
                                  tm(ZCOORD,ZCOORD),tm(ZCOORD,XYZMAG)
        READ(iFile,'(4(E23.16))') tm(XYZMAG,XCOORD),tm(XYZMAG,YCOORD), & 
                                  tm(XYZMAG,ZCOORD),tm(XYZMAG,XYZMAG)        
              
! ------------------------------------------------------------------------------
!       Loop over local patches and store transformation matrix
! ------------------------------------------------------------------------------      
            
        patchLoop: DO iPatch = 1,pGrid%nPatches
          pPatch => pRegion%patches(iPatch)
          
          IF ( pPatch%iPatchGlobal == iPatchGlobal ) THEN
            IF ( global%myProcid == MASTERPROC .AND. &
                 global%verbLevel > VERBOSE_LOW ) THEN
              WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Local patch:', &
                                             iPatch,'.'
            END IF ! global%myProcid
           
            pPatch%transformFlag = .TRUE. ! Set flag for later check
           
            pPatch%tm(XCOORD,XCOORD) = tm(XCOORD,XCOORD)
            pPatch%tm(XCOORD,YCOORD) = tm(XCOORD,YCOORD)
            pPatch%tm(XCOORD,ZCOORD) = tm(XCOORD,ZCOORD)
            pPatch%tm(XCOORD,XYZMAG) = tm(XCOORD,XYZMAG)
            pPatch%tm(YCOORD,XCOORD) = tm(YCOORD,XCOORD)
            pPatch%tm(YCOORD,YCOORD) = tm(YCOORD,YCOORD)
            pPatch%tm(YCOORD,ZCOORD) = tm(YCOORD,ZCOORD)
            pPatch%tm(YCOORD,XYZMAG) = tm(YCOORD,XYZMAG)
            pPatch%tm(ZCOORD,XCOORD) = tm(ZCOORD,XCOORD)
            pPatch%tm(ZCOORD,YCOORD) = tm(ZCOORD,YCOORD)
            pPatch%tm(ZCOORD,ZCOORD) = tm(ZCOORD,ZCOORD)
            pPatch%tm(ZCOORD,XYZMAG) = tm(ZCOORD,XYZMAG)
            pPatch%tm(XYZMAG,XCOORD) = tm(XYZMAG,XCOORD)
            pPatch%tm(XYZMAG,YCOORD) = tm(XYZMAG,YCOORD)
            pPatch%tm(XYZMAG,ZCOORD) = tm(XYZMAG,ZCOORD)
            pPatch%tm(XYZMAG,XYZMAG) = tm(XYZMAG,XYZMAG)
                                                                
            EXIT patchLoop           
          END IF ! pPatch%iPatchGlobal 
        END DO patchLoop
    
! ==============================================================================
!     Invalid section string
! ==============================================================================

      CASE ( '# End' )
        EXIT emptyLoop

! ==============================================================================
!     Default
! ==============================================================================

      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(3X,A)') sectionString
        END IF ! global%verbLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)        
    END SELECT ! TRIM
  END DO emptyLoop
 
! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag        
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Check that periodic patches have transform
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    IF ( pPatch%bcType == BC_PERIODIC ) THEN 
      IF ( pPatch%transformFlag .EQV. .FALSE. ) THEN 
! TEMPORARY
        WRITE(*,*) 'ERROR! Periodic patch without transform!'
        STOP
! END TEMPORARY      
      END IF ! pPatch%transformFlag
    END IF ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading periodicity transforms done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SYPE_ReadTransforms








! ******************************************************************************
!
! Purpose: Determine whether have symmetry and/or periodic patches by peeking
!   into boundary-condition file without actually reading any information.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: 
!   1. Required because may need to know whether have patches before patch data
!      structure is created, so cannot use RFLU_SYPE_HaveSyPePatches. This need
!      occurs in rflupart when converting external grid representations to that
!      used by Rocflu and need to estimate max number of faces, vertices, etc. 
!      Those max values will have to be larger than when do not have sype 
!      patches because virtual cells and vertices are going to be added.
!
! ******************************************************************************

SUBROUTINE RFLU_SYPE_SetSyPePatchesFlag(global)                                   

  USE ModBuildFileNames, ONLY: BuildFileNamePlain

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global   

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName
  CHARACTER(256) :: line
  INTEGER :: errorFlag,iPatch,loopCounter

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_SYPE_SetSyPePatchesFlag',&
  'RFLU_ModSymmetryPeriodic.F90')

! ******************************************************************************
! Set variables
! ******************************************************************************

  loopCounter = 0

  global%syPePatchesFlag = .FALSE.

! ******************************************************************************
! Open file
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.bc',iFileName)

  OPEN(IF_INPUT,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error
    
! ******************************************************************************
! Look for 
! ******************************************************************************

  keyWordLoop: DO
    READ(IF_INPUT,'(A256)',IOSTAT=errorFlag) line

    IF ( errorFlag > 0 ) THEN ! Error occurred
      CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName))
    ELSE IF ( errorFlag < 0 ) THEN ! Encountered end of file
      EXIT keyWordLoop
    END IF ! errorFlag

    SELECT CASE( TRIM(line) )
      CASE ('# BC_PERIODIC')
        global%syPePatchesFlag = .TRUE.
        EXIT keyWordLoop           
      CASE ('# BC_SYMMETRY')
        global%syPePatchesFlag = .TRUE.
        EXIT keyWordLoop        
      CASE ('# END')
        EXIT keyWordLoop
    END SELECT ! TRIM(line)

    loopCounter = loopCounter + 1 
    
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN ! Prevent infinite loop
      CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
    END IF ! loopCounter
  END DO keyWordLoop

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_INPUT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SYPE_SetSyPePatchesFlag







! ******************************************************************************
!
! Purpose: Write transforms between periodic patches.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. Should only be called for serial region.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SYPE_WriteTransforms(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNamePlain  
 
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

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: errorFlag,iFile,iPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SYPE_WriteTransforms',&
  'RFLU_ModSymmetryPeriodic.F90')

  IF ( global%myProcid == MASTERPROC ) THEN
    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing periodicity transforms...'
      
      IF ( global%verbLevel > VERBOSE_LOW ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel
    END IF ! global%verbLevel 
  END IF ! global%myProcid

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
 
  pGrid => pRegion%grid

! ******************************************************************************
! Open file
! ******************************************************************************

  iFile = IF_PTMATRIX

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.ptm',iFileName)

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)   
  global%error = errorFlag        
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Write header
! ******************************************************************************

  sectionString = '# ROCFLU periodic-transformation matrix file'
  WRITE(iFile,'(A)') TRIM(sectionString)  
     
! ******************************************************************************
! Loop over patches, write transformation matrix
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
       
    IF ( pPatch%bcType == BC_PERIODIC ) THEN
      WRITE(iFile,'(A)') '# Patch'
      
      WRITE(iFile,'(I3)') pPatch%iPatchGlobal
      
      WRITE(iFile,'(4(E23.16))') pPatch%tm(XCOORD,XCOORD), & 
                                 pPatch%tm(XCOORD,YCOORD), & 
                                 pPatch%tm(XCOORD,ZCOORD), & 
                                 pPatch%tm(XCOORD,XYZMAG)
      WRITE(iFile,'(4(E23.16))') pPatch%tm(YCOORD,XCOORD), & 
                                 pPatch%tm(YCOORD,YCOORD), & 
                                 pPatch%tm(YCOORD,ZCOORD), & 
                                 pPatch%tm(YCOORD,XYZMAG)
      WRITE(iFile,'(4(E23.16))') pPatch%tm(ZCOORD,XCOORD), & 
                                 pPatch%tm(ZCOORD,YCOORD), & 
                                 pPatch%tm(ZCOORD,ZCOORD), & 
                                 pPatch%tm(ZCOORD,XYZMAG)
      WRITE(iFile,'(4(E23.16))') pPatch%tm(XYZMAG,XCOORD), & 
                                 pPatch%tm(XYZMAG,YCOORD), & 
                                 pPatch%tm(XYZMAG,ZCOORD), & 
                                 pPatch%tm(XYZMAG,XYZMAG)                                                                                                                                                                     
    END IF ! pPatch%bcType    
  END DO ! iPatch

! ******************************************************************************
! Write footer
! ******************************************************************************

  sectionString = '# End'
  WRITE(iFile,'(A)') TRIM(sectionString)  

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag        
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing periodicity transforms done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SYPE_WriteTransforms







END MODULE RFLU_ModSymmetryPeriodic

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModSymmetryPeriodic.F90,v $
! Revision 1.8  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/12/18 02:31:06  haselbac
! Bug fix: Inconsistent IFs in P2VC routines
!
! Revision 1.5  2006/08/18 14:03:15  haselbac
! Added read/write routines for transform matrix
!
! Revision 1.4  2006/04/17 19:57:05  haselbac
! Bug fix: Removed setting of pBorder in destroyign of patch-to-virtual cell list
!
! Revision 1.3  2006/04/12 16:10:05  haselbac
! Increased nLayers for 2nd order virtual cells to 1
!
! Revision 1.2  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.1  2006/03/25 21:38:54  haselbac
! Initial revision
!
! ******************************************************************************

























