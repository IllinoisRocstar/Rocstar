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
! Purpose: Suite of routines to partition region.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModPartitionRegion.F90,v 1.17 2008/12/06 08:45:03 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPartitionRegion

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_PART_AddVirtualCells, &
            RFLU_PART_BuildBorderFaceList, & 
            RFLU_PART_BuildCellLists, & 
            RFLU_PART_BuildPatchLists, &           
            RFLU_PART_BuildVertexData, & 
            RFLU_PART_BuildVertexLists, & 
            RFLU_PART_BuildReg2CellMap, &
            RFLU_PART_CreateCellLists, & 
            RFLU_PART_CreatePatchLists, &  
            RFLU_PART_CreateReg2CellMap, &
            RFLU_PART_DestroyCellData, & 
            RFLU_PART_DestroyBorderFaceList, & 
            RFLU_PART_DestroyCellLists, & 
            RFLU_PART_DestroyPatchLists, & 
            RFLU_PART_DestroyVertexData, &  
            RFLU_PART_PartitionRegion, & 
            RFLU_PART_RenumberVertexLists 

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPartitionRegion.F90,v $ $Revision: 1.17 $'


! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Add virtual cells.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_AddVirtualCells(pRegion,pRegionSerial)                                   

    USE ModSortSearch

    USE RFLU_ModHashTable
    USE RFLU_ModPartitionRegionUtils, ONLY: RFLU_PART_AddVirtualCellsInv1, & 
                                            RFLU_PART_AddVirtualCellsInv2

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

    INTEGER :: errorFlag,i,icg,icg2,icl,ict,iLayer,iLoc,iReg,j,key, &
               nCellsVirt,nCellsVirtMax,nLayers
    INTEGER, DIMENSION(:), ALLOCATABLE :: vc
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegionSerial%global

    CALL RegisterFunction(global,'RFLU_PART_AddVirtualCells',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Adding virtual cells...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid

! ******************************************************************************
!   Build list of virtual cells
! ******************************************************************************

    nCellsVirtMax = pGrid%nCellsMax - pGrid%nCells
    
    ALLOCATE(vc(nCellsVirtMax),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vc')
    END IF ! global%error         
           
    IF ( pRegionSerial%mixtInput%spaceOrder > 1 ) THEN 
      CALL RFLU_PART_AddVirtualCellsInv2(pRegion,pRegionSerial,vc, & 
                                         nCellsVirtMax,nCellsVirt)
    ELSE 
      CALL RFLU_PART_AddVirtualCellsInv1(pRegion,pRegionSerial,vc, &
                                         nCellsVirtMax,nCellsVirt)
    END IF ! pRegionSerial%mixtInput%spaceOrder    
                    
! ******************************************************************************
!   Loop over list of virtual cells, add to connectivity lists according to 
!   cell type
! ******************************************************************************

    DO i = 1,nCellsVirt
      icg = vc(i)

      ict = pGridSerial%cellGlob2Loc(1,icg)
      icl = pGridSerial%cellGlob2Loc(2,icg)

! ==============================================================================
!     Specify connectivity and set local-to-global mapping
! ==============================================================================

      SELECT CASE ( ict )           

! ------------------------------------------------------------------------------
!       Tetrahedra 
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_TET )
          IF ( pGrid%nTetsTot == pGrid%nTetsMax ) THEN 
            global%warnCounter = global%warnCounter + 1

	    IF ( global%verbLevel > VERBOSE_LOW ) THEN
	      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
        	 '*** WARNING *** About to exceed tetrahedra list dimensions.'
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                 '                Increasing list dimensions and continuing.'
	    END IF ! global

            CALL RFLU_PART_RecreateCellList(global,4,pGrid%nTetsMax, & 
                                            pGrid%tet2v,pGrid%tet2CellGlob)
          END IF ! pGrid%nTetsTot

          pGrid%nCellsTot = pGrid%nCellsTot + 1  
          pGrid%nTetsTot  = pGrid%nTetsTot  + 1
            
          pGrid%tet2v(1,pGrid%nTetsTot) = pGridSerial%tet2v(1,icl)
          pGrid%tet2v(2,pGrid%nTetsTot) = pGridSerial%tet2v(2,icl)
          pGrid%tet2v(3,pGrid%nTetsTot) = pGridSerial%tet2v(3,icl)
          pGrid%tet2v(4,pGrid%nTetsTot) = pGridSerial%tet2v(4,icl)

          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_TET
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nTetsTot                

          pGrid%tet2CellGlob(pGrid%nTetsTot) = pGrid%nCellsTot

! ------------------------------------------------------------------------------
!       Hexahedra 
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_HEX )
          IF ( pGrid%nHexsTot == pGrid%nHexsMax ) THEN 
            global%warnCounter = global%warnCounter + 1

	    IF ( global%verbLevel > VERBOSE_LOW ) THEN
	      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
        	 '*** WARNING *** About to exceed hexahedra list dimensions.'
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                 '                Increasing list dimensions and continuing.'
	    END IF ! global

            CALL RFLU_PART_RecreateCellList(global,8,pGrid%nHexsMax, & 
                                            pGrid%hex2v,pGrid%hex2CellGlob)
          END IF ! pGrid%nHexsTot
 
          pGrid%nCellsTot = pGrid%nCellsTot + 1  
          pGrid%nHexsTot  = pGrid%nHexsTot  + 1
  
          pGrid%hex2v(1,pGrid%nHexsTot) = pGridSerial%hex2v(1,icl)
          pGrid%hex2v(2,pGrid%nHexsTot) = pGridSerial%hex2v(2,icl)
          pGrid%hex2v(3,pGrid%nHexsTot) = pGridSerial%hex2v(3,icl)
          pGrid%hex2v(4,pGrid%nHexsTot) = pGridSerial%hex2v(4,icl)
          pGrid%hex2v(5,pGrid%nHexsTot) = pGridSerial%hex2v(5,icl)
          pGrid%hex2v(6,pGrid%nHexsTot) = pGridSerial%hex2v(6,icl)
          pGrid%hex2v(7,pGrid%nHexsTot) = pGridSerial%hex2v(7,icl)
          pGrid%hex2v(8,pGrid%nHexsTot) = pGridSerial%hex2v(8,icl)

          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_HEX
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nHexsTot                

          pGrid%hex2CellGlob(pGrid%nHexsTot) = pGrid%nCellsTot                           

! ------------------------------------------------------------------------------
!       Prisms
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_PRI ) 
          IF ( pGrid%nPrisTot == pGrid%nPrisMax ) THEN 
            global%warnCounter = global%warnCounter + 1

	    IF ( global%verbLevel > VERBOSE_LOW ) THEN
	      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
        	 '*** WARNING *** About to exceed prism list dimensions.'
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                 '                Increasing list dimensions and continuing.'
	    END IF ! global

            CALL RFLU_PART_RecreateCellList(global,6,pGrid%nPrisMax, & 
                                            pGrid%pri2v,pGrid%pri2CellGlob)
          END IF ! pGrid%nPrisTot

          pGrid%nCellsTot = pGrid%nCellsTot + 1              
          pGrid%nPrisTot  = pGrid%nPrisTot  + 1
           
          pGrid%pri2v(1,pGrid%nPrisTot) = pGridSerial%pri2v(1,icl)
          pGrid%pri2v(2,pGrid%nPrisTot) = pGridSerial%pri2v(2,icl)
          pGrid%pri2v(3,pGrid%nPrisTot) = pGridSerial%pri2v(3,icl)
          pGrid%pri2v(4,pGrid%nPrisTot) = pGridSerial%pri2v(4,icl) 
          pGrid%pri2v(5,pGrid%nPrisTot) = pGridSerial%pri2v(5,icl) 
          pGrid%pri2v(6,pGrid%nPrisTot) = pGridSerial%pri2v(6,icl) 

          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_PRI
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nPrisTot                

          pGrid%pri2CellGlob(pGrid%nPrisTot) = pGrid%nCellsTot                                                          

! ------------------------------------------------------------------------------
!       Pyramids
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_PYR ) 
          IF ( pGrid%nPyrsTot == pGrid%nPyrsMax ) THEN 
            global%warnCounter = global%warnCounter + 1

	    IF ( global%verbLevel > VERBOSE_LOW ) THEN
	      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
        	 '*** WARNING *** About to exceed pyramid list dimensions.'
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                 '                Increasing list dimensions and continuing.'
	    END IF ! global

            CALL RFLU_PART_RecreateCellList(global,5,pGrid%nPyrsMax, & 
                                            pGrid%pyr2v,pGrid%pyr2CellGlob)
          END IF ! pGrid%nPyrsTot

          pGrid%nCellsTot = pGrid%nCellsTot + 1              
          pGrid%nPyrsTot  = pGrid%nPyrsTot  + 1

          pGrid%pyr2v(1,pGrid%nPyrsTot) = pGridSerial%pyr2v(1,icl)
          pGrid%pyr2v(2,pGrid%nPyrsTot) = pGridSerial%pyr2v(2,icl)
          pGrid%pyr2v(3,pGrid%nPyrsTot) = pGridSerial%pyr2v(3,icl)
          pGrid%pyr2v(4,pGrid%nPyrsTot) = pGridSerial%pyr2v(4,icl) 
          pGrid%pyr2v(5,pGrid%nPyrsTot) = pGridSerial%pyr2v(5,icl)

          pGrid%cellGlob2Loc(1,pGrid%nCellsTot) = CELL_TYPE_PYR
          pGrid%cellGlob2Loc(2,pGrid%nCellsTot) = pGrid%nPyrsTot                

          pGrid%pyr2CellGlob(pGrid%nPyrsTot) = pGrid%nCellsTot                

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)          
      END SELECT ! pGridSerial%cellGlob2Loc          

! ==============================================================================
!     Update mapping from partitioned cells to serial cells         
! ==============================================================================

      pGrid%pc2sc(pGrid%nCellsTot) = icg      
    END DO ! icl

! ******************************************************************************
!   Deallocate temporary memory for virtual cells
! ******************************************************************************

    DEALLOCATE(vc,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'virtCells')
    END IF ! global%error 
           
! ******************************************************************************
!   Write information about numbers of virtual cells
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
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Adding virtual cells done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_AddVirtualCells







! ******************************************************************************
!
! Purpose: Build border face lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_BuildBorderFaceList(pRegion)                                   

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

    CHARACTER(CHRLEN) :: errorString
    INTEGER :: c1,c2,errorFlag,ifg,ifl,iReg,iReg1,iReg2,nFacesCut
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: avf
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_BuildBorderFaceList',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building border face lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Extract list of ALL actual-virtual faces. NOTE number of cut faces is known
!   from partitioning call.  
! ******************************************************************************

    ALLOCATE(pGrid%avf(3,pGrid%nFacesCut),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%avf')
    END IF ! global%error

    nFacesCut = 0

    DO ifg = 1,pGrid%nFaces + pGrid%nFacesVV
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)
      
      IF ( pGrid%sc2r(c1) /= pGrid%sc2r(c2) ) THEN
        nFacesCut = nFacesCut + 1
       
        pGrid%avf(1,nFacesCut) = ifg
        pGrid%avf(2,nFacesCut) = pGrid%sc2r(c1)
        pGrid%avf(3,nFacesCut) = pGrid%sc2r(c2)                        
      END IF ! pGrid%sc2r
    END DO ! ifg
      
    IF ( nFacesCut /= pGrid%nFacesCut ) THEN 
      WRITE(errorString,'(2(1X,I6))') nFacesCut,pGrid%nFacesCut 
      CALL ErrorStop(global,ERR_NFACESCUT_INVALID,__LINE__,TRIM(errorString))
    END IF ! nFacesCut

! ******************************************************************************
!   Cast list of actual-virtual faces in CSR form so each region can access
!   its actual-virtual faces directly
! ******************************************************************************

! ==============================================================================
!   Allocate memory for CSR access array and count number of cut faces for each
!   region. NOTE this leads to duplicated faces, but more convenient to access
!   that way.
! ==============================================================================
      
    ALLOCATE(pGrid%avfCSRInfo(global%nRegionsLocal),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%avfCSRInfo')
    END IF ! global%error
      
    DO iReg = 1,global%nRegionsLocal 
      pGrid%avfCSRInfo(iReg) = 0
    END DO ! iReg  
      
    DO ifl = 1,pGrid%nFacesCut
      iReg1 = pGrid%avf(2,ifl)
      iReg2 = pGrid%avf(3,ifl)
    
      pGrid%avfCSRInfo(iReg1) = pGrid%avfCSRInfo(iReg1) + 1
      pGrid%avfCSRInfo(iReg2) = pGrid%avfCSRInfo(iReg2) + 1       
    END DO ! ifl  
    
! ==============================================================================
!   Sum number of cut faces for each region (after adding offset) so that can 
!   count down when building list and get CSAR access array pointing to first 
!   entry for each region.
! ==============================================================================    
    
    pGrid%avfCSRInfo(1) = pGrid%avfCSRInfo(1) + 1  
    
    DO iReg = 2,global%nRegionsLocal
      pGrid%avfCSRInfo(iReg) = pGrid%avfCSRInfo(iReg  ) & 
                             + pGrid%avfCSRInfo(iReg-1) 
    END DO ! iReg      

! ==============================================================================
!   Build actual-virtual face list in CSR format
! ==============================================================================
        
    ALLOCATE(pGrid%avfCSR(2*pGrid%nFacesCut),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%avfCSR')
    END IF ! global%error

    DO ifl = 1,pGrid%nFacesCut
      ifg   = pGrid%avf(1,ifl)
      iReg1 = pGrid%avf(2,ifl)
      iReg2 = pGrid%avf(3,ifl)
    
      pGrid%avfCSRInfo(iReg1) = pGrid%avfCSRInfo(iReg1) - 1
      pGrid%avfCSRInfo(iReg2) = pGrid%avfCSRInfo(iReg2) - 1      
        
      pGrid%avfCSR(pGrid%avfCSRInfo(iReg1)) = ifg
      pGrid%avfCSR(pGrid%avfCSRInfo(iReg2)) = ifg       
    END DO ! ifl      

! ******************************************************************************
!   Deallocate original list 
! ******************************************************************************
    
    DEALLOCATE(pGrid%avf,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%avf')
    END IF ! global%error
      
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building border face lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_BuildBorderFaceList








! ******************************************************************************
!
! Purpose: Build cell lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_BuildCellLists(pRegion,pRegionSerial)

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

    INTEGER :: errorFlag,i,iBeg,icg,icg2,icl,icl2,ict,iEnd,iReg,nHexsAct, &
               nHexsVir,nPrisAct,nPrisVir,nPyrsAct,nPyrsVir,nTetsAct,nTetsVir
    TYPE(t_grid), POINTER :: pGrid,pGridSerial 
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_BuildCellLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid
                
! ******************************************************************************
!   Get cell connectivity. NOTE at this stage the cell connectivity is still 
!   given in terms of the vertex numbers of the serial region.
! ******************************************************************************
   
    nTetsAct = 0
    nTetsVir = 0
    
    nHexsAct = 0
    nHexsVir = 0
    
    nPrisAct = 0
    nPrisVir = 0
    
    nPyrsAct = 0
    nPyrsVir = 0 
      
    iReg = pRegion%iRegionGlobal        
    iBeg = pGridSerial%r2pcCSRInfo(iReg)
            
    IF ( iReg /= global%nRegionsLocal ) THEN     
      iEnd = pGridSerial%r2pcCSRInfo(iReg+1)-1 
    ELSE 
      iEnd = pGridSerial%nCellsTot
    END IF ! iReg      
      
    DO i = iBeg,iEnd      
      icg = pGridSerial%r2pcCSR(i) 
      
      ict = pGridSerial%cellGlob2Loc(1,icg)
      icl = pGridSerial%cellGlob2Loc(2,icg)

      SELECT CASE ( ict ) 
      
! ------------------------------------------------------------------------------
!       Tetrahedra
! ------------------------------------------------------------------------------      
      
        CASE ( CELL_TYPE_TET ) 
          IF ( icl <= pGridSerial%nTets ) THEN 
            nTetsAct = nTetsAct + 1                            
            icl2     = nTetsAct                                     
          ELSE 
            nTetsVir = nTetsVir + 1            
            icl2     = nTetsVir + pGrid%nTets            
          END IF ! icl
                    
          pGrid%tet2v(1,icl2) = pGridSerial%tet2v(1,icl)
          pGrid%tet2v(2,icl2) = pGridSerial%tet2v(2,icl)
          pGrid%tet2v(3,icl2) = pGridSerial%tet2v(3,icl)
          pGrid%tet2v(4,icl2) = pGridSerial%tet2v(4,icl)          
          
          icg2 = pGrid%tet2CellGlob(icl2) 

          pGrid%pc2sc(icg2) = icg          

! ------------------------------------------------------------------------------
!       Hexahedra
! ------------------------------------------------------------------------------      
          
        CASE ( CELL_TYPE_HEX ) 
          IF ( icl <= pGridSerial%nHexs ) THEN 
            nHexsAct = nHexsAct + 1 
            icl2     = nHexsAct       
          ELSE 
            nHexsVir = nHexsVir + 1
            icl2     = nHexsVir + pGrid%nHexs           
          END IF ! icl    
          
          pGrid%hex2v(1,icl2) = pGridSerial%hex2v(1,icl)
          pGrid%hex2v(2,icl2) = pGridSerial%hex2v(2,icl)
          pGrid%hex2v(3,icl2) = pGridSerial%hex2v(3,icl)
          pGrid%hex2v(4,icl2) = pGridSerial%hex2v(4,icl) 
          pGrid%hex2v(5,icl2) = pGridSerial%hex2v(5,icl) 
          pGrid%hex2v(6,icl2) = pGridSerial%hex2v(6,icl) 
          pGrid%hex2v(7,icl2) = pGridSerial%hex2v(7,icl) 
          pGrid%hex2v(8,icl2) = pGridSerial%hex2v(8,icl)                                                      

          icg2 = pGrid%hex2CellGlob(icl2) 

          pGrid%pc2sc(icg2) = icg 
          
! ------------------------------------------------------------------------------
!       Prisms
! ------------------------------------------------------------------------------      

        CASE ( CELL_TYPE_PRI ) 
          IF ( icl <= pGridSerial%nPris ) THEN 
            nPrisAct = nPrisAct + 1
            icl2     = nPrisAct        
          ELSE 
            nPrisVir = nPrisVir + 1            
            icl2     = nPrisVir + pGrid%nPris
          END IF ! icl
          
          pGrid%pri2v(1,icl2) = pGridSerial%pri2v(1,icl)
          pGrid%pri2v(2,icl2) = pGridSerial%pri2v(2,icl)
          pGrid%pri2v(3,icl2) = pGridSerial%pri2v(3,icl)
          pGrid%pri2v(4,icl2) = pGridSerial%pri2v(4,icl) 
          pGrid%pri2v(5,icl2) = pGridSerial%pri2v(5,icl) 
          pGrid%pri2v(6,icl2) = pGridSerial%pri2v(6,icl)         

          icg2 = pGrid%pri2CellGlob(icl2) 

          pGrid%pc2sc(icg2) = icg
          
! ------------------------------------------------------------------------------
!       Pyramids
! ------------------------------------------------------------------------------      

        CASE ( CELL_TYPE_PYR ) 
          IF ( icl <= pGridSerial%nPyrs ) THEN 
            nPyrsAct = nPyrsAct + 1 
            icl2     = nPyrsAct           
          ELSE 
            nPyrsVir = nPyrsVir + 1            
            icl2     = nPyrsVir + pGrid%nPyrs
          END IF ! icl     
          
          pGrid%pyr2v(1,icl2) = pGridSerial%pyr2v(1,icl)
          pGrid%pyr2v(2,icl2) = pGridSerial%pyr2v(2,icl)
          pGrid%pyr2v(3,icl2) = pGridSerial%pyr2v(3,icl)
          pGrid%pyr2v(4,icl2) = pGridSerial%pyr2v(4,icl) 
          pGrid%pyr2v(5,icl2) = pGridSerial%pyr2v(5,icl)             

          icg2 = pGrid%pyr2CellGlob(icl2) 

          pGrid%pc2sc(icg2) = icg
          
! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------      

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! ict
    END DO ! i      
      
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_BuildCellLists











! ******************************************************************************
!
! Purpose: Build patch lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region
!
! Output: None.
!
! Notes:
!   1. At the end of this routine, the face lists of the partitioned region are
!      still in terms of the vertex numbers of the serial region. 
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_BuildPatchLists(pRegion,pRegionSerial)

    USE ModSortSearch

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

    INTEGER :: errorFlag,icg,icgMax,icgMin,icg2,icl,ict,ifl,ifl2,iLoc,iPatch, &
               iReg,nBQuadsAct,nBQuadsVir,nBTrisAct,nBTrisVir,offs,v1g,v2g, &
               v3g,v4g
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch,pPatchSerial     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_BuildPatchLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building patch lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid
        
! ******************************************************************************       
!   Loop over serial patches and set variables if that patch exists on this 
!   partition
! ******************************************************************************       

    icgMin = MINVAL(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot))
    icgMax = MAXVAL(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot))    

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)    

      pPatchSerial => pRegionSerial%patches(pPatch%iPatchGlobal)
        
      pPatch%bcType = pPatchSerial%bcType ! Required for sype cases  
        
      nBTrisAct = 0
      nBTrisVir = 0   

      nBQuadsAct = 0
      nBQuadsVir = 0 
      
      offs = pGrid%pbf2sbfCSRInfo(iPatch) - 1     
            
! ==============================================================================
!     Non-virtual patch. Treat patches differently because every cell in 
!     partitioned region must be adjacent to virtual patches, so searching is 
!     not necessary. NOTE that it is not necessary to distinguish between the 
!     two cases for correct running of the code, it is only done for efficiency.
! ==============================================================================

      IF ( pPatchSerial%bcType /= BC_VIRTUAL ) THEN
        DO ifl = 1,pPatchSerial%nBFacesTot        
          icg  = pPatchSerial%bf2c(ifl)                
          iReg = pGridSerial%sc2r(icg)

! ------------------------------------------------------------------------------
!         Find whether cell adjacent to boundary face exists in cell list                
! ------------------------------------------------------------------------------ 

          IF ( icg >= icgMin .AND. icg <= icgMax ) THEN ! Search
            CALL BinarySearchInteger(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot), &
                                     pGrid%nBCellsTot,icg,iLoc)
          ELSE ! No need to search
	    iLoc = ELEMENT_NOT_FOUND				   
          END IF ! icg

! ------- Cell exists, so add face to boundary-face lists ----------------------                  

          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
            icg2 = pGrid%sbc2pc(2,iLoc)                           

            IF ( icg2 <= pGrid%nCells ) THEN ! Actual-boundary face
              IF ( pPatchSerial%bf2v(4,ifl) == VERT_NONE ) THEN                                  
                nBTrisAct = nBTrisAct + 1             

                v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))                            

                pPatch%bTri2v(1,nBTrisAct) = v1g
                pPatch%bTri2v(2,nBTrisAct) = v2g
                pPatch%bTri2v(3,nBTrisAct) = v3g 

                ifl2 = nBTrisAct          

                pGrid%pbf2sbfCSR(offs+ifl2) = ifl
              ELSE 
                nBQuadsAct = nBQuadsAct + 1             

                v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))
                v4g = pPatchSerial%bv(pPatchSerial%bf2v(4,ifl))                                           

                pPatch%bQuad2v(1,nBQuadsAct) = v1g
                pPatch%bQuad2v(2,nBQuadsAct) = v2g
                pPatch%bQuad2v(3,nBQuadsAct) = v3g
                pPatch%bQuad2v(4,nBQuadsAct) = v4g

                ifl2 = pPatch%nBTrisTot + nBQuadsAct

                pGrid%pbf2sbfCSR(offs + ifl2) = ifl              
              END IF ! pPatchSerial
            ELSE ! Virtual-boundary face
              IF ( pPatchSerial%bf2v(4,ifl) == VERT_NONE ) THEN          
                nBTrisVir = nBTrisVir + 1

                v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))                            

                pPatch%bTri2v(1,pPatch%nBTris+nBTrisVir) = v1g
                pPatch%bTri2v(2,pPatch%nBTris+nBTrisVir) = v2g
                pPatch%bTri2v(3,pPatch%nBTris+nBTrisVir) = v3g 

                ifl2 = pPatch%nBTris + nBTrisVir

                pGrid%pbf2sbfCSR(offs+ifl2) = ifl               
              ELSE 
                nBQuadsVir = nBQuadsVir + 1 

                v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))
                v4g = pPatchSerial%bv(pPatchSerial%bf2v(4,ifl))                                           

                pPatch%bQuad2v(1,pPatch%nBQuads+nBQuadsVir) = v1g
                pPatch%bQuad2v(2,pPatch%nBQuads+nBQuadsVir) = v2g
                pPatch%bQuad2v(3,pPatch%nBQuads+nBQuadsVir) = v3g
                pPatch%bQuad2v(4,pPatch%nBQuads+nBQuadsVir) = v4g

                ifl2 = pPatch%nBTrisTot + pPatch%nBQuads + nBQuadsVir

                pGrid%pbf2sbfCSR(offs+ifl2) = ifl                
              END IF ! pPatchSerial          
            END IF ! icg2

! ------- Cell does not exist although it should -------------------------------

          ELSE 
            IF ( iReg == pRegion%iRegionGlobal ) THEN           
              CALL ErrorStop(global,ERR_CELL_NOT_FOUND,__LINE__)
            END IF ! iReg
          END IF ! iLoc          
        END DO ! ifl  
      
! ==============================================================================
!     Virtual patch
! ==============================================================================      

      ELSE 
         
! ------------------------------------------------------------------------------
!       Allocate and build sorted bf2c list along with key so that can get 
!       local boundary face index on serial grid for given global cell index
! ------------------------------------------------------------------------------         
         
        ALLOCATE(pPatchSerial%bf2cSorted(pPatchSerial%nBFacesTot), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                         'pPatchSerial%bf2cSorted')
        END IF ! global%error
             
        ALLOCATE(pPatchSerial%bf2cSortedKeys(pPatchSerial%nBFacesTot), & 
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                         'pPatchSerial%bf2cSortedKeys')
        END IF ! global%error
             
        DO ifl = 1,pPatchSerial%nBFacesTot        
          pPatchSerial%bf2cSorted(ifl)     = pPatchSerial%bf2c(ifl)
          pPatchSerial%bf2cSortedKeys(ifl) = ifl          
        END DO ! ifl
               
        CALL QuickSortIntegerInteger(pPatchSerial%bf2cSorted, & 
                                     pPatchSerial%bf2cSortedKeys, & 
                                     pPatchSerial%nBFacesTot)     
        
! ------------------------------------------------------------------------------
!       Loop over cells in partitioned region. Every single one of them must be 
!       on virtual patch.
! ------------------------------------------------------------------------------        
        
        DO icg = 1,pGrid%nCellsTot
          icg2 = pGrid%pc2sc(icg)
                   
! ------- Check on sorted bf2c list --------------------------------------------          
          
          IF ( pPatchSerial%bf2cSorted(icg2) == icg2 ) THEN           
            ifl = pPatchSerial%bf2cSortedKeys(icg2)

            ict = pGrid%cellGlob2Loc(1,icg)
            icl = pGrid%cellGlob2Loc(2,icg)
            
            SELECT CASE ( ict )
              CASE ( CELL_TYPE_HEX )  
                IF ( icl <= pGrid%nHexs ) THEN ! Actual quad boundary face
                  nBQuadsAct = nBQuadsAct + 1             

                  v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                  v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                  v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))
                  v4g = pPatchSerial%bv(pPatchSerial%bf2v(4,ifl))                                           

                  pPatch%bQuad2v(1,nBQuadsAct) = v1g
                  pPatch%bQuad2v(2,nBQuadsAct) = v2g
                  pPatch%bQuad2v(3,nBQuadsAct) = v3g
                  pPatch%bQuad2v(4,nBQuadsAct) = v4g

                  ifl2 = pPatch%nBTrisTot + nBQuadsAct

                  pGrid%pbf2sbfCSR(offs + ifl2) = ifl                    
                ELSE ! Virtual quad boundary face
                  nBQuadsVir = nBQuadsVir + 1 

                  v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                  v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                  v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))
                  v4g = pPatchSerial%bv(pPatchSerial%bf2v(4,ifl))                                           

                  pPatch%bQuad2v(1,pPatch%nBQuads+nBQuadsVir) = v1g
                  pPatch%bQuad2v(2,pPatch%nBQuads+nBQuadsVir) = v2g
                  pPatch%bQuad2v(3,pPatch%nBQuads+nBQuadsVir) = v3g
                  pPatch%bQuad2v(4,pPatch%nBQuads+nBQuadsVir) = v4g

                  ifl2 = pPatch%nBTrisTot + pPatch%nBQuads + nBQuadsVir

                  pGrid%pbf2sbfCSR(offs+ifl2) = ifl                         
                END IF ! icl 
              CASE ( CELL_TYPE_PRI ) 
                IF ( icl <= pGrid%nPris ) THEN ! Actual tri boundary face
                  nBTrisAct = nBTrisAct + 1             

                  v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                  v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                  v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))                            

                  pPatch%bTri2v(1,nBTrisAct) = v1g
                  pPatch%bTri2v(2,nBTrisAct) = v2g
                  pPatch%bTri2v(3,nBTrisAct) = v3g 

                  ifl2 = nBTrisAct          

                  pGrid%pbf2sbfCSR(offs+ifl2) = ifl                
                ELSE ! Virtual tri boundary face                
                  nBTrisVir = nBTrisVir + 1

                  v1g = pPatchSerial%bv(pPatchSerial%bf2v(1,ifl))
                  v2g = pPatchSerial%bv(pPatchSerial%bf2v(2,ifl))
                  v3g = pPatchSerial%bv(pPatchSerial%bf2v(3,ifl))                            

                  pPatch%bTri2v(1,pPatch%nBTris+nBTrisVir) = v1g
                  pPatch%bTri2v(2,pPatch%nBTris+nBTrisVir) = v2g
                  pPatch%bTri2v(3,pPatch%nBTris+nBTrisVir) = v3g 

                  ifl2 = pPatch%nBTris + nBTrisVir

                  pGrid%pbf2sbfCSR(offs+ifl2) = ifl                                               
                END IF ! icl                               
              CASE DEFAULT
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! ict

! ------- Getting here means that cell is not adjacent to virtual patch --------
            
          ELSE 
            CALL ErrorStop(global,ERR_BF2CSORTED_INVALID,__LINE__)        
          END IF ! pPatchSerial%bf2cSorted 
        END DO ! icg
                
        DEALLOCATE(pPatchSerial%bf2cSorted,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                         'pPatchSerial%bf2cSorted')
        END IF ! global%error        
        
        DEALLOCATE(pPatchSerial%bf2cSortedKeys,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                         'pPatchSerial%bf2cSortedKeys')
        END IF ! global%error        
      END IF ! pPatchSerial%bcType                                           
    END DO ! iPatch 
                    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building patch lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_BuildPatchLists







! ******************************************************************************
!
! Purpose: Build region-to-cell map.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_BuildReg2CellMap(pRegion)

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

    INTEGER :: errorFlag,icg,iReg
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_BuildReg2CellMap',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building region-to-cell map...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid
                
! ******************************************************************************
!   Count cells in each region and set up info array
! ******************************************************************************
    
    DO icg = 1,pGrid%nCellsTot
      iReg = pGrid%sc2r(icg)
      
      pGrid%r2pcCSRInfo(iReg) = pGrid%r2pcCSRInfo(iReg) + 1
    END DO ! icg
     
    pGrid%r2pcCSRInfo(0) = pGrid%r2pcCSRInfo(0) + 1 
     
    DO iReg = 1,global%nRegionsLocal
      pGrid%r2pcCSRInfo(iReg) = pGrid%r2pcCSRInfo(iReg  ) & 
                              + pGrid%r2pcCSRInfo(iReg-1)
    END DO ! iReg
     
! ******************************************************************************
!   Enter cell into mapping array. NOTE loop backwards because stepping back in
!   r2pcCSRInfo array (from last position to first position); that way get same
!   ordering in cell indices as on serial region.
! ******************************************************************************
    
    DO icg = pGrid%nCellsTot,1,-1
      iReg = pGrid%sc2r(icg)
      
      pGrid%r2pcCSRInfo(iReg) = pGrid%r2pcCSRInfo(iReg) - 1
      
      pGrid%r2pcCSR(pGrid%r2pcCSRInfo(iReg)) = icg
    END DO ! icg    
      
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building region-to-cell map done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_BuildReg2CellMap









! ******************************************************************************
!
! Purpose: Build vertex data.
!
! Description: None.
!
! Input:
!   levels	Level data structure
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_BuildVertexData(pRegion,pRegionSerial)

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

    INTEGER :: errorFlag,ivg,ivg2  
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_BuildVertexData',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building vertex data...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    ALLOCATE(pGrid%xyz(XCOORD:ZCOORD,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%xyz')
    END IF ! global%error         

    DO ivg = 1,pGrid%nVertTot
      ivg2 = pGrid%pv2sv(ivg)
      
      pGrid%xyz(XCOORD,ivg) = pGridSerial%xyz(XCOORD,ivg2)
      pGrid%xyz(YCOORD,ivg) = pGridSerial%xyz(YCOORD,ivg2)
      pGrid%xyz(ZCOORD,ivg) = pGridSerial%xyz(ZCOORD,ivg2)            
    END DO ! ivg

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building vertex data done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_BuildVertexData







! ******************************************************************************
!
! Purpose: Build vertex lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_BuildVertexLists(pRegion,pRegionSerial)

    USE ModSortSearch

    USE RFLU_ModHashTable
    USE RFLU_ModRenumberList, ONLY: RFLU_RenumberList2
    USE RFLU_ModRenumberings, ONLY: RFLU_RNMB_CreatePV2SVMap                                                         

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

    INTEGER :: errorFlag,icl,iLoc,ivg,ivgIndx,ivgStat,ivl,key,nVertAct,nVertInt, &
               nVertVir
    INTEGER, DIMENSION(:), ALLOCATABLE :: indx,tempList
    TYPE(t_grid), POINTER :: pGrid,pGridSerial 
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_BuildVertexLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building vertex lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid

! ******************************************************************************
!   Estimate number of vertices, create hash table, and allocate memory for 
!   mapping from vertices in this partition to serial partition
! ******************************************************************************

    pGrid%nVert    = 0
    pGrid%nVertTot = 0
    
! TO DO
!   Must find improved estimate of maximum number of vertices
! END TO DO     
    pGrid%nVertMax = 1.5_RFREAL*8*pGrid%nCellsTot

    CALL RFLU_CreateHashTable(global,pGrid%nVertMax)  

    ALLOCATE(indx(hashTableSize),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'indx')
    END IF ! global%error
               
    CALL RFLU_RNMB_CreatePV2SVMap(pRegion)

    ALLOCATE(pGrid%vertKind(pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%vertKind')
    END IF ! global%error
    
    DO ivg = 1,pGrid%nVertMax
      pGrid%vertKind(ivg) = VERT_NONE 
    END DO ! ivg
                
! ******************************************************************************
!   Determine number of vertices by constructing hash table of global vertex
!   indices. After sorting, this list will be used as the key for renumbering 
!   the vertices in this partition.
! ******************************************************************************
    
! ==============================================================================
!   Actual cells
! ==============================================================================    
    
    DO icl = 1,pGrid%nTets
      DO ivl = 1,4
        CALL RFLU_HashBuildKey(pGrid%tet2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%tet2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)				  
      END DO ! ivl
    END DO ! icl

    DO icl = 1,pGrid%nHexs
      DO ivl = 1,8
        CALL RFLU_HashBuildKey(pGrid%hex2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%hex2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)
      END DO ! ivl
    END DO ! icl

    DO icl = 1,pGrid%nPris
      DO ivl = 1,6
        CALL RFLU_HashBuildKey(pGrid%pri2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%pri2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)
      END DO ! ivl
    END DO ! icl
    
    DO icl = 1,pGrid%nPyrs
      DO ivl = 1,5
        CALL RFLU_HashBuildKey(pGrid%pyr2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%pyr2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)
      END DO ! ivl
    END DO ! icl 
    
    DO ivg = 1,pGrid%nVertTot
      pGrid%vertKind(ivg) = VERT_KIND_A
    END DO ! ivg
    
! ==============================================================================
!   Virtual cells
! ==============================================================================    
        
    DO icl = pGrid%nTets+1,pGrid%nTetsTot
      DO ivl = 1,4
        CALL RFLU_HashBuildKey(pGrid%tet2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%tet2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)
				  
        IF ( ivgStat == HASHTABLE_ENTRYSTATUS_NEW ) THEN 
	  pGrid%vertKind(pGrid%nVertTot) = VERT_KIND_V
	ELSE 
	  IF ( pGrid%vertKind(ivgIndx) == VERT_KIND_A ) THEN 
	    pGrid%vertKind(ivgIndx) = VERT_KIND_AV
	  END IF ! pGrid%vertKind
	END IF ! ivgStat				  	
      END DO ! ivl
    END DO ! icl

    DO icl = pGrid%nHexs+1,pGrid%nHexsTot
      DO ivl = 1,8
        CALL RFLU_HashBuildKey(pGrid%hex2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%hex2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)
				  
        IF ( ivgStat == HASHTABLE_ENTRYSTATUS_NEW ) THEN 
	  pGrid%vertKind(pGrid%nVertTot) = VERT_KIND_V
	ELSE 
	  IF ( pGrid%vertKind(ivgIndx) == VERT_KIND_A ) THEN 
	    pGrid%vertKind(ivgIndx) = VERT_KIND_AV
	  END IF ! pGrid%vertKind
	END IF ! ivgStat				  
      END DO ! ivl
    END DO ! icl

    DO icl = pGrid%nPris+1,pGrid%nPrisTot
      DO ivl = 1,6
        CALL RFLU_HashBuildKey(pGrid%pri2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%pri2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)
				  
        IF ( ivgStat == HASHTABLE_ENTRYSTATUS_NEW ) THEN 
	  pGrid%vertKind(pGrid%nVertTot) = VERT_KIND_V
	ELSE 
	  IF ( pGrid%vertKind(ivgIndx) == VERT_KIND_A ) THEN 
	    pGrid%vertKind(ivgIndx) = VERT_KIND_AV
	  END IF ! pGrid%vertKind
	END IF ! ivgStat				  
      END DO ! ivl
    END DO ! icl
    
    DO icl = pGrid%nPyrs+1,pGrid%nPyrsTot
      DO ivl = 1,5
        CALL RFLU_HashBuildKey(pGrid%pyr2v(ivl,icl:icl),1,key)          
        CALL RFLU_HashVertexFancy(global,key,pGrid%pyr2v(ivl,icl), & 
                                  pGrid%nVertTot,pGrid%pv2sv, & 
				  indx,ivgStat,ivgIndx)
				  
        IF ( ivgStat == HASHTABLE_ENTRYSTATUS_NEW ) THEN 
	  pGrid%vertKind(pGrid%nVertTot) = VERT_KIND_V
	ELSE 
	  IF ( pGrid%vertKind(ivgIndx) == VERT_KIND_A ) THEN 
	    pGrid%vertKind(ivgIndx) = VERT_KIND_AV
	  END IF ! pGrid%vertKind
	END IF ! ivgStat				  
      END DO ! ivl
    END DO ! icl     

! ******************************************************************************
!   Destroy hash table
! ******************************************************************************

    CALL RFLU_DestroyHashTable(global)

    DEALLOCATE(indx,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'indx')
    END IF ! global%error    
                       
! ==============================================================================
!   Count number of vertices. NOTE at this stage, pGrid%nVert only holds the 
!   number of vertices which are EXCLUSIVELY in the region, and does NOT hold
!   the number of vertices on the interface, as it usually does.
! ==============================================================================    
    
    pGrid%nVert    = 0
    pGrid%nVertInt = 0
    
    DO ivg = 1,pGrid%nVertTot
      SELECT CASE ( pGrid%vertKind(ivg) ) 
        CASE ( VERT_KIND_A ) 
          pGrid%nVert = pGrid%nVert + 1
        CASE ( VERT_KIND_V ) 
        
        CASE ( VERT_KIND_AV ) 
          pGrid%nVertInt = pGrid%nVertInt + 1
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pGrid%vertKind
    END DO ! ivg

! ******************************************************************************
!   Copy list and reorder so have actual vertices at the beginning, followed by
!   actual-virtual vertices, followed by virtual vertices.
! ******************************************************************************
    
    ALLOCATE(tempList(pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'tempList')
    END IF ! global%error
        
    nVertAct = 0
    nVertInt = 0
    nVertVir = 0
    
    DO ivg = 1,pGrid%nVertTot
      SELECT CASE ( pGrid%vertKind(ivg) ) 
        CASE ( VERT_KIND_A ) 
          nVertAct = nVertAct + 1
          ivl = nVertAct
          tempList(nVertAct) = pGrid%pv2sv(ivg)
        CASE ( VERT_KIND_V ) 
          nVertVir = nVertVir + 1
          ivl = pGrid%nVert + pGrid%nVertInt + nVertVir
          tempList(ivl) = pGrid%pv2sv(ivg)        
        CASE ( VERT_KIND_AV ) 
          nVertInt = nVertInt + 1
          ivl = pGrid%nVert + nVertInt
          tempList(ivl) = pGrid%pv2sv(ivg)        
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pGrid%vertKind
    END DO ! ivg

! ******************************************************************************
!   Copy sorted list to give mapping from partitioned vertex to corresponding
!   serial vertex 
! ******************************************************************************
    
    DO ivg = 1,pGrid%nVertTot
      pGrid%pv2sv(ivg) = tempList(ivg)   
    END DO ! ivg    
    
    DEALLOCATE(tempList,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'tempList')
    END IF ! global%error        

! ******************************************************************************
!   Set number of actual vertices
! ******************************************************************************
              
    pGrid%nVert = pGrid%nVert + pGrid%nVertInt          
              
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building vertex lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_BuildVertexLists









! ******************************************************************************
!
! Purpose: Create cell lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_CreateCellLists(pRegion,pRegionSerial)

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

    INTEGER :: errorFlag,i,iBeg,icg,icg2,icl,icl2,ict,iEnd,iReg
    TYPE(t_grid), POINTER :: pGrid,pGridSerial 
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_CreateCellLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating cell lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid
                
! ******************************************************************************
!   Determine number of cells for each type
! ******************************************************************************
    
    iReg = pRegion%iRegionGlobal        
    iBeg = pGridSerial%r2pcCSRInfo(iReg)
            
    IF ( iReg /= global%nRegionsLocal ) THEN     
      iEnd = pGridSerial%r2pcCSRInfo(iReg+1)-1 
    ELSE 
      iEnd = pGridSerial%nCellsTot
    END IF ! iReg
         
    DO i = iBeg,iEnd      
      icg = pGridSerial%r2pcCSR(i)
      
      ict = pGridSerial%cellGlob2Loc(1,icg)
      icl = pGridSerial%cellGlob2Loc(2,icg)

      SELECT CASE ( ict ) 
        CASE ( CELL_TYPE_TET ) 
          IF ( icl <= pGridSerial%nTets ) THEN 
            pGrid%nTets    = pGrid%nTets    + 1
            pGrid%nTetsTot = pGrid%nTetsTot + 1            
          ELSE 
            pGrid%nTetsTot = pGrid%nTetsTot + 1            
          END IF ! icl
        CASE ( CELL_TYPE_HEX ) 
          IF ( icl <= pGridSerial%nHexs ) THEN 
            pGrid%nHexs    = pGrid%nHexs    + 1
            pGrid%nHexsTot = pGrid%nHexsTot + 1            
          ELSE 
            pGrid%nHexsTot = pGrid%nHexsTot + 1            
          END IF ! icl        
        CASE ( CELL_TYPE_PRI ) 
          IF ( icl <= pGridSerial%nPris ) THEN 
            pGrid%nPris    = pGrid%nPris    + 1
            pGrid%nPrisTot = pGrid%nPrisTot + 1            
          ELSE 
            pGrid%nPrisTot = pGrid%nPrisTot + 1            
          END IF ! icl
        CASE ( CELL_TYPE_PYR ) 
          IF ( icl <= pGridSerial%nPyrs ) THEN 
            pGrid%nPyrs    = pGrid%nPyrs    + 1
            pGrid%nPyrsTot = pGrid%nPyrsTot + 1            
          ELSE 
            pGrid%nPyrsTot = pGrid%nPyrsTot + 1            
          END IF ! icl        
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! ict
    END DO ! i    

! ******************************************************************************
!   Set maximum dimensions. NOTE must be generous in setting maximum dimensions
!   because have not yet added virtual cells. NOTE also that a region may not
!   contain any actual cells of a given type, but may contain virtual cells of 
!   that type. 
! ******************************************************************************

    pGrid%nTetsMax  = 4*pGrid%nTetsTot
    pGrid%nHexsMax  = 4*pGrid%nHexsTot 
    pGrid%nPrisMax  = 4*pGrid%nPrisTot 
    pGrid%nPyrsMax  = 4*pGrid%nPyrsTot    

    IF ( pGrid%nTetsMax < 100 ) THEN 
      pGrid%nTetsMax = 10*pGridSerial%nTetsTot/global%nRegions
    END IF ! pGrid%nTetsMax

    IF ( pGrid%nHexsMax < 100 ) THEN 
      pGrid%nHexsMax = 10*pGridSerial%nHexsTot/global%nRegions
    END IF ! pGrid%nHexsMax

    IF ( pGrid%nPrisMax < 100 ) THEN 
      pGrid%nPrisMax = 10*pGridSerial%nPrisTot/global%nRegions
    END IF ! pGrid%nPrisMax

    IF ( pGrid%nPyrsMax < 100 ) THEN 
      pGrid%nPyrsMax = 10*pGridSerial%nPyrsTot/global%nRegions
    END IF ! pGrid%nPyrsMax

! ******************************************************************************
!   Set overall dimensions and write statistics
! ******************************************************************************

    pGrid%nCells    = pGrid%nTets    + pGrid%nHexs    & 
                    + pGrid%nPris    + pGrid%nPyrs
    pGrid%nCellsTot = pGrid%nTetsTot + pGrid%nHexsTot & 
                    + pGrid%nPrisTot + pGrid%nPyrsTot
    pGrid%nCellsMax = pGrid%nTetsMax + pGrid%nHexsMax &
                    + pGrid%nPrisMax + pGrid%nPyrsMax
    
    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cell statistics:'
      WRITE(STDOUT,'(A,5X,A,3X,3(1X,I8))') SOLVER_NAME,'Cells:     ', & 
                                           pGrid%nCells,    & 
                                           pGrid%nCellsTot, & 
                                           pGrid%nCellsMax
      WRITE(STDOUT,'(A,7X,A,1X,3(1X,I8))') SOLVER_NAME,'Tetrahedra:', & 
                                           pGrid%nTets,    &
                                           pGrid%nTetsTot, &
                                           pGrid%nTetsMax
      WRITE(STDOUT,'(A,7X,A,1X,3(1X,I8))') SOLVER_NAME,'Hexahedra: ', & 
                                           pGrid%nHexs,    &
                                           pGrid%nHexsTot, &
                                           pGrid%nHexsMax
      WRITE(STDOUT,'(A,7X,A,1X,3(1X,I8))') SOLVER_NAME,'Prisms:    ', & 
                                           pGrid%nPris,    & 
                                           pGrid%nPrisTot, &
                                           pGrid%nPrisMax
      WRITE(STDOUT,'(A,7X,A,1X,3(1X,I8))') SOLVER_NAME,'Pyramids:  ', & 
                                           pGrid%nPyrs,    &
                                           pGrid%nPyrsTot, &
                                           pGrid%nPyrsMax
    END IF ! global%verbLevel

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

! ==============================================================================
!   Connectivity
! ==============================================================================

    IF ( pGrid%nTetsMax > 0 ) THEN 
      ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
      END IF ! global%error      
    ELSE 
      NULLIFY(pGrid%tet2v)
    END IF ! pGrid%nTetsMax

    IF ( pGrid%nHexsMax > 0 ) THEN 
      ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%hex2v)
    END IF ! pGrid%nHexsMax
    
    IF ( pGrid%nPrisMax > 0 ) THEN 
      ALLOCATE(pGrid%pri2v(6,pGrid%nPrisMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pri2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%pri2v)
    END IF ! pGrid%nPrisMax    
 
    IF ( pGrid%nPyrsMax > 0 ) THEN 
      ALLOCATE(pGrid%pyr2v(5,pGrid%nPyrsMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pyr2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%pyr2v)
    END IF ! pGrid%nPyrsMax 

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating cell lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_CreateCellLists








! ******************************************************************************
!
! Purpose: Create patch list and determine number of patches.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_CreatePatchLists(pRegion,pRegionSerial)

    USE ModSortSearch

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

    INTEGER :: errorFlag,icg,icgMax,icgMin,icg2,ifl,iLoc,iPatch,iReg
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch,pPatchSerial     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_CreatePatchLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating patch lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid
        
! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************       
        
    ALLOCATE(pGrid%patchCounter(2,FACE_TYPE_TRI:FACE_TYPE_QUAD, &
                                pGridSerial%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%patchCounter')
    END IF ! global%error 
      
    DO iPatch = 1,pGridSerial%nPatches
      pGrid%patchCounter(1,FACE_TYPE_TRI ,iPatch) = 0 ! Initial value important
      pGrid%patchCounter(2,FACE_TYPE_TRI ,iPatch) = 0 ! Initial value important      
      pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch) = 0 ! Initial value important
      pGrid%patchCounter(2,FACE_TYPE_QUAD,iPatch) = 0 ! Initial value important             
    END DO ! iPatch           

! ******************************************************************************
!   Determine number of patches in this region
! ******************************************************************************       

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining number of patches...' 
    END IF ! global%verbLevel
    
! ==============================================================================
!   Count number of faces on serial patch which are in this region. NOTE need
!   to loop over all serial faces because need to capture both actual and 
!   virtual faces of given partitioned region. So looping only over faces which 
!   are adjacent to cells in a given partitioned region cannot give virtual 
!   cells.
! ==============================================================================
                 
    icgMin = MINVAL(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot))
    icgMax = MAXVAL(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot))	       

    DO iPatch = 1,pGridSerial%nPatches
      pPatchSerial => pRegionSerial%patches(iPatch)
      
! ------------------------------------------------------------------------------
!     For patches which are not virtual, search for cells adjacent to serial 
!     boundary faces in partitioned region. NOTE avoid searching on virtual 
!     patches because for 2d cases, searching on virtual patches makes this 
!     procedure slow (number of virtual faces equals number of cells in serial
!     region; searching both virtual patches actually unnecessary because the 
!     two patches are necessarily partitioned exactly the same). NOTE that the 
!     distinction between virtual and other patches is NOT necessary for correct
!     running of the code, it only improves performance.
! ------------------------------------------------------------------------------      
      
      IF ( pPatchSerial%bcType /= BC_VIRTUAL ) THEN 
        DO ifl = 1,pPatchSerial%nBFacesTot
          icg  = pPatchSerial%bf2c(ifl)                
          iReg = pGridSerial%sc2r(icg)

          IF ( icg >= icgMin .AND. icg <= icgMax ) THEN ! Search
            CALL BinarySearchInteger(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot), &    
                                     pGrid%nBCellsTot,icg,iLoc)
          ELSE ! No need to search
	    iLoc = ELEMENT_NOT_FOUND
	  END IF ! icg				    

          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN
            icg2 = pGrid%sbc2pc(2,iLoc)

            IF ( icg2 <= pGrid%nCells ) THEN                              
              IF ( pPatchSerial%bf2v(4,ifl) == VERT_NONE ) THEN          
                pGrid%patchCounter(1,FACE_TYPE_TRI,iPatch) = & 
                  pGrid%patchCounter(1,FACE_TYPE_TRI,iPatch) + 1
              ELSE 
                pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch) = & 
                  pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch) + 1            
              END IF ! pPatchSerial
            ELSE 
              IF ( pPatchSerial%bf2v(4,ifl) == VERT_NONE ) THEN          
                pGrid%patchCounter(2,FACE_TYPE_TRI,iPatch) = & 
                  pGrid%patchCounter(2,FACE_TYPE_TRI,iPatch) + 1
              ELSE 
                pGrid%patchCounter(2,FACE_TYPE_QUAD,iPatch) = & 
                  pGrid%patchCounter(2,FACE_TYPE_QUAD,iPatch) + 1            
              END IF ! pPatchSerial          
            END IF ! icg2                           
          ELSE 
            IF ( iReg == pRegion%iRegionGlobal ) THEN ! Cell must exist
              CALL ErrorStop(global,ERR_CELL_NOT_FOUND,__LINE__)        
            END IF ! iReg
          END IF ! iLoc
        END DO ! ifl
        
! ------------------------------------------------------------------------------
!     For virtual patches, search is unnecessary. Can set number of faces 
!     directly, see comment above.
! ------------------------------------------------------------------------------        
        
      ELSE 
        pGrid%patchCounter(1,FACE_TYPE_TRI,iPatch)  = pGrid%nPris
        pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch) = pGrid%nHexs 
          
        pGrid%patchCounter(2,FACE_TYPE_TRI,iPatch) = pGrid%nPrisTot & 
                                                   - pGrid%nPris
        pGrid%patchCounter(2,FACE_TYPE_QUAD,iPatch) = pGrid%nHexsTot & 
                                                    - pGrid%nHexs                   
      END IF ! pPatch%bcType
    END DO ! iPatch

! ==============================================================================
!   Number of patches is given by sum of patch with non-zero number of actual 
!   or virtual faces
! ==============================================================================
              
    pGrid%nPatches = 0    
          
    DO iPatch = 1,pGridSerial%nPatches
      IF ( (pGrid%patchCounter(1,FACE_TYPE_TRI ,iPatch) > 0) .OR. &
           (pGrid%patchCounter(2,FACE_TYPE_TRI ,iPatch) > 0) .OR. &
           (pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch) > 0) .OR. &            
           (pGrid%patchCounter(2,FACE_TYPE_QUAD,iPatch) > 0) ) THEN 
        pGrid%nPatches = pGrid%nPatches + 1
      END IF ! pGrid%patchCounter
    END DO ! iPatch

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining number of patches done.' 
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Allocate memory for patch data structure
! ******************************************************************************       

    ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%error           

! ******************************************************************************       
!   Loop over serial patches and set variables if that patch exists on this 
!   partition
! ******************************************************************************       

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining patch dimensions...' 
    END IF ! global%verbLevel
    
    pGrid%nPatches = 0 ! Reset counter
    
    pGrid%nBFaces    = 0 
    pGrid%nBFacesTot = 0 

    DO iPatch = 1,pGridSerial%nPatches
      pPatchSerial => pRegionSerial%patches(iPatch)

! ==============================================================================
!     Patch exists on this partition
! ==============================================================================
    
      IF ( (pGrid%patchCounter(1,FACE_TYPE_TRI ,iPatch) > 0) .OR. &
           (pGrid%patchCounter(2,FACE_TYPE_TRI ,iPatch) > 0) .OR. & 
           (pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch) > 0) .OR. &   
           (pGrid%patchCounter(2,FACE_TYPE_QUAD,iPatch) > 0) ) THEN 
        pGrid%nPatches = pGrid%nPatches + 1

        pPatch => pRegion%patches(pGrid%nPatches)

! ------------------------------------------------------------------------------
!       Set patch variables
! ------------------------------------------------------------------------------
        
        pPatch%iPatchGlobal = iPatch

        pPatch%nBTris    = pGrid%patchCounter(1,FACE_TYPE_TRI,iPatch)
        pPatch%nBTrisTot = pGrid%patchCounter(1,FACE_TYPE_TRI,iPatch) &
                         + pGrid%patchCounter(2,FACE_TYPE_TRI,iPatch)   

        pPatch%nBQuads    = pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch)
        pPatch%nBQuadsTot = pGrid%patchCounter(1,FACE_TYPE_QUAD,iPatch) &
                          + pGrid%patchCounter(2,FACE_TYPE_QUAD,iPatch)   
                   
        pPatch%nBVert    = 0 
        pPatch%nBVertTot = 0            
                   
        pPatch%nBFaces    = pPatch%nBTris    + pPatch%nBQuads
        pPatch%nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot        

        pGrid%nBFaces     = pGrid%nBFaces    + pPatch%nBFaces
        pGrid%nBFacesTot  = pGrid%nBFacesTot + pPatch%nBFacesTot 

        pPatch%nBCellsVirt = 0  

        pPatch%bcCoupled    = pPatchSerial%bcCoupled
        pPatch%movePatchDir = pPatchSerial%movePatchDir
        pPatch%flatFlag     = pPatchSerial%flatFlag

! ------------------------------------------------------------------------------
!       Allocate memory for face lists
! ------------------------------------------------------------------------------

        IF ( pPatch%nBTrisTot > 0 ) THEN 
          ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisTot),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
          END IF ! global%error   
        ELSE 
          NULLIFY(pPatch%bTri2v)
        END IF ! pPatch%nBTrisTot

        IF ( pPatch%nBQuadsTot > 0 ) THEN 
          ALLOCATE(pPatch%bQuad2v(4,pPatch%nBQuadsTot),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bQuad2v')
          END IF ! global%error   
        ELSE 
          NULLIFY(pPatch%bQuad2v)
        END IF ! pPatch%nBQuadTot                                                                       
      END IF ! pGrid%patchCounter                  
    END DO ! iPatch 

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining patch dimensions done.' 
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Deallocate temporary memory 
! ******************************************************************************       

    DEALLOCATE(pGrid%patchCounter,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%patchCounter')
    END IF ! global%error                    
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating patch lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_CreatePatchLists










! ******************************************************************************
!
! Purpose: Create region-to-cell mapping.
!
! Description: None.
!
! Input:
!   levels	Level data structure
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_CreateReg2CellMap(pRegion)

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

    INTEGER :: errorFlag,iReg  
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_CreateReg2CellMap',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating region-to-cell mapping...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    ALLOCATE(pGrid%r2pcCSR(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%r2pcCSR')
    END IF ! global%error         

    ALLOCATE(pGrid%r2pcCSRInfo(0:global%nRegionsLocal),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%r2pcCSRInfo')
    END IF ! global%error
    
    DO iReg = 0,global%nRegionsLocal
      pGrid%r2pcCSRInfo(iReg) = 0 ! Initial value important
    END DO ! iReg
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating region-to-cell mapping done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_CreateReg2CellMap








! ******************************************************************************
!
! Purpose: Destroy border face lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_DestroyBorderFaceList(pRegion)                                   

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

    CALL RegisterFunction(global,'RFLU_PART_DestroyBorderFaceList',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying border face lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
      
    DEALLOCATE(pGrid%avfCSRInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%avfCSRInfo')
    END IF ! global%error
        
    DEALLOCATE(pGrid%avfCSR,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%avfCSR')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying border face lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_DestroyBorderFaceList







! ******************************************************************************
!
! Purpose: Destroy cell data.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_DestroyCellData(pRegion)

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
    TYPE(t_grid), POINTER :: pGrid
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_DestroyCellData',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying cell data...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid
 
! ******************************************************************************
!   Deallocate memory for cell data
! ******************************************************************************

! ==============================================================================
!   Mixture
! ==============================================================================
   
    DEALLOCATE(pRegion%mixt%cv,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cv')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying cell data done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_DestroyCellData








! ******************************************************************************
!
! Purpose: Destroy cell lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_DestroyCellLists(pRegion)

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

    CALL RegisterFunction(global,'RFLU_PART_DestroyCellLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying cell lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid
                    
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    IF ( pGrid%nTetsMax > 0 ) THEN 
      DEALLOCATE(pGrid%tet2v,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%tet2v')
      END IF ! global%error      
    END IF ! pGrid%nTetsMax

    IF ( pGrid%nHexsMax > 0 ) THEN 
      DEALLOCATE(pGrid%hex2v,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%hex2v')
      END IF ! global%error
    END IF ! pGrid%nHexsMax
    
    IF ( pGrid%nPrisMax > 0 ) THEN 
      DEALLOCATE(pGrid%pri2v,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pri2v')
      END IF ! global%error
    END IF ! pGrid%nPrisMax    
 
    IF ( pGrid%nPyrsMax > 0 ) THEN 
      DEALLOCATE(pGrid%pyr2v,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pyr2v')
      END IF ! global%error
    END IF ! pGrid%nPyrsMax 

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying cell lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_DestroyCellLists








! ******************************************************************************
!
! Purpose: Destroy patch list.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_DestroyPatchLists(pRegion)

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

    INTEGER :: errorFlag,iPatch
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch    
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_DestroyPatchLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying patch lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid         

! ******************************************************************************       
!   Deallocate memory for face lists
! ******************************************************************************       

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%nBTrisTot > 0 ) THEN 
        DEALLOCATE(pPatch%bTri2v,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bTri2v')
        END IF ! global%error   
      END IF ! pPatch%nBTrisTot

      IF ( pPatch%nBQuadsTot > 0 ) THEN 
        DEALLOCATE(pPatch%bQuad2v,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bQuad2v')
        END IF ! global%error   
      END IF ! pPatch%nBQuadTot                  
    END DO ! iPatch 
    
! ******************************************************************************
!   Deallocate memory for patch lists
! ******************************************************************************       
         
    DEALLOCATE(pRegion%patches,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%error      
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying patch lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_DestroyPatchLists








! ******************************************************************************
!
! Purpose: Destroy vertex data.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_DestroyVertexData(pRegion)

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
    TYPE(t_grid), POINTER :: pGrid
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_DestroyVertexData',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying vertex data...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    DEALLOCATE(pGrid%xyz,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%xyz')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying vertex data done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_DestroyVertexData










! ******************************************************************************
!
! Purpose: Partition region.
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

  SUBROUTINE RFLU_PART_PartitionRegion(pRegion)

    USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_HaveSyPePatches

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

    INTEGER :: c1,c2,errorFlag,icg,icgBeg,icgEnd,ifg,ifl,iPatch,iReg, & 
               nCellsPerReg,nCellsV2,nFaces,wgtFlag
    INTEGER, DIMENSION(5) :: options
    INTEGER, DIMENSION(:), ALLOCATABLE :: f2cCSR,f2cCSRInfo,vwgt,adjwgt
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_PartitionRegion',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Partitioning region...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Convert interior face list to CSR format
! ******************************************************************************    

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Converting face list to CSR format...' 
    END IF ! global%verbLevel

! ==============================================================================
!   Allocate memory. NOTE need to include virtual-virtual faces so that cases
!   with virtual cells arising from symmetry and periodic patches can be 
!   partitioned properly. 
! ==============================================================================    
         
    nFaces = pGrid%nFaces + pGrid%nFacesVV     
         
    ALLOCATE(f2cCSR(2*nFaces),STAT=errorFlag)    
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2cCSR')
    END IF ! global%error

    ALLOCATE(f2cCSRInfo(pGrid%nCellsTot+1),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2cCSRInfo')
    END IF ! global%error
        
    DO icg = 1,pGrid%nCellsTot+1   
      f2cCSRInfo(icg) = 0 ! Initial value important
    END DO ! icg   
    
! ==============================================================================
!   Build lists 
! ==============================================================================    
    
! ------------------------------------------------------------------------------
!   Compute degree of each cell and sum up to get ending index of range of cell
!   neighbors in CSR list
! ------------------------------------------------------------------------------    
    
    DO ifg = 1,pGrid%nFaces + pGrid%nFacesVV
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)
    
      f2cCSRInfo(c1) = f2cCSRInfo(c1) + 1
      f2cCSRInfo(c2) = f2cCSRInfo(c2) + 1
    END DO ! ifg
    
    f2cCSRInfo(1) = f2cCSRInfo(1) + 1
    
    DO icg = 2,pGrid%nCellsTot
      f2cCSRInfo(icg) = f2cCSRInfo(icg) + f2cCSRInfo(icg-1)  
    END DO ! icg

! ------------------------------------------------------------------------------
!   Build CSR list
! ------------------------------------------------------------------------------    
    
    DO ifg = 1,pGrid%nFaces + pGrid%nFacesVV
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

      f2cCSRInfo(c1) = f2cCSRInfo(c1) - 1
      f2cCSRInfo(c2) = f2cCSRInfo(c2) - 1      
     
      f2cCSR(f2cCSRInfo(c1)) = c2 
      f2cCSR(f2cCSRInfo(c2)) = c1   
    END DO ! ifg    
    
    f2cCSRInfo(pGrid%nCellsTot+1) = 2*nFaces + 1 
   
    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Converting face list to CSR format done.' 
    END IF ! global%verbLevel   
   
! ******************************************************************************
!   Partition region
! ******************************************************************************    
        
    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Calling partitioner...' 
    END IF ! global%verbLevel   
       
! ==============================================================================
!   Partition using METIS. NOTE distinction between the two METIS calls coded
!   as recommended in METIS manual.
! ==============================================================================

    IF ( global%prepPartMode == PARTITION_MODE_PROPER ) THEN
      ALLOCATE(vwgt(1),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
	CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vwgt')
      END IF ! global%error

      ALLOCATE(adjwgt(1),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
	CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'adjwgt')
      END IF ! global%error    

      wgtFlag    = 0 ! No weights on graph
      options(1) = 0 ! Use default settings

      IF ( global%nRegionsLocal < 8 ) THEN 
	CALL METIS_PartGraphRecursive(pGrid%nCellsTot,f2cCSRInfo,f2cCSR,vwgt, &
                                      adjwgt,wgtFlag,1,global%nRegionsLocal, &
                                      options,pGrid%nFacesCut,pGrid%sc2r)
      ELSE 
	CALL METIS_PartGraphKway(pGrid%nCellsTot,f2cCSRInfo,f2cCSR,vwgt, & 
                                 adjwgt,wgtFlag,1,global%nRegionsLocal, & 
                                 options,pGrid%nFacesCut,pGrid%sc2r)      
      END IF ! global%nRegionsLocal

      DEALLOCATE(vwgt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
	CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vwgt')
      END IF ! global%error   

      DEALLOCATE(adjwgt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
	CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'adjwgt')
      END IF ! global%error

! ==============================================================================
!   Impose partition mapping. NOTE need to use only pGrid%nCells to compute
!   nCellsPerReg to avoid problems with virtual cells arising from periodic
!   or symmetry boundaries, otherwise get imbalanced regions. NOTE also need
!   to apportion virtual cells arising from periodic or symmetry boundaries 
!   to first and last regions.
! ==============================================================================

    ELSE IF ( global%prepPartMode == PARTITION_MODE_IMPOSED ) THEN
    
! ------------------------------------------------------------------------------
!     Basic imposed mapping
! ------------------------------------------------------------------------------    
    
      nCellsPerReg = pGrid%nCells/global%nRegionsLocal
 
      DO iReg = 1,global%nRegionsLocal
	icgBeg = nCellsPerReg*(iReg - 1) + 1
	icgEnd = nCellsPerReg* iReg

	IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,5X,I4,2(1X,I9))') SOLVER_NAME,iReg,icgBeg,icgEnd
	END IF ! global%verbLevel   

        DO icg = icgBeg,icgEnd
	  pGrid%sc2r(icg) = iReg
        END DO ! icg
      END DO ! iReg

! ------------------------------------------------------------------------------
!     If have virtual cells in serial region, must be due to periodic or 
!     symmetry boundaries, and hence need to take these into account in 
!     special manner. NOTE code below will only work for periodic boundaries
!     and assumes that virtual cells in serial region are added at opposite
!     ends so can be assigned wholly to the first and last regions without
!     being partitioned themselves.
! ------------------------------------------------------------------------------    
      
      IF ( pGrid%nCells /= pGrid%nCellsTot ) THEN 
        IF ( RFLU_SYPE_HaveSyPePatches(pRegion) .EQV. .TRUE. ) THEN
          IF ( MOD(pGrid%nCellsTot-pGrid%nCells,2) == 0 ) THEN
            nCellsV2 = (pGrid%nCellsTot-pGrid%nCells)/2
           
            DO icg = pGrid%nCells+1,pGrid%nCells+nCellsV2
              pGrid%sc2r(icg) = 1
            END DO ! icg
            
            DO icg = pGrid%nCells+nCellsV2+1,pGrid%nCellsTot
              pGrid%sc2r(icg) = global%nRegionsLocal
            END DO ! icg            
          ELSE 
            CALL ErrorStop(global,ERR_VIRTUALCELLS_NOTDB2,__LINE__)
          END IF ! MOD
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! 
      END IF ! pGrid%nCells  

      pGrid%nFacesCut = 0

      DO ifg = 1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        IF ( pGrid%sc2r(c1) /= pGrid%sc2r(c2) ) THEN 
          pGrid%nFacesCut = pGrid%nFacesCut + 1
        END IF ! pGrid%sc2r
      END DO ! ifg    

! ==============================================================================
!   Default
! ==============================================================================

    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%prepPartMode

    IF ( global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Number of cut faces:', &
                                     pGrid%nFacesCut
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Calling partitioner done.' 
    END IF ! global%verbLevel  
       
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************    
       
    DEALLOCATE(f2cCSRInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2cCSRInfo')
    END IF ! global%error
                
    DEALLOCATE(f2cCSR,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2cCSR')
    END IF ! global%error
            
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Partitioning region done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_PartitionRegion








! ******************************************************************************
!
! Purpose: Recreate cell list.
!
! Description: None.
!
! Input:
!   global		Pointer to global data
!   nVertPerCell	Number of vertices per cell
!   nCellsMax		Maximum number of cells
!   x2v			Connectivity array
!   x2cg		Cell mapping array
!
! Output: 
!   nCellsMax		Increased maximum number of cells
!   x2v			Enlarged connectivity array
!   x2cg		Enlarged cell mapping array
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_RecreateCellList(global,nVertPerCell,nCellsMax,x2v,x2cg)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nVertPerCell
    INTEGER, INTENT(INOUT) :: nCellsMax
    INTEGER, DIMENSION(:), POINTER :: x2cg
    INTEGER, DIMENSION(:,:), POINTER :: x2v
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icl,ivl,nCellsMaxOld    
    INTEGER, DIMENSION(:), ALLOCATABLE:: x2cgTemp
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: x2vTemp      

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_PART_RecreateCellList',&
  'RFLU_ModPartitionRegion.F90')

! ******************************************************************************
!   Increase maximum number of cells
! ******************************************************************************

    nCellsMaxOld =   nCellsMax 
    nCellsMax    = 2*nCellsMax

! ******************************************************************************
!   Copy existing arrays into larger arrays
! ******************************************************************************

! ==============================================================================
!   Connectivity array
! ==============================================================================

    ALLOCATE(x2vTemp(nVertPerCell,nCellsMaxOld),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2vTemp')
    END IF ! global%error
          
    DO icl = 1,nCellsMaxOld
      DO ivl = 1,nVertPerCell
        x2vTemp(ivl,icl) = x2v(ivl,icl)
      END DO ! ivl
    END DO ! icl

    DEALLOCATE(x2v,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2v')
    END IF ! global%error

    ALLOCATE(x2v(nVertPerCell,nCellsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2v')
    END IF ! global%error

    DO icl = 1,nCellsMaxOld
      DO ivl = 1,nVertPerCell
        x2v(ivl,icl) = x2vTemp(ivl,icl)
      END DO ! ivl
    END DO ! icl    

    DEALLOCATE(x2vTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2vTemp')
    END IF ! global%error

! ==============================================================================
!   Cell mapping array
! ==============================================================================

    ALLOCATE(x2cgTemp(nCellsMaxOld),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2cgTemp')
    END IF ! global%error
          
    DO icl = 1,nCellsMaxOld
      x2cgTemp(icl) = x2cg(icl)
    END DO ! icl

    DEALLOCATE(x2cg,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2cg')
    END IF ! global%error

    ALLOCATE(x2cg(nCellsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2cg')
    END IF ! global%error

    DO icl = 1,nCellsMaxOld
      x2cg(icl) = x2cgTemp(icl)
    END DO ! icl    

    DEALLOCATE(x2cgTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2cgTemp')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Renumbering vertex lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_RecreateCellList









! ******************************************************************************
!
! Purpose: Renumber vertex lists.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_RenumberVertexLists(pRegion)

    USE RFLU_ModRenumberList, ONLY: RFLU_RenumberList2                                                         

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

    INTEGER :: errorFlag,icl,iPatch
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_RenumberVertexLists',&
  'RFLU_ModPartitionRegion.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Renumbering vertex lists...' 
    END IF ! global%verbLevel

    IF ( global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid
        
! ******************************************************************************
!   Renumber volume connectivity lists
! ******************************************************************************
    
! ==============================================================================
!   Tetrahedra
! ==============================================================================    
    
    IF ( pGrid%nTetsTot > 0 ) THEN 
      CALL RFLU_RenumberList2(global,4,pGrid%nTetsTot, & 
                              pGrid%tet2v(1:4,1:pGrid%nTetsTot), & 
                              pGrid%nVertTot, &
                              pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                              pGrid%sv2pv(2:2,1:pGrid%nVertTot))
    END IF ! pGrid%nTetsTot
    
! ==============================================================================
!   Hexahedra
! ==============================================================================    
    
    IF ( pGrid%nHexsTot > 0 ) THEN 
      CALL RFLU_RenumberList2(global,8,pGrid%nHexsTot, & 
                              pGrid%hex2v(1:8,1:pGrid%nHexsTot), & 
                              pGrid%nVertTot, &
                              pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                              pGrid%sv2pv(2:2,1:pGrid%nVertTot))   
    END IF ! pGrid%nHexsTot  
          
! ==============================================================================
!   Prisms
! ==============================================================================    

    IF ( pGrid%nPrisTot > 0 ) THEN     
      CALL RFLU_RenumberList2(global,6,pGrid%nPrisTot, & 
                              pGrid%pri2v(1:6,1:pGrid%nPrisTot), & 
                              pGrid%nVertTot, &
                              pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                              pGrid%sv2pv(2:2,1:pGrid%nVertTot))        
    END IF ! pGrid%nPrisTot      
          
! ==============================================================================
!   Pyramids
! ==============================================================================    

    IF ( pGrid%nPyrsTot > 0 ) THEN     
      CALL RFLU_RenumberList2(global,5,pGrid%nPyrsTot, & 
                              pGrid%pyr2v(1:5,1:pGrid%nPyrsTot), & 
                              pGrid%nVertTot, &
                              pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                              pGrid%sv2pv(2:2,1:pGrid%nVertTot))
    END IF ! pGrid%nPyrsTot

! ******************************************************************************
!   Renumber surface connectivity lists
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( pPatch%nBTrisTot > 0 ) THEN                        
        CALL RFLU_RenumberList2(global,3,pPatch%nBTrisTot, & 
                                pPatch%bTri2v(1:3,1:pPatch%nBTrisTot), & 
                                pGrid%nVertTot, &                                                                                     
                                pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                                pGrid%sv2pv(2:2,1:pGrid%nVertTot))
      END IF ! pPatch%nBTrisTot
        
      IF ( pPatch%nBQuadsTot > 0 ) THEN                        
        CALL RFLU_RenumberList2(global,4,pPatch%nBQuadsTot, & 
                                pPatch%bQuad2v(1:4,1:pPatch%nBQuadsTot), & 
                                pGrid%nVertTot, &
                                pGrid%sv2pv(1:1,1:pGrid%nVertTot), &
                                pGrid%sv2pv(2:2,1:pGrid%nVertTot))                     
      END IF ! pPatch%nBQuadsTot
    END DO ! iPatch
          
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Renumbering vertex lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_RenumberVertexLists









! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModPartitionRegion


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModPartitionRegion.F90,v $
! Revision 1.17  2008/12/06 08:45:03  mtcampbe
! Updated license.
!
! Revision 1.16  2008/11/19 22:18:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.15  2006/08/21 19:15:25  haselbac
! Bug fixes for newly added modifications
!
! Revision 1.14  2006/08/21 16:45:24  haselbac
! Bug fix or extension: Balanced partitioning with periodic boundaries
!
! Revision 1.13  2006/04/12 16:11:32  haselbac
! Bug fix in building patch lists, did not loop over nBFacesTot
!
! Revision 1.12  2006/03/25 22:05:48  haselbac
! Substantial changes because of sype patches
!
! Revision 1.11  2005/12/03 19:36:40  haselbac
! Bug fix: Removed hardcoded IF statement on pGrid%nTetsTot
!
! Revision 1.10  2005/08/05 15:28:28  haselbac
! Improved speed of routines for creating and building patch data str
!
! Revision 1.9  2005/07/06 15:55:38  haselbac
! Added imposed partitioning mode
!
! Revision 1.8  2005/07/03 02:42:51  haselbac
! Bug fix: Removed INTENT dataitem from pointer declarations
!
! Revision 1.7  2005/07/01 16:17:59  haselbac
! Changed adding of virtual cells so will not exceed max dims anymore
!
! Revision 1.6  2005/06/13 22:44:09  haselbac
! Bug fix: Initialize movePatchDir
!
! Revision 1.5  2005/05/05 01:44:57  haselbac
! Increased max number of cells for small cases
!
! Revision 1.4  2005/05/04 03:37:18  haselbac
! Bug fix: Added setting of bcCoupled when creating patches
!
! Revision 1.3  2005/04/21 01:39:24  haselbac
! Modified building of vertex lists, creation and building of patch lists
!
! Revision 1.2  2005/04/20 14:47:30  haselbac
! Changed setting of max cell dimensions
!
! Revision 1.1  2005/04/15 15:09:12  haselbac
! Initial revision
!
! Revision 1.6  2005/01/20 14:56:19  haselbac
! Some clean-up, adapted to RNMB changes, added use of sbc2pc mapping
!
! Revision 1.5  2005/01/17 19:49:11  haselbac
! Proper specification of virtual cells, bug fix, clean-up
!
! Revision 1.4  2004/12/29 21:13:11  haselbac
! Substantial changes, creations of virtual cells and bface renumb
!
! Revision 1.3  2004/12/04 03:41:05  haselbac
! Substantial rewrite and expansion
!
! Revision 1.2  2004/11/11 15:11:24  haselbac
! Commented out METIS calls so as not to break compilation on some machines
!
! Revision 1.1  2004/11/08 23:27:01  haselbac
! Initial revision, work in progress
!
! ******************************************************************************
























