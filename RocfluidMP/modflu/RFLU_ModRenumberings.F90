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
! Purpose: Suite of routines to create, write, read, and destroy renumbering 
!   maps.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModRenumberings.F90,v 1.16 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRenumberings

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE ModSortSearch

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_RNMB_BuildSBC2PCMap, & 
            RFLU_RNMB_BuildSC2PCMap, & 
            RFLU_RNMB_BuildSV2PVMap, &  
            RFLU_RNMB_CreatePBF2SBFMap, &             
            RFLU_RNMB_CreatePC2SCMap, & 
            RFLU_RNMB_CreatePV2SVMap, &
            RFLU_RNMB_CreateSC2RMap, & 
            RFLU_RNMB_DestroySBC2PCMap, &            
            RFLU_RNMB_DestroyPC2SCMap, &
            RFLU_RNMB_DestroyPBF2SBFMap, &            
            RFLU_RNMB_DestroyPV2SVMap, &
            RFLU_RNMB_DestroySC2PCMap, &    
            RFLU_RNMB_DestroySC2RMap, &  
            RFLU_RNMB_DestroySV2PVMap, &
            RFLU_RNMB_ReadSC2RMap, &            
            RFLU_RNMB_ReadPxx2SxxMaps, & 
            RFLU_RNMB_WriteSC2RMap, &               
            RFLU_RNMB_WritePxx2SxxMaps
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
        
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModRenumberings.F90,v $ $Revision: 1.16 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Build cell mapping from serial global boundary cell index to 
!   partitioned global cell index.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pRegionSerial       Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_RNMB_BuildSBC2PCMap(pRegion,pRegionSerial)                                   

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

    TYPE(t_global), POINTER :: global     
            
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_BuildSBC2PCMap',&
  'RFLU_ModRenumberings.F90')

! ******************************************************************************
!   Call appropriate routine depending on dimensionality
! ******************************************************************************

    SELECT CASE ( pRegionSerial%mixtInput%dimens ) 
      CASE ( 1,2 ) 
        CALL RFLU_RNMB_BuildSBC2PCMap2D(pRegion)
      CASE ( 3 ) 
        CALL RFLU_RNMB_BuildSBC2PCMap3D(pRegion,pRegionSerial)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegionSerial%mixtInput%dimens
                                                                
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_BuildSBC2PCMap







! ******************************************************************************
!
! Purpose: Build cell mapping from serial global boundary cell index to 
!   partitioned global cell index for 2d grids.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes:
!   1. For 2d grids, every cell is on at least one boundary, so searching like 
!      for 3d case is superfluous. 
!
! ******************************************************************************

  SUBROUTINE RFLU_RNMB_BuildSBC2PCMap2D(pRegion)                                   

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

    INTEGER :: errorFlag,icg,icg2
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global     
            
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_BuildSBC2PCMap2D',&
  'RFLU_ModRenumberings.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sbc2pc mapping...' 
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
!   Allocate memory
! ******************************************************************************

    pGrid%nBCellsTot = pGrid%nCellsTot

    ALLOCATE(pGrid%sbc2pc(2,pGrid%nBCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%sbc2pc')
    END IF ! global%error           
                
! ******************************************************************************
!   Build list 
! ******************************************************************************

    DO icg = 1,pGrid%nCellsTot
      icg2 = pGrid%pc2sc(icg)

      pGrid%sbc2pc(1,icg) = icg2
      pGrid%sbc2pc(2,icg) = icg            
    END DO ! icg

! ******************************************************************************
!   Sort list
! ******************************************************************************

    CALL QuickSortIntegerInteger(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot), & 
                                 pGrid%sbc2pc(2:2,1:pGrid%nBCellsTot), & 
                                 pGrid%nBCellsTot)
                                                                
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sbc2pc mapping done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_BuildSBC2PCMap2D










! ******************************************************************************
!
! Purpose: Build cell mapping from serial global boundary cell index to 
!   partitioned global cell index for 3d grids.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pRegionSerial       Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_RNMB_BuildSBC2PCMap3D(pRegion,pRegionSerial)                                   

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

    INTEGER :: errorFlag,iBcmMax,iBcmMaxLoc,iBcmMin,iBcmMinLoc,icg,icg2,icl, & 
               iLoc,iLocMax,iLocMin,iPc2scMax,iPc2scMaxLoc,iPc2scMin, &
               iPc2scMinLoc               
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: pc2scSorted,sbc2pc
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
    TYPE(t_global), POINTER :: global     
            
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_BuildSBC2PCMap3D',&
  'RFLU_ModRenumberings.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sbc2pc mapping...' 
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

! TO DO 
!   Better estimate needed for pGrid%nBCellsTot instead of pGrid%nCellsTot
! END TO DO 

    ALLOCATE(sbc2pc(2,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'sbc2pc')
    END IF ! global%error           
                
! ******************************************************************************
!   Build temporary list with boundary cells. NOTE cannot loop over boundary 
!   cells on serial grid and check whether in given region because will not 
!   capture virtual cells that way.
! ******************************************************************************

    pGrid%nBCellsTot = 0 

! ==============================================================================
!   Find extrema of entries and their locations so that search can be 
!   restricted.
! ==============================================================================

    iPc2scMin = MINVAL(pGrid%pc2sc(1:pGrid%nCellsTot))
    iPc2scMax = MAXVAL(pGrid%pc2sc(1:pGrid%nCellsTot))

    iBcmMin = pGridSerial%bcm(1)                   ! NOTE Because bcm sorted 
    iBcmMax = pGridSerial%bcm(pGridSerial%nBCells) ! NOTE Because bcm sorted 

! ==============================================================================
!   Search for boundary cells
! ==============================================================================
          
! ------------------------------------------------------------------------------          
!   Range of serial boundary cells exceeds that of partitioned cells. Do not 
!   need to sort bcm array because already sorted.
! ------------------------------------------------------------------------------          
                        
    IF ( iBcmMax - iBcmMin >= iPc2scMax - iPc2scMin ) THEN 
      CALL BinarySearchInteger(pGridSerial%bcm(1:pGridSerial%nBCells), & 
                               pGridSerial%nBCells,iPc2scMin,iLoc,iLocMin)

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        iPc2scMinLoc = iLoc  
      ELSE
        iPc2scMinLoc = iLocMin
      END IF ! iLoc                                                          

      CALL BinarySearchInteger(pGridSerial%bcm(1:pGridSerial%nBCells), & 
                               pGridSerial%nBCells,iPc2scMax,iLoc,iLocMax)

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        iPc2scMaxLoc = iLoc  
      ELSE
        iPc2scMaxLoc = MIN(iLocMax,pGridSerial%nBCells)
      END IF ! iLoc     

      DO icg = 1,pGrid%nCellsTot
        icg2 = pGrid%pc2sc(icg)

        CALL BinarySearchInteger(pGridSerial%bcm(iPc2scMinLoc:iPc2scMaxLoc), & 
                                 iPc2scMaxLoc-iPc2scMinLoc+1,icg2,iLoc)

        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
          pGrid%nBCellsTot = pGrid%nBCellsTot + 1

          sbc2pc(1,pGrid%nBCellsTot) = icg2
          sbc2pc(2,pGrid%nBCellsTot) = icg         
        END IF ! iLoc                            
      END DO ! icg

! ------------------------------------------------------------------------------
!   Range of partitioned cells exceeds that of serial boundary cells. Must sort
!   partitioned cells first to get range which lies within range of serial 
!   boundary cells.
! ------------------------------------------------------------------------------          
                
    ELSE 
      ALLOCATE(pc2scSorted(2,pGrid%nCellsTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pc2scSorted')
      END IF ! global%error
      
      DO icg = 1,pGrid%nCellsTot
        pc2scSorted(1,icg) = pGrid%pc2sc(icg)
        pc2scSorted(2,icg) = icg
      END DO ! icg
      
      CALL QuickSortIntegerInteger(pc2scSorted(1,1:pGrid%nCellsTot), & 
                                   pc2scSorted(2,1:pGrid%nCellsTot), & 
                                   pGrid%nCellsTot)
    
      CALL BinarySearchInteger(pc2scSorted(1,1:pGrid%nCellsTot), & 
                               pGrid%nCellsTot,iBcmMin,iLoc,iLocMin)

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        iBcmMinLoc = iLoc  
      ELSE
        iBcmMinLoc = iLocMin
      END IF ! iLoc                                                          

      CALL BinarySearchInteger(pc2scSorted(1,1:pGrid%nCellsTot), & 
                               pGrid%nCellsTot,iBcmMax,iLoc,iLocMax)

      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        iBcmMaxLoc = iLoc  
      ELSE
        iBcmMaxLoc = MIN(iLocMax,pGrid%nCellsTot) 
      END IF ! iLoc               
            
      DO icl = iBcmMinLoc,iBcmMaxLoc
        icg2 = pc2scSorted(1,icl)
        icg  = pc2scSorted(2,icl)

        CALL BinarySearchInteger(pGridSerial%bcm(1:pGridSerial%nBCells), & 
                                 pGridSerial%nBCells,icg2,iLoc)

        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
          pGrid%nBCellsTot = pGrid%nBCellsTot + 1

          sbc2pc(1,pGrid%nBCellsTot) = icg2
          sbc2pc(2,pGrid%nBCellsTot) = icg         
        END IF ! iLoc                            
      END DO ! icl      
    
      DEALLOCATE(pc2scSorted,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pc2scSorted')
      END IF ! global%error
    END IF ! iBcmMax                              

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Number of boundary cells:', &
                                     pGrid%nBCellsTot
    END IF ! global%verbLevel

! ******************************************************************************
!   Copy and deallocate temporary memory
! ******************************************************************************

    ALLOCATE(pGrid%sbc2pc(2,pGrid%nBCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'sbc2pc')
    END IF ! global%error     

    DO icg = 1,pGrid%nBCellsTot
      pGrid%sbc2pc(1,icg) = sbc2pc(1,icg)
      pGrid%sbc2pc(2,icg) = sbc2pc(2,icg)
    END DO ! icl

    DEALLOCATE(sbc2pc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'sbc2pc')
    END IF ! global%error  

    IF ( pGrid%nBCellsTot > 0 ) THEN
      CALL QuickSortIntegerInteger(pGrid%sbc2pc(1:1,1:pGrid%nBCellsTot), & 
                                   pGrid%sbc2pc(2:2,1:pGrid%nBCellsTot), & 
                                   pGrid%nBCellsTot)
    END IF ! pGrid%nBCellsTot
                                                                
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sbc2pc mapping done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_BuildSBC2PCMap3D










! ******************************************************************************
!
! Purpose: Build cell mapping from serial global cell index to partitioned 
!   global cell index.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   sortFlag		Sorting flag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_RNMB_BuildSC2PCMap(pRegion,sortFlag)                                   

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    LOGICAL, OPTIONAL, INTENT(IN) :: sortFlag
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_BuildSC2PCMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sc2pc mapping...' 
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
!   Allocate and build array
! ******************************************************************************
    
    ALLOCATE(pGrid%sc2pc(2,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%sc2pc')
    END IF ! global%error           
      
    DO icg = 1,pGrid%nCellsTot
      pGrid%sc2pc(1,icg) = pGrid%pc2sc(icg)
      pGrid%sc2pc(2,icg) = icg
    END DO ! icg

! ******************************************************************************
!   Sort array entries so can search. NOTE that regardless of sortFlag 
!   presence and value, always sort sc2pc. If sortFlag is present and FALSE,
!   then sort actual cell entries only.
! ******************************************************************************
      
    IF ( PRESENT(sortFlag) ) THEN
      IF ( sortFlag .EQV. .TRUE. ) THEN 
        CALL QuickSortIntegerInteger(pGrid%sc2pc(1:1,1:pGrid%nCellsTot), & 
                                     pGrid%sc2pc(2:2,1:pGrid%nCellsTot), & 
                                     pGrid%nCellsTot)        
      ELSE 
        CALL QuickSortIntegerInteger(pGrid%sc2pc(1:1,1:pGrid%nCells), & 
                                     pGrid%sc2pc(2:2,1:pGrid%nCells), & 
                                     pGrid%nCells)         
      END IF ! sortFlag 
    ELSE    
      CALL QuickSortIntegerInteger(pGrid%sc2pc(1:1,1:pGrid%nCellsTot), & 
                                   pGrid%sc2pc(2:2,1:pGrid%nCellsTot), & 
                                   pGrid%nCellsTot)
    END IF ! PRESENT(sortFlag)
                                                    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sc2pc mapping done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_BuildSC2PCMap







! ******************************************************************************
!
! Purpose: Build vertex mapping from serial vertex index to partitioned 
!   vertex index.
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

  SUBROUTINE RFLU_RNMB_BuildSV2PVMap(pRegion)                                   

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

    INTEGER :: errorFlag,ivg
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_BuildSV2PVMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sv2pv mapping...' 
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
!   Determine number of cells for each type
! ******************************************************************************
    
    ALLOCATE(pGrid%sv2pv(2,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%sv2pv')
    END IF ! global%error           
      
    DO ivg = 1,pGrid%nVertTot
      pGrid%sv2pv(1,ivg) = pGrid%pv2sv(ivg)
      pGrid%sv2pv(2,ivg) = ivg
    END DO ! ivg
      
    CALL QuickSortIntegerInteger(pGrid%sv2pv(1:1,1:pGrid%nVertTot), & 
                                 pGrid%sv2pv(2:2,1:pGrid%nVertTot), & 
                                 pGrid%nVertTot)
                                                   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building sv2pv mapping done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_BuildSV2PVMap








! ******************************************************************************
!
! Purpose: Create cell-to-region mapping.
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

  SUBROUTINE RFLU_RNMB_CreateSC2RMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_CreateSC2RMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating sc2r mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    ALLOCATE(pGrid%sc2r(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%sc2r')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating sc2r mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_CreateSC2RMap









! ******************************************************************************
!
! Purpose: Create partitioned-cell to serial-cell mapping.
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

  SUBROUTINE RFLU_RNMB_CreatePC2SCMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_CreatePC2SCMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pc2sc mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    ALLOCATE(pGrid%pc2sc(pGrid%nCellsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pc2sc')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pc2sc mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_CreatePC2SCMap









! ******************************************************************************
!
! Purpose: Create partitioned boundary-face to serial boundary-face mapping.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Build CSR access array here for convenience. 
!   2. Can only build CSR access array if patch data structure exists. Need 
!      this check because when building communication lists, grid and hence 
!      patches do not exist, but need to be able to store boundary-face mapping.
!      So this routine still needs to be called, but cannot build CSR access
!      array.
!
! ******************************************************************************

  SUBROUTINE RFLU_RNMB_CreatePBF2SBFMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_CreatePBF2SBFMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pbf2sbf mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for boundary face mapping
! ******************************************************************************
    
    IF ( pGrid%nBFacesTot > 0 ) THEN              
      ALLOCATE(pGrid%pbf2sbfCSR(pGrid%nBFacesTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pbf2sbfCSR')
      END IF ! global%error          
      
      IF ( ASSOCIATED(pRegion%patches) .EQV. .TRUE. ) THEN 
        ALLOCATE(pGrid%pbf2sbfCSRInfo(pGrid%nPatches),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pbf2sbfCSRInfo')
        END IF ! global%error                    
        
        pGrid%pbf2sbfCSRInfo(1) = 1
        
        DO iPatch = 2,pGrid%nPatches
          pPatch => pRegion%patches(iPatch-1)
          
          pGrid%pbf2sbfCSRInfo(iPatch) = pGrid%pbf2sbfCSRInfo(iPatch-1) & 
                                       + pPatch%nBFacesTot                         
        END DO ! iPatch                    
      END IF ! ASSOCIATED              
    ELSE 
      NULLIFY(pGrid%pbf2sbfCSR)
      NULLIFY(pGrid%pbf2sbfCSRInfo)         
    END IF ! pGrid%nBFacesTot        
   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pbf2sbf mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_CreatePBF2SBFMap










! ******************************************************************************
!
! Purpose: Create partitioned-vertex to serial-vertex mapping.
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

  SUBROUTINE RFLU_RNMB_CreatePV2SVMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_CreatePV2SVMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pv2sv mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    ALLOCATE(pGrid%pv2sv(pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pv2sv')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pv2sv mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_CreatePV2SVMap








! ******************************************************************************
!
! Purpose: Destroy partitioned-cell to serial-cell mapping.
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

  SUBROUTINE RFLU_RNMB_DestroyPC2SCMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_DestroyPC2SCMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pc2sc mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    DEALLOCATE(pGrid%pc2sc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pc2sc')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pc2sc mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_DestroyPC2SCMap










! ******************************************************************************
!
! Purpose: Destroy partitioned boundary-face to serial boundary-face mapping.
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

  SUBROUTINE RFLU_RNMB_DestroyPBF2SBFMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_DestroyPBF2SBFMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pbf2sbf mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for boundary face mapping
! ******************************************************************************
      
    IF ( pGrid%nPatches > 0 ) THEN        
      DEALLOCATE(pGrid%pbf2sbfCSR,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pbf2sbfCSR')
      END IF ! global%error                             
    END IF ! pGrid%nPatches

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pbf2sbf mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_DestroyPBF2SBFMap








! ******************************************************************************
!
! Purpose: Destroy partitioned-vertex to serial-vertex mapping.
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

  SUBROUTINE RFLU_RNMB_DestroyPV2SVMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_DestroyPV2SVMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pv2sv mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    DEALLOCATE(pGrid%pv2sv,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pv2sv')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pv2sv mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_DestroyPV2SVMap







! ******************************************************************************
!
! Purpose: Destroy mapping from serial global boundary cell index to 
!   partitioned global cell index.
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

  SUBROUTINE RFLU_RNMB_DestroySBC2PCMap(pRegion)                                   

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

    CALL RegisterFunction(global,'RFLU_RNMB_DestroySBC2PCMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sbc2pc mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid
                
! ******************************************************************************
!   Determine number of cells for each type
! ******************************************************************************
    
    pGrid%nBCellsTot = 0 
    
    DEALLOCATE(pGrid%sbc2pc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%sbc2pc')
    END IF ! global%error           
                                    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sbc2pc mapping done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_DestroySBC2PCMap








! ******************************************************************************
!
! Purpose: Destroy mapping from serial global cell index to partitioned global 
!   cell index.
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

  SUBROUTINE RFLU_RNMB_DestroySC2PCMap(pRegion)                                   

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

    CALL RegisterFunction(global,'RFLU_RNMB_DestroySC2PCMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sc2pc mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid
                
! ******************************************************************************
!   Determine number of cells for each type
! ******************************************************************************
    
    DEALLOCATE(pGrid%sc2pc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%sc2pc')
    END IF ! global%error           
                                    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sc2pc mapping done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_DestroySC2PCMap








! ******************************************************************************
!
! Purpose: Destroy cell-to-region mapping.
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

  SUBROUTINE RFLU_RNMB_DestroySC2RMap(pRegion)

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

    CALL RegisterFunction(global,'RFLU_RNMB_DestroySC2RMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sc2r mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory for cell mapping
! ******************************************************************************
   
    DEALLOCATE(pGrid%sc2r,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%sc2r')
    END IF ! global%error         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sc2r mapping done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_DestroySC2RMap








! ******************************************************************************
!
! Purpose: Destroy mapping from serial vertex index to partitioned vertex index.
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

  SUBROUTINE RFLU_RNMB_DestroySV2PVMap(pRegion)                                   

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

    CALL RegisterFunction(global,'RFLU_RNMB_DestroySV2PVMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sv2pv mapping...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers 
! ******************************************************************************

    pGrid => pRegion%grid
                
! ******************************************************************************
!   Determine number of cells for each type
! ******************************************************************************
    
    DEALLOCATE(pGrid%sv2pv,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%sv2pv')
    END IF ! global%error           
                                    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying sv2pv mapping done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RNMB_DestroySV2PVMap







! ******************************************************************************
!
! Purpose: Read serial cell to region map.
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

  SUBROUTINE RFLU_RNMB_ReadSC2RMap(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

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

    INTEGER :: errorFlag,icg,iFile,ivg,loopCounter,nCellsTot
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_ReadSC2RMap',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading sc2r map...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid 

! ==============================================================================
!   Open file
! ==============================================================================

    iFile = IF_RNMB

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.rnm', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# ROCFLU renumbering file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ==============================================================================
!   Dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    READ(iFile,'(I8)') nCellsTot

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================
    
    IF ( nCellsTot /= pGrid%nCellsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nCellsTot    
    
! ==============================================================================
!   Rest of file
! ==============================================================================

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1
    
      READ(iFile,'(A)') sectionString
    
      SELECT CASE ( TRIM(sectionString) ) 
    
! ------------------------------------------------------------------------------
!       Cell renumbering
! ------------------------------------------------------------------------------

        CASE ( '# Cells' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
          END IF ! global%myProcid
              
          READ(iFile,'(10(I8))') (pGrid%sc2r(icg),icg=1,pGrid%nCellsTot)    
   
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%myProcid

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(3X,A)') sectionString
          END IF ! global%myProcid

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================  

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter  
    END DO ! <empty>
   
! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading sc2r map done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RNMB_ReadSC2RMap









! ******************************************************************************
!
! Purpose: Read maps.
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

  SUBROUTINE RFLU_RNMB_ReadPxx2SxxMaps(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

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

    INTEGER :: errorFlag,icg,iFile,ifl,ivg,loopCounter,nCellsTot,nVertTot
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_ReadPxx2SxxMaps',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading Pxx2Sxx maps...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN     
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ==============================================================================
!   Open file
! ==============================================================================

    iFile = IF_RNMB

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.rnm', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# ROCFLU renumbering file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ==============================================================================
!   Dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    READ(iFile,'(2(I8))') nVertTot,nCellsTot

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================

    IF ( nVertTot /= pGrid%nVertTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nVertTot
    
    IF ( nCellsTot /= pGrid%nCellsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nCellsTot    
    
! ==============================================================================
!   Rest of file
! ==============================================================================

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1
    
      READ(iFile,'(A)') sectionString
    
      SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!       Vertex renumbering
! ------------------------------------------------------------------------------

        CASE ( '# Vertices' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Vertices...'
          END IF ! global%myProcid 
              
          READ(iFile,'(10(I8))') (pGrid%pv2sv(ivg),ivg=1,pGrid%nVertTot)
    
! ------------------------------------------------------------------------------
!       Cell renumbering
! ------------------------------------------------------------------------------

        CASE ( '# Cells' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
          END IF ! global%myProcid 
              
          READ(iFile,'(10(I8))') (pGrid%pc2sc(icg),icg=1,pGrid%nCellsTot)    

! ------------------------------------------------------------------------------
!       Boundary faces
! ------------------------------------------------------------------------------

        CASE ( '# Boundary faces' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary faces...'
          END IF ! global%myProcid 
              
          READ(iFile,'(10(I8))') (pGrid%pbf2sbfCSR(ifl),ifl=1,pGrid%nBFacesTot)
   
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%myProcid  

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(3X,A)') sectionString
          END IF ! global%myProcid      

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================  

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter  
    END DO ! <empty>
   
! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading Pxx2Sxx maps done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RNMB_ReadPxx2SxxMaps








! ******************************************************************************
!
! Purpose: Write serial cell to region map.
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

  SUBROUTINE RFLU_RNMB_WriteSC2RMap(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

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

    INTEGER :: errorFlag,icg,iFile,ivg
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_WriteSC2RMaps',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing sc2r map...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN     
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ==============================================================================
!   Open file
! ==============================================================================

    iFile = IF_RNMB

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.rnm', & 
                            pRegion%iRegionGlobal,iFileName) 

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
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    sectionString = '# ROCFLU renumbering file'  
    WRITE(iFile,'(A)') TRIM(sectionString)  

! ==============================================================================
!   Dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Dimensions...'
    END IF ! global%myProcid

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pGrid%nCellsTot
    
! ==============================================================================
!   Cell-to-region map
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
    END IF ! global%myProcid

    sectionString = '# Cells'  
    WRITE(iFile,'(A)') TRIM(sectionString)      
    WRITE(iFile,'(10(I8))') (pGrid%sc2r(icg),icg=1,pGrid%nCellsTot)    

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%myProcid

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString) 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing sc2r map done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RNMB_WriteSC2RMap








! ******************************************************************************
!
! Purpose: Write maps.
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

  SUBROUTINE RFLU_RNMB_WritePxx2SxxMaps(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

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

    INTEGER :: errorFlag,icg,iFile,ifl,ivg
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RNMB_WritePxx2SxxMaps',&
  'RFLU_ModRenumberings.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Pxx2Sxx maps...'
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid 

! ==============================================================================
!   Open file
! ==============================================================================

    iFile = IF_RNMB

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.rnm', & 
                            pRegion%iRegionGlobal,iFileName) 

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
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    sectionString = '# ROCFLU renumbering file'  
    WRITE(iFile,'(A)') TRIM(sectionString)  

! ==============================================================================
!   Dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Dimensions...'
    END IF ! global%myProcid

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(2(I8))') pGrid%nVertTot,pGrid%nCellsTot

! ==============================================================================
!   Vertex renumbering
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Vertices...'
    END IF ! global%myProcid
    
    sectionString = '# Vertices'
    WRITE(iFile,'(A)') TRIM(sectionString)     
    WRITE(iFile,'(10(I8))') (pGrid%pv2sv(ivg),ivg=1,pGrid%nVertTot)
    
! ==============================================================================
!   Cell renumbering
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
    END IF ! global%myProcid

    sectionString = '# Cells'  
    WRITE(iFile,'(A)') TRIM(sectionString)      
    WRITE(iFile,'(10(I8))') (pGrid%pc2sc(icg),icg=1,pGrid%nCellsTot)    

! ==============================================================================
!   Boundary face renumbering
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary faces...'
    END IF ! global%myProcid

    sectionString = '# Boundary faces'  
    WRITE(iFile,'(A)') TRIM(sectionString)      
    WRITE(iFile,'(10(I8))') (pGrid%pbf2sbfCSR(ifl),ifl=1,pGrid%nBFacesTot)   

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%myProcid

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString) 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Pxx2Sxx maps done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RNMB_WritePxx2SxxMaps






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModRenumberings


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRenumberings.F90,v $
! Revision 1.16  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.13  2007/03/27 00:20:02  haselbac
! Added optional argument to BuildSC2PCMap for new PLAG init
!
! Revision 1.12  2007/02/27 13:06:34  haselbac
! Enabled 1d computations
!
! Revision 1.11  2006/12/15 13:25:43  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.10  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.9  2005/08/05 15:26:48  haselbac
! Added 2d routine for building sb2pc map for efficiency reasons
!
! Revision 1.8  2005/05/17 01:12:06  haselbac
! Bug fix in building sbc2pc map: Improper upper loop index if element not found
!
! Revision 1.7  2005/05/11 13:41:18  haselbac
! Bug fix in building SBC2PC map: Only call QuickSort if have > 0 items
!
! Revision 1.6  2005/04/21 01:37:50  haselbac
! Modified building of sbc2pc to speed it up
!
! Revision 1.5  2005/04/15 15:07:02  haselbac
! Cosmetics only
!
! Revision 1.4  2005/01/20 14:51:18  haselbac
! Added sbc2pc routines, renamed routines consistently
!
! Revision 1.3  2005/01/17 19:53:54  haselbac
! Clean-up
!
! Revision 1.2  2004/12/29 21:09:16  haselbac
! Added boundary faces, changed file extension, cosmetics
!
! Revision 1.1  2004/12/04 03:44:30  haselbac
! Initial revision
!
! ******************************************************************************


























