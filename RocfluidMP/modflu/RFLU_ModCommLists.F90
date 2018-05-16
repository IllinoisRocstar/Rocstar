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
! Purpose: Suite of routines for communication lists.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModCommLists.F90,v 1.14 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModCommLists

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  USE ModBorder, ONLY: t_border

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_COMM_BuildCommLists, & 
            RFLU_COMM_CheckCountBorders, & 
            RFLU_COMM_CountBorders, &
            RFLU_COMM_CountBordersSerial, &              
            RFLU_COMM_CreateBorderCntr, & 
            RFLU_COMM_CreateBorders, &
            RFLU_COMM_CreateCommLists, &
            RFLU_COMM_DestroyBorderCntr, &  
            RFLU_COMM_DestroyBorders, & 
            RFLU_COMM_DestroyCommLists, & 
            RFLU_COMM_GetProcLocRegIds, &  
            RFLU_COMM_ReadCommLists, & 
            RFLU_COMM_WriteCommLists
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModCommLists.F90,v $ $Revision: 1.14 $' 
        
  INTEGER, PARAMETER, PUBLIC :: CREATE_BORDERS_MODE_DIM_KNOWN   = 1, & 
                                CREATE_BORDERS_MODE_DIM_UNKNOWN = 2
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  







! ******************************************************************************
!
! Purpose: Build communication lists. 
!
! Description: None.
!
! Input:
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COMM_BuildCommLists(regions)

! TEMPORARY
    USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_HaveSyPePatches
! END TEMPORARY

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
       
! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global     

! TEMPORARY
    TYPE(t_region), POINTER :: pRegionSerial    
! END TEMPORARY
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'RFLU_COMM_BuildCommLists',&
  'RFLU_ModCommLists.F90')
    
! ******************************************************************************
!   Call routines
! ******************************************************************************

    CALL RFLU_COMM_BuildCommListsCells(regions)

! TEMPORARY - At moment, vertex communication lists do not always get built
!             correctly. The problem appears to arise when the periodic patches
!             get partitioned asymmetrically. It is possible that there are 
!             multiple bugs, as the error messages for first- and second-order
!             runs are not identical. (Example case for which problems occur:
!             tubes case.)
    pRegionSerial => regions(0)

    IF ( RFLU_SYPE_HaveSyPePatches(pRegionSerial) .EQV. .FALSE. ) THEN
      CALL RFLU_COMM_BuildCommListsVert(regions)
    END IF ! RFLU_SYPE_HaveSyPePatches
! END TEMPORARY

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_COMM_BuildCommLists








! ******************************************************************************
!
! Purpose: Build cell communication lists. 
!
! Description: None.
!
! Input:
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COMM_BuildCommListsCells(regions)

    USE ModSortSearch

    USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_GetActualSerialCell

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
       
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString
    INTEGER :: errorFlag,iBorder,iBorder2, &
               icg,icgs,icgs2,icg2,icl,ict,ict2,iLoc,iReg,iReg2,iReg3,nBorders
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: borderInfo
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_grid), POINTER :: pGrid,pGrid2,pGridSerial  
    TYPE(t_global), POINTER :: global     
    TYPE(t_region), POINTER :: pRegion,pRegion2,pRegionSerial    
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'RFLU_COMM_BuildCommListsCells',&
  'RFLU_ModCommLists.F90')

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building cell communication lists...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pRegionSerial => regions(0)
    pGridSerial   => pRegionSerial%grid
    
! ******************************************************************************
!   Determine cells to be received for each region
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Determining cells to be received...'
    END IF ! global%verbLevel

    DO iReg = 1,global%nRegions
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel      

! ==============================================================================
!     Allocate temporary memory and initialize
! ==============================================================================

      ALLOCATE(borderInfo(2,pGrid%nBorders),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderInfo')
      END IF ! global%error 

      DO iBorder = 1,pGrid%nBorders
        borderInfo(1,iBorder) = pGrid%borderCntr(iBorder)
        borderInfo(2,iBorder) = 0  
      END DO ! iBorder          

! ==============================================================================
!     Loop over virtual cells and count number of cells to be received for each 
!     border
! ==============================================================================
            
      DO icg = pGrid%nCells+1,pGrid%nCellsTot
        icgs  = pGrid%pc2sc(icg)      
                
! ------------------------------------------------------------------------------
!       Get region index of cell
! ------------------------------------------------------------------------------                
                
! ----- Actual cell ------------------------------------------------------------                
                
        IF ( icgs <= pGridSerial%nCells ) THEN   
          iReg2 = pGridSerial%sc2r(icgs)

! ----- Virtual cell -----------------------------------------------------------                

        ELSE 
          CALL RFLU_SYPE_GetActualSerialCell(pRegionSerial,icgs,icgs2)
          
          iReg2 = pGridSerial%sc2r(icgs2)
        END IF ! icgs

! ------------------------------------------------------------------------------
!       Search for this region in data structure and increment counter
! ------------------------------------------------------------------------------
                   
        CALL BinarySearchInteger(borderInfo(1:1,1:pGrid%nBorders), & 
                                 pGrid%nBorders,iReg2,iLoc)
                                             
        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
          borderInfo(2,iLoc) = borderInfo(2,iLoc) + 1
        ELSE 
          CALL ErrorStop(global,ERR_REGION_ID_NOT_FOUND,__LINE__)
        END IF ! iLoc           
      END DO ! icg 

! ==============================================================================
!     Copy id of neighboring region and number of cells to be received into 
!     border data structure
! ==============================================================================

      DO iBorder = 1,pGrid%nBorders
        pGrid%borders(iBorder)%iRegionGlobal = borderInfo(1,iBorder)        
        pGrid%borders(iBorder)%nCellsRecv    = borderInfo(2,iBorder)        
        
        ALLOCATE(pGrid%borders(iBorder)%icgRecv(borderInfo(2,iBorder)), & 
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%borders%icgRecv')
        END IF ! global%error
      END DO ! iBorder

! ==============================================================================
!     Build list of cells to be received
! ==============================================================================

      DO iReg2 = 1,pGrid%nBorders
        borderInfo(2,iReg2) = 0 ! NOTE DO NOT reinitialize borderInfo(1,:)     
      END DO ! iReg2 

      DO icg = pGrid%nCells+1,pGrid%nCellsTot
        icgs  = pGrid%pc2sc(icg)      

! ----- Actual cell ------------------------------------------------------------                
                
        IF ( icgs <= pGridSerial%nCells ) THEN   
          iReg2 = pGridSerial%sc2r(icgs)

! ----- Virtual cell -----------------------------------------------------------                

        ELSE 
          CALL RFLU_SYPE_GetActualSerialCell(pRegionSerial,icgs,icgs2)
                     
          iReg2 = pGridSerial%sc2r(icgs2)
        END IF ! icgs

! ------------------------------------------------------------------------------
!       Find entry of this neighboring region in border data structure
! ------------------------------------------------------------------------------
      
        CALL BinarySearchInteger(borderInfo(1:1,1:pGrid%nBorders), & 
                                 pGrid%nBorders,iReg2,iLoc)
            
! ------------------------------------------------------------------------------
!       Increment counter and add cell
! ------------------------------------------------------------------------------

        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
          borderInfo(2,iLoc) = borderInfo(2,iLoc) + 1
          
          pGrid%borders(iLoc)%icgRecv(borderInfo(2,iLoc)) = icg
        ELSE 
          CALL ErrorStop(global,ERR_CELL_NOT_FOUND,__LINE__)                   
        END IF ! iLoc           
      END DO ! icg

! ==============================================================================
!     Deallocate temporary memory
! ==============================================================================

      DEALLOCATE(borderInfo,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderInfo')
      END IF ! global%error
    END DO ! iReg

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Determining cells to be received done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Determine cells to be sent for each region and border of other region
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining cells to be sent...'
    END IF ! global%verbLevel

    DO iReg = 1,global%nRegions
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
          
      IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel

! ==============================================================================
!     Loop over borders
! ==============================================================================          
          
      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)

        iReg2 = pBorder%iRegionGlobal
        
        pRegion2 => regions(iReg2)
        pGrid2   => pRegion2%grid
        
! ------------------------------------------------------------------------------
!       Find index of border of neighboring region whose neighboring region is
!       equal to region currently being looped over
! ------------------------------------------------------------------------------        
        
        innerLoop: DO iBorder2 = 1,pGrid2%nBorders
          IF ( pGrid2%borders(iBorder2)%iRegionGlobal == iReg ) THEN
            pBorder2 => pGrid2%borders(iBorder2)
          
! --------- Set number of cells to be sent and border id -----------------------
           
            pBorder2%iBorder    = iBorder
            pBorder2%nCellsSend = pGrid%borders(iBorder)%nCellsRecv
              
! --------- Allocate memory ----------------------------------------------------              
              
            ALLOCATE(pBorder2%icgSend(pBorder2%nCellsSend),STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN 
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%icgSend')
            END IF ! global%error
            
! --------- Loop over cells to be sent and get global index of cell ------------            
            
            DO icl = 1,pBorder2%nCellsSend
              icg  = pBorder%icgRecv(icl)
              ict  = pGrid%cellGlob2Loc(1,icg)
              icgs = pGrid%pc2sc(icg)
                            
              IF ( icgs > pGridSerial%nCells ) THEN 
                CALL RFLU_SYPE_GetActualSerialCell(pRegionSerial,icgs,icgs2)
                                              
                icgs = icgs2
              END IF ! icgs
                                          
! ----------- Consistency check            
                
              iReg3 = pGridSerial%sc2r(icgs)  
                                 
              IF ( iReg3 /= iReg2 ) THEN 
                WRITE(errorString,'(2(1X,I7.7))') iReg3,iReg2
                CALL ErrorStop(global,ERR_REGION_IDS_INVALID,__LINE__, & 
                               TRIM(errorString))        
              END IF ! iReg3
            
! ----------- Find index of serial cell             
            
              CALL BinarySearchInteger(pGrid2%sc2pc(1:1,1:pGrid2%nCellsTot), &
                                       pGrid2%nCellsTot,icgs,iLoc) 
                                   
! ----------- Enter partitioned cell index into list of cells to be sent             
            
              IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
                icg2 = pGrid2%sc2pc(2,iLoc)
                ict2 = pGrid2%cellGlob2Loc(1,icg2)
                
                IF ( icg2 > pGrid2%nCells ) THEN ! Not actual cell 
                  CALL ErrorStop(global,ERR_CELL_KIND_INVALID,__LINE__)               
                END IF ! icg2
                
                IF ( ict /= ict2 ) THEN 
                  WRITE(errorString,'(2(1X,I2))') ict,ict2
                  CALL ErrorStop(global,ERR_CELL_TYPE_INVALID,__LINE__, & 
                                 TRIM(errorString)) 
                END IF ! ict
 
                pBorder2%icgSend(icl) = icg2
              ELSE         
                CALL ErrorStop(global,ERR_CELL_NOT_FOUND,__LINE__)            
              END IF ! iLoc
            END DO ! icg

            EXIT innerLoop  
          END IF ! pGrid2%borders
        END DO innerLoop     
      END DO ! iBorder                 
    END DO ! iReg

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining cells to be sent done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building cell communication lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_COMM_BuildCommListsCells









! ******************************************************************************
!
! Purpose: Build vertex communication lists. 
!
! Description: None.
!
! Input:
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COMM_BuildCommListsVert(regions)

    USE ModSortSearch

    USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_GetRelatedVertex, & 
                                        RFLU_SYPE_HaveSyPePatches
    USE RFLU_ModTopologyUtils, ONLY: RFLU_BuildCellVertList

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
       
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString
    LOGICAL :: foundFlag,sypeFlag
    INTEGER :: errorFlag,iBorder,iBorderSerial,iBorderSerial2,iBorder2,iLoc, &
               iReg,iReg2,ivg,ivgs,ivgs2,ivg2,ivl,ivl2,ivls2,nVertRecv, & 
               nVertRecvTemp,nVertRecvTempMax,nVertSend,nVertSendTemp, &
               nVertSendTempMax,nVertShared,nVertSharedTempMax,vListDim, &
               vListDimMax
    INTEGER, DIMENSION(:), ALLOCATABLE :: ivgRecvTemp,ivgSendTemp, &
                                          ivgSharedSorted,ivgSharedTemp, &
                                          sortKey,vList
    TYPE(t_border), POINTER :: pBorder,pBorderSerial,pBorderSerial2,pBorder2
    TYPE(t_grid), POINTER :: pGrid,pGrid2,pGridSerial 
    TYPE(t_global), POINTER :: global     
    TYPE(t_region), POINTER :: pRegion,pRegion2,pRegionSerial

! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'RFLU_COMM_BuildCommListsVert',&
  'RFLU_ModCommLists.F90')

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building vertex communication lists...'
    END IF ! global%verbLevel
    
    pRegionSerial => regions(0)
    pGridSerial   => pRegionSerial%grid
    
! ******************************************************************************
!   Check whether have symmetry or periodic boundaries in serial region
! ******************************************************************************

    syPeFlag = RFLU_SYPE_HaveSyPePatches(pRegionSerial)   
        
! ******************************************************************************
!   Determine vertices to be received for each region
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining vertices to be '// &
                               'received and shared...'
    END IF ! global%verbLevel

    DO iReg = 1,global%nRegions
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
      
      IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel      
      
! ==============================================================================
!     Loop over borders
! ==============================================================================

      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)

        IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
          WRITE(STDOUT,'(A,7X,A,1X,I3)') SOLVER_NAME,'Border:',iBorder
        END IF ! global%verbLevel   

! ------------------------------------------------------------------------------
!       Allocate temporary memory for vertex lists
! ------------------------------------------------------------------------------

! TO DO 
!       Need to find better way of estimating maximum number of vertices
!       from given list of cells. Usual estimates likely to be unreliable 
!       because they assume that boundary effects are negligible, which is 
!       almost guaranteed NOT to be the case for virtual cells.
! END TO DO
  
        nVertSendTempMax   = 8*pBorder%nCellsSend
        nVertRecvTempMax   = 8*pBorder%nCellsRecv
        nVertSharedTempMax = MAX(nVertSendTempMax,nVertRecvTempMax)

        ALLOCATE(ivgSendTemp(nVertSendTempMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ivgSendTemp')
        END IF ! global%error

        ALLOCATE(ivgRecvTemp(nVertRecvTempMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ivgRecvTemp')
        END IF ! global%error

        ALLOCATE(ivgSharedTemp(nVertSharedTempMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ivgSharedTemp')
        END IF ! global%error

! ------------------------------------------------------------------------------
!       Build list of vertices from list of cells to be sent and received. NOTE
!       these lists are NOT the vertices which actually need to be sent and 
!       received
! ------------------------------------------------------------------------------

        IF ( pBorder%nCellsSend > 0 ) THEN 
          CALL RFLU_BuildCellVertList(global,pGrid,pBorder%icgSend, & 
                                      pBorder%nCellsSend,ivgSendTemp, & 
                                      nVertSendTempMax,nVertSendTemp)                                      
          CALL QuickSortInteger(ivgSendTemp(1:nVertSendTemp),nVertSendTemp)
        ELSE 
          nVertSendTemp = 0
        END IF ! pBorder%nCellsSend
      
        IF ( pBorder%nCellsRecv > 0 ) THEN 
          CALL RFLU_BuildCellVertList(global,pGrid,pBorder%icgRecv, & 
                                      pBorder%nCellsRecv,ivgRecvTemp, & 
                                      nVertRecvTempMax,nVertRecvTemp)
          CALL QuickSortInteger(ivgRecvTemp(1:nVertRecvTemp),nVertRecvTemp)
        ELSE 
          nVertRecvTemp = 0 
        END IF ! pBorder%nCellsRecv
              
! ------------------------------------------------------------------------------
!       Build list of shared vertices from vertices which are shared among 
!       lists of vertices to be sent and received and list of vertices to be
!       received as those vertices which are not shared. List of vertices to be
!       sent is built below.
! ------------------------------------------------------------------------------

        IF ( nVertSendTemp > 0 .AND. nVertRecvTemp > 0 ) THEN 
          CALL RemoveCommonSortedIntegersFancy(ivgSendTemp(1:nVertSendTemp), & 
                                               nVertSendTemp,nVertSend, & 
                                               ivgRecvTemp(1:nVertRecvTemp), &
                                               nVertRecvTemp,nVertRecv, & 
                                               ivgSharedTemp, &
                                               nVertSharedTempMax, &
                                               nVertShared,errorFlag)
! TEMPORARY
          IF ( errorFlag /= ERR_NONE ) THEN
            WRITE(*,*) 'RemoveCommonSortedIntegersFancy returned error!'
            STOP
          END IF ! errorFlag
! END TEMPORARY
        ELSE 
          nVertShared = 0 

          IF ( nVertSendTemp > 0 ) THEN 
            nVertSend = nVertSendTemp
          ELSE 
            nVertSend = 0
          END IF ! nVertSendTemp

          IF ( nVertRecvTemp > 0 ) THEN 
            nVertRecv = nVertRecvTemp
          ELSE 
            nVertRecv = 0
          END IF ! nVertRecvTemp
        END IF ! nVertSendTemp                                              
      
! ------------------------------------------------------------------------------
!       Deallocate temporary memory and copy vertex lists
! ------------------------------------------------------------------------------

        DEALLOCATE(ivgSendTemp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ivgSendTemp')
        END IF ! global%error

        pBorder%nVertRecv   = nVertRecv
        pBorder%nVertShared = nVertShared 

        IF ( pBorder%nVertRecv > 0 ) THEN
          ALLOCATE(pBorder%ivgRecv(pBorder%nVertRecv),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgRecv')
          END IF ! global%error
        ELSE 
          NULLIFY(pBorder%ivgRecv)
        END IF ! pBorder%nVertRecv

        IF ( pBorder%nVertShared > 0 ) THEN 
          ALLOCATE(pBorder%ivgShared(pBorder%nVertShared),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgShared')
          END IF ! global%error        
        ELSE 
          NULLIFY(pBorder%ivgShared)
        END IF ! pBorder%nVertShared

        DO ivl = 1,pBorder%nVertRecv
          pBorder%ivgRecv(ivl) = ivgRecvTemp(ivl)
        END DO ! ivl

        DO ivl = 1,pBorder%nVertShared
          pBorder%ivgShared(ivl) = ivgSharedTemp(ivl)
        END DO ! ivl

        DEALLOCATE(ivgRecvTemp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ivgRecvTemp')
        END IF ! global%error

        DEALLOCATE(ivgSharedTemp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ivgSharedTemp')
        END IF ! global%error
      END DO ! iBorder
    END DO ! iReg

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Determining vertices to be '// &
                               'received and shared done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Make shared vertex lists consistent
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Ensuring consistency of shared-vertex lists...'
    END IF ! global%verbLevel

    DO iReg = 1,global%nRegions
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      IF ( global%verbLevel >= VERBOSE_HIGH) THEN 
        WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel
      
      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)

        iReg2    = pBorder%iRegionGlobal
        iBorder2 = pBorder%iBorder 
        
        pRegion2 => regions(iReg2)
        pGrid2   => pRegion2%grid
        pBorder2 => pGrid2%borders(iBorder2)

        IF ( pBorder%nVertShared /= pBorder2%nVertShared ) THEN              
          CALL ErrorStop(global,ERR_NVERTSHARED_MISMATCH,__LINE__) 
        END IF ! pBorder%nVertShared

        IF ( (iReg2 < iReg) .OR. & 
             ((iReg2 == iReg) .AND. (iBorder2 < iBorder)) ) THEN 
          DO ivl = 1,pBorder%nVertShared
            ivg  = pBorder%ivgShared(ivl)
            ivgs = pGrid%pv2sv(ivg)

            CALL BinarySearchInteger(pGrid2%sv2pv(1:1,1:pGrid2%nVertTot), &
                                     pGrid2%nVertTot,ivgs,iLoc) 

            IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
              ivg2 = pGrid2%sv2pv(2,iLoc)

              IF ( ivg2 > pGrid2%nVert ) THEN ! Not actual vertex
                CALL ErrorStop(global,ERR_VERTEX_KIND_INVALID,__LINE__)               
              END IF ! ivg2

              pBorder2%ivgShared(ivl) = ivg2
            ELSE         
              IF ( sypeFlag .EQV. .TRUE. ) THEN
                CALL RFLU_SYPE_GetRelatedVertex(pRegionSerial,pRegion2,ivgs,ivg2)

                pBorder2%ivgShared(ivl) = ivg2             
              ELSE 
                CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)            
              END IF ! sypeFlag 
            END IF ! iLoc            
          END DO ! ivl  
        END IF ! iReg2
      END DO ! iBorder
    END DO ! iReg

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Ensuring consistency of '// &
                               'shared-vertex lists done.'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Determine vertices to be sent for each region and border of other region
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Determining vertices to be sent...'
    END IF ! global%verbLevel

    DO iReg = 1,global%nRegions
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
          
      IF ( global%verbLevel >= VERBOSE_HIGH) THEN 
        WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel          
          
! ==============================================================================
!     Loop over borders
! ==============================================================================          
          
      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)

! ------------------------------------------------------------------------------
!       Find index of border of neighboring region whose neighboring region is
!       equal to region currently being looped over
! ------------------------------------------------------------------------------        

        iReg2    = pBorder%iRegionGlobal
        iBorder2 = pBorder%iBorder 
        
        pRegion2 => regions(iReg2)
        pGrid2   => pRegion2%grid
        pBorder2 => pGrid2%borders(iBorder2)
                          
! ----- Set number of vertices to be sent --------------------------------------
           
        pBorder2%nVertSend = pBorder%nVertRecv
              
! ----- Allocate memory --------------------------------------------------------              
              
        IF ( pBorder2%nVertSend > 0 ) THEN        
          ALLOCATE(pBorder2%ivgSend(pBorder2%nVertSend),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgSend')
          END IF ! global%error           
        ELSE 
          NULLIFY(pBorder2%ivgSend)
        END IF ! pBorder2%nVertSend

! ----- Loop over vertices to be sent and get global index of vertex -----------            
            
        DO ivl = 1,pBorder2%nVertSend
          ivg  = pBorder%ivgRecv(ivl)
          ivgs = pGrid%pv2sv(ivg)
                          
! ------- Find index of serial vertex            

          CALL BinarySearchInteger(pGrid2%sv2pv(1:1,1:pGrid2%nVertTot), &
                                   pGrid2%nVertTot,ivgs,iLoc) 

! ------- Enter partitioned vertex index into list of vertices to be sent             
            
          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
            ivg2 = pGrid2%sv2pv(2,iLoc)

            IF ( ivg2 > pGrid2%nVert ) THEN ! Not actual vertex
              CALL ErrorStop(global,ERR_VERTEX_KIND_INVALID,__LINE__)               
            END IF ! ivg2

            pBorder2%ivgSend(ivl) = ivg2
          ELSE         
            IF ( sypeFlag .EQV. .TRUE. ) THEN  
              CALL RFLU_SYPE_GetRelatedVertex(pRegionSerial,pRegion2,ivgs,ivg2)

              pBorder2%ivgSend(ivl) = ivg2
            ELSE
              CALL ErrorStop(global,ERR_VERTEX_NOT_FOUND,__LINE__)            
            END IF ! sypeFlag
          END IF ! iLoc
        END DO ! ivl
      END DO ! iBorder
    END DO ! iReg

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Determining vertices to be sent done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building vertex communication lists done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_COMM_BuildCommListsVert









! ******************************************************************************
!
! Purpose: Check counting of number of borders. 
!
! Description: None.
!
! Input:
!   regions             Regions data
!
! Output: None.
!
! Notes: 
!   1. Needed because when have too many partitions, in particular with higher-
!      order scheme, can have virtual layers extend through a neighboring region
!      into another region. This becomes a problem because in these cases, the 
!      communication pattern may become asymmetric, i.e., region A may need to 
!      receive data from region B, but region B may NOT need to receive data 
!      from region A. Because the routine which counts the number of borders 
!      of a given region determines the number of borders by the number of 
!      regions from which data needs to be received, the count can be 
!      incomplete. For example, in one ACM case, region 25 needed to receive
!      2 cells from region 20, but region 20 did not need to receive any data
!      from region 25. Then the communication pattern is wrong, because region
!      25 would not send anything to 20. This routine checks whether the 
!      communication pattern in symmetric and makes it symmetric if necessary.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_COMM_CheckCountBorders(regions)

    USE ModSortSearch

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
    
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,iBorder,iBorder2,iLoc,iReg,iReg2
    INTEGER, DIMENSION(:), ALLOCATABLE :: borderCntr
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid,pGrid2
    TYPE(t_region), POINTER :: pRegion,pRegion2

! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(1)%global

    CALL RegisterFunction(global,'RFLU_COMM_CheckCountBorders',&
  'RFLU_ModCommLists.F90')

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking counting of borders...'      
    END IF ! global%verbLevel

! ******************************************************************************
!   
! ******************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
    
      IF ( global%verbLevel >= VERBOSE_HIGH) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
      END IF ! global%verbLevel    
    
      DO iBorder = 1,pGrid%nBorders
        iReg2 = pGrid%borderCntr(iBorder)
        
        pRegion2 => regions(iReg2)
        pGrid2   => pRegion2%grid
        
        CALL BinarySearchInteger(pGrid2%borderCntr(1:pGrid2%nBorders), & 
                                 pGrid2%nBorders,iReg,iLoc)
                                
        IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
          global%warnCounter = global%warnCounter + 1
        
          IF ( global%verbLevel >= VERBOSE_NONE ) THEN       
            WRITE(STDOUT,'(A,3X,A,2(1X,A,1X,I5),A)') SOLVER_NAME, &
                  '*** WARNING ***','Border asymmetry between regions', &
                  iReg,'and',iReg2,'.'
            WRITE(STDOUT,'(A,19X,A,1X,I5,A)') SOLVER_NAME, & 
                  'Correcting border data structure for region',iReg2,'.'
          END IF ! global%verbLevel
       
          IF ( pGrid2%nBorders == SIZE(pGrid2%borderCntr) ) THEN 
            ALLOCATE(borderCntr(pGrid2%nBorders),STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN 
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderCntr')
            END IF ! global%error

            DO iBorder2 = 1,pGrid2%nBorders
              borderCntr(iBorder2) = pGrid2%borderCntr(iBorder2)
            END DO ! iBorder2

            CALL RFLU_COMM_DestroyBorderCntr(pRegion2)
            CALL RFLU_COMM_CreateBorderCntr(pRegion2,2*pGrid2%nBorders)

            DO iBorder2 = 1,pGrid2%nBorders
              pGrid2%borderCntr(iBorder2) = borderCntr(iBorder2)
            END DO ! iBorder2 

            DO iBorder2 = pGrid2%nBorders+1,SIZE(pGrid2%borderCntr)
              pGrid2%borderCntr(iBorder2) = 0
            END DO ! iBorder2          

            DEALLOCATE(borderCntr,STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN 
              CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderCntr')
            END IF ! global%error                    
          END IF ! pGrid2%nBorders
          
          pGrid2%nBorders = pGrid2%nBorders + 1
        
          pGrid2%borderCntr(pGrid2%nBorders) = iReg
        
          CALL QuickSortInteger(pGrid2%borderCntr(1:pGrid2%nBorders), & 
                                pGrid2%nBorders)  
        END IF ! iLoc
      END DO ! iBorder
    END DO ! iReg
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking counting of borders done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_CheckCountBorders










! ******************************************************************************
!
! Purpose: Count number of borders
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
    
  SUBROUTINE RFLU_COMM_CountBorders(pRegion,pRegionSerial)

    USE ModSortSearch
    
    USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_GetActualSerialCell, & 
                                        RFLU_SYPE_HaveSyPePatches

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

    LOGICAL :: syPeFlag
    INTEGER :: errorFlag,iBorder,icg,icgs,icgs2,iLoc,iReg,iReg2,nBordersMax
    INTEGER, DIMENSION(:), ALLOCATABLE :: borderCntr
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid,pGridSerial

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_CountBorders',&
  'RFLU_ModCommLists.F90')

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Counting borders...'      
    END IF ! global%verbLevel

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set pointers and initialize 
! ******************************************************************************

    pGrid       => pRegion%grid  
    pGridSerial => pRegionSerial%grid      

! ******************************************************************************
!   Check whether have symmetry or periodic boundaries in serial region
! ******************************************************************************

    syPeFlag = RFLU_SYPE_HaveSyPePatches(pRegionSerial)
    
! ******************************************************************************
!   Loop over virtual cells and increment border counter
! ******************************************************************************

    pGrid%nBorders = 0

    DO icg = pGrid%nCells+1,pGrid%nCellsTot
      icgs = pGrid%pc2sc(icg)    
      iReg = pGridSerial%sc2r(icgs)
            
      IF ( iReg /= pRegion%iRegionGlobal ) THEN 
        IF ( pGrid%nBorders > 0 ) THEN 
          CALL BinarySearchInteger(pGrid%borderCntr(1:pGrid%nBorders), & 
                                   pGrid%nBorders,iReg,iLoc)                                   
        ELSE 
          iLoc = ELEMENT_NOT_FOUND
        END IF ! pGrid%nBorders      
      ELSE 
        IF ( syPeFlag .EQV. .TRUE. ) THEN 
          IF ( icgs > pGridSerial%nCells ) THEN 
            CALL RFLU_SYPE_GetActualSerialCell(pRegionSerial,icgs,icgs2)
                 
            iReg = pGridSerial%sc2r(icgs2)  
                                  
            IF ( pGrid%nBorders > 0 ) THEN 
              CALL BinarySearchInteger(pGrid%borderCntr(1:pGrid%nBorders), & 
                                       pGrid%nBorders,iReg,iLoc)
            ELSE 
              iLoc = ELEMENT_NOT_FOUND
            END IF ! pGrid%nBorders                       
          ELSE 
! TEMPORARY
            WRITE(*,*) 'ERROR!'
            STOP
! END TEMPORARY            
          END IF ! icgs
        ELSE 
! TEMPORARY
          WRITE(*,*) 'ERROR!'
          STOP
! END TEMPORARY                     
        END IF ! syPeFlag
      END IF ! pGrid%nBorders
      
      IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
        IF ( pGrid%nBorders == SIZE(pGrid%borderCntr) ) THEN
          ALLOCATE(borderCntr(pGrid%nBorders),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderCntr')
          END IF ! global%error
          
          DO iBorder = 1,pGrid%nBorders
            borderCntr(iBorder) = pGrid%borderCntr(iBorder)
          END DO ! iBorder
           
          CALL RFLU_COMM_DestroyBorderCntr(pRegion)
          CALL RFLU_COMM_CreateBorderCntr(pRegion,2*pGrid%nBorders)
          
          DO iBorder = 1,pGrid%nBorders
            pGrid%borderCntr(iBorder) = borderCntr(iBorder)
          END DO ! iBorder 
          
          DO iBorder = pGrid%nBorders+1,SIZE(pGrid%borderCntr)
            pGrid%borderCntr(iBorder) = 0
          END DO ! iBorder          

          DEALLOCATE(borderCntr,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderCntr')
          END IF ! global%error                                      
        END IF ! pGrid%nBorders 
           
        pGrid%nBorders = pGrid%nBorders + 1
        
        pGrid%borderCntr(pGrid%nBorders) = iReg
        
        CALL QuickSortInteger(pGrid%borderCntr(1:pGrid%nBorders), & 
                              pGrid%nBorders)  
      END IF ! iLoc
    END DO ! icg    

! ******************************************************************************
!   Check that do not have any virtual cells which are mapped to own region if 
!   do not have any symmetry or periodic patches
! ******************************************************************************

    IF ( syPeFlag .EQV. .FALSE. ) THEN     
      CALL BinarySearchInteger(pGrid%borderCntr(1:pGrid%nBorders), & 
                               pGrid%nBorders,pRegion%iRegionGlobal,iLoc)
      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN    
        CALL ErrorStop(global,ERR_PARTITION_INVALID,__LINE__)  
      END IF ! iLoc
    END IF ! syPeFlag

! ******************************************************************************
!   Write info
! ******************************************************************************
    
    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME,'Number of borders:', & 
                                     pGrid%nBorders
    END IF ! global%verbLevel
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Counting borders done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_CountBorders








! ******************************************************************************
!
! Purpose: Count number of borders for serial region.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_COMM_CountBordersSerial(pRegion)

    USE ModBndPatch, ONLY: t_patch

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

    CALL RegisterFunction(global,'RFLU_COMM_CountBordersSerial',&
  'RFLU_ModCommLists.F90')

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Counting serial borders...'      
    END IF ! global%verbLevel

    IF ( pRegion%iRegionGlobal /= 0 ) THEN 
      CALL ErrorStop(global,ERR_REGION_ID_INVALID,__LINE__)
    END IF ! pRegion%iRegionGlobal
    
! ******************************************************************************
!   Set grid pointers
! ******************************************************************************

    pGrid => pRegion%grid  

! ******************************************************************************
!   Count number of borders
! ******************************************************************************

    pGrid%nBorders = 0

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( (pPatch%bcType == BC_SYMMETRY) .OR. & 
           (pPatch%bcType == BC_PERIODIC) ) THEN 
        pGrid%nBorders = pGrid%nBorders + 1
        
        pPatch%iBorder = pGrid%nBorders
      ELSE 
        pPatch%iBorder = PATCH_IBORDER_DEFAULT
      END IF ! pPatch%bcType
    END DO ! iPatch
        
    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME,'Number of borders:', & 
                                     pGrid%nBorders
    END IF ! global%verbLevel        
        
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Counting serial borders done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_CountBordersSerial






! ******************************************************************************
!
! Purpose: Create border counter.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   nBordersMaxOpt      Desired Size of borderCntr array
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_COMM_CreateBorderCntr(pRegion,nBordersMaxOpt)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN), OPTIONAL :: nBordersMaxOpt
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,iBorder,nBorders,nBordersMax
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_CreateBorderCntr',&
  'RFLU_ModCommLists.F90')
    
! ******************************************************************************
!   Set grid pointer and variables
! ******************************************************************************

    pGrid => pRegion%grid       

    nBordersMax = 10 ! Nominal maximum number of borders

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    IF ( PRESENT(nBordersMaxOpt) ) THEN 
      nBorders = nBordersMaxOpt
    ELSE
      nBorders = nBordersMax
    END IF ! PRESENT

    ALLOCATE(pGrid%borderCntr(nBorders),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%borderCntr')
    END IF ! global%error

    DO iBorder = 1,nBorders
      pGrid%borderCntr(iBorder) = 0
    END DO ! iBorder
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_CreateBorderCntr








! ******************************************************************************
!
! Purpose: Create borders
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   createMode          Creation mode 
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_COMM_CreateBorders(pRegion,createMode)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN), OPTIONAL :: createMode
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_CreateBorders',&
  'RFLU_ModCommLists.F90')

    IF ( global%myProcid == MASTERPROC .AND. &  
         global%verbLevel >= VERBOSE_HIGH ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating borders...'      
    END IF ! global%myProcid 

    IF ( global%myProcid == MASTERPROC .AND. &  
         global%verbLevel >= VERBOSE_HIGH ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid 
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid       

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    IF ( pGrid%nBorders > 0 ) THEN 
      ALLOCATE(pGrid%borders(pGrid%nBorders),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%borders')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%borders)
    END IF ! pGrid%nBorders

! ******************************************************************************
!   Initialize. NOTE in preprocessor do not know dimensions when borders are
!   created (only know the number), but in postprocessor and solver DO know 
!   dimensions. NOTE only initialize iProc and iRegionLocal, will be set later
!   in RFLU_COMM_GetProcLocRegIds.
! ******************************************************************************
    
    IF ( PRESENT(createMode) ) THEN 
      IF ( createMode == CREATE_BORDERS_MODE_DIM_KNOWN ) THEN 
        DO iBorder = 1,pGrid%nBorders
          pBorder => pGrid%borders(iBorder)
      
          pBorder%iRegionGlobal = pGrid%borderInfo(BORDER_INFO_IRGLOB,iBorder)
          pBorder%iBorder       = pGrid%borderInfo(BORDER_INFO_IBORD ,iBorder)
          pBorder%nCellsSend    = pGrid%borderInfo(BORDER_INFO_NCSEND,iBorder)      
          pBorder%nCellsRecv    = pGrid%borderInfo(BORDER_INFO_NCRECV,iBorder)
          pBorder%nVertSend     = pGrid%borderInfo(BORDER_INFO_NVSEND,iBorder)      
          pBorder%nVertRecv     = pGrid%borderInfo(BORDER_INFO_NVRECV,iBorder)
          pBorder%nVertShared   = pGrid%borderInfo(BORDER_INFO_NVSHAR,iBorder)          
          
          pBorder%iProc         = CRAZY_VALUE_INT 
          pBorder%iRegionLocal  = CRAZY_VALUE_INT
        END DO ! iBorder
      ELSE
        DO iBorder = 1,pGrid%nBorders
          pBorder => pGrid%borders(iBorder)
      
          pBorder%iRegionGlobal = CRAZY_VALUE_INT
          pBorder%iBorder       = CRAZY_VALUE_INT

          pBorder%nCellsSend    = 0    
          pBorder%nCellsRecv    = 0
          pBorder%nVertSend     = 0   
          pBorder%nVertRecv     = 0
          pBorder%nVertShared   = 0         
          
          pBorder%iProc         = CRAZY_VALUE_INT 
          pBorder%iRegionLocal  = CRAZY_VALUE_INT
        END DO ! iBorder                  
      END IF ! createMode
    END IF ! PRESENT
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &  
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating borders done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_CreateBorders







! ******************************************************************************
!
! Purpose: Create communication lists.
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
    
  SUBROUTINE RFLU_COMM_CreateCommLists(pRegion)

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

    INTEGER :: errorFlag,iBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_CreateCommLists',&
  'RFLU_ModCommLists.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating communication lists...'      
    END IF ! global%myProcid
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid       

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      ALLOCATE(pBorder%icgSend(pBorder%nCellsSend),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%icgSend')
      END IF ! global%error

      ALLOCATE(pBorder%icgRecv(pBorder%nCellsRecv),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%icgRecv')
      END IF ! global%error       

      ALLOCATE(pBorder%ivgSend(pBorder%nVertSend),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgSend')
      END IF ! global%error

      ALLOCATE(pBorder%ivgRecv(pBorder%nVertRecv),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgRecv')
      END IF ! global%error       

      ALLOCATE(pBorder%ivgShared(pBorder%nVertShared),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%ivgShared')
      END IF ! global%error
    END DO ! iBorder
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating communication lists done.'      
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_CreateCommLists








! ******************************************************************************
!
! Purpose: Destroy border counter.
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
    
  SUBROUTINE RFLU_COMM_DestroyBorderCntr(pRegion)

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

    CALL RegisterFunction(global,'RFLU_COMM_DestroyBorderCntr',&
  'RFLU_ModCommLists.F90')
    
! ******************************************************************************
!   Set grid pointer 
! ******************************************************************************

    pGrid => pRegion%grid       

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(pGrid%borderCntr,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%borderCntr')
    END IF ! global%error
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_DestroyBorderCntr









! ******************************************************************************
!
! Purpose: Destroy borders
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
    
  SUBROUTINE RFLU_COMM_DestroyBorders(pRegion)

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

    INTEGER :: errorFlag,iBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_DestroyBorders',&
  'RFLU_ModCommLists.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN     
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying borders...'      
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid       

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    DEALLOCATE(pGrid%borders,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%borders')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying borders done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_DestroyBorders






! ******************************************************************************
!
! Purpose: Destroy communication lists.
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
    
  SUBROUTINE RFLU_COMM_DestroyCommLists(pRegion)

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

    INTEGER :: errorFlag,iBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_DestroyCommLists',&
  'RFLU_ModCommLists.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying communication lists...'      
    END IF !  global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid       

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      DEALLOCATE(pBorder%icgSend,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%icgSend')
      END IF ! global%error

      DEALLOCATE(pBorder%icgRecv,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%icgRecv')
      END IF ! global%error 

      DEALLOCATE(pBorder%ivgSend,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%ivgSend')
      END IF ! global%error

      DEALLOCATE(pBorder%ivgRecv,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%ivgRecv')
      END IF ! global%error       

      DEALLOCATE(pBorder%ivgShared,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%ivgShared')
      END IF ! global%error
    END DO ! iBorder
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying communication lists done.'      
    END IF !  global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_DestroyCommLists







! ******************************************************************************
!
! Purpose: Get process and local region ids for borders.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine must be called after borders were created and region
!      mapping was read.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_COMM_GetProcLocRegIds(pRegion)

    USE RFLU_ModRegionMapping, ONLY: RFLU_CloseRegionMappingFile, & 
                                     RFLU_GetProcLocRegIds, & 
                                     RFLU_OpenRegionMappingFile

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

    INTEGER :: errorFlag,iBorder,nRegIds
    INTEGER, DIMENSION(:), ALLOCATABLE :: locRegIds,procIds,regIds
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_GetProcLocRegIds',&
  'RFLU_ModCommLists.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting border process ids..'      
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid
    
! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid       

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(regIds(pGrid%nBorders),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'regIds')
    END IF ! global%error

    ALLOCATE(procIds(pGrid%nBorders),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'procIds')
    END IF ! global%error

    ALLOCATE(locRegIds(pGrid%nBorders),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'locRegIds')
    END IF ! global%error
    
! ******************************************************************************
!   Set region ids, get process ids from mapping file, and impose process ids
! ******************************************************************************
    
    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
      
      regIds(iBorder) = pBorder%iRegionGlobal       
    END DO ! iBorder
    
    CALL RFLU_OpenRegionMappingFile(global)
    CALL RFLU_GetProcLocRegIds(global,regIds,pGrid%nBorders,procIds,locRegIds)
    CALL RFLU_CloseRegionMappingFile(global)
    
    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
      
      pBorder%iProc        = procIds(iBorder) - 1 ! NOTE offset      
      pBorder%iRegionLocal = locRegIds(iBorder)
    END DO ! iBorder    

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(regIds,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'regIds')
    END IF ! global%error

    DEALLOCATE(procIds,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'procIds')
    END IF ! global%error
        
    DEALLOCATE(locRegIds,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'locRegIds')
    END IF ! global%error
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting border process ids done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COMM_GetProcLocRegIds






! ******************************************************************************
!
! Purpose: Read communication lists
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

  SUBROUTINE RFLU_COMM_ReadCommLists(pRegion)

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

    INTEGER :: errorFlag,iBorder,icl,iFile,ivl,loopCounter,nBorders, & 
               nCellsRecv,nCellsSend,nVertRecv,nVertSend,nVertShared
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_ReadCommLists',&
  'RFLU_ModCommLists.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading communication lists...'
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Open file
! ******************************************************************************

    iFile = IF_COMM_LISTS

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.com', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# ROCFLU communication lists file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM
    
! ******************************************************************************
!   Dimensions
! ******************************************************************************
  
    pGrid => pRegion%grid      

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM
    
    READ(iFile,'(I8)') nBorders

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================
    
    IF ( nBorders /= pGrid%nBorders ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nBorders    

! ******************************************************************************
!   Rest of file
! ******************************************************************************

    loopCounter = 0

! ==============================================================================  
!   Set up infinite loop
! ==============================================================================  

    DO ! set up infinite loop
      loopCounter = loopCounter + 1
    
      READ(iFile,'(A)') sectionString
    
      SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!       Information
! ------------------------------------------------------------------------------

        CASE ( '# Information' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Information...'
          END IF ! global%myProcid 
              
          DO iBorder = 1,pGrid%nBorders          
            pBorder => pGrid%borders(iBorder)
    
            READ(iFile,'(2(I8))') pBorder%iRegionGlobal,pBorder%iBorder  
          END DO ! iBorder 
    
! ------------------------------------------------------------------------------
!       Cells
! ------------------------------------------------------------------------------

        CASE ( '# Cells' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
          END IF ! global%myProcid 
              
          DO iBorder = 1,pGrid%nBorders          
            pBorder => pGrid%borders(iBorder)
    
            READ(iFile,'(2(I8))') nCellsSend,nCellsRecv

            IF ( nCellsSend /= pBorder%nCellsSend .OR. & 
                 nCellsRecv /= pBorder%nCellsRecv ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nCellsSend              
 
            READ(iFile,'(10(I8))') (pBorder%icgSend(icl), & 
                                    icl=1,pBorder%nCellsSend)
            READ(iFile,'(10(I8))') (pBorder%icgRecv(icl), & 
                                    icl=1,pBorder%nCellsRecv)  
          END DO ! iBorder   

! ------------------------------------------------------------------------------
!       Vertices
! ------------------------------------------------------------------------------

        CASE ( '# Vertices' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Vertices...'
          END IF ! global%myProcid 
              
          DO iBorder = 1,pGrid%nBorders          
            pBorder => pGrid%borders(iBorder)
    
            READ(iFile,'(3(I8))') nVertSend,nVertRecv,nVertShared
                        
            IF ( nVertSend /= pBorder%nVertSend .OR. & 
                 nVertRecv /= pBorder%nVertRecv .OR. &
                 nVertShared /= pBorder%nVertShared ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nVertSend              
 
            READ(iFile,'(10(I8))') (pBorder%ivgSend(ivl), & 
                                    ivl=1,pBorder%nVertSend)
            READ(iFile,'(10(I8))') (pBorder%ivgRecv(ivl), & 
                                    ivl=1,pBorder%nVertRecv) 
            READ(iFile,'(10(I8))') (pBorder%ivgShared(ivl), & 
                                    ivl=1,pBorder%nVertShared)   
          END DO ! iBorder   
   
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_HIGH) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%myProcid      

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel >= VERBOSE_HIGH) THEN  
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

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading communication lists done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_COMM_ReadCommLists







! ******************************************************************************
!
! Purpose: Write communication lists
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

  SUBROUTINE RFLU_COMM_WriteCommLists(pRegion)

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

    INTEGER :: errorFlag,iBorder,icl,iFile,ivl
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COMM_WriteCommLists',&
  'RFLU_ModCommLists.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing communication lists...'
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ******************************************************************************
!   Open file
! ******************************************************************************

    iFile = IF_COMM_LISTS

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.com', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid 

    sectionString = '# ROCFLU communication lists file'  
    WRITE(iFile,'(A)') TRIM(sectionString)  

! ******************************************************************************
!   Dimensions
! ******************************************************************************
  
    pGrid => pRegion%grid  

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Dimensions...'
    END IF ! global%myProcid    
    
    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pGrid%nBorders

! ******************************************************************************
!   Information 
! ******************************************************************************
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Information...'
    END IF ! global%myProcid    
    
    sectionString = '# Information'
    WRITE(iFile,'(A)') TRIM(sectionString)     

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
      WRITE(iFile,'(2(I8))') pBorder%iRegionGlobal,pBorder%iBorder
    END DO ! iBorder

! ******************************************************************************
!   Cells 
! ******************************************************************************
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
    END IF ! global%myProcid    
    
    sectionString = '# Cells'
    WRITE(iFile,'(A)') TRIM(sectionString)     

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
      WRITE(iFile,'( 2(I8))') pBorder%nCellsSend,pBorder%nCellsRecv
      WRITE(iFile,'(10(I8))') (pBorder%icgSend(icl),icl=1,pBorder%nCellsSend)
      WRITE(iFile,'(10(I8))') (pBorder%icgRecv(icl),icl=1,pBorder%nCellsRecv)      
    END DO ! iBorder

! ******************************************************************************
!   Vertices
! ******************************************************************************
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Vertices...'
    END IF ! global%myProcid    
    
    sectionString = '# Vertices'
    WRITE(iFile,'(A)') TRIM(sectionString)     

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
            
      WRITE(iFile,'( 3(I8))') pBorder%nVertSend,pBorder%nVertRecv, &
                              pBorder%nVertShared
      WRITE(iFile,'(10(I8))') (pBorder%ivgSend(ivl),ivl=1,pBorder%nVertSend)
      WRITE(iFile,'(10(I8))') (pBorder%ivgRecv(ivl),ivl=1,pBorder%nVertRecv)
      WRITE(iFile,'(10(I8))') (pBorder%ivgShared(ivl),ivl=1,pBorder%nVertShared)       
    END DO ! iBorder
    
! ******************************************************************************
!   End marker
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%myProcid

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString) 

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

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing communication lists done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_COMM_WriteCommLists


  


! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModCommLists


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModCommLists.F90,v $
! Revision 1.14  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.11  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.10  2006/03/28 21:29:09  haselbac
! Backed out changes for read/writing comm lists
!
! Revision 1.9  2006/03/25 21:51:10  haselbac
! Substantial changes because of sype patches
!
! Revision 1.8  2005/09/19 18:39:59  haselbac
! Bug fix for asymmetric borders: Now use pGrid%borderCntr and check
!
! Revision 1.7  2005/06/20 17:04:45  haselbac
! Bug fix: New way of building shared v lists to deal w odd virtual cell topologies
!
! Revision 1.6  2005/04/18 20:25:51  haselbac
! Bug fix: pGrid used before set
!
! Revision 1.5  2005/04/15 15:06:48  haselbac
! Added routines to build vertex comm lists, more extensions, clean-up
!
! Revision 1.4  2005/01/17 19:53:20  haselbac
! Added proper error treatment
!
! Revision 1.3  2005/01/14 21:21:31  haselbac
! Added routines, init of iProc
!
! Revision 1.2  2004/12/29 21:05:04  haselbac
! Various enhancements
!
! Revision 1.1  2004/12/04 03:44:14  haselbac
! Initial revision
!
! ******************************************************************************





















