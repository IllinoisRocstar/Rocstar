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
!*******************************************************************************
!
! Purpose: Get special cells.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region data
! 
! Output: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: RFLU_GetSpecialCells.F90,v 1.3 2008/12/06 08:44:54 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!*******************************************************************************

SUBROUTINE RFLU_GetSpecialCells(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError  
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModSortSearch

  USE RFLU_ModGrid

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

  CHARACTER :: infoType
  CHARACTER(CHRLEN) :: RCSIdentString  
  INTEGER :: cellIndx,errorFlag,faceIndx,iCellsSpecial,icl,iloc, & 
             nVertPerCell,patchIndx,vertIndx
  INTEGER :: v(8)
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetSpecialCells.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GetSpecialCells', &
                        'RFLU_GetSpecialCells.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting special cells...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal  
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and initialize
! ******************************************************************************

  pGrid => pRegion%grid

  iCellsSpecial = 0 
  pGrid%cellsSpecial(1:NCELLS_SPECIAL_MAX) = 0

! ******************************************************************************
! Get information from user
! ******************************************************************************

  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter information on special cells:' 
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'b  - cell adjacent to boundary face'   
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'c  - single cell'
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'f  - cells adjacent to interior face'
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'v  - cells meeting at vertex'      
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'q  - quit'

! ==============================================================================
! Set up infinite loop
! ==============================================================================

  DO
  
! ------------------------------------------------------------------------------
!   Enter information type
! ------------------------------------------------------------------------------  
   
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information type:'
    READ(STDIN,'(A)') infoType
  
    SELECT CASE ( infoType )
    
! --- cell adjacent to boundary face -------------------------------------------

      CASE ( 'b' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter patch index:'
        READ(STDIN,*,IOSTAT=errorFlag) patchIndx

        IF ( errorFlag /= ERR_NONE ) THEN 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag        
          
        IF ( patchIndx > 0 .AND. patchIndx <= pGrid%nPatches ) THEN
          pPatch => pRegion%patches(patchIndx)
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter face index:'
          READ(STDIN,*,IOSTAT=errorFlag) faceIndx           

          IF ( errorFlag /= ERR_NONE ) THEN 
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                     '*** WARNING *** Invalid input.'
            CYCLE
          END IF ! errorFlag 

          IF ( faceIndx > 0 .AND. faceIndx <= pPatch%nBFacesTot ) THEN
            IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
              CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
            END IF ! iCellsSpecial          
           
            iCellsSpecial = iCellsSpecial + 1          
            pGrid%cellsSpecial(iCellsSpecial) = pPatch%bf2c(faceIndx)
  
            WRITE(STDOUT,'(A,5X,A,1X,I8)') SOLVER_NAME,'Added cell:', & 
                                           pPatch%bf2c(faceIndx)        
          ELSE 
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                     '*** WARNING *** Invalid input.' 
            CYCLE 
          END IF ! faceIndx        
        ELSE 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! patchIndx
    
! --- single cell --------------------------------------------------------------    
    
      CASE ( 'c' )
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter cell index:'
        READ(STDIN,*,IOSTAT=errorFlag) cellIndx
        
        IF ( errorFlag /= ERR_NONE ) THEN 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag        
        
        IF ( cellIndx > 0 .AND. cellIndx <= pGrid%nCellsTot ) THEN 
          IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
            CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
          END IF ! iCellsSpecial         
        
          iCellsSpecial = iCellsSpecial + 1          
          pGrid%cellsSpecial(iCellsSpecial) = cellIndx
        ELSE 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! cellIndx

! --- cells adjacent to interior face ------------------------------------------    

      CASE ( 'f' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter interior face index:'
        READ(STDIN,*,IOSTAT=errorFlag) faceIndx
 
        IF ( errorFlag /= ERR_NONE ) THEN 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag         
        
        IF ( faceIndx > 0 .AND. faceIndx <= pGrid%nFacesTot ) THEN
          IF ( iCellsSpecial == NCELLS_SPECIAL_MAX-1 ) THEN 
            CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
          END IF ! iCellsSpecial         
         
          iCellsSpecial = iCellsSpecial + 1          
          pGrid%cellsSpecial(iCellsSpecial) = pGrid%f2c(1,faceIndx)
                              
          iCellsSpecial = iCellsSpecial + 1          
          pGrid%cellsSpecial(iCellsSpecial) = pGrid%f2c(2,faceIndx) 
          
          WRITE(STDOUT,'(A,5X,A,1X,I8,1X,A,1X,I8)') SOLVER_NAME, & 
                'Added cells:',pGrid%f2c(1,faceIndx),'and', & 
                pGrid%f2c(2,faceIndx)         
        ELSE 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! faceIndx        
      
! --- cells meeting at vertex --------------------------------------------------
 
      CASE ( 'v' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter vertex index:'
        READ(STDIN,*,IOSTAT=errorFlag) vertIndx        
      
        IF ( errorFlag /= ERR_NONE ) THEN 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag               
      
        IF ( vertIndx > 0 .AND. vertIndx <= pGrid%nVertTot ) THEN 

! ------- Tetrahedra

          IF ( pGrid%nTetsTot > 0 ) THEN 
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Tetrahedra...'            

            nVertPerCell = 4

            DO icl = 1,pGrid%nTetsTot
              v(1:nVertPerCell) = pGrid%tet2v(1:nVertPerCell,icl)

              CALL QuickSortInteger(v,nVertPerCell)
              CALL BinarySearchInteger(v,nVertPerCell,vertIndx,iloc)

              IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                  CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                END IF ! iCellsSpecial              
              
                iCellsSpecial = iCellsSpecial + 1          
                pGrid%cellsSpecial(iCellsSpecial) = pGrid%tet2CellGlob(icl)

                WRITE(STDOUT,'(A,7X,A,1X,I8)') SOLVER_NAME,'Added cell:', & 
                                               pGrid%tet2CellGlob(icl)                                 
              END IF ! iloc            
            END DO ! icl        
          END IF ! pGrid%nTetsTot
        
! ------- Hexahedra        
        
          IF ( pGrid%nHexsTot > 0 ) THEN 
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Hexahedra...'             
        
            nVertPerCell = 8        
        
            DO icl = 1,pGrid%nHexsTot
              v(1:nVertPerCell) = pGrid%hex2v(1:nVertPerCell,icl)

              CALL QuickSortInteger(v,nVertPerCell)
              CALL BinarySearchInteger(v,nVertPerCell,vertIndx,iloc)

              IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                  CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                END IF ! iCellsSpecial              

                iCellsSpecial = iCellsSpecial + 1          
                pGrid%cellsSpecial(iCellsSpecial) = pGrid%hex2CellGlob(icl)

                WRITE(STDOUT,'(A,7X,A,1X,I8)') SOLVER_NAME,'Added cell:', & 
                                               pGrid%hex2CellGlob(icl)                                 
              END IF ! iloc            
            END DO ! icl  
          END IF ! pGrid%nPyrsTot                
        
! ------- Prisms        
        
          IF ( pGrid%nPrisTot > 0 ) THEN 
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Prisms...'          

            nVertPerCell = 6

            DO icl = 1,pGrid%nPrisTot
              v(1:nVertPerCell) = pGrid%pri2v(1:nVertPerCell,icl)

              CALL QuickSortInteger(v,nVertPerCell)
              CALL BinarySearchInteger(v,nVertPerCell,vertIndx,iloc)

              IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                  CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                END IF ! iCellsSpecial              

                iCellsSpecial = iCellsSpecial + 1          
                pGrid%cellsSpecial(iCellsSpecial) = pGrid%pri2CellGlob(icl)

                WRITE(STDOUT,'(A,7X,A,1X,I8)') SOLVER_NAME,'Added cell:', & 
                                               pGrid%pri2CellGlob(icl)                                 
              END IF ! iloc            
            END DO ! icl      
          END IF ! pGrid%nPyrsTot                 
        
! ------- Pyramids        
        
          IF ( pGrid%nPyrsTot > 0 ) THEN 
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Pyramids...'          
        
            nVertPerCell = 5        
        
            DO icl = 1,pGrid%nPyrsTot
              v(1:nVertPerCell) = pGrid%pyr2v(1:nVertPerCell,icl)

              CALL QuickSortInteger(v,nVertPerCell)
              CALL BinarySearchInteger(v,nVertPerCell,vertIndx,iloc)

              IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                  CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                END IF ! iCellsSpecial              

                iCellsSpecial = iCellsSpecial + 1          
                pGrid%cellsSpecial(iCellsSpecial) = pGrid%pyr2CellGlob(icl)

                WRITE(STDOUT,'(A,7X,A,1X,I8)') SOLVER_NAME,'Added cell:', & 
                                               pGrid%pyr2CellGlob(icl)                                 
              END IF ! iloc            
            END DO ! icl         
          END IF ! pGrid%nPyrsTot
        
        ELSE 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! vertIndx          
      
! --- quit ---------------------------------------------------------------------      
      
      CASE ( 'q' )
        EXIT
        
! --- default ------------------------------------------------------------------        
        
      CASE DEFAULT 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
        CYCLE 
    END SELECT  
  END DO ! <empty> 

! ==============================================================================
! Set number of special cells
! ==============================================================================

  pGrid%nCellsSpecial = iCellsSpecial

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting special cells done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetSpecialCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetSpecialCells.F90,v $
! Revision 1.3  2008/12/06 08:44:54  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2003/04/01 17:02:46  haselbac
! Initial revision
!
! Revision 1.2  2003/03/20 20:07:19  haselbac
! Modified RegFun call to avoid probs with
! long 'RFLU_GetSpecialCells.F90' names
!
! Revision 1.1  2003/03/15 19:16:54  haselbac
! Initial revision
!
!******************************************************************************








