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
! Purpose: Pick special cells.
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
! ******************************************************************************
!
! $Id: RFLU_PickSpecialCells.F90,v 1.11 2008/12/06 08:45:04 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PickSpecialCells(pRegion)

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

  CHARACTER :: infoType,stencilType
  CHARACTER(CHRLEN) :: RCSIdentString  
  INTEGER :: cellIndx,errorFlag,faceIndx,fnDir,fnDirEnd,iCellsSpecial,icl, &
             iloc,iPatch,nVertPerCell,patchIndx,vertIndx
  INTEGER :: v(8)
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PickSpecialCells.F90,v $ $Revision: 1.11 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PickSpecialCells', &
                        'RFLU_PickSpecialCells.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking special cells...'
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
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'b - cell adjacent to boundary face'   
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'c - single cell'
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'f - cells adjacent to interior face'
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'s - stencil members'
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'v - cells meeting at vertex'      
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'q - quit'

! ******************************************************************************
! Set up infinite loop
! ******************************************************************************

  DO
  
! ==============================================================================
!   Enter information type
! ==============================================================================
   
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information type:'
    READ(STDIN,'(A)') infoType
  
    SELECT CASE ( infoType )

! ------------------------------------------------------------------------------    
!     Cell adjacent to boundary face
! ------------------------------------------------------------------------------

      CASE ( 'b' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter patch index:'
        READ(STDIN,*,IOSTAT=errorFlag) patchIndx

        IF ( errorFlag /= ERR_NONE ) THEN 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag        
          
        IF ( patchIndx > 0 .AND. patchIndx <= pGrid%nPatches ) THEN
          pPatch => pRegion%patches(patchIndx)
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter face index:'
          READ(STDIN,*,IOSTAT=errorFlag) faceIndx           

          IF ( errorFlag /= ERR_NONE ) THEN
            global%warnCounter = global%warnCounter + 1           
           
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
            global%warnCounter = global%warnCounter + 1           
          
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                     '*** WARNING *** Invalid input.' 
            CYCLE 
          END IF ! faceIndx        
        ELSE 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! patchIndx
    
! ------------------------------------------------------------------------------
!     Single cell     
! ------------------------------------------------------------------------------
    
      CASE ( 'c' )
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter cell index:'
        READ(STDIN,*,IOSTAT=errorFlag) cellIndx
        
        IF ( errorFlag /= ERR_NONE ) THEN
          global%warnCounter = global%warnCounter + 1         
         
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
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! cellIndx

! ------------------------------------------------------------------------------
!     Cells adjacent to interior face     
! ------------------------------------------------------------------------------

      CASE ( 'f' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter interior face index:'
        READ(STDIN,*,IOSTAT=errorFlag) faceIndx
 
        IF ( errorFlag /= ERR_NONE ) THEN
          global%warnCounter = global%warnCounter + 1         
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag         
        
        IF ( faceIndx > 0 .AND. faceIndx <= pGrid%nFacesTot ) THEN
          IF ( iCellsSpecial == NCELLS_SPECIAL_MAX-1 ) THEN ! NOTE '-1'
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
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! faceIndx        
      
! ------------------------------------------------------------------------------
!     Stencil members 
! ------------------------------------------------------------------------------

      CASE ( 's' ) 
        WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter type of stencil:' 
        WRITE(STDOUT,'(A,9X,A)') SOLVER_NAME,'b - boundary-face stencil'
        WRITE(STDOUT,'(A,9X,A)') SOLVER_NAME,'c - cell stencil'
        WRITE(STDOUT,'(A,9X,A)') SOLVER_NAME,'f - face stencil'
        WRITE(STDOUT,'(A,9X,A)') SOLVER_NAME,'v - vertex stencil' 
        READ(STDIN,'(A)') stencilType       
       
        SELECT CASE ( stencilType )
        
! ------- Boundary-face stencil ------------------------------------------------

          CASE ( 'b' ) 
            WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter patch index:'
            READ(STDIN,*,IOSTAT=errorFlag) iPatch
            
            IF ( iPatch < 1 .OR. iPatch > pGrid%nPatches ) THEN 
              global%warnCounter = global%warnCounter + 1         

              WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                       '*** WARNING *** Invalid input.' 
              CYCLE               
            END IF ! iPatch
            
            pPatch => pRegion%patches(iPatch)
            
            WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter face index:'
            READ(STDIN,*,IOSTAT=errorFlag) faceIndx                        
        
            IF ( faceIndx > 0 .AND. faceIndx <= pPatch%nBFaces ) THEN            
              WRITE(STDOUT,'(A,9X,A,1X,I3)') SOLVER_NAME, &  
                'Number of stencil members:',pPatch%bf2cs(faceIndx)%nCellMembs
                        
              DO icl = 1,pPatch%bf2cs(faceIndx)%nCellMembs       
                iCellsSpecial = iCellsSpecial + 1
                
                IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                  CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                END IF ! iCellsSpecial                         
                
                pGrid%cellsSpecial(iCellsSpecial) = &
                  pPatch%bf2cs(faceIndx)%cellMembs(icl)
              END DO ! icl
            ELSE 
              global%warnCounter = global%warnCounter + 1         

              WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                       '*** WARNING *** Invalid input.' 
              CYCLE 
            END IF ! cellIndx            
        
! ------- Cell stencil ---------------------------------------------------------        
        
          CASE ( 'c' ) 
            SELECT CASE ( pRegion%mixtInput%stencilDimensCells ) 
              CASE ( 1 ) 
                IF ( ASSOCIATED(pGrid%c2cs1D) .EQV. .FALSE. ) THEN 
                  global%warnCounter = global%warnCounter + 1  

                  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME, & 
                                           '*** WARNING *** Stencil not built.' 
                  CYCLE                     
                END IF ! ASSOCIATED          

                SELECT CASE ( pRegion%mixtInput%dimens ) 
                  CASE ( 1 ) 
                    fnDirEnd = 1
                  CASE ( 2 ) 
                    fnDirEnd = 2
                  CASE ( 3 ) 
                    fnDirEnd = 3
                  CASE DEFAULT ! Defensive coding
                    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                END SELECT ! pRegion%mixtInput%dimensCells

                WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter cell index:'
                READ(STDIN,*,IOSTAT=errorFlag) cellIndx
                
                WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter direction:'
                READ(STDIN,*,IOSTAT=errorFlag) fnDir

                IF ( fnDir < 1 .OR. fnDir > fnDirEnd ) THEN                  
                  global%warnCounter = global%warnCounter + 1         

                  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                           '*** WARNING *** Invalid input.' 
                  CYCLE                 
                END IF ! fnDir

                IF ( cellIndx > 0 .AND. cellIndx <= pGrid%nCellsTot ) THEN
                  IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                    CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                  END IF ! iCellsSpecial         

                  iCellsSpecial = iCellsSpecial + 1          
                  pGrid%cellsSpecial(iCellsSpecial) = cellIndx              

                  WRITE(STDOUT,'(A,9X,A,1X,I3)') SOLVER_NAME, &  
                    'Number of stencil members:', &
                    pGrid%c2cs1D(fnDir,cellIndx)%nCellMembs  

                  DO icl = 1,pGrid%c2cs1D(fnDir,cellIndx)%nCellMembs          
                    iCellsSpecial = iCellsSpecial + 1

                    IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                      CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                    END IF ! iCellsSpecial                         

                    pGrid%cellsSpecial(iCellsSpecial) = &
                      pGrid%c2cs1D(fnDir,cellIndx)%cellMembs(icl)
                  END DO ! icl
                ELSE 
                  global%warnCounter = global%warnCounter + 1         

                  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                           '*** WARNING *** Invalid input.' 
                  CYCLE 
                END IF ! cellIndx
              CASE ( 2,3 )          
                IF ( ASSOCIATED(pGrid%c2cs) .EQV. .FALSE. ) THEN 
                  global%warnCounter = global%warnCounter + 1  

                  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME, & 
                                           '*** WARNING *** Stencil not built.' 
                  CYCLE                     
                END IF ! ASSOCIATED          

                WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter cell index:'
                READ(STDIN,*,IOSTAT=errorFlag) cellIndx

                IF ( cellIndx > 0 .AND. cellIndx <= pGrid%nCellsTot ) THEN
                  IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                    CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                  END IF ! iCellsSpecial         

                  iCellsSpecial = iCellsSpecial + 1          
                  pGrid%cellsSpecial(iCellsSpecial) = cellIndx              

                  WRITE(STDOUT,'(A,9X,A,1X,I3)') SOLVER_NAME, &  
                    'Number of stencil members:',pGrid%c2cs(cellIndx)%nCellMembs  

                  DO icl = 1,pGrid%c2cs(cellIndx)%nCellMembs          
                    iCellsSpecial = iCellsSpecial + 1

                    IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                      CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                    END IF ! iCellsSpecial                         

                    pGrid%cellsSpecial(iCellsSpecial) = &
                      pGrid%c2cs(cellIndx)%cellMembs(icl)
                  END DO ! icl
                ELSE 
                  global%warnCounter = global%warnCounter + 1         

                  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                           '*** WARNING *** Invalid input.' 
                  CYCLE 
                END IF ! cellIndx
              CASE DEFAULT
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! pRegion%mixtInput%stencilDimensCells            
            
! ------- Face stencil ---------------------------------------------------------        

          CASE ( 'f' )
            IF ( ASSOCIATED(pGrid%f2cs) .EQV. .FALSE. ) THEN 
              global%warnCounter = global%warnCounter + 1  
         
              WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME, & 
                                       '*** WARNING *** Stencil not built.' 
              CYCLE                     
            END IF ! ASSOCIATED          
           
            WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter face index:'
            READ(STDIN,*,IOSTAT=errorFlag) faceIndx
            
            IF ( faceIndx > 0 .AND. faceIndx <= pGrid%nFacesTot ) THEN            
              WRITE(STDOUT,'(A,9X,A,1X,I3)') SOLVER_NAME, &  
                'Number of stencil members:',pGrid%f2cs(faceIndx)%nCellMembs
                        
              DO icl = 1,pGrid%f2cs(faceIndx)%nCellMembs          
                iCellsSpecial = iCellsSpecial + 1
                
                IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                  CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                END IF ! iCellsSpecial                         
                
                pGrid%cellsSpecial(iCellsSpecial) = &
                  pGrid%f2cs(faceIndx)%cellMembs(icl)
              END DO ! icl
            ELSE 
              global%warnCounter = global%warnCounter + 1         

              WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                       '*** WARNING *** Invalid input.' 
              CYCLE 
            END IF ! cellIndx

! ------- Vertex stencil --------------------------------------------------------        
          
          CASE ( 'v' ) 
            IF ( ASSOCIATED(pGrid%v2cs) .EQV. .FALSE. ) THEN 
              global%warnCounter = global%warnCounter + 1  
         
              WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME, & 
                                       '*** WARNING *** Stencil not built.' 
              CYCLE                     
            END IF ! ASSOCIATED          
          
            WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'Enter vertex index:'
            READ(STDIN,*,IOSTAT=errorFlag) vertIndx
            
            IF ( vertIndx > 0 .AND. vertIndx <= pGrid%nVertTot ) THEN
              WRITE(STDOUT,'(A,9X,A,1X,I3)') SOLVER_NAME, &  
                'Number of stencil members:',pGrid%v2cs(vertIndx)%nCellMembs 

              DO icl = 1,pGrid%v2cs(vertIndx)%nCellMembs          
                iCellsSpecial = iCellsSpecial + 1
                
                IF ( iCellsSpecial == NCELLS_SPECIAL_MAX ) THEN 
                  CALL ErrorStop(global,ERR_NCELLS_SPECIAL_MAX,__LINE__)
                END IF ! iCellsSpecial                         
                
                pGrid%cellsSpecial(iCellsSpecial) = &
                  pGrid%v2cs(vertIndx)%cellMembs(icl)
              END DO ! icl
            ELSE 
              global%warnCounter = global%warnCounter + 1         

              WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                       '*** WARNING *** Invalid input.' 
              CYCLE 
            END IF ! cellIndx          
          
! ------- Default --------------------------------------------------------------          
          
          CASE DEFAULT
            global%warnCounter = global%warnCounter + 1 
      
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                                     '*** WARNING *** Invalid input.' 
            CYCLE         
        END SELECT ! stencilType 
      
! ------------------------------------------------------------------------------
!     Cells meeting at vertex 
! ------------------------------------------------------------------------------
 
      CASE ( 'v' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter vertex index:'
        READ(STDIN,*,IOSTAT=errorFlag) vertIndx        
      
        IF ( errorFlag /= ERR_NONE ) THEN
          global%warnCounter = global%warnCounter + 1         
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag               
      
        IF ( vertIndx > 0 .AND. vertIndx <= pGrid%nVertTot ) THEN 

! ------- Tetrahedra -----------------------------------------------------------

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
        
! ------- Hexahedra ------------------------------------------------------------       
        
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
        
! ------- Prisms ---------------------------------------------------------------       
        
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
        
! ------- Pyramids -------------------------------------------------------------       
        
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
          global%warnCounter = global%warnCounter + 1 
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! vertIndx          
      
! ------------------------------------------------------------------------------
!     Quit       
! ------------------------------------------------------------------------------
      
      CASE ( 'q' )
        EXIT
        
! ------------------------------------------------------------------------------
!     Default        
! ------------------------------------------------------------------------------
        
      CASE DEFAULT 
        global%warnCounter = global%warnCounter + 1 
      
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
        CYCLE 
    END SELECT  
  END DO ! <empty> 

! ******************************************************************************
! Set number of special cells
! ******************************************************************************

  pGrid%nCellsSpecial = iCellsSpecial

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking special cells done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PickSpecialCells

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PickSpecialCells.F90,v $
! Revision 1.11  2008/12/06 08:45:04  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:18:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2007/02/27 13:18:17  haselbac
! Enabled 1d computations
!
! Revision 1.8  2006/04/07 14:56:41  haselbac
! Adapted to new stencilDimens param
!
! Revision 1.7  2006/01/06 22:18:16  haselbac
! Added treatment of 1d stencils
!
! Revision 1.6  2005/01/10 19:37:42  haselbac
! Added capability of picking boundary-face stencils
!
! Revision 1.5  2004/12/27 23:33:17  haselbac
! Added writing of number of stencil members
!
! Revision 1.4  2004/10/19 19:30:09  haselbac
! Added checks for existence of stencils
!
! Revision 1.3  2004/02/13 03:02:10  haselbac
! Added stencils to selection
!
! Revision 1.2  2003/07/22 02:08:33  haselbac
! Added global%warnCounter
!
! Revision 1.1.1.1  2003/06/04 22:31:20  haselbac
! Initial revision
!
! ******************************************************************************







