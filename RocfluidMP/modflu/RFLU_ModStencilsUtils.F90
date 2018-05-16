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
! Purpose: Suite of utility routines to construct stencils.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModStencilsUtils.F90,v 1.14 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModStencilsUtils

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModStencil, ONLY: t_stencil
  USE ModSortSearch
  USE ModMPI
  
  USE RFLU_ModCellFaceEdgeInfo  
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_AddBFaces, & 
            RFLU_AddCellLayer, &  
            RFLU_AddCellLayer_1D, &
            RFLU_AddCellLayer_1D_G, &           
            RFLU_AddFaceVertNeighbs, &
            RFLU_ComputeStencilSize, &            
            RFLU_ComputeStencilWeights, &             
            RFLU_SortBFaces

       
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModStencilsUtils.F90,v $ $Revision: 1.14 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  


! ******************************************************************************
!
! Purpose: Extend basic stencil by adding boundary faces.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   nBFaceMembsMaxTemp  Maximum allowed number of boundary faces in stencil
!   nCellMembs          Number of cell members in stencil
!   cellMembs           List of cell members of stencil
!
! Output:
!   nBFaceMembs         Number of boundary faces in stencil
!   bFaceMembs          List of boundary faces in stencil
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_AddBFaces(pRegion,nBFaceMembsMaxTemp,nCellMembs,cellMembs, &
                            nBFaceMembs,bFaceMembs)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nBFaceMembsMaxTemp,nCellMembs
    INTEGER, INTENT(OUT) :: nBFaceMembs
    INTEGER, INTENT(IN) :: cellMembs(nCellMembs)
    INTEGER, INTENT(OUT) :: bFaceMembs(2,nBFaceMembsMaxTemp)                 
    TYPE(t_region), POINTER :: pRegion      

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifl,iloc,iPatch,isg,isl      
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid          
    TYPE(t_patch), POINTER :: pPatch      

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    global => pRegion%global

    pGrid => pRegion%grid

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_AddBFaces',&
  'RFLU_ModStencilsUtils.F90')     

! ******************************************************************************
!   Build sorted bf2c lists (and sorting keys) so can be searched later on
! ******************************************************************************        

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( (pPatch%bcType == BC_NOSLIPWALL_HFLUX) .OR. & 
           (pPatch%bcType == BC_NOSLIPWALL_TEMP ) .OR. & 
           (pPatch%bcType == BC_SLIPWALL        ) .OR. & 
           (pPatch%bcType == BC_INJECTION       ) ) THEN 
        IF ( pPatch%nBFaces > 0 ) THEN 
          ALLOCATE(pPatch%bf2cSorted(pPatch%nBFaces),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2cSorted')
          END IF ! global%error

          ALLOCATE(pPatch%bf2cSortedKeys(pPatch%nBFaces),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2cSortedKeys')
          END IF ! global%error

          DO ifl = 1,pPatch%nBFaces
            pPatch%bf2cSorted(ifl) = pPatch%bf2c(ifl)
          END DO ! ifl

          CALL QuickSortInteger(pPatch%bf2cSorted,pPatch%nBFaces)

          DO ifl = 1,pPatch%nBFaces
            CALL BinarySearchInteger(pPatch%bf2cSorted,pPatch%nBFaces, &
                                     pPatch%bf2c(ifl),iloc)

            IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
              pPatch%bf2cSortedKeys(iloc) = ifl
            ELSE 
              CALL ErrorStop(global,ERR_BINARY_SEARCH,__LINE__)
            END IF ! iloc                                   
          END DO ! ifl 
        END IF ! pPatch%nBFaces       
      END IF ! pPatch%bcType        
    END DO ! iPatch 

! ******************************************************************************
!   Initialize list
! ******************************************************************************        

    nBFaceMembs = 0

    DO isl = 1,nBFaceMembsMaxTemp
      bFaceMembs(1,isl) = CRAZY_VALUE_INT  
      bFaceMembs(2,isl) = CRAZY_VALUE_INT                 
    END DO ! isl

! ******************************************************************************
!   Find boundary cells in stencil
! ******************************************************************************        

    islLoop: DO isl = 1,nCellMembs
      isg = cellMembs(isl)

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( (pPatch%bcType == BC_NOSLIPWALL_HFLUX) .OR. & 
             (pPatch%bcType == BC_NOSLIPWALL_TEMP ) .OR. & 
             (pPatch%bcType == BC_SLIPWALL        ) .OR. & 
             (pPatch%bcType == BC_INJECTION       ) ) THEN 
          IF ( pPatch%nBFaces > 0 ) THEN 
            CALL BinarySearchInteger(pPatch%bf2cSorted,pPatch%nBFaces,isg, & 
                                     iloc)

            IF ( iloc /= ELEMENT_NOT_FOUND ) THEN
              IF ( nBFaceMembs < nBFaceMembsMaxTemp ) THEN  
                nBFaceMembs = nBFaceMembs + 1

                bFaceMembs(1,nBFaceMembs) = iPatch
                bFaceMembs(2,nBFaceMembs) = pPatch%bf2cSortedKeys(iloc)
              ELSE 
                EXIT islLoop
              END IF ! nBFaceMembs
            END IF ! iloc
          END IF ! pPatch%nBFaces 
        END IF ! pPatch%bcType                             
      END DO ! iPatch                    
    END DO islLoop          

! ******************************************************************************
!   Deallocate memory for sorted bf2c lists
! ******************************************************************************        

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( (pPatch%bcType == BC_NOSLIPWALL_HFLUX) .OR. & 
           (pPatch%bcType == BC_NOSLIPWALL_TEMP ) .OR. & 
           (pPatch%bcType == BC_SLIPWALL        ) .OR. & 
           (pPatch%bcType == BC_INJECTION       ) ) THEN
        IF ( pPatch%nBFaces > 0 ) THEN 
          DEALLOCATE(pPatch%bf2cSorted,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2cSorted')
          END IF ! global%error

          DEALLOCATE(pPatch%bf2cSortedKeys,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, &
                           'pPatch%bf2cSortedKeys')
          END IF ! global%error 
        END IF ! pPatch%nBFaces
      END IF ! pPatch%bcType       
    END DO ! iPatch 

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_AddBFaces
  
  
  





! ******************************************************************************
!
! Purpose: Extend basic stencil in only ONE direction
!
! Description: Loop over cells in the last layer of stencil and add cells
!  whose normal vector is aligned along a particular direction, i.e. x, y, z
!  Not just add vertex-neighbors because we only want to extend the stencil 
!  in one direction
!
! Input:
!   global              Pointer to global data
!   pGrid               Pointer to grid
!   stencilSizeMax      Maximum allowed size of stencil
!   ixg                 Index of cell for which cells are added
!   degr                Degree of stencil before addition of cell layer
!   x2csBeg             Beginning index for the last layer of cells in stencil
!   x2csEnd             Ending index for the last layer of cells in stencil
!   x2cs                Stencil before addition of cell layer
!   fnDir               Direction of stencil
!
! Output:
!   degr                Degree of stencil after addition of cell layer
!   x2cs                Stencil after addition of cell layer
!
! Notes: None.
!
! ******************************************************************************
 
  SUBROUTINE RFLU_AddCellLayer_1D(global,pGrid,stencilSizeMax,ixg,degr, &
                                  x2csBeg,x2csEnd,x2cs,fnDir)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: fnDir,ixg,stencilSizeMax,x2csBeg,x2csEnd
    INTEGER, INTENT(INOUT) :: degr
    INTEGER, INTENT(INOUT) :: x2cs(stencilSizeMax)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: c1,c2,degrHelp,icg,icg2,icl,ifg,ifl,iloc,iPatch,isl,iv2c,ivg, &
               ivl,nFaces
    INTEGER :: x2csHelp(stencilSizeMax),x2csSort(stencilSizeMax)
    REAL(RFREAL) :: fn

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_AddCellLayer_1D',&
  'RFLU_ModStencilsUtils.F90')

! ******************************************************************************
!   Add layer of cells to existing stencil in 1D, store in separate array
! ******************************************************************************

! ==============================================================================      
!   Initialize
! ==============================================================================

    degrHelp = 0 ! renamed because 'degr' is already used

    DO isl = 1,stencilSizeMax
      x2csHelp(isl) = 0 ! store new layer's members
      x2csSort(isl) = 0 ! duplicated array for x2cs
    END DO ! isl

    DO isl = 1,degr
      x2csSort(isl) = x2cs(isl)
    END DO ! isl

    CALL QuickSortInteger(x2csSort(1:degr),degr)

! ==============================================================================
!   Loop over existing members in the last layer and find candidates for stencil 
! ==============================================================================

    outerLoop: DO isl = x2csBeg,x2csEnd
      icg = x2cs(isl)
      icl = pGrid%cellGlob2Loc(2,icg)

      nFaces = SIZE(pGrid%hex2f,2)

      DO ifl = 1,nFaces
        iPatch = pGrid%hex2f(1,ifl,icl)
        ifg    = pGrid%hex2f(2,ifl,icl)
      
        IF ( iPatch == 0 ) THEN ! Interior face
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          IF ( c1 /= ixg .AND. c2 /= ixg ) THEN
            IF ( c1 == icg ) THEN
              fn =  pGrid%fn(fnDir,ifg)
            ELSE IF ( c2 == icg ) THEN
              fn = -pGrid%fn(fnDir,ifg)
            ELSE ! defensive programming
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! c1

            IF ( ABS(fn) >= 0.999_RFREAL ) THEN 
              IF ( c1 == icg ) THEN 
                icg2 = c2
              ELSE IF ( c2 == icg ) THEN 
                icg2 = c1
              ELSE ! defensive programming
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
              END IF ! c1

              CALL BinarySearchInteger(x2csSort(1:degr),degr,icg2,iloc)
              
              IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                IF ( degrHelp > 0 ) THEN 
                  CALL BinarySearchInteger(x2csHelp(1:degrHelp),degrHelp, &
                                           icg2,iloc)
                  IF ( iloc == ELEMENT_NOT_FOUND ) THEN
                    IF ( degrHelp < stencilSizeMax ) THEN  
                      degrHelp = degrHelp + 1
                      x2csHelp(degrHelp) = icg2
                      IF ( degrHelp > 1 ) THEN
                        CALL QuickSortInteger(x2csHelp(1:degrHelp),degrHelp)
                      END IF ! degrHelp
                    ELSE
                      EXIT outerLoop
                    END IF ! degrHelp
                  END IF ! iloc
                ELSE ! First member                
                  degrHelp = degrHelp + 1
                  x2csHelp(degrHelp) = icg2
                END IF ! degrHelp
              END IF ! iloc
            END IF ! ABS(fn)
          END IF ! c1, c2                 
        ELSE IF ( iPatch > 0 ) THEN ! Boundary face  
! TO DO 
!       Currently do not add boundary faces for constrained reconstruction
! END TO DO 
        ELSE ! Defensive programming
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! iPatch
      END DO ! ifl
    END DO outerLoop

! ******************************************************************************
!   Merge new layer with existing stencil
! ******************************************************************************

    DO isl = 1,degrHelp
      IF ( degr < stencilSizeMax ) THEN 
        degr = degr + 1
        x2cs(degr) = x2csHelp(isl)
      END IF ! degr
    END DO ! isl    

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_AddCellLayer_1D






! ******************************************************************************
!
! Purpose: Extend basic stencil in ONE direction by adding vertex neighbours 
!   of cells based on geometry.
!
! Description: Loop over cells in last layer of stencil and add all cells
!  which share vertices with that layer of cells and are aligned with the
!  given direction.
!
! Input:
!   global              Pointer to global data
!   pGrid               Pointer to grid
!   stencilSizeMax      Maximum allowed size of stencil
!   ixg                 Index of cell for which cells are added
!   degr                Degree of stencil before addition of cell layer
!   x2csBeg             Beginning index for last layer of cells in stencil
!   x2csEnd             Ending index for last layer of cells in stencil
!   x2cs                Stencil before addition of cell layer
!   rc			Location for which stencil is constructed
!   dir			Direction in which stencil is to be constructed
!
! Output:
!   degr                Degree of stencil after addition of cell layer
!   x2cs                Stencil after addition of cell layer
!
! Notes: None.
!
! ******************************************************************************
 
  SUBROUTINE RFLU_AddCellLayer_1D_G(global,pGrid,stencilSizeMax,ixg,degr, &
                                    x2csBeg,x2csEnd,x2cs,rc,dir)

    USE RFLU_ModGeometryTools, ONLY: RFLU_TestVectorCartAxisAligned

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: dir,ixg,stencilSizeMax,x2csBeg,x2csEnd
    INTEGER, INTENT(INOUT) :: degr
    INTEGER, INTENT(INOUT) :: x2cs(stencilSizeMax)
    REAL(RFREAL), INTENT(IN) :: rc(XCOORD:ZCOORD)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: alignFlag
    INTEGER :: degrHelp,icg,icg2,icl,ict,iloc,isl,ivg,ivl,iv2c
    INTEGER :: x2csHelp(stencilSizeMax),x2csSort(stencilSizeMax)
    REAL(RFREAL) :: dr(XCOORD:ZCOORD)

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_AddCellLayer_1D_G',&
  'RFLU_ModStencilsUtils.F90')

! ******************************************************************************
!   Add layer of cells to existing stencil, store in separate array
! ******************************************************************************

! ==============================================================================      
!   Initialize
! ==============================================================================

    degrHelp = 0

    DO isl = 1,stencilSizeMax
      x2csHelp(isl) = 0
      x2csSort(isl) = 0
    END DO ! isl

    DO isl = 1,degr
      x2csSort(isl) = x2cs(isl)
    END DO ! isl

    CALL QuickSortInteger(x2csSort(1:degr),degr)

! ==============================================================================
!   Loop over existing members in last layer and find candidates for stencil 
! ==============================================================================

    outerLoop: DO isl = x2csBeg,x2csEnd                       
      icg = x2cs(isl)
      ict = RFLU_GetGlobalCellType(global,pGrid,icg)
      icl = pGrid%cellGlob2Loc(2,icg)

      SELECT CASE ( ict ) 

! ------------------------------------------------------------------------------
!       Hexahedra
! ------------------------------------------------------------------------------                    

        CASE ( CELL_TYPE_HEX ) 
          DO ivl = 1,8
            ivg = pGrid%hex2v(ivl,icl)

            DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
              icg2 = pGrid%v2c(iv2c)

              IF ( icg2 /= ixg ) THEN 
                CALL BinarySearchInteger(x2csSort(1:degr),degr,icg2,iloc)

                IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                  IF ( degrHelp > 0 ) THEN 
                    CALL BinarySearchInteger(x2csHelp(1:degrHelp),degrHelp, &
                                             icg2,iloc)
                    IF ( iloc == ELEMENT_NOT_FOUND ) THEN
                      IF ( degrHelp < stencilSizeMax ) THEN
                        dr(XCOORD) = pGrid%cofg(XCOORD,icg2) - rc(XCOORD)
                        dr(YCOORD) = pGrid%cofg(YCOORD,icg2) - rc(YCOORD)
                        dr(ZCOORD) = pGrid%cofg(ZCOORD,icg2) - rc(ZCOORD)
                        
                        alignFlag = RFLU_TestVectorCartAxisAligned(global,dr,dir)
                        
                        IF ( alignFlag .EQV. .TRUE. ) THEN                                                 
                          degrHelp = degrHelp + 1
                          x2csHelp(degrHelp) = icg2

                          IF ( degrHelp > 1 ) THEN 
                            CALL QuickSortInteger(x2csHelp(1:degrHelp),degrHelp)
                          END IF ! degrHelp
                        END IF ! alignFlag
                      ELSE
                        EXIT outerLoop
                      END IF ! degrHelp                       
                    END IF ! iloc                          
                  ELSE
                    dr(XCOORD) = pGrid%cofg(XCOORD,icg2) - rc(XCOORD)
                    dr(YCOORD) = pGrid%cofg(YCOORD,icg2) - rc(YCOORD)
                    dr(ZCOORD) = pGrid%cofg(ZCOORD,icg2) - rc(ZCOORD)

                    alignFlag = RFLU_TestVectorCartAxisAligned(global,dr,dir)

                    IF ( alignFlag .EQV. .TRUE. ) THEN
                      degrHelp = degrHelp + 1
                      x2csHelp(degrHelp) = icg2
                    END IF ! alignFlag                       
                  END IF ! degrHelp
                END IF ! iloc
              END IF ! icg2

            END DO ! iv2c
          END DO ! ivl

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------            

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                
      END SELECT ! ict
    END DO outerLoop 

! ******************************************************************************
!   Merge new layer with existing stencil
! ******************************************************************************

    DO isl = 1,degrHelp
      IF ( degr < stencilSizeMax ) THEN 
        degr = degr + 1
        x2cs(degr) = x2csHelp(isl)
      END IF ! degr
    END DO ! isl    

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_AddCellLayer_1D_G



  
  

! ******************************************************************************
!
! Purpose: Extend basic stencil by adding vertex neighbours of cells.
!
! Description: Loop over cells in last layer of stencil and add all cells
!  which share vertices with that layer of cells.
!
! Input:
!   global              Pointer to global data
!   pGrid               Pointer to grid
!   stencilSizeMax      Maximum allowed size of stencil
!   ixg                 Index of cell for which cells are added
!   degr                Degree of stencil before addition of cell layer
!   x2csBeg             Beginning index for last layer of cells in stencil
!   x2csEnd             Ending index for last layer of cells in stencil
!   x2cs                Stencil before addition of cell layer
!
! Output:
!   degr                Degree of stencil after addition of cell layer
!   x2cs                Stencil after addition of cell layer
!
! Notes: None.
!
! ******************************************************************************
 
  SUBROUTINE RFLU_AddCellLayer(global,pGrid,stencilSizeMax,ixg,degr,x2csBeg, &
                               x2csEnd,x2cs)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: ixg,stencilSizeMax,x2csBeg,x2csEnd
    INTEGER, INTENT(INOUT) :: degr
    INTEGER, INTENT(INOUT) :: x2cs(stencilSizeMax)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: degrHelp,icg,icg2,icl,ict,iloc,isl,ivg,ivl,iv2c
    INTEGER :: x2csHelp(stencilSizeMax),x2csSort(stencilSizeMax)

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_AddCellLayer',&
  'RFLU_ModStencilsUtils.F90')

! ******************************************************************************
!   Add layer of cells to existing stencil, store in separate array
! ******************************************************************************

! ==============================================================================      
!   Initialize
! ==============================================================================

    degrHelp = 0

    DO isl = 1,stencilSizeMax
      x2csHelp(isl) = 0
      x2csSort(isl) = 0
    END DO ! isl

    DO isl = 1,degr
      x2csSort(isl) = x2cs(isl)
    END DO ! isl

    CALL QuickSortInteger(x2csSort(1:degr),degr)

! ==============================================================================
!   Loop over existing members in last layer and find candidates for stencil 
! ==============================================================================

    outerLoop: DO isl = x2csBeg,x2csEnd                       
      icg = x2cs(isl)
      ict = RFLU_GetGlobalCellType(global,pGrid,icg)
      icl = pGrid%cellGlob2Loc(2,icg)

      SELECT CASE ( ict ) 

! ------------------------------------------------------------------------------
!       Tetrahedra
! ------------------------------------------------------------------------------        

        CASE ( CELL_TYPE_TET ) 
          DO ivl = 1,4
            ivg = pGrid%tet2v(ivl,icl)

            DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
              icg2 = pGrid%v2c(iv2c)

              IF ( icg2 /= ixg ) THEN 
                CALL BinarySearchInteger(x2csSort(1:degr),degr,icg2,iloc)

                IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                  IF ( degrHelp > 0 ) THEN 
                    CALL BinarySearchInteger(x2csHelp(1:degrHelp),degrHelp, &
                                             icg2,iloc)

                    IF ( iloc == ELEMENT_NOT_FOUND ) THEN
                      IF ( degrHelp < stencilSizeMax ) THEN  
                        degrHelp = degrHelp + 1
                        x2csHelp(degrHelp) = icg2

                        IF ( degrHelp > 1 ) THEN 
                          CALL QuickSortInteger(x2csHelp(1:degrHelp),degrHelp)
                        END IF ! degrHelp
                      ELSE
                        EXIT outerLoop
                      END IF ! degrHelp                       
                    END IF ! iloc                         
                  ELSE
                    degrHelp = degrHelp + 1
                    x2csHelp(degrHelp) = icg2
                  END IF ! degrHelp
                END IF ! iloc
              END IF ! icg2 

            END DO ! iv2c
          END DO ! ivl

! ------------------------------------------------------------------------------
!       Hexahedra
! ------------------------------------------------------------------------------                    

        CASE ( CELL_TYPE_HEX ) 
          DO ivl = 1,8
            ivg = pGrid%hex2v(ivl,icl)

            DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
              icg2 = pGrid%v2c(iv2c)

              IF ( icg2 /= ixg ) THEN 
                CALL BinarySearchInteger(x2csSort(1:degr),degr,icg2,iloc)

                IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                  IF ( degrHelp > 0 ) THEN 
                    CALL BinarySearchInteger(x2csHelp(1:degrHelp),degrHelp, &
                                             icg2,iloc)
                    IF ( iloc == ELEMENT_NOT_FOUND ) THEN
                      IF ( degrHelp < stencilSizeMax ) THEN  
                        degrHelp = degrHelp + 1
                        x2csHelp(degrHelp) = icg2

                        IF ( degrHelp > 1 ) THEN 
                          CALL QuickSortInteger(x2csHelp(1:degrHelp),degrHelp)
                        END IF ! degrHelp
                      ELSE
                        EXIT outerLoop
                      END IF ! degrHelp                       
                    END IF ! iloc                          
                  ELSE
                    degrHelp = degrHelp + 1
                    x2csHelp(degrHelp) = icg2
                  END IF ! degrHelp
                END IF ! iloc
              END IF ! icg2

            END DO ! iv2c
          END DO ! ivl

! ------------------------------------------------------------------------------
!       Prisms
! ------------------------------------------------------------------------------                    

        CASE ( CELL_TYPE_PRI )           
          DO ivl = 1,6
            ivg = pGrid%pri2v(ivl,icl)

            DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
              icg2 = pGrid%v2c(iv2c)

               IF ( icg2 /= ixg ) THEN 
                CALL BinarySearchInteger(x2csSort(1:degr),degr,icg2,iloc)

                IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                  IF ( degrHelp > 0 ) THEN 
                    CALL BinarySearchInteger(x2csHelp(1:degrHelp),degrHelp, &
                                             icg2,iloc)

                    IF ( iloc == ELEMENT_NOT_FOUND ) THEN
                      IF ( degrHelp < stencilSizeMax ) THEN  
                        degrHelp = degrHelp + 1
                        x2csHelp(degrHelp) = icg2

                        IF ( degrHelp > 1 ) THEN 
                          CALL QuickSortInteger(x2csHelp(1:degrHelp),degrHelp)
                        END IF ! degrHelp
                      ELSE
                        EXIT outerLoop
                      END IF ! degrHelp                       
                    END IF ! iloc                          
                  ELSE
                    degrHelp = degrHelp + 1
                    x2csHelp(degrHelp) = icg2
                  END IF ! degrHelp
                END IF ! iloc
              END IF ! icg2

            END DO ! iv2c
          END DO ! ivl  

! ------------------------------------------------------------------------------
!       Pyramids
! ------------------------------------------------------------------------------

        CASE ( CELL_TYPE_PYR ) 
          DO ivl = 1,5
            ivg = pGrid%pyr2v(ivl,icl)

            DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
              icg2 = pGrid%v2c(iv2c)

              IF ( icg2 /= ixg ) THEN 
                CALL BinarySearchInteger(x2csSort(1:degr),degr,icg2,iloc)

                IF ( iloc == ELEMENT_NOT_FOUND ) THEN 
                  IF ( degrHelp > 0 ) THEN 
                    CALL BinarySearchInteger(x2csHelp(1:degrHelp),degrHelp, &
                                             icg2,iloc)

                    IF ( iloc == ELEMENT_NOT_FOUND ) THEN
                      IF ( degrHelp < stencilSizeMax ) THEN  
                        degrHelp = degrHelp + 1
                        x2csHelp(degrHelp) = icg2

                        IF ( degrHelp > 1 ) THEN 
                          CALL QuickSortInteger(x2csHelp(1:degrHelp),degrHelp)
                        END IF ! degrHelp
                      ELSE
                        EXIT outerLoop
                      END IF ! degrHelp                       
                    END IF ! iloc                           
                  ELSE
                    degrHelp = degrHelp + 1
                    x2csHelp(degrHelp) = icg2
                  END IF ! degrHelp
                END IF ! iloc
              END IF ! icg2

            END DO ! iv2c              
          END DO ! ivl

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------            

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                
      END SELECT ! ict
    END DO outerLoop 

! ******************************************************************************
!   Merge new layer with existing stencil
! ******************************************************************************

    DO isl = 1,degrHelp
      IF ( degr < stencilSizeMax ) THEN 
        degr = degr + 1
        x2cs(degr) = x2csHelp(isl)
      END IF ! degr
    END DO ! isl    

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_AddCellLayer






! ******************************************************************************
!
! Purpose: Construct initial stencil by adding cells sharing vertices which
!   define face. 
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pGrid               Pointer to grid
!   stencilSizeMax      Maximum allowed size of stencil
!   f2v                 List of vertices which define face
!
! Output:
!   degr                Number of cells sharing vertices defining face
!   x2cs                List of cells sharing vertices defining face
!
! Notes:
!   1. It is assumed that there are no members in the stencil when this 
!      routine is called. 
!
! ******************************************************************************

  SUBROUTINE RFLU_AddFaceVertNeighbs(global,pGrid,stencilSizeMax,f2v,degr, &
                                     x2cs)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: stencilSizeMax
    INTEGER, INTENT(IN) :: f2v(4)
    INTEGER, INTENT(INOUT) :: degr
    INTEGER, INTENT(INOUT) :: x2cs(stencilSizeMax)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,iloc,ivg,ivl,ivl2

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_AddFaceVertNeighbs',&
  'RFLU_ModStencilsUtils.F90')

! ******************************************************************************
!   Add layer of cells to existing stencil, store in separate array
! ******************************************************************************

! ==============================================================================      
!   Initialize
! ==============================================================================

    degr = 0

! ==============================================================================
!   Loop over vertices of face and add cells  
! ==============================================================================

    outerLoop: DO ivl = 1,4
      ivg = f2v(ivl)

      IF ( ivg /= VERT_NONE ) THEN 

        DO ivl2 = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
          icg = pGrid%v2c(ivl2)

          IF ( degr /= 0 ) THEN 
            CALL BinarySearchInteger(x2cs(1:degr),degr,icg,iloc)

            IF ( iloc == ELEMENT_NOT_FOUND ) THEN
              IF ( degr < stencilSizeMax ) THEN  
                degr = degr + 1
                x2cs(degr) = icg

                IF ( degr > 1 ) THEN 
                  CALL QuickSortInteger(x2cs(1:degr),degr)
                END IF ! degr
              ELSE
                EXIT outerLoop
              END IF ! degr
            END IF ! iloc
          ELSE 
            degr = degr + 1

            x2cs(degr) = icg
          END IF ! degr
        END DO ! ivl2          

      END IF ! ivg
    END DO outerLoop
                                       
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_AddFaceVertNeighbs






! *******************************************************************************
!
! Purpose: Compute stencil size for given order.
!
! Description: None.
!
! Input:
!   global      Global pointer
!   dimens      Dimensionality
!   factor      Safety factor
!   order       Polynomial order
!
! Output: 
!   RFLU_ComputeStencilSize     (Minimum) number of points in stencil
!
! Notes: 
!   1. Include additional support (+1) to allow for interpolation.
!   2. Make larger than necessary to allow for least-squares if order > 1.
!   3. If order = 0, do not include factor. 
!
! ******************************************************************************

  INTEGER FUNCTION RFLU_ComputeStencilSize(global,dimens,factor,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: dimens,factor,order
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeStencilSize',&
  'RFLU_ModStencilsUtils.F90')

! ******************************************************************************
!   Compute stencil size
! ******************************************************************************

    IF ( order > 0 ) THEN 
      SELECT CASE ( dimens ) 
        CASE ( 2 ) 
          RFLU_ComputeStencilSize = factor*(order + 1)*(order + 2)/2
        CASE ( 3 ) 
          RFLU_ComputeStencilSize = factor*(order + 1)*(order + 2)*(order + 3)/6        
        CASE DEFAULT
          CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
      END SELECT ! dimens
    ELSE 
      RFLU_ComputeStencilSize = 1
    END IF ! order

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END FUNCTION RFLU_ComputeStencilSize






! ******************************************************************************
!
! Purpose: Compute stencil weights.
!
! Description: Compute stencil weights in iterative fashion. If weights are
!   rejected, requested order is lowered and the weights are recomputed. Once
!   order is equal to zero, the weights are set instead of computed.
!
! Input:
!   global              Global pointer
!   dimens              Dimensionality
!   wtsMode             Mode of computing weights
!   scalMode            Mode of scaling
!   derivDegree         Degree of derivative (0 = interpolation)
!   order               Order of accuracy of approximation
!   nRows               Number of members in stencil
!   dr                  Relative position vectors of members in stencil
!
! Output: 
!   orderNominal        Order of accuracy of approximation
!   wts                 Stencil weights
!   sCount              Number of singular values
!
! Notes:
!   1. It is important to note that the order is the polynomial order of 
!      approximation, that is, order equal to 1 means a linear polynomial
!      which may be first- or second-order accurate...
!   2. Still need to check for singular systems because although initial 
!      stencils ought to be non-singular, they may become singular with grid
!      motion.
!   3. Although this is a routine to compute weights, it is not located in the
!      module RFLU_ComputeWeights, because it requires the function 
!      RFLU_ComputeStencilSize, and this leads to a circular dependency:
!      RFLU_ModWeights depends on RFLU_ModStencils (because of 
!      RFLU_ComputeStencilSize), and RFLU_ModStencils depends on 
!      RFLU_ModWeights (because of RFLU_ComputeStencilWeights)
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeStencilWeights(global,dimens,wtsMode,scalMode, &
                                        derivDegree,orderNominal,nRows,dr, &
                                        wts,sCount)

    USE ModTools, ONLY: CompFact

    USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: derivDegree,dimens,nRows,scalMode,wtsMode
    INTEGER, INTENT(IN) :: orderNominal
    INTEGER, INTENT(OUT), OPTIONAL :: sCount 
    REAL(RFREAL), INTENT(IN) :: dr(XCOORD:ZCOORD,nRows)
    REAL(RFREAL), INTENT(OUT) :: wts(:,:)
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: RCSIdentString
    INTEGER :: errorFlag,isl,j,nCols,order,p,q,r
    REAL(RFREAL) :: colScal,term,wtsNorm,wtsNormMax,wtsNormMin
    REAL(RFREAL) :: rowScal(nRows)
    REAL(RFREAL) :: drCopy(XCOORD:ZCOORD,nRows)
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,aInv

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeStencilWeights',&
  'RFLU_ModStencilsUtils.F90')

    wtsNormMin = 1.0_RFREAL 
    wtsNormMax = 10.0_RFREAL 

! ******************************************************************************
!   Check that stencil large enough, reduce order if necessary
! ******************************************************************************

    IF ( nRows < RFLU_ComputeStencilSize(global,dimens,1,orderNominal) ) THEN 
      SELECT CASE ( dimens )
        CASE ( 2 ) 
          SELECT CASE ( nRows )
            CASE ( 0:2 ) ! Constant
              order = 0 
            CASE ( 3:5 ) ! Linear
              order = 1
            CASE ( 6:9 ) ! Quadratic
              order = 2
            CASE ( 10:14 ) ! Cubic
              order = 3
            CASE DEFAULT ! More than cubic
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! nRows           
        CASE ( 3 ) 
          SELECT CASE ( nRows )
            CASE ( 0:3 ) ! Constant
              order = 0 
            CASE ( 4:9 ) ! Linear
              order = 1
            CASE ( 10:19 ) ! Quadratic
              order = 2
            CASE ( 20:34 ) ! Cubic
              order = 3
            CASE DEFAULT ! More than cubic
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! nRows
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! dimens
    ELSE 
      order = orderNominal 
    END IF ! nRowsEff

! ******************************************************************************
!   Row scaling
! ******************************************************************************

    IF ( scalMode == COMPWTS_SCAL_NONE ) THEN 
      DO isl = 1,nRows
        rowScal(isl) = 1.0_RFREAL
      END DO ! isl  
    ELSE IF ( scalMode == COMPWTS_SCAL_INVDIST ) THEN 
      DO isl = 1,nRows
        rowScal(isl) = 1.0_RFREAL/SQRT(dr(XCOORD,isl)*dr(XCOORD,isl) & 
                                     + dr(YCOORD,isl)*dr(YCOORD,isl) & 
                                     + dr(ZCOORD,isl)*dr(ZCOORD,isl))
      END DO ! isl
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! scalMode

! ******************************************************************************
!   Column scaling
! ******************************************************************************

    colScal = MAXVAL(ABS(dr(XCOORD:ZCOORD,1:nRows)))

    DO isl = 1,nRows
      drCopy(XCOORD,isl) = dr(XCOORD,isl)/colScal
      drCopy(YCOORD,isl) = dr(YCOORD,isl)/colScal
      drCopy(ZCOORD,isl) = dr(ZCOORD,isl)/colScal                
    END DO ! isl

! ******************************************************************************
!   Compute weights iteratively
! ******************************************************************************

    DO 

! ==============================================================================
!     Order greater than zero, so compute weights
! ==============================================================================

      IF ( order > 0 ) THEN 
        nCols = RFLU_ComputeStencilSize(global,dimens,1,order)        

! ------------------------------------------------------------------------------
!       Allocate temporary memory
! ------------------------------------------------------------------------------

        ALLOCATE(a(nRows,nCols),STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'a')
        END IF ! global%error      

        ALLOCATE(aInv(nCols,nRows),STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'aInv')
        END IF ! global%error      

! ------------------------------------------------------------------------------
!       Build matrix
! ------------------------------------------------------------------------------

        SELECT CASE ( dimens )
          CASE ( 2 ) 
            DO isl = 1,nRows
              a(isl,1) = rowScal(isl)

              j = 2

              DO p = 1,order
                DO q = 0,p
                  term = rowScal(isl)/(CompFact(p-q)*CompFact(q))                 
                  a(isl,j) = term*drCopy(XCOORD,isl)**(p-q) &
                                 *drCopy(YCOORD,isl)**q   
                  j = j + 1
                END DO ! q
              END DO ! p                    
            END DO ! isl             
          CASE ( 3 ) 
            DO isl = 1,nRows
              a(isl,1) = rowScal(isl)

              j = 2

              DO p = 1,order
                DO q = 0,p
                  DO r = 0,q
                    term = rowScal(isl)/(CompFact(p-q)*CompFact(q-r)*CompFact(r))                 
                    a(isl,j) = term*drCopy(XCOORD,isl)**(p-q) &
                                   *drCopy(YCOORD,isl)**(q-r) &
                                   *drCopy(ZCOORD,isl)**r 
                    j = j + 1
                  END DO ! r
                END DO ! q
              END DO ! p                    
            END DO ! isl
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! dimens    

! ------------------------------------------------------------------------------
!       Invert matrix
! ------------------------------------------------------------------------------

        CALL RFLU_InvertMatrixSVD(global,nRows,nCols,a,aInv,sCount)

! ------------------------------------------------------------------------------
!       Extract weights from inverse. NOTE appearance of scalings.
! ------------------------------------------------------------------------------

        IF ( derivDegree == DERIV_DEGREE_0 ) THEN 
          DO isl = 1,nRows
            wts(1,isl) = rowScal(isl)*aInv(1,isl) 
          END DO ! isl
        ELSE IF ( derivDegree == DERIV_DEGREE_1 ) THEN
          SELECT CASE ( dimens )
            CASE ( 2 ) 
              DO isl = 1,nRows
                wts(XCOORD,isl) = rowScal(isl)*aInv(2,isl)/colScal
                wts(YCOORD,isl) = rowScal(isl)*aInv(3,isl)/colScal
                wts(ZCOORD,isl) = 0.0_RFREAL
              END DO ! isl                
            CASE ( 3 ) 
              DO isl = 1,nRows
                wts(XCOORD,isl) = rowScal(isl)*aInv(2,isl)/colScal
                wts(YCOORD,isl) = rowScal(isl)*aInv(3,isl)/colScal
                wts(ZCOORD,isl) = rowScal(isl)*aInv(4,isl)/colScal 
              END DO ! isl   
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)   
          END SELECT ! dimens
        ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! derivDegree

! ------------------------------------------------------------------------------
!       Check whether weights make sense
! ------------------------------------------------------------------------------

        wtsNorm = 0.0_RFREAL

        DO isl = 1,nRows
          wtsNorm = wtsNorm + ABS(aInv(1,isl))
        END DO ! isl

! ------------------------------------------------------------------------------
!       Deallocate temporary memory
! ------------------------------------------------------------------------------

        DEALLOCATE(a,STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'a')
        END IF ! global%error

        DEALLOCATE(aInv,STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'aInv')
        END IF ! global%error

! ------------------------------------------------------------------------------
!       Decide whether to accept or reject weights
! ------------------------------------------------------------------------------

        IF ( sCount == 0 .AND. & 
             wtsNorm > wtsNormMin .AND. wtsNorm < wtsNormMax ) THEN 
          EXIT
        ELSE 
          IF ( wtsMode == COMPWTS_MODE_ADAPT ) THEN 
            order = order - 1                     
          ELSE 
            EXIT
          END IF ! wtsMode
        END IF ! errNorm     

! ==============================================================================
!     Order equal to zero, so set weights 
! ==============================================================================

      ELSE           
        IF ( derivDegree == DERIV_DEGREE_0 ) THEN 
          DO isl = 1,nRows
            wts(1,isl) = 1.0_RFREAL/REAL(nRows,KIND=RFREAL) 
          END DO ! isl
        ELSE IF ( derivDegree == DERIV_DEGREE_1 ) THEN 
          DO isl = 1,nRows
            wts(XCOORD,isl) = 0.0_RFREAL
            wts(YCOORD,isl) = 0.0_RFREAL
            wts(ZCOORD,isl) = 0.0_RFREAL
          END DO ! isl      
        ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! derivDegree         

        EXIT 
      END IF ! order                                
    END DO ! <empty>

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeStencilWeights







! ******************************************************************************
!
! Purpose: Sort boundary faces by increasing distance from location at which
!   gradient is reconstructed.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   xyz                 Coordinates at which gradient is reconstructed
!   nBFaceMembs         Number of boundary face stencil members
!   bFaceMembs          Boundary face stencil members
!
! Output:
!   bFaceMembs          Boundary face stencil members
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SortBFaces(pRegion,xyz,nBFaceMembs,bFaceMembs)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nBFaceMembs
    INTEGER, INTENT(INOUT) :: bFaceMembs(2,nBFaceMembs)
    REAL(RFREAL), INTENT(IN) :: xyz(XCOORD:ZCOORD) 
    TYPE(t_region), POINTER :: pRegion      

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: ifl,iPatch,isl
    INTEGER :: key(nBFaceMembs)
    INTEGER :: bFaceMembsBak(2,nBFaceMembs)      
    REAL(RFREAL) :: dist(nBFaceMembs)      
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid          
    TYPE(t_patch), POINTER :: pPatch      

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    global => pRegion%global
    pGrid  => pRegion%grid

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_SortBFaces',&
  'RFLU_ModStencilsUtils.F90')     

! ******************************************************************************
!   Back up data structure so can copy later on
! ******************************************************************************

    DO isl = 1,nBFaceMembs 
      bFaceMembsBak(1,isl) = bFaceMembs(1,isl)
      bFaceMembsBak(2,isl) = bFaceMembs(2,isl)
    END DO ! isl   

! ******************************************************************************
!   Compute squared distance and build key for sorting
! ******************************************************************************

    DO isl = 1,nBFaceMembs
      iPatch = bFaceMembs(1,isl)
      ifl    = bFaceMembs(2,isl)

      pPatch => pRegion%patches(iPatch)

      key(isl)  = isl
      dist(isl) = (pPatch%fc(XCOORD,ifl) - xyz(XCOORD))**2 & 
                + (pPatch%fc(YCOORD,ifl) - xyz(YCOORD))**2 & 
                + (pPatch%fc(ZCOORD,ifl) - xyz(ZCOORD))**2
    END DO ! isl    
  
! ******************************************************************************
!   Sort according to increasing distance
! ******************************************************************************

    CALL QuickSortRFREALInteger(dist,key,nBFaceMembs)  

    DO isl = 1,nBFaceMembs
      bFaceMembs(1,isl) = bFaceMembsBak(1,key(isl))
      bFaceMembs(2,isl) = bFaceMembsBak(2,key(isl))        
    END DO ! isl        

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)      

  END SUBROUTINE RFLU_SortBFaces





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModStencilsUtils


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModStencilsUtils.F90,v $
! Revision 1.14  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2007/03/06 18:08:37  haselbac
! Added function to add cell layer in 1d based on geometry (not on hex2f)
!
! Revision 1.11  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.10  2006/04/01 15:56:37  haselbac
! Bug fix: Only add bfaces if have actual faces
!
! Revision 1.9  2006/03/25 21:57:59  haselbac
! Changed bf2cSorted list to contain only actual entries (bcos of sype changes)
!
! Revision 1.8  2006/01/06 22:14:15  haselbac
! Added new routine for adding cell layer in 1d
!
! Revision 1.7  2005/10/09 13:27:18  haselbac
! Temporarily disable addition of v bfaces to stencils
!
! Revision 1.6  2005/10/05 14:10:01  haselbac
! Added routines which used to be in RFLU_ModStencils
!
! Revision 1.5  2005/04/15 15:07:05  haselbac
! Added routines to build and destroy c2c stencil in CSR format
!
! Revision 1.4  2004/12/29 21:09:44  haselbac
! Cosmetics only
!
! Revision 1.3  2004/01/22 16:03:59  haselbac
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC 
! and titan
!
! Revision 1.2  2003/12/05 16:53:11  haselbac
! Added cell-to-cell CSR routines for parallel higher-order scheme
!
! Revision 1.1  2003/12/04 03:29:22  haselbac
! Initial revision
!
! ******************************************************************************














