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
! Purpose: Suite of functions to determine information on cell and face 
!   types and kinds, before and after renumbering.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModCellFaceEdgeInfo.F90,v 1.12 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModCellFaceEdgeInfo

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid

  USE RFLU_ModCellMapping

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_GetGlobalCellType, & 
            RFLU_GetGlobalCellKind, &  
            RFLU_GetVirtualCellReg, &  
            RFLU_GetFaceKind, & 
            RFLU_GetEdgeKind
  
  
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModCellFaceEdgeInfo.F90,v $ $Revision: 1.12 $' 
                            
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Determine global cell type (after renumbering).
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pGrid       Pointer to grid data
!   icg         Global cell id
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  INTEGER FUNCTION RFLU_GetGlobalCellType(global,pGrid,icg)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: icg
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************  
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_GetGlobalCellType',&
  'RFLU_ModCellFaceEdgeInfo.F90') 

! ******************************************************************************  
!   Determine global cell type
! ******************************************************************************  

    IF ( icg > 0 ) THEN ! Interior cell
      RFLU_GetGlobalCellType = pGrid%cellGlob2Loc(1,icg) 
    ELSE IF ( icg == CELL_TYPE_BND ) THEN ! Boundary cell
      RFLU_GetGlobalCellType = CELL_TYPE_BND
    ELSE IF ( icg == CELL_TYPE_EXT ) THEN ! Exterior cell
      RFLU_GetGlobalCellType = CELL_TYPE_EXT 
    ELSE   
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! icg 

! ******************************************************************************  
!   End
! ******************************************************************************  

    CALL DeregisterFunction(global)

  END FUNCTION RFLU_GetGlobalCellType








! ******************************************************************************
!
! Purpose: Determine global cell kind (after renumbering)
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pGrid       Pointer to grid data
!   icg         Global cell id
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  INTEGER FUNCTION RFLU_GetGlobalCellKind(global,pGrid,icg)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: icg
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString
    INTEGER :: icl,ict    

! ******************************************************************************  
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_GetGlobalCellKind',&
  'RFLU_ModCellFaceEdgeInfo.F90')

! ******************************************************************************  
!   Get global cell kind
! ******************************************************************************  

    IF ( icg > 0 ) THEN 
      ict = pGrid%cellGlob2Loc(1,icg)
      icl = pGrid%cellGlob2Loc(2,icg)  
    ELSE IF ( icg == CELL_TYPE_BND ) THEN
      ict = CELL_TYPE_BND
    ELSE IF ( icg == CELL_TYPE_EXT ) THEN 
      ict = CELL_TYPE_EXT
    ELSE   
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! icg

    SELECT CASE (ict) 
      CASE ( CELL_TYPE_TET )
        IF ( icl <= pGrid%nTets ) THEN 
          RFLU_GetGlobalCellKind = CELL_KIND_ACTUAL
        ELSE 
          RFLU_GetGlobalCellKind = CELL_KIND_VIRTUAL
        END IF ! icl
      CASE ( CELL_TYPE_HEX )  
        IF ( icl <= pGrid%nHexs ) THEN 
          RFLU_GetGlobalCellKind = CELL_KIND_ACTUAL
        ELSE 
          RFLU_GetGlobalCellKind = CELL_KIND_VIRTUAL
        END IF ! icl
      CASE ( CELL_TYPE_PRI ) 
        IF ( icl <= pGrid%nPris ) THEN 
          RFLU_GetGlobalCellKind = CELL_KIND_ACTUAL
        ELSE 
          RFLU_GetGlobalCellKind = CELL_KIND_VIRTUAL
        END IF ! icl
      CASE ( CELL_TYPE_PYR )         
        IF ( icl <= pGrid%nPyrs ) THEN 
          RFLU_GetGlobalCellKind = CELL_KIND_ACTUAL
        ELSE 
          RFLU_GetGlobalCellKind = CELL_KIND_VIRTUAL
        END IF ! icl
      CASE ( CELL_TYPE_BND )
        RFLU_GetGlobalCellKind = CELL_KIND_BND
      CASE ( CELL_TYPE_EXT ) 
        RFLU_GetGlobalCellKind = CELL_KIND_EXT
      CASE DEFAULT
        WRITE(errorString,'(1X,I3)') ict
        CALL ErrorStop(global,ERR_CELL_TYPE,__LINE__,TRIM(errorString))
    END SELECT ! icType

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END FUNCTION RFLU_GetGlobalCellKind







! ******************************************************************************
!
! Purpose: Determine virtual cell region index.
!
! Description: Loop over cells which are to be received, if cell among these, 
!   then know region from border index.
!
! Input:
!   global      Pointer to global data
!   pGrid       Pointer to grid data
!   icg         Global cell id
!
! Output: None.
!
! Notes: 
!   1. If the cell is not found in any border list, return ELEMENT_NOT_FOUND.
!
! ******************************************************************************
  
  INTEGER FUNCTION RFLU_GetVirtualCellReg(global,pGrid,icg)

    USE ModBorder, ONLY: t_border

    USE ModSortSearch

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: icg
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: iBorder,iLoc
    TYPE(t_border), POINTER :: pBorder

! ******************************************************************************  
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_GetVirtualCellReg',&
  'RFLU_ModCellFaceEdgeInfo.F90') 

! ******************************************************************************  
!   Determine global cell type
! ******************************************************************************  

    RFLU_GetVirtualCellReg = ELEMENT_NOT_FOUND  

    borderLoop: DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      CALL BinarySearchInteger(pBorder%icgRecv,pBorder%nCellsRecv,icg,iLoc)
      
      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN
        RFLU_GetVirtualCellReg = pBorder%iRegionGlobal

        EXIT borderLoop
      END IF ! iLoc
    END DO borderLoop

! ******************************************************************************  
!   End
! ******************************************************************************  

    CALL DeregisterFunction(global)

  END FUNCTION RFLU_GetVirtualCellReg







! ******************************************************************************
!
! Purpose: Determine face kind before renumbering.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   c1k         Kind of cell 1
!   c2k         Kind of cell 2 
!
! Output: None.
!
! Notes: 
!   1. Can only get FACE_KIND_AB if have partitioned case with no dummy cells, 
!      which must not arise in practice but may be of interest for testing 
!      purposes.
!   2. Can only get FACE_KIND_VX if have partitioned case.
!
! ******************************************************************************
  
  INTEGER FUNCTION RFLU_GetFaceKind(global,c1k,c2k)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: c1k,c2k
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString    

! ******************************************************************************  
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_GetFaceKind',&
  'RFLU_ModCellFaceEdgeInfo.F90')

! ******************************************************************************  
!   Get face kind
! ******************************************************************************  

    SELECT CASE (c1k + c2k) 
      CASE ( CELL_KIND_ACTUAL  + CELL_KIND_ACTUAL )
        RFLU_GetFaceKind = FACE_KIND_AA
      CASE ( CELL_KIND_ACTUAL  + CELL_KIND_VIRTUAL ) 
        RFLU_GetFaceKind = FACE_KIND_AV
      CASE ( CELL_KIND_VIRTUAL + CELL_KIND_VIRTUAL )
        RFLU_GetFaceKind = FACE_KIND_VV
      CASE ( CELL_KIND_VIRTUAL + CELL_KIND_BND )
        RFLU_GetFaceKind = FACE_KIND_VB 
      CASE ( CELL_KIND_VIRTUAL + CELL_KIND_EXT ) 
        RFLU_GetFaceKind = FACE_KIND_VX
      CASE ( CELL_KIND_ACTUAL  + CELL_KIND_BND )
        RFLU_GetFaceKind = FACE_KIND_AB           
      CASE DEFAULT
        WRITE(errorString,'(1X,I3)') c1k+c2k
        CALL ErrorStop(global,ERR_FACE_KIND,__LINE__,TRIM(errorString))
    END SELECT ! c1k + c2k

! ******************************************************************************  
!   End
! ******************************************************************************  

    CALL DeregisterFunction(global)    

  END FUNCTION RFLU_GetFaceKind
    
    
    
    
    
! ******************************************************************************
!
! Purpose: Determine edge kind.
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
    
  INTEGER FUNCTION RFLU_GetEdgeKind(global,pGrid,v1,v2)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: v1,v2
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: v1k,v2k    

! ******************************************************************************  
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_GetEdgeKind',&
  'RFLU_ModCellFaceEdgeInfo.F90')

! ******************************************************************************  
!   Get face kind
! ******************************************************************************  

    IF ( v1 <= pGrid%nVert ) THEN 
      v1k = VERT_KIND_ACTUAL
    ELSE 
      v1k = VERT_KIND_VIRTUAL
    END IF ! v1

    IF ( v2 <= pGrid%nVert ) THEN 
      v2k = VERT_KIND_ACTUAL
    ELSE 
      v2k = VERT_KIND_VIRTUAL
    END IF ! v2

    SELECT CASE (v1k + v2k) 
      CASE ( VERT_KIND_ACTUAL + VERT_KIND_ACTUAL )
        RFLU_GetEdgeKind = EDGE_KIND_AA
      CASE ( VERT_KIND_ACTUAL + VERT_KIND_VIRTUAL ) 
        RFLU_GetEdgeKind = EDGE_KIND_AV
      CASE ( VERT_KIND_VIRTUAL + VERT_KIND_VIRTUAL )
        RFLU_GetEdgeKind = EDGE_KIND_VV 
      CASE DEFAULT
        CALL ErrorStop(global,ERR_CELL_TYPE,__LINE__)
    END SELECT ! v1k + v2k

! ******************************************************************************  
!   End
! ******************************************************************************  

    CALL DeregisterFunction(global) 

  END FUNCTION RFLU_GetEdgeKind
    
    
    

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModCellFaceEdgeInfo


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModCellFaceEdgeInfo.F90,v $
! Revision 1.12  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:31  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.9  2005/04/15 15:06:45  haselbac
! Added RFLU_GetVirtualCellReg
!
! Revision 1.8  2004/11/03 17:01:20  haselbac
! Rewrite because of removal of vertex anc cell flags, cosmetics
!                                       
! Revision 1.7  2004/01/22 16:03:58  haselbac                                  
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC   and titan
!
! Revision 1.6  2003/11/20 21:30:32  haselbac                                  
! Added register and deregister calls for better error reporting               
!
! Revision 1.5  2003/11/20 16:40:36  mdbrandy                                  
! Backing out RocfluidMP changes from 11-17-03                                 
!
! Revision 1.2  2003/08/19 22:47:54  haselbac                                  
! Improved error reporting                                                     
!
! Revision 1.1  2003/03/15 18:16:29  haselbac                                  
! Formerly called RFLU_ModCellFaceInfo.F90                                     
!
! Revision 1.5  2003/01/28 16:25:14  haselbac                                  
! Made simpler and more self-consistent                                        
!
! Revision 1.4  2002/10/05 19:02:31  haselbac                                  
! Fixed comment                                                                
!
! Revision 1.3  2002/09/09 15:02:02  haselbac                                  
! global now under regions                                                     
!
! Revision 1.2  2002/07/25 14:58:19  haselbac                                  
! Fixed bug in RFLU_GetCellIndexAfter and RFLU_GetCellKindBefore               
!
! Revision 1.1  2002/06/27 15:48:16  haselbac                                  
! Initial revision                                                             
!
! ******************************************************************************











