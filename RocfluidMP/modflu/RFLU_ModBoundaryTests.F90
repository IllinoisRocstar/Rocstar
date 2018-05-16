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
! Purpose: Collection of routines for testing whether vertices or cells are on
!   boundaries.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModBoundaryTests.F90,v 1.3 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModBoundaryTests

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
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
    RCSIdentString = '$RCSfile: RFLU_ModBoundaryTests.F90,v $ $Revision: 1.3 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: RFLU_TestIsBoundaryCell, & 
            RFLU_TestIsBoundaryVertex

! ==============================================================================
! Private functions
! ==============================================================================



! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Determine whether cell is adjacent to boundaries.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   icg                 Global cell index
!
! Output: 
!   RFLU_TestIsBoundaryCell = .TRUE.     If cell adjacent to boundaries
!   RFLU_TestIsBoundaryCell = .FALSE.    If cell not adjacent to boundaries
!
! Notes: 
!   1. This routine is not efficient if calling it for large numbers of cells.
!      Would need to consider sorting actual boundary-face lists which would
!      avoid the need to sort them (repeatedly) here.
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_TestIsBoundaryCell(pRegion,icg)
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icg
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifl,iLoc,iPatch
  INTEGER, DIMENSION(:), ALLOCATABLE :: bf2cSorted
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TestIsBoundaryCell',&
  'RFLU_ModBoundaryTests.F90')

! ==============================================================================
! Set grid pointer and initialize RFLU_TestIsBoundaryCell
! ==============================================================================

  pGrid => pRegion%grid
  
  RFLU_TestIsBoundaryCell = .FALSE.

! ******************************************************************************
! Check whether point is on boundaries
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
  
    ALLOCATE(bf2cSorted(pPatch%nBFacesTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bf2cSorted')
    END IF ! global%error
    
    DO ifl = 1,pPatch%nBFacesTot ! Explicit loop because of ASCI White
      bf2cSorted(ifl) = pPatch%bf2c(ifl)
    END DO ! ifl
    
    CALL QuickSortInteger(bf2cSorted,pPatch%nBFacesTot)
    
    CALL BinarySearchInteger(bf2cSorted,pPatch%nBFacesTot,icg,iLoc)
    
    IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
      RFLU_TestIsBoundaryCell = .TRUE.
    END IF ! iLoc
    
    DEALLOCATE(bf2cSorted,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bf2cSorted')
    END IF ! global%error   
    
    IF ( RFLU_TestIsBoundaryCell .EQV. .TRUE. ) THEN 
      EXIT
    END IF ! RFLU_TestIsBoundaryCell         
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_TestIsBoundaryCell


   
   
   
! ******************************************************************************
!
! Purpose: Determine whether vertex is on boundaries.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   ivg                 Global vertex index
!
! Output: 
!   RFLU_TestIsBoundaryVertex = .TRUE.     If vertex on boundaries
!   RFLU_TestIsBoundaryVertex = .FALSE.    If vertex not on boundaries
!
! Notes: None.
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_TestIsBoundaryVertex(pRegion,ivg)
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: ivg
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,iLoc,iPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TestIsBoundaryVertex',&
  'RFLU_ModBoundaryTests.F90')

! ==============================================================================
! Set grid pointer and initialize RFLU_TestIsBoundaryVertex
! ==============================================================================

  pGrid => pRegion%grid
  
  RFLU_TestIsBoundaryVertex = .FALSE.

! ******************************************************************************
! Check whether vertex is on boundary
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    CALL BinarySearchInteger(pPatch%bv,pPatch%nBVertTot,ivg,iLoc)
    
    IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
      RFLU_TestIsBoundaryVertex = .TRUE.

      EXIT
    END IF ! RFLU_TestIsBoundaryVertex        
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_TestIsBoundaryVertex
   
   
   
  


END MODULE RFLU_ModBoundaryTests

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModBoundaryTests.F90,v $
! Revision 1.3  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:31  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/08/04 02:59:30  haselbac
! Initial revision
!
! ******************************************************************************










