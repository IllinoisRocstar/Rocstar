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
! Purpose: Destroy grid.
!
! Description: None.
!
! Input:
!   pRegion             Region pointer
!
! Output: None.
!
! Notes: 
!   1. Patch array deallocated at end of routine.
!
! ******************************************************************************
!
! $Id: RFLU_DestroyGrid.F90,v 1.10 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_DestroyGrid(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI
  
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

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_DestroyGrid.F90,v $ $Revision: 1.10 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DestroyGrid',&
  'RFLU_DestroyGrid.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying grid...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Coordinates and vertex flags
! ******************************************************************************

  DEALLOCATE(pGrid%xyz,STAT=errorFlag)
  global%error = errorFlag         
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%xyz')
  END IF ! global%error
  
  NULLIFY(pGrid%xyz)

! ******************************************************************************
! Connectivity
! ******************************************************************************

! ==============================================================================
! Tetrahedra
! ==============================================================================

  IF ( pGrid%nTetsTot > 0 ) THEN 
    DEALLOCATE(pGrid%tet2v,STAT=errorFlag) 
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%tet2v')
    END IF ! global%error
      
    NULLIFY(pGrid%tet2v)   
  END IF ! pGrid%nTetsTot

! ==============================================================================
! Hexahedra
! ==============================================================================

  IF ( pGrid%nHexsTot > 0 ) THEN 
    DEALLOCATE(pGrid%hex2v,STAT=errorFlag) 
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%hex2v')
    END IF ! global%error 
     
    NULLIFY(pGrid%hex2v)      
  END IF ! pGrid%nHexsTot

! ==============================================================================
! Prisms
! ==============================================================================

  IF ( pGrid%nPrisTot > 0 ) THEN 
    DEALLOCATE(pGrid%pri2v,STAT=errorFlag)
    global%error = errorFlag          
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%pri2v')
    END IF ! global%error
     
    NULLIFY(pGrid%pri2v)      
  END IF ! pGrid%nPrisTot               
        
! ==============================================================================
! Pyramids 
! ==============================================================================

  IF ( pGrid%nPyrsTot > 0 ) THEN 
    DEALLOCATE(pGrid%pyr2v,STAT=errorFlag)
    global%error = errorFlag          
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%pyr2v')
    END IF ! global%error
    
    NULLIFY(pGrid%pyr2v)   
  END IF ! pGrid%nPyrsTot

! ******************************************************************************
! Patches
! ******************************************************************************

! ==============================================================================
! Loop over patches
! ==============================================================================

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%nBTrisTot > 0 ) THEN 
      DEALLOCATE(pPatch%bTri2v,STAT=errorFlag)
      global%error = errorFlag             
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%patches%bTri2v')
      END IF ! global%error
      NULLIFY(pPatch%bTri2v)       
    END IF ! pPatch%nBTrisTot

    IF ( pPatch%nBQuadsTot > 0 ) THEN     
      DEALLOCATE(pPatch%bQuad2v,STAT=errorFlag)
      global%error = errorFlag             
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%patch%bQuad2v')
      END IF ! global%error
      NULLIFY(pPatch%bQuad2v)           
    END IF ! pPatch%nBQuadsTot 
    
    IF ( pPatch%nBCellsVirt > 0 ) THEN                                                                          
      DEALLOCATE(pPatch%bvc,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvc')
      END IF ! global%error
      NULLIFY(pPatch%bvc)        
    END IF ! pPatch%bcType                  
  END DO ! iPatch

! ******************************************************************************
! Deallocate actual patch array
! ******************************************************************************

  IF ( pGrid%nPatches > 0 ) THEN
    DEALLOCATE(pRegion%patches,STAT=errorFlag)  
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%patches')
    END IF ! global%error
    NULLIFY(pRegion%patches)   
  END IF ! pGrid%nPatches

! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying grid done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DestroyGrid

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DestroyGrid.F90,v $
! Revision 1.10  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.7  2006/03/25 21:42:18  haselbac
! Changes made bcos of sype boundaries
!
! Revision 1.6  2004/11/03 16:58:55  haselbac
! Removed deallocation of vertex and cell flags
!
! Revision 1.5  2004/10/19 19:24:15  haselbac
! Made consistent with changes in RFLU_CreateGrid.F90, cosmetics
!
! Revision 1.4  2003/12/10 03:58:02  haselbac
! Major bug fix: Fixed dealloc of patches if nPatches = 0
!
! Revision 1.3  2003/12/04 03:23:51  haselbac
! Clean-up
!
! Revision 1.2  2003/03/15 16:53:34  haselbac
! Clean up, bug fixes
!
! Revision 1.1  2003/01/28 15:53:32  haselbac
! Initial revision
!
! ******************************************************************************







