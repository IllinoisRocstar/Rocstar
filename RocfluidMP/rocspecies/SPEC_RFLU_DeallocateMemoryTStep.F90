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
!******************************************************************************
!
! Purpose: Deallocate memory for species related to time stepping.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: SPEC_RFLU_DeallocateMemoryTStep.F90,v 1.6 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_RFLU_DeallocateMemoryTStep(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid  
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input    
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace
   
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag, iPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput  
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_DeallocateMemoryTStep.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_DeallocateMemoryTStep',&
  'SPEC_RFLU_DeallocateMemoryTStep.F90')

! ******************************************************************************
! Set pointers 
! ******************************************************************************

  pGrid      => pRegion%grid
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================     
! Old state vector
! ==============================================================================

  DEALLOCATE(pRegion%spec%cvOld,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%spec%cvOld')
  END IF ! global%error  

! ==============================================================================     
! Residuals 
! ==============================================================================
  
  DEALLOCATE(pRegion%spec%rhs,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%spec%rhs')
  END IF ! global%error  
  
  DEALLOCATE(pRegion%spec%diss,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%spec%diss')
  END IF ! global%error  

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    DEALLOCATE(pRegion%spec%rhsSum,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%spec%rhsSum')
    END IF ! global%error              
  END IF ! global%flowType   

! ==============================================================================
! Gradients 
! ==============================================================================
  
! ------------------------------------------------------------------------------
! Cell gradients
! ------------------------------------------------------------------------------  
  
  IF ( pMixtInput%spaceOrder > 1 ) THEN 
    DEALLOCATE(pRegion%spec%gradCell,STAT=errorFlag)
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%spec%gradCell')
    END IF ! global%error     
  END IF ! pMixtInput%spaceOrder

! ------------------------------------------------------------------------------
! Face gradients
! ------------------------------------------------------------------------------  
  
  IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN 
    DEALLOCATE(pRegion%spec%gradFace,STAT=errorFlag)
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%spec%gradFace')
    END IF ! global%error
  END IF ! pMixtInput%flowModel
  
  IF ( pGrid%nPatches > 0 ) THEN
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        DEALLOCATE(pPatch%spec%gradFace,STAT=errorFlag)
        global%error = errorFlag         
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%spec%gradFace')
        END IF ! global%error                              
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
  END IF ! pGrid%nPatches
  
    
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_DeallocateMemoryTStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_DeallocateMemoryTStep.F90,v $
! Revision 1.6  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:40:28  mparmar
! Moved region%mixt%bGradFace to patch%mixt%gradFace
!
! Revision 1.3  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.2  2004/01/29 22:59:41  haselbac
! Added deallocation of cell and face gradients
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
!******************************************************************************







