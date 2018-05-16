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
! Purpose: Suite for turbulence postprocessing routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: TURB_ModPostFlu.F90,v 1.3 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE TURB_ModPostFlu

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI
 
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: TURB_RFLU_AllocMemPost, &
            TURB_RFLU_AllocMemPostVert, &
            TURB_RFLU_DeallocMemPost, &
            TURB_RFLU_DeallocMemPostVert
        
! ******************************************************************************
! Declarations and definitions
! ****************************************************************************** 
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: TURB_ModPostFlu.F90,v $ $Revision: 1.3 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


! ******************************************************************************
!
! Purpose: Allocate memory of turbulence post quantities.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE TURB_RFLU_AllocMemPost(pRegion)

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
  INTEGER :: nPostv, nTv, nVars, errorFlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'TURB_RFLU_AllocMemPost',&
  'TURB_ModPostFlu.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  nPostv =  1   ! for tv(TV_MIXT_MUET,:)
  nTv    =  pRegion%mixtInput%nTv
  nVars  =  pRegion%turbInput%nOutField
  pGrid  => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Turbulence variables
! ==============================================================================

  IF (nVars > 1) THEN
    ALLOCATE(pRegion%turb%vort(1,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%turb%vort')
    END IF ! global%error
  ENDIF

  IF (nVars > 2) THEN
    ALLOCATE(pRegion%turb%lens(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%turb%lens')
    END IF ! global%error

    nPostv = nPostv + 1
  ENDIF

! ==============================================================================
! Post quantities
! ==============================================================================

  pRegion%turbInput%nPostv = nPostv

  ALLOCATE(pRegion%turb%postv(nPostv,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%turb%postv')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE TURB_RFLU_AllocMemPost



! ******************************************************************************
!
! Purpose: Deallocate memory of turbulence post quantities.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE TURB_RFLU_DeallocMemPost(pRegion)

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
  INTEGER :: nVars, errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'TURB_RFLU_DeallocMemPost',&
  'TURB_ModPostFlu.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  nVars  =  pRegion%turbInput%nOutField

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Turbulence variables
! ==============================================================================

  IF (nVars > 1) THEN
    DEALLOCATE(pRegion%turb%vort,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%turb%vort')
    END IF ! global%error
  ENDIF

  IF (nVars > 2) THEN
    DEALLOCATE(pRegion%turb%lens,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%turb%lens')
    END IF ! global%error
  ENDIF

! ==============================================================================
! Post quantities
! ==============================================================================

  DEALLOCATE(pRegion%turb%postv,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%turb%postv')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE TURB_RFLU_DeallocMemPost



! ******************************************************************************
!
! Purpose: Allocate memory to store interpolated post quantities at vertices.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE TURB_RFLU_AllocMemPostVert(pRegion)

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
  INTEGER :: nPostv, errorFlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'TURB_RFLU_AllocMemPostVert',&
  'TURB_ModPostFlu.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  nPostv =  pRegion%turbInput%nPostv
  pGrid  => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Post quantities
! ==============================================================================

  ALLOCATE(pRegion%turb%postvVert(nPostv,pGrid%nVertTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%turb%postvVert')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE TURB_RFLU_AllocMemPostVert



! ******************************************************************************
!
! Purpose: Deallocate memory of post quantities at vertices.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE TURB_RFLU_DeallocMemPostVert(pRegion)

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
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'TURB_RFLU_DeallocMemPostVert',&
  'TURB_ModPostFlu.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Post quantities
! ==============================================================================

  DEALLOCATE(pRegion%turb%postvVert,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%turb%postvVert')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE TURB_RFLU_DeallocMemPostVert


! ******************************************************************************
! End
! ******************************************************************************
      
END MODULE TURB_ModPostFlu

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_ModPostFlu.F90,v $
! Revision 1.3  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/01/12 09:51:15  wasistho
! initial import
!
!
!
! ******************************************************************************










