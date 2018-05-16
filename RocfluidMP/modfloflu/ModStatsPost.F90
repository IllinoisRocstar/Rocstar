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
! Purpose: Suite for statistics postprocessing routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModStatsPost.F90,v 1.3 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE ModStatsPost

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI
 
  IMPLICIT NONE

  PRIVATE

#ifdef RFLU
  PUBLIC :: STAT_RFLU_AllocMemPost, &
            STAT_RFLU_AllocMemPostVert, &
            STAT_RFLU_DeallocMemPost, &
            STAT_RFLU_DeallocMemPostVert
        
! ******************************************************************************
! Declarations and definitions
! ****************************************************************************** 
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: ModStatsPost.F90,v $ $Revision: 1.3 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


! ******************************************************************************
!
! Purpose: Allocate memory of statistics post quantities.
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

SUBROUTINE STAT_RFLU_AllocMemPost(pRegion)

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
  INTEGER :: nStat, errorFlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'STAT_RFLU_AllocMemPost',&
  'ModStatsPost.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid  => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Mixture statistics
! ==============================================================================

  nStat = global%mixtNStat
  IF (nStat > 0) THEN
    ALLOCATE(pRegion%mixt%tav(nStat,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%tav')
    END IF ! global%error
  ENDIF

#ifdef TURB
! ==============================================================================
! Turbulence statistics
! ==============================================================================

  nStat = global%turbNStat
  IF (nStat > 0) THEN
    ALLOCATE(pRegion%turb%tav(nStat,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%turb%tav')
    END IF ! global%error
  ENDIF
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE STAT_RFLU_AllocMemPost



! ******************************************************************************
!
! Purpose: Deallocate memory of statistics post quantities.
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

SUBROUTINE STAT_RFLU_DeallocMemPost(pRegion)

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
  INTEGER :: nStat, errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'STAT_RFLU_DeallocMemPost',&
  'ModStatsPost.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Mixture statistics
! ==============================================================================

  nStat = global%mixtNStat
  IF (nStat > 0) THEN
    DEALLOCATE(pRegion%mixt%tav,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tav')
    END IF ! global%error
  ENDIF

#ifdef TURB
! ==============================================================================
! Turbulence statistics
! ==============================================================================

  nStat = global%turbNStat
  IF (nStat > 0) THEN
    DEALLOCATE(pRegion%turb%tav,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%turb%tav')
    END IF ! global%error
  ENDIF
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE STAT_RFLU_DeallocMemPost



! ******************************************************************************
!
! Purpose: Allocate memory to store interpolated statistics arrays at vertices.
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

SUBROUTINE STAT_RFLU_AllocMemPostVert(pRegion)

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
  INTEGER :: nStat, errorFlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'STAT_RFLU_AllocMemPostVert',&
  'ModStatsPost.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid  => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Mixture statistics
! ==============================================================================

  nStat = global%mixtNStat
  IF (nStat > 0) THEN
    ALLOCATE(pRegion%mixt%tavVert(nStat,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%tavVert')
    END IF ! global%error
  ENDIF

#ifdef TURB
! ==============================================================================
! Turbulence statistics
! ==============================================================================

  nStat = global%turbNStat
  IF (nStat > 0) THEN
    ALLOCATE(pRegion%turb%tavVert(nStat,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%turb%tavVert')
    END IF ! global%error
  ENDIF
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE STAT_RFLU_AllocMemPostVert



! ******************************************************************************
!
! Purpose: Deallocate memory of statistics arrays at vertices.
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

SUBROUTINE STAT_RFLU_DeallocMemPostVert(pRegion)

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
  INTEGER :: nStat, errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'STAT_RFLU_DeallocMemPostVert',&
  'ModStatsPost.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Mixture statistics
! ==============================================================================

  nStat = global%mixtNStat
  IF (nStat > 0) THEN
    DEALLOCATE(pRegion%mixt%tavVert,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tavVert')
    END IF ! global%error
  ENDIF

#ifdef TURB
! ==============================================================================
! Turbulence statistics
! ==============================================================================

  nStat = global%turbNStat
  IF (nStat > 0) THEN
    DEALLOCATE(pRegion%turb%tavVert,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%turb%tavVert')
    END IF ! global%error
  ENDIF
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE STAT_RFLU_DeallocMemPostVert


! ******************************************************************************
! End
! ******************************************************************************
#endif
      
END MODULE ModStatsPost

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModStatsPost.F90,v $
! Revision 1.3  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/01/12 09:51:49  wasistho
! initial import
!
!
!
! ******************************************************************************










