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
! Purpose: Pick regions by coordinate range.
!
! Description: Check whether a given region is inside a user-specified 
!   bounding box; if yes, flag region as active, otherwise as inactive.
!
! Input: 
!   pRegion		Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. If a given region is not inside the bounding box, then postActiveFlag
!      is set to FALSE.
!
! ******************************************************************************
!
! $Id: RFLU_PickRegionsCoord.F90,v 1.5 2008/12/06 08:45:04 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PickRegionsCoord(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModMPI
  USE ModParameters
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: xMax,xMin,yMax,yMin,zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PickRegionsCoord.F90,v $ $Revision: 1.5 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PickRegionsCoord',&
  'RFLU_PickRegionsCoord.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking region by coordinates...'
  END IF ! global%verbLevel

! ******************************************************************************
! Pick by bounding box
! ******************************************************************************

  pGrid => pRegion%grid

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVertTot))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVertTot))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVertTot))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVertTot))

  IF ( xMin > global%pickXCoordLow .AND. xMax < global%pickXCoordUpp .AND. & 
       yMin > global%pickYCoordLow .AND. yMax < global%pickYCoordUpp .AND. & 
       zMin > global%pickZCoordLow .AND. zMax < global%pickZCoordUpp ) THEN 
    pRegion%postActiveFlag = .TRUE.

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Picked region:', &
                                     pRegion%iRegionGlobal
    END IF ! global%verbLevel
  ELSE 
    pRegion%postActiveFlag = .FALSE.  
    
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I6,1X,A)') SOLVER_NAME,'Region:', &
                                          pRegion%iRegionGlobal,'not picked.'
    END IF ! global%verbLevel    
  END IF ! xMin    
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. & 
       (global%verbLevel > VERBOSE_NONE) ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking region by coordinates done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PickRegionsCoord

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PickRegionsCoord.F90,v $
! Revision 1.5  2008/12/06 08:45:04  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/12/10 23:29:34  haselbac
! Renamed geom post variables, cosmetics
!
! Revision 1.2  2004/10/19 19:30:08  haselbac
! Added output statement, cosmetics
!
! Revision 1.1  2003/08/07 15:14:40  haselbac
! Initial revision
!
! ******************************************************************************







