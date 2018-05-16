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
! Purpose: Write file with special cells.
!
! Description: None.
!
! Input: 
!   pRegion              Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_WritePostInfo.F90,v 1.7 2008/12/06 08:45:04 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_WritePostInfo(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI
  
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
  INTEGER :: ics,iFile,ifs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WritePostInfo.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_WritePostInfo',&
  'RFLU_WritePostInfo.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing post-processor info...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal
  END IF ! global%myProcid

! ******************************************************************************
! Set grid pointer and variables
! ******************************************************************************

  pGrid => pRegion%grid

  iFile = IF_POSTINFO

! ******************************************************************************
! Write to file
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Number of special cells:', &
                                   pGrid%nCellsSpecial
    WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Number of special faces:', &
                                   pGrid%nFacesSpecial                                   
  END IF ! global%myProcid

  WRITE(iFile,'(A,1X,I5.5)') '# Region',pRegion%iRegionGlobal
  WRITE(iFile,'(L1)') pRegion%postActiveFlag
  
  IF ( pRegion%postActiveFlag .EQV. .TRUE. ) THEN   
    WRITE(iFile,*) pGrid%nCellsSpecial

    DO ics = 1,pGrid%nCellsSpecial
      WRITE(iFile,*) pGrid%cellsSpecial(ics)
    END DO ! ics
    
    WRITE(iFile,*) pGrid%nFacesSpecial
    
    DO ifs = 1,pGrid%nFacesSpecial
      WRITE(iFile,*) pGrid%facesSpecial(1:2,ifs)
    END DO ! ifs    
  END IF ! pRegion%postActiveFlag

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing post-processor info done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WritePostInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WritePostInfo.F90,v $
! Revision 1.7  2008/12/06 08:45:04  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/09/27 01:40:14  haselbac
! Added writing of special faces
!
! Revision 1.4  2004/03/23 03:18:46  haselbac
! Changed format statements
!
! Revision 1.3  2004/03/15 21:10:09  haselbac
! Fixed alignment
!
! Revision 1.2  2003/08/07 15:35:50  haselbac
! Changed var name, fixed small bug
!
! Revision 1.1.1.1  2003/06/04 22:31:20  haselbac
! Initial revision
!
! ******************************************************************************







