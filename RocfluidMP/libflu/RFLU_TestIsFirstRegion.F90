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
! Purpose: Determine whether region is first one.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: 
!   RFLU_TestIsFirstRegion = .TRUE.     if region is first region
!   RFLU_TestIsFirstRegion = .FALSE.    if region is not first region
!
! Notes: 
!   1. First region has index 0 for serial runs, and 1 for parallel runs.
!
! ******************************************************************************
!
! $Id: RFLU_TestIsFirstRegion.F90,v 1.3 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_TestIsFirstRegion(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
   
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
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_TestIsFirstRegion.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

! ******************************************************************************
! Test whether is first region
! ******************************************************************************

  RFLU_TestIsFirstRegion = .FALSE. 

  IF ( global%nRegions > 1 ) THEN 
    IF ( pRegion%iRegionGlobal == 1 ) THEN 
      RFLU_TestIsFirstRegion = .TRUE.  
    END IF ! pRegion%iRegionGlobal
  ELSE 
    IF ( pRegion%iRegionGlobal == 0 ) THEN 
      RFLU_TestIsFirstRegion = .TRUE. 
    END IF ! pRegion%iRegionGlobal  
  END IF ! global%nRegions

! ******************************************************************************
! End
! ******************************************************************************

END FUNCTION RFLU_TestIsFirstRegion

!*******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_TestIsFirstRegion.F90,v $
! Revision 1.3  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/10/19 19:23:51  haselbac
! Initial revision
!
!*******************************************************************************






