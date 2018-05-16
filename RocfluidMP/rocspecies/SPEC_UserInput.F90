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
! Purpose: Read user input, store it in the data structure and check.
!
! Description: None.
!
! Input:
!   regions        Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_UserInput.F90,v 1.5 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_UserInput(regions)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

  USE SPEC_ModInterfaces, ONLY: SPEC_CheckUserInput, &
                                SPEC_DerivedInputValues, &
                                SPEC_InitInputValues, &
                                SPEC_ReadInputFile

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_UserInput.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_UserInput',&
  'SPEC_UserInput.F90')

! ******************************************************************************
! Initialize, read, set, and check user input for species
! ******************************************************************************

  CALL SPEC_InitInputValues(regions)
  CALL SPEC_ReadInputFile(regions) ! global%specUsed is set here

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_DerivedInputValues(regions)
    CALL SPEC_CheckUserInput(regions)
  END IF ! global%specUsed

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_UserInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_UserInput.F90,v $
! Revision 1.5  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/11/10 02:39:24  haselbac
! Clean-up
!
! Revision 1.2  2004/07/28 15:29:21  jferry
! created global variable for spec use
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







