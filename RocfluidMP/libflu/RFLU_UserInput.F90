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
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_UserInput.F90,v 1.6 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_UserInput(regions)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
       
  USE ModInterfaces, ONLY: RFLU_CheckUserInput, &
                           RFLU_DerivedInputValues, & 
                           RFLU_InitInputValues, &
                           ReadInputFile     
       
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_UserInput.F90,v $ $Revision: 1.6 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_UserInput',&
  'RFLU_UserInput.F90')

! ******************************************************************************
! Initialize, read input values, set derived values, and check 
! ******************************************************************************

  CALL RFLU_InitInputValues(regions)
  CALL ReadInputFile(regions)
  CALL RFLU_DerivedInputValues(regions)
  CALL RFLU_CheckUserInput(regions)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_UserInput.F90,v $
! Revision 1.6  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/11/10 02:17:44  haselbac
! Clean-up
!
! Revision 1.3  2004/10/19 19:37:56  haselbac
! Cosmetics only
!
! Revision 1.2  2003/11/25 21:02:53  haselbac
! Cosmetic changes only
!
! Revision 1.1  2003/01/28 15:53:32  haselbac
! Initial revision, moved from rocflu
!
! Revision 1.4  2002/09/09 15:51:56  haselbac
! global under regions, adapated interfaces
!
! Revision 1.3  2002/08/18 02:31:32  wasistho
! Added RFLU_CheckUserInput
!
! Revision 1.2  2002/05/04 17:14:02  haselbac
! Added call to RFLU_DerivedInputValues
!
! Revision 1.1  2002/03/26 19:24:49  haselbac
! Initial revision
!
!******************************************************************************







