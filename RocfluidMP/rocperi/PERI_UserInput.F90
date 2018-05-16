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
! Purpose: Read PERI user input and store it in PERI data structure.
!
! Description: PERI simulation parameters are first set to default values.
!              They are then overruled by user input parameters. Based
!              on this, derived PERI variables are defined. Finally,
!              necessary checkings are performed against inconsistencies.
!
! Input: regions = user input in all regions
!
! Output: regions = PERI user input, parameters and derived variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_UserInput.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_UserInput( regions ) ! PUBLIC

  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE PERI_ModInterfaces, ONLY : PERI_CheckParamInput, PERI_DerivedInputValues, &
                                 PERI_InitInputValues, PERI_ReadInputFile
  USE ModError
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_UserInput.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'PERI_UserInput',&
  'PERI_UserInput.F90' )

! initialize parameters --------------------------------------------------

  CALL PERI_InitInputValues( regions )

! read user input

  CALL PERI_ReadInputFile( regions )

  IF (regions(1)%periInput%flowKind /= PERI_FLOW_NONE) THEN

! - set PERI parameters from user input

    CALL PERI_DerivedInputValues( regions )

! - check user input and parameters set in the code

    CALL PERI_CheckParamInput( regions )

  ENDIF  ! periInput%flowModel

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_UserInput.F90,v $
! Revision 1.4  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************







