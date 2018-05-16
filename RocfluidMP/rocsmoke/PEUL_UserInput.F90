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
! Purpose: read user input and store it in the data structure.
!          Check user input.
!
! Description: none.
!
! Input: regions = dimensions and topology (finest grid).
!
! Output: regions = dimensions, topology and user input on all grid levels.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_UserInput.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_UserInput( regions ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters

  USE PEUL_ModInterfaces, ONLY : PEUL_InitInputValues, PEUL_ReadInputFile, &
                                 PEUL_ReadBcInputFile, PEUL_DerivedInputValues
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_patch),  POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_UserInput.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_UserInput',&
  'PEUL_UserInput.F90' )

! begin -----------------------------------------------------------------------

! initialize parameters

  CALL PEUL_InitInputValues( regions )

! read user input (global%peulUsed is set here)

  CALL PEUL_ReadInputFile( regions )

  IF (global%peulUsed) THEN

! - read boundary conditions

    CALL PEUL_ReadBcInputFile( regions )

! - set model & numerical parameters from user input

    CALL PEUL_DerivedInputValues( regions )

  END IF ! peulUsed

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_UserInput.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:10:03  haselbac
! Initial revision after changing case
!
! Revision 1.4  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.3  2003/03/11 16:04:58  jferry
! Created data type for material properties
!
! Revision 1.2  2003/02/12 19:03:42  jferry
! removed unused loop variables
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







