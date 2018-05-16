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
! Purpose: Read RADI user input and store it in RADI data structure.
!
! Description: RADI simulation parameters are first set to default values.
!              They are then overruled by user input parameters. Based
!              on this, derived RADI variables are defined. Finally,
!              necessary checkings are performed against inconsistencies.
!
! Input: regions = user input in all regions
!
! Output: regions = RADI user input, parameters and derived variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_UserInput.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE RADI_UserInput( regions )
#endif
#ifdef RFLU
SUBROUTINE RADI_UserInput
#endif

  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE RADI_ModInterfaces, ONLY : RADI_InitInputValues, RADI_ReadInputFile, &
                                 RADI_DerivedInputValues, RADI_CheckParamInput
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'RADI_UserInput',&
  'RADI_UserInput.F90' )

! initialize parameters --------------------------------------------------

#ifdef RFLO
  CALL RADI_InitInputValues( regions )
#endif
#ifdef RFLU
  CALL RADI_InitInputValues
#endif

! read user input

#ifdef RFLO
  CALL RADI_ReadInputFile( regions )
#endif
#ifdef RFLU
  CALL RADI_ReadInputFile
#endif

! perform global test whether radiation is active in any region

  DO iReg=1,global%nRegions
    IF (regions(iReg)%radiInput%radiModel /= RADI_MODEL_NONE) THEN
      global%radiActive = .TRUE.
    ENDIF
  ENDDO

  IF (global%radiActive) THEN

! - set RADI derived parameters or variables based on user input

#ifdef RFLO
    CALL RADI_DerivedInputValues( regions )
#endif
#ifdef RFLU
    CALL RADI_DerivedInputValues
#endif

! - check user input and parameters set in the code

#ifdef RFLO
    CALL RADI_CheckParamInput( regions )
#endif
#ifdef RFLU
    CALL RADI_CheckParamInput
#endif

  ENDIF  ! radiActive

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_UserInput.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/22 02:35:50  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.2  2003/08/01 22:16:25  wasistho
! prepared rocrad for Genx
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************







