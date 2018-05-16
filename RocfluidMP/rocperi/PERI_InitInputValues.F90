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
! Purpose: Initialize user input parameters for PERI to default values.
!
! Description: User input parameters are set to default before overruled by 
!              user input
!
! Input: regions data
!
! Output: regions = initial/default values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_InitInputValues.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_InitInputValues( regions )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModPeriodic, ONLY   : t_peri_input
  USE ModError
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN)           :: RCSIdentString
  TYPE(t_global), POINTER     :: global
  TYPE(t_peri_input), POINTER :: input

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_InitInputValues.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'PERI_InitInputValues',&
  'PERI_InitInputValues.F90' )

! region related values -------------------------------------------------------

#ifdef RFLO
  DO iReg=1,global%nRegions
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
#endif

    input => regions(iReg)%periInput

    input%flowKind   = PERI_FLOW_NONE
    input%nVar       = 0
    input%split(:)   = 0
    input%pgradType  = 0
    input%meanPgrad  = 0._RFREAL
    input%minjRate   = 0._RFREAL
    input%bulkmFlux  = 0._RFREAL
    input%cprEpsilon = 0._RFREAL
    input%headPres   = 0._RFREAL
    input%headTemp   = 0._RFREAL

  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_InitInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_InitInputValues.F90,v $
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
! Revision 1.3  2003/09/18 01:56:54  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.2  2003/07/17 01:20:08  wasistho
! RFLU, USE ModTurbulence to USE ModPeriodic
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************







