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
! *****************************************************************************
!
! Purpose: Check parameters specified by the user:
!          These checks require knowledge of input to all MP modules
!
! Description: None.
!
! Input:
!   regions        Region data
!
! Output: None.
!
! Notes:
!
! *****************************************************************************
!
! $Id: RFLO_CheckDerivedUserInput.F90,v 1.5 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE RFLO_CheckDerivedUserInput(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY: t_mixt_input
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! *****************************************************************************
! Arguments
! *****************************************************************************

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! *****************************************************************************
! Locals
! *****************************************************************************

  INTEGER :: iReg
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLO_CheckDerivedUserInput.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLO_CheckDerivedUserInput',&
  'RFLO_CheckDerivedUserInput.F90')

! *****************************************************************************
! Check region related data
! *****************************************************************************

  DO iReg=1,global%nRegions
    pMixtInput => regions(iReg)%mixtInput

! TEMPORARY: Deactivated while checking RK3 Stability
#ifdef PLAG
! =============================================================================
!   Check for consistency between RK scheme and PLAG module
! =============================================================================

!    IF ( (global%plagUsed .EQV. .TRUE.) .AND. &
!         (global%rkScheme /= RK_SCHEME_3_WRAY) ) THEN
!      CALL ErrorStop(global,ERR_RK_SCHEME_INVALID,__LINE__, &
!                     'PLAG requires RK3.')
!    END IF ! global
#endif
! END TEMPORARY: Deactivated while checking RK3 Stability

#ifndef TURB
! =============================================================================
!   Check for consistency between flow model and state of TURB module
! =============================================================================

    IF ( pMixtInput%flowModel == FLOW_NAVST .AND. &
         pMixtInput%turbModel /= TURB_MODEL_NONE ) THEN
      CALL ErrorStop(global,ERR_TURB_MODULE,__LINE__, &
                     'Turbulence is on, compile TURB')
    END IF ! pMixtInput
#endif    

! =============================================================================
!   Check for valid input for viscosity model
! =============================================================================

    IF (  pMixtInput%computeTv .AND. &
        ( pMixtInput%viscModel <  VISC_SUTHR .OR.  &
          pMixtInput%viscModel >  VISC_ANTIB )     ) THEN
      CALL ErrorStop(global,ERR_UNKNOWN_VISCMODEL,__LINE__)
    END IF ! pMixtInput
  END DO ! iReg

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_CheckDerivedUserInput

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CheckDerivedUserInput.F90,v $
! Revision 1.5  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/12/06 21:07:57  fnajjar
! Deactivated consistency check for PLAG and RK scheme while testing stability
!
! Revision 1.2  2004/12/04 02:39:10  haselbac
! Bug fix: Missing brackets in logical expression
!
! Revision 1.1  2004/11/30 20:12:48  fnajjar
! Initial import
!
! *****************************************************************************







