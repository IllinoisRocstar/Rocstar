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
! $Id: RFLU_CheckDerivedUserInput.F90,v 1.6 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE RFLU_CheckDerivedUserInput(regions)

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
    '$RCSfile: RFLU_CheckDerivedUserInput.F90,v $ $Revision: 1.6 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_CheckDerivedUserInput',&
  'RFLU_CheckDerivedUserInput.F90')

! *****************************************************************************
! Check region related data
! *****************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    pMixtInput => regions(iReg)%mixtInput

#ifdef PLAG
! =============================================================================
!   Check for consistency between RK scheme and PLAG module
! =============================================================================

    IF ( (global%plagUsed .EQV. .TRUE.) .AND. &
         (global%rkScheme /= RK_SCHEME_3_WRAY) ) THEN
      CALL ErrorStop(global,ERR_RK_SCHEME_INVALID,__LINE__, &
                     'PLAG requires RK3.')
    END IF ! global

! =============================================================================
!   Check for consistency between PLAG module and flag for converting 
!   Lagrangian to Eulerian field
! =============================================================================
 
    IF ( global%plagUsed .EQV. .FALSE. ) THEN
      IF ( global%postLag2EulFlag .EQV. .TRUE. ) THEN 
        IF ( iReg == LBOUND(regions,1) ) THEN
          global%warnCounter = global%warnCounter + 1

          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'*** WARNING *** Invalid '// &
                                   'input for Eulerian postprocessing flag.' 
          WRITE(STDOUT,'(A,20X,A)') SOLVER_NAME,'Setting flag to false.'
        END IF ! iReg

        global%postLag2EulFlag = .FALSE.
      END IF ! global%postLag2EulFlag
    END IF ! global%plagUsed
#endif

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

END SUBROUTINE RFLU_CheckDerivedUserInput

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckDerivedUserInput.F90,v $
! Revision 1.6  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/12/10 16:54:03  haselbac
! Added check for postLag2EulFlag
!
! Revision 1.3  2004/12/04 02:33:12  haselbac
! Bug fix: Missing brackets in logical expression
!
! Revision 1.2  2004/11/30 20:07:48  fnajjar
! Added error trap for running PLAG with none RK3 scheme
!
! Revision 1.1  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! *****************************************************************************







