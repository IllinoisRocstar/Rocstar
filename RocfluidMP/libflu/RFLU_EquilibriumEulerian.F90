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
! Purpose: Add Equilibrium Eulerian corrections
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_EquilibriumEulerian.F90,v 1.5 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_EquilibriumEulerian(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: RFLU_FinishSD
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_EqEulCorr, &
                                  SPEC_RFLU_SetEEv
#endif

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

  INTEGER :: iSpec
  TYPE(t_global), POINTER :: global
#ifdef SPEC
  TYPE(t_spec_type), POINTER :: pSpecType
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EquilibriumEulerian',&
  'RFLU_EquilibriumEulerian.F90')

! ******************************************************************************
! Call various subroutines to finish computation, if necessary
! ******************************************************************************

#ifdef SPEC
  IF ( pRegion%mixtInput%indSd == 1 ) THEN
    CALL RFLU_FinishSD(pRegion)

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)    
    
      IF ( pSpecType%velocityMethod == SPEC_METHV_EQEUL ) THEN
        CALL SPEC_RFLU_SetEEv(pRegion,iSpec)
        CALL SPEC_EqEulCorr(pRegion,iSpec)
      END IF ! pRegion%specInput
    END DO ! iSpec
  END IF ! pRegion%mixtInput%indSd
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EquilibriumEulerian

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_EquilibriumEulerian.F90,v $
! Revision 1.5  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/11/27 01:48:31  haselbac
! Adapted to changes in EE implementation
!
! Revision 1.2  2004/08/02 13:57:28  haselbac
! Bug fix: Missing ifdefs, cosmetics
!
! Revision 1.1  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! ******************************************************************************







