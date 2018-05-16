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
! Purpose: Wrapper for writing information on flow solution.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_PrintFlowInfoWrapper.F90,v 1.4 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_PrintFlowInfoWrapper(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_PrintFlowInfo

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_PrintFlowInfo
#endif

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Local variables
! ==============================================================================

  TYPE(t_global), POINTER :: global
  CHARACTER(CHRLEN) :: RCSIdentString

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PrintFlowInfoWrapper.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PrintFlowInfoWrapper',&
  'RFLU_PrintFlowInfoWrapper.F90')

! ******************************************************************************
! Print flow info
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

  CALL RFLU_PrintFlowInfo(pRegion)

! ==============================================================================
! Physical modules
! ==============================================================================

#ifdef SPEC
  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_PrintFlowInfo(pRegion)
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PrintFlowInfoWrapper


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_PrintFlowInfoWrapper.F90,v $
!   Revision 1.4  2008/12/06 08:44:12  mtcampbe
!   Updated license.
!
!   Revision 1.3  2008/11/19 22:17:25  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.2  2004/07/28 15:29:19  jferry
!   created global variable for spec use
!
!   Revision 1.1  2003/11/25 21:02:58  haselbac
!   Initial revision
!
! ******************************************************************************







