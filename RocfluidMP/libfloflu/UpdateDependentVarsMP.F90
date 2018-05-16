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
! Purpose: Update dependent variables.
!
! Description: None.
!
! Input:
!   region        Region data
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: UpdateDependentVarsMP.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE UpdateDependentVarsMP(region)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
#ifdef RFLU
  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons
#endif

#ifdef RFLO
  USE ModInterfaces, ONLY: RFLO_GetDimensDummy,RFLO_GetCellOffset

#include "Indexing.h"
#endif

  USE ModInterfaces, ONLY: MixtureProperties

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_nonCvUpdate
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_UpdateDependentVars
#endif

  IMPLICIT NONE

! *****************************************************************************
! Declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ibc,iec
#ifdef RFLO
  INTEGER :: iCOff,idcbeg,idcend,ijCOff,iLev,jdcbeg,jdcend,kdcbeg,kdcend
#endif
  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: UpdateDependentVarsMP.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction(global,'UpdateDependentVarsMP',&
  'UpdateDependentVarsMP.F90')

! *****************************************************************************
! Set variables
! *****************************************************************************

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy(region,iLev,idcbeg,idcend,jdcbeg,jdcend,kdcbeg, &
                           kdcend)
  CALL RFLO_GetCellOffset(region,iLev,iCOff,ijCOff)
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot
#endif

! *****************************************************************************
! Update dependent variables
! *****************************************************************************

! =============================================================================
! Mixture. NOTE here the last parameter MUST be TRUE, otherwise the gas
! properties do not get initialized correctly.
! =============================================================================

  CALL MixtureProperties(region,ibc,iec,.TRUE.)

#ifdef RFLU
#ifdef PLAG
! =============================================================================
! Particles
! =============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    pRegion => region%pRegion
    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)
    CALL PLAG_nonCvUpdate(pRegion)
    CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)
  END IF ! plagUsed
#endif
#endif

#ifdef SPEC
! =============================================================================
! Species
! =============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    pRegion => region%pRegion
    CALL SPEC_UpdateDependentVars(pRegion,ibc,iec)
  END IF ! global%specUsed
#endif

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE UpdateDependentVarsMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateDependentVarsMP.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/04/15 15:06:04  haselbac
! Adapted interface to SPEC_UpdateDependentVars
!
! Revision 1.1  2004/12/01 16:51:59  haselbac
! Initial revision after changing case
!
! Revision 1.6  2004/11/29 17:15:20  wasistho
! use ModInterfacesSpecies
!
! Revision 1.5  2004/11/14 19:37:18  haselbac
! Changed interfaces to PLAG_nonCvUpdate and SPEC_UpdateDependentVars
!
! Revision 1.4  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.3  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/02/26 21:01:50  haselbac
! Added PLAG support
!
! Revision 1.1  2004/01/29 22:52:33  haselbac
! Initial revision
!
!******************************************************************************







