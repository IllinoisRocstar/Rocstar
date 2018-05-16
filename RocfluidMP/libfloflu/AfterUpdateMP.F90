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
! Purpose: check and/or modify conserved variable fields for RocfluidMP
!
! Description: none.
!
! Input: pRegion = data of current region,
!
! Output: pRegion%levels%*%cv = modified solution
!
! Notes: none.
!
!******************************************************************************
!
! $Id: AfterUpdateMP.F90,v 1.6 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
!******************************************************************************

SUBROUTINE AfterUpdateMP( pRegion,istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
#ifdef INRT
  USE INRT_ModParameters
#endif

#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_CheckValidity
#endif

#ifdef RFLU
  USE ModInterfaces, ONLY : RFLU_CheckPositivityWrapper, &
                            RFLU_CheckValidityWrapper,   &
                            RFLU_EnforceBoundsWrapper
#endif

#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_SetParticleTemp,       &
                                    INRT_VaporEnergyConversion, &
                                    INRT_BurnStatusUpdate
#ifdef RFLO
  USE ModInterfacesInteract, ONLY : INRT_TwoDimAverage
#endif
#endif

#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_EnforcePositivity
#endif

#ifdef PLAG
  USE PLAG_ModCheckVars, ONLY: PLAG_CheckValidity
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER    :: pRegion
  INTEGER,        INTENT(IN) :: istage

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  LOGICAL :: finalStage

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: AfterUpdateMP.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction( global,'AfterUpdateMP',&
  'AfterUpdateMP.F90' )

  finalStage = (istage == global%nrkSteps)

! check positivity ------------------------------------------------------------

#ifdef RFLO
  IF (finalStage) &
    CALL RFLO_CheckValidity(pRegion)
#ifdef PLAG
  IF (finalStage) &
    CALL PLAG_CheckValidity(pRegion)
#endif
#endif

#ifdef RFLU
  CALL RFLU_CheckValidityWrapper(pRegion)
  CALL RFLU_EnforceBoundsWrapper(pRegion)
  CALL RFLU_CheckPositivityWrapper(pRegion)
#endif

#ifdef INRT
  IF (global%inrtUsed .AND. finalStage) THEN

    IF (pRegion%inrtInput%inrts(INRT_TYPE_BURNING)%used) THEN
      CALL INRT_SetParticleTemp(pRegion)
      CALL INRT_VaporEnergyConversion(pRegion)
    END IF ! burning Used

#ifdef RFLO
    IF (pRegion%inrtInput%twoDAverage /= 0) THEN
      CALL INRT_TwoDimAverage(pRegion)
    END IF
#endif

    IF (pRegion%inrtInput%inrts(INRT_TYPE_BURNING)%used) THEN
      CALL INRT_BurnStatusUpdate(pRegion)
    END IF ! burning used

  END IF ! inrtUsed
#endif

#ifdef PEUL
  IF (global%peulUsed .AND. finalStage) THEN
    CALL PEUL_EnforcePositivity(pRegion)
  END IF ! peulUsed
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE AfterUpdateMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: AfterUpdateMP.F90,v $
! Revision 1.6  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/03/27 15:17:26  haselbac
! Fixed screwed-up last check-in
!
! Revision 1.3  2006/03/26 20:21:07  haselbac
! Changed to wrappers bcos of GL model
!
! Revision 1.2  2005/12/01 21:50:11  fnajjar
! Added call to PLAG_CheckValidity
!
! Revision 1.1  2004/12/01 16:47:41  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/08/04 00:26:49  wasistho
! call rflo_checkvalidity only at finalstage
!
! Revision 1.4  2004/07/26 19:02:36  wasistho
! add RFLO_CheckValidity
!
! Revision 1.3  2004/03/25 21:14:20  jferry
! changed AfterUpdate to call most subroutines only after final RK stage
!
! Revision 1.2  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
!******************************************************************************







