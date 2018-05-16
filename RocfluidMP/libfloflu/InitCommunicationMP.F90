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
! Purpose: initiate communication for RocfluidMP framework.
!
! Description: none.
!
! Input: regions = data of all regions,
!        iReg    = current region,
!        istage  = current RK stage.
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: InitCommunicationMP.F90,v 1.4 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE InitCommunicationMP( regions,iReg,istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
#ifdef RFLO
  USE RFLO_ModBoundaryConditions, ONLY : RFLO_BoundaryConditionsSend
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_BoundaryConditionsSend
#endif
#ifdef TURB
#ifdef RFLO
  USE ModInterfacesTurbulence, ONLY : TURB_RFLO_RansBndConditionsSend
#endif
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg, istage

! ... local variables
  REAL(RFREAL) :: ark(5), time, subtime, subdt

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'InitCommunicationMP',&
  'InitCommunicationMP.F90' )

! compute value of substage physical time and time-step =======================

  ark(:) = regions(iReg)%mixtInput%ark(:)

  IF (istage == 1) THEN
    subdt   = ark(1)*global%dtMin
    subtime = global%currentTime + subdt
  ELSE IF (istage == global%nrkSteps) THEN
    subdt   = (1.0_RFREAL - ark(global%nrkSteps-1))*global%dtMin
    subtime = global%currentTime + global%dtMin
  ELSE
    subdt   = (ark(istage) - ark(istage - 1))*global%dtMin
    subtime = global%currentTime + ark(istage)*global%dtMin
  ENDIF ! istage

! send conservative variables to other processors -----------------------------

#ifdef RFLO
  CALL RFLO_BoundaryConditionsSend( regions,iReg )

#ifdef PEUL
  IF ( global%peulUsed ) &
    CALL PEUL_BoundaryConditionsSend( regions,iReg )
#endif

#ifdef TURB
  IF (regions(iReg)%mixtInput%flowModel == FLOW_NAVST .AND. &
      regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) &
    CALL TURB_RFLO_RansBndConditionsSend( regions,iReg )
#endif
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE InitCommunicationMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: InitCommunicationMP.F90,v $
! Revision 1.4  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/28 22:49:25  wasistho
! moved RFLO_Bcond* and RFLO_BoundaryCond* routines into RFLO_ModBoundaryConditions
!
! Revision 1.1  2004/12/01 16:48:39  haselbac
! Initial revision after changing case
!
! Revision 1.14  2004/03/27 03:02:24  wasistho
! added ifdef RFLO within ifdef TURB
!
! Revision 1.13  2004/03/11 03:32:39  wasistho
! changed rocturb nomenclature
!
! Revision 1.12  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.11  2003/11/25 21:01:40  haselbac
! Removed RFLU routines, now in RFLU_UpdateDummyCells
!
! Revision 1.10  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/10/03 20:41:09  haselbac
! Removed call to UpdateTbc - wrong place
!
! Revision 1.6  2003/10/03 20:12:22  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.5  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.4  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/09 14:17:23  fnajjar
! Added PEUL_BoundaryConditionsSend for MPI-based rocsmoke
!
! Revision 1.1  2003/03/28 19:45:59  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************







