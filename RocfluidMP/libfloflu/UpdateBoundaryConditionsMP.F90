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
! Purpose: update boundary conditions for RocfluidMP framework.
!
! Description: none.
!
! Input: regions = data of all regions,
!        istage  = current RK stage.
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: UpdateBoundaryConditionsMP.F90,v 1.6 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE UpdateBoundaryConditionsMP( regions,istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
#ifdef RFLO
  USE RFLO_ModBoundaryConditions, ONLY : RFLO_BoundaryConditionsSet, &
        RFLO_BoundaryConditionsRecv
  USE ModInterfaces, ONLY : RFLO_ClearSendRequests, &
        RFLO_GetBoundaryValues, RFLO_SendBoundaryValues,     &
        RFLO_SendBoundaryValuesAlpha, UpdateTbc

#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_RFLO_RansBndConditionsSet, &
        TURB_RFLO_RansBndConditionsRecv, TURB_RFLO_RansClearSendRequests
#endif
!#ifdef PERI
!  USE PERI_ModHybridDES, ONLY : PERI_CoMeanCorrection
!#endif
#endif

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_PatchUpdateWrapper, &
                                      PLAG_CECellsWrapper
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_BoundaryConditionsSet, &
        PEUL_BoundaryConditionsRecv, PEUL_ClearSendRequests
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN) :: istage

! ... loop variables
  INTEGER :: iReg

! ... local variables
  REAL(RFREAL) :: ark(5), trk(5), time, subdt, subtime

  TYPE(t_global), POINTER :: global

#ifdef GENX
  DOUBLE PRECISION :: alpha
#endif

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'UpdateBoundaryConditions',&
  'UpdateBoundaryConditionsMP.F90' )

! get time-stepping coefficients ==============================================

  ark(:) = regions(1)%mixtInput%ark(:)
  trk(:) = regions(1)%mixtInput%trk(:)

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

! update boundary conditions ------------------------------------------------

#ifdef RFLO
#ifdef GENX
  time  = global%currentTime + global%dtMin*trk(istage)
  alpha = (time-global%timeStamp)/global%dTimeSystem

! - send fluids density at interface
  CALL COM_call_function( global%genxHandleBc,2,alpha,1 )

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      IF (regions(iReg)%mixtInput%externalBc) THEN
        CALL RFLO_SendBoundaryValuesAlpha( regions(iReg) )
      ENDIF
    ENDIF
  ENDDO

! - get BC values at the interface; set BC values in dummy cells

  CALL COM_call_function( global%genxHandleBc,2,alpha,2 )
#endif

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      IF (regions(iReg)%mixtInput%externalBc) THEN
        CALL RFLO_GetBoundaryValues( regions(iReg) )
      ENDIF
      CALL UpdateTbc( regions(iReg),subtime,subdt,istage==global%nrkSteps )
      CALL RFLO_BoundaryConditionsSet( regions,iReg )
#ifdef PEUL
      IF (global%peulUsed) &
        CALL PEUL_BoundaryConditionsSet( regions,iReg )
#endif
#ifdef TURB
      IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
          (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) &
        CALL TURB_RFLO_RansBndConditionsSet( regions,iReg )
#endif
    ENDIF
  ENDDO

! - receive variables from other processors; send BC data to external driver --

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      CALL RFLO_BoundaryConditionsRecv( regions,iReg )
#ifdef PEUL
      IF (global%peulUsed) &
        CALL PEUL_BoundaryConditionsRecv( regions,iReg )
#endif
#ifdef TURB
      IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
          (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) &
        CALL TURB_RFLO_RansBndConditionsRecv( regions,iReg )
#endif
    ENDIF
  ENDDO

! - clear send requests -------------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      CALL RFLO_ClearSendRequests( regions,iReg,.false. )
#ifdef PEUL
      IF (global%peulUsed) &
        CALL PEUL_ClearSendRequests( regions,iReg )
#endif
#ifdef TURB
      IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
          (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) &
        CALL TURB_RFLO_RansClearSendRequests( regions,iReg )
#endif
    ENDIF
  ENDDO

#ifdef PLAG
  CALL PLAG_CECellsWrapper( regions )

  CALL PLAG_PatchUpdateWrapper( regions )
#endif

! - post-BC operations of physical modules ------------------------------------

#ifdef PERI
!  DO iReg=1,global%nRegions
!    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
!        regions(iReg)%active==ACTIVE) THEN              ! on my processor
!      IF (regions(iReg)%periInput%flowKind /= OFF) THEN
!        CALL PERI_CoMeanCorrection( regions(iReg) ) 
!      ENDIF
!    ENDIF
!  ENDDO
#endif

#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE UpdateBoundaryConditionsMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateBoundaryConditionsMP.F90,v $
! Revision 1.6  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/04/17 05:21:22  wasistho
! commented PERI_coMeanCorrection
!
! Revision 1.3  2005/04/06 02:16:59  wasistho
! mv call to PERI_CoMeanCorrection to UpdateBoundaryConditionsMP
!
! Revision 1.2  2004/12/28 22:49:37  wasistho
! moved RFLO_Bcond* and RFLO_BoundaryCond* routines into RFLO_ModBoundaryConditions
!
! Revision 1.1  2004/12/01 16:51:56  haselbac
! Initial revision after changing case
!
! Revision 1.21  2004/03/27 03:02:36  wasistho
! added ifdef RFLO within ifdef TURB
!
! Revision 1.20  2004/03/11 03:32:25  wasistho
! changed rocturb nomenclature
!
! Revision 1.19  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.18  2004/02/26 21:01:48  haselbac
! Removed istage argument from call to PLAG_PatchUpdateWrapper
!
! Revision 1.17  2004/02/14 01:51:08  fnajjar
! Removed WRITE statements surrounding PLAG_CECellsWrapper
!
! Revision 1.16  2004/02/02 22:48:37  haselbac
! Added ifdef RFLO - temporary measure
!
! Revision 1.15  2004/01/26 23:38:18  fnajjar
! Reinstated call to PLAG_CECellsWrapper
!
! Revision 1.14  2004/01/16 21:14:26  fnajjar
! Deactivated call to PLAG_CECellsWrapper as routine not complete
!
! Revision 1.13  2003/11/21 22:33:33  fnajjar
! Turned off comments
!
! Revision 1.12  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/11/12 21:19:11  fnajjar
! Added Corner-Edge cell routine
!
! Revision 1.8  2003/10/03 20:13:12  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.7  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.6  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.4  2003/05/09 17:01:03  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.3  2003/04/10 01:22:41  jblazek
! Got rid of pRegion in ViscousFluxesMP.
!
! Revision 1.2  2003/04/09 14:19:29  fnajjar
! Added routines to receive and clear requests for MPI-based rocsmoke
!
! Revision 1.1  2003/03/28 19:48:05  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************







