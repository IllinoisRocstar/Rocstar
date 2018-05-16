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
! Purpose: set boundary conditions or exchange data between adjacent
!          regions being on same processor for smoke.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%peul = updated smoke values (cv,dv,tv)
!                                     in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_BoundaryConditionsSet.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_BoundaryConditionsSet( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE PEUL_ModInterfaces, ONLY : PEUL_BcondInflow,PEUL_BcondInjection, &
                                 PEUL_BcondOutflow,PEUL_BcondSlipWall, &
                                 PEUL_BcondSymmetry,PEUL_ExchangeDummyConf
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER    :: regions(:)
  INTEGER,        INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, nPatches, bcType, iRegSrc, iPatchSrc

  TYPE(t_patch),  POINTER :: patch, patchSrc
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_BoundaryConditionsSet.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PEUL_BoundaryConditionsSet',&
  'PEUL_BoundaryConditionsSet.F90' )

! begin -----------------------------------------------------------------------

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patch%bcType
    iRegSrc   = patch%srcRegion
    iPatchSrc = patch%srcPatch

! - inflow

    SELECT CASE(bcType)

    CASE (BC_INFLOW:BC_INFLOW+BC_RANGE)
      CALL PEUL_BcondInflow( regions(iReg),patch )

! - outflow

    CASE (BC_OUTFLOW:BC_OUTFLOW+BC_RANGE)
      CALL PEUL_BcondOutflow( regions(iReg),patch )

! - conforming region interface

    CASE (BC_REGIONCONF:BC_REGIONCONF+BC_RANGE)
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid == global%myProcid) THEN
        CALL PEUL_ExchangeDummyConf( regions(iReg),regions(iRegSrc), &
                                     patch,patchSrc )
      ENDIF

! - non-conforming region interface (integer)

!    CASE (BC_REGIONINT:BC_REGIONINT+BC_RANGE)
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      IF (regions(iRegSrc)%procid == global%myProcid) THEN
!        CALL PEUL_ExchangeDummyInt( regions(iReg),regions(iRegSrc), &
!                                    patch,patchSrc )
!      ENDIF

! - non-conforming region interface (irregular)

!    CASE (BC_REGNONCONF:BC_REGNONCONF+BC_RANGE)
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      IF (regions(iRegSrc)%procid == global%myProcid) THEN
!        CALL PEUL_ExchangeDummyIreg( regions(iReg),regions(iRegSrc), &
!                                     patch,patchSrc )
!      ENDIF

! - slip wall

    CASE (BC_SLIPWALL:BC_SLIPWALL+BC_RANGE)
      CALL PEUL_BcondSlipWall( regions(iReg),patch )

! - noslip wall

    CASE (BC_NOSLIPWALL:BC_NOSLIPWALL+BC_RANGE)
      CALL PEUL_BcondSlipWall( regions(iReg),patch )

! - far field

    CASE (BC_FARFIELD:BC_FARFIELD+BC_RANGE)
      CALL PEUL_BcondInflow( regions(iReg),patch )

! - injection

    CASE (BC_INJECTION:BC_INJECTION+BC_RANGE)
      CALL PEUL_BcondInjection( regions(iReg),patch )

! - symmetry

    CASE (BC_SYMMETRY:BC_SYMMETRY+BC_RANGE)
      CALL PEUL_BcondSymmetry( regions(iReg),patch )

! - translational periodicity

!    CASE (BC_TRA_PERI:BC_TRA_PERI+BC_RANGE)
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      IF (regions(iRegSrc)%procid == global%myProcid) THEN
!        CALL PEUL_ExchangeDummyConf( regions(iReg),regions(iRegSrc), &
!                                     patch,patchSrc )
!      ENDIF

! - rotational periodicity

!    CASE (BC_ROT_PERI:BC_ROT_PERI+BC_RANGE)
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      CALL PEUL_BcondRotatPeriod( regions(iReg),regions(iRegSrc), &
!                                  patch,patchSrc )

! - other

    CASE DEFAULT
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )

    END SELECT ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_BoundaryConditionsSet

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_BoundaryConditionsSet.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:23  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.4  2003/04/09 15:09:56  jferry
! added slip wall boundary conditions
!
! Revision 1.3  2003/04/09 14:23:36  fnajjar
! Activated PEUL_exchangeDummyConf for Multi-region non-MPI capabilities
!
! Revision 1.2  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







