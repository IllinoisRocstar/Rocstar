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
! Purpose: send data to dummy cells of adjacent regions being located
!          on another processor specific to PEUL module.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: data to other processors.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_BoundaryConditionsSend.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_BoundaryConditionsSend( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  USE PEUL_ModInterfaces, ONLY : PEUL_SendDummyConf,PEUL_SendCornerEdgeCells
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
    '$RCSfile: PEUL_BoundaryConditionsSend.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PEUL_BoundaryConditionsSend',&
  'PEUL_BoundaryConditionsSend.F90' )

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! send data for edge and corner cells -----------------------------------------

  CALL PEUL_SendCornerEdgeCells( regions,iReg )

! loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patch%bcType
    iRegSrc   = patch%srcRegion
    iPatchSrc = patch%srcPatch

! - conforming region interface

    IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
        CALL PEUL_SendDummyConf( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - non-conforming region interface (integer)

!    ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
!        CALL PEUL_SendDummyInt( regions(iReg),regions(iRegSrc),patch )
!      ENDIF

! - non-conforming region interface (irregular)

!    ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
!        CALL PEUL_SendDummyIreg( regions(iReg),regions(iRegSrc),patch )
!      ENDIF

! - translational periodicity

!    ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
!        CALL PEUL_SendDummyConf( regions(iReg),regions(iRegSrc),patch )
!      ENDIF

! - rotational periodicity

!    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
!      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
!
!      CALL PEUL_BcondRotatPeriod( regions(iReg),regions(iRegSrc), &
!                                  patch,patchSrc )

    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_BoundaryConditionsSend

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_BoundaryConditionsSend.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:22  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
! Revision 1.1  2003/04/09 14:33:15  fnajjar
! Initial Import of MPI-based rocsmoke
!
!******************************************************************************







