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
! Purpose: receive data for dummy cells from adjacent regions being
!          on another processor specific to PEUL Module.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%peul%cv =
!           updated smoke concentrations (cv) in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_BoundaryConditionsRecv.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_BoundaryConditionsRecv( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  USE PEUL_ModInterfaces, ONLY : PEUL_ReceiveDummyVals, &
        PEUL_SetCornerEdgeCells, PEUL_ExchangeCornerEdgeCells, &
        PEUL_ReceiveCornerEdgeCells, PEUL_CorrectCornerEdgeCells
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
    '$RCSfile: PEUL_BoundaryConditionsRecv.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PEUL_BoundaryConditionsRecv', 'PEUL_BoundaryConditionsRecv.F90' )

! begin -----------------------------------------------------------------------

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! receive data (regular cells) from other processors --------------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patch%bcType
    iRegSrc   = patch%srcRegion
    iPatchSrc = patch%srcPatch

! - region interface, periodic boundary

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
        (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
        (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
      IF (regions(iRegSrc)%procid /= regions(iReg)%global%myProcid) THEN
        patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

        CALL PEUL_ReceiveDummyVals( regions(iReg),regions(iRegSrc), &
                                    patch,patchSrc )
      ENDIF
    ENDIF

  ENDDO   ! iPatch

! copy/exchange data for edge and corner cells --------------------------------

  CALL PEUL_SetCornerEdgeCells( regions(iReg) )

  DO iPatch=1,nPatches
    patch => regions(iReg)%levels(iLev)%patches(iPatch)
    bcType = patch%bcType
    IF (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE) THEN
      CALL PEUL_CorrectCornerEdgeCells( regions(iReg),patch,bcType )
    ENDIF
  ENDDO   ! iPatch

  CALL PEUL_ExchangeCornerEdgeCells( regions,iReg )
  CALL PEUL_ReceiveCornerEdgeCells(  regions,iReg )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_BoundaryConditionsRecv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_BoundaryConditionsRecv.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:20  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
! Revision 1.1  2003/04/09 14:32:53  fnajjar
! Initial Import of MPI-based rocsmoke
!
!******************************************************************************







