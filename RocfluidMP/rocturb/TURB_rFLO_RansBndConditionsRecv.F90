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
! Purpose: receive RaNS data for dummy cells from adjacent regions being
!          on another processor.
!
! Description: none.
!
! Input: regions = RaNS data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%turb = updated flow RaNS cv
!                                     in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_rFLO_RansBndConditionsRecv.F90,v 1.4 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_RansBndConditionsRecv( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE TURB_ModInterfaces, ONLY : TURB_FloRansRecvDummyVals, &
               TURB_FloRansCorrCornEdgeCells, TURB_FloRansExchCornEdgeCells, &
               TURB_FloRansRecvCornEdgeCells, TURB_FloRansSetCornEdgeCells 
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  INTEGER :: iLev, nPatches, bcType, iRegSrc, iPatchSrc

  TYPE(t_patch), POINTER :: patch, patchSrc

!******************************************************************************

  CALL RegisterFunction( regions(iReg)%global, &
                         'TURB_RFLO_RansBndConditionsRecv', 'TURB_rFLO_RansBndConditionsRecv.F90' )

  IF (regions(iReg)%turbInput%modelClass /= MODEL_RANS) GOTO 999

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

        CALL TURB_FloRansRecvDummyVals( regions(iReg),regions(iRegSrc), &
                                        patch,patchSrc )
      ENDIF
    ENDIF

  ENDDO   ! iPatch

! copy/exchange data for edge and corner cells --------------------------------

! set corner/edge cells outcommented due to extrapolation in tcv update
!  CALL TURB_FloRansSetCornEdgeCells( regions(iReg) )  

  DO iPatch=1,nPatches
    patch => regions(iReg)%levels(iLev)%patches(iPatch)
    bcType = patch%bcType
    IF ((bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
        (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE .AND. &
         regions(iReg)%mixtInput%flowModel==FLOW_NAVST) .OR. &
        (bcType>=BC_SYMMETRY   .AND. bcType<=BC_SYMMETRY  +BC_RANGE)) THEN
      CALL TURB_FloRansCorrCornEdgeCells( regions(iReg),patch,bcType )
    ENDIF
  ENDDO   ! iPatch

  CALL TURB_FloRansExchCornEdgeCells( regions,iReg )
  CALL TURB_FloRansRecvCornEdgeCells( regions,iReg )

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( regions(iReg)%global )

END SUBROUTINE TURB_RFLO_RansBndConditionsRecv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_RansBndConditionsRecv.F90,v $
! Revision 1.4  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/11 03:26:33  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2004/01/23 00:32:43  wasistho
! added receiving RaNS edge/corner variables
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







