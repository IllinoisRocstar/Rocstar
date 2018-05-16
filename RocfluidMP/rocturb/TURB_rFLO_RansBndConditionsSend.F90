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
! Purpose: send RaNS data to dummy cells of adjacent regions being located
!          on another processor.
!
! Description: none.
!
! Input: regions = RaNS data of all regions
!        iReg    = index of current region.
!
! Output: RaNS data to other processors.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_rFLO_RansBndConditionsSend.F90,v 1.4 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_RansBndConditionsSend( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE TURB_ModInterfaces, ONLY : TURB_FloRansSendDummyConf, &
                                 TURB_FloRansSendCornEdgeCells
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

  TYPE(t_patch), POINTER  :: patch, patchSrc
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'TURB_RFLO_RansBndConditionsSend',&
  'TURB_rFLO_RansBndConditionsSend.F90' )

  IF (regions(iReg)%turbInput%modelClass /= MODEL_RANS) GOTO 999

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! send data for edge and corner cells -----------------------------------------

  CALL TURB_FloRansSendCornEdgeCells( regions,iReg )

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
        CALL TURB_FloRansSendDummyConf( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - non-conforming region interface (integer)

    ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
!        CALL TURB_FloRansSendDummyInt( regions(iReg),regions(iRegSrc),patch )
!      ENDIF
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'non-conforming integer bc is not ready yet for RaNS.' )

! - non-conforming region interface (irregular)

    ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
!        CALL TURB_FloRansSendDummyIreg( regions(iReg),regions(iRegSrc),patch )
!      ENDIF
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'non-conforming irregular bc is not ready yet for RaNS.' )

! - translational periodicity

    ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
        CALL TURB_FloRansSendDummyConf( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - rotational periodicity

    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      CALL TURB_FloRansSendDummyRot( regions(iReg),regions(iRegSrc), &
!                                     patch,patchSrc )
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'rotational periodic bc is not ready yet for RaNS.' )

    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RFLO_RansBndConditionsSend

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_RansBndConditionsSend.F90,v $
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
! Revision 1.2  2004/01/23 00:32:07  wasistho
! added sending RaNS edge/corner variables
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







