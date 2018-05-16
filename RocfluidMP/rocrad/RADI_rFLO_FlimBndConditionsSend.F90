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
! Purpose: send FLD radiation data to dummy cells of adjacent regions being 
!          located on another processor.
!
! Description: none.
!
! Input: regions = FLD radiation data of all regions
!        iReg    = index of current region.
!
! Output: FLD radiation data to other processors.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_rFLO_FlimBndConditionsSend.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_RFLO_FlimBndConditionsSend( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE RADI_ModInterfaces, ONLY : RADI_FloFlimSendDummyConf, &
                                 RADI_FloFlimSendCornEdgeCells
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
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

  CALL RegisterFunction( global,'RADI_RFLO_FlimBndConditionsSend',&
  'RADI_rFLO_FlimBndConditionsSend.F90' )

  IF (regions(iReg)%radiInput%radiModel /= RADI_MODEL_FLDTRAN) GOTO 999

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! send data for edge and corner cells -----------------------------------------

  CALL RADI_FloFlimSendCornEdgeCells( regions,iReg )

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
        CALL RADI_FloFlimSendDummyConf( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - non-conforming region interface (integer)

    ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
!        CALL RADI_FloFlimSendDummyInt( regions(iReg),regions(iRegSrc),patch )
!      ENDIF
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'non-conforming integer bc is not ready yet for RADI.' )

! - non-conforming region interface (irregular)

    ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
!        CALL RADI_FloFlimSendDummyIreg( regions(iReg),regions(iRegSrc),patch )
!      ENDIF
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'non-conforming irregular bc is not ready yet for RADI.' )

! - translational periodicity

    ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
        CALL RADI_FloFlimSendDummyConf( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - rotational periodicity

    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      CALL RADI_FloFlimSendDummyRot( regions(iReg),regions(iRegSrc), &
!                                     patch,patchSrc )
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'rotational periodic bc is not ready yet for RADI.' )

    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_RFLO_FlimBndConditionsSend

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_rFLO_FlimBndConditionsSend.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







