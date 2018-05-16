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
! Purpose: set RaNS boundary conditions or exchange data between adjacent
!          regions being on same processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%turb%cv = updated cv of RaNS eqs.
!                                        in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_rFLO_RansBndConditionsSet.F90,v 1.4 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_RansBndConditionsSet( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_FloRansExchangeDummyConf, &
                    TURB_FloRansBcondInflow,     TURB_FloRansBcondInjection, &
                    TURB_FloRansBcondNoslipWall, TURB_FloRansBcondSymmetry, &
                    TURB_FloRansBcondZeroGrad
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

  CALL RegisterFunction( global,'TURB_RFLO_RansBndConditionsSet',&
  'TURB_rFLO_RansBndConditionsSet.F90' )

  IF (regions(iReg)%turbInput%modelClass /= MODEL_RANS) GOTO 999

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

    IF (bcType>=BC_INFLOW .AND. bcType<=BC_INFLOW+BC_RANGE) THEN
      CALL TURB_FloRansBcondInflow( regions(iReg),patch )

! - outflow, slipwall, farfield

    ELSE IF ((bcType>=BC_OUTFLOW  .AND. bcType<=BC_OUTFLOW+BC_RANGE)  .OR. &
             (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) .OR. &
             (bcType>=BC_FARFIELD .AND. bcType<=BC_FARFIELD+BC_RANGE)) THEN
      CALL TURB_FloRansBcondZeroGrad( regions(iReg),patch )

! - conforming region interface

    ELSE IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid == global%myProcid) THEN
        CALL TURB_FloRansExchangeDummyConf( regions(iReg),regions(iRegSrc), &
                                            patch,patchSrc )
      ENDIF

! - non-conforming region interface (integer)

    ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      IF (regions(iRegSrc)%procid == global%myProcid) THEN
!        CALL TURB_FloRansExchangeDummyInt( regions(iReg),regions(iRegSrc), &
!                                           patch,patchSrc )
!      ENDIF
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'non-conforming integer bc is not ready yet for RaNS.' )

! - non-conforming region interface (irregular)

    ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      IF (regions(iRegSrc)%procid == global%myProcid) THEN
!        CALL TURB_FloRansExchangeDummyIreg( regions(iReg),regions(iRegSrc), &
!                                            patch,patchSrc )
!      ENDIF
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'non-conforming irregular bc is not ready yet for RaNS.' )

! - noslip wall

    ELSE IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
      CALL TURB_FloRansBcondNoslipWall( regions(iReg),patch )

! - injection

    ELSE IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
      CALL TURB_FloRansBcondInjection( regions(iReg),patch )

! - symmetry

    ELSE IF (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE) THEN
      CALL TURB_FloRansBcondSymmetry( regions(iReg),patch )

! - translational periodicity

    ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid == global%myProcid) THEN
        CALL TURB_FloRansExchangeDummyConf( regions(iReg),regions(iRegSrc), &
                                            patch,patchSrc )
      ENDIF

! - rotational periodicity

    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

!      CALL TURB_FloRansBcondRotatPeriod( regions(iReg),regions(iRegSrc), &
!                                         patch,patchSrc )
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__, &
                     'rotational periodic bc is not ready yet for RaNS.' )

    ELSE
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )
    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RFLO_RansBndConditionsSet

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_RansBndConditionsSet.F90,v $
! Revision 1.4  2008/12/06 08:44:45  mtcampbe
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
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







