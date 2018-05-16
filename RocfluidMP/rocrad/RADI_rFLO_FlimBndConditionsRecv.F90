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
! Purpose: receive FLD radiation data for dummy cells from adjacent regions 
!          being on another processor.
!
! Description: none.
!
! Input: regions = radiation data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%radi = updated radiation cv
!                                     in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_rFLO_FlimBndConditionsRecv.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_RFLO_FlimBndConditionsRecv( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE RADI_ModInterfaces, ONLY : RADI_FloFlimRecvDummyVals, &
               RADI_FloFlimCorrCornEdgeCells, RADI_FloFlimExchCornEdgeCells, &
               RADI_FloFlimRecvCornEdgeCells, RADI_FloFlimSetCornEdgeCells 
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

  TYPE(t_patch), POINTER :: patch, patchSrc

!******************************************************************************

  CALL RegisterFunction( regions(iReg)%global, &
                         'RADI_RFLO_FlimBndConditionsRecv', 'RADI_rFLO_FlimBndConditionsRecv.F90' )

  IF (regions(iReg)%radiInput%radiModel /= RADI_MODEL_FLDTRAN) GOTO 999

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

        CALL RADI_FloFlimRecvDummyVals( regions(iReg),regions(iRegSrc), &
                                        patch,patchSrc )
      ENDIF
    ENDIF

  ENDDO   ! iPatch

! copy/exchange data for edge and corner cells --------------------------------

! set corner/edge cells outcommented due to extrapolation in tcv update
!  CALL RADI_FloFlimSetCornEdgeCells( regions(iReg) )  

  DO iPatch=1,nPatches
    patch => regions(iReg)%levels(iLev)%patches(iPatch)
    bcType = patch%bcType
    IF ((bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
        (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE .AND. &
         regions(iReg)%mixtInput%flowModel==FLOW_NAVST) .OR. &
        (bcType>=BC_SYMMETRY   .AND. bcType<=BC_SYMMETRY  +BC_RANGE)) THEN
      CALL RADI_FloFlimCorrCornEdgeCells( regions(iReg),patch,bcType )
    ENDIF
  ENDDO   ! iPatch

  CALL RADI_FloFlimExchCornEdgeCells( regions,iReg )
  CALL RADI_FloFlimRecvCornEdgeCells( regions,iReg )

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( regions(iReg)%global )

END SUBROUTINE RADI_RFLO_FlimBndConditionsRecv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_rFLO_FlimBndConditionsRecv.F90,v $
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







