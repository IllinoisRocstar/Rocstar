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
! Purpose: check BC data specified by the user for all modules.
!
! Description: none.
!
! Input: regions = input parameters for all grid regions.
!
! Output: none.
!
! Notes: also allocates TBC arrays.
!
!******************************************************************************
!
! $Id: RFLO_CheckBcInput.F90,v 1.6 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CheckBcInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModMixture,    ONLY : t_mixt_input
  USE ModBndPatch,   ONLY : t_patch, t_bcvalues
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, iLev, fType

! ... local variables
  CHARACTER(CHRLEN) :: msg

  INTEGER :: errorFlag

  TYPE(t_global),     POINTER :: global
  TYPE(t_mixt_input), POINTER :: input
  TYPE(t_patch),      POINTER :: patch
  TYPE(t_bcvalues),   POINTER :: bc

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_CheckBcInput',&
  'RFLO_CheckBcInput.F90' )

! check if all BCs set (for all active modules), and allocate TBCs

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      input => regions(iReg)%mixtInput
      DO iLev=1,regions(iReg)%nGridLevels
        DO iPatch=1,regions(iReg)%nPatches
          patch => regions(iReg)%levels(iLev)%patches(iPatch)
          DO fType = 1,FTYPE_MAX

            SELECT CASE(fType)
            CASE (FTYPE_MIXT)
              bc => patch%mixt
            CASE (FTYPE_TURB)
              IF (input%turbModel == TURB_MODEL_NONE) CYCLE
              bc => patch%turb
            CASE (FTYPE_PLAG)
              CYCLE
            CASE (FTYPE_PEUL)
              IF (global%peulUsed .eqv. .false.) CYCLE
              bc => patch%peul
            CASE (FTYPE_SPEC)
              IF (input%gasModel == GAS_MODEL_TCPERF) CYCLE
              bc => patch%spec
            CASE (FTYPE_RADI)
              IF (input%radiUsed .eqv. .false.) CYCLE
              bc => patch%valRadi
            END SELECT

            IF (bc%bcSet .eqv. .false.) THEN
              WRITE(msg,1000) iReg,iLev,iPatch,patch%bcType,fType
              CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__,TRIM(msg) )
            ENDIF

            IF (bc%nData > 0) THEN
              ALLOCATE( bc%tbcs(bc%nData),stat=errorFlag )
              global%error = errorFlag
              IF (global%error /= 0) &
                CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
              bc%tbcs(:)%tbcType = TBC_NONE
            ENDIF

          ENDDO ! fType
        ENDDO ! iPatch
      ENDDO ! iLev
    ENDIF ! region active and on my processor
  ENDDO ! iReg

! finalize

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', level ',I1,', patch ',I3,', bcType ',I3, &
            ', fType ',I3)

END SUBROUTINE RFLO_CheckBcInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CheckBcInput.F90,v $
! Revision 1.6  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2008/08/13 17:25:33  mtcampbe
! Fixed IF statements for booleans to use .eqv. .false. instead of .NOT.
! (beats me)
!
! Revision 1.3  2006/08/19 15:39:30  mparmar
! Renamed patch variables
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.4  2004/03/05 22:09:02  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.3  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.2  2003/02/26 23:19:48  jferry
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.1  2003/02/11 22:49:53  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.1  2003/02/11 22:30:20  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
!******************************************************************************







