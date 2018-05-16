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
! Purpose: find source patches for inter-region communication.
!
! Description: none.
!
! Input: regions = dimensions and topology of all regions.
!
! Output: regions%levels(1)%patches = source patches.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_FindSourcePatches.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_FindSourcePatches( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModParameters
  USE ModError
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, iPatchSrc

! ... local variables
  CHARACTER(CHRLEN) :: msg

  INTEGER :: bcType, iRegSrc, lbs, l1bs, l1es, l2bs, l2es

  TYPE(t_patch), POINTER  :: patch, patchSrc
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_FindSourcePatches',&
  'RFLO_FindSourcePatches.F90' )

! loop over regions and patches

  DO iReg=1,global%nRegions
    DO iPatch=1,regions(iReg)%nPatches

      patch  => regions(iReg)%levels(1)%patches(iPatch)
      bcType = patch%bcType

      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
          (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
        iRegSrc = patch%srcRegion
        lbs     = patch%srcLbound
        l1bs    = MIN(patch%srcL1beg,patch%srcL1end)
        l1es    = MAX(patch%srcL1beg,patch%srcL1end)
        l2bs    = MIN(patch%srcL2beg,patch%srcL2end)
        l2es    = MAX(patch%srcL2beg,patch%srcL2end)

! ----- search for corresponding source patch on source region

        patch%srcPatch = -999

        DO iPatchSrc=1,regions(iRegSrc)%nPatches
          IF (iRegSrc/=iReg .OR. iPatchSrc/=iPatch) THEN
            patchSrc => regions(iRegSrc)%levels(1)%patches(iPatchSrc)
            IF (patchSrc%bcType==bcType .AND. &
                patchSrc%lbound==lbs    .AND. &
                patchSrc%l1beg ==l1bs   .AND. &
                patchSrc%l1end ==l1es   .AND. &
                patchSrc%l2beg ==l2bs   .AND. &
                patchSrc%l2end ==l2es) THEN       ! OK, patch found
              patch%srcPatch = iPatchSrc
            ENDIF
          ENDIF
        ENDDO    ! iPatchSrc

        IF (patch%srcPatch < 0) THEN
          WRITE(msg,1000) iReg,iPatch,bcType
          CALL ErrorStop( global,ERR_PATCH_NOSOURCE,__LINE__,msg )
        ENDIF
      ENDIF   ! bcType

    ENDDO     ! iPatch
  ENDDO       ! iReg

! finalize

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', patch ',I3,', BC type ',I3)

END SUBROUTINE RFLO_FindSourcePatches

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_FindSourcePatches.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.9  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.8  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.7  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.5  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.3  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.2  2001/12/19 23:09:22  jblazek
! Added routines to read grid and solution.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************







