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
! Purpose: update RaNS variables in dummy cells at injection patch
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%turb%cv = RaNS variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_floRansBcondInjection.F90,v 1.6 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansBcondInjection( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, l

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, lbound, n1, n2, i2d, distrib, nOff
  INTEGER :: iCOff, ijCOff, ijkC, ijkC1, ijkC2, ijkD, nCv

  REAL(RFREAL) :: mRate
  REAL(RFREAL), POINTER :: tcv(:,:), vals(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_FloRansBcondInjection',&
  'TURB_floRansBcondInjection.F90' )

! get dimensions and pointers

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  nOff    =  ABS(patch%l1end-patch%l1beg) + 1
  distrib =  patch%mixt%distrib
  nCv     =  region%turbInput%nCv
  tcv     => region%levels(iLev)%turb%cv
  vals    => patch%mixt%vals

! loop over all cells of a patch

  IF ((region%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN

    DO idum=1,region%nDumCells
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkD  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
            ijkC  = &
            IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
            ijkC1 = &
            IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)
            ijkC2 = &
            IndIJK(i+(idum-1)*idir,j+(idum-1)*jdir,k+(idum-1)*kdir,iCOff,ijCOff)

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
            ENDIF
            i2d   = distrib * IndIJ(n1,n2,nOff)
            mRate = vals(BCDAT_INJECT_MFRATE,i2d)

            IF (mRate > 0._RFREAL) THEN  ! burning surface
              IF (idum == 1) THEN
                DO l=1,nCv
                  tcv(l,ijkD) = tcv(l,ijkC)
                ENDDO
              ELSE
                DO l=1,nCv
                  tcv(l,ijkD) = 2._RFREAL*tcv(l,ijkC) - tcv(l,ijkC1)
                ENDDO
              ENDIF
            ELSE                       ! non burning surface, as ns-wall
              DO l=1,nCv
                tcv(l,ijkD) = - tcv(l,ijkC2)
              ENDDO
            ENDIF
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDDO        ! idum

  ELSE

    CALL ErrorStop( region%global,ERR_TURB_RANSINPUT,__LINE__, &
                   'injection bc is not ready yet for RaNS model selected' )
  ENDIF

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_FloRansBcondInjection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansBcondInjection.F90,v $
! Revision 1.6  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:40:44  mparmar
! Renamed patch variables
!
! Revision 1.3  2005/03/09 06:36:53  wasistho
! incorporated HDESSA
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2004/01/21 03:44:04  wasistho
! fixed distrib and ijkC, C1, C2, D
!
! Revision 1.1  2003/10/07 02:17:02  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







