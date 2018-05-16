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
! Purpose: update values in dummy cells for smoke at injection boundary patch.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%peul = smoke variables in dummy cells.
!
! Notes: if the mass flow rate is zero or negative, the corresponding
!        cell face is treated like solid wall.
!
!******************************************************************************
!
! $Id: PEUL_BcondInjection.F90,v 1.4 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_BcondInjection( region,patch )

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  TYPE(t_patch),  INTENT(IN)    :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, n1, n2

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, iLev
  INTEGER :: iCOff, ijCOff, nOff, ijkC, ijkC1, ijkD, mixtI2d, peulI2d
  INTEGER :: mixtDistrib, peulDistrib, dumExtrapol, nPtypes, lobound

  REAL(RFREAL) :: rhogas, mRate

  REAL(RFREAL), POINTER :: mixtCv(:,:), peulCv(:,:), vals(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_BcondInjection.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_BcondInjection',&
  'PEUL_BcondInjection.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers

  iLev   = region%currLevel
  lobound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  nOff        = ABS(patch%l1end-patch%l1beg) + 1
  dumExtrapol = patch%mixt%switches(BCSWI_INJECT_EXTRAP)
  mixtDistrib = patch%mixt%distrib
  peulDistrib = patch%peul%distrib

  mixtCv => region%levels(iLev)%mixt%cv
  peulCv => region%levels(iLev)%peul%cv
  vals   => patch%peul%vals

  nPtypes = region%peulInput%nPtypes
  IF (region%levels(iLev)%peul%nCv /= nPtypes) &
    CALL ErrorStop( global,ERR_PEUL_NPMISMATCH,__LINE__ )
  IF (LBOUND(vals,1) /= 1 .OR. UBOUND(vals,1) /= nPtypes) &
    CALL ErrorStop( global,ERR_PEUL_NPMISMATCH,__LINE__ )

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
          ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)
          ijkD  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)

          SELECT CASE(lobound)
          CASE (1:2)
            n1 = j - jbeg
            n2 = k - kbeg
          CASE (3:4)
            n1 = k - kbeg
            n2 = i - ibeg
          CASE(5:6)
            n1 = i - ibeg
            n2 = j - jbeg
          END SELECT ! lobound

          mixtI2d = mixtDistrib * IndIJ(n1,n2,nOff)
          peulI2d = peulDistrib * IndIJ(n1,n2,nOff)

          mRate   = vals(BCDAT_INJECT_MFRATE,mixtI2d)

! ------- surface burning and 1st dummy cell

          IF (mRate>0._RFREAL .AND. idum==1) THEN

            IF (dumExtrapol == EXTRAPOL_CONST) THEN
              rhogas = mixtCv(CV_MIXT_DENS,ijkD)
              peulCv(:,ijkD) = rhogas*vals(:,peulI2d)

            ELSE
              rhogas = 0.5_RFREAL*(mixtCv(CV_MIXT_DENS,ijkC) + &
                                   mixtCv(CV_MIXT_DENS,ijkD))
              peulCv(:,ijkD) = 2._RFREAL*rhogas*vals(:,peulI2d) - peulCv(:,ijkC)
            END IF ! dumExtrapol

! ------- not burning or dummy cell>1

          ELSE ! treat like slip wall

            IF (dumExtrapol == EXTRAPOL_CONST) THEN
              peulCv(:,ijkD) = peulCv(:,ijkC)
            ELSE
              peulCv(:,ijkD) = 2._RFREAL*peulCv(:,ijkC) - peulCv(:,ijkC1)
            END IF ! dumExtrapol

          END IF ! mRate .AND. idum

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_BcondInjection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_BcondInjection.F90,v $
! Revision 1.4  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:18  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/12/01 21:09:15  haselbac
! Initial revision after changing case
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/09 15:11:31  jferry
! removed dependence on temporary storage structures
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







