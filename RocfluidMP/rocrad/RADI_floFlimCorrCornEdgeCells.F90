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
! Purpose: correct values in the corner and edge cells at diffuse, refracting,
!          etc boundaries.
!
! Description: none.
!
! Input: region    = region dimensions, user input
!        patch     = current patch
!        radBcType = type of boundary condition.
!
! Output: radiation Energy in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_floFlimCorrCornEdgeCells.F90,v 1.4 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimCorrCornEdgeCells( region,patch,radBcType )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetDimensPhys, RFLO_GetDimensDummy
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: radBcType

! ... loop variables
  INTEGER :: idum, i, j, k, l, n1, n2, iside

! ... local variables
  LOGICAL :: doDummy(4)

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ibegp, iendp, jbegp, jendp, kbegp, kendp, idir, jdir, kdir
  INTEGER :: ibeg(4), iend(4), jbeg(4), jend(4), kbeg(4), kend(4), nsides
  INTEGER :: iLev, nCv, lbound, iCOff, ijCOff, nOff, distrib
  INTEGER :: ijkC, ijkD, i2d

  REAL(RFREAL)          :: mRate
  REAL(RFREAL), POINTER :: rcv(:,:), vals(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FloFlimCorrCornEdgeCells', &
                         'RADI_floFlimCorrCornEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound
  nCv    = region%radiInput%nCv

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibegp,iendp,jbegp,jendp,kbegp,kendp )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib

  rcv  => region%levels(iLev)%radi%cv
  vals => patch%mixt%vals

! set loop indices ------------------------------------------------------------

  doDummy(:) = .false.
  nsides     = 0

  IF      (lbound==1 .OR. lbound==2) THEN
    ibeg(:) = ibegp
    iend(:) = iendp
    IF (jbegp == jpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jdcbeg
      jend(nsides)    = jpcbeg - 1
      IF (kbegp == kpcbeg) THEN
        kbeg(nsides)  = kdcbeg
      ELSE
        kbeg(nsides)  = kbegp
      ENDIF
      IF (kendp == kpcend) THEN
        kend(nsides)  = kdcend
      ELSE
        kend(nsides)  = kendp
      ENDIF
    ENDIF
    IF (jendp == jpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jpcend + 1
      jend(nsides)    = jdcend
      IF (kbegp == kpcbeg) THEN
        kbeg(nsides)  = kdcbeg
      ELSE
        kbeg(nsides)  = kbegp
      ENDIF
      IF (kendp == kpcend) THEN
        kend(nsides)  = kdcend
      ELSE
        kend(nsides)  = kendp
      ENDIF
    ENDIF
    IF (kbegp == kpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jbegp
      jend(nsides)    = jendp
      kbeg(nsides)    = kdcbeg
      kend(nsides)    = kpcbeg - 1
    ENDIF
    IF (kendp == kpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jbegp
      jend(nsides)    = jendp
      kbeg(nsides)    = kpcend + 1
      kend(nsides)    = kdcend
    ENDIF

  ELSE IF (lbound==3 .OR. lbound==4) THEN
    jbeg(:) = jbegp
    jend(:) = jendp
    IF (kbegp == kpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kdcbeg
      kend(nsides)    = kpcbeg - 1
      IF (ibegp == ipcbeg) THEN
        ibeg(nsides)  = idcbeg
      ELSE
        ibeg(nsides)  = ibegp
      ENDIF
      IF (iendp == ipcend) THEN
        iend(nsides)  = idcend
      ELSE
        iend(nsides)  = iendp
      ENDIF
    ENDIF
    IF (kendp == kpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kpcend + 1
      kend(nsides)    = kdcend
      IF (ibegp == ipcbeg) THEN
        ibeg(nsides)  = idcbeg
      ELSE
        ibeg(nsides)  = ibegp
      ENDIF
      IF (iendp == ipcend) THEN
        iend(nsides)  = idcend
      ELSE
        iend(nsides)  = iendp
      ENDIF
    ENDIF
    IF (ibegp == ipcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kbegp
      kend(nsides)    = kendp
      ibeg(nsides)    = idcbeg
      iend(nsides)    = ipcbeg - 1
    ENDIF
    IF (iendp == ipcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kbegp
      kend(nsides)    = kendp
      ibeg(nsides)    = ipcend + 1
      iend(nsides)    = idcend
    ENDIF

  ELSE IF (lbound==5 .OR. lbound==6) THEN
    kbeg(:) = kbegp
    kend(:) = kendp
    IF (ibegp == ipcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = idcbeg
      iend(nsides)    = ipcbeg - 1
      IF (jbegp == jpcbeg) THEN
        jbeg(nsides)  = jdcbeg
      ELSE
        jbeg(nsides)  = jbegp
      ENDIF
      IF (jendp == jpcend) THEN
        jend(nsides)  = jdcend
      ELSE
        jend(nsides)  = jendp
      ENDIF
    ENDIF
    IF (iendp == ipcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = ipcend + 1
      iend(nsides)    = idcend
      IF (jbegp == jpcbeg) THEN
        jbeg(nsides)  = jdcbeg
      ELSE
        jbeg(nsides)  = jbegp
      ENDIF
      IF (jendp == jpcend) THEN
        jend(nsides)  = jdcend
      ELSE
        jend(nsides)  = jendp
      ENDIF
    ENDIF
    IF (jbegp == jpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = ibegp
      iend(nsides)    = iendp
      jbeg(nsides)    = jdcbeg
      jend(nsides)    = jpcbeg - 1
    ENDIF
    IF (jendp == jpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = ibegp
      iend(nsides)    = iendp
      jbeg(nsides)    = jpcend + 1
      jend(nsides)    = jdcend
    ENDIF
  ENDIF

! loop over all cells of a patch ----------------------------------------------

  DO iside=1,nsides
    IF (doDummy(iside)) THEN

      DO idum=1,region%nDumCells
        DO k=kbeg(iside),kend(iside)
          DO j=jbeg(iside),jend(iside)
            DO i=ibeg(iside),iend(iside)
              ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
              ijkC = IndIJK(i+(idum-1)*idir,j+(idum-1)*jdir,k+(idum-1)*kdir,iCOff,ijCOff)

! ----------- noslip wall

              IF (radBcType>=RADI_BC_DIFFUS .AND. &
                  radBcType<=RADI_BC_DIFFUS+RADI_BC_RANGE) THEN
                DO l=1,nCv
                  rcv(l,ijkD) = -rcv(l,ijkC)
                ENDDO

! ----------- symmetry boundary

              ELSE IF (radBcType>=RADI_BC_REFRAC .AND. &
                       radBcType<=RADI_BC_REFRAC+RADI_BC_RANGE) THEN
                DO l=1,nCv
                  rcv(l,ijkD) =  rcv(l,ijkC)
                ENDDO

! ----------- non-burning injection boundary (like noslip wall)

              ELSE IF (radBcType>=RADI_BC_CYCLIC .AND. &
                       radBcType<=RADI_BC_CYCLIC+RADI_BC_RANGE) THEN
                IF      (lbound==1 .OR. lbound==2) THEN
                  n1 = MAX(MIN(j,jendp),jbegp) - jbegp
                  n2 = MAX(MIN(k,kendp),kbegp) - kbegp
                ELSE IF (lbound==3 .OR. lbound==4) THEN
                  n1 = MAX(MIN(k,kendp),kbegp) - kbegp
                  n2 = MAX(MIN(i,iendp),ibegp) - ibegp
                ELSE IF (lbound==5 .OR. lbound==6) THEN
                  n1 = MAX(MIN(i,iendp),ibegp) - ibegp
                  n2 = MAX(MIN(j,jendp),jbegp) - jbegp
                ENDIF
                i2d   = distrib * IndIJ(n1,n2,nOff)
                mRate = vals(BCDAT_INJECT_MFRATE,i2d)
                IF (mRate <= 0._RFREAL) THEN
                  DO l=1,nCv
                    rcv(l,ijkD) = -rcv(l,ijkC)
                  ENDDO
                ENDIF
              ENDIF

            ENDDO  ! i
          ENDDO    ! j
        ENDDO      ! k
      ENDDO        ! idum

    ENDIF  ! doDummy
  ENDDO    ! iside

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_FloFlimCorrCornEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimCorrCornEdgeCells.F90,v $
! Revision 1.4  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:13  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







