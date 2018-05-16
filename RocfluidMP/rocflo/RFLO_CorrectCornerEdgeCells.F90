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
! Purpose: correct values in the corner and edge cells at noslip (viscous)
!          walls, symmetry boundaries, and non-burning injection boundaries.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch
!        bcType = type of boundary condition.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_CorrectCornerEdgeCells.F90,v 1.5 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CorrectCornerEdgeCells( region,patch,bcType )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, MixtureProperties, MixtPerf_R_M, MixtPerf_G_CpR, &
        MixtPerf_D_PRT, MixtPerf_Eo_DGPUVW, RFLO_GetDimensPhys, &
        RFLO_GetDimensDummy, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER        :: bcType
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, n1, n2, iside

! ... local variables
  LOGICAL :: doDummy(4), moveGrid

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ibegp, iendp, jbegp, jendp, kbegp, kendp, idir, jdir, kdir
  INTEGER :: ibeg(4), iend(4), jbeg(4), jend(4), kbeg(4), kend(4), nsides
  INTEGER :: iLev, lbound, iCOff, ijCOff, nOff, distrib, bcOpt
  INTEGER :: gasModel, ijkC, ijkD, i2d, iNOff, ijNOff, ijkN(4)
  INTEGER :: inode, jnode, knode

  REAL(RFREAL)          :: dt, mRate, sgn(3), dxn(4), dyn(4), dzn(4)
  REAL(RFREAL)          :: uSurf, vSurf, wSurf
  REAL(RFREAL), POINTER :: cv(:,:), vals(:,:), xyz(:,:), xyzOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CorrectCornerEdgeCells', &
                         'RFLO_CorrectCornerEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibegp,iendp,jbegp,jendp,kbegp,kendp )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib
  gasModel = region%mixtInput%gasModel
  moveGrid  = region%mixtInput%moveGrid
  dt        = region%global%dtMin

  cv   => region%levels(iLev)%mixt%cv
  vals => patch%mixt%vals
  xyz  => region%levels(iLev)%grid%xyz
  IF (moveGrid) xyzOld => region%levels(iLev)%gridOld%xyz

! to take the correct surface coordinates

  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

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

! settings to mirror velocity components (symmetry) ---------------------------

  IF      (lbound==1 .OR. lbound==2) THEN
    sgn(1) = -1._RFREAL
    sgn(2) = +1._RFREAL
    sgn(3) = +1._RFREAL
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    sgn(1) = +1._RFREAL
    sgn(2) = -1._RFREAL
    sgn(3) = +1._RFREAL
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    sgn(1) = +1._RFREAL
    sgn(2) = +1._RFREAL
    sgn(3) = -1._RFREAL
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

              IF (bcType>=BC_NOSLIPWALL .AND. &
                  bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
                IF (moveGrid) THEN
                  ijkN(1) = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
                  IF      (lbound==1 .OR. lbound==2) THEN
                    ijkN(2) = IndIJK(i+inode,j+jnode+1,k+knode  ,iNOff,ijNOff)
                    ijkN(3) = IndIJK(i+inode,j+jnode+1,k+knode+1,iNOff,ijNOff)
                    ijkN(4) = IndIJK(i+inode,j+jnode  ,k+knode+1,iNOff,ijNOff)
                  ELSE IF (lbound==3 .OR. lbound==4) THEN
                    ijkN(2) = IndIJK(i+inode  ,j+jnode,k+knode+1,iNOff,ijNOff)
                    ijkN(3) = IndIJK(i+inode+1,j+jnode,k+knode+1,iNOff,ijNOff)
                    ijkN(4) = IndIJK(i+inode+1,j+jnode,k+knode  ,iNOff,ijNOff)
                  ELSE IF (lbound==5 .OR. lbound==6) THEN
                    ijkN(2) = IndIJK(i+inode+1,j+jnode  ,k+knode,iNOff,ijNOff)
                    ijkN(3) = IndIJK(i+inode+1,j+jnode+1,k+knode,iNOff,ijNOff)
                    ijkN(4) = IndIJK(i+inode  ,j+jnode+1,k+knode,iNOff,ijNOff)
                  ENDIF
                  dxn(1) = xyz(XCOORD,ijkN(1)) - xyzOld(XCOORD,ijkN(1))
                  dyn(1) = xyz(YCOORD,ijkN(1)) - xyzOld(YCOORD,ijkN(1))
                  dzn(1) = xyz(ZCOORD,ijkN(1)) - xyzOld(ZCOORD,ijkN(1))
                  dxn(2) = xyz(XCOORD,ijkN(2)) - xyzOld(XCOORD,ijkN(2))
                  dyn(2) = xyz(YCOORD,ijkN(2)) - xyzOld(YCOORD,ijkN(2))
                  dzn(2) = xyz(ZCOORD,ijkN(2)) - xyzOld(ZCOORD,ijkN(2))
                  dxn(3) = xyz(XCOORD,ijkN(3)) - xyzOld(XCOORD,ijkN(3))
                  dyn(3) = xyz(YCOORD,ijkN(3)) - xyzOld(YCOORD,ijkN(3))
                  dzn(3) = xyz(ZCOORD,ijkN(3)) - xyzOld(ZCOORD,ijkN(3))
                  dxn(4) = xyz(XCOORD,ijkN(4)) - xyzOld(XCOORD,ijkN(4))
                  dyn(4) = xyz(YCOORD,ijkN(4)) - xyzOld(YCOORD,ijkN(4))
                  dzn(4) = xyz(ZCOORD,ijkN(4)) - xyzOld(ZCOORD,ijkN(4))
                  uSurf  = 0.25_RFREAL*(dxn(1)+dxn(2)+dxn(3)+dxn(4))/dt
                  vSurf  = 0.25_RFREAL*(dyn(1)+dyn(2)+dyn(3)+dyn(4))/dt
                  wSurf  = 0.25_RFREAL*(dzn(1)+dzn(2)+dzn(3)+dzn(4))/dt
                ELSE
                  uSurf = 0._RFREAL
                  vSurf = 0._RFREAL
                  wSurf = 0._RFREAL
                ENDIF
                cv(CV_MIXT_DENS,ijkD) = cv(CV_MIXT_DENS,ijkC)
                cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*uSurf - cv(CV_MIXT_XMOM,ijkC)
                cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*vSurf - cv(CV_MIXT_YMOM,ijkC)
                cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*wSurf - cv(CV_MIXT_ZMOM,ijkC)
                cv(CV_MIXT_ENER,ijkD) = cv(CV_MIXT_ENER,ijkC)

! ----------- symmetry boundary

              ELSE IF (bcType>=BC_SYMMETRY .AND. &
                       bcType<=BC_SYMMETRY+BC_RANGE) THEN
                cv(CV_MIXT_DENS,ijkD) =        cv(CV_MIXT_DENS,ijkC)
                cv(CV_MIXT_XMOM,ijkD) = sgn(1)*cv(CV_MIXT_XMOM,ijkC)
                cv(CV_MIXT_YMOM,ijkD) = sgn(2)*cv(CV_MIXT_YMOM,ijkC)
                cv(CV_MIXT_ZMOM,ijkD) = sgn(3)*cv(CV_MIXT_ZMOM,ijkC)
                cv(CV_MIXT_ENER,ijkD) =        cv(CV_MIXT_ENER,ijkC)

! ----------- non-burning injection boundary (like noslip wall)

              ELSE IF (bcType>=BC_INJECTION .AND. &
                       bcType<=BC_INJECTION+BC_RANGE) THEN
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
                  cv(CV_MIXT_DENS,ijkD) =  cv(CV_MIXT_DENS,ijkC)
                  cv(CV_MIXT_XMOM,ijkD) = -cv(CV_MIXT_XMOM,ijkC)
                  cv(CV_MIXT_YMOM,ijkD) = -cv(CV_MIXT_YMOM,ijkC)
                  cv(CV_MIXT_ZMOM,ijkD) = -cv(CV_MIXT_ZMOM,ijkC)
                  cv(CV_MIXT_ENER,ijkD) =  cv(CV_MIXT_ENER,ijkC)
                ENDIF
              ENDIF

              IF (gasModel == GAS_MODEL_TCPERF) THEN
                CALL MixtureProperties( region,ijkD,ijkD,.false. )
              ELSE
                CALL MixtureProperties( region,ijkD,ijkD,.true.  )
              ENDIF

            ENDDO  ! i
          ENDDO    ! j
        ENDDO      ! k
      ENDDO        ! idum

    ENDIF  ! doDummy
  ENDDO    ! iside

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CorrectCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CorrectCornerEdgeCells.F90,v $
! Revision 1.5  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:39:33  mparmar
! Renamed patch variables
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.1  2003/05/20 20:46:57  jblazek
! Values in edge & corner cells now corrected at noslip and symmetry walls.
!
!******************************************************************************







