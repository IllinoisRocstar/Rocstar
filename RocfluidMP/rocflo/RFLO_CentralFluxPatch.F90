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
! Purpose: compute central convective fluxes through a patch
!          by using an average of variables.
!
! Description: none.
!
! Input: region = data of current region
!        patch  = current patch.
!
! Output: region%levels%mixt%rhs = convective fluxes added to the residual.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_CentralFluxPatch.F90,v 1.6 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CentralFluxPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, BcondInjectionPerf
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, n1, n2

! ... local variables
  INTEGER :: iLev, lbound, bcType, distrib, flowModel, gasModel, indCp, indMol
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, i2d, nOff
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkCD, ijkCB0, ijkCB1, ijkNB
  INTEGER :: inode, jnode, knode, indSvel, acId0, acId1, aeroCoeff, j2d

  REAL(RFREAL)          :: sgn, rhoa, rhoua, rhova, rhowa, rhoea, rhoVrel(3)
  REAL(RFREAL)          :: pa, mRate, tBurn, rgas, dS, sf(3)
  REAL(RFREAL)          :: uinj, vinj, winj, vcont, sv, ac0, ac1
  REAL(RFREAL)          :: pRef, rRef, vRef, rCpRef
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), rhs(:,:)
  REAL(RFREAL), POINTER :: avgCo(:,:), sFace(:,:), sVel(:), vals(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CentralFluxPatch',&
  'RFLO_CentralFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  bcType    = patch%bcType
  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib
  flowModel = region%mixtInput%flowModel
  gasModel  = region%mixtInput%gasModel
  indCp     = region%levels(iLev)%mixt%indCp
  indMol    = region%levels(iLev)%mixt%indMol
  indSvel   = region%levels(iLev)%grid%indSvel

  cv   => region%levels(iLev)%mixt%cv
  dv   => region%levels(iLev)%mixt%dv
  gv   => region%levels(iLev)%mixt%gv
  rhs  => region%levels(iLev)%mixt%rhs
  vals => patch%mixt%vals

  aeroCoeff = region%global%aeroCoeffs
  pRef      = region%global%refPressure
  rRef      = region%global%refDensity
  vRef      = region%global%refVelocity

  rCpRef = 2.0_RFREAL/(rRef*vRef*vRef)

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  acId0 = 2
  acId1 = 1
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
    acId0 = 1
    acId1 = 2
  ENDIF

! get the appropriate face vector and grid speed

  IF (lbound==1 .OR. lbound==2) THEN
    avgCo => region%levels(iLev)%grid%c2fCoI
    sFace => region%levels(iLev)%grid%si
    sVel  => region%levels(iLev)%grid%siVel
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    avgCo => region%levels(iLev)%grid%c2fCoJ
    sFace => region%levels(iLev)%grid%sj
    sVel  => region%levels(iLev)%grid%sjVel
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    avgCo => region%levels(iLev)%grid%c2fCoK
    sFace => region%levels(iLev)%grid%sk
    sVel  => region%levels(iLev)%grid%skVel
  ENDIF

! stationary grid -------------------------------------------------------------
! slip wall

  IF (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) THEN

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCB1 = IndIJK(i+idir ,j+jdir ,k+kdir ,iCOff,ijCOff)  ! interior
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          pa     = 0.5_RFREAL*(3._RFREAL*dv(DV_MIXT_PRES,ijkCB0)- &
                                         dv(DV_MIXT_PRES,ijkCB1))
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
          sv     = sgn*sVel(ijkNB*indSvel)

          rhs(CV_MIXT_XMOM,ijkCB0) = rhs(CV_MIXT_XMOM,ijkCB0) + pa*sf(1)
          rhs(CV_MIXT_YMOM,ijkCB0) = rhs(CV_MIXT_YMOM,ijkCB0) + pa*sf(2)
          rhs(CV_MIXT_ZMOM,ijkCB0) = rhs(CV_MIXT_ZMOM,ijkCB0) + pa*sf(3)
          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) + pa*sv

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
          j2d  = aeroCoeff * IndIJ(n1,n2,nOff)
          patch%cp(j2d) = rCpRef*(pa - pRef)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! noslip wall

  ELSE IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          pa     = dv(DV_MIXT_PRES,ijkCB0)
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
          sv     = sgn*sVel(ijkNB*indSvel)

          rhs(CV_MIXT_XMOM,ijkCB0) = rhs(CV_MIXT_XMOM,ijkCB0) + pa*sf(1)
          rhs(CV_MIXT_YMOM,ijkCB0) = rhs(CV_MIXT_YMOM,ijkCB0) + pa*sf(2)
          rhs(CV_MIXT_ZMOM,ijkCB0) = rhs(CV_MIXT_ZMOM,ijkCB0) + pa*sf(3)
          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) + pa*sv

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
          j2d  = aeroCoeff * IndIJ(n1,n2,nOff)
          patch%cp(j2d) = rCpRef*(pa - pRef)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! symmetry boundary

  ELSE IF (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE) THEN

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          pa     = 0.5_RFREAL*(dv(DV_MIXT_PRES,ijkCB0)+dv(DV_MIXT_PRES,ijkCD))
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
          sv     = sgn*sVel(ijkNB*indSvel)

          rhs(CV_MIXT_XMOM,ijkCB0) = rhs(CV_MIXT_XMOM,ijkCB0) + pa*sf(1)
          rhs(CV_MIXT_YMOM,ijkCB0) = rhs(CV_MIXT_YMOM,ijkCB0) + pa*sf(2)
          rhs(CV_MIXT_ZMOM,ijkCB0) = rhs(CV_MIXT_ZMOM,ijkCB0) + pa*sf(3)
          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) + pa*sv

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
          j2d  = aeroCoeff * IndIJ(n1,n2,nOff)
          patch%cp(j2d) = rCpRef*(pa - pRef)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! injection boundary (as wall if mass flow rate <= 0)

  ELSE IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCB1 = IndIJK(i+idir ,j+jdir ,k+kdir ,iCOff,ijCOff)  ! interior
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
          sv     = sgn*sVel(ijkNB*indSvel)

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
          i2d        = distrib * IndIJ(n1,n2,nOff)
          mRate      = vals(BCDAT_INJECT_MFRATE,i2d)
          tBurn      = vals(BCDAT_INJECT_TEMP  ,i2d)
          rhoVrel(1) = vals(BCDAT_INJECT_RFVFU ,i2d)
          rhoVrel(2) = vals(BCDAT_INJECT_RFVFV ,i2d)
          rhoVrel(3) = vals(BCDAT_INJECT_RFVFW ,i2d)

          IF (mRate > 0._RFREAL) THEN        ! surface burning
            dS = SQRT(sf(1)*sf(1)+sf(2)*sf(2)+sf(3)*sf(3))
            IF (gasModel == GAS_MODEL_TCPERF) THEN
              CALL BcondInjectionPerf( distrib,mRate,tBurn,rhoVrel, &
                                       sf(1)/dS,sf(2)/dS,sf(3)/dS, &
                                       gv(GV_MIXT_CP  ,ijkCB0*indCp ), &
                                       gv(GV_MIXT_MOL ,ijkCB0*indMol), &
                                       dv(DV_MIXT_PRES,ijkCB0       ), &
                                       rhoa,rhoua,rhova,rhowa,rhoea,pa, &
                                       uinj,vinj,winj )
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,__LINE__ )
            ENDIF
            vcont = uinj*sf(1) + vinj*sf(2) + winj*sf(3)

          ELSE                               ! not burning - slip/noslip wall
            rhoa     = 0._RFREAL
            rhoua    = 0._RFREAL
            rhova    = 0._RFREAL
            rhowa    = 0._RFREAL
            rhoea    = 0._RFREAL
            vcont    = 0._RFREAL
            IF (flowModel == FLOW_EULER) THEN
              pa = 0.5_RFREAL*(3._RFREAL*dv(DV_MIXT_PRES,ijkCB0)- &
                                         dv(DV_MIXT_PRES,ijkCB1))
            ELSE
              pa = dv(DV_MIXT_PRES,ijkCB0)
            ENDIF
          ENDIF

          rhs(CV_MIXT_DENS,ijkCB0) = rhs(CV_MIXT_DENS,ijkCB0) + vcont*rhoa
          rhs(CV_MIXT_XMOM,ijkCB0) = rhs(CV_MIXT_XMOM,ijkCB0) + vcont*rhoua + &
                                                                pa*sf(1)
          rhs(CV_MIXT_YMOM,ijkCB0) = rhs(CV_MIXT_YMOM,ijkCB0) + vcont*rhova + &
                                                                pa*sf(2)
          rhs(CV_MIXT_ZMOM,ijkCB0) = rhs(CV_MIXT_ZMOM,ijkCB0) + vcont*rhowa + &
                                                                pa*sf(3)
          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) + &
                                     vcont*(rhoea+pa) + pa*sv

          j2d  = aeroCoeff * IndIJ(n1,n2,nOff)
          patch%cp(j2d) = rCpRef*(pa - pRef)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! non-conforming region interface

  ELSE IF (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) THEN

  ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN

! everything else

  ELSE

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          ac0    = avgCo(acId0,ijkNB)
          ac1    = avgCo(acId1,ijkNB)

          rhoa   = ac0*cv(CV_MIXT_DENS,ijkCB0)+ac1*cv(CV_MIXT_DENS,ijkCD)
          rhoua  = ac0*cv(CV_MIXT_XMOM,ijkCB0)+ac1*cv(CV_MIXT_XMOM,ijkCD)
          rhova  = ac0*cv(CV_MIXT_YMOM,ijkCB0)+ac1*cv(CV_MIXT_YMOM,ijkCD)
          rhowa  = ac0*cv(CV_MIXT_ZMOM,ijkCB0)+ac1*cv(CV_MIXT_ZMOM,ijkCD)
          rhoea  = ac0*cv(CV_MIXT_ENER,ijkCB0)+ac1*cv(CV_MIXT_ENER,ijkCD)
          pa     = ac0*dv(DV_MIXT_PRES,ijkCB0)+ac1*dv(DV_MIXT_PRES,ijkCD)
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
          sv     = sgn*sVel(ijkNB*indSvel)
          vcont  = (rhoua*sf(1)+rhova*sf(2)+rhowa*sf(3))/rhoa - sv

          rhs(CV_MIXT_DENS,ijkCB0) = rhs(CV_MIXT_DENS,ijkCB0) + vcont*rhoa
          rhs(CV_MIXT_XMOM,ijkCB0) = rhs(CV_MIXT_XMOM,ijkCB0) + &
                                     vcont*rhoua + pa*sf(1)
          rhs(CV_MIXT_YMOM,ijkCB0) = rhs(CV_MIXT_YMOM,ijkCB0) + &
                                     vcont*rhova + pa*sf(2)
          rhs(CV_MIXT_ZMOM,ijkCB0) = rhs(CV_MIXT_ZMOM,ijkCB0) + &
                                     vcont*rhowa + pa*sf(3)
          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) + &
                                     vcont*(rhoea+pa) + sv*pa

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
          j2d  = aeroCoeff * IndIJ(n1,n2,nOff)
          patch%cp(j2d) = rCpRef*(pa - pRef)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  ENDIF         ! bcType

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CentralFluxPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CentralFluxPatch.F90,v $
! Revision 1.6  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:39:28  mparmar
! Renamed patch variables
!
! Revision 1.3  2006/03/13 03:41:29  wasistho
! defined aero coeffs
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.25  2004/08/02 21:56:27  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.24  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.20  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.19  2003/05/29 17:28:43  jblazek
! Implemented Roe scheme.
!
! Revision 1.18  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.17  2003/04/09 15:07:01  jferry
! removed temporary rocsmoke storage structures
!
! Revision 1.16  2003/02/11 22:53:19  jferry
! Initial import of Rocsmoke
!
! Revision 1.15  2002/12/06 23:06:08  jblazek
! Grid speeds added only at free boundaries.
!
! Revision 1.14  2002/10/03 22:00:30  jblazek
! Removed grid speeds from fluxes at burning boundaries.
!
! Revision 1.13  2002/10/03 21:25:39  jblazek
! Changed init. of burnig boundaries for GenX.
!
! Revision 1.12  2002/09/27 00:29:04  jblazek
! Changed injection BC for moving boundaries.
!
! Revision 1.11  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.10  2002/08/30 18:25:55  jblazek
! Forgot to multiply grid speeds by face area ...
!
! Revision 1.9  2002/08/29 23:25:54  jblazek
! Added support for moving grids.
!
! Revision 1.8  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.7  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.6  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.5  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.4  2002/02/06 00:15:39  jblazek
! Improved injection BC. Added pointers to gradients.
!
! Revision 1.3  2002/02/04 15:30:25  jblazek
! Added injection boundary condition.
!
! Revision 1.2  2002/02/01 22:17:38  jblazek
! Change addressing of face vectors at block boundaries.
!
! Revision 1.1  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
!******************************************************************************







