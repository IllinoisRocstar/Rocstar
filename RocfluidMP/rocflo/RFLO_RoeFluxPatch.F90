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
! Purpose: compute convective fluxes through a patch
!          by using Roe`s upwind scheme.
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
! $Id: RFLO_RoeFluxPatch.F90,v 1.6 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_RoeFluxPatch( region,patch )

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
  INTEGER :: iCOff, ijCOff, ijkCD, ijkCD1, ijkCD2, ijkCB0, ijkCB1, ijkNB
  INTEGER :: iNOff, ijNOff, inode, jnode, knode, indSvel, spaceOrder
  INTEGER :: aeroCoeff, j2d

  REAL(RFREAL) :: sgn, rhoa, rhoua, rhova, rhowa, rhoea, rhoVrel(3), pa, &
                  mRate, tBurn, rgas, dS, sv, limFac, limFac3, rVolRef, &
                  eps2n, uinj, vinj, winj, vcont, rhl, rhr, qsl, qsr, vola, &
                  gam, ggm1, rl, rr, ul, ur, vl, vr, wl, wr, pl, pr, &
                  hl, hr, qsrl, qsrr, eps2(3), dVar(5), dVarm(5), dVarp(5), &
                  deltl(5), deltr(5), sf(3)
  REAL(RFREAL) :: pRef, rRef, vRef, rCpRef
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), rhs(:,:), vol(:)
  REAL(RFREAL), POINTER :: sFace(:,:), sVel(:), vals(:,:)

! ... functions
  REAL(RFREAL) :: MUSCL3, af, bf, eps

!******************************************************************************
! limiter function

  MUSCL3(af,bf,eps) = (bf*(2._RFREAL*af*af+eps)+af*(bf*bf+2._RFREAL*eps))/ &
                      (2._RFREAL*af*af+2._RFREAL*bf*bf-af*bf+ &
                       3._RFREAL*eps+1.E-30_RFREAL)

  CALL RegisterFunction( region%global,'RFLO_RoeFluxPatch',&
  'RFLO_RoeFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  bcType     = patch%bcType
  nOff       = ABS(patch%l1end-patch%l1beg) + 1
  distrib    = patch%mixt%distrib
  flowModel  = region%mixtInput%flowModel
  gasModel  = region%mixtInput%gasModel
  spaceOrder = region%mixtInput%spaceOrder
  indCp      = region%levels(iLev)%mixt%indCp
  indMol     = region%levels(iLev)%mixt%indMol
  indSvel    = region%levels(iLev)%grid%indSvel
  limFac     = region%mixtInput%limFac

  cv   => region%levels(iLev)%mixt%cv
  dv   => region%levels(iLev)%mixt%dv
  gv   => region%levels(iLev)%mixt%gv
  vol  => region%levels(iLev)%grid%vol
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
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector and grid speed

  IF (lbound==1 .OR. lbound==2) THEN
    sFace => region%levels(iLev)%grid%si
    sVel  => region%levels(iLev)%grid%siVel
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    sFace => region%levels(iLev)%grid%sj
    sVel  => region%levels(iLev)%grid%sjVel
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    sFace => region%levels(iLev)%grid%sk
    sVel  => region%levels(iLev)%grid%skVel
  ENDIF

! normalise epsilon^2 for all limited variables (rho, u, v, w, p) -------------

  limFac3 = limFac*limFac*limFac
  rVolRef = 1._RFREAL/region%global%limVolRef**1.5_RFREAL
  eps2(1) = limFac3*region%global%limRef(1)*region%global%limRef(1)*rVolRef
  eps2(2) = limFac3*region%global%limRef(2)*region%global%limRef(2)*rVolRef
  eps2(3) = limFac3*region%global%limRef(3)*region%global%limRef(3)*rVolRef

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

    IF (spaceOrder == DISCR_ORDER_1) THEN
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
            ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
            ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            sf(1)  = sgn*sFace(XCOORD,ijkNB)
            sf(2)  = sgn*sFace(YCOORD,ijkNB)
            sf(3)  = sgn*sFace(ZCOORD,ijkNB)
            sv     = sgn*sVel(ijkNB*indSvel)
            rhl    = dv(DV_MIXT_PRES,ijkCD) + cv(CV_MIXT_ENER,ijkCD)
            qsl    = dv(DV_MIXT_UVEL,ijkCD)*sf(1) + &
                     dv(DV_MIXT_VVEL,ijkCD)*sf(2) + &
                     dv(DV_MIXT_WVEL,ijkCD)*sf(3) - sv
            rhr    = dv(DV_MIXT_PRES,ijkCB0) + cv(CV_MIXT_ENER,ijkCB0)
            qsr    = dv(DV_MIXT_UVEL,ijkCB0)*sf(1) + &
                     dv(DV_MIXT_VVEL,ijkCB0)*sf(2) + &
                     dv(DV_MIXT_WVEL,ijkCB0)*sf(3) - sv
            pa     = 0.5_RFREAL*(dv(DV_MIXT_PRES,ijkCB0)+dv(DV_MIXT_PRES,ijkCD))
            rhs(CV_MIXT_DENS,ijkCB0) = rhs(CV_MIXT_DENS,ijkCB0) + &
                                       0.5_RFREAL*(qsl*cv(CV_MIXT_DENS,ijkCD)+ &
                                                   qsr*cv(CV_MIXT_DENS,ijkCB0))
            rhs(CV_MIXT_XMOM,ijkCB0) = rhs(CV_MIXT_XMOM,ijkCB0) + pa*sf(1) + &
                                       0.5_RFREAL*(qsl*cv(CV_MIXT_XMOM,ijkCD)+ &
                                                   qsr*cv(CV_MIXT_XMOM,ijkCB0))
            rhs(CV_MIXT_YMOM,ijkCB0) = rhs(CV_MIXT_YMOM,ijkCB0) + pa*sf(2) + &
                                       0.5_RFREAL*(qsl*cv(CV_MIXT_YMOM,ijkCD)+ &
                                                   qsr*cv(CV_MIXT_YMOM,ijkCB0))
            rhs(CV_MIXT_ZMOM,ijkCB0) = rhs(CV_MIXT_ZMOM,ijkCB0) + pa*sf(3) + &
                                       0.5_RFREAL*(qsl*cv(CV_MIXT_ZMOM,ijkCD)+ &
                                                   qsr*cv(CV_MIXT_ZMOM,ijkCB0))
            rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) + &
                                       0.5_RFREAL*(qsl*rhl+qsr*rhr) + sv*pa

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

          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k
    ELSE IF (spaceOrder == DISCR_ORDER_2) THEN
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkCB0 = IndIJK(i       ,j        ,k      ,iCOff,ijCOff)  ! boundary
            ijkCB1 = IndIJK(i+  idir,j+  jdir,k+  kdir,iCOff,ijCOff)  ! interior
            ijkCD1 = IndIJK(i-  idir,j-  jdir,k-  kdir,iCOff,ijCOff)  ! dummy 1
            ijkCD2 = IndIJK(i-2*idir,j-2*jdir,k-2*kdir,iCOff,ijCOff)  ! dummy 2
            ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            sf(1)  = sgn*sFace(XCOORD,ijkNB)
            sf(2)  = sgn*sFace(YCOORD,ijkNB)
            sf(3)  = sgn*sFace(ZCOORD,ijkNB)
            sv     = sgn*sVel(ijkNB*indSvel)
            vola   = (0.5_RFREAL*(vol(ijkCB0)+vol(ijkCD1)))**1.5_RFREAL

            dVar(1)  = cv(CV_MIXT_DENS,ijkCB0) - cv(CV_MIXT_DENS,ijkCD1)
            dVar(2)  = dv(DV_MIXT_UVEL,ijkCB0) - dv(DV_MIXT_UVEL,ijkCD1)
            dVar(3)  = dv(DV_MIXT_VVEL,ijkCB0) - dv(DV_MIXT_VVEL,ijkCD1)
            dVar(4)  = dv(DV_MIXT_WVEL,ijkCB0) - dv(DV_MIXT_WVEL,ijkCD1)
            dVar(5)  = dv(DV_MIXT_PRES,ijkCB0) - dv(DV_MIXT_PRES,ijkCD1)

            dVarm(1) = cv(CV_MIXT_DENS,ijkCD1) - cv(CV_MIXT_DENS,ijkCD2)
            dVarm(2) = dv(DV_MIXT_UVEL,ijkCD1) - dv(DV_MIXT_UVEL,ijkCD2)
            dVarm(3) = dv(DV_MIXT_VVEL,ijkCD1) - dv(DV_MIXT_VVEL,ijkCD2)
            dVarm(4) = dv(DV_MIXT_WVEL,ijkCD1) - dv(DV_MIXT_WVEL,ijkCD2)
            dVarm(5) = dv(DV_MIXT_PRES,ijkCD1) - dv(DV_MIXT_PRES,ijkCD2)

            dVarp(1) = cv(CV_MIXT_DENS,ijkCB1) - cv(CV_MIXT_DENS,ijkCB0)
            dVarp(2) = dv(DV_MIXT_UVEL,ijkCB1) - dv(DV_MIXT_UVEL,ijkCB0)
            dVarp(3) = dv(DV_MIXT_VVEL,ijkCB1) - dv(DV_MIXT_VVEL,ijkCB0)
            dVarp(4) = dv(DV_MIXT_WVEL,ijkCB1) - dv(DV_MIXT_WVEL,ijkCB0)
            dVarp(5) = dv(DV_MIXT_PRES,ijkCB1) - dv(DV_MIXT_PRES,ijkCB0)

            eps2n    = eps2(1)*vola
            deltl(1) = 0.5_RFREAL*MUSCL3(dVar(1) ,dVarm(1),eps2n)
            deltr(1) = 0.5_RFREAL*MUSCL3(dVarp(1),dVar(1) ,eps2n)
            eps2n    = eps2(2)*vola
            deltl(2) = 0.5_RFREAL*MUSCL3(dVar(2) ,dVarm(2),eps2n)
            deltr(2) = 0.5_RFREAL*MUSCL3(dVarp(2),dVar(2) ,eps2n)
            deltl(3) = 0.5_RFREAL*MUSCL3(dVar(3) ,dVarm(3),eps2n)
            deltr(3) = 0.5_RFREAL*MUSCL3(dVarp(3),dVar(3) ,eps2n)
            deltl(4) = 0.5_RFREAL*MUSCL3(dVar(4) ,dVarm(4),eps2n)
            deltr(4) = 0.5_RFREAL*MUSCL3(dVarp(4),dVar(4) ,eps2n)
            eps2n    = eps2(3)*vola
            deltl(5) = 0.5_RFREAL*MUSCL3(dVar(5) ,dVarm(5),eps2n)
            deltr(5) = 0.5_RFREAL*MUSCL3(dVarp(5),dVar(5) ,eps2n)

            rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkCD1*indMol)
            gam   = gv(GV_MIXT_CP,ijkCD1*indCp)/(gv(GV_MIXT_CP,ijkCD1*indCp)-rgas)
            ggm1  = gam/(gam-1._RFREAL)
            rl    = cv(CV_MIXT_DENS,ijkCD1) + deltl(1)
            ul    = dv(DV_MIXT_UVEL,ijkCD1) + deltl(2)
            vl    = dv(DV_MIXT_VVEL,ijkCD1) + deltl(3)
            wl    = dv(DV_MIXT_WVEL,ijkCD1) + deltl(4)
            pl    = dv(DV_MIXT_PRES,ijkCD1) + deltl(5)
            hl    = ggm1*pl/rl + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)
            qsrl  = (ul*sf(1)+vl*sf(2)+ wl*sf(3)-sv)*rl

            rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkCB0*indMol)
            gam   = gv(GV_MIXT_CP,ijkCB0*indCp)/(gv(GV_MIXT_CP,ijkCB0*indCp)-rgas)
            ggm1  = gam/(gam-1._RFREAL)
            rr    = cv(CV_MIXT_DENS,ijkCB0) - deltr(1)
            ur    = dv(DV_MIXT_UVEL,ijkCB0) - deltr(2)
            vr    = dv(DV_MIXT_VVEL,ijkCB0) - deltr(3)
            wr    = dv(DV_MIXT_WVEL,ijkCB0) - deltr(4)
            pr    = dv(DV_MIXT_PRES,ijkCB0) - deltr(5)
            hr    = ggm1*pr/rr + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)
            qsrr  = (ur*sf(1)+vr*sf(2)+wr*sf(3)-sv)*rr
            pa    = 0.5_RFREAL*(pl+pr)

            rhs(CV_MIXT_DENS,ijkCB0) = rhs(CV_MIXT_DENS,ijkCB0) + &
                                       0.5_RFREAL*(qsrl   +qsrr   )
            rhs(CV_MIXT_XMOM,ijkCB0) = rhs(CV_MIXT_XMOM,ijkCB0) + pa*sf(1) + &
                                       0.5_RFREAL*(qsrl*ul+qsrr*ur)
            rhs(CV_MIXT_YMOM,ijkCB0) = rhs(CV_MIXT_YMOM,ijkCB0) + pa*sf(2) + &
                                       0.5_RFREAL*(qsrl*vl+qsrr*vr)
            rhs(CV_MIXT_ZMOM,ijkCB0) = rhs(CV_MIXT_ZMOM,ijkCB0) + pa*sf(3) + &
                                       0.5_RFREAL*(qsrl*wl+qsrr*wr)
            rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) + &
                                       0.5_RFREAL*(qsrl*hl+qsrr*hr) + sv*pa

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

          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k
    ENDIF       ! spaceOrder

  ENDIF         ! bcType

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_RoeFluxPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_RoeFluxPatch.F90,v $
! Revision 1.6  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:39:42  mparmar
! Renamed patch variables
!
! Revision 1.3  2006/03/13 03:41:49  wasistho
! defined aero coeffs
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.1  2003/05/29 17:28:43  jblazek
! Implemented Roe scheme.
!
!******************************************************************************







