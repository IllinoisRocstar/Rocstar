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
! Purpose: compute numerical dissipation based on 2nd-order
!          Roe`s upwind scheme.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%diss = dissipative fluxes.
!
! Notes: uses MUSCL scheme with kappa=1/3.
!
!******************************************************************************
!
! $Id: RFLO_RoeDissipSecond.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_RoeDissipSecond( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkCm1, ijkCm2, ijkCp1, ijkN
  INTEGER :: iLev, indCp, indMol, indSvel

  REAL(RFREAL) :: beta5, dS, nx, ny, nz, sVel, rgas, gaml, gamr, ggm1, gam1, &
                  rl, ul, vl, wl, pl, hl, rr, ur, vr, wr, pr, hr, rav, dd, &
                  dd1, uav, vav, wav, hav, q2a, c2a, cav, uvw, du, eabs1, &
                  eabs2, eabs5, h1, h2, h3, h4, h5, epsEntr, delta, limFac, &
                  limFac3, rvolRef, eps2(3), eps2n, vola, dVar(5), dVarm(5), &
                  dVarp(5), deltl(5), deltr(5), fd(5)
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), diss(:,:), &
                           vol(:), si(:,:), sj(:,:), sk(:,:), &
                           siVel(:), sjVel(:), skVel(:)

! ... functions
  REAL(RFREAL) :: MUSCL3, af, bf, eps

!******************************************************************************
! limiter function

  MUSCL3(af,bf,eps) = (bf*(2._RFREAL*af*af+eps)+af*(bf*bf+2._RFREAL*eps))/ &
                      (2._RFREAL*af*af+2._RFREAL*bf*bf-af*bf+ &
                       3._RFREAL*eps+1.E-30_RFREAL)

  CALL RegisterFunction( region%global,'RFLO_RoeDissipSecond',&
  'RFLO_RoeDissipSecond.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv    => region%levels(iLev)%mixt%cv
  dv    => region%levels(iLev)%mixt%dv
  gv    => region%levels(iLev)%mixt%gv
  vol   => region%levels(iLev)%grid%vol
  si    => region%levels(iLev)%grid%si
  sj    => region%levels(iLev)%grid%sj
  sk    => region%levels(iLev)%grid%sk
  siVel => region%levels(iLev)%grid%siVel
  sjVel => region%levels(iLev)%grid%sjVel
  skVel => region%levels(iLev)%grid%skVel
  diss  => region%levels(iLev)%mixt%diss

! get coefficients ------------------------------------------------------------

  beta5   = 0.5_RFREAL*region%mixtInput%betrk(region%irkStep)
  epsEntr = region%mixtInput%epsEntr
  limFac  = region%mixtInput%limFac
  indSvel = region%levels(iLev)%grid%indSvel
  indCp   = region%levels(iLev)%mixt%indCp
  indMol  = region%levels(iLev)%mixt%indMol

! normalise epsilon^2 for all limited variables (rho, u, v, w, p) -------------

  limFac3 = limFac*limFac*limFac
  rVolRef = 1._RFREAL/region%global%limVolRef**1.5_RFREAL
  eps2(1) = limFac3*region%global%limRef(1)*region%global%limRef(1)*rVolRef
  eps2(2) = limFac3*region%global%limRef(2)*region%global%limRef(2)*rVolRef
  eps2(3) = limFac3*region%global%limRef(3)*region%global%limRef(3)*rVolRef

! dissipation in i-direction --------------------------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend+1
        ijkC0  = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkCm1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkCm2 = IndIJK(i-2,j,k,iCOff,ijCOff)
        ijkCp1 = IndIJK(i+1,j,k,iCOff,ijCOff)
        ijkN   = IndIJK(i  ,j,k,iNOff,ijNOff)
        dS     = SQRT(si(XCOORD,ijkN)*si(XCOORD,ijkN)+ &
                      si(YCOORD,ijkN)*si(YCOORD,ijkN)+ &
                      si(ZCOORD,ijkN)*si(ZCOORD,ijkN))
        nx     = si(XCOORD,ijkN)/dS
        ny     = si(YCOORD,ijkN)/dS
        nz     = si(ZCOORD,ijkN)/dS
        sVel   = siVel(ijkN*indSvel)/dS
        vola   = (0.5_RFREAL*(vol(ijkC0)+vol(ijkCm1)))**1.5_RFREAL
        dS     = dS*beta5

! ----- limited extrapolations

        dVar(1)  = cv(CV_MIXT_DENS,ijkC0 ) - cv(CV_MIXT_DENS,ijkCm1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0 ) - dv(DV_MIXT_UVEL,ijkCm1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0 ) - dv(DV_MIXT_VVEL,ijkCm1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0 ) - dv(DV_MIXT_WVEL,ijkCm1)
        dVar(5)  = dv(DV_MIXT_PRES,ijkC0 ) - dv(DV_MIXT_PRES,ijkCm1)

        dVarm(1) = cv(CV_MIXT_DENS,ijkCm1) - cv(CV_MIXT_DENS,ijkCm2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkCm1) - dv(DV_MIXT_UVEL,ijkCm2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkCm1) - dv(DV_MIXT_VVEL,ijkCm2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkCm1) - dv(DV_MIXT_WVEL,ijkCm2)
        dVarm(5) = dv(DV_MIXT_PRES,ijkCm1) - dv(DV_MIXT_PRES,ijkCm2)

        dVarp(1) = cv(CV_MIXT_DENS,ijkCp1) - cv(CV_MIXT_DENS,ijkC0 )
        dVarp(2) = dv(DV_MIXT_UVEL,ijkCp1) - dv(DV_MIXT_UVEL,ijkC0 )
        dVarp(3) = dv(DV_MIXT_VVEL,ijkCp1) - dv(DV_MIXT_VVEL,ijkC0 )
        dVarp(4) = dv(DV_MIXT_WVEL,ijkCp1) - dv(DV_MIXT_WVEL,ijkC0 )
        dVarp(5) = dv(DV_MIXT_PRES,ijkCp1) - dv(DV_MIXT_PRES,ijkC0 )

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

! ----- left and right states

        rl   = cv(CV_MIXT_DENS,ijkCm1) + deltl(1)
        ul   = dv(DV_MIXT_UVEL,ijkCm1) + deltl(2)
        vl   = dv(DV_MIXT_VVEL,ijkCm1) + deltl(3)
        wl   = dv(DV_MIXT_WVEL,ijkCm1) + deltl(4)
        pl   = dv(DV_MIXT_PRES,ijkCm1) + deltl(5)
        rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkCm1*indMol)
        gaml = gv(GV_MIXT_CP,ijkCm1*indCp)/(gv(GV_MIXT_CP,ijkCm1*indCp)-rgas)
        ggm1 = gaml/(gaml-1._RFREAL)
        hl   = ggm1*pl/rl + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)

        rr   = cv(CV_MIXT_DENS,ijkC0) - deltr(1)
        ur   = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr   = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr   = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        pr   = dv(DV_MIXT_PRES,ijkC0) - deltr(5)
        rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC0*indMol)
        gamr = gv(GV_MIXT_CP,ijkC0*indCp)/(gv(GV_MIXT_CP,ijkC0*indCp)-rgas)
        ggm1 = gamr/(gamr-1._RFREAL)
        hr   = ggm1*pr/rr + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)

! ----- Roe`s average

        rav   = SQRT(rl*rr)
        dd    = rav/rl
        dd1   = 1._RFREAL/(1._RFREAL+dd)
        uav   = (ul+dd*ur)*dd1
        vav   = (vl+dd*vr)*dd1
        wav   = (wl+dd*wr)*dd1
        hav   = (hl+dd*hr)*dd1
        q2a   = 0.5_RFREAL*(uav*uav+vav*vav+wav*wav)
        gam1  = 0.5_RFREAL*(gaml+gamr) - 1._RFREAL
        c2a   = gam1*(hav-q2a)
        cav   = SQRT(c2a)
        uvw   = uav*nx + vav*ny + wav*nz - sVel
        du    = (ur-ul)*nx + (vr-vl)*ny + (wr-wl)*nz

! ----- eigenvalues

        h1    = ABS(uvw - cav)
        h2    = ABS(uvw)
        eabs5 = ABS(uvw + cav)
        delta = epsEntr*eabs5
        eabs1 = Entropy_corr2( h1,delta )
        eabs2 = Entropy_corr2( h2,delta )

! ----- upwind fluxes

        h1 = rav*cav*du
        h2 = eabs1*(pr-pl - h1)/(2._RFREAL*c2a)
        h3 = eabs2*(rr-rl - (pr-pl)/c2a)
        h4 = eabs2*rav
        h5 = eabs5*(pr-pl + h1)/(2._RFREAL*c2a)

        fd(1)  = h2 + h3 + h5
        fd(2)  = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + &
                 h5*(uav+cav*nx)
        fd(3)  = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + &
                 h5*(vav+cav*ny)
        fd(4)  = h2*(wav-cav*nz) + h3*wav + h4*(wr-wl-du*nz) + &
                 h5*(wav+cav*nz)
        fd(5)  = h2*(hav-cav*uvw) + h3*q2a + h4*(uav*(ur-ul)+ &
                 vav*(vr-vl)+wav*(wr-wl)-uvw*du) + &
                 h5*(hav+cav*uvw)

        fd(1)  = fd(1)*dS
        fd(2)  = fd(2)*dS
        fd(3)  = fd(3)*dS
        fd(4)  = fd(4)*dS
        fd(5)  = fd(5)*dS

        diss(CV_MIXT_DENS,ijkC0 ) = diss(CV_MIXT_DENS,ijkC0 ) - fd(1)
        diss(CV_MIXT_XMOM,ijkC0 ) = diss(CV_MIXT_XMOM,ijkC0 ) - fd(2)
        diss(CV_MIXT_YMOM,ijkC0 ) = diss(CV_MIXT_YMOM,ijkC0 ) - fd(3)
        diss(CV_MIXT_ZMOM,ijkC0 ) = diss(CV_MIXT_ZMOM,ijkC0 ) - fd(4)
        diss(CV_MIXT_ENER,ijkC0 ) = diss(CV_MIXT_ENER,ijkC0 ) - fd(5)

        diss(CV_MIXT_DENS,ijkCm1) = diss(CV_MIXT_DENS,ijkCm1) + fd(1)
        diss(CV_MIXT_XMOM,ijkCm1) = diss(CV_MIXT_XMOM,ijkCm1) + fd(2)
        diss(CV_MIXT_YMOM,ijkCm1) = diss(CV_MIXT_YMOM,ijkCm1) + fd(3)
        diss(CV_MIXT_ZMOM,ijkCm1) = diss(CV_MIXT_ZMOM,ijkCm1) + fd(4)
        diss(CV_MIXT_ENER,ijkCm1) = diss(CV_MIXT_ENER,ijkCm1) + fd(5)
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! dissipation in j-direction --------------------------------------------------

  DO k=kpcbeg,kpcend
    DO i=ipcbeg,ipcend
      DO j=jpcbeg,jpcend+1
        ijkC0  = IndIJK(i,j  ,k,iCOff,ijCOff)
        ijkCm1 = IndIJK(i,j-1,k,iCOff,ijCOff)
        ijkCm2 = IndIJK(i,j-2,k,iCOff,ijCOff)
        ijkCp1 = IndIJK(i,j+1,k,iCOff,ijCOff)
        ijkN   = IndIJK(i,j  ,k,iNOff,ijNOff)
        dS     = SQRT(sj(XCOORD,ijkN)*sj(XCOORD,ijkN)+ &
                      sj(YCOORD,ijkN)*sj(YCOORD,ijkN)+ &
                      sj(ZCOORD,ijkN)*sj(ZCOORD,ijkN))
        nx     = sj(XCOORD,ijkN)/dS
        ny     = sj(YCOORD,ijkN)/dS
        nz     = sj(ZCOORD,ijkN)/dS
        sVel   = sjVel(ijkN*indSvel)/dS
        vola   = (0.5_RFREAL*(vol(ijkC0)+vol(ijkCm1)))**1.5_RFREAL
        dS     = dS*beta5

! ----- limited extrapolations

        dVar(1)  = cv(CV_MIXT_DENS,ijkC0 ) - cv(CV_MIXT_DENS,ijkCm1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0 ) - dv(DV_MIXT_UVEL,ijkCm1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0 ) - dv(DV_MIXT_VVEL,ijkCm1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0 ) - dv(DV_MIXT_WVEL,ijkCm1)
        dVar(5)  = dv(DV_MIXT_PRES,ijkC0 ) - dv(DV_MIXT_PRES,ijkCm1)

        dVarm(1) = cv(CV_MIXT_DENS,ijkCm1) - cv(CV_MIXT_DENS,ijkCm2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkCm1) - dv(DV_MIXT_UVEL,ijkCm2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkCm1) - dv(DV_MIXT_VVEL,ijkCm2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkCm1) - dv(DV_MIXT_WVEL,ijkCm2)
        dVarm(5) = dv(DV_MIXT_PRES,ijkCm1) - dv(DV_MIXT_PRES,ijkCm2)

        dVarp(1) = cv(CV_MIXT_DENS,ijkCp1) - cv(CV_MIXT_DENS,ijkC0 )
        dVarp(2) = dv(DV_MIXT_UVEL,ijkCp1) - dv(DV_MIXT_UVEL,ijkC0 )
        dVarp(3) = dv(DV_MIXT_VVEL,ijkCp1) - dv(DV_MIXT_VVEL,ijkC0 )
        dVarp(4) = dv(DV_MIXT_WVEL,ijkCp1) - dv(DV_MIXT_WVEL,ijkC0 )
        dVarp(5) = dv(DV_MIXT_PRES,ijkCp1) - dv(DV_MIXT_PRES,ijkC0 )

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

! ----- left and right states

        rl   = cv(CV_MIXT_DENS,ijkCm1) + deltl(1)
        ul   = dv(DV_MIXT_UVEL,ijkCm1) + deltl(2)
        vl   = dv(DV_MIXT_VVEL,ijkCm1) + deltl(3)
        wl   = dv(DV_MIXT_WVEL,ijkCm1) + deltl(4)
        pl   = dv(DV_MIXT_PRES,ijkCm1) + deltl(5)
        rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkCm1*indMol)
        gaml = gv(GV_MIXT_CP,ijkCm1*indCp)/(gv(GV_MIXT_CP,ijkCm1*indCp)-rgas)
        ggm1 = gaml/(gaml-1._RFREAL)
        hl   = ggm1*pl/rl + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)

        rr   = cv(CV_MIXT_DENS,ijkC0) - deltr(1)
        ur   = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr   = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr   = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        pr   = dv(DV_MIXT_PRES,ijkC0) - deltr(5)
        rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC0*indMol)
        gamr = gv(GV_MIXT_CP,ijkC0*indCp)/(gv(GV_MIXT_CP,ijkC0*indCp)-rgas)
        ggm1 = gamr/(gamr-1._RFREAL)
        hr   = ggm1*pr/rr + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)

! ----- Roe`s average

        rav   = SQRT(rl*rr)
        dd    = rav/rl
        dd1   = 1._RFREAL/(1._RFREAL+dd)
        uav   = (ul+dd*ur)*dd1
        vav   = (vl+dd*vr)*dd1
        wav   = (wl+dd*wr)*dd1
        hav   = (hl+dd*hr)*dd1
        q2a   = 0.5_RFREAL*(uav*uav+vav*vav+wav*wav)
        gam1  = 0.5_RFREAL*(gaml+gamr) - 1._RFREAL
        c2a   = gam1*(hav-q2a)
        cav   = SQRT(c2a)
        uvw   = uav*nx + vav*ny + wav*nz - sVel
        du    = (ur-ul)*nx + (vr-vl)*ny + (wr-wl)*nz

! ----- eigenvalues

        h1    = ABS(uvw - cav)
        h2    = ABS(uvw)
        eabs5 = ABS(uvw + cav)
        delta = epsEntr*eabs5
        eabs1 = Entropy_corr2( h1,delta )
        eabs2 = Entropy_corr2( h2,delta )

! ----- upwind fluxes

        h1 = rav*cav*du
        h2 = eabs1*(pr-pl - h1)/(2._RFREAL*c2a)
        h3 = eabs2*(rr-rl - (pr-pl)/c2a)
        h4 = eabs2*rav
        h5 = eabs5*(pr-pl + h1)/(2._RFREAL*c2a)

        fd(1)  = h2 + h3 + h5
        fd(2)  = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + &
                 h5*(uav+cav*nx)
        fd(3)  = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + &
                 h5*(vav+cav*ny)
        fd(4)  = h2*(wav-cav*nz) + h3*wav + h4*(wr-wl-du*nz) + &
                 h5*(wav+cav*nz)
        fd(5)  = h2*(hav-cav*uvw) + h3*q2a + h4*(uav*(ur-ul)+ &
                 vav*(vr-vl)+wav*(wr-wl)-uvw*du) + &
                 h5*(hav+cav*uvw)

        fd(1)  = fd(1)*dS
        fd(2)  = fd(2)*dS
        fd(3)  = fd(3)*dS
        fd(4)  = fd(4)*dS
        fd(5)  = fd(5)*dS

        diss(CV_MIXT_DENS,ijkC0 ) = diss(CV_MIXT_DENS,ijkC0 ) - fd(1)
        diss(CV_MIXT_XMOM,ijkC0 ) = diss(CV_MIXT_XMOM,ijkC0 ) - fd(2)
        diss(CV_MIXT_YMOM,ijkC0 ) = diss(CV_MIXT_YMOM,ijkC0 ) - fd(3)
        diss(CV_MIXT_ZMOM,ijkC0 ) = diss(CV_MIXT_ZMOM,ijkC0 ) - fd(4)
        diss(CV_MIXT_ENER,ijkC0 ) = diss(CV_MIXT_ENER,ijkC0 ) - fd(5)

        diss(CV_MIXT_DENS,ijkCm1) = diss(CV_MIXT_DENS,ijkCm1) + fd(1)
        diss(CV_MIXT_XMOM,ijkCm1) = diss(CV_MIXT_XMOM,ijkCm1) + fd(2)
        diss(CV_MIXT_YMOM,ijkCm1) = diss(CV_MIXT_YMOM,ijkCm1) + fd(3)
        diss(CV_MIXT_ZMOM,ijkCm1) = diss(CV_MIXT_ZMOM,ijkCm1) + fd(4)
        diss(CV_MIXT_ENER,ijkCm1) = diss(CV_MIXT_ENER,ijkCm1) + fd(5)
      ENDDO  ! j
    ENDDO    ! i
  ENDDO      ! k

! dissipation in k-direction --------------------------------------------------

  DO j=jpcbeg,jpcend
    DO i=ipcbeg,ipcend
      DO k=kpcbeg,kpcend+1
        ijkC0  = IndIJK(i,j,k  ,iCOff,ijCOff)
        ijkCm1 = IndIJK(i,j,k-1,iCOff,ijCOff)
        ijkCm2 = IndIJK(i,j,k-2,iCOff,ijCOff)
        ijkCp1 = IndIJK(i,j,k+1,iCOff,ijCOff)
        ijkN   = IndIJK(i,j,k  ,iNOff,ijNOff)
        dS     = SQRT(sk(XCOORD,ijkN)*sk(XCOORD,ijkN)+ &
                      sk(YCOORD,ijkN)*sk(YCOORD,ijkN)+ &
                      sk(ZCOORD,ijkN)*sk(ZCOORD,ijkN))
        nx     = sk(XCOORD,ijkN)/dS
        ny     = sk(YCOORD,ijkN)/dS
        nz     = sk(ZCOORD,ijkN)/dS
        sVel   = skVel(ijkN*indSvel)/dS
        vola   = (0.5_RFREAL*(vol(ijkC0)+vol(ijkCm1)))**1.5_RFREAL
        dS     = dS*beta5

! ----- limited extrapolations

        dVar(1)  = cv(CV_MIXT_DENS,ijkC0 ) - cv(CV_MIXT_DENS,ijkCm1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0 ) - dv(DV_MIXT_UVEL,ijkCm1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0 ) - dv(DV_MIXT_VVEL,ijkCm1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0 ) - dv(DV_MIXT_WVEL,ijkCm1)
        dVar(5)  = dv(DV_MIXT_PRES,ijkC0 ) - dv(DV_MIXT_PRES,ijkCm1)

        dVarm(1) = cv(CV_MIXT_DENS,ijkCm1) - cv(CV_MIXT_DENS,ijkCm2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkCm1) - dv(DV_MIXT_UVEL,ijkCm2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkCm1) - dv(DV_MIXT_VVEL,ijkCm2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkCm1) - dv(DV_MIXT_WVEL,ijkCm2)
        dVarm(5) = dv(DV_MIXT_PRES,ijkCm1) - dv(DV_MIXT_PRES,ijkCm2)

        dVarp(1) = cv(CV_MIXT_DENS,ijkCp1) - cv(CV_MIXT_DENS,ijkC0 )
        dVarp(2) = dv(DV_MIXT_UVEL,ijkCp1) - dv(DV_MIXT_UVEL,ijkC0 )
        dVarp(3) = dv(DV_MIXT_VVEL,ijkCp1) - dv(DV_MIXT_VVEL,ijkC0 )
        dVarp(4) = dv(DV_MIXT_WVEL,ijkCp1) - dv(DV_MIXT_WVEL,ijkC0 )
        dVarp(5) = dv(DV_MIXT_PRES,ijkCp1) - dv(DV_MIXT_PRES,ijkC0 )

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

! ----- left and right states

        rl   = cv(CV_MIXT_DENS,ijkCm1) + deltl(1)
        ul   = dv(DV_MIXT_UVEL,ijkCm1) + deltl(2)
        vl   = dv(DV_MIXT_VVEL,ijkCm1) + deltl(3)
        wl   = dv(DV_MIXT_WVEL,ijkCm1) + deltl(4)
        pl   = dv(DV_MIXT_PRES,ijkCm1) + deltl(5)
        rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkCm1*indMol)
        gaml = gv(GV_MIXT_CP,ijkCm1*indCp)/(gv(GV_MIXT_CP,ijkCm1*indCp)-rgas)
        ggm1 = gaml/(gaml-1._RFREAL)
        hl   = ggm1*pl/rl + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)

        rr   = cv(CV_MIXT_DENS,ijkC0) - deltr(1)
        ur   = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr   = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr   = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        pr   = dv(DV_MIXT_PRES,ijkC0) - deltr(5)
        rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC0*indMol)
        gamr = gv(GV_MIXT_CP,ijkC0*indCp)/(gv(GV_MIXT_CP,ijkC0*indCp)-rgas)
        ggm1 = gamr/(gamr-1._RFREAL)
        hr   = ggm1*pr/rr + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)

! ----- Roe`s average

        rav   = SQRT(rl*rr)
        dd    = rav/rl
        dd1   = 1._RFREAL/(1._RFREAL+dd)
        uav   = (ul+dd*ur)*dd1
        vav   = (vl+dd*vr)*dd1
        wav   = (wl+dd*wr)*dd1
        hav   = (hl+dd*hr)*dd1
        q2a   = 0.5_RFREAL*(uav*uav+vav*vav+wav*wav)
        gam1  = 0.5_RFREAL*(gaml+gamr) - 1._RFREAL
        c2a   = gam1*(hav-q2a)
        cav   = SQRT(c2a)
        uvw   = uav*nx + vav*ny + wav*nz - sVel
        du    = (ur-ul)*nx + (vr-vl)*ny + (wr-wl)*nz

! ----- eigenvalues

        h1    = ABS(uvw - cav)
        h2    = ABS(uvw)
        eabs5 = ABS(uvw + cav)
        delta = epsEntr*eabs5
        eabs1 = Entropy_corr2( h1,delta )
        eabs2 = Entropy_corr2( h2,delta )

! ----- upwind fluxes

        h1 = rav*cav*du
        h2 = eabs1*(pr-pl - h1)/(2._RFREAL*c2a)
        h3 = eabs2*(rr-rl - (pr-pl)/c2a)
        h4 = eabs2*rav
        h5 = eabs5*(pr-pl + h1)/(2._RFREAL*c2a)

        fd(1)  = h2 + h3 + h5
        fd(2)  = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + &
                 h5*(uav+cav*nx)
        fd(3)  = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + &
                 h5*(vav+cav*ny)
        fd(4)  = h2*(wav-cav*nz) + h3*wav + h4*(wr-wl-du*nz) + &
                 h5*(wav+cav*nz)
        fd(5)  = h2*(hav-cav*uvw) + h3*q2a + h4*(uav*(ur-ul)+ &
                 vav*(vr-vl)+wav*(wr-wl)-uvw*du) + &
                 h5*(hav+cav*uvw)

        fd(1)  = fd(1)*dS
        fd(2)  = fd(2)*dS
        fd(3)  = fd(3)*dS
        fd(4)  = fd(4)*dS
        fd(5)  = fd(5)*dS

        diss(CV_MIXT_DENS,ijkC0 ) = diss(CV_MIXT_DENS,ijkC0 ) - fd(1)
        diss(CV_MIXT_XMOM,ijkC0 ) = diss(CV_MIXT_XMOM,ijkC0 ) - fd(2)
        diss(CV_MIXT_YMOM,ijkC0 ) = diss(CV_MIXT_YMOM,ijkC0 ) - fd(3)
        diss(CV_MIXT_ZMOM,ijkC0 ) = diss(CV_MIXT_ZMOM,ijkC0 ) - fd(4)
        diss(CV_MIXT_ENER,ijkC0 ) = diss(CV_MIXT_ENER,ijkC0 ) - fd(5)

        diss(CV_MIXT_DENS,ijkCm1) = diss(CV_MIXT_DENS,ijkCm1) + fd(1)
        diss(CV_MIXT_XMOM,ijkCm1) = diss(CV_MIXT_XMOM,ijkCm1) + fd(2)
        diss(CV_MIXT_YMOM,ijkCm1) = diss(CV_MIXT_YMOM,ijkCm1) + fd(3)
        diss(CV_MIXT_ZMOM,ijkCm1) = diss(CV_MIXT_ZMOM,ijkCm1) + fd(4)
        diss(CV_MIXT_ENER,ijkCm1) = diss(CV_MIXT_ENER,ijkCm1) + fd(5)
      ENDDO  ! k
    ENDDO    ! i
  ENDDO      ! j

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

! #############################################################################

  CONTAINS

    REAL(RFREAL) FUNCTION Entropy_corr2( z,d )

      REAL(RFREAL) :: z, d

      IF (z > d) THEN
        Entropy_corr2 = z
      ELSE
        Entropy_corr2 = 0.5_RFREAL*(z*z+d*d)/d
      ENDIF

    END FUNCTION Entropy_corr2

END SUBROUTINE RFLO_RoeDissipSecond

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_RoeDissipSecond.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
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







