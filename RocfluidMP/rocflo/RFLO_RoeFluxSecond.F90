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
! Purpose: compute convective fluxes based on 2nd-order Roe upwind scheme.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%rhs = convective fluxes added to the residual.
!
! Notes: uses MUSCL scheme with kappa=1/3.
!
!******************************************************************************
!
! $Id: RFLO_RoeFluxSecond.F90,v 1.3 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_RoeFluxSecond( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset, RFLO_RoeFluxPatch
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iPatch

! ... local variables
  INTEGER :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkC1, ijkC2, ijkC3, ijkN
  INTEGER :: indCp, indMol, indSvel

  REAL(RFREAL) :: limFac, limFac3, rVolRef, vola, eps2(3), dVar(5), dVarm(5), &
                  dVarp(5), rhl, rhr, qsl, qsr, pav, deltl(5), deltr(5), fc(5),&
                  eps2n, rgas, gam, ggm1, rl, rr, ul, ur, vl, vr, wl, wr, &
                  pl, pr, hl, hr, qsrl, qsrr
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), rhs(:,:), si(:,:), sj(:,:), sk(:,:)
  REAL(RFREAL), POINTER :: gv(:,:), vol(:), siVel(:), sjVel(:), skVel(:)

! ... functions
  REAL(RFREAL) :: MUSCL3, af, bf, eps

!******************************************************************************
! limiter function

  MUSCL3(af,bf,eps) = (bf*(2._RFREAL*af*af+eps)+af*(bf*bf+2._RFREAL*eps))/ &
                      (2._RFREAL*af*af+2._RFREAL*bf*bf-af*bf+ &
                       3._RFREAL*eps+1.E-30_RFREAL)

  CALL RegisterFunction( region%global,'RFLO_RoeFluxSecond',&
  'RFLO_RoeFluxSecond.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv
  gv     => region%levels(iLev)%mixt%gv
  rhs    => region%levels(iLev)%mixt%rhs
  vol    => region%levels(iLev)%grid%vol
  si     => region%levels(iLev)%grid%si
  sj     => region%levels(iLev)%grid%sj
  sk     => region%levels(iLev)%grid%sk
  siVel  => region%levels(iLev)%grid%siVel
  sjVel  => region%levels(iLev)%grid%sjVel
  skVel  => region%levels(iLev)%grid%skVel
  indSvel = region%levels(iLev)%grid%indSvel
  indCp   = region%levels(iLev)%mixt%indCp
  indMol  = region%levels(iLev)%mixt%indMol
  limFac  = region%mixtInput%limFac

! normalise epsilon^2 for all limited variables (rho, u, v, w, p) -------------

  limFac3 = limFac*limFac*limFac
  rVolRef = 1._RFREAL/region%global%limVolRef**1.5_RFREAL
  eps2(1) = limFac3*region%global%limRef(1)*region%global%limRef(1)*rVolRef
  eps2(2) = limFac3*region%global%limRef(2)*region%global%limRef(2)*rVolRef
  eps2(3) = limFac3*region%global%limRef(3)*region%global%limRef(3)*rVolRef

! flux in i-direction (except through boundary) -------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg+1,ipcend
        ijkC0 = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i-2,j,k,iCOff,ijCOff)
        ijkC3 = IndIJK(i+1,j,k,iCOff,ijCOff)
        ijkN  = IndIJK(i  ,j,k,iNOff,ijNOff)

        vola = (0.5_RFREAL*(vol(ijkC0)+vol(ijkC1)))**1.5_RFREAL

! ----- limited extrapolations

        dVar(1)  = cv(CV_MIXT_DENS,ijkC0) - cv(CV_MIXT_DENS,ijkC1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0) - dv(DV_MIXT_UVEL,ijkC1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0) - dv(DV_MIXT_VVEL,ijkC1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0) - dv(DV_MIXT_WVEL,ijkC1)
        dVar(5)  = dv(DV_MIXT_PRES,ijkC0) - dv(DV_MIXT_PRES,ijkC1)

        dVarm(1) = cv(CV_MIXT_DENS,ijkC1) - cv(CV_MIXT_DENS,ijkC2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkC1) - dv(DV_MIXT_UVEL,ijkC2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkC1) - dv(DV_MIXT_VVEL,ijkC2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkC1) - dv(DV_MIXT_WVEL,ijkC2)
        dVarm(5) = dv(DV_MIXT_PRES,ijkC1) - dv(DV_MIXT_PRES,ijkC2)

        dVarp(1) = cv(CV_MIXT_DENS,ijkC3) - cv(CV_MIXT_DENS,ijkC0)
        dVarp(2) = dv(DV_MIXT_UVEL,ijkC3) - dv(DV_MIXT_UVEL,ijkC0)
        dVarp(3) = dv(DV_MIXT_VVEL,ijkC3) - dv(DV_MIXT_VVEL,ijkC0)
        dVarp(4) = dv(DV_MIXT_WVEL,ijkC3) - dv(DV_MIXT_WVEL,ijkC0)
        dVarp(5) = dv(DV_MIXT_PRES,ijkC3) - dv(DV_MIXT_PRES,ijkC0)

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

        rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC1*indMol)
        gam   = gv(GV_MIXT_CP,ijkC1*indCp)/(gv(GV_MIXT_CP,ijkC1*indCp)-rgas)
        ggm1  = gam/(gam-1._RFREAL)
        rl    = cv(CV_MIXT_DENS,ijkC1) + deltl(1)
        ul    = dv(DV_MIXT_UVEL,ijkC1) + deltl(2)
        vl    = dv(DV_MIXT_VVEL,ijkC1) + deltl(3)
        wl    = dv(DV_MIXT_WVEL,ijkC1) + deltl(4)
        pl    = dv(DV_MIXT_PRES,ijkC1) + deltl(5)
        hl    = ggm1*pl/rl + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)
        qsrl  = (ul*si(XCOORD,ijkN)+vl*si(YCOORD,ijkN)+ &
                 wl*si(ZCOORD,ijkN)-siVel(ijkN*indSvel))*rl

        rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC0*indMol)
        gam   = gv(GV_MIXT_CP,ijkC0*indCp)/(gv(GV_MIXT_CP,ijkC0*indCp)-rgas)
        ggm1  = gam/(gam-1._RFREAL)
        rr    = cv(CV_MIXT_DENS,ijkC0) - deltr(1)
        ur    = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr    = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr    = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        pr    = dv(DV_MIXT_PRES,ijkC0) - deltr(5)
        hr    = ggm1*pr/rr + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)
        qsrr  = (ur*si(XCOORD,ijkN)+vr*si(YCOORD,ijkN)+ &
                 wr*si(ZCOORD,ijkN)-siVel(ijkN*indSvel))*rr

! ----- fluxes

        pav   = 0.5_RFREAL*(pl+pr)
        fc(1) = 0.5_RFREAL*(qsrl   +qsrr   )
        fc(2) = 0.5_RFREAL*(qsrl*ul+qsrr*ur) + si(XCOORD,ijkN)*pav
        fc(3) = 0.5_RFREAL*(qsrl*vl+qsrr*vr) + si(YCOORD,ijkN)*pav
        fc(4) = 0.5_RFREAL*(qsrl*wl+qsrr*wr) + si(ZCOORD,ijkN)*pav
        fc(5) = 0.5_RFREAL*(qsrl*hl+qsrr*hr) + siVel(ijkN*indSvel)*pav

        rhs(CV_MIXT_DENS,ijkC0) = rhs(CV_MIXT_DENS,ijkC0) + fc(1)
        rhs(CV_MIXT_XMOM,ijkC0) = rhs(CV_MIXT_XMOM,ijkC0) + fc(2)
        rhs(CV_MIXT_YMOM,ijkC0) = rhs(CV_MIXT_YMOM,ijkC0) + fc(3)
        rhs(CV_MIXT_ZMOM,ijkC0) = rhs(CV_MIXT_ZMOM,ijkC0) + fc(4)
        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) + fc(5)

        rhs(CV_MIXT_DENS,ijkC1) = rhs(CV_MIXT_DENS,ijkC1) - fc(1)
        rhs(CV_MIXT_XMOM,ijkC1) = rhs(CV_MIXT_XMOM,ijkC1) - fc(2)
        rhs(CV_MIXT_YMOM,ijkC1) = rhs(CV_MIXT_YMOM,ijkC1) - fc(3)
        rhs(CV_MIXT_ZMOM,ijkC1) = rhs(CV_MIXT_ZMOM,ijkC1) - fc(4)
        rhs(CV_MIXT_ENER,ijkC1) = rhs(CV_MIXT_ENER,ijkC1) - fc(5)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! flux in j-direction (except through boundary) -------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg+1,jpcend
      DO i=ipcbeg,ipcend
        ijkC0 = IndIJK(i,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j-1,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j-2,k,iCOff,ijCOff)
        ijkC3 = IndIJK(i,j+1,k,iCOff,ijCOff)
        ijkN  = IndIJK(i,j  ,k,iNOff,ijNOff)

        vola = (0.5_RFREAL*(vol(ijkC0)+vol(ijkC1)))**1.5_RFREAL

! ----- limited extrapolations

        dVar(1)  = cv(CV_MIXT_DENS,ijkC0) - cv(CV_MIXT_DENS,ijkC1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0) - dv(DV_MIXT_UVEL,ijkC1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0) - dv(DV_MIXT_VVEL,ijkC1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0) - dv(DV_MIXT_WVEL,ijkC1)
        dVar(5)  = dv(DV_MIXT_PRES,ijkC0) - dv(DV_MIXT_PRES,ijkC1)

        dVarm(1) = cv(CV_MIXT_DENS,ijkC1) - cv(CV_MIXT_DENS,ijkC2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkC1) - dv(DV_MIXT_UVEL,ijkC2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkC1) - dv(DV_MIXT_VVEL,ijkC2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkC1) - dv(DV_MIXT_WVEL,ijkC2)
        dVarm(5) = dv(DV_MIXT_PRES,ijkC1) - dv(DV_MIXT_PRES,ijkC2)

        dVarp(1) = cv(CV_MIXT_DENS,ijkC3) - cv(CV_MIXT_DENS,ijkC0)
        dVarp(2) = dv(DV_MIXT_UVEL,ijkC3) - dv(DV_MIXT_UVEL,ijkC0)
        dVarp(3) = dv(DV_MIXT_VVEL,ijkC3) - dv(DV_MIXT_VVEL,ijkC0)
        dVarp(4) = dv(DV_MIXT_WVEL,ijkC3) - dv(DV_MIXT_WVEL,ijkC0)
        dVarp(5) = dv(DV_MIXT_PRES,ijkC3) - dv(DV_MIXT_PRES,ijkC0)

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

        rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC1*indMol)
        gam   = gv(GV_MIXT_CP,ijkC1*indCp)/(gv(GV_MIXT_CP,ijkC1*indCp)-rgas)
        ggm1  = gam/(gam-1._RFREAL)
        rl    = cv(CV_MIXT_DENS,ijkC1) + deltl(1)
        ul    = dv(DV_MIXT_UVEL,ijkC1) + deltl(2)
        vl    = dv(DV_MIXT_VVEL,ijkC1) + deltl(3)
        wl    = dv(DV_MIXT_WVEL,ijkC1) + deltl(4)
        pl    = dv(DV_MIXT_PRES,ijkC1) + deltl(5)
        hl    = ggm1*pl/rl + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)
        qsrl  = (ul*sj(XCOORD,ijkN)+vl*sj(YCOORD,ijkN)+ &
                 wl*sj(ZCOORD,ijkN)-sjVel(ijkN*indSvel))*rl

        rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC0*indMol)
        gam   = gv(GV_MIXT_CP,ijkC0*indCp)/(gv(GV_MIXT_CP,ijkC0*indCp)-rgas)
        ggm1  = gam/(gam-1._RFREAL)
        rr    = cv(CV_MIXT_DENS,ijkC0) - deltr(1)
        ur    = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr    = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr    = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        pr    = dv(DV_MIXT_PRES,ijkC0) - deltr(5)
        hr    = ggm1*pr/rr + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)
        qsrr  = (ur*sj(XCOORD,ijkN)+vr*sj(YCOORD,ijkN)+ &
                 wr*sj(ZCOORD,ijkN)-sjVel(ijkN*indSvel))*rr

! ----- fluxes

        pav   = 0.5_RFREAL*(pl+pr)
        fc(1) = 0.5_RFREAL*(qsrl   +qsrr   )
        fc(2) = 0.5_RFREAL*(qsrl*ul+qsrr*ur) + sj(XCOORD,ijkN)*pav
        fc(3) = 0.5_RFREAL*(qsrl*vl+qsrr*vr) + sj(YCOORD,ijkN)*pav
        fc(4) = 0.5_RFREAL*(qsrl*wl+qsrr*wr) + sj(ZCOORD,ijkN)*pav
        fc(5) = 0.5_RFREAL*(qsrl*hl+qsrr*hr) + sjVel(ijkN*indSvel)*pav

        rhs(CV_MIXT_DENS,ijkC0) = rhs(CV_MIXT_DENS,ijkC0) + fc(1)
        rhs(CV_MIXT_XMOM,ijkC0) = rhs(CV_MIXT_XMOM,ijkC0) + fc(2)
        rhs(CV_MIXT_YMOM,ijkC0) = rhs(CV_MIXT_YMOM,ijkC0) + fc(3)
        rhs(CV_MIXT_ZMOM,ijkC0) = rhs(CV_MIXT_ZMOM,ijkC0) + fc(4)
        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) + fc(5)

        rhs(CV_MIXT_DENS,ijkC1) = rhs(CV_MIXT_DENS,ijkC1) - fc(1)
        rhs(CV_MIXT_XMOM,ijkC1) = rhs(CV_MIXT_XMOM,ijkC1) - fc(2)
        rhs(CV_MIXT_YMOM,ijkC1) = rhs(CV_MIXT_YMOM,ijkC1) - fc(3)
        rhs(CV_MIXT_ZMOM,ijkC1) = rhs(CV_MIXT_ZMOM,ijkC1) - fc(4)
        rhs(CV_MIXT_ENER,ijkC1) = rhs(CV_MIXT_ENER,ijkC1) - fc(5)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! flux in k-direction (except through boundary) -------------------------------

  DO k=kpcbeg+1,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        ijkC0 = IndIJK(i,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j,k-1,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j,k-2,iCOff,ijCOff)
        ijkC3 = IndIJK(i,j,k+1,iCOff,ijCOff)
        ijkN  = IndIJK(i,j,k  ,iNOff,ijNOff)

        vola = (0.5_RFREAL*(vol(ijkC0)+vol(ijkC1)))**1.5_RFREAL

! ----- limited extrapolations

        dVar(1)  = cv(CV_MIXT_DENS,ijkC0) - cv(CV_MIXT_DENS,ijkC1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0) - dv(DV_MIXT_UVEL,ijkC1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0) - dv(DV_MIXT_VVEL,ijkC1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0) - dv(DV_MIXT_WVEL,ijkC1)
        dVar(5)  = dv(DV_MIXT_PRES,ijkC0) - dv(DV_MIXT_PRES,ijkC1)

        dVarm(1) = cv(CV_MIXT_DENS,ijkC1) - cv(CV_MIXT_DENS,ijkC2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkC1) - dv(DV_MIXT_UVEL,ijkC2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkC1) - dv(DV_MIXT_VVEL,ijkC2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkC1) - dv(DV_MIXT_WVEL,ijkC2)
        dVarm(5) = dv(DV_MIXT_PRES,ijkC1) - dv(DV_MIXT_PRES,ijkC2)

        dVarp(1) = cv(CV_MIXT_DENS,ijkC3) - cv(CV_MIXT_DENS,ijkC0)
        dVarp(2) = dv(DV_MIXT_UVEL,ijkC3) - dv(DV_MIXT_UVEL,ijkC0)
        dVarp(3) = dv(DV_MIXT_VVEL,ijkC3) - dv(DV_MIXT_VVEL,ijkC0)
        dVarp(4) = dv(DV_MIXT_WVEL,ijkC3) - dv(DV_MIXT_WVEL,ijkC0)
        dVarp(5) = dv(DV_MIXT_PRES,ijkC3) - dv(DV_MIXT_PRES,ijkC0)

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

        rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC1*indMol)
        gam   = gv(GV_MIXT_CP,ijkC1*indCp)/(gv(GV_MIXT_CP,ijkC1*indCp)-rgas)
        ggm1  = gam/(gam-1._RFREAL)
        rl    = cv(CV_MIXT_DENS,ijkC1) + deltl(1)
        ul    = dv(DV_MIXT_UVEL,ijkC1) + deltl(2)
        vl    = dv(DV_MIXT_VVEL,ijkC1) + deltl(3)
        wl    = dv(DV_MIXT_WVEL,ijkC1) + deltl(4)
        pl    = dv(DV_MIXT_PRES,ijkC1) + deltl(5)
        hl    = ggm1*pl/rl + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)
        qsrl  = (ul*sk(XCOORD,ijkN)+vl*sk(YCOORD,ijkN)+ &
                 wl*sk(ZCOORD,ijkN)-skVel(ijkN*indSvel))*rl

        rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC0*indMol)
        gam   = gv(GV_MIXT_CP,ijkC0*indCp)/(gv(GV_MIXT_CP,ijkC0*indCp)-rgas)
        ggm1  = gam/(gam-1._RFREAL)
        rr    = cv(CV_MIXT_DENS,ijkC0) - deltr(1)
        ur    = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr    = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr    = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        pr    = dv(DV_MIXT_PRES,ijkC0) - deltr(5)
        hr    = ggm1*pr/rr + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)
        qsrr  = (ur*sk(XCOORD,ijkN)+vr*sk(YCOORD,ijkN)+ &
                 wr*sk(ZCOORD,ijkN)-skVel(ijkN*indSvel))*rr

! ----- fluxes

        pav   = 0.5_RFREAL*(pl+pr)
        fc(1) = 0.5_RFREAL*(qsrl   +qsrr   )
        fc(2) = 0.5_RFREAL*(qsrl*ul+qsrr*ur) + sk(XCOORD,ijkN)*pav
        fc(3) = 0.5_RFREAL*(qsrl*vl+qsrr*vr) + sk(YCOORD,ijkN)*pav
        fc(4) = 0.5_RFREAL*(qsrl*wl+qsrr*wr) + sk(ZCOORD,ijkN)*pav
        fc(5) = 0.5_RFREAL*(qsrl*hl+qsrr*hr) + skVel(ijkN*indSvel)*pav

        rhs(CV_MIXT_DENS,ijkC0) = rhs(CV_MIXT_DENS,ijkC0) + fc(1)
        rhs(CV_MIXT_XMOM,ijkC0) = rhs(CV_MIXT_XMOM,ijkC0) + fc(2)
        rhs(CV_MIXT_YMOM,ijkC0) = rhs(CV_MIXT_YMOM,ijkC0) + fc(3)
        rhs(CV_MIXT_ZMOM,ijkC0) = rhs(CV_MIXT_ZMOM,ijkC0) + fc(4)
        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) + fc(5)

        rhs(CV_MIXT_DENS,ijkC1) = rhs(CV_MIXT_DENS,ijkC1) - fc(1)
        rhs(CV_MIXT_XMOM,ijkC1) = rhs(CV_MIXT_XMOM,ijkC1) - fc(2)
        rhs(CV_MIXT_YMOM,ijkC1) = rhs(CV_MIXT_YMOM,ijkC1) - fc(3)
        rhs(CV_MIXT_ZMOM,ijkC1) = rhs(CV_MIXT_ZMOM,ijkC1) - fc(4)
        rhs(CV_MIXT_ENER,ijkC1) = rhs(CV_MIXT_ENER,ijkC1) - fc(5)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! fluxes through boundaries ---------------------------------------------------

  DO iPatch=1,region%nPatches
    CALL RFLO_RoeFluxPatch( region,region%levels(iLev)%patches(iPatch) )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_RoeFluxSecond

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_RoeFluxSecond.F90,v $
! Revision 1.3  2008/12/06 08:44:28  mtcampbe
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







