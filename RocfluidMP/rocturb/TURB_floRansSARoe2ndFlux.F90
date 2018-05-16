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
! Purpose: compute SA convective fluxes based on 2nd-order Roe upwind scheme.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%turb%rhs = convective fluxes added to the residual.
!
! Notes: uses MUSCL scheme with kappa=1/3.
!
!******************************************************************************
!
! $Id: TURB_floRansSARoe2ndFlux.F90,v 1.4 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansSARoe2ndFlux( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloRansSARoeFluxPatch
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iPatch

! ... local variables
  INTEGER :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkC1, ijkC2, ijkC3, ijkN
  INTEGER :: indSvel

  REAL(RFREAL) :: limFac, limFac3, rVolRef, vola, eps2(3), dVar(4), dVarm(4), &
                  dVarp(4), qsl, qsr, deltl(4), deltr(4), fc, eps2n, &
                  rl, rr, ul, ur, vl, vr, wl, wr, qsrl, qsrr
  REAL(RFREAL), POINTER :: dv(:,:), si(:,:), sj(:,:), sk(:,:), vol(:)
  REAL(RFREAL), POINTER :: tcv(:,:), trhs(:,:), siVel(:), sjVel(:), skVel(:)

! ... functions
  REAL(RFREAL) :: MUSCL3, af, bf, eps

!******************************************************************************
! limiter function

  MUSCL3(af,bf,eps) = (bf*(2._RFREAL*af*af+eps)+af*(bf*bf+2._RFREAL*eps))/ &
                      (2._RFREAL*af*af+2._RFREAL*bf*bf-af*bf+ &
                       3._RFREAL*eps+1.E-30_RFREAL)

  CALL RegisterFunction( region%global,'TURB_FloRansSARoe2ndFlux',&
  'TURB_floRansSARoe2ndFlux.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  dv     => region%levels(iLev)%mixt%dv
  tcv    => region%levels(iLev)%turb%cv
  trhs   => region%levels(iLev)%turb%rhs
  vol    => region%levels(iLev)%grid%vol
  si     => region%levels(iLev)%grid%si
  sj     => region%levels(iLev)%grid%sj
  sk     => region%levels(iLev)%grid%sk
  siVel  => region%levels(iLev)%grid%siVel
  sjVel  => region%levels(iLev)%grid%sjVel
  skVel  => region%levels(iLev)%grid%skVel
  indSvel = region%levels(iLev)%grid%indSvel
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

        dVar(1)  = tcv(CV_SA_NUTIL,ijkC0) - tcv(CV_SA_NUTIL,ijkC1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0) - dv(DV_MIXT_UVEL,ijkC1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0) - dv(DV_MIXT_VVEL,ijkC1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0) - dv(DV_MIXT_WVEL,ijkC1)

        dVarm(1) = tcv(CV_SA_NUTIL,ijkC1) - tcv(CV_SA_NUTIL,ijkC2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkC1) - dv(DV_MIXT_UVEL,ijkC2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkC1) - dv(DV_MIXT_VVEL,ijkC2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkC1) - dv(DV_MIXT_WVEL,ijkC2)

        dVarp(1) = tcv(CV_SA_NUTIL,ijkC3) - tcv(CV_SA_NUTIL,ijkC0)
        dVarp(2) = dv(DV_MIXT_UVEL,ijkC3) - dv(DV_MIXT_UVEL,ijkC0)
        dVarp(3) = dv(DV_MIXT_VVEL,ijkC3) - dv(DV_MIXT_VVEL,ijkC0)
        dVarp(4) = dv(DV_MIXT_WVEL,ijkC3) - dv(DV_MIXT_WVEL,ijkC0)

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

! ----- left and right states

        rl    = tcv(CV_SA_NUTIL,ijkC1) + deltl(1)
        ul    = dv(DV_MIXT_UVEL,ijkC1) + deltl(2)
        vl    = dv(DV_MIXT_VVEL,ijkC1) + deltl(3)
        wl    = dv(DV_MIXT_WVEL,ijkC1) + deltl(4)
        qsrl  = (ul*si(XCOORD,ijkN)+vl*si(YCOORD,ijkN)+ &
                 wl*si(ZCOORD,ijkN)-siVel(ijkN*indSvel))*rl

        rr    = tcv(CV_SA_NUTIL,ijkC0) - deltr(1)
        ur    = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr    = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr    = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        qsrr  = (ur*si(XCOORD,ijkN)+vr*si(YCOORD,ijkN)+ &
                 wr*si(ZCOORD,ijkN)-siVel(ijkN*indSvel))*rr

! ----- fluxes

        fc    = 0.5_RFREAL*(qsrl   +qsrr   )
        trhs(CV_SA_NUTIL,ijkC0) = trhs(CV_SA_NUTIL,ijkC0) + fc
        trhs(CV_SA_NUTIL,ijkC1) = trhs(CV_SA_NUTIL,ijkC1) - fc
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

        dVar(1)  = tcv(CV_SA_NUTIL,ijkC0) - tcv(CV_SA_NUTIL,ijkC1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0) - dv(DV_MIXT_UVEL,ijkC1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0) - dv(DV_MIXT_VVEL,ijkC1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0) - dv(DV_MIXT_WVEL,ijkC1)

        dVarm(1) = tcv(CV_SA_NUTIL,ijkC1) - tcv(CV_SA_NUTIL,ijkC2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkC1) - dv(DV_MIXT_UVEL,ijkC2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkC1) - dv(DV_MIXT_VVEL,ijkC2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkC1) - dv(DV_MIXT_WVEL,ijkC2)

        dVarp(1) = tcv(CV_SA_NUTIL,ijkC3) - tcv(CV_SA_NUTIL,ijkC0)
        dVarp(2) = dv(DV_MIXT_UVEL,ijkC3) - dv(DV_MIXT_UVEL,ijkC0)
        dVarp(3) = dv(DV_MIXT_VVEL,ijkC3) - dv(DV_MIXT_VVEL,ijkC0)
        dVarp(4) = dv(DV_MIXT_WVEL,ijkC3) - dv(DV_MIXT_WVEL,ijkC0)

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

! ----- left and right states

        rl    = tcv(CV_SA_NUTIL,ijkC1) + deltl(1)
        ul    = dv(DV_MIXT_UVEL,ijkC1) + deltl(2)
        vl    = dv(DV_MIXT_VVEL,ijkC1) + deltl(3)
        wl    = dv(DV_MIXT_WVEL,ijkC1) + deltl(4)
        qsrl  = (ul*sj(XCOORD,ijkN)+vl*sj(YCOORD,ijkN)+ &
                 wl*sj(ZCOORD,ijkN)-sjVel(ijkN*indSvel))*rl

        rr    = tcv(CV_SA_NUTIL,ijkC0) - deltr(1)
        ur    = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr    = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr    = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        qsrr  = (ur*sj(XCOORD,ijkN)+vr*sj(YCOORD,ijkN)+ &
                 wr*sj(ZCOORD,ijkN)-sjVel(ijkN*indSvel))*rr

! ----- fluxes

        fc    = 0.5_RFREAL*(qsrl   +qsrr   )
        trhs(CV_SA_NUTIL,ijkC0) = trhs(CV_SA_NUTIL,ijkC0) + fc
        trhs(CV_SA_NUTIL,ijkC1) = trhs(CV_SA_NUTIL,ijkC1) - fc
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

        dVar(1)  = tcv(CV_SA_NUTIL,ijkC0) - tcv(CV_SA_NUTIL,ijkC1)
        dVar(2)  = dv(DV_MIXT_UVEL,ijkC0) - dv(DV_MIXT_UVEL,ijkC1)
        dVar(3)  = dv(DV_MIXT_VVEL,ijkC0) - dv(DV_MIXT_VVEL,ijkC1)
        dVar(4)  = dv(DV_MIXT_WVEL,ijkC0) - dv(DV_MIXT_WVEL,ijkC1)

        dVarm(1) = tcv(CV_SA_NUTIL,ijkC1) - tcv(CV_SA_NUTIL,ijkC2)
        dVarm(2) = dv(DV_MIXT_UVEL,ijkC1) - dv(DV_MIXT_UVEL,ijkC2)
        dVarm(3) = dv(DV_MIXT_VVEL,ijkC1) - dv(DV_MIXT_VVEL,ijkC2)
        dVarm(4) = dv(DV_MIXT_WVEL,ijkC1) - dv(DV_MIXT_WVEL,ijkC2)

        dVarp(1) = tcv(CV_SA_NUTIL,ijkC3) - tcv(CV_SA_NUTIL,ijkC0)
        dVarp(2) = dv(DV_MIXT_UVEL,ijkC3) - dv(DV_MIXT_UVEL,ijkC0)
        dVarp(3) = dv(DV_MIXT_VVEL,ijkC3) - dv(DV_MIXT_VVEL,ijkC0)
        dVarp(4) = dv(DV_MIXT_WVEL,ijkC3) - dv(DV_MIXT_WVEL,ijkC0)

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

! ----- left and right states

        rl    = tcv(CV_SA_NUTIL,ijkC1) + deltl(1)
        ul    = dv(DV_MIXT_UVEL,ijkC1) + deltl(2)
        vl    = dv(DV_MIXT_VVEL,ijkC1) + deltl(3)
        wl    = dv(DV_MIXT_WVEL,ijkC1) + deltl(4)
        qsrl  = (ul*sk(XCOORD,ijkN)+vl*sk(YCOORD,ijkN)+ &
                 wl*sk(ZCOORD,ijkN)-skVel(ijkN*indSvel))*rl

        rr    = tcv(CV_SA_NUTIL,ijkC0) - deltr(1)
        ur    = dv(DV_MIXT_UVEL,ijkC0) - deltr(2)
        vr    = dv(DV_MIXT_VVEL,ijkC0) - deltr(3)
        wr    = dv(DV_MIXT_WVEL,ijkC0) - deltr(4)
        qsrr  = (ur*sk(XCOORD,ijkN)+vr*sk(YCOORD,ijkN)+ &
                 wr*sk(ZCOORD,ijkN)-skVel(ijkN*indSvel))*rr

! ----- fluxes

        fc    = 0.5_RFREAL*(qsrl   +qsrr   )
        trhs(CV_SA_NUTIL,ijkC0) = trhs(CV_SA_NUTIL,ijkC0) + fc
        trhs(CV_SA_NUTIL,ijkC1) = trhs(CV_SA_NUTIL,ijkC1) - fc
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! fluxes through boundaries ---------------------------------------------------

  DO iPatch=1,region%nPatches
    CALL TURB_FloRansSARoeFluxPatch( region, &
                                     region%levels(iLev)%patches(iPatch) )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_FloRansSARoe2ndFlux

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansSARoe2ndFlux.F90,v $
! Revision 1.4  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/10/27 23:10:26  wasistho
! established RaNS 2nd order upwind scheme
!
! Revision 1.1  2003/10/27 04:54:38  wasistho
! added RaNS upwind schemes
!
!
!******************************************************************************







