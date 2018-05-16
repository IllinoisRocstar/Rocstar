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
! Purpose: compute convective fluxes based on 1st-order Roe upwind scheme.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%rhs = convective fluxes added to the residual.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_RoeFluxFirst.F90,v 1.3 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_RoeFluxFirst( region )

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
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkC1, ijkN, indSvel

  REAL(RFREAL)          :: rhl, rhr, qsl, qsr, pav, fc(5)
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), rhs(:,:), si(:,:), sj(:,:), sk(:,:)
  REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_RoeFluxFirst',&
  'RFLO_RoeFluxFirst.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv
  rhs    => region%levels(iLev)%mixt%rhs
  si     => region%levels(iLev)%grid%si
  sj     => region%levels(iLev)%grid%sj
  sk     => region%levels(iLev)%grid%sk
  siVel  => region%levels(iLev)%grid%siVel
  sjVel  => region%levels(iLev)%grid%sjVel
  skVel  => region%levels(iLev)%grid%skVel
  indSvel = region%levels(iLev)%grid%indSvel

! flux in i-direction (except through boundary) -------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg+1,ipcend
        ijkC0 = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkN  = IndIJK(i  ,j,k,iNOff,ijNOff)

        rhl   = dv(DV_MIXT_PRES,ijkC1) + cv(CV_MIXT_ENER,ijkC1)
        qsl   = dv(DV_MIXT_UVEL,ijkC1)*si(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC1)*si(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC1)*si(ZCOORD,ijkN) - siVel(ijkN*indSvel)
        rhr   = dv(DV_MIXT_PRES,ijkC0) + cv(CV_MIXT_ENER,ijkC0)
        qsr   = dv(DV_MIXT_UVEL,ijkC0)*si(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC0)*si(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC0)*si(ZCOORD,ijkN) - siVel(ijkN*indSvel)
        pav   = 0.5_RFREAL*(dv(DV_MIXT_PRES,ijkC0)+dv(DV_MIXT_PRES,ijkC1))

        fc(1) = 0.5_RFREAL*(qsl*cv(CV_MIXT_DENS,ijkC1)+ &
                            qsr*cv(CV_MIXT_DENS,ijkC0))
        fc(2) = 0.5_RFREAL*(qsl*cv(CV_MIXT_XMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_XMOM,ijkC0)) + pav*si(XCOORD,ijkN)
        fc(3) = 0.5_RFREAL*(qsl*cv(CV_MIXT_YMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_YMOM,ijkC0)) + pav*si(YCOORD,ijkN)
        fc(4) = 0.5_RFREAL*(qsl*cv(CV_MIXT_ZMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_ZMOM,ijkC0)) + pav*si(ZCOORD,ijkN)
        fc(5) = 0.5_RFREAL*(qsl*rhl+qsr*rhr) + siVel(ijkN*indSvel)*pav

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
        ijkN  = IndIJK(i,j  ,k,iNOff,ijNOff)

        rhl   = dv(DV_MIXT_PRES,ijkC1) + cv(CV_MIXT_ENER,ijkC1)
        qsl   = dv(DV_MIXT_UVEL,ijkC1)*sj(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC1)*sj(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC1)*sj(ZCOORD,ijkN) - sjVel(ijkN*indSvel)
        rhr   = dv(DV_MIXT_PRES,ijkC0) + cv(CV_MIXT_ENER,ijkC0)
        qsr   = dv(DV_MIXT_UVEL,ijkC0)*sj(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC0)*sj(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC0)*sj(ZCOORD,ijkN) - sjVel(ijkN*indSvel)
        pav   = 0.5_RFREAL*(dv(DV_MIXT_PRES,ijkC0)+dv(DV_MIXT_PRES,ijkC1))

        fc(1) = 0.5_RFREAL*(qsl*cv(CV_MIXT_DENS,ijkC1)+ &
                            qsr*cv(CV_MIXT_DENS,ijkC0))
        fc(2) = 0.5_RFREAL*(qsl*cv(CV_MIXT_XMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_XMOM,ijkC0)) + pav*sj(XCOORD,ijkN)
        fc(3) = 0.5_RFREAL*(qsl*cv(CV_MIXT_YMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_YMOM,ijkC0)) + pav*sj(YCOORD,ijkN)
        fc(4) = 0.5_RFREAL*(qsl*cv(CV_MIXT_ZMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_ZMOM,ijkC0)) + pav*sj(ZCOORD,ijkN)
        fc(5) = 0.5_RFREAL*(qsl*rhl+qsr*rhr) + sjVel(ijkN*indSvel)*pav

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
        ijkN  = IndIJK(i,j,k  ,iNOff,ijNOff)

        rhl   = dv(DV_MIXT_PRES,ijkC1) + cv(CV_MIXT_ENER,ijkC1)
        qsl   = dv(DV_MIXT_UVEL,ijkC1)*sk(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC1)*sk(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC1)*sk(ZCOORD,ijkN) - skVel(ijkN*indSvel)
        rhr   = dv(DV_MIXT_PRES,ijkC0) + cv(CV_MIXT_ENER,ijkC0)
        qsr   = dv(DV_MIXT_UVEL,ijkC0)*sk(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC0)*sk(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC0)*sk(ZCOORD,ijkN) - skVel(ijkN*indSvel)
        pav   = 0.5_RFREAL*(dv(DV_MIXT_PRES,ijkC0)+dv(DV_MIXT_PRES,ijkC1))

        fc(1) = 0.5_RFREAL*(qsl*cv(CV_MIXT_DENS,ijkC1)+ &
                            qsr*cv(CV_MIXT_DENS,ijkC0))
        fc(2) = 0.5_RFREAL*(qsl*cv(CV_MIXT_XMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_XMOM,ijkC0)) + pav*sk(XCOORD,ijkN)
        fc(3) = 0.5_RFREAL*(qsl*cv(CV_MIXT_YMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_YMOM,ijkC0)) + pav*sk(YCOORD,ijkN)
        fc(4) = 0.5_RFREAL*(qsl*cv(CV_MIXT_ZMOM,ijkC1)+ &
                            qsr*cv(CV_MIXT_ZMOM,ijkC0)) + pav*sk(ZCOORD,ijkN)
        fc(5) = 0.5_RFREAL*(qsl*rhl+qsr*rhr) + skVel(ijkN*indSvel)*pav

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

END SUBROUTINE RFLO_RoeFluxFirst

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_RoeFluxFirst.F90,v $
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







