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
! Purpose: compute SA convective fluxes based on 1st-order Roe upwind scheme.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%turb%rhs = convective fluxes added to the residual.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_floRansSARoe1stFlux.F90,v 1.4 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansSARoe1stFlux( region )

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
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkC1, ijkN, indSvel

  REAL(RFREAL)          :: qsl, qsr, fc
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), si(:,:), sj(:,:), sk(:,:)
  REAL(RFREAL), POINTER :: tcv(:,:), trhs(:,:), siVel(:), sjVel(:), skVel(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_FloRansSARoe1stFlux',&
  'TURB_floRansSARoe1stFlux.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv
  tcv    => region%levels(iLev)%turb%cv
  trhs   => region%levels(iLev)%turb%rhs
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

        qsl   = dv(DV_MIXT_UVEL,ijkC1)*si(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC1)*si(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC1)*si(ZCOORD,ijkN) - siVel(ijkN*indSvel)
        qsr   = dv(DV_MIXT_UVEL,ijkC0)*si(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC0)*si(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC0)*si(ZCOORD,ijkN) - siVel(ijkN*indSvel)

        fc    = 0.5_RFREAL*(qsl*tcv(CV_SA_NUTIL,ijkC1)+ &
                            qsr*tcv(CV_SA_NUTIL,ijkC0))

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
        ijkN  = IndIJK(i,j  ,k,iNOff,ijNOff)

        qsl   = dv(DV_MIXT_UVEL,ijkC1)*sj(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC1)*sj(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC1)*sj(ZCOORD,ijkN) - sjVel(ijkN*indSvel)
        qsr   = dv(DV_MIXT_UVEL,ijkC0)*sj(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC0)*sj(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC0)*sj(ZCOORD,ijkN) - sjVel(ijkN*indSvel)

        fc    = 0.5_RFREAL*(qsl*tcv(CV_SA_NUTIL,ijkC1)+ &
                            qsr*tcv(CV_SA_NUTIL,ijkC0))

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
        ijkN  = IndIJK(i,j,k  ,iNOff,ijNOff)

        qsl   = dv(DV_MIXT_UVEL,ijkC1)*sk(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC1)*sk(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC1)*sk(ZCOORD,ijkN) - skVel(ijkN*indSvel)
        qsr   = dv(DV_MIXT_UVEL,ijkC0)*sk(XCOORD,ijkN) + &
                dv(DV_MIXT_VVEL,ijkC0)*sk(YCOORD,ijkN) + &
                dv(DV_MIXT_WVEL,ijkC0)*sk(ZCOORD,ijkN) - skVel(ijkN*indSvel)

        fc    = 0.5_RFREAL*(qsl*tcv(CV_SA_NUTIL,ijkC1)+ &
                            qsr*tcv(CV_SA_NUTIL,ijkC0))

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

END SUBROUTINE TURB_FloRansSARoe1stFlux

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansSARoe1stFlux.F90,v $
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
! Revision 1.1  2003/10/27 04:54:38  wasistho
! added RaNS upwind schemes
!
!
!******************************************************************************







