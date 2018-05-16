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
! Purpose: compute central convective flux FLD eq. by average of variables.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%radi%rhs = convective fluxes added to FLD residual.
!
! Notes: grid speeds si,j,kVel already contain face area
!
!******************************************************************************
!
! $Id: RADI_floFlimCentFlux.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimCentFlux( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE RADI_ModInterfaces, ONLY : RADI_FloFlimCentFluxPatch
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iPatch

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkC1, ijkN, indSvel

  REAL(RFREAL)          :: radEa, ua, va, wa, vcont, sVel, fc
  REAL(RFREAL), POINTER :: dv(:,:), rcv(:,:), rrhs(:,:)
  REAL(RFREAL), POINTER :: aci(:,:), acj(:,:), ack(:,:)
  REAL(RFREAL), POINTER :: si(:,:), sj(:,:), sk(:,:)
  REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FloFlimCentFlux',&
  'RADI_floFlimCentFlux.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  dv     => region%levels(iLev)%mixt%dv
  rcv    => region%levels(iLev)%radi%cv
  rrhs   => region%levels(iLev)%radi%rhs
  aci    => region%levels(iLev)%grid%c2fCoI
  acj    => region%levels(iLev)%grid%c2fCoJ
  ack    => region%levels(iLev)%grid%c2fCoK
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

        radEa  = aci(2,ijkN)*rcv(CV_RADI_ENER,ijkC0)+ &
                 aci(1,ijkN)*rcv(CV_RADI_ENER,ijkC1)
        ua     = aci(2,ijkN)*dv(DV_MIXT_UVEL,ijkC0)+ &
                 aci(1,ijkN)*dv(DV_MIXT_UVEL,ijkC1)
        va     = aci(2,ijkN)*dv(DV_MIXT_VVEL,ijkC0)+ &
                 aci(1,ijkN)*dv(DV_MIXT_VVEL,ijkC1)
        wa     = aci(2,ijkN)*dv(DV_MIXT_WVEL,ijkC0)+ &
                 aci(1,ijkN)*dv(DV_MIXT_WVEL,ijkC1)

        sVel   = siVel(ijkN*indSvel)
        vcont  = ua*si(XCOORD,ijkN)+va*si(YCOORD,ijkN)+wa*si(ZCOORD,ijkN)-sVel

        fc = vcont*radEa

        rrhs(CV_RADI_ENER,ijkC0) = rrhs(CV_RADI_ENER,ijkC0) + fc
        rrhs(CV_RADI_ENER,ijkC1) = rrhs(CV_RADI_ENER,ijkC1) - fc
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

        radEa  = acj(2,ijkN)*rcv(CV_RADI_ENER,ijkC0)+ &
                 acj(1,ijkN)*rcv(CV_RADI_ENER,ijkC1)
        ua     = acj(2,ijkN)*dv(DV_MIXT_UVEL,ijkC0)+ &
                 acj(1,ijkN)*dv(DV_MIXT_UVEL,ijkC1)
        va     = acj(2,ijkN)*dv(DV_MIXT_VVEL,ijkC0)+ &
                 acj(1,ijkN)*dv(DV_MIXT_VVEL,ijkC1)
        wa     = acj(2,ijkN)*dv(DV_MIXT_WVEL,ijkC0)+ &
                 acj(1,ijkN)*dv(DV_MIXT_WVEL,ijkC1)

        sVel   = sjVel(ijkN*indSvel)
        vcont  = ua*sj(XCOORD,ijkN)+va*sj(YCOORD,ijkN)+wa*sj(ZCOORD,ijkN)-sVel

        fc = vcont*radEa

        rrhs(CV_RADI_ENER,ijkC0) = rrhs(CV_RADI_ENER,ijkC0) + fc
        rrhs(CV_RADI_ENER,ijkC1) = rrhs(CV_RADI_ENER,ijkC1) - fc
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

        radEa  = ack(2,ijkN)*rcv(CV_RADI_ENER,ijkC0)+ &
                 ack(1,ijkN)*rcv(CV_RADI_ENER,ijkC1)
        ua     = ack(2,ijkN)*dv(DV_MIXT_UVEL,ijkC0)+ &
                 ack(1,ijkN)*dv(DV_MIXT_UVEL,ijkC1)
        va     = ack(2,ijkN)*dv(DV_MIXT_VVEL,ijkC0)+ &
                 ack(1,ijkN)*dv(DV_MIXT_VVEL,ijkC1)
        wa     = ack(2,ijkN)*dv(DV_MIXT_WVEL,ijkC0)+ &
                 ack(1,ijkN)*dv(DV_MIXT_WVEL,ijkC1)

        sVel   = skVel(ijkN*indSvel)
        vcont  = ua*sk(XCOORD,ijkN)+va*sk(YCOORD,ijkN)+wa*sk(ZCOORD,ijkN)-sVel

        fc = vcont*radEa

        rrhs(CV_RADI_ENER,ijkC0) = rrhs(CV_RADI_ENER,ijkC0) + fc
        rrhs(CV_RADI_ENER,ijkC1) = rrhs(CV_RADI_ENER,ijkC1) - fc
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! fluxes through boundaries ---------------------------------------------------

  DO iPatch=1,region%nPatches
    CALL RADI_FloFlimCentFluxPatch( region, &
                                    region%levels(iLev)%patches(iPatch) )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_FloFlimCentFlux

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimCentFlux.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







