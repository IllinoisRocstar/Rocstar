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
! Purpose: compute central convective flux of SA eq. by average of variables.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%turb%rhs = convective fluxes added to SA residual.
!
! Notes: grid speeds si,j,kVel already contain face area
!
!******************************************************************************
!
! $Id: TURB_floRansSACentFlux.F90,v 1.5 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansSACentFlux( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloRansSACentFluxPatch
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

  REAL(RFREAL)          :: rnutila, ua, va, wa, vcont, sVel, fc
  REAL(RFREAL), POINTER :: dv(:,:), tcv(:,:), trhs(:,:)
  REAL(RFREAL), POINTER :: aci(:,:), acj(:,:), ack(:,:)
  REAL(RFREAL), POINTER :: si(:,:), sj(:,:), sk(:,:)
  REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_FloRansSACentFlux',&
  'TURB_floRansSACentFlux.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  dv     => region%levels(iLev)%mixt%dv
  tcv    => region%levels(iLev)%turb%cv
  trhs   => region%levels(iLev)%turb%rhs
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

        rnutila= aci(2,ijkN)*tcv(CV_SA_NUTIL,ijkC0)+ &
                 aci(1,ijkN)*tcv(CV_SA_NUTIL,ijkC1)
        ua     = aci(2,ijkN)*dv(DV_MIXT_UVEL,ijkC0)+ &
                 aci(1,ijkN)*dv(DV_MIXT_UVEL,ijkC1)
        va     = aci(2,ijkN)*dv(DV_MIXT_VVEL,ijkC0)+ &
                 aci(1,ijkN)*dv(DV_MIXT_VVEL,ijkC1)
        wa     = aci(2,ijkN)*dv(DV_MIXT_WVEL,ijkC0)+ &
                 aci(1,ijkN)*dv(DV_MIXT_WVEL,ijkC1)

        sVel   = siVel(ijkN*indSvel)
        vcont  = ua*si(XCOORD,ijkN)+va*si(YCOORD,ijkN)+wa*si(ZCOORD,ijkN)-sVel

        fc = vcont*rnutila

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

        rnutila= acj(2,ijkN)*tcv(CV_SA_NUTIL,ijkC0)+ &
                 acj(1,ijkN)*tcv(CV_SA_NUTIL,ijkC1)
        ua     = acj(2,ijkN)*dv(DV_MIXT_UVEL,ijkC0)+ &
                 acj(1,ijkN)*dv(DV_MIXT_UVEL,ijkC1)
        va     = acj(2,ijkN)*dv(DV_MIXT_VVEL,ijkC0)+ &
                 acj(1,ijkN)*dv(DV_MIXT_VVEL,ijkC1)
        wa     = acj(2,ijkN)*dv(DV_MIXT_WVEL,ijkC0)+ &
                 acj(1,ijkN)*dv(DV_MIXT_WVEL,ijkC1)

        sVel   = sjVel(ijkN*indSvel)
        vcont  = ua*sj(XCOORD,ijkN)+va*sj(YCOORD,ijkN)+wa*sj(ZCOORD,ijkN)-sVel

        fc = vcont*rnutila

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

        rnutila= ack(2,ijkN)*tcv(CV_SA_NUTIL,ijkC0)+ &
                 ack(1,ijkN)*tcv(CV_SA_NUTIL,ijkC1)
        ua     = ack(2,ijkN)*dv(DV_MIXT_UVEL,ijkC0)+ &
                 ack(1,ijkN)*dv(DV_MIXT_UVEL,ijkC1)
        va     = ack(2,ijkN)*dv(DV_MIXT_VVEL,ijkC0)+ &
                 ack(1,ijkN)*dv(DV_MIXT_VVEL,ijkC1)
        wa     = ack(2,ijkN)*dv(DV_MIXT_WVEL,ijkC0)+ &
                 ack(1,ijkN)*dv(DV_MIXT_WVEL,ijkC1)

        sVel   = skVel(ijkN*indSvel)
        vcont  = ua*sk(XCOORD,ijkN)+va*sk(YCOORD,ijkN)+wa*sk(ZCOORD,ijkN)-sVel

        fc = vcont*rnutila

        trhs(CV_SA_NUTIL,ijkC0) = trhs(CV_SA_NUTIL,ijkC0) + fc
        trhs(CV_SA_NUTIL,ijkC1) = trhs(CV_SA_NUTIL,ijkC1) - fc
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! fluxes through boundaries ---------------------------------------------------

  DO iPatch=1,region%nPatches
    CALL TURB_FloRansSACentFluxPatch( region, &
                                      region%levels(iLev)%patches(iPatch) )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_FloRansSACentFlux

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansSACentFlux.F90,v $
! Revision 1.5  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/08/02 21:55:07  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.3  2003/10/21 22:08:41  wasistho
! modified grid speed as it already contain face area
!
! Revision 1.2  2003/10/20 20:27:30  wasistho
! made consistent with compressible SA formulation
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







