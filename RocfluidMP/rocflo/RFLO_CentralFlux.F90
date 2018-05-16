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
! Purpose: compute central convective fluxes by average of variables.
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
! $Id: RFLO_CentralFlux.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CentralFlux( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
        RFLO_GetNodeOffset, RFLO_CentralFluxPatch
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

  REAL(RFREAL)          :: rhoa, rhoua, rhova, rhowa, rhoea, pa, vcont, fc(5)
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), rhs(:,:), si(:,:), sj(:,:), sk(:,:)
  REAL(RFREAL), POINTER :: aci(:,:), acj(:,:), ack(:,:)
  REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CentralFlux',&
  'RFLO_CentralFlux.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv
  rhs    => region%levels(iLev)%mixt%rhs
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

        rhoa  = aci(2,ijkN)*cv(CV_MIXT_DENS,ijkC0)+ &
                aci(1,ijkN)*cv(CV_MIXT_DENS,ijkC1)
        rhoua = aci(2,ijkN)*cv(CV_MIXT_XMOM,ijkC0)+ &
                aci(1,ijkN)*cv(CV_MIXT_XMOM,ijkC1)
        rhova = aci(2,ijkN)*cv(CV_MIXT_YMOM,ijkC0)+ &
                aci(1,ijkN)*cv(CV_MIXT_YMOM,ijkC1)
        rhowa = aci(2,ijkN)*cv(CV_MIXT_ZMOM,ijkC0)+ &
                aci(1,ijkN)*cv(CV_MIXT_ZMOM,ijkC1)
        rhoea = aci(2,ijkN)*cv(CV_MIXT_ENER,ijkC0)+ &
                aci(1,ijkN)*cv(CV_MIXT_ENER,ijkC1)
        pa    = aci(2,ijkN)*dv(DV_MIXT_PRES,ijkC0)+ &
                aci(1,ijkN)*dv(DV_MIXT_PRES,ijkC1)
        vcont = (rhoua*si(XCOORD,ijkN)+rhova*si(YCOORD,ijkN)+&
                 rhowa*si(ZCOORD,ijkN))/rhoa - siVel(ijkN*indSvel)

        fc(1) = vcont*rhoa
        fc(2) = vcont*rhoua + pa*si(XCOORD,ijkN)
        fc(3) = vcont*rhova + pa*si(YCOORD,ijkN)
        fc(4) = vcont*rhowa + pa*si(ZCOORD,ijkN)
        fc(5) = vcont*(rhoea+pa) + siVel(ijkN*indSvel)*pa

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

        rhoa  = acj(2,ijkN)*cv(CV_MIXT_DENS,ijkC0)+ &
                acj(1,ijkN)*cv(CV_MIXT_DENS,ijkC1)
        rhoua = acj(2,ijkN)*cv(CV_MIXT_XMOM,ijkC0)+ &
                acj(1,ijkN)*cv(CV_MIXT_XMOM,ijkC1)
        rhova = acj(2,ijkN)*cv(CV_MIXT_YMOM,ijkC0)+ &
                acj(1,ijkN)*cv(CV_MIXT_YMOM,ijkC1)
        rhowa = acj(2,ijkN)*cv(CV_MIXT_ZMOM,ijkC0)+ &
                acj(1,ijkN)*cv(CV_MIXT_ZMOM,ijkC1)
        rhoea = acj(2,ijkN)*cv(CV_MIXT_ENER,ijkC0)+ &
                acj(1,ijkN)*cv(CV_MIXT_ENER,ijkC1)
        pa    = acj(2,ijkN)*dv(DV_MIXT_PRES,ijkC0)+ &
                acj(1,ijkN)*dv(DV_MIXT_PRES,ijkC1)
        vcont = (rhoua*sj(XCOORD,ijkN)+rhova*sj(YCOORD,ijkN)+&
                 rhowa*sj(ZCOORD,ijkN))/rhoa - sjVel(ijkN*indSvel)

        fc(1) = vcont*rhoa
        fc(2) = vcont*rhoua + pa*sj(XCOORD,ijkN)
        fc(3) = vcont*rhova + pa*sj(YCOORD,ijkN)
        fc(4) = vcont*rhowa + pa*sj(ZCOORD,ijkN)
        fc(5) = vcont*(rhoea+pa) + sjVel(ijkN*indSvel)*pa

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

        rhoa  = ack(2,ijkN)*cv(CV_MIXT_DENS,ijkC0)+ &
                ack(1,ijkN)*cv(CV_MIXT_DENS,ijkC1)
        rhoua = ack(2,ijkN)*cv(CV_MIXT_XMOM,ijkC0)+ &
                ack(1,ijkN)*cv(CV_MIXT_XMOM,ijkC1)
        rhova = ack(2,ijkN)*cv(CV_MIXT_YMOM,ijkC0)+ &
                ack(1,ijkN)*cv(CV_MIXT_YMOM,ijkC1)
        rhowa = ack(2,ijkN)*cv(CV_MIXT_ZMOM,ijkC0)+ &
                ack(1,ijkN)*cv(CV_MIXT_ZMOM,ijkC1)
        rhoea = ack(2,ijkN)*cv(CV_MIXT_ENER,ijkC0)+ &
                ack(1,ijkN)*cv(CV_MIXT_ENER,ijkC1)
        pa    = ack(2,ijkN)*dv(DV_MIXT_PRES,ijkC0)+ &
                ack(1,ijkN)*dv(DV_MIXT_PRES,ijkC1)
        vcont = (rhoua*sk(XCOORD,ijkN)+rhova*sk(YCOORD,ijkN)+&
                 rhowa*sk(ZCOORD,ijkN))/rhoa - skVel(ijkN*indSvel)

        fc(1) = vcont*rhoa
        fc(2) = vcont*rhoua + pa*sk(XCOORD,ijkN)
        fc(3) = vcont*rhova + pa*sk(YCOORD,ijkN)
        fc(4) = vcont*rhowa + pa*sk(ZCOORD,ijkN)
        fc(5) = vcont*(rhoea+pa) + skVel(ijkN*indSvel)*pa

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
    CALL RFLO_CentralFluxPatch( region,region%levels(iLev)%patches(iPatch) )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CentralFlux

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CentralFlux.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.14  2004/08/02 21:56:23  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.13  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.8  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.7  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/08/30 18:25:55  jblazek
! Forgot to multiply grid speeds by face area ...
!
! Revision 1.5  2002/08/29 23:25:54  jblazek
! Added support for moving grids.
!
! Revision 1.4  2002/03/18 23:11:32  jblazek
! Finished multiblock and MPI.
!
! Revision 1.3  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.1  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
!******************************************************************************







