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
! Purpose: compute central convective fluxes for smoke by average of variables.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%peul%rhs = convective fluxes added to the residual.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_CentralFlux.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_CentralFlux( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  USE ModInterfaces,      ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                                 RFLO_GetNodeOffset
  USE PEUL_ModInterfaces, ONLY : PEUL_CentralFluxPatch
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i,j,k,iPatch,ipt

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend,nPtypes
  INTEGER :: iLev,iCOff,ijCOff,iNOff,ijNOff,ijkC0,ijkC1,ijkN,indSvel

  REAL(RFREAL) :: sRhoa,gRhoa,gRhoua,gRhova,gRhowa,vcont,dS,sVel,fc
  REAL(RFREAL), POINTER :: sCv(:,:),gCv(:,:),sRhs(:,:)
  REAL(RFREAL), POINTER :: si(:,:),sj(:,:),sk(:,:)
  REAL(RFREAL), POINTER :: siVel(:),sjVel(:),skVel(:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_CentralFlux.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_CentralFlux',&
  'PEUL_CentralFlux.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  gCv    => region%levels(iLev)%mixt%cv
  sCv    => region%levels(iLev)%peul%cv
  sRhs   => region%levels(iLev)%peul%rhs
  si     => region%levels(iLev)%grid%si
  sj     => region%levels(iLev)%grid%sj
  sk     => region%levels(iLev)%grid%sk
  siVel  => region%levels(iLev)%grid%siVel
  sjVel  => region%levels(iLev)%grid%sjVel
  skVel  => region%levels(iLev)%grid%skVel
  indSvel = region%levels(iLev)%grid%indSvel

  nPtypes = region%peulInput%nPtypes
  IF (nPtypes /= region%levels(iLev)%peul%nCv) &
    CALL ErrorStop( global,ERR_PEUL_NPMISMATCH,__LINE__ )

! flux in i-direction (except through boundary) -------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg+1,ipcend
        ijkC0 = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkN  = IndIJK(i  ,j,k,iNOff,ijNOff)

        dS    = SQRT(si(XCOORD,ijkN)*si(XCOORD,ijkN)+ &
                     si(YCOORD,ijkN)*si(YCOORD,ijkN)+ &
                     si(ZCOORD,ijkN)*si(ZCOORD,ijkN))
        sVel  = siVel(ijkN*indSvel)*dS

        gRhoa  = 0.5_RFREAL*(gCv(CV_MIXT_DENS,ijkC0)+gCv(CV_MIXT_DENS,ijkC1))
        gRhoua = 0.5_RFREAL*(gCv(CV_MIXT_XMOM,ijkC0)+gCv(CV_MIXT_XMOM,ijkC1))
        gRhova = 0.5_RFREAL*(gCv(CV_MIXT_YMOM,ijkC0)+gCv(CV_MIXT_YMOM,ijkC1))
        gRhowa = 0.5_RFREAL*(gCv(CV_MIXT_ZMOM,ijkC0)+gCv(CV_MIXT_ZMOM,ijkC1))

        vcont = (gRhoua*si(XCOORD,ijkN)+gRhova*si(YCOORD,ijkN)+&
                 gRhowa*si(ZCOORD,ijkN))/gRhoa - sVel

        DO ipt=1,nPtypes
          sRhoa = 0.5_RFREAL*(sCv(ipt,ijkC0)+sCv(ipt,ijkC1))
          fc    = vcont*sRhoa
          sRhs(ipt,ijkC0) = sRhs(ipt,ijkC0) + fc
          sRhs(ipt,ijkC1) = sRhs(ipt,ijkC1) - fc
        ENDDO ! ipt

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

        dS    = SQRT(sj(XCOORD,ijkN)*sj(XCOORD,ijkN)+ &
                     sj(YCOORD,ijkN)*sj(YCOORD,ijkN)+ &
                     sj(ZCOORD,ijkN)*sj(ZCOORD,ijkN))
        sVel  = sjVel(ijkN*indSvel)*dS

        gRhoa  = 0.5_RFREAL*(gCv(CV_MIXT_DENS,ijkC0)+gCv(CV_MIXT_DENS,ijkC1))
        gRhoua = 0.5_RFREAL*(gCv(CV_MIXT_XMOM,ijkC0)+gCv(CV_MIXT_XMOM,ijkC1))
        gRhova = 0.5_RFREAL*(gCv(CV_MIXT_YMOM,ijkC0)+gCv(CV_MIXT_YMOM,ijkC1))
        gRhowa = 0.5_RFREAL*(gCv(CV_MIXT_ZMOM,ijkC0)+gCv(CV_MIXT_ZMOM,ijkC1))

        vcont = (gRhoua*sj(XCOORD,ijkN)+gRhova*sj(YCOORD,ijkN)+&
                 gRhowa*sj(ZCOORD,ijkN))/gRhoa - sVel

        DO ipt=1,nPtypes
          sRhoa = 0.5_RFREAL*(sCv(ipt,ijkC0)+sCv(ipt,ijkC1))
          fc    = vcont*sRhoa
          sRhs(ipt,ijkC0) = sRhs(ipt,ijkC0) + fc
          sRhs(ipt,ijkC1) = sRhs(ipt,ijkC1) - fc
        ENDDO ! ipt

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

        dS    = SQRT(sk(XCOORD,ijkN)*sk(XCOORD,ijkN)+ &
                     sk(YCOORD,ijkN)*sk(YCOORD,ijkN)+ &
                     sk(ZCOORD,ijkN)*sk(ZCOORD,ijkN))
        sVel  = skVel(ijkN*indSvel)*dS

        gRhoa  = 0.5_RFREAL*(gCv(CV_MIXT_DENS,ijkC0)+gCv(CV_MIXT_DENS,ijkC1))
        gRhoua = 0.5_RFREAL*(gCv(CV_MIXT_XMOM,ijkC0)+gCv(CV_MIXT_XMOM,ijkC1))
        gRhova = 0.5_RFREAL*(gCv(CV_MIXT_YMOM,ijkC0)+gCv(CV_MIXT_YMOM,ijkC1))
        gRhowa = 0.5_RFREAL*(gCv(CV_MIXT_ZMOM,ijkC0)+gCv(CV_MIXT_ZMOM,ijkC1))

        vcont = (gRhoua*sk(XCOORD,ijkN)+gRhova*sk(YCOORD,ijkN)+&
                 gRhowa*sk(ZCOORD,ijkN))/gRhoa - sVel

        DO ipt=1,nPtypes
          sRhoa = 0.5_RFREAL*(sCv(ipt,ijkC0)+sCv(ipt,ijkC1))
          fc    = vcont*sRhoa
          sRhs(ipt,ijkC0) = sRhs(ipt,ijkC0) + fc
          sRhs(ipt,ijkC1) = sRhs(ipt,ijkC1) - fc
        ENDDO ! ipt

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! fluxes through boundaries ---------------------------------------------------

  DO iPatch=1,region%nPatches
    CALL PEUL_CentralFluxPatch( region,region%levels(iLev)%patches(iPatch) )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_CentralFlux

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_CentralFlux.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:26  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/07/28 15:42:13  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.6  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2003/09/25 15:42:57  jferry
! Added mixture source terms due to active smoke (for Detangle interaction)
!
! Revision 1.4  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.3  2003/05/01 22:57:24  jferry
! substituted macro for IndIJK
!
! Revision 1.2  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







