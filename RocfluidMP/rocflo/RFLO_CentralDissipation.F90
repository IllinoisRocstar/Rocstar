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
! Purpose: compute artificial numerical dissipation (JST type).
!
! Description: the dissipation consists of a blend of 2nd- and 4th-order
!              differences (undivided).
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%diss = dissipative fluxes.
!
! Notes: the TVD type of pressure switch is taken from:
!        Turkel, E.; Swanson, R.C.; Vatsa, V.N.; White, J.A.: Multigrid
!        for Hypersonic Viscous Two- and Three-Dimensional Flows. AIAA
!        Paper 91-1572, 1991.
!
!******************************************************************************
!
! $Id: RFLO_CentralDissipation.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CentralDissipation( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, ii, jj, kk

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, ijkC0, ijkCm1, ijkCp1, ijkCp2, pSwitchType

  REAL(RFREAL)          :: beta, vis2, vis4, eval, pmax, eps2, eps4, fd(5)
  REAL(RFREAL)          :: pSwitchOmega, pTvd, pSum
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), diss(:,:), srad(:,:), dp(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CentralDissipation',&
  'RFLO_CentralDissipation.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv   => region%levels(iLev)%mixt%cv
  dv   => region%levels(iLev)%mixt%dv
  diss => region%levels(iLev)%mixt%diss
  srad => region%levels(iLev)%mixt%srad
  dp   => region%work1D

! get coefficients and switch type --------------------------------------------

  beta         = region%mixtInput%betrk(region%irkStep)
  vis2         = beta*region%mixtInput%vis2
  vis4         = beta*region%mixtInput%vis4
  pSwitchType  = region%mixtInput%pSwitchType
  pSwitchOmega = region%mixtInput%pSwitchOmega

! dissipation in i-direction --------------------------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend

! --- pressure switch

      ii = 0
      IF (pSwitchType == PSWITCH_STD) THEN
        DO i=ipcbeg-1,ipcend+1
          ii     = ii + 1
          ijkC0  = IndIJK(i  ,j,k,iCOff,ijCOff)
          ijkCm1 = IndIJK(i-1,j,k,iCOff,ijCOff)
          ijkCp1 = IndIJK(i+1,j,k,iCOff,ijCOff)
          dp(ii) = ABS((          dv(DV_MIXT_PRES,ijkCp1)- &
                        2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                  dv(DV_MIXT_PRES,ijkCm1))/ &
                       (          dv(DV_MIXT_PRES,ijkCp1)+ &
                        2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                  dv(DV_MIXT_PRES,ijkCm1)))
        ENDDO
      ELSE
        DO i=ipcbeg-1,ipcend+1
          ii     = ii + 1
          ijkC0  = IndIJK(i  ,j,k,iCOff,ijCOff)
          ijkCm1 = IndIJK(i-1,j,k,iCOff,ijCOff)
          ijkCp1 = IndIJK(i+1,j,k,iCOff,ijCOff)
          pSum   =           dv(DV_MIXT_PRES,ijkCp1)+ &
                   2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                             dv(DV_MIXT_PRES,ijkCm1)
          pTvd   = ABS(dv(DV_MIXT_PRES,ijkCp1)-dv(DV_MIXT_PRES,ijkC0 ))+ &
                   ABS(dv(DV_MIXT_PRES,ijkC0 )-dv(DV_MIXT_PRES,ijkCm1))
          dp(ii) = ABS(          dv(DV_MIXT_PRES,ijkCp1)- &
                       2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                 dv(DV_MIXT_PRES,ijkCm1))/ &
                      ((1._RFREAL-pSwitchOmega)*pTvd+pSwitchOmega*pSum)
        ENDDO
      ENDIF

! --- dissipative fluxes at I+1/2

      ii = 0
      DO i=ipcbeg-1,ipcend
        ii     = ii + 1
        ijkC0  = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkCm1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkCp1 = IndIJK(i+1,j,k,iCOff,ijCOff)
        ijkCp2 = IndIJK(i+2,j,k,iCOff,ijCOff)
        eval   = 0.5_RFREAL*(srad(ICOORD,ijkC0)+srad(ICOORD,ijkCp1) + &
                             MAX(srad(JCOORD,ijkC0)+srad(JCOORD,ijkCp1), &
                                 srad(KCOORD,ijkC0)+srad(KCOORD,ijkCp1)))
        pmax   = MAX(dp(ii),dp(ii+1))
        eps2   = eval*vis2*pmax
        eps4   = eval*vis4
        eps4   = DIM(eps4,eps2)
        fd(1)  = eps2*(cv(CV_MIXT_DENS,ijkCp1)-cv(CV_MIXT_DENS,ijkC0)) - &
                 eps4*(cv(CV_MIXT_DENS,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_DENS,ijkCp1)- &
                                  cv(CV_MIXT_DENS,ijkC0 ))- &
                       cv(CV_MIXT_DENS,ijkCm1))
        fd(2)  = eps2*(cv(CV_MIXT_XMOM,ijkCp1)-cv(CV_MIXT_XMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_XMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_XMOM,ijkCp1)- &
                                  cv(CV_MIXT_XMOM,ijkC0 ))- &
                       cv(CV_MIXT_XMOM,ijkCm1))
        fd(3)  = eps2*(cv(CV_MIXT_YMOM,ijkCp1)-cv(CV_MIXT_YMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_YMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_YMOM,ijkCp1)- &
                                  cv(CV_MIXT_YMOM,ijkC0 ))- &
                       cv(CV_MIXT_YMOM,ijkCm1))
        fd(4)  = eps2*(cv(CV_MIXT_ZMOM,ijkCp1)-cv(CV_MIXT_ZMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_ZMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_ZMOM,ijkCp1)- &
                                  cv(CV_MIXT_ZMOM,ijkC0 ))- &
                       cv(CV_MIXT_ZMOM,ijkCm1))
        fd(5)  = eps2*(cv(CV_MIXT_ENER,ijkCp1)-cv(CV_MIXT_ENER,ijkC0)) - &
                 eps4*(cv(CV_MIXT_ENER,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_ENER,ijkCp1)- &
                                  cv(CV_MIXT_ENER,ijkC0 ))- &
                       cv(CV_MIXT_ENER,ijkCm1))

        diss(CV_MIXT_DENS,ijkC0 ) = diss(CV_MIXT_DENS,ijkC0 ) + fd(1)
        diss(CV_MIXT_XMOM,ijkC0 ) = diss(CV_MIXT_XMOM,ijkC0 ) + fd(2)
        diss(CV_MIXT_YMOM,ijkC0 ) = diss(CV_MIXT_YMOM,ijkC0 ) + fd(3)
        diss(CV_MIXT_ZMOM,ijkC0 ) = diss(CV_MIXT_ZMOM,ijkC0 ) + fd(4)
        diss(CV_MIXT_ENER,ijkC0 ) = diss(CV_MIXT_ENER,ijkC0 ) + fd(5)

        diss(CV_MIXT_DENS,ijkCp1) = diss(CV_MIXT_DENS,ijkCp1) - fd(1)
        diss(CV_MIXT_XMOM,ijkCp1) = diss(CV_MIXT_XMOM,ijkCp1) - fd(2)
        diss(CV_MIXT_YMOM,ijkCp1) = diss(CV_MIXT_YMOM,ijkCp1) - fd(3)
        diss(CV_MIXT_ZMOM,ijkCp1) = diss(CV_MIXT_ZMOM,ijkCp1) - fd(4)
        diss(CV_MIXT_ENER,ijkCp1) = diss(CV_MIXT_ENER,ijkCp1) - fd(5)
      ENDDO  ! i

    ENDDO    ! j
  ENDDO      ! k

! dissipation in j-direction --------------------------------------------------

  DO k=kpcbeg,kpcend
    DO i=ipcbeg,ipcend

! --- pressure switch

      jj = 0
      IF (pSwitchType == PSWITCH_STD) THEN
        DO j=jpcbeg-1,jpcend+1
          jj     = jj + 1
          ijkC0  = IndIJK(i,j  ,k,iCOff,ijCOff)
          ijkCm1 = IndIJK(i,j-1,k,iCOff,ijCOff)
          ijkCp1 = IndIJK(i,j+1,k,iCOff,ijCOff)
          dp(jj) = ABS((          dv(DV_MIXT_PRES,ijkCp1)- &
                        2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                  dv(DV_MIXT_PRES,ijkCm1))/ &
                       (          dv(DV_MIXT_PRES,ijkCp1)+ &
                        2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                  dv(DV_MIXT_PRES,ijkCm1)))
        ENDDO
      ELSE
        DO j=jpcbeg-1,jpcend+1
          jj     = jj + 1
          ijkC0  = IndIJK(i,j  ,k,iCOff,ijCOff)
          ijkCm1 = IndIJK(i,j-1,k,iCOff,ijCOff)
          ijkCp1 = IndIJK(i,j+1,k,iCOff,ijCOff)
          pSum   =           dv(DV_MIXT_PRES,ijkCp1)+ &
                   2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                             dv(DV_MIXT_PRES,ijkCm1)
          pTvd   = ABS(dv(DV_MIXT_PRES,ijkCp1)-dv(DV_MIXT_PRES,ijkC0 ))+ &
                   ABS(dv(DV_MIXT_PRES,ijkC0 )-dv(DV_MIXT_PRES,ijkCm1))
          dp(jj) = ABS(          dv(DV_MIXT_PRES,ijkCp1)- &
                       2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                 dv(DV_MIXT_PRES,ijkCm1))/ &
                      ((1._RFREAL-pSwitchOmega)*pTvd+pSwitchOmega*pSum)
        ENDDO
      ENDIF

! --- dissipative fluxes at J+1/2

      jj = 0
      DO j=jpcbeg-1,jpcend
        jj     = jj + 1
        ijkC0  = IndIJK(i,j  ,k,iCOff,ijCOff)
        ijkCm1 = IndIJK(i,j-1,k,iCOff,ijCOff)
        ijkCp1 = IndIJK(i,j+1,k,iCOff,ijCOff)
        ijkCp2 = IndIJK(i,j+2,k,iCOff,ijCOff)
        eval   = 0.5_RFREAL*(srad(JCOORD,ijkC0)+srad(JCOORD,ijkCp1) + &
                             MAX(srad(ICOORD,ijkC0)+srad(ICOORD,ijkCp1), &
                                 srad(KCOORD,ijkC0)+srad(KCOORD,ijkCp1)))
        pmax   = MAX(dp(jj),dp(jj+1))
        eps2   = eval*vis2*pmax
        eps4   = eval*vis4
        eps4   = DIM(eps4,eps2)
        fd(1)  = eps2*(cv(CV_MIXT_DENS,ijkCp1)-cv(CV_MIXT_DENS,ijkC0)) - &
                 eps4*(cv(CV_MIXT_DENS,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_DENS,ijkCp1)- &
                                  cv(CV_MIXT_DENS,ijkC0 ))- &
                       cv(CV_MIXT_DENS,ijkCm1))
        fd(2)  = eps2*(cv(CV_MIXT_XMOM,ijkCp1)-cv(CV_MIXT_XMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_XMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_XMOM,ijkCp1)- &
                                  cv(CV_MIXT_XMOM,ijkC0 ))- &
                       cv(CV_MIXT_XMOM,ijkCm1))
        fd(3)  = eps2*(cv(CV_MIXT_YMOM,ijkCp1)-cv(CV_MIXT_YMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_YMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_YMOM,ijkCp1)- &
                                  cv(CV_MIXT_YMOM,ijkC0 ))- &
                       cv(CV_MIXT_YMOM,ijkCm1))
        fd(4)  = eps2*(cv(CV_MIXT_ZMOM,ijkCp1)-cv(CV_MIXT_ZMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_ZMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_ZMOM,ijkCp1)- &
                                  cv(CV_MIXT_ZMOM,ijkC0 ))- &
                       cv(CV_MIXT_ZMOM,ijkCm1))
        fd(5)  = eps2*(cv(CV_MIXT_ENER,ijkCp1)-cv(CV_MIXT_ENER,ijkC0)) - &
                 eps4*(cv(CV_MIXT_ENER,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_ENER,ijkCp1)- &
                                  cv(CV_MIXT_ENER,ijkC0 ))- &
                       cv(CV_MIXT_ENER,ijkCm1))

        diss(CV_MIXT_DENS,ijkC0 ) = diss(CV_MIXT_DENS,ijkC0 ) + fd(1)
        diss(CV_MIXT_XMOM,ijkC0 ) = diss(CV_MIXT_XMOM,ijkC0 ) + fd(2)
        diss(CV_MIXT_YMOM,ijkC0 ) = diss(CV_MIXT_YMOM,ijkC0 ) + fd(3)
        diss(CV_MIXT_ZMOM,ijkC0 ) = diss(CV_MIXT_ZMOM,ijkC0 ) + fd(4)
        diss(CV_MIXT_ENER,ijkC0 ) = diss(CV_MIXT_ENER,ijkC0 ) + fd(5)

        diss(CV_MIXT_DENS,ijkCp1) = diss(CV_MIXT_DENS,ijkCp1) - fd(1)
        diss(CV_MIXT_XMOM,ijkCp1) = diss(CV_MIXT_XMOM,ijkCp1) - fd(2)
        diss(CV_MIXT_YMOM,ijkCp1) = diss(CV_MIXT_YMOM,ijkCp1) - fd(3)
        diss(CV_MIXT_ZMOM,ijkCp1) = diss(CV_MIXT_ZMOM,ijkCp1) - fd(4)
        diss(CV_MIXT_ENER,ijkCp1) = diss(CV_MIXT_ENER,ijkCp1) - fd(5)
      ENDDO  ! j

    ENDDO    ! i
  ENDDO      ! k

! dissipation in k-direction --------------------------------------------------

  DO j=jpcbeg,jpcend
    DO i=ipcbeg,ipcend

! --- pressure switch

      kk = 0
      IF (pSwitchType == PSWITCH_STD) THEN
        DO k=kpcbeg-1,kpcend+1
          kk     = kk + 1
          ijkC0  = IndIJK(i,j,k  ,iCOff,ijCOff)
          ijkCm1 = IndIJK(i,j,k-1,iCOff,ijCOff)
          ijkCp1 = IndIJK(i,j,k+1,iCOff,ijCOff)
          dp(kk) = ABS((          dv(DV_MIXT_PRES,ijkCp1)- &
                        2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                  dv(DV_MIXT_PRES,ijkCm1))/ &
                       (          dv(DV_MIXT_PRES,ijkCp1)+ &
                        2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                  dv(DV_MIXT_PRES,ijkCm1)))
        ENDDO
      ELSE
        DO k=kpcbeg-1,kpcend+1
          kk     = kk + 1
          ijkC0  = IndIJK(i,j,k  ,iCOff,ijCOff)
          ijkCm1 = IndIJK(i,j,k-1,iCOff,ijCOff)
          ijkCp1 = IndIJK(i,j,k+1,iCOff,ijCOff)
          pSum   =           dv(DV_MIXT_PRES,ijkCp1)+ &
                   2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                             dv(DV_MIXT_PRES,ijkCm1)
          pTvd   = ABS(dv(DV_MIXT_PRES,ijkCp1)-dv(DV_MIXT_PRES,ijkC0 ))+ &
                   ABS(dv(DV_MIXT_PRES,ijkC0 )-dv(DV_MIXT_PRES,ijkCm1))
          dp(kk) = ABS(          dv(DV_MIXT_PRES,ijkCp1)- &
                       2._RFREAL*dv(DV_MIXT_PRES,ijkC0 )+ &
                                 dv(DV_MIXT_PRES,ijkCm1))/ &
                      ((1._RFREAL-pSwitchOmega)*pTvd+pSwitchOmega*pSum)
        ENDDO
      ENDIF

! --- dissipative fluxes at K+1/2

      kk = 0
      DO k=kpcbeg-1,kpcend
        kk     = kk + 1
        ijkC0  = IndIJK(i,j,k  ,iCOff,ijCOff)
        ijkCm1 = IndIJK(i,j,k-1,iCOff,ijCOff)
        ijkCp1 = IndIJK(i,j,k+1,iCOff,ijCOff)
        ijkCp2 = IndIJK(i,j,k+2,iCOff,ijCOff)
        eval   = 0.5_RFREAL*(srad(KCOORD,ijkC0)+srad(KCOORD,ijkCp1) + &
                             MAX(srad(ICOORD,ijkC0)+srad(ICOORD,ijkCp1), &
                                 srad(JCOORD,ijkC0)+srad(JCOORD,ijkCp1)))
        pmax   = MAX(dp(kk),dp(kk+1))
        eps2   = eval*vis2*pmax
        eps4   = eval*vis4
        eps4   = DIM(eps4,eps2)
        fd(1)  = eps2*(cv(CV_MIXT_DENS,ijkCp1)-cv(CV_MIXT_DENS,ijkC0)) - &
                 eps4*(cv(CV_MIXT_DENS,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_DENS,ijkCp1)- &
                                  cv(CV_MIXT_DENS,ijkC0 ))- &
                       cv(CV_MIXT_DENS,ijkCm1))
        fd(2)  = eps2*(cv(CV_MIXT_XMOM,ijkCp1)-cv(CV_MIXT_XMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_XMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_XMOM,ijkCp1)- &
                                  cv(CV_MIXT_XMOM,ijkC0 ))- &
                       cv(CV_MIXT_XMOM,ijkCm1))
        fd(3)  = eps2*(cv(CV_MIXT_YMOM,ijkCp1)-cv(CV_MIXT_YMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_YMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_YMOM,ijkCp1)- &
                                  cv(CV_MIXT_YMOM,ijkC0 ))- &
                       cv(CV_MIXT_YMOM,ijkCm1))
        fd(4)  = eps2*(cv(CV_MIXT_ZMOM,ijkCp1)-cv(CV_MIXT_ZMOM,ijkC0)) - &
                 eps4*(cv(CV_MIXT_ZMOM,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_ZMOM,ijkCp1)- &
                                  cv(CV_MIXT_ZMOM,ijkC0 ))- &
                       cv(CV_MIXT_ZMOM,ijkCm1))
        fd(5)  = eps2*(cv(CV_MIXT_ENER,ijkCp1)-cv(CV_MIXT_ENER,ijkC0)) - &
                 eps4*(cv(CV_MIXT_ENER,ijkCp2)- &
                       3._RFREAL*(cv(CV_MIXT_ENER,ijkCp1)- &
                                  cv(CV_MIXT_ENER,ijkC0 ))- &
                       cv(CV_MIXT_ENER,ijkCm1))

        diss(CV_MIXT_DENS,ijkC0 ) = diss(CV_MIXT_DENS,ijkC0 ) + fd(1)
        diss(CV_MIXT_XMOM,ijkC0 ) = diss(CV_MIXT_XMOM,ijkC0 ) + fd(2)
        diss(CV_MIXT_YMOM,ijkC0 ) = diss(CV_MIXT_YMOM,ijkC0 ) + fd(3)
        diss(CV_MIXT_ZMOM,ijkC0 ) = diss(CV_MIXT_ZMOM,ijkC0 ) + fd(4)
        diss(CV_MIXT_ENER,ijkC0 ) = diss(CV_MIXT_ENER,ijkC0 ) + fd(5)

        diss(CV_MIXT_DENS,ijkCp1) = diss(CV_MIXT_DENS,ijkCp1) - fd(1)
        diss(CV_MIXT_XMOM,ijkCp1) = diss(CV_MIXT_XMOM,ijkCp1) - fd(2)
        diss(CV_MIXT_YMOM,ijkCp1) = diss(CV_MIXT_YMOM,ijkCp1) - fd(3)
        diss(CV_MIXT_ZMOM,ijkCp1) = diss(CV_MIXT_ZMOM,ijkCp1) - fd(4)
        diss(CV_MIXT_ENER,ijkCp1) = diss(CV_MIXT_ENER,ijkCp1) - fd(5)
      ENDDO  ! k

    ENDDO    ! i
  ENDDO      ! j

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CentralDissipation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CentralDissipation.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.10  2003/11/20 16:40:38  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/07/25 00:36:48  jblazek
! Option for TVD type pressure switch.
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







