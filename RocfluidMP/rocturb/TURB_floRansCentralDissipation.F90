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
! Purpose: compute artificial numerical dissipation for turbulence RaNS 
!          equation, if any.
!
! Description: the dissipation is 4th-order only, with no 2nd-order part
!
! Input: region = data of current region.
!
! Output: region%levels%turb%diss = RaNS dissipative fluxes.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_floRansCentralDissipation.F90,v 1.4 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansCentralDissipation( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i, j, k, ii, jj, kk, iCv

! ... local variables

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: nCv, iLev, iCOff, ijCOff, ijkC0, ijkCm1, ijkCp1, ijkCp2

  REAL(RFREAL) :: beta, eval, eps2, eps4, pmax, fd, vis2, vis4
  REAL(RFREAL), POINTER :: cv(:,:), diss(:,:), srad(:,:), dp(:), dv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_FloRansCentralDissipation',&
  'TURB_floRansCentralDissipation.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  nCv  =  region%turbInput%nCv
  cv   => region%levels(iLev)%turb%cv
  diss => region%levels(iLev)%turb%diss
  srad => region%levels(iLev)%turb%srad
  dp   => region%work1D

  dv   => region%levels(iLev)%mixt%dv

  beta =  region%mixtInput%betrk(region%irkStep)
  vis2 =  beta*region%turbInput%vis2
  vis4 =  beta*region%turbInput%vis4

! dissipation in i-direction ------------------------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend

! --- pressure switch

      ii = 0
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
        DO iCv = 1,nCv
          fd   = eps2*(  cv(iCv,ijkCp1) - cv(iCv,ijkC0)) + &
                 eps4*( (cv(iCv,ijkCm1) - cv(iCv,ijkCp2)) + &
              3._RFREAL*(cv(iCv,ijkCp1) - cv(iCv,ijkC0 )) )

          diss(iCv,ijkC0 ) = diss(iCv,ijkC0 ) + fd
          diss(iCv,ijkCp1) = diss(iCv,ijkCp1) - fd
        ENDDO ! iCv

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! dissipation in j-direction ------------------------------------------------

  DO k=kpcbeg,kpcend
    DO i=ipcbeg,ipcend

! --- pressure switch

      jj = 0
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
        DO iCv = 1,nCv
          fd   = eps2*(  cv(iCv,ijkCp1) - cv(iCv,ijkC0)) + &
                 eps4*( (cv(iCv,ijkCm1) - cv(iCv,ijkCp2)) + &
              3._RFREAL*(cv(iCv,ijkCp1) - cv(iCv,ijkC0 )) )

          diss(iCv,ijkC0 ) = diss(iCv,ijkC0 ) + fd
          diss(iCv,ijkCp1) = diss(iCv,ijkCp1) - fd
        ENDDO ! iCv

      ENDDO   ! j
    ENDDO     ! i
  ENDDO       ! k

! dissipation in k-direction --------------------------------------------------

  DO j=jpcbeg,jpcend
    DO i=ipcbeg,ipcend

! --- pressure switch

      kk = 0
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
        DO iCv = 1,nCv
          fd   = eps2*(  cv(iCv,ijkCp1) - cv(iCv,ijkC0)) + &
                 eps4*( (cv(iCv,ijkCm1) - cv(iCv,ijkCp2)) + &
              3._RFREAL*(cv(iCv,ijkCp1) - cv(iCv,ijkC0 )) )

          diss(iCv,ijkC0 ) = diss(iCv,ijkC0 ) + fd
          diss(iCv,ijkCp1) = diss(iCv,ijkCp1) - fd
        ENDDO ! iCv

      ENDDO   ! k
    ENDDO     ! i
  ENDDO       ! j

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloRansCentralDissipation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansCentralDissipation.F90,v $
! Revision 1.4  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.5  2003/10/27 04:51:37  wasistho
! added RaNS upwind schemes
!
! Revision 1.4  2003/10/15 22:06:07  wasistho
! use turb srad instead of mixture
!
! Revision 1.3  2003/10/15 03:41:21  wasistho
! added 2nd order dissipation coeff. k2
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







