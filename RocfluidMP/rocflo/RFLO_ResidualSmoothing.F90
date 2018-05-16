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
! Purpose: conduct implicit residual smoothing (central form).
!
! Description: none.
!
! Input: region = region data (residual, smoothing coeff., temporary array)
!
! Output: region%levels%mixt%rhs = smoothed residual.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ResidualSmoothing.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ResidualSmoothing( region )

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
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, off1d, ijkC, ijkC1, ijDum, ijDum1
  INTEGER :: im1, jm1, km1, ip1, jp1, kp1

  REAL(RFREAL)          :: t
  REAL(RFREAL), POINTER :: rhs(:,:), epsIrs(:,:), d(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ResidualSmoothing',&
  'RFLO_ResidualSmoothing.F90' )

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  epsIrs => region%levels(iLev)%mixt%epsIrs
  rhs    => region%levels(iLev)%mixt%rhs
  d      => region%work1D

! solution of tridiagonal system in I-direction

  DO k=kpcbeg,kpcend

    i     = ipcbeg - 1
    off1d = ipcend + 1
    DO j=jpcbeg,jpcend
      ijDum    = IndIJ(i,j,off1d) + 1
      d(ijDum) = 0._RFREAL
    ENDDO

    DO i=ipcbeg,ipcend
      im1 = i - 1
      DO j=jpcbeg,jpcend
        ijkC     = IndIJK(i  ,j,k,iCoff,ijCOff)
        ijkC1    = IndIJK(im1,j,k,iCoff,ijCOff)
        ijDum    = IndIJ(i  ,j,off1d) + 1
        ijDum1   = IndIJ(im1,j,off1d) + 1
        t        = 1._RFREAL/(1._RFREAL+2._RFREAL*epsIrs(ICOORD,ijkC)- &
                              epsIrs(ICOORD,ijkC)*d(ijDum1))
        d(ijDum) = t*epsIrs(ICOORD,ijkC)

        rhs(CV_MIXT_DENS,ijkC) = t*(rhs(CV_MIXT_DENS,ijkC) + &
                                    epsIrs(ICOORD,ijkC)*rhs(CV_MIXT_DENS,ijkC1))
        rhs(CV_MIXT_XMOM,ijkC) = t*(rhs(CV_MIXT_XMOM,ijkC) + &
                                    epsIrs(ICOORD,ijkC)*rhs(CV_MIXT_XMOM,ijkC1))
        rhs(CV_MIXT_YMOM,ijkC) = t*(rhs(CV_MIXT_YMOM,ijkC) + &
                                    epsIrs(ICOORD,ijkC)*rhs(CV_MIXT_YMOM,ijkC1))
        rhs(CV_MIXT_ZMOM,ijkC) = t*(rhs(CV_MIXT_ZMOM,ijkC) + &
                                    epsIrs(ICOORD,ijkC)*rhs(CV_MIXT_ZMOM,ijkC1))
        rhs(CV_MIXT_ENER,ijkC) = t*(rhs(CV_MIXT_ENER,ijkC) + &
                                    epsIrs(ICOORD,ijkC)*rhs(CV_MIXT_ENER,ijkC1))
      ENDDO
    ENDDO

    DO i=ipcend-1,ipcbeg,-1
      ip1 = i + 1
      DO j=jpcbeg,jpcend
        ijkC  = IndIJK(i  ,j,k,iCoff,ijCOff)
        ijkC1 = IndIJK(ip1,j,k,iCoff,ijCOff)
        ijDum = IndIJ(i,j,off1d) + 1

        rhs(CV_MIXT_DENS,ijkC) = rhs(CV_MIXT_DENS,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_DENS,ijkC1)
        rhs(CV_MIXT_XMOM,ijkC) = rhs(CV_MIXT_XMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_XMOM,ijkC1)
        rhs(CV_MIXT_YMOM,ijkC) = rhs(CV_MIXT_YMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_YMOM,ijkC1)
        rhs(CV_MIXT_ZMOM,ijkC) = rhs(CV_MIXT_ZMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_ZMOM,ijkC1)
        rhs(CV_MIXT_ENER,ijkC) = rhs(CV_MIXT_ENER,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_ENER,ijkC1)
      ENDDO
    ENDDO

  ENDDO    ! k

! solution of tridiagonal system in J-direction

  DO k=kpcbeg,kpcend

    j     = jpcbeg - 1
    off1d = jpcend + 1
    DO i=ipcbeg,ipcend
      ijDum    = IndIJ(j,i,off1d) + 1
      d(ijDum) = 0._RFREAL
    ENDDO

    DO j=jpcbeg,jpcend
      jm1 = j - 1
      DO i=ipcbeg,ipcend
        ijkC     = IndIJK(i,j  ,k,iCoff,ijCOff)
        ijkC1    = IndIJK(i,jm1,k,iCoff,ijCOff)
        ijDum    = IndIJ(j  ,i,off1d) + 1
        ijDum1   = IndIJ(jm1,i,off1d) + 1
        t        = 1._RFREAL/(1._RFREAL+2._RFREAL*epsIrs(JCOORD,ijkC)- &
                              epsIrs(JCOORD,ijkC)*d(ijDum1))
        d(ijDum) = t*epsIrs(JCOORD,ijkC)

        rhs(CV_MIXT_DENS,ijkC) = t*(rhs(CV_MIXT_DENS,ijkC) + &
                                    epsIrs(JCOORD,ijkC)*rhs(CV_MIXT_DENS,ijkC1))
        rhs(CV_MIXT_XMOM,ijkC) = t*(rhs(CV_MIXT_XMOM,ijkC) + &
                                    epsIrs(JCOORD,ijkC)*rhs(CV_MIXT_XMOM,ijkC1))
        rhs(CV_MIXT_YMOM,ijkC) = t*(rhs(CV_MIXT_YMOM,ijkC) + &
                                    epsIrs(JCOORD,ijkC)*rhs(CV_MIXT_YMOM,ijkC1))
        rhs(CV_MIXT_ZMOM,ijkC) = t*(rhs(CV_MIXT_ZMOM,ijkC) + &
                                    epsIrs(JCOORD,ijkC)*rhs(CV_MIXT_ZMOM,ijkC1))
        rhs(CV_MIXT_ENER,ijkC) = t*(rhs(CV_MIXT_ENER,ijkC) + &
                                    epsIrs(JCOORD,ijkC)*rhs(CV_MIXT_ENER,ijkC1))
      ENDDO
    ENDDO

    DO j=jpcend-1,jpcbeg,-1
      jp1 = j + 1
      DO i=ipcbeg,ipcend
        ijkC  = IndIJK(i,j  ,k,iCoff,ijCOff)
        ijkC1 = IndIJK(i,jp1,k,iCoff,ijCOff)
        ijDum = IndIJ(j,i,off1d) + 1

        rhs(CV_MIXT_DENS,ijkC) = rhs(CV_MIXT_DENS,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_DENS,ijkC1)
        rhs(CV_MIXT_XMOM,ijkC) = rhs(CV_MIXT_XMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_XMOM,ijkC1)
        rhs(CV_MIXT_YMOM,ijkC) = rhs(CV_MIXT_YMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_YMOM,ijkC1)
        rhs(CV_MIXT_ZMOM,ijkC) = rhs(CV_MIXT_ZMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_ZMOM,ijkC1)
        rhs(CV_MIXT_ENER,ijkC) = rhs(CV_MIXT_ENER,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_ENER,ijkC1)
      ENDDO
    ENDDO

  ENDDO    ! k

! solution of tridiagonal system in K-direction

  DO j=jpcbeg,jpcend

    k     = kpcbeg - 1
    off1d = kpcend + 1
    DO i=ipcbeg,ipcend
      ijDum    = IndIJ(k,i,off1d) + 1
      d(ijDum) = 0._RFREAL
    ENDDO

    DO k=kpcbeg,kpcend
      km1 = k - 1
      DO i=ipcbeg,ipcend
        ijkC     = IndIJK(i,j,k  ,iCoff,ijCOff)
        ijkC1    = IndIJK(i,j,km1,iCoff,ijCOff)
        ijDum    = IndIJ(k  ,i,off1d) + 1
        ijDum1   = IndIJ(km1,i,off1d) + 1
        t        = 1._RFREAL/(1._RFREAL+2._RFREAL*epsIrs(KCOORD,ijkC)- &
                              epsIrs(KCOORD,ijkC)*d(ijDum1))
        d(ijDum) = t*epsIrs(KCOORD,ijkC)

        rhs(CV_MIXT_DENS,ijkC) = t*(rhs(CV_MIXT_DENS,ijkC) + &
                                    epsIrs(KCOORD,ijkC)*rhs(CV_MIXT_DENS,ijkC1))
        rhs(CV_MIXT_XMOM,ijkC) = t*(rhs(CV_MIXT_XMOM,ijkC) + &
                                    epsIrs(KCOORD,ijkC)*rhs(CV_MIXT_XMOM,ijkC1))
        rhs(CV_MIXT_YMOM,ijkC) = t*(rhs(CV_MIXT_YMOM,ijkC) + &
                                    epsIrs(KCOORD,ijkC)*rhs(CV_MIXT_YMOM,ijkC1))
        rhs(CV_MIXT_ZMOM,ijkC) = t*(rhs(CV_MIXT_ZMOM,ijkC) + &
                                    epsIrs(KCOORD,ijkC)*rhs(CV_MIXT_ZMOM,ijkC1))
        rhs(CV_MIXT_ENER,ijkC) = t*(rhs(CV_MIXT_ENER,ijkC) + &
                                    epsIrs(KCOORD,ijkC)*rhs(CV_MIXT_ENER,ijkC1))
      ENDDO
    ENDDO

    DO k=kpcend-1,kpcbeg,-1
      kp1 = k + 1
      DO i=ipcbeg,ipcend
        ijkC  = IndIJK(i,j,k  ,iCoff,ijCOff)
        ijkC1 = IndIJK(i,j,kp1,iCoff,ijCOff)
        ijDum = IndIJ(k,i,off1d) + 1

        rhs(CV_MIXT_DENS,ijkC) = rhs(CV_MIXT_DENS,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_DENS,ijkC1)
        rhs(CV_MIXT_XMOM,ijkC) = rhs(CV_MIXT_XMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_XMOM,ijkC1)
        rhs(CV_MIXT_YMOM,ijkC) = rhs(CV_MIXT_YMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_YMOM,ijkC1)
        rhs(CV_MIXT_ZMOM,ijkC) = rhs(CV_MIXT_ZMOM,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_ZMOM,ijkC1)
        rhs(CV_MIXT_ENER,ijkC) = rhs(CV_MIXT_ENER,ijkC) + &
                                 d(ijDum)*rhs(CV_MIXT_ENER,ijkC1)
      ENDDO
    ENDDO

  ENDDO    ! j

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ResidualSmoothing

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ResidualSmoothing.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.10  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.5  2003/01/10 17:58:43  jblazek
! Added missing explicit interfaces.
!
! Revision 1.4  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.1  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
!******************************************************************************







