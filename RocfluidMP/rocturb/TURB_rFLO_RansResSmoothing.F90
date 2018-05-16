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
! Purpose: conduct implicit residual smoothing (central form) for RaNS class
!          of turbulence
!
! Description: none.
!
! Input: region = region data (residual, smoothing coeff., temporary array)
!
! Output: region%levels%turb%rhs = smoothed residual.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_rFLO_RansResSmoothing.F90,v 1.3 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_RansResSmoothing( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i, j, k, iCv

! ... local variables

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, off1d, ijkC, ijkC1, ijDum, ijDum1
  INTEGER :: im1, jm1, km1, ip1, jp1, kp1, nCv

  REAL(RFREAL)          :: t
  REAL(RFREAL), POINTER :: rhs(:,:), epsIrs(:,:), d(:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_RFLO_RansResSmoothing',&
  'TURB_rFLO_RansResSmoothing.F90' )

  IF (region%turbInput%modelClass /= MODEL_RANS) GOTO 999

! begin -----------------------------------------------------------------------

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  epsIrs => region%levels(iLev)%turb%epsIrs
  rhs    => region%levels(iLev)%turb%rhs
  d      => region%work1D

  nCv = region%turbInput%nCv

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

        DO iCv=1,nCv
          rhs(iCv,ijkC) = t*(rhs(iCv,ijkC) + epsIrs(ICOORD,ijkC)*rhs(iCv,ijkC1))
        ENDDO

      ENDDO
    ENDDO

    DO i=ipcend-1,ipcbeg,-1
      ip1 = i + 1
      DO j=jpcbeg,jpcend
        ijkC  = IndIJK(i  ,j,k,iCoff,ijCOff)
        ijkC1 = IndIJK(ip1,j,k,iCoff,ijCOff)
        ijDum = IndIJ(i,j,off1d) + 1

        DO iCv=1,nCv
          rhs(iCv,ijkC) = rhs(iCv,ijkC) + d(ijDum)*rhs(iCv,ijkC1)
        ENDDO

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

        DO iCv=1,nCv
          rhs(iCv,ijkC) = t*(rhs(iCv,ijkC) + epsIrs(ICOORD,ijkC)*rhs(iCv,ijkC1))
        ENDDO

      ENDDO
    ENDDO

    DO j=jpcend-1,jpcbeg,-1
      jp1 = j + 1
      DO i=ipcbeg,ipcend
        ijkC  = IndIJK(i,j  ,k,iCoff,ijCOff)
        ijkC1 = IndIJK(i,jp1,k,iCoff,ijCOff)
        ijDum = IndIJ(j,i,off1d) + 1

        DO iCv=1,nCv
          rhs(iCv,ijkC) = rhs(iCv,ijkC) + d(ijDum)*rhs(iCv,ijkC1)
        ENDDO

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

        DO iCv=1,nCv
          rhs(iCv,ijkC) = t*(rhs(iCv,ijkC) + epsIrs(ICOORD,ijkC)*rhs(iCv,ijkC1))
        ENDDO

      ENDDO
    ENDDO

    DO k=kpcend-1,kpcbeg,-1
      kp1 = k + 1
      DO i=ipcbeg,ipcend
        ijkC  = IndIJK(i,j,k  ,iCoff,ijCOff)
        ijkC1 = IndIJK(i,j,kp1,iCoff,ijCOff)
        ijDum = IndIJ(k,i,off1d) + 1

        DO iCv=1,nCv
          rhs(iCv,ijkC) = rhs(iCv,ijkC) + d(ijDum)*rhs(iCv,ijkC1)
        ENDDO

      ENDDO
    ENDDO

  ENDDO    ! j

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RFLO_RansResSmoothing

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_RansResSmoothing.F90,v $
! Revision 1.3  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/03/11 03:26:33  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







