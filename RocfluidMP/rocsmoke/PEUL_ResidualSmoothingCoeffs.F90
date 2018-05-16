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
! Purpose: calculate coefficients for the implicit residual smoothing for smoke
!
! Description: none.
!
! Input: regions%levels%peul%srad = convective spectral radii for smoke.
!
! Output: regions%levels%peul%epsIrs = smoothing coefficients.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ResidualSmoothingCoeffs.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_residualSmoothingCoeffs( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, ijkC

  REAL(RFREAL) :: cflRat, smoocf, radij, radjk, radik, radji, radkj, radki
  REAL(RFREAL) :: psi, ex, ey, ez
  REAL(RFREAL), POINTER :: srad(:,:), epsIrs(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_ResidualSmoothingCoeffs.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_residualSmoothingCoeffs',&
  'PEUL_ResidualSmoothingCoeffs.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  srad   => region%levels(iLev)%peul%srad
  epsIrs => region%levels(iLev)%peul%epsIrs

  smoocf = region%peulInput%smoocf

! loop over physical cells

  cflRat = SQRT(1._RFREAL+4._RFREAL*smoocf)
  psi    = 0.0625_RFREAL

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        ijkC  = IndIJK(i,j,k,iCoff,ijCOff)
        radij = srad(ICOORD,ijkC)/srad(JCOORD,ijkC)
        radjk = srad(JCOORD,ijkC)/srad(KCOORD,ijkC)
        radik = srad(ICOORD,ijkC)/srad(KCOORD,ijkC)
        radji = 1._RFREAL/radij
        radkj = 1._RFREAL/radjk
        radki = 1._RFREAL/radik

        ex = 0.25_RFREAL*((cflRat/(1._RFREAL+psi*(radji+radki)))**2-1._RFREAL)
        ey = 0.25_RFREAL*((cflRat/(1._RFREAL+psi*(radij+radkj)))**2-1._RFREAL)
        ez = 0.25_RFREAL*((cflRat/(1._RFREAL+psi*(radik+radjk)))**2-1._RFREAL)
        ex = MIN(smoocf,ex)
        ey = MIN(smoocf,ey)
        ez = MIN(smoocf,ez)

        epsIrs(ICOORD,ijkC) = MAX(0._RFREAL,ex)
        epsIrs(JCOORD,ijkC) = MAX(0._RFREAL,ey)
        epsIrs(KCOORD,ijkC) = MAX(0._RFREAL,ez)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_residualSmoothingCoeffs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ResidualSmoothingCoeffs.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:53  haselbac
! Initial revision after changing case
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/09 15:12:04  jferry
! miscellaneous stylistic changes
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







