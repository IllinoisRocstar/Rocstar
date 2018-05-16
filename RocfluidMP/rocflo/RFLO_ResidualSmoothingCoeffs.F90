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
! Purpose: calculate coefficients for the implicit residual smoothing.
!
! Description: none.
!
! Input: regions%levels%mixt%srad = convective spectral radii.
!
! Output: regions%levels%mixt%epsIrs = smoothing coefficients.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ResidualSmoothingCoeffs.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_residualSmoothingCoeffs( region )

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
  INTEGER :: iLev, iCOff, ijCOff, ijkC

  REAL(RFREAL) :: cflRat, smoocf, radij, radjk, radik, radji, radkj, radki
  REAL(RFREAL) :: psi, ex, ey, ez
  REAL(RFREAL), POINTER :: srad(:,:), epsIrs(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_residualSmoothingCoeffs',&
  'RFLO_ResidualSmoothingCoeffs.F90' )

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  srad   => region%levels(iLev)%mixt%srad
  epsIrs => region%levels(iLev)%mixt%epsIrs

  smoocf = region%mixtInput%smoocf

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

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_residualSmoothingCoeffs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ResidualSmoothingCoeffs.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.8  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
!******************************************************************************







