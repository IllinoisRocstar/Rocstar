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
! Purpose: calculate coefficients for the implicit residual smoothing for 
!          RaNS class of turbulence
!
! Description: none.
!
! Input: region%levels%turb%srad = convective spectral radii for smoke.
!
! Output: region%levels%turb%epsIrs = smoothing coefficients.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_rFLO_RansResSmoothingCoeff.F90,v 1.3 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_RansResSmoothingCoeff( region ) ! PUBLIC

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

  global => region%global

  CALL RegisterFunction( global,'TURB_RFLO_RansResSmoothingCoeff',&
  'TURB_rFLO_RansResSmoothingCoeff.F90' )

  IF (region%turbInput%modelClass /= MODEL_RANS) GOTO 999

! begin -----------------------------------------------------------------------

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  srad   => region%levels(iLev)%turb%srad
  epsIrs => region%levels(iLev)%turb%epsIrs

  smoocf = region%turbInput%smoocf

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

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RFLO_RansResSmoothingCoeff

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_RansResSmoothingCoeff.F90,v $
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
! Revision 1.2  2003/10/20 00:37:06  wasistho
! goto 999 before register function
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







