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
!          FLD radiation model.
!
! Description: none.
!
! Input: region%levels%radi%srad = convective spectral radii for smoke.
!
! Output: region%levels%radi%epsIrs = smoothing coefficients.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_rFLO_FlimResSmoothingCoeff.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_RFLO_FlimResSmoothingCoeff( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
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

  CALL RegisterFunction( global,'RADI_RFLO_FlimResSmoothingCoeff',&
  'RADI_rFLO_FlimResSmoothingCoeff.F90' )

  IF (region%radiInput%radiModel /= RADI_MODEL_FLDTRAN) GOTO 999

! begin -----------------------------------------------------------------------

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  srad   => region%levels(iLev)%radi%srad
  epsIrs => region%levels(iLev)%radi%epsIrs

  smoocf = region%radiInput%smoocf

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

END SUBROUTINE RADI_RFLO_FlimResSmoothingCoeff

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_rFLO_FlimResSmoothingCoeff.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







