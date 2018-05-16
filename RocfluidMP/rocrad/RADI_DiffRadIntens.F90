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
! Purpose: compute diffusion approximation (Rosseland) radiation intensity
!          at specified directions.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%radi%radInt(:,n) = n-directions radiation intensities
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_DiffRadIntens.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_DiffRadIntens( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, nAng, iCOff, ijCOff, ijkC

  REAL(RFREAL)          :: stBoltz, pi, rad, avgFac
  REAL(RFREAL)          :: coefc, tempc, fact, rati, adir(3), thet, phi, sgrad
  REAL(RFREAL), POINTER :: dv(:,:), coef(:,:), radInt(:,:), wvInt(:,:)
  REAL(RFREAL), POINTER :: angles(:,:), goFact(:)

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_DiffRadIntens',&
  'RADI_DiffRadIntens.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel
  nAng = region%radiInput%nAng

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  dv     => region%levels(iLev)%mixt%dv
  goFact => region%levels(iLev)%radi%goFact
  coef   => region%levels(iLev)%radi%radCoef
  radInt => region%levels(iLev)%radi%radInt
  wvInt  => region%levels(iLev)%radi%wvInt
  angles => region%radiInput%angles

  stBoltz = region%radiInput%stBoltz
  pi      = global%pi
  rad     = global%rad
  avgFac  = 1._RFREAL/6._RFREAL

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC = IndIJK(i  ,j,k,iCOff,ijCOff)

        coefc = coef(ijkC,RADI_COEFF_EXTINCT)
        tempc = dv(DV_MIXT_TEMP,ijkC)
        fact  = goFact(ijkC)*stBoltz*tempc**3/pi
        rati  = 4._RFREAL/coefc

! ----- wvInt contains grad(T)n+s+e+w+f+b
        wvInt(:,ijkC) = avgFac*wvInt(:,ijkC)

! ----- get intensity in all given directions
        DO l = 1,nAng
          thet    = angles(l,RADI_ANGLE_POLAR)
          phi     = angles(l,RADI_ANGLE_AZIMU)
          adir(1) = SIN(thet)*COS(phi)
          adir(2) = SIN(thet)*SIN(phi)
          adir(3) = COS(thet)
          sgrad   = adir(1)*wvInt(XCOORD,ijkC) + &
                    adir(2)*wvInt(YCOORD,ijkC) + &
                    adir(3)*wvInt(ZCOORD,ijkC)
          radInt(l,ijkC) = fact*( tempc - rati*sgrad)
        ENDDO ! l

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_DiffRadIntens

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_DiffRadIntens.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.3  2003/07/30 22:24:43  wasistho
! enter part and smoke data into radiation
!
! Revision 1.2  2003/07/23 03:14:22  wasistho
! cured baby illness
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!******************************************************************************







