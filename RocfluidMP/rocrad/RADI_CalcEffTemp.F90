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
! Purpose: compute the effective temperature as a compossite temperature of 
!          participating components (gas, particles, smoke, etc).
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%radi%dv = radiation effective temperature as dv comp.
!
! Notes: none. 
!
!******************************************************************************
!
! $Id: RADI_CalcEffTemp.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_CalcEffTemp( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff, ijkC

  REAL(RFREAL), POINTER :: dv(:,:), rdv(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_CalcEffTemp',&
  'RADI_CalcEffTemp.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  dv  => region%levels(iLev)%mixt%dv
  rdv => region%levels(iLev)%radi%dv

! get effective temperature --------------------------------------------------

  IF (region%radiInput%media == RADI_MEDIA_ARTIF) THEN
    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkC = IndIJK(i  ,j,k,iCOff,ijCOff)
          rdv(DV_RADI_TEFF,ijkC) = dv(DV_MIXT_TEMP,ijkC)
        ENDDO   ! i
      ENDDO   ! j
    ENDDO   ! k
  ELSE    ! real media
    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkC = IndIJK(i  ,j,k,iCOff,ijCOff)
          rdv(DV_RADI_TEFF,ijkC) = dv(DV_MIXT_TEMP,ijkC)
        ENDDO   ! i
      ENDDO   ! j
    ENDDO   ! k
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_CalcEffTemp

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_CalcEffTemp.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.1  2004/09/18 18:02:07  wasistho
! initial import for Limited Flux Diffusion radiation
!
!
!
!******************************************************************************







