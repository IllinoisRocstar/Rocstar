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
! Purpose: compute convective fluxes for rocsmoke, add them to dissipation
!          in order to obtain the residual.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%peul%rhs = convective fluxes + dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ConvectiveFluxes.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ConvectiveFluxes( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
#ifdef RFLO
  USE ModInterfaces,      ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE PEUL_ModInterfaces, ONLY : PEUL_CentralFlux

#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic, iCv

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif
  INTEGER :: ibc, iec, spaceDiscr, spaceOrder

  REAL(RFREAL), POINTER :: diss(:,:), rhs(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_ConvectiveFluxes.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction( region%global,'PEUL_ConvectiveFluxes',&
  'PEUL_ConvectiveFluxes.F90' )

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                            kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  diss => region%levels(iLev)%peul%diss
  rhs  => region%levels(iLev)%peul%rhs
#endif

  spaceDiscr = region%mixtInput%spaceDiscr
  spaceOrder = region%mixtInput%spaceOrder

! initialize residual (set = dissipation) -------------------------------------

  DO ic=ibc,iec
    DO iCv=1,region%levels(iLev)%peul%nCv
      rhs(iCv,ic) = -diss(iCv,ic)
    ENDDO
  ENDDO

! 2nd-order central scheme ----------------------------------------------------

#ifdef RFLO
  IF (spaceDiscr==DISCR_CEN_SCAL .AND. &
      (spaceOrder==DISCR_ORDER_1 .OR. spaceOrder==DISCR_ORDER_2)) THEN
    CALL PEUL_CentralFlux( region )
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE PEUL_ConvectiveFluxes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ConvectiveFluxes.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:29  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!
!******************************************************************************







