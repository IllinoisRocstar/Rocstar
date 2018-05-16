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
! Purpose: compute convective fluxes for RaNS class turbulence model, depending 
!          on model selected, and add them to dissipation in order to obtain 
!          the residual.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%turb%rhs = turb convective fluxes + turb dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansConvectiveFluxes.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansConvectiveFluxes( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLO
  USE ModInterfaces,      ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloRansSACentFlux, &
                                 TURB_FloRansSARoe1stFlux, &
                                 TURB_FloRansSARoe2ndFlux

#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic, idx

! ... local variables
  INTEGER :: ibc, iec, discrType, discrOrder, turbModel, idxbeg, idxend
  REAL(RFREAL), POINTER :: tdiss(:,:), trhs(:,:)

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_RansConvectiveFluxes',&
  'TURB_RansConvectiveFluxes.F90' )

  IF (region%turbInput%modelClass /= MODEL_RANS) GOTO 999

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                            kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  tdiss => region%levels(iLev)%turb%diss
  trhs  => region%levels(iLev)%turb%rhs
#endif
#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot
  
  tdiss => region%turb%diss
  trhs  => region%turb%rhs  
#endif

  turbModel  = region%mixtInput%turbModel
  discrType  = region%turbInput%spaceDiscr
  discrOrder = region%turbInput%spaceOrder

! select start and end index of 1st dimension depending on RaNS model selected

  IF ((turbModel == TURB_MODEL_SA) .OR. &
      (turbModel == TURB_MODEL_DESSA) .OR. &
      (turbModel == TURB_MODEL_HDESSA)) THEN
    idxbeg = CV_SA_NUTIL
    idxend = CV_SA_NUTIL
  ENDIF

! initialize residual (set = dissipation) -------------------------------------

  DO ic=ibc,iec
    DO idx=idxbeg,idxend
      trhs(idx,ic) = -tdiss(idx,ic)
    ENDDO
  ENDDO

! 2nd-order central scheme ----------------------------------------------------

#ifdef RFLO
  IF ( discrType ==RANS_DISCR_CEN) THEN
    IF (turbModel == TURB_MODEL_SA .OR. turbModel == TURB_MODEL_DESSA .OR. &
        turbModel == TURB_MODEL_HDESSA) THEN
      IF (discrOrder==RANS_DISCR_ORD1 .OR. discrOrder==RANS_DISCR_ORD2) THEN
        CALL TURB_FloRansSACentFlux( region )
      ENDIF
    ENDIF
  ENDIF
#endif

  IF (discrType == RANS_DISCR_UPW) THEN
    IF (turbModel == TURB_MODEL_SA .OR. turbModel == TURB_MODEL_DESSA .OR. &
        turbModel == TURB_MODEL_HDESSA) THEN
      IF (discrOrder==RANS_DISCR_ORD1) THEN
#ifdef RFLO
        CALL TURB_FloRansSARoe1stFlux( region )
#endif
#ifdef RFLU
!        CALL TURB_FluRansSARoe1stFlux( region )
#endif
      ELSEIF (discrOrder==RANS_DISCR_ORD2) THEN
#ifdef RFLO
        CALL TURB_FloRansSARoe2ndFlux( region )
#endif
#ifdef RFLU
!        CALL TURB_FluRansSARoe2ndFlux( region )
#endif
      ENDIF
    ENDIF
  ENDIF

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_RansConvectiveFluxes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansConvectiveFluxes.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/03/09 06:35:24  wasistho
! incorporated HDESSA
!
! Revision 1.3  2004/03/20 00:29:17  wasistho
! set turb_rflo_ransNumericalDiss to turb_ransNumerical..
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.2  2003/10/27 04:51:50  wasistho
! added RaNS upwind schemes
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







