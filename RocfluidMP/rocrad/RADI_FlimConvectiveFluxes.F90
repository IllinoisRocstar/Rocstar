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
! Purpose: compute convective fluxes for FLDTRAN radiation model, and 
!          add them to dissipation in order to obtain the residual.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%radi%rhs = radi convective fluxes + radi dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimConvectiveFluxes.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimConvectiveFluxes( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
#ifdef RFLO
  USE ModInterfaces,      ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE RADI_ModInterfaces, ONLY : RADI_FloFlimCentFlux

#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic, idx

! ... local variables
  INTEGER :: ibc, iec, discrType, discrOrder, idxbeg, idxend
  REAL(RFREAL), POINTER :: diss(:,:), rhs(:,:)

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FlimConvectiveFluxes',&
  'RADI_FlimConvectiveFluxes.F90' )

  IF (region%radiInput%radiModel /= RADI_MODEL_FLDTRAN) GOTO 999

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                            kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  diss => region%levels(iLev)%radi%diss
  rhs  => region%levels(iLev)%radi%rhs
#endif
#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot
  
  diss => region%radi%diss
  rhs  => region%radi%rhs  
#endif

  discrType  = region%radiInput%spaceDiscr
  discrOrder = region%radiInput%spaceOrder

! select start and end index of 1st dimension

  idxbeg = CV_RADI_ENER
  idxend = CV_RADI_ENER

! initialize residual (set = dissipation) -------------------------------------

  DO ic=ibc,iec
    DO idx=idxbeg,idxend
      rhs(idx,ic) = -diss(idx,ic)
    ENDDO
  ENDDO

! 2nd-order central scheme ----------------------------------------------------

#ifdef RFLO
  IF (discrType ==FLD_DISCR_CEN) THEN
    IF (discrOrder==FLD_DISCR_ORD1 .OR. discrOrder==FLD_DISCR_ORD2) THEN
      CALL RADI_FloFlimCentFlux( region )
    ENDIF
  ENDIF
#endif

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_FlimConvectiveFluxes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimConvectiveFluxes.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







