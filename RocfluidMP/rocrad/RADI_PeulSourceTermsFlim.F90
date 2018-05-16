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
! Purpose: added radiation source terms to the rhs of PEUL Energy eq.
!
! Description: see formulation of FLD (Howell et al. JCP 184 (2003) 53-78)
!
! Input: region = data of current region,
!
! Output: region%levels%peul%rhs, radiaton source terms added to PEUL rhs
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_PeulSourceTermsFlim.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_PeulSourceTermsFlim( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, ijkC, iCon
  
! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: ibc, iec, nConstit
  REAL(RFREAL), POINTER :: dv(:,:), rhs(:,:), vol(:)
  REAL(RFREAL), POINTER :: rcv(:,:), coef(:,:)
  REAL(RFREAL) :: stBoltz, fourStb, planckAlp, bTermAlp, sTerm

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff,ijCOff
#endif

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_PeulSourceTermsFlim',&
  'RADI_PeulSourceTermsFlim.F90' )

! get coefficients and parameters ---------------------------------------------

  stBoltz = region%radiInput%stBoltz
  fourStb = 4._RFREAL*stBoltz

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev  = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

#ifdef PEUL
  dv    => region%levels(iLev)%peul%dv
  rhs   => region%levels(iLev)%peul%rhs
#endif
  rcv   => region%levels(iLev)%radi%cv
!  coef  => region%levels(iLev)%radi%peulCoef
  vol   => region%levels(iLev)%grid%vol

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
#endif
#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  dv    => region%peul%dv
  rhs   => region%peul%rhs
  rcv   => region%radi%cv
!  coef  => region%radi%peulCoef
  vol   => region%grid%vol

  DO ijkC = ibc,iec
#endif
! ----- src term = bbSourceAlp-planckAlp*c*Er ! bbSourceAlp=planckAlp*bTermAlp

!        planckAlp = coef(ijkC,RADI_COEFF_PLANCK)
!        bTermAlp  = fourStb*dv(DV_PEUL_TEMP,ijkC)**4
!        sTerm = bTermAlp - dv(DV_MIXT_SOUN,ijkC)*rcv(CV_RADI_ENER,ijkC)
!        sTerm = planckAlp*sTerm

!        rhs(CV_PEUL_ENER,ijkC) = rhs(CV_PEUL_ENER,ijkC) - vol(ijkC)*sTerm

#ifdef RFLO
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
#endif
#ifdef RFLU
  ENDDO       ! ijkC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_PeulSourceTermsFlim

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_PeulSourceTermsFlim.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/02/15 19:38:22  wasistho
! put peul within ifdef
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







