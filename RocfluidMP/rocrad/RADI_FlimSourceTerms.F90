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
! Purpose: added source terms of FLD transport equation to the rhs.
!
! Description: see formulation of FLD (Howell et al. JCP 184 (2003) 53-78)
!
! Input: region = data of current region,
!
! Output: region%levels%radi%rhs , source terms added to rhs
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimSourceTerms.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimSourceTerms( region )

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
  REAL(RFREAL), POINTER :: dv(:,:), vol(:)
  REAL(RFREAL), POINTER :: rcv(:,:), rrhs(:,:)
  REAL(RFREAL) :: stBoltz, sTerm, mPlanck

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff,ijCOff
#endif

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_FlimSourceTerms',&
  'RADI_FlimSourceTerms.F90' )

! get coefficients and parameters ---------------------------------------------

  stBoltz = region%radiInput%stBoltz

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev  = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  dv     => region%levels(iLev)%mixt%dv
  rcv    => region%levels(iLev)%radi%cv
  rrhs   => region%levels(iLev)%radi%rhs
  vol    => region%levels(iLev)%grid%vol

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
#endif
#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  dv     => region%mixt%dv
  rcv    => region%radi%cv
  rrhs   => region%radi%rhs
  vol    => region%grid%vol

  DO ijkC = ibc,iec
#endif
        sTerm = 0._RFREAL
!        DO iCon = 1, nConstit 
!          sTerm  = sTerm + volFrac(iCon,ijkC)*bbSource(iCon,ijkC)
!        ENDDO
        sTerm = sTerm - mPlanck*dv(DV_MIXT_SOUN,ijkC)*rcv(CV_RADI_ENER,ijkC)

        rrhs(CV_RADI_ENER,ijkC) = rrhs(CV_RADI_ENER,ijkC) + vol(ijkC)*sTerm

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

END SUBROUTINE RADI_FlimSourceTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimSourceTerms.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
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







