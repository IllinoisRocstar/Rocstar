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
! Purpose: add mean pressure gradient term to rhs
!
! Description: mean pressure gradient is added as source term of axial
!              momentum to force the flow and periodicity
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%rhs = channel pressure gradient added to residual
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_CnlForceTerm.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_CnlForceTerm( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, ijkC0

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff
  REAL(RFREAL), POINTER :: dv(:,:)
#endif
#ifdef RFLU
  REAL(RFREAL), POINTER :: cv(:,:)
#endif

  REAL(RFREAL), POINTER :: rhs(:,:), vol(:)
  REAL(RFREAL) :: refMeanPgrad

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_CnlForceTerm.F90,v $'

  CALL RegisterFunction( region%global,'PERI_CnlForceTerm',&
  'PERI_CnlForceTerm.F90' )

! get dimensions, pointers and parameters -------------------------------------

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  dv   => region%levels(iLev)%mixt%dv
  rhs  => region%levels(iLev)%mixt%rhs
  vol  => region%levels(iLev)%grid%vol
#endif
#ifdef RFLU
  cv   => region%mixt%cv
  rhs  => region%mixt%rhs
  vol  => region%grid%vol
#endif

  refMeanPgrad = region%periInput%meanPgrad

#ifdef RFLO
  DO k = kpcbeg, kpcend
    DO j = jpcbeg, jpcend
      DO i = ipcbeg, ipcend
        ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)
#endif
#ifdef RFLU
  DO ijkC0 = 1, region%grid%nCells
#endif
        rhs(CV_MIXT_XMOM,ijkC0) = rhs(CV_MIXT_XMOM,ijkC0) + vol(ijkC0)*refMeanPgrad
#ifdef RFLO
        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) + &
                                  vol(ijkC0)*refMeanPgrad*dv(DV_MIXT_UVEL,ijkC0)
      ENDDO ! i
    ENDDO   ! j
  ENDDO     ! k
#endif
#ifdef RFLU
        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) + vol(ijkC0)* &
                                  refMeanPgrad*cv(CV_MIXT_XMOM,ijkC0)/ &
                                  cv(CV_MIXT_DENS,ijkC0)
  ENDDO     ! ijkC0
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE PERI_CnlForceTerm

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_CnlForceTerm.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/17 20:03:05  wasistho
! compiled with RFLU
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.3  2004/01/22 03:59:22  wasistho
! add dpdx*u1 to rhs(energy-eq.)
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!******************************************************************************







