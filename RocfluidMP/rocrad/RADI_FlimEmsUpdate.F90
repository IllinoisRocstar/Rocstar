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
! Purpose: update FLD radiation solution for steady flow or dual time-stepping
!          using explicit multistage scheme
!
! Description: none.
!
! Input: region = data of current region,
!
! Output: region%levels%radi%cv = new solution after each stage.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimEmsUpdate.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimEmsUpdate( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset, &
                            RFLO_ExtrapIntCellVec

#include "Indexing.h"
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic, idx
  
! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: iStage, ibc, iec, idxb, idxe
  REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:), rhs(:,:)
  REAL(RFREAL), POINTER :: vol(:), dt(:)
  REAL(RFREAL)          :: ark(5), fac, adtv, cfl

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_FlimEmsUpdate',&
  'RADI_FlimEmsUpdate.F90' )

! get parameters, dimensions and pointers -------------------------------------

  iStage = region%irkStep
  ark(:) = region%mixtInput%ark(:)
  cfl    = region%mixtInput%cfl

  IF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
    idxb = CV_RADI_ENER
    idxe = CV_RADI_ENER
  ENDIF

#ifdef RFLO
  iLev  = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  cv     => region%levels(iLev)%radi%cv
  cvOld  => region%levels(iLev)%radi%cvOld
  rhs    => region%levels(iLev)%radi%rhs
  vol    => region%levels(iLev)%grid%vol
  dt     => region%levels(iLev)%dt
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  cv     => region%radi%cv
  cvOld  => region%radi%cvOld
  rhs    => region%radi%rhs
  vol    => region%grid%vol
  dt     => region%dt
#endif

! update FLD radiation solution ------------------------------------------------

  fac = ark(iStage)*cfl
  DO ic = ibc,iec
    adtv  = fac*dt(ic)/vol(ic)
    DO idx = idxb,idxe
      rhs(idx,ic) = adtv*rhs(idx,ic)
    ENDDO
  ENDDO

  IF (global%solverType == SOLV_IMPLICIT) THEN
    fac = 1.5_RFREAL*ark(iStage)*cfl/global%dtMin
    DO ic = ibc,iec
      adtv = 1._RFREAL/(1._RFREAL+fac*dt(ic))
      DO idx = idxb,idxe
        cv(idx,ic) = cvOld(idx,ic) - adtv*rhs(idx,ic)
      ENDDO
    ENDDO
  ELSE
    DO ic=ibc,iec
      DO idx = idxb,idxe
        cv(idx,ic) = cvOld(idx,ic) - rhs(idx,ic)
      ENDDO
    ENDDO
  ENDIF

! extrapolate solution to dummy cells

#ifdef RFLO
  CALL RFLO_ExtrapIntCellVec( region,idxb,idxe,cv )
#endif
#ifdef RFLU
!  CALL RFLU_ExtrapIntCellVec( region,idxb,idxe,cv )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FlimEmsUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimEmsUpdate.F90,v $
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








