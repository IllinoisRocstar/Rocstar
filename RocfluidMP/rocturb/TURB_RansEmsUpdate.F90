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
! Purpose: update RaNS/DES solution for steady flow or dual time-stepping
!          using explicit multistage scheme
!
! Description: none.
!
! Input: region = data of current region,
!
! Output: region%levels%turb%cv = new solution after each stage.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansEmsUpdate.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansEmsUpdate( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloExtrapIntCellVec

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
  REAL(RFREAL), POINTER :: tcv(:,:), tcvOld(:,:), trhs(:,:), dsterm(:,:)
  REAL(RFREAL), POINTER :: vol(:), dt(:)
  REAL(RFREAL)          :: ark(5), fac, adtv, cfl

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_RansEmsUpdate',&
  'TURB_RansEmsUpdate.F90' )

! get parameters, dimensions and pointers -------------------------------------

  iStage = region%irkStep
  ark(:) = region%mixtInput%ark(:)
  cfl    = region%mixtInput%cfl

  IF (region%mixtInput%turbModel == TURB_MODEL_SA .OR. &
      region%mixtInput%turbModel == TURB_MODEL_DESSA .OR. &
      region%mixtInput%turbModel == TURB_MODEL_HDESSA) THEN
    idxb = CV_SA_NUTIL
    idxe = CV_SA_NUTIL
  ENDIF

#ifdef RFLO
  iLev  = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  tcv     => region%levels(iLev)%turb%cv
  tcvOld  => region%levels(iLev)%turb%cvOld
  trhs    => region%levels(iLev)%turb%rhs
  dsterm  => region%levels(iLev)%turb%dsterm
  vol     => region%levels(iLev)%grid%vol
  dt      => region%levels(iLev)%dt
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  tcv     => region%turb%cv
  tcvOld  => region%turb%cvOld
  trhs    => region%turb%rhs
  dsterm  => region%turb%dsterm
  vol     => region%grid%vol
  dt      => region%dt
#endif

! update RaNS solution -------------------------------------------------------

  fac = ark(iStage)*cfl
  DO ic = ibc,iec
    adtv  = fac*dt(ic)/vol(ic)
    DO idx = idxb,idxe
      trhs(idx,ic) = adtv*trhs(idx,ic)
    ENDDO
  ENDDO

  IF (global%solverType == SOLV_IMPLICIT) THEN
    fac = 1.5_RFREAL*ark(iStage)*cfl/global%dtMin
    DO ic = ibc,iec
      adtv = 1._RFREAL/(1._RFREAL+fac*dt(ic))
      DO idx = idxb,idxe
        tcv(idx,ic) = tcvOld(idx,ic) - adtv*trhs(idx,ic)
      ENDDO
    ENDDO
  ELSE
    DO ic=ibc,iec
      DO idx = idxb,idxe
        tcv(idx,ic) = tcvOld(idx,ic) - trhs(idx,ic)
      ENDDO
    ENDDO
  ENDIF

! clipping for stability

  IF (region%mixtInput%turbModel == TURB_MODEL_SA .OR. &
      region%mixtInput%turbModel == TURB_MODEL_DESSA .OR. &
      region%mixtInput%turbModel == TURB_MODEL_HDESSA) THEN
    DO ic=ibc,iec
      tcv(CV_SA_NUTIL,ic) = MAX( tcv(CV_SA_NUTIL,ic),REAL_SMALL )
    ENDDO !ic
  ENDIF

! extrapolate solution to dummy cells

#ifdef RFLO
  CALL TURB_FloExtrapIntCellVec( region,idxb,idxe,tcv )
#endif
#ifdef RFLU
!  CALL TURB_FluExtrapIntCellVec( region,idxb,idxe,tcv )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RansEmsUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansEmsUpdate.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/03/09 06:35:35  wasistho
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
! Revision 1.4  2004/02/12 03:46:24  wasistho
! filled in RaNS lengthscale in dummy cells
!
! Revision 1.3  2004/01/22 04:00:26  wasistho
! replaced dtMod by dt(ic)
!
! Revision 1.2  2003/10/21 20:32:03  wasistho
! added dt relaxation in steady flow due to RANS source term
!
! Revision 1.1  2003/10/16 20:20:03  wasistho
! installed RaNS in steady state flow (Exp.Mult.Stg)
!
!
!******************************************************************************








