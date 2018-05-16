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
! Purpose: utility functions for the dual time-steping.
!
! Description: file contains the following subroutines:
!
!  - DualTstInit    = initialize flow solutions at all time levels
!  - DualTstPredict = guess start solution for subiterations
!  - DualTstSterm   = compute the source term
!  - DualTstShift   = shift time levels.
!
! Input: region    = current region
!        regions   = all regions
!        timeLevel = time level of conserved variables to be stored (0-2).
!
! Output: region(s) = updated values for region(s).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_DualTstUtil.F90,v 1.6 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_DualTstInit( regions,timeLevel )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
#ifdef TURB
  USE TURB_ModParameters
#endif
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: timeLevel

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, ic, id

! ... local variables
  INTEGER :: iLev, ibc, iec
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend, iCOff, ijCOff

  REAL(RFREAL), POINTER :: cv(:,:), cvn(:,:)

  TYPE(t_global), POINTER :: global

#ifdef TURB
  LOGICAL :: turbUsed
#endif

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_DualTstInit',&
  'RFLO_DualTstUtil.F90' )

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE) THEN             ! on my processor

      iLev = regions(iReg)%currLevel
      CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                                jdcbeg,jdcend,kdcbeg,kdcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
      ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
      iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

      cv => regions(iReg)%levels(iLev)%mixt%cv
      IF      (timeLevel == 0) THEN
        cvn => regions(iReg)%levels(iLev)%mixt%cvn
      ELSE IF (timeLevel == 1) THEN
        cvn => regions(iReg)%levels(iLev)%mixt%cvn1
      ELSE
        cvn => regions(iReg)%levels(iLev)%mixt%cvn2
      ENDIF

      DO ic=ibc,iec
        cvn(CV_MIXT_DENS,ic) = cv(CV_MIXT_DENS,ic)
        cvn(CV_MIXT_XMOM,ic) = cv(CV_MIXT_XMOM,ic)
        cvn(CV_MIXT_YMOM,ic) = cv(CV_MIXT_YMOM,ic)
        cvn(CV_MIXT_ZMOM,ic) = cv(CV_MIXT_ZMOM,ic)
        cvn(CV_MIXT_ENER,ic) = cv(CV_MIXT_ENER,ic)
      ENDDO

! --- turbulence part --------------------------------------------------------

#ifdef TURB
      turbUsed = (regions(iReg)%mixtInput%flowModel == FLOW_NAVST .AND. &
                  regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)

      IF (turbUsed) THEN
        IF (regions(iReg)%turbInput%modelClass == MODEL_RANS .AND. &
            regions(iReg)%turbInput%nCv > 0) THEN

          cv => regions(iReg)%levels(iLev)%turb%cv
          IF      (timeLevel == 0) THEN
            cvn => regions(iReg)%levels(iLev)%turb%cvn
          ELSE IF (timeLevel == 1) THEN
            cvn => regions(iReg)%levels(iLev)%turb%cvn1
          ELSE
            cvn => regions(iReg)%levels(iLev)%turb%cvn2
          ENDIF

          DO ic=ibc,iec
            DO id=1,regions(iReg)%turbInput%nCv
              cvn(id,ic) = cv(id,ic)
            ENDDO
          ENDDO
        ENDIF ! rans
    ENDIF     ! turbused
#endif

    ENDIF     ! region on this processor and active
  ENDDO       ! iReg

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_DualTstInit

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_DualTstPredict( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset, &
                            RFLO_GetDimensPhys, MixtureProperties
  USE ModError
  USE ModParameters
#ifdef TURB
  USE TURB_ModParameters
#endif
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, ic, id

! ... local variables
  INTEGER :: iLev, ibc, iec
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend, iCOff, ijCOff
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend

  REAL(RFREAL), POINTER :: cv(:,:), cvn(:,:), cvn1(:,:), cvn2(:,:)

  TYPE(t_global), POINTER :: global

#ifdef TURB
  LOGICAL :: turbUsed
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_DualTstPredict',&
  'RFLO_DualTstUtil.F90' )

  iLev = region%currLevel
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  cv   => region%levels(iLev)%mixt%cv
  cvn  => region%levels(iLev)%mixt%cvn
  cvn1 => region%levels(iLev)%mixt%cvn1
  cvn2 => region%levels(iLev)%mixt%cvn2

! cv(pseudo) = cv(n) + (3cv(n) - 4cv(n-1) + cv(n-2))/2 -----------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        ic = IndIJK(i,j,k,iCOff,ijCOff)
        cv(CV_MIXT_DENS,ic) = cvn(CV_MIXT_DENS,ic) + &
                              0.5_RFREAL*(3._RFREAL*cvn (CV_MIXT_DENS,ic)-&
                                          4._RFREAL*cvn1(CV_MIXT_DENS,ic)+&
                                                    cvn2(CV_MIXT_DENS,ic))
        cv(CV_MIXT_XMOM,ic) = cvn(CV_MIXT_XMOM,ic) + &
                              0.5_RFREAL*(3._RFREAL*cvn (CV_MIXT_XMOM,ic)-&
                                          4._RFREAL*cvn1(CV_MIXT_XMOM,ic)+&
                                                    cvn2(CV_MIXT_XMOM,ic))
        cv(CV_MIXT_YMOM,ic) = cvn(CV_MIXT_YMOM,ic) + &
                              0.5_RFREAL*(3._RFREAL*cvn (CV_MIXT_YMOM,ic)-&
                                          4._RFREAL*cvn1(CV_MIXT_YMOM,ic)+&
                                                    cvn2(CV_MIXT_YMOM,ic))
        cv(CV_MIXT_ZMOM,ic) = cvn(CV_MIXT_ZMOM,ic) + &
                              0.5_RFREAL*(3._RFREAL*cvn (CV_MIXT_ZMOM,ic)-&
                                          4._RFREAL*cvn1(CV_MIXT_ZMOM,ic)+&
                                                    cvn2(CV_MIXT_ZMOM,ic))
        cv(CV_MIXT_ENER,ic) = cvn(CV_MIXT_ENER,ic) + &
                              0.5_RFREAL*(3._RFREAL*cvn (CV_MIXT_ENER,ic)-&
                                          4._RFREAL*cvn1(CV_MIXT_ENER,ic)+&
                                                    cvn2(CV_MIXT_ENER,ic))
      ENDDO
    ENDDO
  ENDDO

  IF (region%mixtInput%gasModel == GAS_MODEL_TCPERF) THEN
    CALL MixtureProperties( region,ibc,iec,.false. )
  ELSE
    CALL MixtureProperties( region,ibc,iec,.true.  )
  ENDIF

! turbulence part ------------------------------------------------------------

#ifdef TURB
  turbUsed = (region%mixtInput%flowModel == FLOW_NAVST .AND. &
              region%mixtInput%turbModel /= TURB_MODEL_NONE)

  IF (turbUsed) THEN
    IF (region%turbInput%modelClass == MODEL_RANS .AND. &
        region%turbInput%nCv > 0) THEN

      cv   => region%levels(iLev)%turb%cv
      cvn  => region%levels(iLev)%turb%cvn
      cvn1 => region%levels(iLev)%turb%cvn1
      cvn2 => region%levels(iLev)%turb%cvn2

      DO k=kpcbeg,kpcend
        DO j=jpcbeg,jpcend
          DO i=ipcbeg,ipcend
            ic = IndIJK(i,j,k,iCOff,ijCOff)

            DO id=1,region%turbInput%nCv
              cv(id,ic) = cvn(id,ic) + &
                          0.5_RFREAL*(3._RFREAL*cvn (id,ic)-&
                                      4._RFREAL*cvn1(id,ic)+&
                                                cvn2(id,ic))
            ENDDO ! id
          ENDDO
        ENDDO
      ENDDO
    ENDIF  ! rans
  ENDIF    ! turbused
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_DualTstPredict

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_DualTstSterm( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
#ifdef TURB
  USE TURB_ModParameters
#endif
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic, id

! ... local variables
  INTEGER :: iLev, ibc, iec
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend, iCOff, ijCOff

  REAL(RFREAL)          :: rdt2, rdt5
  REAL(RFREAL), POINTER :: cvn(:,:), cvn1(:,:), sDual(:,:)
  REAL(RFREAL), POINTER :: vol(:), volOld(:)

  TYPE(t_global), POINTER :: global

#ifdef TURB
  LOGICAL :: turbUsed
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_DualTstSterm',&
  'RFLO_DualTstUtil.F90' )

  iLev = region%currLevel
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  cvn    => region%levels(iLev)%mixt%cvn
  cvn1   => region%levels(iLev)%mixt%cvn1
  sDual  => region%levels(iLev)%mixt%sDual
  vol    => region%levels(iLev)%grid%vol
  rdt2   =  2._RFREAL/global%dtMin
  rdt5   =  0.5_RFREAL/global%dtMin

! compute dual tst source term Q* --------------------------------------------

  IF (region%mixtInput%moveGrid) THEN
    volOld => region%levels(iLev)%gridOld%vol
    DO ic=ibc,iec
      sDual(CV_MIXT_DENS,ic) = rdt2*vol(ic)   *cvn (CV_MIXT_DENS,ic) - &
                               rdt5*volOld(ic)*cvn1(CV_MIXT_DENS,ic)
      sDual(CV_MIXT_XMOM,ic) = rdt2*vol(ic)   *cvn (CV_MIXT_XMOM,ic) - &
                               rdt5*volOld(ic)*cvn1(CV_MIXT_XMOM,ic)
      sDual(CV_MIXT_YMOM,ic) = rdt2*vol(ic)   *cvn (CV_MIXT_YMOM,ic) - &
                               rdt5*volOld(ic)*cvn1(CV_MIXT_YMOM,ic)
      sDual(CV_MIXT_ZMOM,ic) = rdt2*vol(ic)   *cvn (CV_MIXT_ZMOM,ic) - &
                               rdt5*volOld(ic)*cvn1(CV_MIXT_ZMOM,ic)
      sDual(CV_MIXT_ENER,ic) = rdt2*vol(ic)   *cvn (CV_MIXT_ENER,ic) - &
                               rdt5*volOld(ic)*cvn1(CV_MIXT_ENER,ic)
    ENDDO
  ELSE
    DO ic=ibc,iec
      sDual(CV_MIXT_DENS,ic) = vol(ic)*(rdt2*cvn (CV_MIXT_DENS,ic)- &
                                        rdt5*cvn1(CV_MIXT_DENS,ic))
      sDual(CV_MIXT_XMOM,ic) = vol(ic)*(rdt2*cvn (CV_MIXT_XMOM,ic)- &
                                        rdt5*cvn1(CV_MIXT_XMOM,ic))
      sDual(CV_MIXT_YMOM,ic) = vol(ic)*(rdt2*cvn (CV_MIXT_YMOM,ic)- &
                                        rdt5*cvn1(CV_MIXT_YMOM,ic))
      sDual(CV_MIXT_ZMOM,ic) = vol(ic)*(rdt2*cvn (CV_MIXT_ZMOM,ic)- &
                                        rdt5*cvn1(CV_MIXT_ZMOM,ic))
      sDual(CV_MIXT_ENER,ic) = vol(ic)*(rdt2*cvn (CV_MIXT_ENER,ic)- &
                                        rdt5*cvn1(CV_MIXT_ENER,ic))
    ENDDO
  ENDIF

! turbulence part ------------------------------------------------------------

#ifdef TURB
  turbUsed = (region%mixtInput%flowModel == FLOW_NAVST .AND. &
              region%mixtInput%turbModel /= TURB_MODEL_NONE)

  IF (turbUsed) THEN
    IF (region%turbInput%modelClass == MODEL_RANS .AND. &
        region%turbInput%nCv > 0) THEN

      cvn    => region%levels(iLev)%turb%cvn
      cvn1   => region%levels(iLev)%turb%cvn1
      sDual  => region%levels(iLev)%turb%sDual

      IF (region%mixtInput%moveGrid) THEN
        volOld => region%levels(iLev)%gridOld%vol
        DO ic=ibc,iec
          DO id=1,region%turbInput%nCv
            sDual(id,ic) = rdt2*vol(ic)*cvn (id,ic) - &
                           rdt5*volOld(ic)*cvn1(id,ic)
          ENDDO
        ENDDO
      ELSE
        DO ic=ibc,iec
          DO id=1,region%turbInput%nCv
            sDual(id,ic) = vol(ic)*(rdt2*cvn (id,ic)- &
                                    rdt5*cvn1(id,ic))
          ENDDO
        ENDDO
      ENDIF  ! movegrid
    ENDIF    ! rans
  ENDIF      ! turbused
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_DualTstSterm

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_DualTstShift( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
#ifdef TURB
  USE TURB_ModParameters
#endif
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic, id

! ... local variables
  INTEGER :: iLev, ibc, iec
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend, iCOff, ijCOff

  REAL(RFREAL), POINTER :: cv(:,:), cvn(:,:), cvn1(:,:), cvn2(:,:)

  TYPE(t_global), POINTER :: global

#ifdef TURB
  LOGICAL :: turbUsed
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_DualTstShift',&
  'RFLO_DualTstUtil.F90' )

  iLev = region%currLevel
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  cv   => region%levels(iLev)%mixt%cv
  cvn  => region%levels(iLev)%mixt%cvn
  cvn1 => region%levels(iLev)%mixt%cvn1
  cvn2 => region%levels(iLev)%mixt%cvn2

  DO ic=ibc,iec
    cvn2(CV_MIXT_DENS,ic) = cvn1(CV_MIXT_DENS,ic)
    cvn2(CV_MIXT_XMOM,ic) = cvn1(CV_MIXT_XMOM,ic)
    cvn2(CV_MIXT_YMOM,ic) = cvn1(CV_MIXT_YMOM,ic)
    cvn2(CV_MIXT_ZMOM,ic) = cvn1(CV_MIXT_ZMOM,ic)
    cvn2(CV_MIXT_ENER,ic) = cvn1(CV_MIXT_ENER,ic)
    cvn1(CV_MIXT_DENS,ic) = cvn (CV_MIXT_DENS,ic)
    cvn1(CV_MIXT_XMOM,ic) = cvn (CV_MIXT_XMOM,ic)
    cvn1(CV_MIXT_YMOM,ic) = cvn (CV_MIXT_YMOM,ic)
    cvn1(CV_MIXT_ZMOM,ic) = cvn (CV_MIXT_ZMOM,ic)
    cvn1(CV_MIXT_ENER,ic) = cvn (CV_MIXT_ENER,ic)
    cvn (CV_MIXT_DENS,ic) = cv  (CV_MIXT_DENS,ic)
    cvn (CV_MIXT_XMOM,ic) = cv  (CV_MIXT_XMOM,ic)
    cvn (CV_MIXT_YMOM,ic) = cv  (CV_MIXT_YMOM,ic)
    cvn (CV_MIXT_ZMOM,ic) = cv  (CV_MIXT_ZMOM,ic)
    cvn (CV_MIXT_ENER,ic) = cv  (CV_MIXT_ENER,ic)
  ENDDO

! turbulence part ------------------------------------------------------------

#ifdef TURB
  turbUsed = (region%mixtInput%flowModel == FLOW_NAVST .AND. &
              region%mixtInput%turbModel /= TURB_MODEL_NONE)

  IF (turbUsed) THEN
    IF (region%turbInput%modelClass == MODEL_RANS .AND. &
        region%turbInput%nCv > 0) THEN

      cv   => region%levels(iLev)%turb%cv
      cvn  => region%levels(iLev)%turb%cvn
      cvn1 => region%levels(iLev)%turb%cvn1
      cvn2 => region%levels(iLev)%turb%cvn2

      DO ic=ibc,iec
        DO id=1,region%turbInput%nCv
          cvn2(id,ic) = cvn1(id,ic)
          cvn1(id,ic) = cvn (id,ic)
          cvn (id,ic) = cv  (id,ic)
        ENDDO
      ENDDO
    ENDIF    ! rans
  ENDIF      ! turbused
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_DualTstShift

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_DualTstUtil.F90,v $
! Revision 1.6  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.3  2004/12/09 22:16:34  wasistho
! added data turbulence
!
! Revision 1.2  2004/12/04 07:22:58  wasistho
! finish up dual tst
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/07/08 21:21:37  jblazek
! Modified start up procedure for dual-time stepping.
!
! Revision 1.1  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
!******************************************************************************










