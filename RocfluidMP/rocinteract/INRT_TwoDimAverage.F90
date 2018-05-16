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
! Purpose: Average the cv variables for fluid and smoke in some direction
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%cv
!         region%levels%peul%cv
!         region%levels%plag%cv
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_TwoDimAverage.F90,v 1.4 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_TwoDimAverage( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE INRT_ModParameters
#ifdef PLAG
  USE PLAG_ModParameters
#endif

  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i,j,k,iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend
  INTEGER :: iLev,iCOff,ijCOff,ijkC0
  INTEGER :: nPcls

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: mixtCv,peulCv,plagCv

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_TwoDimAverage.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_TwoDimAverage',&
  'INRT_TwoDimAverage.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  mixtCv => region%levels(iLev)%mixt%cv
#ifdef PEUL
  peulCv => region%levels(iLev)%peul%cv
#endif
#ifdef PLAG
  plagCv => region%levels(iLev)%plag%cv
#endif

  SELECT CASE (region%inrtInput%twoDAverage)

  CASE (0)
    CONTINUE ! do not average

  CASE (1)

! - zero z-momentum of gas ----------------------------------------------------

    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        DO i=ipcbeg,ipcend

          ijkC0 = IndIJK(i,j,k,iCOff,ijCOff)
          mixtCv(CV_MIXT_ZMOM,ijkC0) = 0._RFREAL

        END DO   ! i
      END DO     ! j
    END DO       ! k

! - zero z-momentum of Lagrangian particles -----------------------------------

#ifdef PLAG
    nPcls = 0
    IF (global%plagUsed) nPcls = region%levels(iLev)%plag%nPcls
    IF (nPcls > 0) THEN

      DO iPcls = 1,nPcls
        plagCv(CV_PLAG_ZMOM,iPcls) = 0._RFREAL
      END DO ! iPcls

    END IF ! nPcls
#endif

    CALL k_average(mixtCv,CV_MIXT_NEQS)

#ifdef PEUL
    IF (global%peulUsed) THEN
      CALL k_average(peulCv,region%levels(iLev)%peul%nCv)
    END IF ! peulUsed
#endif

  CASE DEFAULT
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)

  END SELECT ! twoDAverage

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

CONTAINS

! average in k-direction ------------------------------------------------------

  SUBROUTINE k_average(cv,nCv)

    REAL(RFREAL), POINTER    :: cv(:,:)
    INTEGER,      INTENT(IN) :: nCv

    INTEGER      :: iCv
    REAL(RFREAL) :: csum,vsum
    REAL(RFREAL), POINTER :: vol(:)

    vol => region%levels(iLev)%grid%vol

    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        DO iCv = 1,nCv

          csum = 0._RFREAL
          vsum = 0._RFREAL
          DO k=kpcbeg,kpcend
            ijkC0 = IndIJK(i,j,k,iCOff,ijCOff)
            csum = csum + vol(ijkC0) * cv(iCv,ijkC0)
            vsum = vsum + vol(ijkC0)
          END DO ! k
          csum = csum / vsum

          DO k=kpcbeg,kpcend
            ijkC0 = IndIJK(i,j,k,iCOff,ijCOff)
            cv(iCv,ijkC0) = csum
          END DO ! k

        END DO   ! iCv
      END DO     ! i
    END DO       ! j

  END SUBROUTINE k_average

END SUBROUTINE INRT_TwoDimAverage

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_TwoDimAverage.F90,v $
! Revision 1.4  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/02/15 20:17:03  wasistho
! put peul and plag within ifdef
!
! Revision 1.1  2004/12/01 21:56:47  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
!******************************************************************************







