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
! Purpose: Check if there are negative values of concentration, and, when
!          indicated, report and possibly fix these values.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%peul%cv
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_EnforcePositivity.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_EnforcePositivity( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i,j,k,ipt

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend,nPtypes
  INTEGER :: iLev,iCOff,ijCOff,ijkC0
  INTEGER :: ncnt1,ncnt2,iloc1,jloc1,kloc1,iloc2,jloc2,kloc2

  LOGICAL :: clip,firstType

  REAL(RFREAL) :: minc,maxc,maxConc
  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_EnforcePositivity.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_EnforcePositivity',&
  'PEUL_EnforcePositivity.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv => region%levels(iLev)%peul%cv

  nPtypes = region%peulInput%nPtypes

  firstType = .true.

  DO ipt=1,nPtypes

    maxConc = region%peulInput%ptypes(ipt)%maxConc

    minc =  1.E20_RFREAL
    maxc = -1.E20_RFREAL
    ncnt1 = 0
    ncnt2 = 0
    iloc1 = CRAZY_VALUE_INT
    jloc1 = CRAZY_VALUE_INT
    kloc1 = CRAZY_VALUE_INT
    iloc2 = CRAZY_VALUE_INT
    jloc2 = CRAZY_VALUE_INT
    kloc2 = CRAZY_VALUE_INT

    clip = .false.

    SELECT CASE (region%peulInput%ptypes(ipt)%clipModel)

    CASE (PEUL_CLIP_MODEL_USED)
      clip = .true.

    END SELECT ! clipModel

    IF (region%peulInput%ptypes(ipt)%negReport == PEUL_NEG_REPORT_NONE &
        .AND. .NOT.clip) CYCLE

    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        DO i=ipcbeg,ipcend
          ijkC0 = IndIJK(i,j,k,iCOff,ijCOff)

          IF (cv(ipt,ijkC0) < minc) THEN
            minc = cv(ipt,ijkC0)
            iloc1 = i
            jloc1 = j
            kloc1 = k
          ENDIF ! minc

          IF (cv(ipt,ijkC0) > maxc) THEN
            maxc = cv(ipt,ijkC0)
            iloc2 = i
            jloc2 = j
            kloc2 = k
          ENDIF ! maxc

          IF (cv(ipt,ijkC0) < 0._RFREAL) THEN
            ncnt1 = ncnt1 + 1
            IF (clip) cv(ipt,ijkC0) = 0._RFREAL
          ENDIF ! cv(ipt,ijkC0)

          IF (cv(ipt,ijkC0) > maxConc) THEN
            ncnt2 = ncnt2 + 1
            IF (clip) cv(ipt,ijkC0) = maxConc
          ENDIF ! cv(ipt,ijkC0)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

    IF (ncnt1 > 0 .OR. ncnt2 > 0) THEN

      SELECT CASE (region%peulInput%ptypes(ipt)%negReport)

      CASE (PEUL_NEG_REPORT_USED)

        IF (firstType) THEN

          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                'Negative smoke concentration detected!'

          IF ( global%flowType == FLOW_UNSTEADY ) THEN
            WRITE(STDOUT,'(A,3X,A,1X,ES12.5,A,I1)') SOLVER_NAME, &
              'Current time and stage:',global%currentTime,', ',region%irkStep
          ELSE
            WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, &
              'Current iteration number:', global%currentIter
          ENDIF ! flowType

          WRITE(STDOUT,'(A,2X,12(1X,A))') SOLVER_NAME,    &
                                               ' Region', &
                                                 '  ipt', &
                                               ' # negs', &
                                       '    Least value', &
                                                 '    i', &
                                                 '    j', &
                                                 '    k', &
                                               ' # bigs', &
                                       ' Greatest value', &
                                                 '    i', &
                                                 '    j', &
                                                 '    k'
          firstType = .false.

        ENDIF ! firstType

        WRITE(STDOUT,'(A,I10,I6,I8,ES16.6,3I6,I8,ES16.6,3I6)') &
          SOLVER_NAME,region%iRegionGlobal,ipt,ncnt1,minc,iloc1,jloc1,kloc1, &
          ncnt2,maxc,iloc2,jloc2,kloc2

      END SELECT ! negReport

    ENDIF ! ncnt1 > 0 .OR. ncnt2 > 0

  ENDDO ! ipt

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_EnforcePositivity

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_EnforcePositivity.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:33  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/05/03 15:09:05  jferry
! added ability to clip smoke values greater than maximal space packing
!
! Revision 1.2  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/03/02 21:44:52  jferry
! Added clipping options
!
!******************************************************************************







