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
! Purpose: register turbulence variables with GenX.
!
! Description: none.
!
! Input: regions = patches and region (volume) variables
!
! Output: to Roccom via RFLO_initGenxInterface.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_InitGenxInterface.F90,v 1.7 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_InitGenxInterface( regions,wins,winv )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  CHARACTER(CHRLEN) :: wins, winv
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
#ifdef STATS
  CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
  INTEGER :: iStat
#endif
  INTEGER :: iLev, ibc, iec, pid, errorFlag, ilb
  INTEGER :: iCOff, ijCOff, idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_InitGenxInterface',&
  'TURB_InitGenxInterface.F90' )

! input data (currently none) -------------------------------------------------

! output global data

  CALL COM_new_dataitem( TRIM(winv)//'.esg1Sum' ,'w',COM_DOUBLE,1,'J/(m^3s)' )
  CALL COM_new_dataitem( TRIM(winv)//'.esg4Sum' ,'w',COM_DOUBLE,1,'J/(m^3s)' )

! output surface data (currently none)
 
! output volume (for visualization) and restart data 

  CALL COM_new_dataitem( TRIM(winv)//'.mut' ,'e',COM_DOUBLE,1,'kg/(ms)' )
  CALL COM_new_dataitem( TRIM(winv)//'.lens','e',COM_DOUBLE,1,'m' )
  CALL COM_new_dataitem( TRIM(winv)//'.vort','e',COM_DOUBLE,3,'1/s' )

! statistics

#ifdef STATS
  IF ((global%flowType == FLOW_UNSTEADY) .AND. (global%doStat == ACTIVE)) THEN
    IF (global%turbNStat > 0) THEN
      statNm => global%turbStatNm
      DO iStat=1,global%turbNStat

        CALL COM_new_dataitem( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),'e', &
                                COM_DOUBLE,1,TRIM(statNm(1,2,iStat)) )
      ENDDO
    ENDIF
  ENDIF
#endif

! store pointers to variables, loop over all regions --------------------------

! global data

  CALL COM_set_array( TRIM(winv)//'.esg1Sum' ,0, global%esg1Sum )
  CALL COM_set_array( TRIM(winv)//'.esg4Sum' ,0, global%esg4Sum )

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      iLev = regions(iReg)%currLevel

! --- volume data

      pid = iReg*REGOFF
      
      CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                                jdcbeg,jdcend,kdcbeg,kdcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
      ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
      iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

      IF (regions(iReg)%turbInput%modelClass /= MODEL_RANS) THEN
        ALLOCATE( regions(iReg)%levels(iLev)%turb%lens(ibc:iec), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
      ENDIF
      IF (regions(iReg)%mixtInput%turbModel == TURB_MODEL_NONE) THEN
        ALLOCATE( regions(iReg)%levels(iLev)%turb%vort(3,ibc:iec), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
      ENDIF

      regions(iReg)%levels(iLev)%mixt%tv   = 0._RFREAL
      regions(iReg)%levels(iLev)%turb%lens = 0._RFREAL
      regions(iReg)%levels(iLev)%turb%vort = 0._RFREAL

      ilb = LBOUND(regions(iReg)%levels(iLev)%mixt%tv,2)

      CALL COM_set_array( TRIM(winv)//'.lens' ,pid, &
           regions(iReg)%levels(iLev)%turb%lens )

      CALL COM_set_array( TRIM(winv)//'.mut',pid, &
           regions(iReg)%levels(iLev)%mixt%tv(3,ilb),4)

      CALL COM_set_array( TRIM(winv)//'.vort',pid, &
           regions(iReg)%levels(iLev)%turb%vort(1,ilb))

! --- statistics

#ifdef STATS
      IF ((global%flowType==FLOW_UNSTEADY) .AND. (global%doStat==ACTIVE)) THEN
        IF (global%turbNStat > 0) THEN
          DO iStat=1,global%turbNStat
            CALL COM_set_array( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)), pid, &
              regions(iReg)%levels(iLev)%turb%tav(iStat,ilb), global%turbNStat)
          ENDDO
        ENDIF
      ENDIF
#endif

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_InitGenxInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_InitGenxInterface.F90,v $
! Revision 1.7  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2005/12/07 20:03:43  wasistho
! removed attemp to store actual tav i.o. accumulated
!
! Revision 1.4  2005/12/07 04:44:06  wasistho
! modified statistics treatment
!
! Revision 1.3  2005/12/07 02:23:50  wasistho
! added integrTime with eps
!
! Revision 1.2  2005/12/06 21:53:05  wasistho
! devided and multiply tav with integrTime
!
! Revision 1.1  2004/12/01 21:24:00  haselbac
! Initial revision after changing case
!
! Revision 1.13  2004/06/30 04:06:22  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.12  2004/06/29 23:53:08  wasistho
! migrated to Roccom-3
!
! Revision 1.11  2004/06/07 23:05:36  wasistho
! provide Genx statistics names, units, and anytime-activation
!
! Revision 1.10  2004/05/12 20:05:36  wasistho
! added USE TURB_ModParameters
!
! Revision 1.9  2004/05/04 20:39:41  wasistho
! added RaNS/DES length scale variables: lens
!
! Revision 1.8  2004/03/02 00:07:51  wasistho
! added global variables esg1Sum and esg4Sum for turb.(exact) restart
!
! Revision 1.7  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/09/16 01:00:23  wasistho
! composed vorticities to a 3-components vector
!
! Revision 1.3  2003/09/10 21:11:37  jblazek
! Setting variables to zero.
!
! Revision 1.2  2003/08/14 20:06:58  jblazek
! Corrected bug associated with radiation flux qr.
!
! Revision 1.1  2003/08/09 02:09:26  wasistho
! added TURB and RADI_initGenxInterface
!
!******************************************************************************







