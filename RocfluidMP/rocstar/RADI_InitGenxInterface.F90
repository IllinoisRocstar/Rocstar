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
! Purpose: register radiation variables with GenX.
!
! Description: none.
!
! Input: regions = region (volume) variables
!
! Output: to Roccom via RFLO_initGenxInterface.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_InitGenxInterface.F90,v 1.3 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_InitGenxInterface( regions,wins,winv )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  CHARACTER(CHRLEN) :: wins, winv
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER :: iLev, ibc, iec, pid, icount, errorFlag, ilb
  INTEGER :: iCOff, ijCOff, idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RADI_InitGenxInterface',&
  'RADI_InitGenxInterface.F90' )

! input data (currently none) -------------------------------------------------

! output volume (for visualization) and/or restart data (currently none) ------

  CALL COM_new_dataitem( TRIM(winv)//'.xcof','e',COM_DOUBLE,1,'1/m')

! store pointers to variables, loop over all regions --------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      iLev   = regions(iReg)%currLevel
      icount = 0

! --- volume data

      pid = iReg*REGOFF
      
      CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                                jdcbeg,jdcend,kdcbeg,kdcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
      ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
      iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

      IF (.NOT.regions(iReg)%mixtInput%radiUsed) THEN
        ALLOCATE( regions(iReg)%levels(iLev)%radi%radCoef(ibc:iec,2), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
      ENDIF

      regions(iReg)%levels(iLev)%radi%radCoef = 0._RFREAL

      CALL COM_set_array( TRIM(winv)//'.xcof',pid, &
                          regions(iReg)%levels(iLev)%radi%radCoef )

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_InitGenxInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_InitGenxInterface.F90,v $
! Revision 1.3  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:23:46  haselbac
! Initial revision after changing case
!
! Revision 1.8  2004/06/30 04:05:43  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.7  2004/06/29 23:53:15  wasistho
! migrated to Roccom-3
!
! Revision 1.6  2003/11/20 16:40:33  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
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







