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
! $Id: PEUL_InitGenxInterface.F90,v 1.3 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_InitGenxInterface( regions,wins,winv )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  CHARACTER(CHRLEN) :: wins, winv
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, ipt

! ... local variables
  INTEGER, PARAMETER :: ASCII_ZERO = 48   ! char representation of zero
  INTEGER :: iLev, pid, errorFlag, ilb, nPtypes, nPtypesMax

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_InitGenxInterface',&
  'PEUL_InitRocstarInterface.F90' )

! input data (currently none) -------------------------------------------------

! output surface data (currently none)

! output volume (for visualization) and restart data

  nPtypesMax = 0
  IF ( global%peulUsed ) THEN
    DO iReg=1,global%nRegions
      nPtypesMax = MAX(nPtypesMax,regions(iReg)%peulInput%nPtypes)
    ENDDO ! iReg
  ENDIF ! peulUsed

  DO ipt=1,nPtypesMax
    CALL COM_new_dataitem( TRIM(winv)//'.peul'//CHAR(ipt+ASCII_ZERO),'e', &
      COM_DOUBLE,1,'kg/(m^3)' )
  ENDDO ! ipt

! store pointers to variables, loop over all regions --------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
        global%peulUsed ) THEN                          ! and smoke used

      iLev = regions(iReg)%currLevel

! --- volume data

      pid = iReg*REGOFF

      ilb = LBOUND(regions(iReg)%levels(iLev)%peul%cv,2)

      nPtypes = regions(iReg)%peulInput%nPtypes

      DO ipt=1,nPtypes
        CALL COM_set_array( TRIM(winv)//'.peul'//CHAR(ipt+ASCII_ZERO), &
          pid,regions(iReg)%levels(iLev)%peul%cv(ipt,ilb),nPtypes)
      ENDDO ! ipt

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_InitGenxInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_InitGenxInterface.F90,v $
! Revision 1.3  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:23:43  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/06/30 04:05:23  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.4  2004/06/29 23:52:26  wasistho
! migrated to Roccom-3
!
! Revision 1.3  2004/03/05 22:08:58  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/01/21 16:34:45  fnajjar
! Fixed semantics of ipt+'0' to bypass Frost compiler error
!
! Revision 1.1  2003/11/21 22:21:23  fnajjar
! Initial import of Rocsmoke GenX interfaces
!
!******************************************************************************







