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
! Purpose: dummy flow solver for GenX.
!
! Description: none.
!
! Input: globalGenx   = pointer to global data
!        dTimeSystem  = system time step
!        genxHandleBc = handle for BC update
!        genxHandleGm = handle for geometry update.
!
! Output: nothing useful.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_FlowSolverDummy.F90,v 1.3 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_FlowSolverDummy( globalGenx,timeSystem,dTimeSystem, &
                                 genxHandleBc,genxHandleGm )

  USE ModDataTypes
  USE ModRocstar, ONLY       : t_globalGenx
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  INTEGER, INTENT(in) :: genxHandleBc, genxHandleGm

  DOUBLE PRECISION, INTENT(in) :: dTimeSystem, timeSystem

  TYPE(t_globalGenx), POINTER :: globalGenx

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER :: iLev

  REAL(RFREAL), PARAMETER :: pressure=1.E+6  ! ### set ASAP ###

  TYPE(t_global), POINTER  :: global
  TYPE (t_region), POINTER :: regions(:)

!******************************************************************************
! initialize some global variables

  global  => globalGenx%global
  regions => globalGenx%regions

  global%genxHandleBc = genxHandleBc
  global%genxHandleGm = genxHandleGm
  global%dTimeSystem  = dTimeSystem
  global%currentTime  = timeSystem

  CALL RegisterFunction( global,'RFLO_FlowSolverDummy',&
  'RFLO_FlowSolverDummy.F90' )

! fill outgoing buffers

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE) THEN             ! on my processor

      iLev = regions(iReg)%currLevel
      regions(iReg)%levels(iLev)%mixt%dv(DV_MIXT_PRES,:) = pressure

      CALL RFLO_SendBoundaryValues( regions(iReg),.false. )
    ENDIF
  ENDDO

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_FlowSolverDummy

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_FlowSolverDummy.F90,v $
! Revision 1.3  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:23:48  haselbac
! Initial revision after changing case
!
! Revision 1.3  2003/05/15 02:57:00  jblazek
! Inlined index function.
!
! Revision 1.2  2002/10/25 19:09:32  jblazek
! Provided possibility to set pressure in dummy flow solver.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







