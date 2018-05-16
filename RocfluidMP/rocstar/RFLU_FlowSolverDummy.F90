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
! Purpose: Dummy flow solver for GenX.
!
! Description: None.
!
! Input: 
!   globalGenx   	Pointer to global data
!   timeSystem		System time
!   dTimeSystem		System time step
!   genxHandleBc	Handle for BC update
!   genxHandleGm	Handle for geometry update.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_FlowSolverDummy.F90,v 1.4 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_FlowSolverDummy(globalGenx,timeSystem,dTimeSystem, & 
                                genxHandleBc,genxHandleGm)

  USE ModDataTypes
  USE ModRocstar, ONLY: t_globalGenx
  USE ModDataStruct, ONLY: t_level
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN) :: genxHandleBc,genxHandleGm
  DOUBLE PRECISION, INTENT(IN) :: dTimeSystem,timeSystem
  TYPE(t_globalGenx), POINTER :: globalGenx

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: pLevel

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_FlowSolverDummy.F90,v $ $Revision: 1.4 $'

! initialize some global variables

  global => globalGenx%global
  pLevel => globalGenx%levels(1)

  global%genxHandleBc = genxHandleBc
  global%genxHandleGm = genxHandleGm
  global%dTimeSystem  = dTimeSystem

! start time stepping

  CALL RegisterFunction(global,'RFLU_FlowSolverDummy',&
  'RFLU_FlowSolverDummy.F90')

  CALL ErrorStop(global,ERR_EXTERNAL_FUNCT,__LINE__)

! finalize

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_FlowSolverDummy

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_FlowSolverDummy.F90,v $
! Revision 1.4  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2002/10/17 19:55:55  haselbac
! Added timeSystem as second argument
!
! Revision 1.1  2002/10/05 18:28:18  haselbac
! Initial revision
!
!******************************************************************************







