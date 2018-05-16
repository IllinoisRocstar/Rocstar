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
! Purpose: main driver of ROCFLO-MP.
!
! Description: gets the case name + verbosity level and calls the flow solver.
!
! Input: case name & verbosity level from the list of arguments
!
! Output: accurate solution.
!
! Notes: none
!
!******************************************************************************
!
! $Id: Main.F90,v 1.5 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
!******************************************************************************

PROGRAM Main

  USE ModDataTypes
  USE ModError
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_level,t_region
  USE ModParameters
  USE ModInterfaces, ONLY : RFLO_InitFlowSolver, RFLO_FlowSolver, &
                            RFLO_EndFlowSolver
  IMPLICIT NONE

! ... local variables
  CHARACTER(CHRLEN) :: casename, verbosity

  INTEGER :: verbLevel, dIterSystem, errorFlag

  REAL(RFREAL) :: dTimeSystem

  TYPE(t_region), POINTER :: regions(:)
  TYPE(t_global), POINTER :: global

!******************************************************************************
! read case name, verbosity level

  CALL GETARG(1,casename)
  CALL GETARG(2,verbosity)

  IF (LEN_TRIM(casename) == 0 .OR. LEN_TRIM(verbosity) == 0) THEN
    WRITE(STDOUT,'(A)')      SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Usage: rflomp <casename> <verbosity>'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'       verbosity = 0 - no output'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                 = 1 - moderate output'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                 = 2 - output all'
    WRITE(STDOUT,'(A)')      SOLVER_NAME
    STOP
  ENDIF   ! LEN_TRIM

  READ(verbosity,*) verbLevel

! initialize ------------------------------------------------------------------

  ALLOCATE( global,STAT=errorFlag )
  IF (errorFlag /= ERR_NONE) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME//' ERROR - pointer allocation failed.'
    STOP
  ENDIF

  CALL RFLO_InitFlowSolver( casename,verbLevel,global,regions )

! call the flow solver --------------------------------------------------------

  dTimeSystem = global%maxTime - global%currentTime
  dIterSystem = global%MaxIter - global%currentIter

  CALL RFLO_FlowSolver( dTimeSystem,dIterSystem,regions )

! shut down -------------------------------------------------------------------

  CALL RFLO_EndFlowSolver( regions )

END PROGRAM Main

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Main.F90,v $
! Revision 1.5  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/05/03 03:07:27  haselbac
! Removed RFLU stuff
!
! Revision 1.2  2005/04/15 15:07:46  haselbac
! Cosmetics only
!
! Revision 1.1  2004/12/01 21:29:33  haselbac
! Initial revision after changing case
!
! Revision 1.9  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2002/10/05 19:40:34  haselbac
! GENX integration of RFLU
!
! Revision 1.4  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.3  2002/09/09 16:47:11  haselbac
! Added RFLU calls
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/08/15 19:43:11  jblazek
! Prepared solver to be optionally linked with an external driver.
!
! Revision 1.11  2002/07/25 16:18:02  haselbac
! Fix problem completely, cpp apparently confused by double backslashes
!
! Revision 1.10  2002/07/25 16:04:18  haselbac
! Corrected output, leads to error on popovich
!
! Revision 1.9  2002/07/25 14:47:23  haselbac
! Added call to RFLU_InitGlobal
!
! Revision 1.8  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.7  2002/06/17 13:32:55  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.6  2002/04/04 19:41:09  jblazek
! Moved dTime and dIter test from main to flowSolver.
!
! Revision 1.5  2002/03/14 19:02:29  haselbac
! Modified argument list of RFLU calls
!
! Revision 1.4  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
! Revision 1.3  2002/01/02 16:00:03  jblazek
! Added input for multigrid parameters.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************






