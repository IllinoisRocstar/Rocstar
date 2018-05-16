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
! ******************************************************************************
!
! Purpose: Main driver of ROCFLU-MP.
!
! Description: None.
!
! Input: 
!   casename    Case name
!   verbLevel   Verbosity
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rflump.F90,v 1.7 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rflump(caseString,verbLevel)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModParameters

  USE ModInterfaces, ONLY: RFLU_EndFlowSolver, &
                           RFLU_FlowSolver, & 
                           RFLU_InitFlowSolver

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString
  INTEGER, INTENT(IN) :: verbLevel
  
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: casename
  INTEGER :: dIterSystem,errorFlag
  REAL(RFREAL) :: dTimeSystem
  TYPE(t_level), POINTER :: levels(:)  
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start, allocate global pointer
! ******************************************************************************

  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME//' ERROR - pointer allocation failed.'
    STOP
  END IF ! global

  casename = caseString(1:LEN(caseString))

! ******************************************************************************
! Call initialization, solver, and  finalization
! ******************************************************************************

  CALL RFLU_InitFlowSolver(casename,verbLevel,global,levels)

  dTimeSystem = global%maxTime - global%currentTime
  dIterSystem = global%MaxIter - global%currentIter

  CALL RFLU_FlowSolver(dTimeSystem,dIterSystem,levels)

  CALL RFLU_EndFlowSolver(levels)

! ******************************************************************************
! Deallocate global pointer
! ******************************************************************************

  DEALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME//' ERROR - pointer deallocation failed.'
    STOP
  END IF ! global

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE rflump

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rflump.F90,v $
! Revision 1.7  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.4  2005/09/14 15:59:33  haselbac
! Minor clean-up
!
! Revision 1.3  2005/09/13 20:46:01  mtcampbe
! Moved profiling calls to (Init)(End)Flowsolver
!
! Revision 1.2  2005/07/07 22:45:14  haselbac
! Added profiling calls
!
! Revision 1.1  2005/05/03 02:55:45  haselbac
! Initial revision
!
! ******************************************************************************






