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
! Purpose: Print and write convergence history.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PrintWriteConvergence.F90,v 1.4 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PrintWriteConvergence(global)

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI

  USE ModTools, ONLY: MakeNonZero

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  REAL(RFREAL) :: resNorm
  REAL(RFREAL), DIMENSION(6) :: localVals,globalVals

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_PrintWriteConvergence',&
  'RFLU_PrintWriteConvergence.F90')

! ******************************************************************************
! Reduce data
! ******************************************************************************

  localVals(1) = global%forceX
  localVals(2) = global%forceY
  localVals(3) = global%forceZ
  localVals(4) = global%massIn
  localVals(5) = global%massOut
  localVals(6) = global%stopRun

  CALL MPI_Allreduce(localVals,globalVals,6,MPI_RFREAL,MPI_SUM,global%mpiComm, &
                     errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag

  global%forceX  = globalVals(1)
  global%forceY  = globalVals(2)
  global%forceZ  = globalVals(3)
  global%massIn  = globalVals(4)
  global%massOut = globalVals(5)
  global%stopRun = globalVals(6)

! ******************************************************************************
! Print and write global data
! ******************************************************************************

! ==============================================================================
! Steady flow 
! ==============================================================================

  IF ( (global%flowType == FLOW_STEADY) .AND. & 
       (global%myProcid == MASTERPROC) ) THEN
    IF ( global%currentIter == 1 ) THEN
      resNorm = 1.0_RFREAL
    ELSE
      resNorm = global%residual/MakeNonZero(global%resInit)
      resNorm = MakeNonZero(resNorm)
    END IF ! global%currentIter

    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,I6,1PE13.4,5E13.4)') SOLVER_NAME, & 
        global%currentIter,LOG10(resNorm), &
        global%forceX,global%forceY,global%forceZ, &
        global%massIn,global%massOut
    END IF ! global%verbLevel

    WRITE(IF_CONVER,'(I6,1PE13.4,5E13.4)') & 
      global%currentIter,LOG10(resNorm), &
      global%forceX,global%forceY,global%forceZ, &
      global%massIn,global%massOut

! ==============================================================================
! Unsteady flow 
! ==============================================================================

  ELSE IF ( (global%flowType == FLOW_UNSTEADY) .AND. & 
            (global%myProcid==MASTERPROC) ) THEN
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,1PE12.5,6E13.4)') SOLVER_NAME, & 
        global%currentTime,global%dtMin, &
        global%forceX,global%forceY,global%forceZ, &
        global%massIn,global%massOut
    END IF ! global%verbLevel

    WRITE(IF_CONVER,'(1PE12.5,6E13.4)') & 
      global%currentTime,global%dtMin, &
      global%forceX,global%forceY,global%forceZ, &
      global%massIn,global%massOut
  END IF ! global%flowType

! ******************************************************************************
! Start
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PrintWriteConvergence

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintWriteConvergence.F90,v $
! Revision 1.4  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.1  2005/04/15 15:07:10  haselbac
! Initial revision
!
! ******************************************************************************







