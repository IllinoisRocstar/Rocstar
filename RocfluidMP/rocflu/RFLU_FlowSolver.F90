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
! Purpose: Flow solver of Rocflu, essentially wrapper around time-stepping
!   routine.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_FlowSolver.F90,v 1.24 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

#ifdef GENX
SUBROUTINE RFLU_FlowSolver(globalGenx,timeSystem,dTimeSystem,genxHandleBc, & 
                           genxHandleGm)
#else
SUBROUTINE RFLU_FlowSolver(dTimeSystem,dIterSystem,levels)
#endif

  USE ModDataTypes
#ifdef GENX
  USE ModRocstar, ONLY: t_globalGenx
#endif  
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  
#ifdef GENX
  USE RFLU_ModGeometry
  USE RFLU_ModRocstarTools, ONLY: RFLU_GENX_InitBFLAG
#endif   

#ifdef PETSC
  USE RFLU_ModNewtonKrylov, ONLY: RFLU_NK_TimeStepping
#endif
  
  USE ModInterfaces, ONLY: RFLU_TimeStepping
#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_GetBoundaryValues,RFLU_PutBoundaryValuesAlpha
#endif

  IMPLICIT NONE

! ******************************************************************************
! Arguments
! ******************************************************************************

#ifdef GENX
  INTEGER :: dIterSystem
  INTEGER, INTENT(IN) :: genxHandleBc,genxHandleGm
  DOUBLE PRECISION, INTENT(IN) :: dTimeSystem,timeSystem
  TYPE(t_globalGenx), POINTER :: globalGenx
#else
  INTEGER, INTENT(IN) :: dIterSystem
  REAL(RFREAL), INTENT(IN) :: dTimeSystem
  TYPE(t_level), POINTER :: levels(:)  
#endif

! ******************************************************************************
! Locals
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString
#ifdef GENX  
  INTEGER :: iReg
  LOGICAL :: corrFlag  
#endif
  TYPE(t_region), POINTER :: pRegion,regions(:)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_FlowSolver.F90,v $ $Revision: 1.24 $'

#ifdef GENX
  global  => globalGenx%global
  regions => globalGenx%levels(1)%regions

  IF ( (global%currentTime-timeSystem) > 1.0E-10_RFREAL ) THEN
    corrFlag = .TRUE.
  ELSE 
    corrFlag = .FALSE.    
  END IF ! global 

  global%currentTime  = timeSystem
  global%timeStamp    = timeSystem
  
  global%genxHandleBc = genxHandleBc
  global%genxHandleGm = genxHandleGm
  dIterSystem         = 0
#else
  global  => levels(1)%regions(1)%global
  regions => levels(1)%regions  
#endif

  global%dTimeSystem = dTimeSystem

  CALL RegisterFunction(global,'RFLU_FlowSolver',&
  'RFLU_FlowSolver.F90')

#ifdef ROCPROF
  CALL FPROFILER_BEGINS("RFLU::FlowSolver")
#endif

! ******************************************************************************
! Start time stepping
! ******************************************************************************
                          
! ==============================================================================
! Write header for convergence history
! ==============================================================================

  IF ( global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME
    IF ( global%flowType == FLOW_STEADY ) THEN 
      WRITE(STDOUT,1000) SOLVER_NAME,SOLVER_NAME
    ELSE IF ( global%flowType == FLOW_UNSTEADY ) THEN 
      WRITE(STDOUT,1010) SOLVER_NAME,SOLVER_NAME
    END IF ! global%flowType
  END IF ! global%myProcid

#ifdef GENX
! ==============================================================================
! Get current boundary condition values for computing time steps
! ==============================================================================

  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)
    CALL COM_call_function(global%genxHandleBc,2,0.0_RFREAL,1)
    CALL RFLU_PutBoundaryValuesAlpha(pRegion)
    CALL COM_call_function(global%genxHandleBc,2,0.0_RFREAL,2)
    CALL RFLU_GENX_InitBFLAG(pRegion)
    CALL RFLU_GetBoundaryValues(pRegion)
  END DO ! iReg 

! ==============================================================================
! Update geometry if starting corrector step
! ==============================================================================

  IF ( corrFlag .EQV. .TRUE. ) THEN  
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Starting corrector step.'
    END IF ! global

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      CALL RFLU_BuildGeometry(pRegion)    
    END DO ! iReg 
  END IF ! corrFlag
#endif

! ******************************************************************************
! Call time-stepping routines
! ******************************************************************************

  IF (global%solverType == SOLV_EXPLICIT) THEN
    IF (global%cycleType == MGCYCLE_NO) THEN
      CALL RFLU_TimeStepping(dTimeSystem,dIterSystem,regions)
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%cycleType
  ELSE 
#ifdef PETSC
    CALL RFLU_NK_TimeStepping(dTimeSystem,dIterSystem,regions)
#endif
  ENDIF ! global%solverType

! ******************************************************************************
! End
! ****************************************************************************** 
  
#ifdef ROCPROF
  CALL FPROFILER_ENDS("RFLU::FlowSolver")
#endif
  
  CALL DeregisterFunction(global)

1000 FORMAT(A,2X,' iter',4X,'res-norm',5X,'force-x',6X,'force-y',6X,'force-z', &
            6X,'mass-in',6X,'mass-out',/,A,1X,84('-'))
1010 FORMAT(A,2X,' time',10X,'delta-t',6X,'force-x',6X,'force-y',6X,'force-z', &
            6X,'mass-in',6X,'mass-out'/,A,1X,90('-'))

END SUBROUTINE RFLU_FlowSolver

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_FlowSolver.F90,v $
! Revision 1.24  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.23  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.22  2007/04/20 16:07:49  mtcampbe
! Updating for burnout support function RFLU_GENX_InitBFLAG
!
! Revision 1.21  2005/09/14 15:59:33  haselbac
! Minor clean-up
!
! Revision 1.20  2005/09/13 20:49:39  mtcampbe
! Added profiling calls
!
! Revision 1.19  2005/08/03 18:30:32  hdewey2
! Add IF for solverType
!
! Revision 1.18  2005/08/02 18:26:14  hdewey2
! Added NK capability
!
! Revision 1.17  2004/10/19 19:29:17  haselbac
! Cosmetics only
!
! Revision 1.16  2003/06/20 22:34:58  haselbac
! Cosmetic changes
!
! Revision 1.15  2003/05/09 17:01:03  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.14  2003/03/05 20:39:31  jiao
! ACH: Added calls to get correct data at every time step inside/outside of PC
!
! Revision 1.13  2003/02/24 18:05:38  haselbac
! Bug fix and clean-up
!
! Revision 1.12  2003/02/24 17:25:20  haselbac
! Add geometry computation for PC iterations within GENX
!
! Revision 1.11  2003/02/24 14:50:33  haselbac
! Bug fix: Added missing initialization of timeStamp
!
! Revision 1.10  2002/10/19 22:22:23  haselbac
! Removed RFLU_GetBValues - not needed here with proper calls
!
! Revision 1.9  2002/10/19 16:13:19  haselbac
! Removed include for Roccom, cosmetic changes to output
!
! Revision 1.8  2002/10/17 20:04:52  haselbac
! Added timeSystem to argument list (GENX)
!
! Revision 1.7  2002/10/17 14:12:59  haselbac
! Added RFLU_GetBValues for proper restart (discussion with Jim J.)
!
! Revision 1.6  2002/10/05 19:21:30  haselbac
! GENX integration, some cosmetics
!
! Revision 1.5  2002/09/09 15:49:58  haselbac
! global now under regions
!
! Revision 1.4  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.3  2002/05/04 17:09:00  haselbac
! Uncommented writing of convergence file
!
! Revision 1.2  2002/04/11 19:02:21  haselbac
! Cosmetic changes and some preparation work
!
! Revision 1.1  2002/03/14 19:12:00  haselbac
! Initial revision
!
! ******************************************************************************







