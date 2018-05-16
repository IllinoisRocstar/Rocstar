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
! Purpose: solve the discretized governing equations.
!
! Description: none.
!
! Input: dTimeSystem = time step to run the solver
!        dIterSystem = no. of iterations to run the solver
!        regions     = dimension, BC`s and flow variables of all regions.
!
! Output: accurate solution.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_FlowSolver.F90,v 1.5 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef GENX
SUBROUTINE RFLO_FlowSolver( globalGenx,timeSystem,dTimeSystem,genxHandleBc, &
                            genxHandleGm )
#else
SUBROUTINE RFLO_FlowSolver( dTimeSystem,dIterSystem,regions )
#endif

  USE ModDataTypes
#ifdef GENX
  USE ModRocstar, ONLY       : t_globalGenx
#endif
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_TimeStepping, RFLO_Multigrid, &
        RFLO_GetBoundaryValues, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_SendBoundaryValues, &
        RFLO_DualTimeStepping, RFLO_DualMultigrid
  USE RFLO_ModGridMetrics, ONLY: RFLO_CalcGridMetrics
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE
#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ... parameters
#ifdef GENX
  INTEGER, INTENT(in) :: genxHandleBc, genxHandleGm

  DOUBLE PRECISION, INTENT(in) :: timeSystem, dTimeSystem

  TYPE(t_globalGenx), POINTER :: globalGenx
#else
  REAL(RFREAL) :: dTimeSystem
#endif

  INTEGER :: dIterSystem

  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: msg

  REAL(RFREAL) :: timerBeg, timerEnd, timerLoc, timerGlob

  TYPE(t_global), POINTER :: global

!******************************************************************************
! initialize some global variables

#ifdef GENX
  global  => globalGenx%global
  regions => globalGenx%regions

  global%genxHandleBc = genxHandleBc
  global%genxHandleGm = genxHandleGm
  dIterSystem         = 0
  IF ((global%currentTime-timeSystem) > 1.E-9_RFREAL) THEN
    global%predCorrIter = .true.
  ELSE
    global%predCorrIter = .false.
  ENDIF
  global%currentTime = timeSystem
  global%timeStamp   = timeSystem
#else
  global => regions(1)%global
#endif

  global%dTimeSystem = dTimeSystem

  CALL RegisterFunction( global,'RFLO_FlowSolver',&
  'RFLO_FlowSolver.F90' )

#ifdef GENX
! restore geometry if predictor-corrector iteration

  IF (global%predCorrIter) THEN
    IF (global%myProcid == MASTERPROC) WRITE(STDOUT,'(A)')  &
      SOLVER_NAME//' Restoring geometry (PC iteration) ...'
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor
        CALL RFLO_GenerateCoarseGrids( regions(iReg) )   ! coarsen finest grid
        CALL RFLO_CopyGeometryDummy( regions(iReg) )     ! copy to dummy nodes
        CALL RFLO_ExtrapolateGeometry( regions(iReg) )   ! extrapolate
      ENDIF     ! region on this processor and active
    ENDDO       ! iReg
    CALL RFLO_ExchangeGeometry( regions )                ! exchange geometry
    CALL RFLO_CalcGridMetrics( regions )
  ENDIF         ! predCorrIter==true

! get BC data from GenX at time=0

  CALL COM_call_function( genxHandleBc,2,0._RFREAL,1 )
  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE) THEN             ! on my processor
      CALL RFLO_SendBoundaryValues( regions(iReg),.false. )
    ENDIF     ! region on this processor and active
  ENDDO       ! iReg
  CALL COM_call_function( genxHandleBc,2,0._RFREAL,2 )

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE          .AND. &  ! on my processor
        regions(iReg)%mixtInput%externalBc) THEN       ! external BC
      CALL RFLO_GetBoundaryValues( regions(iReg) )
    ENDIF
  ENDDO
#endif

! start time stepping

  IF (global%flowType == FLOW_UNSTEADY) THEN
    IF (dTimeSystem <= 0._RFREAL) THEN
      WRITE(msg,1000) global%currentTime,global%maxTime
      CALL ErrorStop( global,ERR_DTIME_NEGATIVE,__LINE__,msg )
    ENDIF
  ELSE
    IF (dIterSystem <= 0) THEN
      WRITE(msg,1005) global%currentIter,global%maxIter
      CALL ErrorStop( global,ERR_DITER_NEGATIVE,__LINE__,msg )
    ENDIF
  ENDIF

! call time-stepping routines

#ifndef GENX
#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  timerBeg = MPI_Wtime()
#endif
#endif

  IF (global%flowType == FLOW_UNSTEADY) THEN
    IF (global%solverType == SOLV_EXPLICIT) THEN
      CALL RFLO_TimeStepping( dTimeSystem,dIterSystem,regions )
    ELSE
      IF (global%cycleType == MGCYCLE_NO) THEN
        CALL RFLO_DualTimeStepping( dTimeSystem,regions )
      ELSE
        CALL RFLO_DualMultigrid( dTimeSystem,regions )
      ENDIF
    ENDIF
  ELSE
    IF (global%cycleType == MGCYCLE_NO) THEN
      CALL RFLO_TimeStepping( dTimeSystem,dIterSystem,regions )
    ELSE
      CALL RFLO_Multigrid( dIterSystem,regions )
    ENDIF
  ENDIF

#ifndef GENX
#ifdef MPI
  timerEnd = MPI_Wtime()
  timerLoc = timerEnd - timerBeg
  CALL MPI_Reduce( timerLoc,timerGlob,1,MPI_RFREAL,MPI_SUM,MASTERPROC, &
                   global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

  IF (global%myProcid == MASTERPROC) &
    WRITE(STDOUT,1020) SOLVER_NAME,timerGlob/REAL(global%nProcAlloc), &
                       global%nProcAlloc
#endif
#endif

! finalize

  CALL DeregisterFunction( global )

! formats

1000 FORMAT('Current time is= ',1PE12.5,' but max. time is= ',E12.5)
1005 FORMAT('Current iteration is= ',I6,' but max. iteration is= ',I6)
1020 FORMAT(/,A,' Elapsed time for this run: ',1PE12.5,' sec. (average over ', &
            I5,' processors)')

END SUBROUTINE RFLO_FlowSolver

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_FlowSolver.F90,v $
! Revision 1.5  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/03/04 04:29:14  wasistho
! moved calcGridMetrics to a rocflo module
!
! Revision 1.2  2005/06/03 13:10:09  rfiedler
! Reduce "Restoring geometry" messages.  Use MASTERPROC.
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.47  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.43  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.42  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.41  2003/05/09 17:01:03  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.40  2003/04/04 21:05:00  jblazek
! Corrected bug in dumping out the solution.
!
! Revision 1.39  2003/03/06 18:29:24  jiao
! Jiri: Another silly change for GenX.
!
! Revision 1.38  2003/03/04 21:56:47  jblazek
! Corrected bug for predictor-corrector iterations.
!
! Revision 1.37  2003/02/21 22:43:00  jblazek
! Corrected timeStamp for GenX.
!
! Revision 1.36  2003/02/07 00:07:03  jblazek
! Slight change of the predictor-corrector check.
!
! Revision 1.35  2003/02/06 23:55:22  jblazek
! Added check for predictor-corrector iterations in GenX.
!
! Revision 1.34  2002/10/25 18:36:47  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.33  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.32  2002/10/18 16:49:20  jblazek
! Changed parameter lists to some GenX routines.
!
! Revision 1.31  2002/10/17 06:48:37  jiao
! Added call to Rocman at time 0..
!
! Revision 1.30  2002/10/16 18:56:14  jblazek
! Within GenX, BC data at t>0 are obtained in GetFlowSolution.
!
! Revision 1.29  2002/10/16 18:30:38  jblazek
! Within GenX, BC data at t=0 are updated in FlowSolver before calling
! the time-stepping routine.
!
! Revision 1.28  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.27  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.26  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.25  2002/04/17 22:45:20  jblazek
! Pressure forces calculated also for injection boundaries.
!
! Revision 1.24  2002/04/12 17:36:23  jblazek
! Added timer.
!
! Revision 1.23  2002/04/04 19:41:09  jblazek
! Moved dTime and dIter test from main to flowSolver.
!
! Revision 1.22  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
! Revision 1.21  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.20  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.19  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.18  2002/02/04 15:30:25  jblazek
! Added injection boundary condition.
!
! Revision 1.17  2002/02/01 22:17:38  jblazek
! Change addressing of face vectors at block boundaries.
!
! Revision 1.16  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.15  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.14  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.13  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.12  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
! Revision 1.11  2002/01/12 00:02:49  jblazek
! Added postprocessor.
!
! Revision 1.10  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.9  2002/01/10 18:21:29  jblazek
! Added iteration number and initial residual to solution file.
!
! Revision 1.8  2002/01/10 00:02:07  jblazek
! Added calculation of mixture properties.
!
! Revision 1.7  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.6  2002/01/02 16:04:20  jblazek
! Added routines to generate geometry for dummy cells.
!
! Revision 1.5  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.4  2001/12/19 23:09:22  jblazek
! Added routines to read grid and solution.
!
! Revision 1.3  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.2  2001/12/08 00:18:42  jblazek
! Added routines to read BC input file.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************







