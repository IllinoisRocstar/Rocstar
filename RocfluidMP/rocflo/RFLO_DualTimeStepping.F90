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
! Purpose: integrate the governing equations in time using
!          dual time-stepping; move/regenerate the grid.
!
! Description: none.
!
! Input: dTimeSystem = total solution time
!        regions     = data for all grid regions
!
! Output: regions = flow variables and grid for all grid regions.
!
! Notes: scheme also applicable for small time steps.
!
!******************************************************************************
!
! $Id: RFLO_DualTimeStepping.F90,v 1.30 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_DualTimeStepping( dTimeSystem,regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_MoveGridBlocks, RFLO_MoveGridGlobal, &
        RFLO_NewGrid, RFLO_DualTstInit, RFLO_DualTstPredict, RFLO_DualTstSterm,&
        RFLO_DualTstShift, ExplicitMultistage, RFLO_CalcMassFlow, &
        RFLO_CalcForces, WriteProbe, RFLO_WriteGrid, RFLO_WriteSolution, &
        WriteConvergence, RFLO_ResidualNorm, RFLO_SendBoundaryValues, &
        RFLO_CalcThrust, WriteThrust, RFLO_ComputeIntegralValues, &
        RFLO_TimeStepInviscid, RFLO_TimeStepViscous, RFLO_MinimumTimeStep, &
        RungeKuttaMP, RFLO_SetMstageCoeffs, RFLO_WriteRandomState

  USE RFLO_ModForcesMoments,    ONLY : RFLO_ComputePatchForceMomCo, &
                                       RFLO_ComputeIntegralForceMomCo, &
                                       RFLO_WriteIntegralForceMomCo
  USE RFLO_ModMoveGridFrame,    ONLY : RFLO_MoveGridFrame
  USE RFLO_ModMoveGridElliptGlo,ONLY : RFLO_MoveGridElliptGlo
  USE RFLO_ModMoveGridElliptFra,ONLY : RFLO_MoveGridElliptFra
  USE RFLO_ModVolMeshSmoothing, ONLY : RFLO_MoveGridVms
  USE RFLO_ModPatchAeroCoeffs,  ONLY : RFLO_WritePatchAeroCoeffs
  USE RFLO_ModRestartInfo,      ONLY : RFLO_WriteRestartInfo
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_WriteSolution
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_WriteSolution
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_RFLO_WriteSolution
#endif
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY : SPEC_WriteSolution
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_CalcMetrics, &
                                      TURB_RFLO_WriteSolution
#endif
#ifdef PERI
  USE PERI_ModHybridDES, ONLY : PERI_CoMeanCorrection 
#endif
#ifdef STATS
  USE ModInterfaces, ONLY : RFLO_WriteStat
  USE ModStatsRoutines, ONLY : GetStatistics
  USE RFLO_ModStatsBoundaryConditions, ONLY : RFLO_StatBoundaryConditionsSet
#endif
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: dTimeSystem

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, subIter

! ... local variables
  INTEGER, SAVE :: iter=0
  INTEGER, ALLOCATABLE :: timeScheme(:)

  LOGICAL :: stopExists, finished, ftermNew, residFterm, &
             doPrint, doWrite, doProbe, doThrust, moveGrid

  REAL(RFREAL), SAVE :: timePrint=0._RFREAL, timeWrite =0._RFREAL, &
                        timeProbe=0._RFREAL, timeThrust=0._RFREAL
  REAL(RFREAL) :: time, subTime, timeBc, alphaBc, residRat
  REAL(RFREAL), ALLOCATABLE :: cfl(:), smoocf(:)
  REAL(RFREAL), POINTER :: cv(:,:), cvn(:,:), cvn1(:,:), cvn2(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_DualTimeStepping',&
  'RFLO_DualTimeStepping.F90' )

! initialize parameters -------------------------------------------------------

#ifdef GENX
  IF (global%predCorrIter) THEN
    iter = 0
    global%predCorrIter = .false.
  ENDIF
#endif

  time         = 0._RFREAL

  IF (global%dtFixed) THEN
    IF (global%currentTime <= EPSILON( 1._RFREAL )) THEN  ! more accurate
      global%dtMin = global%dtImposed                     ! but slower
    ENDIF
  ELSE
    global%dtMin = global%dtImposed                       ! for more speed
  ENDIF

  finished       = .false.   ! run not finished yet
  global%stopRun = 0._RFREAL

! no multigrid here

  ftermNew   = .false.       ! no new forcing term
  residFterm = .false.       ! do not add forcing term to residual

! got some moving grid?

  moveGrid = .false.
  DO iReg=1,global%nRegions
    IF (regions(iReg)%mixtInput%moveGrid) moveGrid = .true.
  ENDDO

! initialize solutions --------------------------------------------------------

  IF (iter == 0) THEN

    global%dualTstSource = .false.    ! ordinary Runge-Kutta here

! - solution at level n-2

    CALL RFLO_DualTstInit( regions,2 )

! - store scheme, CFL, epsIRS (overwritten for Runge-Kutta!)

    ALLOCATE( timeScheme(global%nRegions) )
    ALLOCATE( cfl       (global%nRegions) )
    ALLOCATE( smoocf    (global%nRegions) )
    DO iReg=1,global%nRegions                            ! commented since all
!      IF (regions(iReg)%procid==global%myProcid .AND. & ! procs need trk from
!          regions(iReg)%active==ACTIVE) THEN            ! ireg=1 (see ExpMS)
        timeScheme(iReg) = regions(iReg)%mixtInput%timeScheme
        cfl(iReg)        = regions(iReg)%mixtInput%cfl
        smoocf(iReg)     = regions(iReg)%mixtInput%smoocf
        regions(iReg)%mixtInput%timeScheme = TST_STD4RK
        regions(iReg)%mixtInput%cfl        = 3._RFREAL
        CALL RFLO_SetMstageCoeffs( global,regions(iReg)%mixtInput,global%nrkSteps )
!      ENDIF     ! region on this processor and active
    ENDDO       ! iReg

! - solutions at levels n-1 and n

    subTime = 0._RFREAL

    DO subIter=1,0,-1
      global%dtMin = 1.E+30_RFREAL
      DO iReg=1,global%nRegions
        IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
            regions(iReg)%active==ACTIVE) THEN              ! on my processor
          IF (regions(iReg)%mixtInput%flowModel == FLOW_EULER) THEN
            CALL RFLO_TimeStepInviscid( regions(iReg) )
          ELSE
            CALL RFLO_TimeStepViscous( regions(iReg) )
          ENDIF
        ENDIF     ! region on this processor and active
      ENDDO       ! iReg
      CALL RFLO_MinimumTimeStep( regions )
      IF (time+global%dtMin > dTimeSystem) THEN   ! do not run over max. time
        global%dtMin = dTimeSystem - time
      ENDIF

      IF (moveGrid) CALL RFLO_MoveGridBlocks( regions )

      CALL RungeKuttaMP( regions )
      CALL RFLO_DualTstInit( regions,subIter )

      global%currentTime = global%currentTime + global%dtMin
      time               = time + global%dtMin
      subTime            = subTime + global%dtMin
    ENDDO  ! subIter

! - restore scheme, CFL, epsIRS again

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor
        regions(iReg)%mixtInput%timeScheme = timeScheme(iReg)
        regions(iReg)%mixtInput%cfl        = cfl(iReg)
        regions(iReg)%mixtInput%smoocf     = smoocf(iReg)
        CALL RFLO_SetMstageCoeffs( global,regions(iReg)%mixtInput,global%nrkSteps )
      ENDIF     ! region on this processor and active
    ENDDO       ! iReg
    DEALLOCATE( timeScheme )
    DEALLOCATE( cfl        )
    DEALLOCATE( smoocf     )

! - substract subTime from the first implicit time step

    IF (global%dtImposed <= subTime) THEN
      global%dtMin = 0._RFREAL
      DO
        global%dtMin = global%dtMin + global%dtImposed
        IF (global%dtMin > subTime) EXIT
      ENDDO
      global%dtMin = global%dtMin - subTime
    ELSE
      global%dtMin = global%dtImposed - subTime
    ENDIF
    global%dualTstSource = .true.

  ENDIF   ! iter=0

! time steps ------------------------------------------------------------------

  DO

    IF (global%dtFixed) THEN

! --- Fixed dt for speed (default)
      IF (iter > 0) global%dtMin = global%dtImposed

    ELSE

! --- Adjustable dt for higher accuracy, but slower
      residRat = global%residual/(global%resInit*global%tolSubIter)
      IF (iter == 0)  THEN
        global%dtMin = 0.01_RFREAL * global%dtMin
      ELSEIF (iter == 1) THEN
        IF (global%currentTime < global%dtImposed) & ! even more accrt bt slwer
          global%dtMin = global%dtImposed
      ELSEIF (iter > 1) THEN
        IF (residRat > 1._RFREAL) THEN
          global%dtMin = 0.75_RFREAL*global%dtMin
        ELSE
          IF (global%dtMin < global%dtImposed) &
            global%dtMin = MIN( 1.25_RFREAL*global%dtMin, global%dtImposed )
        ENDIF  ! residRat
      ENDIF    ! iter

    ENDIF      ! dtFixed

    IF (time+global%dtMin > dTimeSystem) THEN   ! do not run over max. time
      global%dtMin = dTimeSystem - time
      finished     = .true.
    ENDIF

! - move grid

    IF (moveGrid) THEN
      IF (global%moveGridScheme == MOVEGRID_BLOCKS) THEN
        CALL RFLO_MoveGridBlocks( regions )
      ELSEIF (global%moveGridScheme == MOVEGRID_GLOBAL) THEN
        CALL RFLO_MoveGridGlobal( regions )
      ELSEIF (global%moveGridScheme == MOVEGRID_FRAME .OR. &
              global%moveGridScheme == MOVEGRID_FOMS) THEN
        CALL RFLO_MoveGridFrame( regions )
      ELSEIF (global%moveGridScheme == MOVEGRID_ELGLOBAL) THEN
        CALL RFLO_MoveGridElliptGlo( regions )
      ELSEIF (global%moveGridScheme == MOVEGRID_ELFRAME) THEN
        CALL RFLO_MoveGridElliptFra( regions )
      ELSE
        CALL RFLO_MoveGridVms( regions )
      ENDIF
    ENDIF

! - new grid needed?

    IF (moveGrid) THEN
      CALL RFLO_NewGrid( regions )
    ENDIF

! - recompute metrics of MP modules due to moving grid

    IF (moveGrid) THEN
#ifdef TURB
      IF (global%turbActive) CALL TURB_CalcMetrics( regions, 0 )
#endif
    ENDIF

! - compute source term for implicit scheme; guess start solution

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor
        CALL RFLO_DualTstSterm( regions(iReg) )
        IF (global%predictSol .OR. iter==0) &
          CALL RFLO_DualTstPredict( regions(iReg) )
      ENDIF     ! region on this processor and active
    ENDDO       ! iReg

! - new solution (do subiterations)

    subIter         = 0
    global%flowType = FLOW_STEADY

    DO
      subIter            = subIter + 1
      global%currentIter = subIter
      CALL ExplicitMultistage( regions,ftermNew,residFterm )
      CALL RFLO_ResidualNorm( regions )
      IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) THEN
        WRITE(STDOUT,'(A,1X,2I6,1PE13.4,E13.4)') &
          SOLVER_NAME,iter,subIter,global%residual/global%resInit,global%residual
      ENDIF
      IF (subIter==global%maxSubIter .OR. &
          global%residual/global%resInit<=global%tolSubIter .OR. &
          global%residual<100._RFREAL*EPSILON(1.0_RFREAL)) EXIT
    ENDDO
    IF (global%myProcid == MASTERPROC .AND. global%verbLevel>=VERBOSE_LOW) THEN
      WRITE(STDOUT,'(A,1X,I6,1PE13.4,E13.4)') &
        SOLVER_NAME,subIter,global%residual/global%resInit,global%residual
    ENDIF

    global%flowType = FLOW_UNSTEADY

! - shift time levels

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor
        CALL RFLO_DualTstShift( regions(iReg) )
      ENDIF     ! region on this processor and active
    ENDDO       ! iReg

#ifdef PERI
! - hybrid (guided) DES

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor
        IF ( regions(iReg)%periInput%flowKind /= OFF ) THEN
          CALL PERI_CoMeanCorrection( regions(iReg) )
        ENDIF ! flowKind 
      ENDIF   ! region on this processor and active
    ENDDO     ! iReg
#endif

#ifdef STATS
! - get statistics

    CALL GetStatistics( regions )
#endif

! - update physical time

    global%currentTime = global%currentTime + global%dtMin
    time = time + global%dtMin
    iter = iter + 1

! - print/write convergence, data, probe or thrust?

    IF (iter == 1) THEN
      timePrint  = global%timeStamp + global%printTime
      timeWrite  = global%timeStamp + global%writeTime
      timeProbe  = global%timeStamp + global%probeSaveTime
      timeThrust = global%timeStamp + global%thrustSaveTime
    ENDIF

    doPrint  = .false.
    doWrite  = .false.
    doProbe  = .false.
    doThrust = .false.
    IF (ABS(timePrint-global%currentTime)<global%dtMin/10._RFREAL .OR. &
        timePrint<global%currentTime .OR. iter==1) THEN
      doPrint = .true.
      IF (iter > 1) timePrint = timePrint + global%printTime
    ENDIF
    IF (ABS(timeWrite-global%currentTime)<global%dtMin/10._RFREAL .OR. &
        timeWrite<global%currentTime) THEN
      doWrite   = .true.
      timeWrite = timeWrite + global%writeTime
    ENDIF
    IF (ABS(timeProbe-global%currentTime)<global%dtMin/10._RFREAL .OR. &
        timeProbe<global%currentTime .OR. iter==1) THEN
      doProbe = .true.
      IF (iter > 1) timeProbe = timeProbe + global%probeSaveTime
    ENDIF
    IF (ABS(timeThrust-global%currentTime)<global%dtMin/10._RFREAL .OR. &
        timeThrust<global%currentTime .OR. iter==1) THEN
      doThrust = .true.
      IF (iter > 1) timeThrust = timeThrust + global%thrustSaveTime
    ENDIF

! - check for stop file

#ifndef GENX
    INQUIRE(file="STOP",exist=stopExists)
    IF (stopExists) global%stopRun = 1.1_RFREAL
#endif

! - check for end of time stepping

    IF (time >= dTimeSystem) finished = .true.

! - compute forces, mass flow & thrust;
! - store probe data, thrust and flow solution

    global%forceX      = 0._RFREAL
    global%forceY      = 0._RFREAL
    global%forceZ      = 0._RFREAL
    global%massIn      = 0._RFREAL
    global%massOut     = 0._RFREAL
    global%thrustMom   = 0._RFREAL
    global%thrustPress = 0._RFREAL

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor

! ----- store probe, aero-coeffs. and thrust data
#ifdef GENX
        IF (global%nProbes>0 .AND. doProbe) &
          CALL WriteProbe( regions,iReg )

        IF (global%aeroCoeffs==ACTIVE .AND. doProbe) &
          CALL RFLO_ComputePatchForceMomCo( regions(iReg) )

        IF (global%thrustType/=THRUST_NONE .AND. doThrust) &
          CALL RFLO_CalcThrust( regions(iReg) )
#else
        IF (global%nProbes>0 .AND. (doProbe.OR.finished)) &
          CALL WriteProbe( regions,iReg )

        IF (global%aeroCoeffs==ACTIVE .AND. (doProbe.OR.finished)) &
          CALL RFLO_ComputePatchForceMomCo( regions(iReg) )

        IF (global%thrustType/=THRUST_NONE .AND. (doThrust.OR.finished)) &
          CALL RFLO_CalcThrust( regions(iReg) )
#endif

! ----- get forces, mass flow
        IF (doPrint .OR. finished) THEN
          CALL RFLO_CalcMassFlow( regions(iReg) )
          CALL RFLO_CalcForces( regions(iReg) )
        ENDIF
      ENDIF     ! region on this processor and active
    ENDDO       ! iReg

! - compute and write global aero-coeffs. to file
#ifdef GENX
    IF (global%aeroCoeffs==ACTIVE .AND. doProbe) THEN
      CALL RFLO_ComputeIntegralForceMomCo( global )
      CALL RFLO_WriteIntegralForceMomCo( global )
    ENDIF
#else
    IF (global%aeroCoeffs==ACTIVE .AND. (doProbe.OR.finished)) THEN
      CALL RFLO_ComputeIntegralForceMomCo( global )
      CALL RFLO_WriteIntegralForceMomCo( global )
    ENDIF
#endif

! - write thrust to file
#ifdef GENX
    IF (global%thrustType/=THRUST_NONE .AND. doThrust) &
      CALL WriteThrust( global )
#else
    IF (global%thrustType/=THRUST_NONE .AND. (doThrust.OR.finished)) &
      CALL WriteThrust( global )
    IF (moveGrid) CALL RFLO_ComputeIntegralValues( regions )
#endif

! - write convergence history (file & screen)

    IF (doPrint .OR. finished) THEN
#ifdef MASS
      CALL RFLO_ComputeIntegralValues( regions )
      WRITE(STDOUT,*) 'Total mass = ',global%totalMass
#endif
      CALL WriteConvergence( global )
#ifndef GENX
      IF (global%stopRun > 1._RFREAL) THEN
        finished = .true.
      ENDIF
#endif
    ENDIF

#ifndef GENX
! - store flow field (and grid if moving)

    IF (doWrite .OR. finished) THEN
      IF (moveGrid) THEN
        IF (global%myProcid==MASTERPROC .AND. &
            global%verbLevel/=VERBOSE_NONE) THEN
          WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Saving grid ...'
        ENDIF
        CALL RFLO_WriteGrid( regions )
      ENDIF
      IF (global%myProcid==MASTERPROC .AND. &
          global%verbLevel/=VERBOSE_NONE) THEN
        WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Saving flow solution ...'
      ENDIF
      IF (global%myProcid==MASTERPROC .AND. &
          global%verbLevel>=VERBOSE_HIGH) THEN
        WRITE(STDOUT,'(A)')   SOLVER_NAME//'   - mixture'
      ENDIF
      CALL RFLO_WriteSolution( regions )
      CALL RFLO_WriteRandomState( regions )

      IF (global%aeroCoeffs == ACTIVE) THEN
        IF (global%myProcid==MASTERPROC .AND. &
            global%verbLevel>=VERBOSE_HIGH) THEN
          WRITE(STDOUT,'(A)') SOLVER_NAME//'   - patch ac'
        ENDIF
        CALL RFLO_WritePatchAeroCoeffs( regions )
      ENDIF
#ifdef STATS
      IF (global%doStat==ACTIVE) THEN
        IF (global%myProcid==MASTERPROC .AND. &
            global%verbLevel/=VERBOSE_NONE) THEN
          WRITE(STDOUT,'(/,A)') SOLVER_NAME,' Saving statistics ...'
        ENDIF
        CALL RFLO_StatBoundaryConditionsSet( regions )
        CALL RFLO_WriteStat( regions )
      ENDIF
#endif
#ifdef PLAG
      CALL PLAG_WriteSolution( regions )
#endif
#ifdef PEUL
      CALL PEUL_WriteSolution( regions )
#endif
#ifdef RADI
      IF (global%radiActive) THEN
        IF (global%myProcid==MASTERPROC .AND. &
            global%verbLevel/=VERBOSE_NONE) THEN
          WRITE(STDOUT,'(/,A)') SOLVER_NAME,' Saving radiation solution ...'
        ENDIF
        CALL RADI_RFLO_WriteSolution( regions )
      ENDIF
#endif
#ifdef SPEC
      CALL SPEC_WriteSolution( regions )
#endif
#ifdef TURB
      IF (global%turbActive) THEN
        IF (global%myProcid==MASTERPROC .AND. &
            global%verbLevel/=VERBOSE_NONE) THEN
          WRITE(STDOUT,'(/,A)') SOLVER_NAME,' Saving turbulence solution ...'
        ENDIF
        CALL TURB_RFLO_WriteSolution( regions )
      ENDIF
#endif
      CALL RFLO_WriteRestartInfo( global )

    ENDIF ! doWrite
#endif

! - run finished?

    IF (finished) THEN
      DO iReg=1,global%nRegions
        IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
            regions(iReg)%active==ACTIVE) THEN             ! on my processor
          CALL RFLO_SendBoundaryValues( regions(iReg),.false. )
        ENDIF     ! region on this processor and active
      ENDDO       ! iReg
      global%timeStamp = global%currentTime
      EXIT
    ENDIF

  ENDDO  ! loop over physical time

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_DualTimeStepping

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_DualTimeStepping.F90,v $
! Revision 1.30  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.29  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.28  2006/05/09 23:32:06  wasistho
! added fixed-or-free dt options in implicit dtst
!
! Revision 1.27  2006/05/05 17:35:27  wasistho
! commented back region-split for parallel, all procs need trk
!
! Revision 1.26  2006/04/15 00:23:28  wasistho
! set constant dt as default timestepping
!
! Revision 1.25  2006/03/24 23:29:01  wasistho
! added forceMomentCoeffs computation and output
!
! Revision 1.24  2006/03/22 03:03:32  wasistho
! added call to RFLO_WritePatchAeroCoeffs
!
! Revision 1.23  2006/03/06 08:07:45  wasistho
! set to more agressive stepping
!
! Revision 1.22  2006/03/02 01:26:40  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.21  2006/02/11 03:55:58  wasistho
! added MOVEGRID_EPDE
!
! Revision 1.20  2006/02/01 20:10:13  wasistho
! added WriteRestartInfo
!
! Revision 1.19  2006/01/28 03:12:45  wasistho
! set adjustable timestep to more accurate option
!
! Revision 1.18  2006/01/27 07:32:04  wasistho
! adjustable timestep for stability
!
! Revision 1.17  2006/01/26 08:25:23  wasistho
! bug fixed: invert the dtMin adjuster
!
! Revision 1.16  2006/01/26 08:10:58  wasistho
! relate dtMin to previous residual
!
! Revision 1.15  2006/01/13 00:08:28  wasistho
! bound global%dtMin within max.time
!
! Revision 1.14  2005/11/07 19:49:11  wasistho
! added MOVEGRID_FOMS
!
! Revision 1.13  2005/06/16 03:52:40  wasistho
! activated RFLO_ModStatsBc
!
! Revision 1.12  2005/06/02 03:21:25  wasistho
! shuffle MoveGridVms with MoveGridFrame
!
! Revision 1.11  2005/05/28 08:08:28  wasistho
! added moveGridFrame
!
! Revision 1.10  2005/05/21 07:08:18  wasistho
! backout RFLO_ModStatsBoundaryConditions temporarily
!
! Revision 1.9  2005/05/21 01:43:54  wasistho
! added rflo_StatBcSet
!
! Revision 1.8  2005/05/21 00:18:52  wasistho
! added moveGridVms
!
! Revision 1.7  2005/04/17 05:09:50  wasistho
! mv guided DES treatment to before statistics
!
! Revision 1.6  2005/03/10 02:03:23  wasistho
! test calling PERI_coMeanCorrection from DualTsT i.o.EMS
!
! Revision 1.5  2005/02/26 04:05:38  wasistho
! added RFLO_ComputeIntegralValues
!
! Revision 1.4  2004/12/28 20:27:12  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.3  2004/12/15 09:20:14  wasistho
! fixed dual tst for Rocstar
!
! Revision 1.2  2004/12/09 22:16:55  wasistho
! added data turbulence Metric computation
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.17  2004/11/29 17:15:38  wasistho
! use ModInterfacesStatistics
!
! Revision 1.16  2004/11/17 16:30:25  haselbac
! Adapted interface of RFLO_SetMStageCoeffs
!
! Revision 1.15  2004/09/23 03:50:04  wasistho
! changed RADI_WriteSol.. to RADI_RFLO_WriteSol..
!
! Revision 1.14  2004/03/11 03:31:33  wasistho
! changed rocturb nomenclature
!
! Revision 1.13  2004/02/11 03:28:09  wasistho
! get rid of argument numVar in TURB_WriteSolution
!
! Revision 1.12  2004/02/07 01:12:56  wasistho
! modified TURB_WriteSolution
!
! Revision 1.11  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.10  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/09/26 21:44:28  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.5  2003/08/28 20:35:46  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.4  2003/08/15 21:15:17  jblazek
! Corrected bug in output of thrust in case of GENX.
!
! Revision 1.3  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.2  2003/07/08 21:21:37  jblazek
! Modified start up procedure for dual-time stepping.
!
! Revision 1.1  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
!******************************************************************************







