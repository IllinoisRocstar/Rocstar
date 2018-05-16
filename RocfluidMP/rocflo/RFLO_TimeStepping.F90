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
! Purpose: integrate the governing equations in time;
!          move/regenerate the grid.
!
! Description: none.
!
! Input: dTimeSystem = total solution time (unsteady flow)
!        dIterSystem = total number of iterations (steady flow)
!        regions     = data for all grid regions
!
! Output: regions = flow variables and grid for all grid regions.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_TimeStepping.F90,v 1.25 2009/08/27 14:04:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_TimeStepping( dTimeSystem,dIterSystem,regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_MoveGridBlocks, RFLO_MoveGridGlobal, &
        RFLO_NewGrid, RungeKuttaMP, RFLO_TimeStepInviscid, &
        RFLO_TimeStepViscous, RFLO_MinimumTimeStep, ExplicitMultistage, &
        RFLO_CalcMassFlow, RFLO_CalcForces, WriteProbe, RFLO_WriteGrid, &
        RFLO_WriteSolution, WriteConvergence, RFLO_InterpolToFinerLevel, &
        RFLO_ResidualNorm, RFLO_SendBoundaryValues, RFLO_CalcThrust, &
        WriteThrust, RFLO_ComputeIntegralValues, DescaleGridSpeeds, &
        RFLO_WriteRandomState
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
  USE ModInterfacesLagrangian, ONLY : PLAG_RFLO_SetMetrics, PLAG_WriteSolution
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_WriteSolution, PEUL_SpectralRadii
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_RFLO_WriteSolution
#endif
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY : SPEC_WriteSolution
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_CalcMetrics, &
                         TURB_RFLO_WriteSolution, TURB_RFLO_RansSpectralRadii
#endif
#ifdef PERI
  USE PERI_ModHybridDES, ONLY : PERI_CoMeanCorrection 
#endif
#ifdef STATS
  USE ModStatsRoutines, ONLY : GetStatistics, StatWriteMP
#endif
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: dIterSystem

  REAL(RFREAL) :: dTimeSystem

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER, SAVE :: iter = 0

  LOGICAL :: stopExists, finished, ftermNew, residFterm, &
             doPrint, doWrite, doProbe, doThrust, moveGrid

  REAL(RFREAL), SAVE :: timePrint=0._RFREAL, timeWrite =0._RFREAL, &
                        timeProbe=0._RFREAL, timeThrust=0._RFREAL
  REAL(RFREAL) :: time, totalMass

  TYPE(t_global), POINTER :: global

  INTEGER :: original_format

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_TimeStepping',&
  'RFLO_TimeStepping.F90' )

! initialize ------------------------------------------------------------------

#ifdef GENX
  IF (global%predCorrIter) THEN
    iter = 0
    global%predCorrIter = .false.
  ENDIF
! write header for convergence history

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME
    IF (global%flowType == FLOW_STEADY) &
      WRITE(STDOUT,1010) SOLVER_NAME,SOLVER_NAME
    IF (global%flowType == FLOW_UNSTEADY) &
      WRITE(STDOUT,1015) SOLVER_NAME,SOLVER_NAME
  ENDIF
#endif

  time = 0._RFREAL

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

! time steps / iterations -----------------------------------------------------

  DO

! - min. time step (unsteady flow)

    IF (global%flowType == FLOW_UNSTEADY) THEN
      global%dtMin = 1.E+30_RFREAL
      DO iReg=1,global%nRegions
        IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
            regions(iReg)%active==ACTIVE) THEN              ! on my processor
          IF (regions(iReg)%mixtInput%flowModel == FLOW_EULER) THEN
            CALL RFLO_TimeStepInviscid( regions(iReg) )
          ELSE
            CALL RFLO_TimeStepViscous( regions(iReg) )
          ENDIF
#ifdef PEUL
          IF (global%peulUsed) &
            CALL PEUL_SpectralRadii( regions(iReg) )
#endif
#ifdef TURB
          IF (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) &
            CALL TURB_RFLO_RansSpectralRadii( regions(iReg) )
#endif
        ENDIF     ! region on this processor and active
      ENDDO       ! iReg
      CALL RFLO_MinimumTimeStep( regions )
      IF (time+global%dtMin > dTimeSystem) THEN   ! do not run over max. time
        global%dtMin = dTimeSystem - time
        finished     = .true.
      ENDIF
    ENDIF

! - move grid

    IF (global%flowType==FLOW_UNSTEADY .AND. moveGrid) THEN
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

    IF (global%flowType==FLOW_UNSTEADY .AND. moveGrid) THEN
      CALL RFLO_NewGrid( regions )
    ENDIF

! - recompute metrics of MP modules due to moving grid
#ifndef PROPONLY
    IF (global%flowType==FLOW_UNSTEADY .AND. moveGrid) THEN
#ifdef TURB
      IF (global%turbActive) CALL TURB_CalcMetrics( regions, 0 )
#endif
#ifdef PLAG
      CALL PLAG_RFLO_SetMetrics( regions )
#endif
    ENDIF

! - new solution

    IF (global%flowType == FLOW_UNSTEADY) THEN
      CALL RungeKuttaMP( regions )
      DO iReg=1,global%nRegions
        IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active,
            regions(iReg)%active==ACTIVE .AND. &           ! on my processor
            regions(iReg)%mixtInput%moveGrid) THEN         ! and moving
          CALL DescaleGridSpeeds( regions(iReg) )
        ENDIF     ! region on this processor and active
      ENDDO       ! iReg
    ELSE
      IF (global%solverType == SOLV_EXPLICIT) THEN
        CALL ExplicitMultistage( regions,ftermNew,residFterm )
      ELSE
        ! implicit scheme
        CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__ )
      ENDIF
    ENDIF

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
#endif ! PROPONLY

! - update time / iteration number

    IF (global%flowType == FLOW_UNSTEADY) THEN
      global%currentTime = global%currentTime + global%dtMin
      time = time + global%dtMin
      iter = iter + 1
    ELSE
      global%currentIter = global%currentIter + 1
      iter = iter + 1
    ENDIF

#ifndef PROPONLY
! - print/write convergence, data, probe or thrust?

    IF (global%flowType==FLOW_UNSTEADY .AND. iter==1) THEN
      timePrint  = global%timeStamp + global%printTime
      timeWrite  = global%timeStamp + global%writeTime
      timeProbe  = global%timeStamp + global%probeSaveTime
      timeThrust = global%timeStamp + global%thrustSaveTime
!      IF (global%myProcid==MASTERPROC .AND. &
!           global%verbLevel/=VERBOSE_NONE) THEN
!         WRITE(STDOUT,*) SOLVER_NAME//'S:timeWrite,timeStamp,writeTime =',timeWrite,global%timeStamp, &
!              global%writeTime
!      ENDIF

    ENDIF

    IF (global%flowType == FLOW_UNSTEADY) THEN
      doPrint  = .false.
      doWrite  = .false.
      doProbe  = .false.
      doThrust = .false.
 !     IF (global%myProcid==MASTERPROC .AND. &
 !          global%verbLevel/=VERBOSE_NONE) THEN
 !        WRITE(STDOUT,*) SOLVER_NAME//'timeWrite,timeStamp,writeTime =',timeWrite,global%timeStamp, &
 !             global%writeTime
 !        WRITE(STDOUT,*) SOLVER_NAME//'currentTime =',global%currentTime
 !     ENDIF
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
    ELSE
      doPrint  = (MOD(global%currentIter,global%printIter     ) == 0)
      doWrite  = (MOD(global%currentIter,global%writeIter     ) == 0)
      doProbe  = (MOD(global%currentIter,global%probeSaveIter ) == 0)
      doThrust = (MOD(global%currentIter,global%thrustSaveIter) == 0)
    ENDIF
#endif ! PROPONLY
! - check for stop file

#ifndef GENX
    INQUIRE(file="STOP",exist=stopExists)
    IF (stopExists) global%stopRun = 1.1_RFREAL
#endif

! - check for end of time stepping

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (time>=dTimeSystem) finished = .true.
    ELSE
#ifndef PROPONLY
      CALL RFLO_ResidualNorm( regions )
      IF (iter==dIterSystem .OR. &
          global%residual/global%resInit<=global%resTol) finished = .true.
#endif
    ENDIF

! - compute forces, mass flow & thrust;
! - store probe data, thrust and flow solution

    global%forceX      = 0._RFREAL
    global%forceY      = 0._RFREAL
    global%forceZ      = 0._RFREAL
    global%massIn      = 0._RFREAL
    global%massOut     = 0._RFREAL
    global%thrustMom   = 0._RFREAL
    global%thrustPress = 0._RFREAL
    totalMass          = 0._RFREAL

#ifndef PROPONLY 
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor

! ----- store probe, aero-coeffs. and thrust data
#ifdef GENX
        IF (global%nProbes>0 .AND. (doProbe .eqv. .true.)) &
          CALL WriteProbe( regions,iReg )

        IF (global%aeroCoeffs==ACTIVE .AND. (doProbe .eqv. .true.)) &
          CALL RFLO_ComputePatchForceMomCo( regions(iReg) )

        IF (global%thrustType/=THRUST_NONE .AND. (doThrust .eqv. .true.)) &
          CALL RFLO_CalcThrust( regions(iReg) )
#else
        IF (global%nProbes>0 .AND. ((doProbe.eqv..true.).OR.(finished.eqv..true.))) &
          CALL WriteProbe( regions,iReg )

        IF (global%aeroCoeffs==ACTIVE .AND. (doProbe.OR.finished)) &
          CALL RFLO_ComputePatchForceMomCo( regions(iReg) )

        IF (global%thrustType/=THRUST_NONE .AND. (doThrust.OR.finished)) &
          CALL RFLO_CalcThrust( regions(iReg) )
#endif

! ----- get forces, mass flow
        IF ((doPrint .eqv. .true.) .OR. (finished .eqv. .true.)) THEN
          CALL RFLO_CalcMassFlow( regions(iReg) )
          CALL RFLO_CalcForces( regions(iReg) )
        ENDIF
      ENDIF     ! region on this processor and active
    ENDDO       ! iReg

! - compute and write global aero-coeffs. to file
#ifdef GENX
    IF (global%aeroCoeffs==ACTIVE .AND. (doProbe .eqv. .true.)) THEN
      CALL RFLO_ComputeIntegralForceMomCo( global )
      CALL RFLO_WriteIntegralForceMomCo( global )
    ENDIF
#else
    IF (global%aeroCoeffs==ACTIVE .AND. ((doProbe .eqv. .true.).OR.(finished.eqv..true.))) THEN
      CALL RFLO_ComputeIntegralForceMomCo( global )
      CALL RFLO_WriteIntegralForceMomCo( global )
    ENDIF
#endif

! - write thrust to file
#ifdef GENX
    IF (global%thrustType/=THRUST_NONE .AND. (doThrust .eqv. .true.)) &
      CALL WriteThrust( global )
#else
    IF (global%thrustType/=THRUST_NONE .AND. ((doThrust .eqv. .true.).OR.(finished .eqv. .true.))) &
      CALL WriteThrust( global )
#endif

! - write convergence history (file & screen)

    IF ((doPrint .eqv. .true.) .OR. (finished .eqv. .true.)) THEN
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

    IF ((doWrite .eqv. .true.) .OR. (finished .eqv. .true.)) THEN
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
          WRITE(STDOUT,'(A)')   SOLVER_NAME//'   - patch ac'
        ENDIF
        CALL RFLO_WritePatchAeroCoeffs( regions )
      ENDIF
#ifdef STATS
      IF (global%doStat==ACTIVE) THEN
        IF (global%myProcid==MASTERPROC .AND. &
            global%verbLevel/=VERBOSE_NONE) THEN
          WRITE(STDOUT,'(/,A)') SOLVER_NAME,' Saving statistics ...'
        ENDIF
        CALL StatWriteMP( regions )
      ENDIF
#endif
#endif ! ndef GENX

#ifdef GENX
#ifdef NATIVE_MP_IO
    IF (doWrite .eqv. .true.) THEN
      original_format = global%solutFormat
      global%gridFormat  = FORMAT_ASCII
      global%solutFormat = FORMAT_ASCII
#ifdef PLAG
        IF (global%myProcid==MASTERPROC .AND. &
            global%verbLevel/=VERBOSE_NONE) THEN
          WRITE(STDOUT,'(/,A)') SOLVER_NAME,' Saving Lagrangian particle solution ...'
        ENDIF
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
      global%gridFormat  = original_format
      global%solutFormat = original_format
   ENDIF
#endif ! NATIVE_MP_IO
#endif ! GENX

#ifndef GENX
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
#endif ! Not GENX
#endif ! NOT PROPONLY
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

! - check if to proceed to next finer grid

    IF (global%flowType==FLOW_STEADY .AND. MOD(iter,global%refineIter)==0) THEN
      DO iReg=1,global%nRegions
        regions(iReg)%currLevel = regions(iReg)%currLevel - 1
        CALL RFLO_InterpolToFinerLevel( regions(iReg) )
      ENDDO
    ENDIF

  ENDDO  ! loop over time / iter

! formats

1010 FORMAT(A,'  iter',4X,'res-norm',5X,'force-x',6X,'force-y',6X,'force-z', &
            6X,'mass-in',6X,'mass-out',/,A,1X,84('-'))
1015 FORMAT(A,'  time',10X,'delta-t',6X,'force-x',6X,'force-y',6X,'force-z', &
            6X,'mass-in',6X,'mass-out'/,A,1X,90('-'))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_TimeStepping

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_TimeStepping.F90,v $
! Revision 1.25  2009/08/27 14:04:52  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.24  2009/08/12 04:15:58  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.23  2009/04/07 14:53:57  mtcampbe
! Switch to FORMAT_ASCII for native mp io
!
! Revision 1.22  2009/03/05 13:00:23  mtcampbe
! Added NATIVEMPIO and PROPONLY options for Rocflo
!
! Revision 1.21  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.20  2008/11/19 22:17:39  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.19  2006/03/24 23:29:23  wasistho
! added forceMomentCoeffs computation and output
!
! Revision 1.18  2006/03/22 03:03:25  wasistho
! added call to RFLO_WritePatchAeroCoeffs
!
! Revision 1.17  2006/03/02 01:27:58  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.16  2006/02/11 03:55:44  wasistho
! added MOVEGRID_EPDE
!
! Revision 1.15  2006/02/01 20:02:18  wasistho
! added WriteRestartInfo
!
! Revision 1.14  2005/11/11 07:32:47  wasistho
! removed commented call to moveGridGlobal
!
! Revision 1.13  2005/11/11 07:27:05  wasistho
! removed obsolete file
!
! Revision 1.12  2005/10/27 05:12:41  wasistho
! added moveGridFoms
!
! Revision 1.11  2005/06/02 03:21:17  wasistho
! shuffle MoveGridVms with MoveGridFrame
!
! Revision 1.10  2005/05/28 08:08:18  wasistho
! added moveGridFrame
!
! Revision 1.9  2005/05/21 01:43:17  wasistho
! added statWriteMP
!
! Revision 1.8  2005/05/21 00:18:33  wasistho
! added moveGridVms
!
! Revision 1.7  2005/04/17 05:09:58  wasistho
! mv guided DES treatment to before statistics
!
! Revision 1.6  2005/03/11 04:24:23  wasistho
! added mean correction for PERI flows with hybrid DES
!
! Revision 1.5  2005/02/26 04:05:29  wasistho
! added RFLO_ComputeIntegralValues
!
! Revision 1.4  2005/02/16 14:43:08  fnajjar
! Included PLAG call to communicate statistics buffers during IO stage
!
! Revision 1.3  2005/01/08 20:38:39  fnajjar
! Added PLAG_RFLO_WriteStat call for IO
!
! Revision 1.2  2004/12/28 20:27:52  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.58  2004/11/29 17:15:57  wasistho
! use ModInterfacesStatistics
!
! Revision 1.57  2004/09/23 03:49:29  wasistho
! changed RADI_WriteSol.. to RADI_RFLO_WriteSol..
!
! Revision 1.56  2004/03/11 03:31:09  wasistho
! changed rocturb nomenclature
!
! Revision 1.55  2004/03/05 22:09:02  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.54  2004/02/11 03:27:57  wasistho
! get rid of argument numVar in TURB_WriteSolution
!
! Revision 1.53  2004/02/07 01:12:43  wasistho
! modified TURB_WriteSolution
!
! Revision 1.51  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.50  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.46  2003/11/12 21:21:06  fnajjar
! Added Corner-Edge cells routine to communicate metrics for PLAG
!
! Revision 1.45  2003/10/15 03:38:41  wasistho
! added call to turbulence spectralRadii routine
!
! Revision 1.44  2003/10/03 20:18:57  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.43  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.42  2003/09/26 21:44:28  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.41  2003/08/28 20:35:34  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.40  2003/08/15 21:15:17  jblazek
! Corrected bug in output of thrust in case of GENX.
!
! Revision 1.39  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.38  2003/08/01 22:14:08  wasistho
! radiWrite/turbWrite to radiActive/turbActive
!
! Revision 1.37  2003/07/22 02:57:42  wasistho
! prepare more accurate rocturb restart
!
! Revision 1.36  2003/07/17 01:03:28  wasistho
! initial activation rocrad
!
! Revision 1.35  2003/07/08 21:21:37  jblazek
! Modified start up procedure for dual-time stepping.
!
! Revision 1.34  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.33  2003/06/02 17:12:01  jblazek
! Added computation of thrust.
!
! Revision 1.32  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.31  2003/04/10 01:22:41  jblazek
! Got rid of pRegion in ViscousFluxesMP.
!
! Revision 1.30  2003/04/04 21:05:00  jblazek
! Corrected bug in dumping out the solution.
!
! Revision 1.29  2003/03/28 19:35:06  fnajjar
! include RungeKuttaMP routine
!
! Revision 1.28  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.27  2003/02/11 22:53:19  jferry
! Initial import of Rocsmoke
!
! Revision 1.26  2003/02/05 21:07:30  jblazek
! Coordinated stop of a run works now for MPI.
!
! Revision 1.25  2003/01/23 17:48:53  jblazek
! Changed algorithm to dump convergence, solution and probe data.
!
! Revision 1.24  2003/01/10 17:58:43  jblazek
! Added missing explicit interfaces.
!
! Revision 1.23  2002/11/02 01:58:14  wasistho
! Added TURB statistics
!
! Revision 1.22  2002/10/02 22:21:59  jiao
! Debugged GenX restart.
!
! Revision 1.21  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.20  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.19  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.18  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.17  2002/07/18 22:51:57  jblazek
! Reduce time step to match max. physical time.
!
! Revision 1.16  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.15  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.14  2002/06/18 00:33:57  wasistho
! Added prefix SOLVER NAME to satistics STDOutput
!
! Revision 1.13  2002/06/14 21:38:45  wasistho
! Added time avg statistics
!
! Revision 1.12  2002/06/12 21:56:29  jblazek
! Added read/write solution for physical modules.
!
! Revision 1.11  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.10  2002/04/11 21:10:27  jblazek
! Set correct time when writing grid only for Tecplot.
!
! Revision 1.9  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.8  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.7  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.6  2002/02/01 21:04:26  jblazek
! Streamlined time stepping routine.
!
! Revision 1.5  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/23 23:36:37  jblazek
! All blocks passed to time integration routines.
!
! Revision 1.2  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
!******************************************************************************







