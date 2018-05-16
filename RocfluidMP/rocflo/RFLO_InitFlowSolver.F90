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
! Purpose: initialize all global variables and regions.
!
! Description: none.
!
! Input: casename  = name of the case
!        verbLevel = verbosity level (VERBOSE_NONE/LOW/HIGH).
!
! Output: global  = global variables
!         regions = dimensions and initial values for all regions
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_InitFlowSolver.F90,v 1.18 2009/08/12 04:15:58 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef GENX
SUBROUTINE RFLO_InitFlowSolver( globalGenx,initialTime,communicator, &
                                genxHandle,inSurf,inVol,obtain_attribute )
#else
SUBROUTINE RFLO_InitFlowSolver( casename,verbLevel,global,regions )
#endif

  USE ModDataTypes
#ifdef GENX
  USE ModRocstar, ONLY       : t_globalGenx
  USE ModInterfaces, ONLY : RFLO_InitGenxInterface
#endif
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_ReadRegionTopology, RFLO_FindSourcePatches, &
        RFLO_GetUserInput, RFLO_CheckRegionFaces, RFLO_FindSourceRegions, &
        RFLO_DoMemoryAllocation, RFLO_GetGeometry, &
        RFLO_InitGridProcedures, RFLO_GetFlowSolution, &
        RFLO_OpenConverFile, RFLO_OpenProbeFile, RFLO_CheckMinimumCells, &
        RFLO_SendBoundaryValues, RFLO_RandomInit, RFLO_FindThrustPatches, &
        RFLO_OpenThrustFile, BuildVersionString
  USE RFLO_ModForcesMoments, ONLY: RFLO_OpenForceMomCoFile
  USE RFLO_ModGridMetrics,   ONLY: RFLO_CalcGridMetrics
  USE RFLO_ModRestartInfo,   ONLY: RFLO_ReadRestartInfo
  USE RFLO_ModDegenerateCornEdge, ONLY: RFLO_WriteDegeneratEC, &
                                        RFLO_MarkDegeneratVert

#ifdef STATS
  USE ModStatsRoutines,        ONLY : StatMapping, InitStatistics, &
                                      StatBuildVersionString
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_CalcMetrics, &
                                      TURB_BuildVersionString
#endif
#ifdef RADI
  USE ModInterfacesRadiation,  ONLY : RADI_BuildVersionString
#endif
#ifdef PERI
  USE ModInterfacesPeriodic,   ONLY : PERI_BuildVersionString
#endif
#ifdef SPEC
  USE ModInterfacesSpecies,    ONLY : SPEC_BuildVersionString
#endif
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_RFLO_SetMetrics, &
                                      PLAG_BuildVersionString
#endif
#ifdef PEUL
  USE ModInterfacesEulerian,   ONLY : PEUL_BuildVersionString
#endif
#ifdef INRT
  USE ModInterfacesInteract,   ONLY : INRT_BuildVersionString
#endif
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE
#ifdef GENX
  INCLUDE "comf90.h"
#endif

! ... parameters
#ifdef GENX
  CHARACTER(*), INTENT(in) :: inSurf, inVol
  DOUBLE PRECISION, INTENT(in) :: initialTime
  INTEGER, INTENT(in)  :: communicator, genxHandle, obtain_attribute

  TYPE(t_globalGenx), POINTER :: globalGenx
#else
  CHARACTER(*) :: casename

  INTEGER :: verbLevel
#endif
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN) :: msg, versionString, headerString, fname

  INTEGER :: headerWidth
  INTEGER :: solver, error, margin, versionWidth, iReg, errorFlag

  LOGICAL :: fileExists,dummyLogical

  INTEGER :: dummy

!******************************************************************************
! initialize some global variables --------------------------------------------

#ifdef GENX
  global => globalGenx%global

  global%timeStamp   = initialTime
  global%currentTime = initialTime
  global%currentIter = -1           ! no iteration
#else
  global%winName     = ''
  global%casename    = casename
  global%verbLevel   = verbLevel
  global%inDir       = './'
  global%outDir      = './'
  global%currentTime = -1._RFREAL   ! no physical time set
  global%currentIter = -1           ! no iteration
#endif

  global%nFunTree = 0
  CALL RegisterFunction( global,'RFLO_InitFlowSolver',&
  'RFLO_InitFlowSolver.F90' )

  global%nProcAlloc = 1
  global%myProcid   = MASTERPROC    ! default process number (if not MPI)
  global%mpierr     = ERR_NONE
  global%error      = ERR_NONE
 
  global%pi            = 4._RFREAL*ATAN(1._RFREAL)
  global%rad           = global%pi/180._RFREAL
  global%calcCellCtr   = .FALSE.
  global%calcFaceCtr   = .FALSE.

! global grid motion

  global%moveGridNbour = 6

#ifdef GENX
! read GenX control file ------------------------------------------------------

  fname = TRIM(global%winName)//'/RocfloControl.txt'

  OPEN(IF_CONTROL,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  READ(IF_CONTROL,'(A)',iostat=errorFlag) global%casename
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

  READ(IF_CONTROL,*,iostat=errorFlag) global%verbLevel
  global%error = errorFlag
  IF (global%error /= 0) global%verbLevel = 1

  READ(IF_CONTROL,'(A)',iostat=errorFlag) global%inDir
  global%error = errorFlag
  IF (global%error /= 0) THEN
    global%inDir = TRIM(global%winName)//'/'
  ELSE
    IF (global%inDir(LEN_TRIM(global%inDir):LEN_TRIM(global%inDir)) /= '/') &
      global%inDir = TRIM(global%inDir)//'/'
  ENDIF

  READ(IF_CONTROL,'(A)',iostat=errorFlag) global%outDir
  global%error = errorFlag
  IF (global%error /= 0) THEN
    global%outDir = TRIM(global%winName)//'/'
  ELSE
    IF (global%outDir(LEN_TRIM(global%outDir):LEN_TRIM(global%outDir)) /= '/') &
      global%outDir = TRIM(global%outDir)//'/'
  ENDIF
  global%error = ERR_NONE
  CLOSE(IF_CONTROL)
#endif

! start up solver and MPI -----------------------------------------------------

#ifdef MPI
#ifndef GENX
  global%mpiComm = MPI_COMM_WORLD
  CALL MPI_Init( global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  global%mpiTagMax = 32768
#else
  global%mpiComm = communicator
#endif
  dummy = 0
  dummyLogical = .true.
  CALL MPI_Attr_get(MPI_COMM_WORLD,MPI_TAG_UB,global%mpiTagMax,dummyLogical,dummy)
!  IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE) THEN
!    WRITE(STDOUT,*) SOLVER_NAME//'  Maximum MPI Tag: ',global%mpiTagMax
!  ENDIF
  CALL MPI_Comm_Size( global%mpiComm,global%nProcAlloc,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

  CALL MPI_Comm_Rank( global%mpiComm,global%myProcid,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE) THEN
    CALL BuildVersionString( versionString )
    headerWidth  = 53
    headerString = ' '
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    headerString(1:1) = '*'
    headerString(headerWidth:headerWidth) = '*'

    WRITE(STDOUT,'(/,A)') SOLVER_NAME//'  *****************************************************'
    WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
    WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                     RocfloMP                      *'
    WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
!    WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
    WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (C) 2015 Illinois Rocstar LLC.       *'
    WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'
  ENDIF

! write out messages from conditional compilation -----------------------------

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_LOW) THEN
    headerWidth  = 32
    headerString = ' '

#ifdef CHECK_GRAD
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' -----------------------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//'                    WARNING                     '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' Compiled to check for gradients of u,v,w and T '
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' -----------------------------------------------'
#endif

#ifdef STATS
    CALL StatBuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' Compiled with Statistics module '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif

#ifdef PLAG
    CALL PLAG_BuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' Compiled with Lagrangian module '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif

#ifdef PEUL
    CALL PEUL_BuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//'  Compiled with Eulerian module  '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif

#ifdef RADI
    CALL RADI_BuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//'  Compiled with Radiation module '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif

#ifdef SPEC
    CALL SPEC_BuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//'   Compiled with Species module  '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif

#ifdef TURB
    CALL TURB_BuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' Compiled with Turbulence module '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif

#ifdef INRT
    CALL INRT_BuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' Compiled with Interaction module'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif

#ifdef PERI
    CALL PERI_BuildVersionString( versionString )
    versionWidth = LEN_TRIM(versionString)
    margin       = (headerWidth-versionWidth)/2
    headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' --------------------------------'
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//'  Compiled with Periodic module  '
    WRITE(STDOUT,'(A)'  ) SOLVER_NAME//' '//TRIM(headerString)
    WRITE(STDOUT,'(A,/)') SOLVER_NAME//' --------------------------------'
#endif
  ENDIF

! check for stop file - delete it if there ------------------------------------

  INQUIRE(file="STOP",exist=fileExists)
  IF (fileExists) THEN
#ifdef GENX
    errorFlag = COM_call_system( "rm -f STOP")
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_SYSTEM_COMMAND,__LINE__,'rm -f STOP' )
#else
    CALL System( 'rm -f STOP' )
#endif
  ENDIF

! read region topology --------------------------------------------------------

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_LOW) THEN
    WRITE(STDOUT,'(A,A,A,/)') SOLVER_NAME//' Case <', &
      TRIM(global%casename),'> running'
  ENDIF
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Reading region topology ...'
  ENDIF

#ifdef GENX
  CALL RFLO_ReadRegionTopology( global,regions )
  globalGenx%regions => regions
#else
  CALL RFLO_ReadRegionTopology( global,regions )
#endif

#ifdef MPI
  global%nProcAlloc = MIN(global%nProcAlloc,global%nRegions)
#endif

! find source patches on adjacent regions

  CALL RFLO_FindSourcePatches( regions )

! get user parameters ---------------------------------------------------------

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_LOW) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Reading user input ...'

  CALL RFLO_GetUserInput( regions )

! initialize random number generator ------------------------------------------

  CALL RFLO_RandomInit( regions )

! read restart info -----------------------------------------------------------

#ifndef GENX
  CALL RFLO_ReadRestartInfo( global )
#endif

! find out solver type (Euler/Navier-Stokes) for GenX -------------------------

#ifdef GENX
  solver = 0
  DO iReg=1,global%nRegions
    IF (regions(iReg)%mixtInput%flowModel == FLOW_NAVST) solver = 1
  ENDDO
#endif

! check if BCs defined for all boundary faces ---------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Checking BCs of regions ...'

  CALL RFLO_CheckRegionFaces( regions )

! check if enough cells to provide data to adjacent region

  CALL RFLO_CheckMinimumCells( regions )

! find source regions for edge & corner cells ---------------------------------
#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Searching source regions for edge & corner cells ...'

  CALL RFLO_FindSourceRegions( regions )

! allocate memory -------------------------------------------------------------

  
#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Allocating memory ...'

  CALL RFLO_DoMemoryAllocation( regions )

! write degenerate edges & corners if exist, and mark them --------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%degenrtEc .AND. global%myProcid==MASTERPROC) THEN
    IF (global%verbLevel>=VERBOSE_HIGH) &
      WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Write degenerated edges and corners into file ...'

    CALL RFLO_WriteDegeneratEC( regions )
  ENDIF
  CALL RFLO_MarkDegeneratVert( regions )

#ifdef STATS
! statistics mapping (must be done before initGenxInterface) ------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
    IF (global%flowType == FLOW_UNSTEADY .AND. global%doStat==ACTIVE) THEN
       IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
            WRITE(STDOUT,'(A)') SOLVER_NAME// &
            ' Doing stat mapping ...'
       CALL StatMapping( global )
    ENDIF
#endif

! get initial grid data -------------------------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Reading grid file, exchanging geometry ...'

  CALL RFLO_GetGeometry( regions,1 )

#ifdef GENX
! initial grid procedures -----------------------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Initializing Grid Procedures ...'
  CALL RFLO_InitGridProcedures( regions )

! init. interface; restore grid and flow solution if restart ------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Preparing GenX interface ...'

  CALL RFLO_InitGenxInterface( regions,genxHandle,solver,inSurf,inVol, &
                               obtain_attribute )

! exchange geometry between regions -------------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Getting geometry ...'
  CALL RFLO_GetGeometry( regions,0 )
#endif

! calculate & check grid metrics ----------------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Calculating & checking grid metrics ...'

  CALL RFLO_CalcGridMetrics( regions )

! calculate metrics of MP modules

#ifdef TURB
#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Calculating & checking turbulence metrics ...'
  IF (global%turbActive) CALL TURB_CalcMetrics( regions, 1 )
#endif

#ifdef PLAG
#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Setting metrics for lagrangian particles ...'
  CALL PLAG_RFLO_SetMetrics( regions )
#endif

! read flow solution ----------------------------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Reading flow solution ...'

  CALL RFLO_GetFlowSolution( regions )

! find patches for thrust integration -----------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%thrustType /= THRUST_NONE) THEN
    IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
      WRITE(STDOUT,'(A)') SOLVER_NAME//' Searching for thrust patches ...'
#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE) THEN             ! on my processor
        CALL RFLO_FindThrustPatches( regions(iReg),iReg )
      ENDIF
    ENDDO
  ENDIF

#ifdef STATS
! statistics initialization ---------------------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Initializing statistics ...'
  CALL InitStatistics( regions )
#endif

! open files for convergence history and time evolution data ------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Opening native output files ...'
  CALL RFLO_OpenConverFile( global )

  CALL RFLO_OpenProbeFile( regions )

  IF (global%thrustType /= THRUST_NONE) &
    CALL RFLO_OpenThrustFile( global )

  IF (global%aeroCoeffs == ACTIVE .AND. global%myProcid==MASTERPROC) &
    CALL RFLO_OpenForceMomCoFile( global )
#ifndef GENX
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE) &
    WRITE(STDOUT,'(//,A,/,A,/,A)') SOLVER_NAME//'  Time stepping', &
                                   SOLVER_NAME//'  =============', &
                                   SOLVER_NAME

! write header for convergence history

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME
    IF (global%flowType == FLOW_STEADY) &
      WRITE(STDOUT,1010) SOLVER_NAME,SOLVER_NAME
    IF (global%flowType == FLOW_UNSTEADY) &
      WRITE(STDOUT,1015) SOLVER_NAME,SOLVER_NAME
  ENDIF
#endif
#ifdef GENX
! send initial data to GenX ---------------------------------------------------

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' Sending boundary values to Rocstar ...'
  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE) THEN             ! on my processor
      CALL RFLO_SendBoundaryValues( regions(iReg),global%timeStamp<=0)
    ENDIF     ! region on this processor and active
  ENDDO       ! iReg
#endif

#ifdef MPI
      CALL MPI_Barrier( global%mpiComm,global%mpierr )
#endif
  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//' All processors initialized ...'
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! formats

1010 FORMAT(A,'  iter',4X,'res-norm',5X,'force-x',6X,'force-y',6X,'force-z', &
            6X,'mass-in',6X,'mass-out',/,A,1X,84('-'))
1015 FORMAT(A,'  time',10X,'delta-t',6X,'force-x',6X,'force-y',6X,'force-z', &
            6X,'mass-in',6X,'mass-out'/,A,1X,90('-'))

END SUBROUTINE RFLO_InitFlowSolver

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_InitFlowSolver.F90,v $
! Revision 1.18  2009/08/12 04:15:58  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.17  2009/03/02 00:19:35  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.16  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2006/03/25 01:14:30  wasistho
! open forceMomCoeff file only by masterproc
!
! Revision 1.13  2006/03/24 23:30:54  wasistho
! added RFLO_OpenForceMomCoFile
!
! Revision 1.12  2006/03/04 04:29:08  wasistho
! moved calcGridMetrics to a rocflo module
!
! Revision 1.11  2006/02/01 20:01:59  wasistho
! added ReadRestartInfo
!
! Revision 1.10  2005/11/11 07:19:52  wasistho
! added RFLO_MarkDegeneratVert
!
! Revision 1.9  2005/10/20 06:50:25  wasistho
! initialize calcFaceCtr
!
! Revision 1.8  2005/06/28 08:51:22  rfiedler
! Remove local currentTime; use timeStamp in place of currentTime to open probes.
!
! Revision 1.7  2005/06/23 01:38:25  wasistho
! initialize global%moveGridNbour to 6
!
! Revision 1.6  2005/05/28 21:24:43  wasistho
! move second getGeometry within ifdef GENX
!
! Revision 1.5  2005/05/28 08:09:04  wasistho
! activate RFLO_InitGridProcedures
!
! Revision 1.4  2005/05/27 08:06:16  wasistho
! allow genx read initial grid
!
! Revision 1.3  2004/12/28 20:27:22  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.2  2004/12/01 00:20:58  wasistho
! added phys. modules version strings
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.44  2004/11/29 17:15:46  wasistho
! use ModInterfacesStatistics
!
! Revision 1.43  2004/08/21 00:32:59  wasistho
! added write file degenerated edge/corners
!
! Revision 1.42  2004/07/23 23:26:47  wasistho
! separate system command btw genx and standalone
!
! Revision 1.41  2004/07/06 23:24:00  jiao
! Jiao: Changed to call RFLO_SendBoundaryValues only at time 0.
!
! Revision 1.40  2004/06/29 23:57:17  wasistho
! migrated to Roccom-3
!
! Revision 1.39  2004/06/07 23:09:35  wasistho
! moved statistics mapping from initStatistics to before initGenxInterfaces
!
! Revision 1.38  2004/02/04 22:29:43  wasistho
! add integer argument isInit to TURB_calcMetrics
!
! Revision 1.37  2003/12/07 04:52:21  jiao
! Changed the call to RFLO_ReadRegionTopology to work with PGI compilers.
!
! Revision 1.36  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.35  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.31  2003/11/12 21:21:06  fnajjar
! Added Corner-Edge cells routine to communicate metrics for PLAG
!
! Revision 1.30  2003/10/03 20:18:43  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.29  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.28  2003/06/02 17:12:01  jblazek
! Added computation of thrust.
!
! Revision 1.27  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.26  2003/04/15 20:57:53  jiao
! Jiri: Added close statement for the genx-control file.
!
! Revision 1.25  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.24  2003/02/17 19:31:10  jferry
! Implemented portable random number generator ModRandom
!
! Revision 1.23  2003/02/14 22:32:37  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.22  2003/02/03 19:20:47  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.21  2002/10/16 18:30:38  jblazek
! Within GenX, BC data at t=0 are updated in FlowSolver before calling
! the time-stepping routine.
!
! Revision 1.20  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.19  2002/10/02 22:21:59  jiao
! Debugged GenX restart.
!
! Revision 1.18  2002/09/27 22:42:41  jferry
! removed PFEU reference
!
! Revision 1.17  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.16  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.15  2002/07/24 17:30:06  wasistho
! Added Rocturb header
!
! Revision 1.14  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.13  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.12  2002/06/18 00:34:32  wasistho
! Added prefix SOLVER NAME to satistics STDOutput
!
! Revision 1.11  2002/06/17 16:09:17  wasistho
! Added STATS compilation message
!
! Revision 1.10  2002/06/14 21:38:45  wasistho
! Added time avg statistics
!
! Revision 1.9  2002/06/13 23:06:20  jblazek
! Added version string.
!
! Revision 1.8  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.7  2002/04/12 17:36:23  jblazek
! Added timer.
!
! Revision 1.6  2002/04/03 02:28:52  jblazek
! Added x,y,z location to probe file header.
!
! Revision 1.5  2002/04/01 19:36:08  jblazek
! Added routine to clear send requests.
!
! Revision 1.4  2002/03/30 00:50:49  jblazek
! Cleaned up with flint.
!
! Revision 1.3  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.2  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.1  2002/02/25 22:36:53  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************







