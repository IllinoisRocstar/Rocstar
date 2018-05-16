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
! Purpose: write out ROCFLO-MP`s grid and/or solution for visualization.
!
! Description: currently supported formats are:
!              - generic binary
!              - TECPLOT
!
! Input: case name from the list of arguments
!
! Output: to file.
!
! Notes: the output is collected in one file, but the regions are processed
!        separately to save memory.
!
!******************************************************************************
!
! $Id: POST_Main.F90,v 1.5 2008/12/06 08:44:49 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

PROGRAM ROCFLO_Post

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE POST_ModInterfaces, ONLY : RFLO_ReadRegionTopology, &
        RFLO_InitInputValues, ReadInputFile, RFLO_DerivedInputValues, &
        RFLO_GetDimensPhysNodes, RFLO_GetDimensDummy, &
        RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, RFLO_GetCellOffset, &
        RFLO_ReadGridRegion, RFLO_ReadSolutionRegion, RFLO_ReadStatRegion, &
        MixtureProperties, WriteTecplotAscii, WriteTecplotBinary, WriteGeneric,&
        RFLO_GenerateCoarseGrids, RFLO_CopyGeometryDummy, BuildVersionString, &
        PrintPostInput
  USE ModMPI
  USE ModParameters
#ifdef TURB
  USE ModTurbulence, ONLY : t_turb
  USE POST_ModInterfaces, ONLY : TURB_ReadInputFile, &
                        TURB_DerivedInputValues, TURB_RansSAGetEddyVis, &
                        TURB_RFLO_ReadSolutionRegion
  USE TURB_ModParameters
#endif
  IMPLICIT NONE

#include "Indexing.h"

! ... loop variables
  INTEGER :: iReg, iLev

! ... local variables
  CHARACTER(CHRLEN) :: casename, verbosity, msg, versionString, headerString

  INTEGER :: ipc, jpc, kpc, ibc, iec, ibn, ien, mixtNStatTec, turbNStatTec
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 53

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)
  TYPE(t_grid)  , POINTER :: grid
  TYPE(t_mixt)  , POINTER :: mixt
#ifdef TURB
  TYPE(t_turb)  , POINTER :: turb
#endif
#ifdef STATS
  LOGICAL :: statsActive
#endif

!******************************************************************************

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'ROCFLO_Post', &
                        'POST_Main.F90' )

! initialize global parameters ------------------------------------------------

  global%verbLevel = VERBOSE_NONE

  global%flowType    = FLOW_STEADY  ! stationary flow
  global%currentTime = -1._RFREAL   ! no physical time set
  global%currentIter = -1           ! no iteration

  global%inDir  = './'              ! directory path
  global%outDir = './'

  global%nProcAlloc = 1
  global%myProcid   = MASTERPROC    ! default process number (not an MPI code)
  global%mpierr     = ERR_NONE
  global%error      = ERR_NONE

  global%pi  = 4._RFREAL*ATAN(1._RFREAL)
  global%rad = global%pi/180._RFREAL

! print header ----------------------------------------------------------------

#ifdef MPI
  CALL MPI_Init( global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

  CALL BuildVersionString( versionString )

  headerString = ' '
  versionWidth = LEN_TRIM(versionString)
  margin       = (headerWidth-versionWidth)/2
  headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
  headerString(1:1) = '*'
  headerString(headerWidth:headerWidth) = '*'

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//'  *****************************************************'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        ROCFLO-MP: Solution Postprocessing         *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        ==================================         *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! print required input and read argument list ---------------------------------

  WRITE(STDOUT,'(/,A,/,A,/,18(A,/))') &
    SOLVER_NAME//' Required Rocflo-Post input:', &
    SOLVER_NAME, &
    SOLVER_NAME//'   <casename>         =>  command line input', &
    SOLVER_NAME//'   <verbosity>        =>  command line input', &
    SOLVER_NAME//'   grid level         => .inp file: # MULTIGRID: START', &
    SOLVER_NAME//'   unsteadiness       => .inp file: # TIMESTEP: FLOWTYPE', &
    SOLVER_NAME//'   turbulence         => .inp file: # TURBULENCE: TURBMODEL', &
    SOLVER_NAME//'                      => .inp file: # TURBULENCE: OUTPUTNUMBER', &
    SOLVER_NAME//'   gas properties     => .inp file: # REFERENCE', &
    SOLVER_NAME//'                      => .inp file: # VISCMODEL', &
    SOLVER_NAME//'                      => .inp file: # MATERIAL (if MP is active)', &
    SOLVER_NAME//'   plot type          => .inp file: # POST: PLTTYPE', &
    SOLVER_NAME//'   time (unsteady)    => .inp file: # POST: TIME', &
    SOLVER_NAME//'   iteration (steady) => .inp file: # POST: ITER', &
    SOLVER_NAME//'   output plot format => .inp file: # POST: OUTFORMAT', &
    SOLVER_NAME//'   include statistics => .inp file: # POST: STATSFLAG', &
    SOLVER_NAME//'   include turbulence => .inp file: # POST: TURBFLAG', &
    SOLVER_NAME//'   include particles  => .inp file: # POST: PLAGFLAG', &
    SOLVER_NAME//'   include radiation  => .inp file: # POST: RADIFLAG', &
    SOLVER_NAME//'   include species    => .inp file: # POST: SPECFLAG'

  CALL GETARG(1,casename)
  CALL GETARG(2,verbosity)

  IF (LEN_TRIM(casename)==0 .OR. &
      LEN_TRIM(verbosity)==0) THEN
    WRITE(STDOUT,'(/,A,/)') &
      SOLVER_NAME//' Usage: rflopost <casename> <verbosity>'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF

  READ(casename ,*) global%casename
  READ(verbosity,*) global%verbLevel

! read region topology --------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading region topology ...'

  CALL RFLO_ReadRegionTopology( global,regions )

! get user parameters ---------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading user input ...'

  CALL RFLO_InitInputValues( regions )
  CALL ReadInputFile( regions )

! copy times and iteration to corresponding post values

  global%currentTime = global%postTime
  global%timeStamp   = global%postTime
  global%currentIter = global%postIter
 
! check grid level and obtain region grid size

  DO iReg=1,global%nRegions
    regions(iReg)%startLevel = global%startLevel
    regions(iReg)%currLevel  = global%startLevel
    IF (regions(iReg)%nGridLevels < regions(iReg)%currLevel) THEN
      WRITE(msg,1000) SOLVER_NAME,iReg,global%startLevel
      CALL ErrorStop( global,ERR_GRID_LEVEL,__LINE__,msg )
    ENDIF
    DO iLev=2,regions(iReg)%nGridLevels
      ipc = regions(iReg)%levels(iLev-1)%grid%ipc
      jpc = regions(iReg)%levels(iLev-1)%grid%jpc
      kpc = regions(iReg)%levels(iLev-1)%grid%kpc
      regions(iReg)%levels(iLev)%grid%ipc = ipc/2
      regions(iReg)%levels(iLev)%grid%jpc = jpc/2
      regions(iReg)%levels(iLev)%grid%kpc = kpc/2
    ENDDO
  ENDDO

! continue collecting remaining required user input

  CALL RFLO_DerivedInputValues( regions )
#ifdef TURB
  CALL TURB_ReadInputFile( regions )
  CALL TURB_DerivedInputValues( regions )
#endif

#ifdef STATS
  statsActive = (global%postStatsFlag .AND. &
                (global%flowType == FLOW_UNSTEADY) .AND. &
                (global%doStat == ACTIVE))
#endif

! check for TECPLOT library ---------------------------------------------------

#ifdef NO_TECPLOT
  IF (global%postOutFmt == PLOT_FMT_TECPLOT) THEN
    WRITE(STDOUT,'(/,A,/)') SOLVER_NAME// &
                            ' Sorry, not linked to TECPLOT library.'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF
#endif

! loop over regions -----------------------------------------------------------

  IF (global%postPlotType == PLOT_GRID_ONLY) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Reading grid and writing plot file ...'
  ELSE
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Reading grid+solution and writing plot file ...'
  ENDIF

  DO iReg=1,global%nRegions

! - allocate memory for the grid (all grid levels)

    DO iLev=1,regions(iReg)%nGridLevels
      CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                     jdnbeg,jdnend,kdnbeg,kdnend )
      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
      ibn  =  IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
      ien  =  IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
      grid => regions(iReg)%levels(iLev)%grid
      ALLOCATE( grid%xyz(3,ibn:ien),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ENDDO

! - allocate memory for the rest (current grid level)

    iLev =  regions(iReg)%currLevel
    grid => regions(iReg)%levels(iLev)%grid
    mixt => regions(iReg)%levels(iLev)%mixt

    CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
    iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

    CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                  jpnbeg,jpnend,kpnbeg,kpnend )
    CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
    ibn = IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff)
    ien = IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff)

    IF (regions(iReg)%mixtInput%moveGrid) THEN
      ALLOCATE( grid%siVel(ibn:ien),stat=errorFlag )
      ALLOCATE( grid%sjVel(ibn:ien),stat=errorFlag )
      ALLOCATE( grid%skVel(ibn:ien),stat=errorFlag )
    ENDIF
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    IF (global%postPlotType == PLOT_GRID_FLOW) THEN
      ALLOCATE( mixt%cv(5,ibc:iec) ,stat=errorFlag )
      ALLOCATE( mixt%dv(mixt%nDv,ibc:iec),stat=errorFlag )
      IF (regions(iReg)%mixtInput%computeTv) THEN
        ALLOCATE( mixt%tv(mixt%nTv,ibc:iec),stat=errorFlag )
      ENDIF
      IF (regions(iReg)%mixtInput%gasModel == GAS_MODEL_TCPERF) THEN
        ALLOCATE( mixt%gv(mixt%nGv,0:1),stat=errorFlag )
      ELSE
        ALLOCATE( mixt%gv(mixt%nGv,ibc:iec),stat=errorFlag )
      ENDIF
#ifdef TURB
      IF (global%postTurbFlag .AND. &
          regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
        turb => regions(iReg)%levels(iLev)%turb
        IF (regions(iReg)%turbInput%nOutField > 1) THEN
          ALLOCATE( turb%vort(XCOORD:ZCOORD,ibc:iec),stat=errorFlag )
        ENDIF
        IF (regions(iReg)%turbInput%modelClass == MODEL_RANS) THEN
          ALLOCATE( turb%cv(regions(iReg)%turbInput%nCv,ibc:iec),stat=errorFlag )
          ALLOCATE( turb%lens(ibc:iec),stat=errorFlag )
        ENDIF
      ENDIF
#endif
#ifdef STATS
      IF (statsActive) THEN
        IF (global%mixtNStat > 0) THEN
          IF (global%mixtNStat < NSTATS_TEC_MIXT) &
            CALL ErrorStop( global,ERR_STATS_TECPLOT,__LINE__, &
                           'mixtNStat < NSTATS_TEC_MIXT' )
          ALLOCATE( mixt%tav(global%mixtNStat,ibc:iec),stat=errorFlag )
        ENDIF
#ifdef TURB
        IF (global%turbNStat > 0) THEN
          IF (global%turbNStat < NSTATS_TEC_TURB) &
            CALL ErrorStop( global,ERR_STATS_TECPLOT,__LINE__, &
                           'turbNStat < NSTATS_TEC_TURB' )
          IF (global%mixtNStat <=0) &
            CALL ErrorStop( global,ERR_STATS_TECPLOT,__LINE__, &
                           'Tecplot expect mixtNStat>0 if turbNStat>0' )
          ALLOCATE( turb%tav(global%turbNStat,ibc:iec),stat=errorFlag )
        ENDIF
#endif
      ENDIF
#endif
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ENDIF

! - read grid

    CALL RFLO_ReadGridRegion( iReg,regions )
    CALL RFLO_GenerateCoarseGrids( regions(iReg) )
    CALL RFLO_CopyGeometryDummy( regions(iReg) )

! - read solution, calc. mixture properties

    IF (global%postPlotType == PLOT_GRID_FLOW) THEN
      CALL RFLO_ReadSolutionRegion( iReg,regions )
      CALL MixtureProperties( regions(iReg),ibc,iec,.true. )
#ifdef TURB
      IF (global%postTurbFlag .AND. &
          regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
        IF (regions(iReg)%turbInput%modelClass == MODEL_LES) THEN
          CALL TURB_RFLO_ReadSolutionRegion( iReg,regions )
        ELSEIF ((regions(iReg)%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
                (regions(iReg)%mixtInput%turbModel == TURB_MODEL_DESSA)) THEN
          CALL TURB_RFLO_ReadSolutionRegion( iReg,regions )
          CALL TURB_RansSAGetEddyVis( regions(iReg) )
        ELSE
        ENDIF
      ENDIF
#endif
#ifdef STATS
      IF (statsActive) THEN
        IF (global%mixtNStat > 0 .OR. global%turbNStat > 0) &
          CALL RFLO_ReadStatRegion( iReg,regions )
      ENDIF
#endif
    ENDIF

! - print on screen actual input selected

    IF (global%verbLevel >= VERBOSE_LOW) CALL PrintPostInput( regions(iReg) )
    IF (global%verbLevel == VERBOSE_LOW) WRITE(STDOUT,'(A,I5.5)') &
                                         SOLVER_NAME//'   - region ',iReg

! - write data

    IF (global%postOutFmt == PLOT_FMT_GENERIC) THEN

! --- write data to generic file

      CALL WriteGeneric( iReg,regions(iReg) )

    ELSE IF (global%postOutFmt == PLOT_FMT_TECPLOT) THEN

! --- write data to TECPLOT file (binary)

#ifndef NO_TECPLOT
      CALL WriteTecplotBinary( iReg,regions(iReg) )
#endif

    ELSE IF (global%postOutFmt == PLOT_FMT_TECASCII) THEN

! --- write data to TECPLOT file (ASCII)

      CALL WriteTecplotAscii( iReg,regions(iReg) )

    ENDIF

! - deallocate memory

    DO iLev=1,regions(iReg)%nGridLevels
      grid => regions(iReg)%levels(iLev)%grid
      DEALLOCATE( grid%xyz,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
    ENDDO

    iLev = regions(iReg)%currLevel
    IF (regions(iReg)%mixtInput%moveGrid) THEN
      grid => regions(iReg)%levels(iLev)%grid
      DEALLOCATE( grid%siVel,stat=errorFlag )
      DEALLOCATE( grid%sjVel,stat=errorFlag )
      DEALLOCATE( grid%skVel,stat=errorFlag )
    ENDIF
    IF (global%postPlotType == PLOT_GRID_FLOW) THEN
      DEALLOCATE( mixt%cv ,stat=errorFlag )
      DEALLOCATE( mixt%dv ,stat=errorFlag )
      DEALLOCATE( mixt%tv ,stat=errorFlag )
      DEALLOCATE( mixt%gv ,stat=errorFlag )
    ENDIF
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

#ifdef TURB
    IF (global%postTurbFlag .AND. &
        regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN

      IF (ASSOCIATED(turb%vort)) DEALLOCATE( turb%vort,stat=errorFlag ) !general
      IF (ASSOCIATED(turb%cv))   DEALLOCATE( turb%cv  ,stat=errorFlag ) !RANS
      IF (ASSOCIATED(turb%lens)) DEALLOCATE( turb%lens,stat=errorFlag ) !RANS
    ENDIF
#endif

#ifdef STATS
    IF (statsActive) THEN
      IF (ASSOCIATED( mixt%tav )) DEALLOCATE( mixt%tav,stat=errorFlag )
#ifdef TURB
      IF (ASSOCIATED( turb%tav )) DEALLOCATE( turb%tav,stat=errorFlag )
#endif
    ENDIF
#endif

  ENDDO   ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif

1000 FORMAT(A,' Region ',I5,', grid level= ',I2,'.')

END PROGRAM ROCFLO_Post

!******************************************************************************
!
! RCS Revision history:
!
! $Log: POST_Main.F90,v $
! Revision 1.5  2008/12/06 08:44:49  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/10/31 21:09:38  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.2  2004/12/03 03:20:04  wasistho
! rflo_modinterfacespost to post_modinterfaces
!
! Revision 1.1  2004/12/03 02:03:16  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:32:01  wasistho
! lower to upper case
!
! Revision 1.39  2004/11/09 10:49:30  wasistho
! added statistics to rflopost
!
! Revision 1.38  2004/08/27 04:00:41  wasistho
! fixed dimensions of relative vels when grid moves
!
! Revision 1.37  2004/08/26 01:15:43  wasistho
! more complete screen-display of required input parameters
!
! Revision 1.36  2004/08/26 00:44:30  wasistho
! bug fixed interconnection turbFlag and turbModel
!
! Revision 1.35  2004/08/10 17:40:04  wasistho
! do not compile writeTecplotBinary if there is no tecplot library
!
! Revision 1.34  2004/07/28 01:50:17  wasistho
! added print input
!
! Revision 1.33  2004/07/24 03:50:15  wasistho
! use postSection instead of command line input
!
! Revision 1.32  2004/03/11 03:36:07  wasistho
! changed rocturb nomenclature
!
! Revision 1.31  2004/03/05 21:16:17  wasistho
! fixed usage to include turbulence
!
! Revision 1.30  2004/03/03 23:55:40  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.29  2004/02/11 03:25:33  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.28  2004/02/07 01:18:37  wasistho
! added turbulence related results in rocflo post processing
!
! Revision 1.27  2003/09/19 22:38:11  jblazek
! Added turbulence input.
!
! Revision 1.26  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.25  2003/03/20 22:23:47  haselbac
! Renamed ModInterfaces
!
! Revision 1.24  2003/03/20 19:41:26  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.23  2003/03/20 19:34:37  haselbac
! Modified RegFun call to avoid probs with long 'POST_Main.F90' names
!
! Revision 1.22  2002/10/19 00:40:31  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.21  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.20  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.19  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.18  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.17  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.16  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.15  2002/08/07 20:51:26  jblazek
! Enhanced for moving grids.
!
! Revision 1.14  2002/07/20 00:42:05  jblazek
! Added ASCII Tecplot format.
!
! Revision 1.13  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.12  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.11  2002/06/14 17:16:41  jblazek
! Added version string.
!
! Revision 1.10  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.9  2002/04/11 21:10:27  jblazek
! Set correct time when writing grid only for Tecplot.
!
! Revision 1.8  2002/03/18 21:31:54  jblazek
! Made codes compatible with MPI.
!
! Revision 1.7  2002/02/22 20:30:39  jblazek
! Changed generic format. Enhanced Tecplot title.
!
! Revision 1.6  2002/02/22 00:05:44  jblazek
! Changed TECPLOT link option.
!
! Revision 1.4  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.3  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.2  2002/01/12 00:02:49  jblazek
! Added postprocessor.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************








