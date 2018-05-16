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
! Purpose: write out ROCPART solution for visualization.
!
! Description: currently supported formats are:
!              - TECPLOT ASCII
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
! $Id: PLAG_PostProcessing.F90,v 1.5 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

PROGRAM ROCPART_Post

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE PLAG_ModInterfacesPost, ONLY : BuildVersionString, ReadInputFile,  &
                                     RFLO_DerivedInputValues,            &
                                     RFLO_ReadRegionTopology,            &
                                     RFLO_InitInputValues,               &
                                     PLAG_AllocateMemoryPost,            &
                                     PLAG_BinSortNozzleInlet,            &
                                     PLAG_BinSortSpatialDist,            &
				     PLAG_DeallocateMemoryPost,          &
                                     PLAG_CalcDerivedVariables,          &
                                     PLAG_IntrpMixtProperties,           &
				     PLAG_ProcessEulerField,             &
                                     PLAG_ReadSolutionFilePost,          &
                                     PLAG_UserInput,                     &
                                     PLAG_WriteTecplotAscii

#ifdef STATS  
  USE PLAG_ModStats, ONLY: PLAG_CreateStat, PLAG_DestroyStat
  USE PLAG_ModInterfacesPost, ONLY : PLAG_ReadStatPost, &
                                     PLAG_WriteStatTecAscii
#endif

  USE ModMPI
  USE ModParameters
  USE PLAG_ModParameters

  IMPLICIT NONE

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString, stamp, outFmt, msg, &
                       verbosity, versionString, headerString

  INTEGER :: iLev, outputFormat
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER :: nPclsSum
  INTEGER :: iRegBin
  INTEGER, PARAMETER :: headerWidth = 53
  INTEGER, PARAMETER :: NSTATS_TEC_PLAG = 10

#ifdef STATS
  LOGICAL :: statsActive
#endif

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)
  TYPE(t_plag)  , POINTER :: pPlag

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_PostProcessing.F90,v $ $Revision: 1.5 $'

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'ROCPART_Post', 'PLAG_PostProcessing.F90' )

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

  global%startLevel = 1
  
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
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        ROCPART: Solution Postprocessing           *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        ================================           *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! read argument list ----------------------------------------------------------

  CALL GETARG(1,global%casename)
  CALL GETARG(2,stamp)
  CALL GETARG(3,outFmt)
  CALL GETARG(4,verbosity)

  IF (LEN_TRIM(global%casename)==0 .OR. &
      LEN_TRIM(stamp)==0           .OR. &
      LEN_TRIM(verbosity)==0       .OR. &
      LEN_TRIM(outFmt)==0) THEN
    WRITE(STDOUT,'(/,A,/,A,/,9(A,/))') &
      SOLVER_NAME//' Usage: plagpost <casename> <time> <format> <verbosity>', &
      SOLVER_NAME, &
      SOLVER_NAME//'        time      = time ', &
      SOLVER_NAME//' ', &
      SOLVER_NAME//'        format    = 3 - Tecplot ASCII',&
      SOLVER_NAME//' ', &
      SOLVER_NAME//'        verbosity = 0-2 '
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF

  READ(outFmt ,*) outputFormat
  
  IF (outputFormat <= 1) THEN
    outputFormat = PLOT_FMT_GENERIC
  ELSE IF (outputFormat == 2) THEN
    outputFormat = PLOT_FMT_TECPLOT
  ELSE
    outputFormat = PLOT_FMT_TECASCII
  ENDIF

! check for TECPLOT library

#ifdef NO_TECPLOT
  IF (outputFormat == PLOT_FMT_TECPLOT) THEN
    WRITE(STDOUT,'(/,A,/)') SOLVER_NAME// &
                            ' Sorry, not linked to TECPLOT library.'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF
#endif

  READ(verbosity,*) global%verbLevel

! read region topology --------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading region topology ...'

  CALL RFLO_ReadRegionTopology( global,regions )

  DO iReg=1,global%nRegions
    regions(iReg)%startLevel = global%startLevel
    regions(iReg)%currLevel  = global%startLevel
    IF (regions(iReg)%nGridLevels < regions(iReg)%currLevel) THEN
      WRITE(msg,1000) SOLVER_NAME,iReg,global%startLevel
      CALL ErrorStop( global,ERR_GRID_LEVEL,__LINE__,msg )
    ENDIF
  ENDDO

! get user parameters ---------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading user input ...'

  CALL RFLO_InitInputValues( regions )
  CALL ReadInputFile( regions )
  CALL PLAG_UserInput( regions )

  IF (global%flowType == FLOW_STEADY) THEN
    WRITE(STDOUT,'(/,A,/)') SOLVER_NAME// &
                            ' Sorry, unable to run Rocpart Post-processing ', &
                            ' tool for steady state case.'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP   
  ELSE
    READ(stamp,*) global%timeStamp
    global%currentTime = global%timeStamp
  ENDIF

  CALL RFLO_DerivedInputValues( regions )

#ifndef PLAG
  WRITE(STDOUT,'(/,A,/)') SOLVER_NAME// &
                          ' Sorry, unable to run Rocpart Post-processing ', &
                          ' tool since code was not compiled with PLAG.'
#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif
  STOP
#endif

#ifdef STATS
  statsActive = ( (global%plagUsed .EQV. .TRUE.)     .AND. &
                  (global%flowType == FLOW_UNSTEADY) .AND. &
                  (global%doStat == ACTIVE)                )
#endif

! allocate memory -------------------------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Reading ROCPART solution and writing plot file ...'

  DO iReg=1,global%nRegions
    CALL PLAG_AllocateMemoryPost( regions(iReg), iReg )
  ENDDO   ! iReg
  
! read solution ---------------------------------------------------------------

  CALL PLAG_ReadSolutionFilePost( regions )

! initialize cumulative sum ---------------------------------------------------

  nPclsSum = 0

! write data ------------------------------------------------------------------

  DO iReg=1,global%nRegions

    WRITE(STDOUT,'(A,I5.5)') SOLVER_NAME//'   - region ',iReg
    
    iLev =  regions(iReg)%currLevel

! - compute cumulative number of particles in all regions

    nPclsSum = nPclsSum + regions(iReg)%levels(iLev)%plag%nPcls 

    IF ( iReg == global%nRegions ) &
      WRITE(STDOUT,'(A,I8.8)') 'Total Number of Particles = ',nPclsSum  
        
! - calculate derived variables

    CALL PLAG_calcDerivedVariables( regions(iReg) )

! - write data 

    SELECT CASE ( outputFormat )

! -- write data to TECPLOT file (ASCII)

      CASE ( PLOT_FMT_TECASCII ) 
        CALL PLAG_WriteTecplotAscii( iReg,iLev,regions(iReg) )

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! outputFormat
  ENDDO   ! iReg

! write bined particle data at nozzle inlet (ASCII format)

  DO iReg=1,global%nRegions
    IF ( iReg == 41 .OR. iReg == 52 ) THEN
      iRegBin = iReg
      CALL PLAG_BinSortNozzleInlet( iReg,iLev,regions(iReg), iRegBin )
      CALL PLAG_BinSortSpatialDist( iReg,iLev,regions(iReg), iRegBin )
    END IF ! iReg 
  ENDDO   ! iReg

! process Eulerian grid-based data for PLAG -----------------------------------

  DO iReg=1,global%nRegions
    CALL PLAG_ProcessEulerField( regions,iReg,nPclsSum )    
  ENDDO   ! iReg
  
! deallocate memory -----------------------------------------------------------

  DO iReg=1,global%nRegions
    CALL PLAG_DeallocateMemoryPost( regions(iReg), iReg )
  ENDDO   ! iReg

! ******************************************************************************
! Save currentTime as global value is clobbered by RFLO_ReadGridRegion
! ******************************************************************************
  
  global%currentTime = global%timeStamp

! process Eulerian grid-based statistics data for PLAG ------------------------

#ifdef STATS
  IF ( statsActive .EQV. .TRUE. ) THEN 

! - trap error if expected number of statistics is less then set value  -------

    IF ( global%plagNStat > 0 ) THEN
      IF ( global%plagNStat < NSTATS_TEC_PLAG ) &
        CALL ErrorStop( global,ERR_STATS_TECPLOT,__LINE__, &
                        'plagNStat < NSTATS_TEC_PLAG' )
    ENDIF ! plagNStat > 0

! - write tecplot statistics file ---------------------------------------------
    
    DO iReg=1,global%nRegions

! -- allocate memory ----------------------------------------------------------

      IF ( iReg== 1 .AND. global%myProcid == MASTERPROC .AND. &                              
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Allocating memory for statistics...'
      END IF ! global%verbLevel
      
      DO iLev=1,regions(iReg)%nGridLevels
        pPlag => regions(iReg)%levels(iLev)%plag
        CALL PLAG_CreateStat( regions(iReg), pPlag )
      END DO ! iLev

      IF ( iReg == global%nRegions .AND. global%myProcid == MASTERPROC .AND. &                              
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Allocating memory for statistics done...'
      END IF ! global%verbLevel
      
! -- read statistics file ------------------------------------------------------
      
      IF ( iReg== 1 .AND. global%myProcid == MASTERPROC .AND. &                              
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading statistics solution file for PLAG...'
      END IF ! global%verbLevel
      
      CALL PLAG_ReadStatPost( regions,iReg )    

      IF ( iReg == global%nRegions .AND. global%myProcid == MASTERPROC .AND. &                              
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading statistics solution file for PLAG done...'
      END IF ! global%verbLevel
      
! -- write data to file --------------------------------------------------------

      WRITE(STDOUT,'(A,I5.5)') SOLVER_NAME//'   Statistics - region ',iReg

      SELECT CASE ( outputFormat )

! --- TECPLOT file (ASCII) -----------------------------------------------------

        CASE ( PLOT_FMT_TECASCII ) 
          CALL PLAG_WriteStatTecAscii( regions,iReg )

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! outputFormat
 
! -- deallocate memory ---------------------------------------------------------

      IF ( iReg== 1 .AND. global%myProcid == MASTERPROC .AND. &                              
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Deallocating memory for statistics...'
      END IF ! global%verbLevel
      
      DO iLev=1,regions(iReg)%nGridLevels
        pPlag => regions(iReg)%levels(iLev)%plag
        CALL PLAG_DestroyStat( regions(iReg), pPlag )
      END DO ! iLev

      IF ( iReg == global%nRegions .AND. global%myProcid == MASTERPROC .AND. &                              
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Deallocating memory for statistics done...'
      END IF ! global%verbLevel
        
    ENDDO   ! iReg


  ENDIF ! statsActive
#endif
   
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif

1000 FORMAT(A,' Region ',I5,', grid level= ',I2,'.')

END PROGRAM ROCPART_Post

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PostProcessing.F90,v $
! Revision 1.5  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/02/16 23:44:16  fnajjar
! Moved statistics-specific routines inside ifdef construct
!
! Revision 1.2  2005/02/16 14:49:23  fnajjar
! Added infrastructure to write Tecplot-based statistics file
!
! Revision 1.1  2004/12/01 22:00:46  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2004/11/17 22:13:17  fnajjar
! Cosmetic changes, removed support for enhanced Tecplot files and added call for eulerian data processing
!
! Revision 1.7  2004/11/13 21:58:48  fnajjar
! Added spatial binary sort and streamlined calls to binary sorts
!
! Revision 1.6  2004/06/29 14:06:19  fnajjar
! Removed call to interpolation of mixture properties since mixture data not allocated
!
! Revision 1.5  2004/05/24 14:24:18  fnajjar
! Included interpolation for mixture and binning routine
!
! Revision 1.4  2004/03/20 23:47:52  fnajjar
! Updated executable name in error trapping from rplagpost to plagpost
!
! Revision 1.3  2003/07/30 23:21:53  fnajjar
! Included cumulative sum of particles in IO
!
! Revision 1.2  2003/05/28 13:57:36  fnajjar
! Removed IndIJK as being an obsolete option
!
! Revision 1.1.1.1  2003/05/06 16:14:38  fnajjar
! Import of postprocessing tool for Rocpart
!
!******************************************************************************







