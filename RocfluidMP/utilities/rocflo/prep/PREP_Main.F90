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
! Purpose: provide intitial guess for ROCFLO-MP.
!
! Description: the initialization is done sepoarately for each region.
!              Flow conditions can vary between the regions.
!
! Input: case name from the list of arguments
!
! Output: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: PREP_Main.F90,v 1.7 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

PROGRAM ROCFLO_Init

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef GENX
  USE PREP_ModInterfaces, ONLY : GenxInitSolution, GenxWriteSolution
#endif
  USE PREP_ModInterfaces, ONLY : RFLO_ReadRegionTopology, GetGrid, &
           InitInputValues, ReadInputFile, ReadBcInputFile, &
           CheckBcValidity, PrintPrepInput, InitializeFlowField, &
           RFLO_WriteSolutionRegion, BuildVersionString
  USE PREP_ModBcDistribution, ONLY : BcDistributionFiles
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ... loop variables
  INTEGER :: iReg, iLev

! ... local variables
  CHARACTER(CHRLEN) :: casename, verbosity, msg, versionString, headerString
#ifdef GENX
  CHARACTER(CHRLEN) :: wins, winv
#endif

  INTEGER :: ipc, jpc, kpc, gridLevel
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 53

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)

!******************************************************************************

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'ROCFLO_Init',&
  'PREP_Main.F90' )

! initialize global parameters ------------------------------------------------

  global%verbLevel = VERBOSE_NONE

  global%flowType    = FLOW_STEADY  ! stationary flow
  global%currentTime = 0._RFREAL    ! no physical time set
  global%currentIter = 0            ! no iteration yet
  global%resInit     = 1._RFREAL

  global%inDir  = './'              ! directory path
  global%outDir = './'

  global%nProcAlloc = 1
  global%myProcid   = MASTERPROC    ! default process number (if not MPI)
  global%mpierr     = ERR_NONE
  global%error      = ERR_NONE

  global%pi  = 4._RFREAL*ATAN(1._RFREAL)
  global%rad = global%pi/180._RFREAL

! print header ----------------------------------------------------------------

!#ifdef MPI
!  CALL MPI_Init( global%mpierr )
!  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
!#endif

  CALL BuildVersionString( versionString )

  headerString = ' '
  versionWidth = LEN_TRIM(versionString)
  margin       = (headerWidth-versionWidth)/2
  headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
  headerString(1:1) = '*'
  headerString(headerWidth:headerWidth) = '*'

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//'  *****************************************************'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        ROCFLO-MP: Solution Initialization         *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        ==================================         *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! print required input and read argument list ---------------------------------

  WRITE(STDOUT,'(/,A,/,A,/,8(A,/))') &
    SOLVER_NAME//' Required Rocflo-Prep input:', &
    SOLVER_NAME, &
    SOLVER_NAME//'   <casename>        =>  command line input', &
    SOLVER_NAME//'   <verbosity>       =>  command line input', &
    SOLVER_NAME//'   grid level        => .inp file: # MULTIGRID: START', &
    SOLVER_NAME//'   input grd format  => .inp file: # FORMATS: GRID', &
    SOLVER_NAME//'   output sol format => .inp file: # FORMATS: SOLUTION', &
    SOLVER_NAME//'   unsteadiness      => .inp file: # TIMESTEP: FLOWTYPE', &
    SOLVER_NAME//'   spec. heat ratio  => .inp file: # REFERENCE: GAMMA', &
    SOLVER_NAME//'   initial solution  => .inp file: # INITFLOW: all'

  CALL GETARG(1,casename)
  CALL GETARG(2,verbosity)

  IF (LEN_TRIM(casename)==0 .OR. &
      LEN_TRIM(verbosity)==0) THEN
    WRITE(STDOUT,'(/,A,/)') &
      SOLVER_NAME//' Usage: rfloprep <casename> <verbosity>'
!#ifdef MPI
!    CALL MPI_Finalize( global%mpierr )
!#endif
    STOP
  ENDIF

  READ(casename ,*) global%casename
  READ(verbosity,*) global%verbLevel

! initial process for Genx ----------------------------------------------------

#ifdef GENX
  CALL COM_Init

! load Rocout module
  CALL COM_set_verbose( 1 )
  CALL SimOUT_load_module( 'OUT')
#endif

! read, check and print user input --------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading region topology ...'

  CALL RFLO_ReadRegionTopology( global,regions )

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Reading user input file  ...'

  CALL InitInputValues( regions )

  CALL ReadInputFile( regions )

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Checking BC consistency  ...'

  CALL ReadBcInputFile( regions )

  CALL CheckBcValidity( regions )

  IF (global%verbLevel >= VERBOSE_LOW) CALL PrintPrepInput( regions )

! check grid level and obtain volume grid -------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Reading grid  ...'

  gridLevel = global%startLevel

  DO iReg=1,global%nRegions
    IF (regions(iReg)%nGridLevels < gridLevel) THEN
      WRITE(msg,1000) SOLVER_NAME,iReg,gridLevel
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
    
    CALL GetGrid( gridLevel,iReg,regions )
  ENDDO

! create bc distribution files if applicable ----------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Write BC distribution to files ...'

  CALL BcDistributionFiles( regions )

! check for Genx format -------------------------------------------------------

  IF (global%solutFormat==2) THEN
#ifndef GENX
    CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__, &
         'compile with GENX=1 for format=2 (Genx HDF format)' )
#endif
  ENDIF

! generate initial guess (regionwise) -----------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Generating and storing solution ...'

  DO iReg=1,global%nregions
    WRITE(STDOUT,'(A,I5.5)') SOLVER_NAME//'   - region ',iReg

    IF (global%solutFormat==2) THEN
#ifdef GENX
      CALL GenxInitSolution( gridLevel,iReg,regions,wins,winv )
      CALL GenxWriteSolution( gridLevel,iReg,regions(iReg),wins,winv )
#endif
    ELSE
      CALL InitializeFlowField( gridLevel,regions(iReg) )
      CALL RFLO_WriteSolutionRegion( iReg,regions )

      DEALLOCATE( regions(iReg)%levels(gridLevel)%mixt%cv,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
    ENDIF   !  solFormat
  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

!#ifdef MPI
!  CALL MPI_Finalize( global%mpierr )
!#endif

1000 FORMAT(A,' Region ',I5,', grid level= ',I2,'.')

END PROGRAM ROCFLO_Init

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_Main.F90,v $
! Revision 1.7  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2005/05/03 03:20:01  wasistho
! enabled modified cyl.Taylor inflow profile
!
! Revision 1.4  2005/05/02 18:07:54  wasistho
! added cylindrical Taylor inflow profile capability
!
! Revision 1.3  2005/04/29 03:31:45  wasistho
! added distribution bc file generator
!
! Revision 1.2  2004/12/03 03:29:42  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.30  2004/11/16 07:16:32  jiao
! Updated to call Rocout_load_module directly instead of calling through
! COM_load_module, to allow static linking.
!
! Revision 1.29  2004/11/14 17:23:09  jiao
! Removed conditional compilation for CHARM.
!
! Revision 1.28  2004/11/14 05:42:16  jiao
! Changed to use TCHARM_GETARG when CHARM is defined.
!
! Revision 1.27  2004/09/08 05:18:30  wasistho
! commented ifdef MPI, as it is always one processor operation anyway
!
! Revision 1.26  2004/07/28 01:47:04  wasistho
! cosmetics and cleanup
!
! Revision 1.25  2004/07/27 20:28:13  wasistho
! added readBcInputFile and checkBcValidity
!
! Revision 1.24  2004/07/27 03:34:47  wasistho
! add printPrepInput
!
! Revision 1.23  2004/07/23 04:32:18  wasistho
! Genx: readin from Rocin, standalone: read .inp file i.o. command line input
!
! Revision 1.22  2004/06/30 00:00:19  wasistho
! migrated to Roccom-3
!
! Revision 1.21  2003/05/25 18:11:30  jiao
! Added support for Charm.
!
! Revision 1.20  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.19  2003/04/10 03:44:38  jblazek
! Merged .ini file into .inp file.
!
! Revision 1.18  2003/03/20 23:31:50  jiao
! ACH: Changed ModInterfaces to RFLO_ModInterfacesPrep.
!
! Revision 1.17  2003/03/20 19:44:22  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.16  2003/03/20 19:35:43  haselbac
! Modified RegFun call to avoid probs with long 'PREP_Main.F90' names
!
! Revision 1.15  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.14  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.13  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.12  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.11  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.10  2002/06/14 17:35:11  jblazek
! Added version string.
!
! Revision 1.9  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.8  2002/03/18 21:31:55  jblazek
! Made codes compatible with MPI.
!
! Revision 1.7  2002/02/21 23:25:07  jblazek
! Blocks renamed as regions.
!
! Revision 1.6  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.5  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.4  2002/01/10 18:21:30  jblazek
! Added iteration number and initial residual to solution file.
!
! Revision 1.3  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.2  2002/01/02 15:57:08  jblazek
! Added flow initialization.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************








