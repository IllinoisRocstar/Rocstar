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
! Purpose: conversion of structured grid into its equivalent hex unstructured.
!
! Description: the conversion is from RFLO format to RFLU format.
!
! Input: case name (and verbosity) from screen input, others from .inp file.
!
! Output: [case name].grda_00000
!
! Notes: none
!
!******************************************************************************
!
! $Id: TFLU_Main.F90,v 1.6 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

PROGRAM ROCFLO_toFlu

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE TFLU_ModInterfaces, ONLY : BuildVersionString, &
                               RFLO_ReadRegionTopology, ReadInputFile, &
                               PrintTofluInput, ConvertFlo2FluMesh, &
                               ConvertFlo2FluPatch, GetBndVertType, &
                               CorrectNedges, WriteFluCellMap, &
                               WriteFluGrid, WriteFluDimens
!  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... loop variables
  INTEGER :: iReg, iLev

! ... local variables
  CHARACTER(CHRLEN) :: casename, verbosity, msg, versionString, headerString

  INTEGER :: ipc, jpc, kpc, ni, nj, nk, gridLevel
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER :: nPatches, nVerts, nHexs, nBfMax, nBnMax
  INTEGER, PARAMETER :: headerWidth = 53

  TYPE(t_global), POINTER :: global
  TYPE(t_grid) , POINTER  :: grid
  TYPE(t_region), POINTER :: regions(:)

!******************************************************************************

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'ROCFLO_Init',&
  'TFLU_Main.F90' )

! initialize global parameters ------------------------------------------------

  global%verbLevel = VERBOSE_NONE

  global%flowType    = FLOW_STEADY  ! stationary flow
  global%currentTime = 0._RFREAL    ! no physical time set
  global%currentIter = 0            ! no iteration yet
  global%resInit     = 1._RFREAL

  global%inDir  = './'              ! directory path
  global%outDir = './'

  global%nProcAlloc = 1
  global%myProcid   = 0             ! default process number (if not MPI)
  global%mpierr     = ERR_NONE
  global%error      = ERR_NONE

  global%startLevel  = 1 
  global%gridFormat  = FORMAT_ASCII 
  global%solutFormat = FORMAT_ASCII 

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
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *          ROCFLO-MP: Conversion to ROCFLU          *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *          ===============================          *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! print required input and read argument list ---------------------------------

  WRITE(STDOUT,'(/,A,/,A,/,8(A,/))') &
    SOLVER_NAME//' Required Rocflo2Flu input:', &
    SOLVER_NAME, &
    SOLVER_NAME//'   <casename>        =>  command line input', &
    SOLVER_NAME//'   <verbosity>       =>  command line input', &
    SOLVER_NAME//'   grid level        => .inp file: # MULTIGRID: START', &
    SOLVER_NAME//'   Rocflo grd format => .inp file: # FORMATS: GRID', &
    SOLVER_NAME//'   Rocflu grd format => .inp file: # FORMATS: SOLUTION'

  CALL GETARG(1,casename)
  CALL GETARG(2,verbosity)

  IF (LEN_TRIM(casename)==0 .OR. &
      LEN_TRIM(verbosity)==0) THEN
    WRITE(STDOUT,'(/,A,/)') &
      SOLVER_NAME//' Usage: rflo2flu <casename> <verbosity>'
!#ifdef MPI
!    CALL MPI_Finalize( global%mpierr )
!#endif
    STOP
  ENDIF

  READ(casename ,*) global%casename
  READ(verbosity,*) global%verbLevel

! read, check and print user input -------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading region topology ...'

  CALL RFLO_ReadRegionTopology( global,regions )

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Reading user input file  ...'

  CALL ReadInputFile( regions )

  IF (global%verbLevel >= VERBOSE_LOW) CALL PrintTofluInput( regions )

  gridLevel = global%startLevel

! check grid level and obtain region grid size --------------------------------

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
  ENDDO

! start conversion multi regions Rocflo to single region Rocflu --------------
! obtain sizes of global unstructured grid arrays at whole domain, all patches

  global%tofluNPatches = 0
  global%tofluNHexs    = 0
  global%tofluNVerts   = 0
  global%tofluNbfMax   = 0
  global%tofluNbnMax   = 0
  global%tofluNFaces   = 0
  global%tofluNEdges   = 0

  DO iReg=1,global%nregions
    regions(iReg)%currLevel = gridLevel
    CALL ConvertFlo2FluMesh( 0,iReg,regions )
  ENDDO   ! iReg

! allocate unstructured coordinates and connectivities -----------------------

  nPatches = global%tofluNPatches
  nVerts   = global%tofluNVerts
  nHexs    = global%tofluNHexs
  nBfMax   = global%tofluNbfMax
  nBnMax   = global%tofluNbnMax

  ALLOCATE( global%tofluXyz(3,nVerts),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ALLOCATE( global%tofluNbVerts(nPatches),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ALLOCATE( global%tofluNbFaces(nPatches),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ALLOCATE( global%tofluHex2v(8,nHexs),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ALLOCATE( global%tofluQuad2v(4,nBfMax,nPatches),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ALLOCATE( global%tofluBLoc2g(nBnMax,nPatches),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ALLOCATE( global%tofluIq(nPatches),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! generate and store Rocflu files --------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Generating and storing Rocflu files ...'

  global%tofluNPatches   = 0
  global%tofluNHexs      = 0
  global%tofluNVerts     = 0
  global%tofluNbVerts(:) = 0
  global%tofluNbFaces(:) = 0
  global%tofluIq(:)      = 0
  global%tofluMaxBind    = 0

  DO iReg=1,global%nregions
    WRITE(STDOUT,'(A,I5.5)') SOLVER_NAME//'   - region ',iReg

    iLev =  regions(iReg)%currLevel
    grid => regions(iReg)%levels(iLev)%grid

    CALL ConvertFlo2FluMesh( 1,iReg,regions )

    CALL ConvertFlo2FluPatch( iReg,regions )

! - note, only grid%xyz is deallocated, grid%tofluLoc2g is still needed

    DEALLOCATE( grid%xyz,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ENDDO   ! iReg

! correct number of edges per global patch

  ALLOCATE( global%tofluBType(6,global%tofluMaxBind),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  global%tofluBType(:,:) = 0  

  DO iReg=1,global%nregions
    CALL GetBndVertType( iReg,regions )
  ENDDO

  CALL CorrectNedges( global )

! write Rocflu files

  CALL WriteFluGrid( global )
  CALL WriteFluDimens( global )
  CALL WriteFluCellMap( global )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

!#ifdef MPI
!  CALL MPI_Finalize( global%mpierr )
!#endif

1000 FORMAT(A,' Region ',I5,', grid level= ',I2,'.')

END PROGRAM ROCFLO_toFlu

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_Main.F90,v $
! Revision 1.6  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/03/07 03:23:32  wasistho
! enabled serial execution
!
! Revision 1.3  2005/12/21 22:38:02  wasistho
! added writeFluCellMap
!
! Revision 1.2  2004/12/03 03:44:10  wasistho
! rflo_modinterfacestoflu to tflu_modinterfaces
!
! Revision 1.1  2004/12/03 02:59:30  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:58:20  wasistho
! lower to upper case
!
! Revision 1.3  2004/08/18 02:14:21  wasistho
! removed explicit dimensions in initializing tofluBType
!
! Revision 1.2  2004/08/18 02:10:20  wasistho
! added new routines to create dimension file
!
! Revision 1.1.1.1  2004/08/17 01:41:39  wasistho
! initial checkin
!
!
!******************************************************************************







