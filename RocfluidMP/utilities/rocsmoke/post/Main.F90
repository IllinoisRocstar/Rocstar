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
! Purpose: write out ROCFLO-MP`s grid and/or solution for visualization,
!          with smoke fields included
!
! Description: currently supported formats are:
!              - TECPLOT Ascii
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
! $Id: Main.F90,v 1.4 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

PROGRAM PEUL_ROCFLO_Post

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModPartEul, ONLY    : t_peul
  USE RFLO_ModInterfacesPost, ONLY : RFLO_ReadRegionTopology, &
        RFLO_InitInputValues, ReadInputFile, RFLO_DerivedInputValues, &
        RFLO_GetDimensDummyNodes, RFLO_GetDimensDummy, RFLO_GetNodeOffset, &
        RFLO_GetCellOffset, RFLO_ReadGridRegion, RFLO_ReadSolutionRegion, &
        MixtureProperties, WriteTecplotAscii, &
        RFLO_GenerateCoarseGrids, RFLO_CopyGeometryDummy, BuildVersionString
   USE RFLO_ModInterfacesPost, ONLY : PEUL_ReadSolutionRegion
#ifdef TURB
   USE RFLO_ModInterfacesPost, ONLY : TURB_ReadInputFile, &
                                      TURB_DerivedInputValues
#endif
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... loop variables
  INTEGER :: iReg, iLev

! ... local variables
  CHARACTER(CHRLEN) :: pltType, level, stamp, outFmt, msg, &
                       versionString, headerString, nPeulStr

  INTEGER :: plotType, outputFormat, nPeul
  INTEGER :: ipc, jpc, kpc, ibc, iec, ibn, ien
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 53

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)
  TYPE(t_grid)  , POINTER :: grid
  TYPE(t_mixt)  , POINTER :: mixt
  TYPE(t_peul)  , POINTER :: peul

!******************************************************************************

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'PEUL_ROCFLO_Post', &
                        'Main.F90' )

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
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        ROCSMOKE: Solution Postprocessing          *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *        =================================          *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! read argument list ----------------------------------------------------------

  CALL GETARG(1,global%casename)
  CALL GETARG(2,pltType)
  CALL GETARG(3,level)
  CALL GETARG(4,stamp)
  CALL GETARG(5,outFmt)
  CALL GETARG(6,nPeulStr)

  IF (LEN_TRIM(global%casename)==0 .OR. &
      LEN_TRIM(pltType)==0         .OR. &
      LEN_TRIM(level)==0           .OR. &
      LEN_TRIM(stamp)==0           .OR. &
      LEN_TRIM(outFmt)==0          .OR. &
      LEN_TRIM(nPeulStr)==0 ) THEN
    WRITE(STDOUT,'(/,A,/,A,/,9(A,/))') &
      SOLVER_NAME//' Usage: peulpost <casename> <type> <level> <time/iter> <format> <nPeul>', &
      SOLVER_NAME, &
      SOLVER_NAME//'        type      = 1 - grid + smoke only', &
      SOLVER_NAME//'                  = 2 - grid + fluid + smoke fields', &
      SOLVER_NAME//' ', &
      SOLVER_NAME//'        level     = grid level (>0)', &
      SOLVER_NAME//'        time/iter = time or iteration number', &
      SOLVER_NAME//' ', &
      SOLVER_NAME//'        format    = 3 - Tecplot ASCII', &
      SOLVER_NAME//'        nPeul     = number of smoke fields included'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF

  READ(pltType ,*) plotType
  READ(level   ,*) global%startLevel
  READ(outFmt  ,*) outputFormat
  READ(nPeulStr,*) nPeul

  IF (plotType <= 1) THEN
    plotType = PLOT_GRID_ONLY
  ELSE
    plotType = PLOT_GRID_FLOW
  ENDIF

  IF (outputFormat == 3) THEN
    outputFormat = PLOT_FMT_TECASCII
  ELSE
    WRITE(STDOUT,'(/,A,/)') SOLVER_NAME// &
                            ' Sorry, binary output not yet supported.'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF

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
    DO iLev=2,regions(iReg)%nGridLevels
      ipc = regions(iReg)%levels(iLev-1)%grid%ipc
      jpc = regions(iReg)%levels(iLev-1)%grid%jpc
      kpc = regions(iReg)%levels(iLev-1)%grid%kpc
      regions(iReg)%levels(iLev)%grid%ipc = ipc/2
      regions(iReg)%levels(iLev)%grid%jpc = jpc/2
      regions(iReg)%levels(iLev)%grid%kpc = kpc/2
    ENDDO
  ENDDO

! get user parameters ---------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading user input ...'

  CALL RFLO_InitInputValues( regions )
  CALL ReadInputFile( regions )
  CALL RFLO_DerivedInputValues( regions )
#ifdef TURB
  CALL TURB_ReadInputFile( regions )
  CALL TURB_DerivedInputValues( regions )
#endif

  IF (global%flowType == FLOW_STEADY) THEN
    READ(stamp,*) global%currentIter
  ELSE
    READ(stamp,*) global%timeStamp
    global%currentTime = global%timeStamp
  ENDIF

! loop over regions -----------------------------------------------------------

  IF (plotType == PLOT_GRID_ONLY) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Reading grid + smoke and writing plot file ...'
  ELSE
    WRITE(STDOUT,'(A)') SOLVER_NAME// &
                        ' Reading grid + fluid + smoke and writing plot file ...'
  ENDIF

  DO iReg=1,global%nRegions

    WRITE(STDOUT,'(A,I5.5)') SOLVER_NAME//'   - region ',iReg

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
    peul => regions(iReg)%levels(iLev)%peul

    CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
    iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

    CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

    IF (regions(iReg)%mixtInput%moveGrid) THEN
      ALLOCATE( grid%siVel(ibn:ien),stat=errorFlag )
      ALLOCATE( grid%sjVel(ibn:ien),stat=errorFlag )
      ALLOCATE( grid%skVel(ibn:ien),stat=errorFlag )
    ENDIF
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    IF (plotType == PLOT_GRID_FLOW) THEN
      ALLOCATE( mixt%cv(5,ibc:iec),stat=errorFlag )
      ALLOCATE( mixt%dv(mixt%nDv,ibc:iec),stat=errorFlag )
      IF (regions(iReg)%mixtInput%computeTv) THEN
        ALLOCATE( mixt%tv(mixt%nTv,ibc:iec),stat=errorFlag )
      ENDIF
      IF (regions(iReg)%mixtInput%gasModel == GAS_MODEL_TCPERF) THEN
        ALLOCATE( mixt%gv(mixt%nGv,0:1),stat=errorFlag )
      ELSE
        ALLOCATE( mixt%gv(mixt%nGv,ibc:iec),stat=errorFlag )
      ENDIF
    ENDIF

    peul%nCv = nPeul
    regions(iReg)%peulInput%nPtypes = nPeul
    IF (nPeul > 0) THEN
      ALLOCATE( peul%cv(nPeul,ibc:iec),stat=errorFlag )
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - read grid

    CALL RFLO_ReadGridRegion( iReg,regions )
    CALL RFLO_GenerateCoarseGrids( regions(iReg) )
    CALL RFLO_CopyGeometryDummy( regions(iReg) )

! - read solution, calc. mixture properties

    IF (plotType == PLOT_GRID_FLOW) THEN
      CALL RFLO_ReadSolutionRegion( iReg,regions )
      CALL MixtureProperties( regions(iReg),ibc,iec,.true. )
    ENDIF

! - read smoke solution

    IF (nPeul > 0) THEN
      WRITE(*,*) ' Entering PEUL_ReadSolutionRegion: iReg = ', iReg
      CALL PEUL_ReadSolutionRegion( iReg,regions )
      WRITE(*,*) ' Exiting  PEUL_ReadSolutionRegion: iReg = ', iReg
    ENDIF

    IF (outputFormat == PLOT_FMT_TECASCII) THEN

! --- write data to TECPLOT file (ASCII)

      CALL WriteTecplotAscii( iReg,iLev,plotType,regions(iReg) )

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
    IF (plotType == PLOT_GRID_FLOW) THEN
      DEALLOCATE( mixt%cv ,stat=errorFlag )
      DEALLOCATE( mixt%dv ,stat=errorFlag )
      DEALLOCATE( mixt%tv ,stat=errorFlag )
      DEALLOCATE( mixt%gv ,stat=errorFlag )
    ENDIF
    IF (nPeul > 0) DEALLOCATE( peul%cv ,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ENDDO   ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif

1000 FORMAT(A,' Region ',I5,', grid level= ',I2,'.')

END PROGRAM PEUL_ROCFLO_Post

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Main.F90,v $
! Revision 1.4  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:39  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/12/01 22:29:09  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2004/03/03 23:55:42  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.2  2003/09/26 22:51:04  jferry
! changed header printed out
!
! Revision 1.1  2003/09/25 15:40:22  jferry
! Implented Rocsmoke post-processing
!
!
!******************************************************************************







