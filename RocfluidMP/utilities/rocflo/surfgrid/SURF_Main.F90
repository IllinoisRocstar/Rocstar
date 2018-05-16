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
! Purpose: provide grid coordinates for all surfaces which interact
!          with GenX.
!
! Description: none.
!
! Input: case name from the list of arguments
!
! Output: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: SURF_Main.F90,v 1.4 2008/12/06 08:44:52 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

#ifdef CHARM
SUBROUTINE MPI_Main
#else
PROGRAM ROCFLO_Surf
#endif

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE SURF_ModInterfaces, ONLY : BuildVersionString, &
          RFLO_ReadRegionTopology, RFLO_ReadBcInputFile, &
          RFLO_CopyTopologyLevels, RFLO_ReadGridRegion, &
          RFLO_GenerateCoarseGrids, RFLO_CopyGeometryDummy, &
          CountInteractingPatches, WriteSurfaceGrid
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... loop variables
  INTEGER :: iReg, iLev

! ... local variables
  CHARACTER(CHRLEN) :: level, gridFormat, msg, versionString, headerString

  INTEGER :: ipc, jpc, kpc, ibn, ien, nInteract
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, iNOff, ijNOff
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 53

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)

!******************************************************************************

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'ROCFLO_Surf',&
  'SURF_Main.F90' )

! initialize global parameters ------------------------------------------------

  global%verbLevel = VERBOSE_NONE

  global%flowType    = FLOW_STEADY  ! stationary flow
  global%currentTime = 0._RFREAL    ! no physical time set
  global%timeStamp   = 0._RFREAL
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
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *         ROCFLO-MP: Surface Grid for GenX          *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *         ================================          *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! read argument list ----------------------------------------------------------

  CALL GETARG(1,global%casename)
  CALL GETARG(2,level)
  CALL GETARG(3,gridFormat)

  IF (LEN_TRIM(global%casename)==0 .OR. &
      LEN_TRIM(level)==0           .OR. &
      LEN_TRIM(gridFormat)==0) THEN
    WRITE(STDOUT,'(/,A,/,5(A,/))')  &
      SOLVER_NAME//' Usage: rflosurf <casename> <level> <grid>', &
      SOLVER_NAME, &
      SOLVER_NAME//'        level = grid level (>0)', &
      SOLVER_NAME, &
      SOLVER_NAME//'        grid  = 0 - ASCII format', &
      SOLVER_NAME//'              = 1 - binary format'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF

  READ(level     ,*) global%startLevel
  READ(gridFormat,*) global%gridFormat

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

! read boundary conditions

  CALL RFLO_ReadBcInputFile( regions )

! copy topology and BCs to all grid levels

  CALL RFLO_CopyTopologyLevels( regions )

! count number of interacting patches

  CALL CountInteractingPatches( regions,nInteract )

! output surface grid (regionwise) --------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Generating and storing surface grid ...'

! open file

  OPEN(IF_PLOT,file=TRIM(global%casename)//'.im',form='formatted', &
       status='unknown',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__ )

  WRITE(IF_PLOT,*) nInteract,'  2'

! loop over all regions

  DO iReg=1,global%nregions
    WRITE(STDOUT,'(A,I5.5)') SOLVER_NAME//'   - region ',iReg

    regions(iReg)%currLevel = global%startLevel

! - allocate memory for grid

    DO iLev=1,regions(iReg)%nGridLevels
      CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                     jdnbeg,jdnend,kdnbeg,kdnend )
      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
      ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
      ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
      ALLOCATE( regions(iReg)%levels(iLev)%grid%xyz(3,ibn:ien),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ENDDO

! - read grid

    CALL RFLO_ReadGridRegion( iReg,regions )
    CALL RFLO_GenerateCoarseGrids( regions(iReg) )
    CALL RFLO_CopyGeometryDummy( regions(iReg) )

! - write out surface grid

    iLev = global%startLevel

    CALL WriteSurfaceGrid( iReg,regions(iReg) )

! - deallocate memory

    DO iLev=1,regions(iReg)%nGridLevels
      DEALLOCATE( regions(iReg)%levels(iLev)%grid%xyz,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
    ENDDO
  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  CLOSE(IF_PLOT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__ )

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif

1000 FORMAT(A,' Region ',I5,', grid level= ',I2,'.')

#ifdef CHARM
END SUBROUTINE MPI_Main
#else
END PROGRAM ROCFLO_Surf
#endif

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SURF_Main.F90,v $
! Revision 1.4  2008/12/06 08:44:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:35:43  wasistho
! rflo_modinterfacessurf to surf_modinterfaces
!
! Revision 1.1  2004/12/03 02:47:00  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:49:09  wasistho
! lower to upper case
!
! Revision 1.5  2003/05/25 18:11:30  jiao
! Added support for Charm.
!
! Revision 1.4  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.3  2003/03/20 22:35:02  haselbac
! Renamed ModInterfaces
!
! Revision 1.2  2003/03/20 19:48:09  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.1  2002/10/19 00:40:31  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
!******************************************************************************







