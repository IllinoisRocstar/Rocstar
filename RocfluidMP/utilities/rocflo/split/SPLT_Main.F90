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
! Purpose: split up single region into multiple grids and generate
!          corresponding topology file.
!
! Description: region can be splitted either in i-, j-, or in k-direction.
!
! Input: file names, split direction, number of grids.
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: SPLT_Main.F90,v 1.4 2008/12/06 08:44:51 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

PROGRAM ROCFLO_Split

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE SPLT_ModInterfaces, ONLY : RFLO_ReadRegionTopology, &
           RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, RFLO_ReadGridRegion, &
           SplitGrid, SplitTopology, RFLO_WriteGridRegion, &
           RFLO_WriteRegionTopology, BuildVersionString
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: casenameOld, casenameNew, versionString, headerString

  INTEGER :: gridFormat, splitDirection, nSplits
  INTEGER :: strLen, ibn, ien, iNOff, ijNOff
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 53

  TYPE(t_global), POINTER :: globalOld, globalNew
  TYPE(t_region), POINTER :: regionsOld(:), regionsNew(:)

!******************************************************************************

  ALLOCATE( globalOld )
  ALLOCATE( globalNew )

  globalOld%nFunTree = 0
  CALL RegisterFunction( globalOld,'ROCFLO_Split',&
  'SPLT_Main.F90' )

! initialize global parameters ------------------------------------------------

  globalOld%verbLevel = VERBOSE_NONE

  globalOld%flowType    = FLOW_STEADY  ! stationary flow
  globalOld%currentTime = -1._RFREAL   ! no physical time set
  globalOld%currentIter = -1           ! no iteration

  globalOld%inDir  = './'              ! directory path
  globalOld%outDir = './'

  globalOld%nProcAlloc = 1
  globalOld%myProcid   = MASTERPROC    ! default process number (not an MPI code)
  globalOld%mpierr     = ERR_NONE
  globalOld%error      = ERR_NONE

  globalOld%pi  = 4._RFREAL*ATAN(1._RFREAL)
  globalOld%rad = globalOld%pi/180._RFREAL

! print header ----------------------------------------------------------------

#ifdef MPI
  CALL MPI_Init( globalOld%mpierr )
  IF (globalOld%mpierr /=0 ) &
    CALL ErrorStop( globalOld,ERR_MPI_TROUBLE,__LINE__ )
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
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *             ROCFLO-MP: Grid Split Up              *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *             ========================              *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! read user input -------------------------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Case name of old grid:'
  READ(STDIN,'(A)') casenameOld

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Case name of new grid:'
  READ(STDIN,'(A)') casenameNew

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Grid format (ASCII=1, binary=2):'
  READ(STDIN,*) gridFormat

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Split direction (i=1, j=2, k=3):'
  READ(STDIN,*) splitDirection

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Number of new regions:'
  READ(STDIN,*) nSplits

! figure out grid format

  IF (gridFormat <= 1) THEN
    globalOld%gridFormat = FORMAT_ASCII
  ELSE
    globalOld%gridFormat = FORMAT_BINARY
  ENDIF

! read old topology -----------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading old topology ...'

  globalOld%casename = casenameOld

  CALL RFLO_ReadRegionTopology( globalOld,regionsOld )

  IF (globalOld%nRegions > 1) THEN
    WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Tool works only for SINGLE region.'
#ifdef MPI
    CALL MPI_Finalize( globalOld%mpierr )
#endif
    STOP
  ENDIF

! read old grid ---------------------------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Reading old grid ...'

! allocate memory

  regionsOld(1)%nDumCells = 0
  CALL RFLO_GetDimensDummyNodes( regionsOld(1),1,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( regionsOld(1),1,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  ALLOCATE( regionsOld(1)%levels(1)%grid%xyz(3,ibn:ien),stat=errorFlag )
  globalOld%error = errorFlag
  IF (globalOld%error /= 0) CALL ErrorStop( globalOld,ERR_ALLOCATE,__LINE__ )

! read grid

  CALL RFLO_ReadGridRegion( 1,regionsOld )

! split the grid --------------------------------------------------------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Splitting the grid ...'

  globalNew = globalOld

  globalNew%casename = casenameNew
  globalNew%nRegions = nSplits

  ALLOCATE( regionsNew(globalNew%nRegions),stat=errorFlag )
  globalOld%error = errorFlag
  IF (globalOld%error /= 0) CALL ErrorStop( globalOld,ERR_ALLOCATE,__LINE__ )

! set global within regionsNew = globalNew

  DO iReg=1,globalNew%nRegions
    regionsNew(iReg)%global => globalNew
  ENDDO

! do the splitting

  CALL SplitGrid( splitDirection,regionsOld,regionsNew )

  CALL SplitTopology( splitDirection,regionsOld,regionsNew )

! write new grid and topology

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Writing data ...'

  DO iReg=1,globalNew%nRegions
    CALL RFLO_WriteGridRegion( iReg,regionsNew )
  ENDDO

  CALL RFLO_WriteRegionTopology( regionsNew )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( globalOld )

  WRITE(STDOUT,'(/,A,/)') SOLVER_NAME//' Finished.'

#ifdef MPI
  CALL MPI_Finalize( globalOld%mpierr )
#endif

END PROGRAM ROCFLO_Split

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPLT_Main.F90,v $
! Revision 1.4  2008/12/06 08:44:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:33:15  wasistho
! rflo_modinterfacessplit to splt_modinterfaces
!
! Revision 1.1  2004/12/03 02:41:15  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:45:30  wasistho
! lower to upper case
!
! Revision 1.7  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.6  2003/03/20 22:32:37  haselbac
! Renamed ModInterfaces
!
! Revision 1.5  2003/03/20 19:45:54  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.4  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.3  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/07/12 21:50:08  jblazek
! Added tool to split single grid into multiple regions.
!
!******************************************************************************







