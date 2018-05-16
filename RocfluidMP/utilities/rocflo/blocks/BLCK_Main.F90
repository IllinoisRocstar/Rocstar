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
! Purpose: write out numbers of blocks being contained in a given
!          box in the physical space.
!
! Description: none.
!
! Input: case name from the list of arguments
!
! Output: to screen.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: BLCK_Main.F90,v 1.3 2008/12/06 08:44:46 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

#ifdef CHARM
SUBROUTINE MPI_Main
#else
PROGRAM ROCFLO_Blocks
#endif

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE BLCK_ModInterfaces, ONLY : RFLO_ReadRegionTopology, &
        RFLO_GetDimensDummyNodes, RFLO_GetDimensPhysNodes, &
        RFLO_GetNodeOffset, RFLO_ReadGridRegion, BuildVersionString
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... loop variables
  INTEGER :: iReg, i, j, k, ijkN

! ... local variables
  CHARACTER(CHRLEN) :: gridFmt, msg, xmin, xmax, ymin, ymax, &
                       zmin, zmax, versionString, headerString

  INTEGER :: ibn, ien, iNOff, ijNOff, regNum, regNumBeg, regNumEnd
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 53

  REAL(RFREAL) :: xminBox, xmaxBox, yminBox, ymaxBox, zminBox, zmaxBox
  REAL(RFREAL), POINTER :: xyz(:,:)
  
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)

!******************************************************************************

  ALLOCATE( global )

  global%nFunTree = 0
  CALL RegisterFunction( global,'ROCFLO_Blocks',&
  'BLCK_Main.F90' )

! initialize global parameters ------------------------------------------------

  global%verbLevel = VERBOSE_NONE

  global%flowType    = FLOW_STEADY  ! stationary flow
  global%currentTime = -1._RFREAL   ! no physical time set
  global%currentIter = -1           ! no iteration
  global%startLevel  = 1            ! always finest grid only

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
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *         ROCFLO-MP: Block Counting Utility         *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *         =================================         *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') SOLVER_NAME//'  *                                                   *'
  WRITE(STDOUT,'(A,/)') SOLVER_NAME//'  *****************************************************'

! read user input ----------------------------------------------------------

  CALL GETARG(1,global%casename)
  CALL GETARG(2,gridFmt)
  CALL GETARG(3,xmin)
  CALL GETARG(4,xmax)
  CALL GETARG(5,ymin)
  CALL GETARG(6,ymax)
  CALL GETARG(7,zmin)
  CALL GETARG(8,zmax)

  IF (LEN_TRIM(global%casename)==0 .OR. &
      LEN_TRIM(gridFmt)==0         .OR. &
      LEN_TRIM(xmin)==0            .OR. &
      LEN_TRIM(xmax)==0            .OR. &
      LEN_TRIM(ymin)==0            .OR. &
      LEN_TRIM(ymax)==0            .OR. &
      LEN_TRIM(zmin)==0            .OR. &
      LEN_TRIM(zmax)==0) THEN
    WRITE(STDOUT,'(/,A,/,A,/,2(A,/))') &
      SOLVER_NAME//' Usage: rfloblocks <casename> <format> <xmin xmax> <ymin ymax> <zmin zmax>', &
      SOLVER_NAME, &
      SOLVER_NAME//'        format = 0 - ASCII grid', &
      SOLVER_NAME//'                 1 - binary grid'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF

  READ(gridFmt ,*) global%gridFormat
  READ(xmin    ,*) xminBox
  READ(xmax    ,*) xmaxBox
  READ(ymin    ,*) yminBox
  READ(ymax    ,*) ymaxBox
  READ(zmin    ,*) zminBox
  READ(zmax    ,*) zmaxBox

  IF (global%gridFormat <= 0) THEN
    global%gridFormat = FORMAT_ASCII
  ELSE
    global%gridFormat = FORMAT_BINARY
  ENDIF

! read region topology --------------------------------------------------------

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Reading region topology ...'

  CALL RFLO_ReadRegionTopology( global,regions )

  DO iReg=1,global%nRegions
    regions(iReg)%startLevel = global%startLevel
    regions(iReg)%currLevel  = global%startLevel
  ENDDO

! loop over regions -----------------------------------------------------------

  WRITE(STDOUT,'(/,A,/)') SOLVER_NAME//' Searching for blocks ...'

  regNumBeg = -99
  regNumEnd = -99

  DO iReg=1,global%nRegions

! - allocate memory for the grid (all grid levels)

    CALL RFLO_GetDimensDummyNodes( regions(iReg),1,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetDimensPhysNodes( regions(iReg),1,ipnbeg,ipnend, &
                                  jpnbeg,jpnend,kpnbeg,kpnend )
    CALL RFLO_GetNodeOffset( regions(iReg),1,iNOff,ijNOff )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
    ALLOCATE( regions(iReg)%levels(1)%grid%xyz(3,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - read grid

    CALL RFLO_ReadGridRegion( iReg,regions )

! - find regions contained in the box

    xyz => regions(iReg)%levels(1)%grid%xyz

    regNum = 0
    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend
        DO i=ipnbeg,ipnend
          ijkN = IndIJK(i,j,k,iNOff,ijNOff)
          IF ((xyz(XCOORD,ijkN)>xminBox .AND. xyz(XCOORD,ijkN)<xmaxBox) .AND. &
              (xyz(YCOORD,ijkN)>yminBox .AND. xyz(YCOORD,ijkN)<ymaxBox) .AND. &
              (xyz(ZCOORD,ijkN)>zminBox .AND. xyz(ZCOORD,ijkN)<zmaxBox)) THEN
            regNum = iReg
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    IF (regNum == regNumEnd+1) THEN
      regNumEnd = regNum
    ELSE
      IF (regNum/=0 .AND. regNumEnd<0) THEN
        regNumBeg = regNum
        regNumEnd = regNum
      ELSE IF (regNum/=0 .AND. regNumEnd>0) THEN
        WRITE(STDOUT,1000) SOLVER_NAME,regNumBeg,regNumEnd
        regNumBeg = regNum
        regNumEnd = regNum
      ENDIF
    ENDIF

! - deallocate memory

    DEALLOCATE( regions(iReg)%levels(1)%grid%xyz,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  ENDDO   ! iReg

  IF (regNumBeg>0 .AND. regNumEnd>0) THEN
    WRITE(STDOUT,1000) SOLVER_NAME,regNumBeg,regNumEnd
  ELSE
    WRITE(STDOUT,1005) SOLVER_NAME
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') SOLVER_NAME//' Finished.'

#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif

1000 FORMAT(A,' From block ',I5,' to ',I5)
1005 FORMAT(A,' No blocks found.')

#ifdef CHARM
END SUBROUTINE MPI_Main
#else
END PROGRAM ROCFLO_Blocks
#endif

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BLCK_Main.F90,v $
! Revision 1.3  2008/12/06 08:44:46  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/03 03:13:06  wasistho
! bloc to blck
!
! Revision 1.1  2004/12/03 01:55:08  wasistho
! add prefix
!
! Revision 1.1  2004/12/03 00:26:14  wasistho
! lower to upper case
!
! Revision 1.6  2003/05/25 18:11:30  jiao
! Added support for Charm.
!
! Revision 1.5  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/20 22:20:21  haselbac
! Renamed ModInterfaces
!
! Revision 1.3  2003/03/20 19:39:01  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.2  2003/03/20 19:33:08  haselbac
! Modified RegFun call to avoid probs with long 'BLCK_Main.F90' names
!
! Revision 1.1  2002/12/20 19:38:02  jblazek
! Added tool to count blocks in a box.
!
!******************************************************************************








