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
! Purpose: read File for Stream 2 Data based on Rocflo Computations.
!
! Description: none.
!
! Input:
!
! Output: Memory location with data for Stream 2
!
! ISSUE: Where do you read the Analytical and Experimental Data on Nodes or Cells
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ReadFileStream2Comput.F90,v 1.4 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ReadFileStream2Comput ( regionsS1, regionsS2 )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModInterfaces, ONLY : RFLO_ReadRegionTopology, RFLO_InitInputValues, &
        ReadInputFile, RFLO_DerivedInputValues, RFLO_GetDimensDummyNodes, &
        RFLO_GetDimensDummy, RFLO_GetNodeOffset, RFLO_GetCellOffset, &
        RFLO_ReadGridRegion, RFLO_ReadSolutionRegion, MixtureProperties
  USE ModMPI
  USE ModParameters
  USE RVAV_ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModDataStruct
  IMPLICIT NONE

#include "Indexing.h"

! ... parameter variables
  TYPE (t_region), POINTER :: regionsS1(:), regionsS2(:)

! ... loop variables
  INTEGER :: iReg, iLev
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: msg, fname

  INTEGER :: iNodes, jNodes, kNodes, nVars, iVars
  INTEGER :: ipc, jpc, kpc, ibc, iec, ibn, ien
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: nRegionsS2
  INTEGER :: ijkN0, ijkC0
  INTEGER :: nGridLevels, ivar, errorFlag

  REAL(RFREAL), DIMENSION(3)              :: xyzMin, xyzMax
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: cvMin , cvMax
  REAL(RFREAL), POINTER                   :: cv(:,:), xyz(:,:)

  TYPE(t_grid),   POINTER :: grid
  TYPE(t_mixt),   POINTER :: mixt
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regionsS1(1)%global

  CALL RegisterFunction( global, 'RVAV_ReadFileStream2Comput',&
  'RVAV_ReadFileStream2Comput.F90' )

! read region topology --------------------------------------------------------

  IF ( global%verbLevel/=VERBOSE_NONE ) &
    WRITE(STDOUT,'(/,A)') 'RFLO Reading region topology for Stream 2...'

  global%casename = TRIM(globalRVAV%casename)
  CALL RFLO_ReadRegionTopology( global, regionsS2 )

  IF ( global%verbLevel/=VERBOSE_NONE ) &
    WRITE(STDOUT,'(/,A)') 'RFLO Reading user input for Stream 2...'

  CALL RFLO_InitInputValues( regionsS2 )
  CALL ReadInputFile( regionsS2 )
  CALL RFLO_DerivedInputValues( regionsS2 )

  DO iReg=1,global%nRegions
    regionsS2(iReg)%startLevel = global%startLevel
    regionsS2(iReg)%currLevel  = global%startLevel

    regionsS2(iReg)%mixtInput%computeTv = regionsS1(iReg)%mixtInput%computeTv
    regionsS2(iReg)%mixtInput%flowModel = regionsS1(iReg)%mixtInput%flowModel
    regionsS2(iReg)%mixtInput%gasModel = regionsS1(iReg)%mixtInput%gasModel

    IF (regionsS2(iReg)%nGridLevels < regionsS2(iReg)%currLevel) THEN
      WRITE(msg,1000) iReg,global%startLevel
      CALL ErrorStop( global, ERR_GRID_LEVEL,__LINE__,msg )
    ENDIF

  ENDDO ! iReg

! loop over regions -----------------------------------------------------------
  IF ( global%verbLevel/=VERBOSE_NONE ) &
    WRITE(STDOUT,'(/,A)') 'Reading grid and solution from Stream 2 - COMPUTED...'

! - reload casename for grid and solution files

  global%casename = TRIM(globalRVAV%casename)//'_s2'

  DO iReg=1,global%nRegions

    WRITE(STDOUT,'(A,I5.5)') '  - region ',iReg

    iLev =  regionsS2(iReg)%currLevel
    grid => regionsS2(iReg)%levels(iLev)%grid
    mixt => regionsS2(iReg)%levels(iLev)%mixt

! - allocate memory

    CALL RFLO_GetDimensDummy( regionsS2(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regionsS2(iReg),iLev,iCOff,ijCOff )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
    iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

    CALL RFLO_GetDimensDummyNodes( regionsS2(iReg),iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( regionsS2(iReg),iLev,iNOff,ijNOff )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

    globalRVAV%iCOffS2  = iCOff
    globalRVAV%ijCOffS2 = ijCOff

    IF ( global%verbLevel/=VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(/,A,2(I10,2X))') 'ibn-ien ',ibn,ien
      WRITE(STDOUT,'(A,3(I8,2X))')   'idnbeg,jdnbeg,kdnbeg',idnbeg,jdnbeg,kdnbeg
      WRITE(STDOUT,'(A,3(I8,2X))')   'idnend,jdnend,kdnend',idnend,jdnend,kdnend

      WRITE(STDOUT,'(/,A,2(I10,2X))') 'ibc-iec ',ibc,iec
      WRITE(STDOUT,'(A,3(I8,2X))')   'idcbeg,jdcbeg,kdcbeg',idcbeg,jdcbeg,kdcbeg
      WRITE(STDOUT,'(A,3(I8,2X))')   'idcend,jdcend,kdcend',idcend,jdcend,kdcend
      WRITE(STDOUT,'(/,A,2(I8,2X))') 'iCOffS2-ijCOffS2 ',globalRVAV%iCOffS2,globalRVAV%ijCOffS2

      WRITE(STDOUT,'(A,I5)')'nTv = ',mixt%nTv
      WRITE(STDOUT,'(A,I5)')'nGv = ',mixt%nGv
      WRITE(STDOUT,'(A,I5)')'nDV = ',mixt%nDv
      WRITE(STDOUT,'(A,L1)')'computeTv = ',regionsS2(iReg)%mixtInput%computeTv
      WRITE(STDOUT,'(A,I5)')'flowModel = ',regionsS2(iReg)%mixtInput%flowModel
      WRITE(STDOUT,'(A,I5)')'gasModel = ',regionsS2(iReg)%mixtInput%gasModel
    ENDIF ! verbLevel

    ALLOCATE( grid%xyz(3,ibn:ien),stat=errorFlag )
    IF (regionsS2(iReg)%mixtInput%moveGrid) THEN
      ALLOCATE( grid%siVel(ibn:ien),stat=errorFlag )
      ALLOCATE( grid%sjVel(ibn:ien),stat=errorFlag )
      ALLOCATE( grid%skVel(ibn:ien),stat=errorFlag )
    ENDIF
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

    ALLOCATE( mixt%cv(5,ibc:iec) ,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

    ALLOCATE( mixt%dv(mixt%nDv,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

    IF (regionsS2(iReg)%mixtInput%computeTv) THEN
      ALLOCATE( mixt%tv(mixt%nTv,ibc:iec),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
    ENDIF

    IF (regionsS2(iReg)%mixtInput%gasModel == GAS_MODEL_TCPERF) THEN
      ALLOCATE( mixt%gv(mixt%nGv,0:1),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
    ELSE
      ALLOCATE( mixt%gv(mixt%nGv,ibc:iec),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
    ENDIF

! - initialize variables

    grid%xyz = 0.0_RFREAL
    mixt%cv  = 0.0_RFREAL

! - read grid

    CALL RFLO_ReadGridRegion( iReg,regionsS2 )

! - read solution, calc. mixture properties

    CALL RFLO_ReadSolutionRegion( iReg,regionsS2 )

    CALL MixtureProperties( regionsS2(iReg),ibc,iec,.true. )

! - Write Min-Max Values

    IF (global%verbLevel /= VERBOSE_NONE) THEN

      WRITE(STDOUT,'(/,A,2(I5,2X))')'ibn-ien ',ibn,ien
      WRITE(STDOUT,'(A,3(I5,2X))')'idnbeg,jdnbeg,kdnbeg',idnbeg,jdnbeg,kdnbeg
      WRITE(STDOUT,'(A,3(I5,2X))')'idnend,jdnend,kdnend',idnend,jdnend,kdnend

      xyzMin = +1.0E+30_RFREAL
      xyzMax = -1.0E+30_RFREAL

      DO i=ibn,ien
        xyzMin(1) = MIN(grid%xyz(1,i), xyzMin(1))
        xyzMin(2) = MIN(grid%xyz(2,i), xyzMin(2))
        xyzMin(3) = MIN(grid%xyz(3,i), xyzMin(3))

        xyzMax(1) = MAX(grid%xyz(1,i), xyzMax(1))
        xyzMax(2) = MAX(grid%xyz(2,i), xyzMax(2))
        xyzMax(3) = MAX(grid%xyz(3,i), xyzMax(3))
      ENDDO ! i

      DO k=kdnbeg+3,kdnbeg+3
        DO j=jdnbeg+3,jdnbeg+3
          DO i=idnbeg,idnend
            ijkN0 = IndIJK(i,j,k,iNOff,ijNOff)
            WRITE(STDOUT,'(I5,3(3X,E12.5))')i,grid%xyz(1,ijkN0),grid%xyz(2,ijkN0),grid%xyz(3,ijkN0)
          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k

      WRITE(STDOUT,'(/,A,2(E12.5,2X))') 'Min-Max of X     ', xyzMin(1),xyzMax(1)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of Y     ', xyzMin(2),xyzMax(2)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of Z     ', xyzMin(3),xyzMax(3)

      WRITE(STDOUT,'(/,A,2(I5,2X))') 'ibc-iec ',ibc,iec
      WRITE(STDOUT,'(A,3(I5,2X))')   'idcbeg,jdcbeg,kdcbeg',idcbeg,jdcbeg,kdcbeg
      WRITE(STDOUT,'(A,3(I5,2X))')   'idcend,jdcend,kdcend',idcend,jdcend,kdcend

      cvMin = +1.0E+30_RFREAL
      cvMax = -1.0E+30_RFREAL

      DO i=ibc,iec
        cvMin(1) = MIN(mixt%cv(1,i), cvMin(1))
        cvMin(2) = MIN(mixt%cv(2,i), cvMin(2))
        cvMin(3) = MIN(mixt%cv(3,i), cvMin(3))
        cvMin(4) = MIN(mixt%cv(4,i), cvMin(4))
        cvMin(5) = MIN(mixt%cv(5,i), cvMin(5))

        cvMax(1) = MAX(mixt%cv(1,i), cvMax(1))
        cvMax(2) = MAX(mixt%cv(2,i), cvMax(2))
        cvMax(3) = MAX(mixt%cv(3,i), cvMax(3))
        cvMax(4) = MAX(mixt%cv(4,i), cvMax(4))
        cvMax(5) = MAX(mixt%cv(5,i), cvMax(5))
      ENDDO ! i

      DO k=kdcbeg,kdcbeg
        DO j=jdcbeg,jdcbeg
          DO i=idcbeg,idcend
            ijkC0 = IndIJK(i,j,k,iCOff,ijCOff)
            WRITE(STDOUT,'(/,A,I5,3(E12.5,2X))') i,mixt%cv(1,ijkC0),mixt%cv(2,ijkC0),mixt%cv(3,ijkC0)
          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k

      WRITE(STDOUT,'(/,A,2(E12.5,2X))') 'Min-Max of CV(1) ', cvMin(1),cvMax(1)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of CV(2) ', cvMin(2),cvMax(2)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of CV(3) ', cvMin(3),cvMax(3)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of CV(4) ', cvMin(4),cvMax(4)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of CV(5) ', cvMin(5),cvMax(5)

    ENDIF ! verbLevel

  ENDDO   ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', grid level= ',I2,'.')

END SUBROUTINE RVAV_readFileStream2Comput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ReadFileStream2Comput.F90,v $
! Revision 1.4  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:39  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/12/01 22:46:35  fnajjar
! Initial revision after changing case
!
! Revision 1.11  2004/03/03 23:55:42  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.10  2003/11/20 16:40:41  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.5  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.4  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.3  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.2  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.1  2002/07/16 22:34:46  f-najjar
! Initial Import
!
!******************************************************************************







