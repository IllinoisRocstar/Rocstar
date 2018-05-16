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
! Purpose: read in x-,y-,z-coordinates of grid nodes.
!
! Description: the following grid formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions = dimensions of all regions.
!
! Output: region%levels%grid%xyz = grid coordinates
!         global%currentTime     = physical time
!
! Notes: additionaly, physical time is read in, which is stored at
!        the beginning of the file (meaningful only in the case of unsteady
!        flow and dynamic grids). It is the finest grid, which is read in.
!
!******************************************************************************
!
! $Id: RFLO_ReadGrid.F90,v 1.5 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadGrid( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal, &
                            RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k, n

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iRegFile, ipc, jpc, kpc, nDim
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iNOff, ijNOff, ijkN, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  LOGICAL :: moveGrid

  REAL(RFREAL), POINTER     :: xyz(:,:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), xyzFile(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadGrid',&
  'RFLO_ReadGrid.F90' )

! allocate fixed-size temporary data arrays -----------------------------------

  ALLOCATE( ivar(4,1),stat=errorFlag )
  ALLOCATE( rvar(1,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open grid file (only master proc.) ------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    moveGrid = .false.
    DO iReg=1,global%nRegions
      IF (regions(iReg)%mixtInput%moveGrid) moveGrid = .true.
    ENDDO

! - unsteady flow

    IF (global%flowType==FLOW_UNSTEADY .AND. &
        moveGrid .AND. global%timeStamp>0._RFREAL) THEN
#ifdef GENX
      IF (global%gridFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.grda'
        OPEN(IF_GRID,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%gridFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.grdb'
        OPEN(IF_GRID,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
#else
      IF (global%gridFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.grda_', &
                                   global%timeStamp
        OPEN(IF_GRID,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%gridFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.grdb_', &
                                   global%timeStamp
        OPEN(IF_GRID,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
#endif

! - steady flow

    ELSE
      IF (global%gridFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.grda'
        OPEN(IF_GRID,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%gridFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.grdb'
        OPEN(IF_GRID,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  ENDIF   ! MASTERPROC

! read & broadcast time in file -----------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_ReadDataFileReal( global,IF_GRID,global%gridFormat,1,1,rvar )
#ifndef GENX
    global%currentTime = rvar(1,1)
#endif
  ENDIF

#ifdef MPI
  CALL MPI_Bcast( global%currentTime,1,MPI_RFREAL,MASTERPROC, &
                  global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! read grid data --------------------------------------------------------------

  DO iReg=1,global%nRegions

! - read region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      CALL RFLO_ReadDataFileInt( global,IF_GRID,global%gridFormat,4,1,ivar )
      iRegFile = ivar(1,1)
      ipc      = ivar(2,1)
      jpc      = ivar(3,1)
      kpc      = ivar(4,1)

      IF (iRegFile /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__,'File: '//TRIM(fname) )
      IF ((ipc /= regions(iReg)%levels(1)%grid%ipc) .OR. &
          (jpc /= regions(iReg)%levels(1)%grid%jpc) .OR. &
          (kpc /= regions(iReg)%levels(1)%grid%kpc)) THEN
        WRITE(msg,1000) iReg,ipc,jpc,kpc
        CALL ErrorStop( global,ERR_GRID_DIMENSIONS,__LINE__,msg )
      ENDIF
    ENDIF

! - master reads & sends data, others receive them

    IF (global%myProcid == MASTERPROC) THEN

      nDim = (ipc+1)*(jpc+1)*(kpc+1)
      ALLOCATE( xyzFile(3,nDim),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL RFLO_ReadDataFileReal( global,IF_GRID,global%gridFormat,3,nDim, &
                                  xyzFile )
#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        CALL MPI_Send( xyzFile,3*nDim,MPI_RFREAL,regions(iReg)%procid,iReg, &
                       global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif

    ELSE   ! not the master

      IF (regions(iReg)%procid == global%myProcid) THEN
        nDim = (regions(iReg)%levels(1)%grid%ipc+1)* &
               (regions(iReg)%levels(1)%grid%jpc+1)* &
               (regions(iReg)%levels(1)%grid%kpc+1)
        ALLOCATE( xyzFile(3,nDim),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#ifdef MPI
        CALL MPI_Recv( xyzFile,3*nDim,MPI_RFREAL,MASTERPROC,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
      ENDIF

    ENDIF

! - copy grid into data structure

    IF (regions(iReg)%procid == global%myProcid) THEN
      CALL RFLO_GetDimensPhysNodes( regions(iReg),1,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
      CALL RFLO_GetNodeOffset( regions(iReg),1,iNOff,ijNOff )
      xyz => regions(iReg)%levels(1)%grid%xyz

      n = 0
      DO k=kpnbeg,kpnend
        DO j=jpnbeg,jpnend
          DO i=ipnbeg,ipnend
            n    = n + 1
            ijkN = IndIJK(i,j,k,iNOff,ijNOff)
            xyz(XCOORD,ijkN) = xyzFile(1,n)
            xyz(YCOORD,ijkN) = xyzFile(2,n)
            xyz(ZCOORD,ijkN) = xyzFile(3,n)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (ALLOCATED(xyzFile)) THEN
      DEALLOCATE( xyzFile,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
    ENDIF

  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(IF_GRID,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', ipc= ',I6,', jpc= ',I6,', kpc= ',I6,'.')

END SUBROUTINE RFLO_ReadGrid

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadGrid.F90,v $
! Revision 1.5  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/06/28 21:19:24  wasistho
! if GENX don't read currenttime
!
! Revision 1.2  2005/05/27 08:04:58  wasistho
! allow genx read initial grid
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.17  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.14  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.13  2002/10/23 18:43:24  jblazek
! Changed temporary pointer arrays into allocatable arrays
! in grid and solution I/O routines.
!
! Revision 1.12  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.11  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.10  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.9  2002/08/29 21:52:21  jblazek
! Added I/O of grid speeds.
!
! Revision 1.8  2002/06/07 16:40:36  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.7  2002/03/18 21:56:39  jblazek
! Finished multiblock and MPI.
!
! Revision 1.6  2002/02/21 23:25:04  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.4  2001/12/22 00:09:37  jblazek
! Added routines to store grid and solution.
!
! Revision 1.3  2001/12/19 23:09:20  jblazek
! Added routines to read grid and solution.
!
! Revision 1.2  2001/12/08 00:18:41  jblazek
! Added routines to read BC input file.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************







