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
! Purpose: read in time averaged statistics of the mixture.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary
!
! Input: regions = dimensions of all regions
!
! Output: region%levels%mixt%tav = time avg mixture variables (current grid level)
!         region%levels%turb%tav = time avg TURB variables (current grid level)
!         global%integrTime      = integrated averaging time
!
! Notes: time averaged solution is read in only for the current grid level;
!        it is also read in for all dummy cells
!
!******************************************************************************
!
! $Id: RFLO_ReadStat.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadStat( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k, l, n, ind

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, nDim, iOff, ijOff, ijk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
  INTEGER :: mixtVarId(2,regions(1)%global%mixtNStat)
  INTEGER :: turbVarId(2,regions(1)%global%turbNStat)
  INTEGER :: nTavgVar, maxNStat, allNStat, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:), jvar(:,:)

  REAL(RFREAL), POINTER     :: mixttav(:,:), turbtav(:,:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), tavFile(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadStat',&
  'RFLO_ReadStat.F90' )

! allocate temporary data arrays --------------------------------------------

  ALLOCATE( ivar(5,1),stat=errorFlag )
  ALLOCATE( rvar(2,1),stat=errorFlag )
  maxNStat = global%mixtNStat
#ifdef TURB
  IF (global%turbNStat>0) maxNStat = MAX( global%mixtNStat,global%turbNStat )
#endif  
  ALLOCATE( jvar(maxNStat+1,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open statistics file (only master proc.) ----------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    IF (global%solutFormat == FORMAT_ASCII) THEN
      WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.stata_', &
                                 global%timeStamp
      OPEN(IF_STAT,file=fname,form='formatted',status='old',iostat=errorFlag)
    ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
      WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.statb_', &
                                 global%timeStamp
      OPEN(IF_STAT,file=fname,form='unformatted',status='old',iostat=errorFlag)
    ELSE 
      CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  ENDIF   ! MASTERPROC

! read & broadcast current and integrated time in file, and stats ID ---------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_ReadDataFileReal( global,IF_STAT,global%solutFormat,2,1,rvar )
  ENDIF

#ifdef MPI
  CALL MPI_Bcast( rvar,2,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

  IF (global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL) THEN
    IF (global%currentTime /= rvar(1,1)) THEN
      WRITE(msg,1000) rvar(1,1),global%currentTime
      CALL ErrorStop( global,ERR_TIME_SOLUTION,__LINE__,msg//' File: '//TRIM(fname) )
    ENDIF
  ELSE
    global%currentTime = rvar(1,1)
  ENDIF
  global%integrTime = rvar(2,1)

  IF (global%myProcid == MASTERPROC) THEN

! - mixture statistics NSTATS and ID

    IF (global%mixtNStat > 0) THEN
      CALL RFLO_ReadDataFileInt( global,IF_STAT,global%solutFormat, &
                                 global%mixtNStat+1,1,jvar )
      nTavgVar  = jvar(1,1)
      IF (nTavgVar /= global%mixtNStat) THEN
        CALL ErrorStop( global,ERR_STATS_RESTART,__LINE__ )
      ENDIF

      mixtVarId(1,:) = jvar(2:global%mixtNStat+1,1)
      mixtVarId(2,:) = mod(mixtVarId(1,:),10)
      mixtVarId(1,:) = (mixtVarId(1,:)-mixtVarId(2,:))/10

      DO ind=1,2
        DO l=1,global%mixtNStat
          IF (mixtVarId(ind,l) /= global%mixtStatId(ind,l)) &
          CALL ErrorStop( global,ERR_STATS_RESTART,__LINE__ )
        ENDDO 
      ENDDO
    ENDIF

! - turbulence statistics NSTATS and ID

#ifdef TURB
    IF (global%turbNStat > 0) THEN
      CALL RFLO_ReadDataFileInt( global,IF_STAT,global%solutFormat, &
                                 global%turbNStat+1,1,jvar )
      nTavgVar  = jvar(1,1)
      IF (nTavgVar /= global%turbNStat) THEN
        CALL ErrorStop( global,ERR_STATS_RESTART,__LINE__ )
      ENDIF

      turbVarId(1,:) = jvar(2:global%turbNStat+1,1)
      turbVarId(2,:) = mod(turbVarId(1,:),10)
      turbVarId(1,:) = (turbVarId(1,:)-turbVarId(2,:))/10

      DO ind=1,2
        DO l=1,global%turbNStat
          IF (turbVarId(ind,l) /= global%turbStatId(ind,l)) &
          CALL ErrorStop( global,ERR_STATS_RESTART,__LINE__ )
        ENDDO 
      ENDDO
    ENDIF
#endif
  ENDIF

! read statistics data from all regions ------------------------------------------

  allNStat = global%mixtNStat+global%turbNStat

  DO iReg=1,global%nRegions

! - get dimensions and pointers

    iLev = regions(iReg)%currLevel
    CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regions(iReg),iLev,iOff,ijOff )
    ijkBeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
    ijkEnd = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
    nDim   = ijkEnd - ijkBeg + 1

! - read region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      CALL RFLO_ReadDataFileInt( global,IF_STAT,global%solutFormat,5,1,ivar )
      iRegFile  = ivar(1,1)
      ipc       = ivar(2,1)
      jpc       = ivar(3,1)
      kpc       = ivar(4,1)
      nDumCells = ivar(5,1)
      
      IF (iRegFile /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__,'File: '//TRIM(fname) )
      IF ((ipc /= regions(iReg)%levels(iLev)%grid%ipc) .OR. &
          (jpc /= regions(iReg)%levels(iLev)%grid%jpc) .OR. &
          (kpc /= regions(iReg)%levels(iLev)%grid%kpc)) THEN
        WRITE(msg,1005) iReg,ipc,jpc,kpc
        CALL ErrorStop( global,ERR_GRID_DIMENSIONS,__LINE__,msg )
      ENDIF
      IF (nDumCells /= regions(iReg)%nDumCells) THEN
        WRITE(msg,1010) iReg,nDumCells,regions(iReg)%nDumCells
        CALL ErrorStop( global,ERR_GRID_DUMCELLS,__LINE__,msg )
      ENDIF

! - master reads & sends data, others receive them

      ALLOCATE( tavFile(allNStat,nDim),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL RFLO_ReadDataFileReal( global,IF_STAT,global%solutFormat, &
                                  allNStat,nDim,tavFile )
#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        CALL MPI_Send( tavFile,allNStat*nDim,MPI_RFREAL, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
    ELSE   ! not the master

      IF (regions(iReg)%procid == global%myProcid) THEN
        ALLOCATE( tavFile(allNStat,nDim),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#ifdef MPI
        CALL MPI_Recv( tavFile,allNStat*nDim,MPI_RFREAL,MASTERPROC, &
                       iReg,global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
      ENDIF

    ENDIF  ! global%myProcid

! - copy statistics into data structure

    IF (regions(iReg)%procid == global%myProcid) THEN
      IF (global%mixtNStat > 0) mixttav => regions(iReg)%levels(iLev)%mixt%tav
      n = 0
      DO k=kdcbeg,kdcend
        DO j=jdcbeg,jdcend
          DO i=idcbeg,idcend
            n   = n + 1
            ijk = IndIJK(i,j,k,iOff,ijOff)
            DO l=1,global%mixtNStat
              mixttav(l,ijk) = tavFile(l,n)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
#ifdef TURB
      IF (global%turbNStat > 0) THEN
        IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
            (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) THEN
          turbtav => regions(iReg)%levels(iLev)%turb%tav
          n = 0
          DO k=kdcbeg,kdcend
            DO j=jdcbeg,jdcend
              DO i=idcbeg,idcend
                n   = n + 1
                ijk = IndIJK(i,j,k,iOff,ijOff)
                DO l=1,global%turbNStat
                  turbtav(l,ijk) = tavFile(l+global%mixtNStat,n)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
#endif
    ENDIF

    IF (ALLOCATED(tavFile)) THEN
      DEALLOCATE( tavFile,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
    ENDIF

  ENDDO     ! iReg

! finalize -----------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(IF_STAT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

  CALL DeregisterFunction( global )

1000 FORMAT('Time in file is= ',1PE12.5,' but it should be= ',E12.5,'.')
1005 FORMAT('Region ',I5,', ipc= ',I6,', jpc= ',I6,', kpc= ',I6,'.')
1010 FORMAT('Region ',I5,', # dummy cells=',I2,' but should be= ',I1)

END SUBROUTINE RFLO_ReadStat

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadStat.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.14  2004/11/09 10:51:48  wasistho
! provide maximum size for jvar
!
! Revision 1.13  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.10  2003/08/11 20:12:51  wasistho
! put safety in reading turbulence statistics
!
! Revision 1.9  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.8  2002/11/04 18:41:57  wasistho
! Modified statistics restart
!
! Revision 1.7  2002/11/02 01:47:23  wasistho
! Added TURB statistics
!
! Revision 1.6  2002/10/23 18:43:24  jblazek
! Changed temporary pointer arrays into allocatable arrays
! in grid and solution I/O routines.
!
! Revision 1.5  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.4  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.3  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/07/22 15:44:08  wasistho
! Cleaned-up conforming Coding Rule
!
! Revision 1.1  2002/06/14 20:53:05  wasistho
! add time avg statistics
!
!******************************************************************************







