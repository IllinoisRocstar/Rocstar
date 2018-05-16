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
! Purpose: write time averaged solution to file for mixture
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary
!
! Input: regions              = dimensions and cons. variables of all regions
!        global%currentTime   = physical time
!        global%integrTime    = integrated time during time averaging process
!        mixtNStat, mixtStatId= number of mixture statistics variables and IDs 
!        mixttav              = time averaged mixture variables
!        turbNStat, turbStatId= number of TURB statistics variables and IDs 
!        turbtav              = time averaged TURB variables
!
! Output: to file
!
! Notes: solution is stored only for the current grid level; it is also
!        stored for all dummy cells; all regions are written into one file
!
!******************************************************************************
!
! $Id: RFLO_WriteStat.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_WriteStat( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_WriteDataFileInt, RFLO_WriteDataFileReal, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k, l, n

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, nDim, iOff, ijOff, ijk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
  INTEGER :: maxNStat, allNStat, errorFlag
  INTEGER :: mixtStatId(regions(1)%global%mixtNStat)
  INTEGER :: turbStatId(regions(1)%global%turbNStat)
  INTEGER, ALLOCATABLE :: ivar(:,:), jvar(:,:)

  REAL(RFREAL), POINTER     :: mixttav(:,:), turbtav(:,:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), tavFile(:,:)
  LOGICAL                   :: zeroTurb

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_WriteStat',&
  'RFLO_WriteStat.F90' )

! allocate temporary data arrays ---------------------------------------------

  ALLOCATE( ivar(5,1),stat=errorFlag )
  ALLOCATE( rvar(2,1),stat=errorFlag )
  maxNStat = MAX( global%mixtNStat,global%turbNStat )
  ALLOCATE( jvar(maxNStat+1,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open statistics file (only master proc.) -----------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    IF (global%solutFormat == FORMAT_ASCII) THEN
      WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.stata_', &
                                 global%currentTime
      OPEN(IF_STAT,file=fname,form='formatted',status='unknown', &
           iostat=errorFlag)
    ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
      WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.statb_', &
                                 global%currentTime
      OPEN(IF_STAT,file=fname,form='unformatted',status='unknown', &
           iostat=errorFlag)
    ELSE
      CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
    ENDIF
    rvar(1,1) = global%currentTime
    rvar(2,1) = global%integrTime

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! write current and integrated time to file ----------------------------------

    CALL RFLO_WriteDataFileReal( global,IF_STAT,global%solutFormat,2,1,rvar )

! mixture statistics NSTATS and ID

    IF (global%mixtNStat > 0) THEN
      jvar(1,1)    = global%mixtNStat
      mixtStatId(:)= global%mixtStatId(1,:)*10 + global%mixtStatId(2,:)
      jvar(2:global%mixtNStat+1,1) = mixtStatId(1:global%mixtNStat)
      CALL RFLO_WriteDataFileInt( global,IF_STAT,global%solutFormat, &
                                  global%mixtNStat+1,1,jvar )
    ENDIF

! turbulence statistics NSTATS and ID

#ifdef TURB
    IF (global%turbNStat > 0) THEN
      jvar(1,1)    = global%turbNStat
      turbStatId(:)= global%turbStatId(1,:)*10 + global%turbStatId(2,:)
      jvar(2:global%turbNStat+1,1) = turbStatId(1:global%turbNStat)
      CALL RFLO_WriteDataFileInt( global,IF_STAT,global%solutFormat, &
                                  global%turbNStat+1,1,jvar )
    ENDIF
#endif
  ENDIF

! write statistics data ------------------------------------------------------

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

! - allocate memory for data field

    IF (regions(iReg)%procid==global%myProcid .OR. &
        global%myProcid==MASTERPROC) THEN
      ALLOCATE( tavFile(allNStat,nDim),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ENDIF

! - copy statistics into data structure

    zeroTurb = .false.

    IF (regions(iReg)%procid == global%myProcid) THEN
      IF (global%mixtNStat > 0) mixttav => regions(iReg)%levels(iLev)%mixt%tav
#ifdef TURB
      IF (global%turbNStat > 0) THEN
        IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
            (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) THEN
          turbtav => regions(iReg)%levels(iLev)%turb%tav
        ELSE
          ALLOCATE( turbtav(global%turbNStat,nDim) )
          turbtav = 0._RFREAL
          zeroTurb = .true.
        ENDIF
      ENDIF
#endif
      n = 0
      DO k=kdcbeg,kdcend
        DO j=jdcbeg,jdcend
          DO i=idcbeg,idcend
            n   = n + 1
            ijk = IndIJK(i,j,k,iOff,ijOff)
            DO l=1,global%mixtNStat
              tavFile(l,n) = mixttav(l,ijk)
            ENDDO
#ifdef TURB
            DO l=1,global%turbNStat
              tavFile(global%mixtNStat+l,n) = turbtav(l,ijk)
            ENDDO
#endif
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (zeroTurb) DEALLOCATE( turbtav )

! - write region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      ivar(1,1) = iReg
      ivar(2,1) = regions(iReg)%levels(iLev)%grid%ipc
      ivar(3,1) = regions(iReg)%levels(iLev)%grid%jpc
      ivar(4,1) = regions(iReg)%levels(ilev)%grid%kpc
      ivar(5,1) = regions(iReg)%nDumCells
      CALL RFLO_WriteDataFileInt( global,IF_STAT,global%solutFormat,5,1,ivar )
    ENDIF

! - master receives and writes data, others send them

    IF (global%myProcid == MASTERPROC) THEN
#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        CALL MPI_Recv( tavFile,allNStat*nDim,MPI_RFREAL, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
      CALL RFLO_WriteDataFileReal( global,IF_STAT,global%solutFormat, &
                                   allNStat,nDim,tavFile )

    ELSE   ! not the master
#ifdef MPI
      IF (regions(iReg)%procid == global%myProcid) THEN
        CALL MPI_Send( tavFile,allNStat*nDim,MPI_RFREAL,MASTERPROC, &
                       iReg,global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
    ENDIF  ! global%myProcid

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

END SUBROUTINE RFLO_WriteStat

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_WriteStat.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.13  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.10  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.9  2002/12/12 03:39:37  wasistho
! facilitate the possibility of NO TURB statistics
!
! Revision 1.8  2002/11/04 18:42:08  wasistho
! Modified statistics restart
!
! Revision 1.7  2002/11/02 01:47:48  wasistho
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
! Revision 1.2  2002/07/22 15:44:24  wasistho
! Cleaned-up conforming Coding Rule
!
! Revision 1.1  2002/06/14 20:53:05  wasistho
! add time avg statistics
!
!******************************************************************************







