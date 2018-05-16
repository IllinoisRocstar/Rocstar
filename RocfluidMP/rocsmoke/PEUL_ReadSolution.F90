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
! Purpose: read in flow solution (for Eulerian particles).
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions = dimensions of all regions.
!
! Output: region%levels%peul%cv = conservative variables (current grid level)
!         global%currentTime    = physical time
!         global%peulResInit    = initial residual.

! Notes: solution is read in only for the current grid level;
!        solution is also read in for all dummy cells.
!
!******************************************************************************
!
! $Id: PEUL_ReadSolution.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReadSolution( regions ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset, &
                            RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE ModInterfacesEulerian, ONLY: PEUL_EnforcePositivity
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPtype, i, j, k, n

! ... local variables
  CHARACTER(2*CHRLEN+22) :: fname,fhead
  CHARACTER(CHRLEN)      :: RCSIdentString, msg

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, iOff, ijOff, ijk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff
  INTEGER :: nDimC, nDimN, nCv, errorFlag
  INTEGER, POINTER :: ivar(:,:)

  REAL(RFREAL), POINTER :: rvar(:,:), cv(:,:), cvFile(:,:)

  LOGICAL :: fileExists

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_ReadSolution.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_ReadSolution',&
  'PEUL_ReadSolution.F90' )

! begin -----------------------------------------------------------------------

  IF (.NOT. global%peulUsed) GOTO 9

! allocate fixed-size temporary data arrays -----------------------------------

  ALLOCATE( ivar(5,1),stat=errorFlag )
  ALLOCATE( rvar(2,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    fhead = TRIM(global%inDir)//TRIM(global%casename)

    SELECT CASE(global%solutFormat)
    CASE (FORMAT_ASCII)
      fhead = TRIM(fhead)//'.peul_sola_'
    CASE (FORMAT_BINARY)
      fhead = TRIM(fhead)//'.peul_solb_'
    CASE DEFAULT
      CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
    END SELECT

    IF (global%flowType == FLOW_UNSTEADY) THEN
      WRITE(fname,'(A,ES11.5)') TRIM(fhead), global%timeStamp
    ELSE
      WRITE(fname,'(A,I6.6)')   TRIM(fhead), global%currentIter
    ENDIF

    INQUIRE(file=fname,exist=fileExists)

    IF (fileExists) THEN

      SELECT CASE(global%solutFormat)
      CASE (FORMAT_ASCII)
        OPEN(IF_PEUL_SOLUT,file=fname,form=  'formatted',status='old', &
            iostat=errorFlag)
      CASE (FORMAT_BINARY)
        OPEN(IF_PEUL_SOLUT,file=fname,form='unformatted',status='old', &
            iostat=errorFlag)
      END SELECT
      global%error = errorFlag
      IF (global%error /= 0) &
        CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

    ENDIF

  ENDIF   ! MASTERPROC

! read & broadcast time and initial residual in file --------------------------

  IF (global%myProcid == MASTERPROC) THEN
    IF (fileExists) THEN
      CALL RFLO_ReadDataFileReal( global,IF_PEUL_SOLUT,global%solutFormat,&
                                  2,1,rvar )
    ELSE
      rvar(1,1) = 0._RFREAL
      rvar(2,1) = 1._RFREAL
    ENDIF
  ENDIF

#ifdef MPI
  CALL MPI_Bcast( rvar,2,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

  IF (global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL) THEN
    IF (global%currentTime /= rvar(1,1) .AND. &
        .NOT.regions(1)%peulInput%constInit) THEN
      WRITE(msg,1000) rvar(1,1),global%currentTime
      CALL ErrorStop( global,ERR_TIME_SOLUTION,__LINE__,msg//' File: '// &
        TRIM(fname) )
    ENDIF
  ENDIF
  global%peulResInit = rvar(2,1)

! read solution data ----------------------------------------------------------

  DO iReg=1,global%nRegions

! - get dimensions and pointers

    iLev = regions(iReg)%currLevel
    CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regions(iReg),iLev,iOff,ijOff )
    ijkBeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
    ijkEnd = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
    nDimC  = ijkEnd - ijkBeg + 1

    CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                  jpnbeg,jpnend,kpnbeg,kpnend )
    CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
    nDimN = (regions(iReg)%levels(iLev)%grid%ipc+1)* &
            (regions(iReg)%levels(iLev)%grid%jpc+1)* &
            (regions(iReg)%levels(iLev)%grid%kpc+1)

    nCv = regions(iReg)%levels(iLev)%peul%nCv

! - read region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      IF (fileExists) THEN
        CALL RFLO_ReadDataFileInt( global,IF_PEUL_SOLUT,global%solutFormat, &
                                   5,1,ivar )
        iRegFile  = ivar(1,1)
        ipc       = ivar(2,1)
        jpc       = ivar(3,1)
        kpc       = ivar(4,1)
        nDumCells = ivar(5,1)
      ELSE
        iRegFile  = iReg
        ipc       = regions(iReg)%levels(iLev)%grid%ipc
        jpc       = regions(iReg)%levels(iLev)%grid%jpc
        kpc       = regions(iReg)%levels(iLev)%grid%kpc
        nDumCells = regions(iReg)%nDumCells
      ENDIF

      IF (iRegFile /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__, &
          'File: '//TRIM(fname) )
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
    ENDIF

! - master reads & sends data, others receive them

    IF (global%myProcid == MASTERPROC) THEN

      ALLOCATE( cvFile(nCv,nDimC),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      IF (fileExists) THEN
        CALL RFLO_ReadDataFileReal( global,IF_PEUL_SOLUT,global%solutFormat, &
                                    nCv,nDimC,cvFile )
      ELSE
        DO iPtype = 1,nCv
          cvFile(iPtype,:) = regions(iReg)%peulInput%ptypes(iPtype)%initc
        END DO ! iPtype
      ENDIF

#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        CALL MPI_Send( cvFile,nCv*nDimC,MPI_RFREAL, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,global%mpierr )
        IF (global%mpierr /=0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

      ENDIF
#endif

    ELSE    ! not the master, but the cringing, obsequious slave

      IF (regions(iReg)%procid == global%myProcid) THEN
        ALLOCATE( cvFile(nCv,nDimC),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

#ifdef MPI
        CALL MPI_Recv( cvFile,nCv*nDimC,MPI_RFREAL,MASTERPROC,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

#endif
      ENDIF

    ENDIF   ! MASTERPROC

! - copy solution into data structure

    IF (regions(iReg)%procid == global%myProcid) THEN
      cv => regions(iReg)%levels(iLev)%peul%cv
      n  = 0
      DO k=kdcbeg,kdcend
        DO j=jdcbeg,jdcend
          DO i=idcbeg,idcend
            n   = n + 1
            ijk = IndIJK(i,j,k,iOff,ijOff)
            cv(:,ijk) = cvFile(:,n)
          ENDDO
        ENDDO
      ENDDO
      CALL PEUL_EnforcePositivity(regions(iReg))
    ENDIF      ! global%myProcid

    IF (ASSOCIATED(cvFile)) THEN
      DEALLOCATE( cvFile,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
      NULLIFY( cvFile )
    ENDIF

  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    IF (fileExists) CLOSE(IF_PEUL_SOLUT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

9 CONTINUE
  CALL DeregisterFunction( global )

1000 FORMAT('Time in file is= ',ES12.5,' but it should be= ',ES12.5,'.')
1005 FORMAT('Region ',I5,', ipc= ',I6,', jpc= ',I6,', kpc= ',I6,'.')
1010 FORMAT('Region ',I5,', # dummy cells=',I2,' but should be= ',I1)

END SUBROUTINE PEUL_ReadSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReadSolution.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:48  haselbac
! Initial revision after changing case
!
! Revision 1.10  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.9  2004/03/02 21:44:52  jferry
! Added clipping options
!
! Revision 1.8  2003/09/26 22:52:12  jferry
! changed file number for read/write of rocsmoke solutions
!
! Revision 1.7  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.6  2003/04/14 16:33:52  jferry
! added option to initialize to constant for t > 0
!
! Revision 1.5  2003/04/09 16:05:18  jferry
! changed naming convention for solution file
!
! Revision 1.4  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.3  2003/03/11 16:04:57  jferry
! Created data type for material properties
!
! Revision 1.2  2003/02/12 23:34:48  jferry
! Replaced [io]stat=global%error with local errorFlag
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







