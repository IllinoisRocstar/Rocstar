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
! Purpose: read in turbulence solution.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions = dimensions of all regions.
!
! Output: 
!        globalClass LES  : region%levels%mixt%tv = turbulent viscosity
!        globalClass RANS : no yet defined
!        global%currentTime    = physical time
!        global%resInit        = initial residual
!
! Notes: solution is read in only for the current grid level;
!        solution is also read in for all dummy cells.
!
!******************************************************************************
!
! $Id: TURB_rFLO_ReadSolution.F90,v 1.7 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_ReadSolution( regions ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k, n

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg, timeString
  LOGICAL                :: fileExists

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, nField, iOff, ijOff, ijk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
  INTEGER :: nDimC, nRvar, nOutSol, globalClass, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL), POINTER     :: tv(:,:), tcv(:,:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), solFile(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_RFLO_ReadSolution',&
  'TURB_rFLO_ReadSolution.F90' )

! allocate fixed-size temporary data arrays -----------------------------------

  nRvar = 4

  ALLOCATE( ivar(6,1),stat=errorFlag )
  ALLOCATE( rvar(nRvar,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! copy time to string ---------------------------------------------------------

  IF (global%flowType == FLOW_UNSTEADY) THEN
    WRITE(timeString,'(1PE11.5)') global%timeStamp
  ELSE
    WRITE(timeString,'(1PE11.5)') 0._RFREAL
  ENDIF

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

! - unsteady flow

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.turba_', &
                                   global%timeStamp
        INQUIRE(FILE=fname,EXIST=fileExists)
        IF (fileExists .EQV. .true.) &
          OPEN(IF_SOLUT,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.turbb_', &
                                   global%timeStamp
        INQUIRE(FILE=fname,EXIST=fileExists)
        IF (fileExists .EQV. .true.) &
          OPEN(IF_SOLUT,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF

! - steady flow

    ELSE
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.turba_', &
                                global%currentIter
        INQUIRE(FILE=fname,EXIST=fileExists)
        IF (fileExists .EQV. .true.) &
          OPEN(IF_SOLUT,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.turbb_', &
                                global%currentIter
        INQUIRE(FILE=fname,EXIST=fileExists)
        IF (fileExists .EQV. .true.) &
          OPEN(IF_SOLUT,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
    ENDIF

    IF (fileExists .EQV. .true.) THEN
      global%error = errorFlag
      IF (global%error /= 0) &
        CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

      IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE)  &
        WRITE(STDOUT,'(A)') SOLVER_NAME//'   - turbulence'
    ELSE
      IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE)  &
        WRITE(STDOUT,'(A)') SOLVER_NAME//' Start new turbulence ...'
    ENDIF

  ENDIF   ! MASTERPROC

#ifdef MPI
  CALL MPI_Bcast( fileExists,1,MPI_LOGICAL,MASTERPROC,global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

  IF (.NOT. (fileExists .eqv. .true.)) THEN
    GOTO 999
  ENDIF

! read & broadcast time and initial residual in file --------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_ReadDataFileReal( global,IF_SOLUT,global%solutFormat,nRvar,1,rvar )
  ENDIF

#ifdef MPI
  CALL MPI_Bcast( rvar,nRvar,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
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
  global%resInit      = rvar(2,1)
  global%esg1Sum      = rvar(3,1)
  global%esg4Sum      = rvar(4,1)

! read solution data ----------------------------------------------------------

! first define no.of output solution by subsquent selection (order matters)

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid == global%myProcid) THEN
      IF (regions(iReg)%turbInput%modelClass == MODEL_LES) THEN
        globalClass = MODEL_LES
      ENDIF
    ENDIF
  ENDDO

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid == global%myProcid) THEN
      IF ((regions(iReg)%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
        globalClass = MODEL_RANS
      ENDIF
    ENDIF
  ENDDO

  DO iReg=1,global%nRegions

! - get dimensions and pointers

    iLev = regions(iReg)%currLevel
    CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regions(iReg),iLev,iOff,ijOff )
    ijkBeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
    ijkEnd = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
    nDimC  = ijkEnd - ijkBeg + 1

! - read region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      CALL RFLO_ReadDataFileInt( global,IF_SOLUT,global%solutFormat,6,1,ivar )
      iRegFile  = ivar(1,1)
      ipc       = ivar(2,1)
      jpc       = ivar(3,1)
      kpc       = ivar(4,1)
      nDumCells = ivar(5,1)
      nField    = ivar(6,1)

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
    ENDIF

! - master reads & sends data, others receive them

    IF (global%myProcid == MASTERPROC) THEN

      nOutSol = nField
      ALLOCATE( solFile(nOutSol,nDimC),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL RFLO_ReadDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                  nOutSol,nDimC,solFile )

#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        CALL MPI_Send( nOutSol,1,MPI_INTEGER, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        CALL MPI_Send( solFile,nOutSol*nDimC,MPI_RFREAL, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif

    ELSE   ! not the master

      IF (regions(iReg)%procid == global%myProcid) THEN
#ifdef MPI
        CALL MPI_Recv( nOutSol,1,MPI_INTEGER,MASTERPROC,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
        ALLOCATE( solFile(nOutSol,nDimC),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

#ifdef MPI
        CALL MPI_Recv( solFile,nOutSol*nDimC,MPI_RFREAL,MASTERPROC,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
      ENDIF

    ENDIF

! - copy solution into data structure

    IF (regions(iReg)%procid == global%myProcid) THEN
      IF (regions(iReg)%turbInput%modelClass == MODEL_LES) THEN
        tv => regions(iReg)%levels(iLev)%mixt%tv
        n  = 0
        DO k=kdcbeg,kdcend
          DO j=jdcbeg,jdcend
            DO i=idcbeg,idcend
              n   = n + 1
              ijk = IndIJK(i,j,k,iOff,ijOff)
              tv(TV_MIXT_MUET,ijk) = solFile(1,n)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF ((regions(iReg)%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
        tcv => regions(iReg)%levels(iLev)%turb%cv
        n  = 0
        DO k=kdcbeg,kdcend
          DO j=jdcbeg,jdcend
            DO i=idcbeg,idcend
              n   = n + 1
              ijk = IndIJK(i,j,k,iOff,ijOff)
              tcv(CV_SA_NUTIL,ijk) = solFile(1,n)
            ENDDO
          ENDDO
        ENDDO
      ENDIF    ! turbModel/modelClass
    ENDIF      ! global%myProcid

    IF (ALLOCATED(solFile).eqv..true.) THEN
      DEALLOCATE( solFile,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
    ENDIF

  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(IF_SOLUT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

999 CONTINUE

  CALL DeregisterFunction( global )

1000 FORMAT('Time in file is= ',1PE12.5,' but it should be= ',E12.5,'.')
1005 FORMAT('Region ',I5,', ipc= ',I6,', jpc= ',I6,', kpc= ',I6,'.')
1010 FORMAT('Region ',I5,', # dummy cells=',I2,' but should be= ',I1)

END SUBROUTINE TURB_RFLO_ReadSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_ReadSolution.F90,v $
! Revision 1.7  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.6  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/02/07 06:04:18  wasistho
! bug fixed: broadcast fileExists before used
!
! Revision 1.3  2006/02/06 06:13:28  wasistho
! enable starting turbulence from laminar run
!
! Revision 1.2  2005/03/09 06:37:23  wasistho
! incorporated HDESSA
!
! Revision 1.1  2004/03/11 03:26:34  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.6  2004/02/26 21:28:07  wasistho
! added esg1Sum and esg4Sum to Real heading for restart
!
! Revision 1.5  2004/02/11 03:25:01  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.4  2003/10/09 22:50:02  wasistho
! mv call to TURB_RansSAGetEddyVis from readSolution to initSolution
!
! Revision 1.3  2003/10/07 02:06:41  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.2  2003/09/05 21:50:07  wasistho
! removed ifdef PERI
!
! Revision 1.1  2003/07/22 03:01:12  wasistho
! prepare more accurate rocturb restart
!
!
!******************************************************************************









