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
! Purpose: write turbulence instantaneous solution to file
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions            = dimensions and turb. variables of all regions
!        global%currentTime = physical time
!        global%resInit     = initial residual.
!
! Output: to file.
!
! Notes: solution is stored only for the current grid level; it is also
!        stored for all dummy cells. All regions are written into one file.
!
!******************************************************************************
!
! $Id: TURB_rFLO_WriteSolution.F90,v 1.5 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_WriteSolution( regions ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_WriteDataFileInt, RFLO_WriteDataFileReal, &
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
  TYPE(t_global), POINTER :: global

  CHARACTER(2*CHRLEN+17)  :: fname

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, iOff, ijOff, ijk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
  INTEGER :: nDimC, nRvar, nOutSol, globalClass, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL), POINTER     :: tv(:,:), tcv(:,:), vort(:,:), lens(:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), solFile(:,:)

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_RFLO_WriteSolution',&
  'TURB_rFLO_WriteSolution.F90' )

! allocate fixed-size temporary data arrays -----------------------------------

  nRvar = 4

  ALLOCATE( ivar(6,1),stat=errorFlag )
  ALLOCATE( rvar(nRvar,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

! - unsteady flow

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.turba_', &
                                   global%currentTime
        OPEN(IF_SOLUT,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.turbb_', &
                                   global%currentTime
        OPEN(IF_SOLUT,file=fname,form='unformatted',status='unknown', &
             iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
      rvar(1,1) = global%currentTime
      rvar(2,1) = 1._RFREAL

! - steady flow

    ELSE
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%outDir)//TRIM(global%casename)//'.turba_', &
                                global%currentIter
        OPEN(IF_SOLUT,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%outDir)//TRIM(global%casename)//'.turbb_', &
                                global%currentIter
        OPEN(IF_SOLUT,file=fname,form='unformatted',status='unknown', &
             iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
      rvar(1,1) = 0._RFREAL
      rvar(2,1) = global%resInit
    ENDIF
    rvar(3,1) = global%esg1Sum
    rvar(4,1) = global%esg4Sum

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  ENDIF   ! MASTERPROC

! write time and initial residual to file --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat,nRvar,1,rvar )
  ENDIF

! write solution data ----------------------------------------------------------

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

! - allocate memory for data field and initiate it to zero

    IF (regions(iReg)%procid==global%myProcid .OR. &
        global%myProcid==MASTERPROC) THEN
      nOutSol = regions(iReg)%turbInput%nOutField 
      ALLOCATE( solFile(nOutSol,nDimc),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
      solFile = 0._RFREAL
    ENDIF

! - copy solution(s), depending on model selected, into data structure

    IF (regions(iReg)%procid == global%myProcid) THEN
      IF (regions(iReg)%turbInput%modelClass == MODEL_LES) THEN
        tv => regions(iReg)%levels(iLev)%mixt%tv
        n = 0
        DO k=kdcbeg,kdcend
          DO j=jdcbeg,jdcend
            DO i=idcbeg,idcend
              n   = n + 1
              ijk = IndIJK(i,j,k,iOff,ijOff)
              solFile(1,n) = tv(TV_MIXT_MUET,ijk)
            ENDDO
          ENDDO
        ENDDO
        IF (nOutSol == 2) THEN
          vort => regions(iReg)%levels(iLev)%turb%vort
          n = 0
          DO k=kdcbeg,kdcend
            DO j=jdcbeg,jdcend
              DO i=idcbeg,idcend
                n   = n + 1
                ijk = IndIJK(i,j,k,iOff,ijOff)
                solFile(2,n) = SQRT( vort(XCOORD,ijk)*vort(XCOORD,ijk) + &
                                     vort(YCOORD,ijk)*vort(YCOORD,ijk) + &
                                     vort(ZCOORD,ijk)*vort(ZCOORD,ijk) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF  ! nOutSol
      ENDIF
      IF ((regions(iReg)%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
        tcv  => regions(iReg)%levels(iLev)%turb%cv
        n = 0
        DO k=kdcbeg,kdcend
          DO j=jdcbeg,jdcend
            DO i=idcbeg,idcend
              n   = n + 1
              ijk = IndIJK(i,j,k,iOff,ijOff)
              solFile(1,n) = tcv(CV_SA_NUTIL,ijk)
            ENDDO
          ENDDO
        ENDDO
        IF (nOutSol == 2) THEN
          vort => regions(iReg)%levels(iLev)%turb%vort
          n = 0
          DO k=kdcbeg,kdcend
            DO j=jdcbeg,jdcend
              DO i=idcbeg,idcend
                n   = n + 1
                ijk = IndIJK(i,j,k,iOff,ijOff)
                solFile(2,n) = SQRT( vort(XCOORD,ijk)*vort(XCOORD,ijk) + &
                                     vort(YCOORD,ijk)*vort(YCOORD,ijk) + &
                                     vort(ZCOORD,ijk)*vort(ZCOORD,ijk) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF  ! nOutSol
        IF (nOutSol == 3) THEN
          vort => regions(iReg)%levels(iLev)%turb%vort
          lens => regions(iReg)%levels(iLev)%turb%lens
          n = 0
          DO k=kdcbeg,kdcend
            DO j=jdcbeg,jdcend
              DO i=idcbeg,idcend
                n   = n + 1
                ijk = IndIJK(i,j,k,iOff,ijOff)
                solFile(2,n) = SQRT( vort(XCOORD,ijk)*vort(XCOORD,ijk) + &
                                     vort(YCOORD,ijk)*vort(YCOORD,ijk) + &
                                     vort(ZCOORD,ijk)*vort(ZCOORD,ijk) )
                solFile(3,n) = lens(ijk)
              ENDDO
            ENDDO
          ENDDO
        ENDIF  ! nOutSol
      ENDIF    ! turbModel/modelClass
    ENDIF      ! global%myProcid

! - write region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      ivar(1,1) = iReg
      ivar(2,1) = regions(iReg)%levels(iLev)%grid%ipc
      ivar(3,1) = regions(iReg)%levels(iLev)%grid%jpc
      ivar(4,1) = regions(iReg)%levels(ilev)%grid%kpc
      ivar(5,1) = regions(iReg)%nDumCells
      ivar(6,1) = nOutSol
      CALL RFLO_WriteDataFileInt( global,IF_SOLUT,global%solutFormat,6,1,ivar )
    ENDIF

! - master receives and writes data, others send them

    IF (global%myProcid == MASTERPROC) THEN
#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        CALL MPI_Recv( solFile,nOutSol*nDimC,MPI_RFREAL, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
      CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat,nOutSol, &
                                   nDimC,solFile )
    ELSE   ! not the master
#ifdef MPI
      IF (regions(iReg)%procid == global%myProcid) THEN
        CALL MPI_Send( solFile,nOutSol*nDimC,MPI_RFREAL, &
                       MASTERPROC,iReg,global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
    ENDIF

    IF (ALLOCATED(solFile) .eqv. .true.) THEN
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

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RFLO_WriteSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_WriteSolution.F90,v $
! Revision 1.5  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.4  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/03/09 06:37:50  wasistho
! incorporated HDESSA
!
! Revision 1.1  2004/03/11 03:26:34  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.5  2004/02/26 21:27:30  wasistho
! added esg1Sum and esg4Sum to Real heading for restart
!
! Revision 1.4  2004/02/11 03:24:50  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.3  2004/02/07 01:13:26  wasistho
! modified TURB_WriteSolution
!
! Revision 1.2  2003/10/07 02:08:12  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.1  2003/07/22 03:01:12  wasistho
! prepare more accurate rocturb restart
!
!
!******************************************************************************









