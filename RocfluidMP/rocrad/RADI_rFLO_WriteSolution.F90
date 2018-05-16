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
! Purpose: write radiation solution, extinction coefficients and intensities, 
!          to file
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions            = dimensions and cons. variables of all regions
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
! $Id: RADI_rFLO_WriteSolution.F90,v 1.4 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_RFLO_WriteSolution( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_WriteDataFileInt, RFLO_WriteDataFileReal, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k, l, n

! ... local variables
  TYPE(t_global), POINTER :: global

  CHARACTER(2*CHRLEN+17)  :: fname

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells
  INTEGER :: iCOff, ijCOff, ijkC, iNOff, ijNOff, ijkN, ijkNi, ijkNj, ijkNk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
  INTEGER :: nAng, nDimC, nRvar, nSolComp, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL), POINTER     :: radInt(:,:), radCoef(:,:), qri(:),qrj(:),qrk(:)
  REAL(RFREAL), POINTER     :: rcv(:,:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), solFile(:,:)

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RADI_RFLO_WriteSolution',&
  'RADI_rFLO_WriteSolution.F90' )

! allocate fixed-size temporary data arrays -----------------------------------

  nRvar = 2

  ALLOCATE( ivar(5,1),stat=errorFlag )
  ALLOCATE( rvar(nRvar,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

! - unsteady flow

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.rada_', &
                                   global%currentTime
        OPEN(IF_SOLUT,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%outDir)//TRIM(global%casename)//'.radb_', &
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
        WRITE(fname,'(A,I6.6)') TRIM(global%outDir)//TRIM(global%casename)//'.rada_', &
                                global%currentIter
        OPEN(IF_SOLUT,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%outDir)//TRIM(global%casename)//'.radb_', &
                                global%currentIter
        OPEN(IF_SOLUT,file=fname,form='unformatted',status='unknown', &
             iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF
      rvar(1,1) = 0._RFREAL
      rvar(2,1) = global%resInit
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  ENDIF   ! MASTERPROC

! write time and initial residual to file -------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat,nRvar,1,rvar )
  ENDIF

! write solution data ---------------------------------------------------------

  DO iReg=1,global%nRegions

! - get dimensions and pointers

    iLev = regions(iReg)%currLevel
    CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
    CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
    ijkBeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
    ijkEnd = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
    nDimC  = ijkEnd - ijkBeg + 1

! - allocate memory for data field

    nAng     = regions(iReg)%radiInput%nAng
    nSolComp = 3+RADI_COEFF_NCOMP+nAng

    IF (regions(iReg)%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
      nSolComp = nSolComp+1
    ENDIF

    IF (regions(iReg)%procid==global%myProcid .OR. &
        global%myProcid==MASTERPROC) THEN
      ALLOCATE( solFile(nSolComp,nDimc),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ENDIF

! - copy solution into data structure

    IF (regions(iReg)%procid == global%myProcid) THEN
      IF (regions(iReg)%mixtInput%radiUsed) THEN
        qri     => regions(iReg)%levels(iLev)%radi%qri
        qrj     => regions(iReg)%levels(iLev)%radi%qrj
        qrk     => regions(iReg)%levels(iLev)%radi%qrk
        radInt  => regions(iReg)%levels(iLev)%radi%radInt
        radCoef => regions(iReg)%levels(iLev)%radi%radCoef
        n = 0
        IF (regions(iReg)%radiInput%radiModel /= RADI_MODEL_FLDTRAN) THEN
          DO k=kdcbeg,kdcend
            DO j=jdcbeg,jdcend
              DO i=idcbeg,idcend
                n    = n + 1
                ijkC = IndIJK(i,j,k,iCOff,ijCOff)
                ijkN = IndIJK(i,j,k,iNOff,ijNOff)
                ijkNi= IndIJK(i+1,j,k,iNOff,ijNOff)
                ijkNj= IndIJK(i,j+1,k,iNOff,ijNOff)
                ijkNk= IndIJK(i,j,k+1,iNOff,ijNOff)
                solFile(1,n) = 0.5_RFREAL*(qri(ijkN)+qri(ijkNi))
                solFile(2,n) = 0.5_RFREAL*(qrj(ijkN)+qrj(ijkNj))
                solFile(3,n) = 0.5_RFREAL*(qrk(ijkN)+qrk(ijkNk))
                solFile(4,n) = radCoef(ijkC,RADI_COEFF_EXTINCT)
                solFile(5,n) = radCoef(ijkC,RADI_COEFF_SCATTER)
                DO l = 1,nAng
                  solFile(5+l,n) = radInt(l,ijkC)
                ENDDO  ! l
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k
        ELSE
          rcv => regions(iReg)%levels(iLev)%radi%cv

          DO k=kdcbeg,kdcend
            DO j=jdcbeg,jdcend
              DO i=idcbeg,idcend
                n    = n + 1
                ijkC = IndIJK(i,j,k,iCOff,ijCOff)
                ijkN = IndIJK(i,j,k,iNOff,ijNOff)
                ijkNi= IndIJK(i+1,j,k,iNOff,ijNOff)
                ijkNj= IndIJK(i,j+1,k,iNOff,ijNOff)
                ijkNk= IndIJK(i,j,k+1,iNOff,ijNOff)
                solFile(1,n) = rcv(CV_RADI_ENER,ijkC)
                solFile(2,n) = 0.5_RFREAL*(qri(ijkN)+qri(ijkNi))
                solFile(3,n) = 0.5_RFREAL*(qrj(ijkN)+qrj(ijkNj))
                solFile(4,n) = 0.5_RFREAL*(qrk(ijkN)+qrk(ijkNk))
                solFile(5,n) = radCoef(ijkC,RADI_COEFF_EXTINCT)
                solFile(6,n) = radCoef(ijkC,RADI_COEFF_SCATTER)
                DO l = 1,nAng
                  solFile(6+l,n) = radInt(l,ijkC)
                ENDDO  ! l
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k
        ENDIF  ! radiModel
      ELSE  ! radiUsed
        solFile = 0._RFREAL
      ENDIF
    ENDIF      ! global%myProcid

! - write region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      ivar(1,1) = iReg
      ivar(2,1) = regions(iReg)%levels(iLev)%grid%ipc
      ivar(3,1) = regions(iReg)%levels(iLev)%grid%jpc
      ivar(4,1) = regions(iReg)%levels(ilev)%grid%kpc
      ivar(5,1) = regions(iReg)%nDumCells
      CALL RFLO_WriteDataFileInt( global,IF_SOLUT,global%solutFormat,5,1,ivar )
    ENDIF

! - master receives and writes data, others send them

    IF (global%myProcid == MASTERPROC) THEN
#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        CALL MPI_Recv( solFile, nSolComp*nDimC,MPI_RFREAL, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
      CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                   nSolComp,nDimC,solFile )
    ELSE   ! not the master
#ifdef MPI
      IF (regions(iReg)%procid == global%myProcid) THEN
        CALL MPI_Send( solFile,nSolComp*nDimC,MPI_RFREAL, &
                       MASTERPROC,iReg,global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
    ENDIF

    IF (ALLOCATED(solFile)) THEN
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

END SUBROUTINE RADI_RFLO_WriteSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_rFLO_WriteSolution.F90,v $
! Revision 1.4  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/23 03:51:12  wasistho
! changed RADI_WriteSol.. to RADI_RFLO_WriteSol..
!
! Revision 1.1  2004/09/22 02:35:50  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.4  2003/08/11 21:52:31  wasistho
! added cell avg. radiation fluxes to output solution
!
! Revision 1.3  2003/07/30 22:24:08  wasistho
! enter part and smoke data into radiation
!
! Revision 1.2  2003/07/22 03:05:58  wasistho
! include logical write-parameter
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!******************************************************************************









