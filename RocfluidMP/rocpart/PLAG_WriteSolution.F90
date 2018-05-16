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
! Purpose: write in PLAG solution.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions = dimensions of all regions.
!
! Output: region%levels%plag%cv        = conservative variables (current grid
!                                        level)
!         region%levels%plag%aiv       = auxilliary integer variables
!         region%levels%plag%arv       = auxilliary real variables
!         region%levels%tile%cv        = conservative variables
!
! Notes: only unsteady solution file format is supported.
!
!******************************************************************************
!
! $Id: PLAG_WriteSolution.F90,v 1.6 2009/10/26 00:19:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_WriteSolution( regions )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_WriteDataFileInt, RFLO_WriteDataFileReal
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: i, iCont, iPatch, iReg

! ... local variables
  CHARACTER(CHRLEN+17) :: fname
  CHARACTER(CHRLEN)    :: RCSIdentString, msg, timeString

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER, PARAMETER :: ONE = 1
  INTEGER :: bcType, errorFlag, iLev, iRegFile, n, n1, n2, nAiv, nArv, &
             nCont, nCv, nCvTile, nDimP, nDimT, nDvTile, nPcls, nIdNumberP
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: nDimPlag
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: nextIdNumber
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: aivFile, ivar
  INTEGER, POINTER,     DIMENSION(:)   :: pCvPlagMass, pCvTileMass
  INTEGER, POINTER,     DIMENSION(:,:) :: pAiv

  REAL(RFREAL), ALLOCATABLE, DIMENSION(:,:) :: arvFile, cvFile, dvFile, rvar
  REAL(RFREAL), POINTER,     DIMENSION(:,:) :: pArv, pCv, pCvTile, pDvTile

  TYPE(t_patch), POINTER     :: pPatch
  TYPE(t_plag),   POINTER    :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global), POINTER    :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_WriteSolution.F90,v $ $Revision: 1.6 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_WriteSolution',&
  'PLAG_WriteSolution.F90' )

  IF (.NOT. (global%plagUsed .eqv. .true.)) THEN
     IF (global%myProcid == MASTERPROC) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Warning: PLAG seems to be OFF.'
     ENDIF
     GOTO 999
  ENDIF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing PLAG solution file...'
  END IF ! global%verbLevel

! allocate fixed-size temporary data arrays -----------------------------------

  ALLOCATE( ivar(3,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'ivar' )
  END IF ! global%error

  ALLOCATE( rvar(1,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'rvar' )
  END IF ! global%error

  ALLOCATE( nDimPlag(global%nRegions),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'nDimPlag' )
  END IF ! global%error

  ALLOCATE( nextIdNumber(global%nRegions),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'nextIdNumber' )
  END IF ! global%error

  nDimPlag = 0
  nextIdNumber = 0

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    IF (global%solutFormat == FORMAT_ASCII) THEN
      WRITE(fname,'(A,1PE11.5)') &
       TRIM(global%outDir)//TRIM(global%casename)//'.plag_sola_',global%currentTime
      OPEN(IF_SOLUT,file=fname,form='formatted',status='unknown',iostat=errorFlag)
    ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
      WRITE(fname,'(A,1PE11.5)') &
       TRIM(global%outDir)//TRIM(global%casename)//'.plag_solb_', global%currentTime
      OPEN(IF_SOLUT,file=fname,form='unformatted',status='unknown',iostat=errorFlag)
    ENDIF ! global%solutFormat

    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )
    END IF ! global%error

  ENDIF   ! MASTERPROC
  rvar(1,1) = global%currentTime

! write current time to file -------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat,1,1,rvar )
  ENDIF

! communicate number of particles to master processor -------------------------
!  master receives and writes data, others send them

  DO iReg=1,global%nRegions

! - get dimensions

    iLev  = regions(iReg)%currLevel

    IF ( global%myProcid == MASTERPROC ) THEN
      nDimPlag(iReg) = regions(iReg)%levels(iLev)%plag%nPcls

#ifdef MPI
      IF ( regions(iReg)%procid /= MASTERPROC ) THEN
        CALL MPI_Recv( nDimPlag(iReg),1,MPI_INTEGER, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF ( global%mpierr /= ERR_NONE ) THEN
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        ENDIF ! global%mpierr
      ENDIF ! regions(iReg)%procid
#endif

    ELSE   ! not the master
      nDimP = regions(iReg)%levels(iLev)%plag%nPcls

#ifdef MPI
      IF ( regions(iReg)%procid == global%myProcid ) THEN
        CALL MPI_Send( nDimP,1,MPI_INTEGER,MASTERPROC,iReg, &
                       global%mpiComm,global%mpierr )
        IF ( global%mpierr /= ERR_NONE ) THEN
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        ENDIF ! global%mpierr
      ENDIF ! regions(iReg)%procid
#endif

    ENDIF ! global%myProcid

  END DO !iReg

! broadcast nDimPlag ----------------------------------------------------------

#ifdef MPI
  CALL MPI_Bcast( nDimPlag,global%nRegions,MPI_INTEGER,MASTERPROC, &
                  global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) THEN
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  END IF ! global%mpierr
#endif

! communicate next Id number to master processor ------------------------------
!  master receives and writes data, others send them

  DO iReg=1,global%nRegions

! - get dimensions

    iLev  = regions(iReg)%currLevel

    IF ( global%myProcid == MASTERPROC ) THEN
      nextIdNumber(iReg) = regions(iReg)%levels(iLev)%plag%nextIdNumber

#ifdef MPI
      IF ( regions(iReg)%procid /= MASTERPROC ) THEN
        CALL MPI_Recv( nextIdNumber(iReg),1,MPI_INTEGER, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF ( global%mpierr /= ERR_NONE ) THEN
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        ENDIF ! global%mpierr
      ENDIF ! regions(iReg)%procid
#endif

    ELSE   ! not the master
      nIdNumberP = regions(iReg)%levels(iLev)%plag%nextIdNumber

#ifdef MPI
      IF ( regions(iReg)%procid == global%myProcid ) THEN
        CALL MPI_Send( nIdNumberP,1,MPI_INTEGER,MASTERPROC,iReg, &
                       global%mpiComm,global%mpierr )
        IF ( global%mpierr /= ERR_NONE ) THEN
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        ENDIF ! global%mpierr
      ENDIF ! regions(iReg)%procid
#endif

    ENDIF ! global%myProcid

  END DO !iReg

! broadcast nextIdNumber ------------------------------------------------------

#ifdef MPI
  CALL MPI_Bcast( nextIdNumber,global%nRegions,MPI_INTEGER,MASTERPROC, &
                  global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) THEN
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  END IF ! global%mpierr
#endif

! write PLAG solution data ----------------------------------------------------

  DO iReg=1,global%nRegions

! - get dimensions ------------------------------------------------------------

    iLev  = regions(iReg)%currLevel
    nCont = regions(iReg)%plagInput%nCont
    nAiv  = AIV_PLAG_LAST
    nArv  = ARV_PLAG_LAST
    nCv   = CV_PLAG_LAST+nCont
    nDimP = nDimPlag(iReg)

! - set pointers --------------------------------------------------------------

    pPlag => regions(iReg)%levels(iLev)%plag
    pAiv  => pPlag%aiv
    pArv  => pPlag%arv
    pCv   => pPlag%cv
    pCvPlagMass => pPlag%cvPlagMass

! - write region number and dimensions (only master) --------------------------

    IF ( global%myProcid == MASTERPROC ) THEN
      ivar(1,1) = iReg
      ivar(2,1) = nDimPlag(iReg)
      ivar(3,1) = nextIdNumber(iReg)
      CALL RFLO_WriteDataFileInt( global,IF_SOLUT,global%solutFormat,3,1,ivar )
    ENDIF ! MASTERPROC

! - activate for number of particles greater than zero ------------------------

    SELECT CASE (nDimP)
      CASE (ONE:)

! -- allocate memory for data field : aivFile ---------------------------------

        IF ( regions(iReg)%procid == global%myProcid .OR. &
             global%myProcid==MASTERPROC                ) THEN
          ALLOCATE( aivFile(nAiv,nDimP),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'aivFile' )
          ENDIF ! global%error
        ENDIF

! -- copy solution into data structure ----------------------------------------

        IF ( regions(iReg)%procid == global%myProcid ) THEN
          n=0
          DO i=1, nDimP
            n = n+1
            aivFile(1,n) = pAiv(AIV_PLAG_PIDINI,i)
            aivFile(2,n) = pAiv(AIV_PLAG_REGINI,i)
            aivFile(3,n) = pAiv(AIV_PLAG_REGCRT,i)
            aivFile(4,n) = pAiv(AIV_PLAG_ICELLS,i)
            aivFile(5,n) = pAiv(AIV_PLAG_INDEXI,i)
            aivFile(6,n) = pAiv(AIV_PLAG_INDEXJ,i)
            aivFile(7,n) = pAiv(AIV_PLAG_INDEXK,i)
            aivFile(8,n) = pAiv(AIV_PLAG_BURNSTAT,i)
            aivFile(9,n) = pAiv(AIV_PLAG_STATUS,i)
          ENDDO ! i
        END IF !regions(iReg)%procid

! -- master receives and writes data, others send them ------------------------

        IF ( global%myProcid == MASTERPROC ) THEN
#ifdef MPI
          IF ( regions(iReg)%procid /= MASTERPROC ) THEN
            CALL MPI_Recv( aivFile,nAiv*nDimP,MPI_INTEGER, &
                           regions(iReg)%procid,iReg, &
                           global%mpiComm,status,global%mpierr )
            IF ( global%mpierr /= ERR_NONE ) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            ENDIF ! global%mpierr
          ENDIF ! regions(iReg)%procid
#endif
          CALL RFLO_WriteDataFileInt( global,IF_SOLUT,global%solutFormat, &
                                      nAiv,nDimP,aivFile )

        ELSE   ! not the master
#ifdef MPI
          IF (regions(iReg)%procid == global%myProcid) THEN
            CALL MPI_Send( aivFile,nAiv*nDimP,MPI_INTEGER,MASTERPROC,iReg, &
                           global%mpiComm,global%mpierr )
            IF ( global%mpierr /= ERR_NONE ) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            ENDIF ! global%mpierr
          ENDIF ! regions(iReg)%procid
#endif
        ENDIF ! global%myProcid

! -- deallocate memory for data field : aivFile -------------------------------

        IF ( ALLOCATED(aivFile) ) THEN
          DEALLOCATE( aivFile,stat=errorFlag )
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'aivFile' )
          ENDIF ! global%error
        ENDIF ! aivFile

! -- allocate memory for data field : arvFile ---------------------------------

        IF ( regions(iReg)%procid == global%myProcid .OR. &
             global%myProcid==MASTERPROC                ) THEN
          ALLOCATE( arvFile(nArv,nDimP),stat=errorFlag )
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'arvFile' )
          ENDIF ! global%error
        ENDIF  ! regions(iReg)%procid

! -- copy solution into data structure ----------------------------------------

        IF ( regions(iReg)%procid == global%myProcid ) THEN
          n=0
          DO i=1, nDimP
            n = n+1
            arvFile(1,n) = pArv(ARV_PLAG_SPLOAD,i)
            arvFile(2,n) = pArv(ARV_PLAG_DISTOT,i)
          ENDDO ! i
        ENDIF ! regions(iReg)%procid

! -- master receives and writes data, others send them ------------------------

        IF ( global%myProcid == MASTERPROC ) THEN
#ifdef MPI
          IF ( regions(iReg)%procid /= MASTERPROC ) THEN
            CALL MPI_Recv( arvFile,nArv*nDimP,MPI_RFREAL, &
                           regions(iReg)%procid,iReg, &
                           global%mpiComm,status,global%mpierr )
            IF ( global%mpierr /= ERR_NONE ) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            ENDIF ! global%mpierr
          ENDIF ! regions(iReg)%procid
#endif
          CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                       nArv,nDimP,arvFile )

        ELSE   ! not the master
#ifdef MPI
          IF ( regions(iReg)%procid == global%myProcid ) THEN
            CALL MPI_Send( arvFile,nArv*nDimP,MPI_RFREAL,MASTERPROC,iReg, &
                           global%mpiComm,global%mpierr )
            IF ( global%mpierr /= ERR_NONE ) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            ENDIF ! global%mpierr
          ENDIF ! regions(iReg)%procid
#endif
        ENDIF ! global%myProcid

! -- deallocate memory for data field : arvFile -------------------------------

        IF ( ALLOCATED(arvFile) ) THEN
          DEALLOCATE( arvFile,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'arvFile' )
          ENDIF ! global%error
        ENDIF ! arvFile

! -- allocate memory for data field : cvFile ----------------------------------

        IF ( regions(iReg)%procid == global%myProcid .OR. &
             global%myProcid==MASTERPROC                ) THEN
          ALLOCATE( cvFile(nCv,nDimP),stat=errorFlag )
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
          ENDIF ! global%error
        ENDIF ! regions(iReg)%procid

! -- copy solution into data structure ----------------------------------------

        IF ( regions(iReg)%procid == global%myProcid ) THEN
          n=0
          DO i=1, nDimP
           n = n+1
           cvFile(1,n) = pCv(CV_PLAG_XMOM,i)
           cvFile(2,n) = pCv(CV_PLAG_YMOM,i)
           cvFile(3,n) = pCv(CV_PLAG_ZMOM,i)
           cvFile(4,n) = pCv(CV_PLAG_ENER,i)
           cvFile(5,n) = pCv(CV_PLAG_XPOS,i)
           cvFile(6,n) = pCv(CV_PLAG_YPOS,i)
           cvFile(7,n) = pCv(CV_PLAG_ZPOS,i)
           cvFile(8,n) = pCv(CV_PLAG_ENERVAPOR,i)
           DO iCont = 1, nCont
            cvFile(CV_PLAG_LAST+iCont,n) = pCv(pCvPlagMass(iCont),i)
           ENDDO ! iCont
         ENDDO ! i
        ENDIF !regions(iReg)%procid

! -- master receives and writes data, others send them ------------------------

        IF ( global%myProcid == MASTERPROC ) THEN
#ifdef MPI
          IF ( regions(iReg)%procid /= MASTERPROC ) THEN
            CALL MPI_Recv( cvFile,nCv*nDimP,MPI_RFREAL, &
                           regions(iReg)%procid,iReg, &
                           global%mpiComm,status,global%mpierr )
            IF ( global%mpierr /= ERR_NONE ) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            ENDIF ! global%mpierr
          ENDIF ! regions(iReg)%procid
#endif
          CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                       nCv,nDimP,cvFile )

        ELSE   ! not the master
#ifdef MPI
          IF ( regions(iReg)%procid == global%myProcid ) THEN
            CALL MPI_Send( cvFile,nCv*nDimP,MPI_RFREAL,MASTERPROC,iReg, &
                           global%mpiComm,global%mpierr )
            IF ( global%mpierr /= ERR_NONE ) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            ENDIF ! global%mpierr
          ENDIF ! regions(iReg)%procid
#endif
        ENDIF ! global%myProcid

! -- deallocate memory for data field : cvFile --------------------------------

        IF ( ALLOCATED(cvFile) ) THEN
          DEALLOCATE( cvFile,stat=errorFlag )
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'cvFile' )
          ENDIF ! global%error
        ENDIF ! cvFile

    END SELECT ! nDimP
  END DO !iReg

! write tile solution data ----------------------------------------------------

  DO iReg=1,global%nRegions

! - get dimensions ------------------------------------------------------------

    iLev  = regions(iReg)%currLevel
    nCvTile   = CV_TILE_LAST+nCont
    nDvTile   = 3

! - loop over all patches -----------------------------------------------------

    DO iPatch=1,regions(iReg)%nPatches
      pPatch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType =  pPatch%bcType

      IF ( bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE ) THEN

! -- get tile dimensions ------------------------------------------------------

        n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1
        n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
        nDimT   = n1*n2

! -- set pointers -------------------------------------------------------------

        pTilePlag   => pPatch%tilePlag
        pCvTile     => pTilePlag%cv
        pCvTileMass => pTilePlag%cvTileMass
        pDvTile     => pTilePlag%dv

! -- write region number and dimensions (only master) -------------------------

        IF ( global%myProcid == MASTERPROC ) THEN
          ivar(1,1) = iReg
          ivar(2,1) = nDimT
          CALL RFLO_WriteDataFileInt( global,IF_SOLUT,global%solutFormat, &
                                      2,1,ivar )
        ENDIF ! global%myProcid

! -- activate for number of tiles greater than zero ---------------------------

        SELECT CASE (nDimT)
          CASE (ONE:)

! --- allocate memory for data field : cvFile ---------------------------------

            IF ( regions(iReg)%procid == global%myProcid .OR. &
                 global%myProcid==MASTERPROC                ) THEN
              ALLOCATE( cvFile(nCvTile,nDimT),stat=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN
                CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
              ENDIF ! global%error
            ENDIF ! regions(iReg)%procid

! --- copy solution into file data structure ----------------------------------

            IF ( regions(iReg)%procid == global%myProcid ) THEN
              n=0
              DO i=1, nDimT
                n = n+1
                cvFile(1,n) = pCvTile(CV_TILE_MOMNRM,i)
                cvFile(2,n) = pCvTile(CV_TILE_ENER,  i)
                DO iCont = 1, nCont
                  cvFile(2+iCont,n) = pCvTile(pCvTileMass(iCont),i)
                ENDDO ! iCont
              ENDDO ! i
            END IF !regions(iReg)%procid

! --- master receives and writes data, others send them -----------------------

            IF ( global%myProcid == MASTERPROC ) THEN
#ifdef MPI
              IF ( regions(iReg)%procid /= MASTERPROC ) THEN
                CALL MPI_Recv( cvFile,nCvTile*nDimT,MPI_RFREAL, &
                               regions(iReg)%procid,iReg, &
                               global%mpiComm,status,global%mpierr )
                IF ( global%mpierr /= ERR_NONE ) THEN
                  CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                ENDIF ! global%mpierr
              ENDIF ! regions(iReg)%procid
#endif
              CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                           nCvTile,nDimT,cvFile )

            ELSE   ! not the master
#ifdef MPI
              IF ( regions(iReg)%procid == global%myProcid ) THEN
                CALL MPI_Send( cvFile,nCvTile*nDimT,MPI_RFREAL,MASTERPROC,iReg, &
                               global%mpiComm,global%mpierr )
                IF ( global%mpierr /= ERR_NONE ) THEN
                  CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                ENDIF ! global%mpierr
              ENDIF ! regions(iReg)%procid
#endif
            ENDIF ! global%myProcid

! --- deallocate memory for data field : cvFile -------------------------------

            IF ( ALLOCATED(cvFile) ) THEN
              DEALLOCATE( cvFile,stat=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN
                CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'cvFile' )
              ENDIF ! global%error
            ENDIF ! cvFile

! --- allocate memory for data field : dvFile ---------------------------------

            IF ( regions(iReg)%procid == global%myProcid .OR. &
                 global%myProcid==MASTERPROC                ) THEN
              ALLOCATE( dvFile(nDvTile,nDimT),stat=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN
                CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
              ENDIF ! global%error
            ENDIF ! regions(iReg)%procid

! --- copy solution into file data structure ----------------------------------

            IF ( regions(iReg)%procid == global%myProcid ) THEN
              n=0
              DO i=1, nDimT
                n = n+1
                dvFile(1,n) = pDvTile(DV_TILE_COUNTDOWN,i)
                dvFile(2,n) = pDvTile(DV_TILE_DIAM    ,i)
                dvFile(3,n) = pDvTile(DV_TILE_SPLOAD  ,i)
              ENDDO ! i
            ENDIF !regions(iReg)%procid

! --- master receives and writes data, others send them -----------------------

            IF ( global%myProcid == MASTERPROC ) THEN
#ifdef MPI
              IF ( regions(iReg)%procid /= MASTERPROC ) THEN
                CALL MPI_Recv( dvFile,nDvTile*nDimT,MPI_RFREAL, &
                               regions(iReg)%procid,iReg, &
                               global%mpiComm,status,global%mpierr )
                IF ( global%mpierr /= ERR_NONE ) THEN
                  CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                ENDIF ! global%mpierr
              ENDIF ! regions(iReg)%procid
#endif
              CALL RFLO_WriteDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                           nDvTile,nDimT,dvFile )

            ELSE   ! not the master
#ifdef MPI
              IF ( regions(iReg)%procid == global%myProcid ) THEN
                CALL MPI_Send( dvFile,nDvTile*nDimT,MPI_RFREAL,MASTERPROC,iReg, &
                               global%mpiComm,global%mpierr )
                IF ( global%mpierr /= ERR_NONE ) THEN
                  CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                ENDIF ! global%mpierr
              ENDIF ! regions(iReg)%procid
#endif
            ENDIF ! global%myProcid

! --- deallocate memory for data field : dvFile -------------------------------

            IF ( ALLOCATED(dvFile) ) THEN
              DEALLOCATE( dvFile,stat=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN
                CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'dvFile' )
              ENDIF ! global%error
            ENDIF ! dvFile

         END SELECT ! nDimT
      ENDIF ! bcType
    ENDDO ! iPatch

  ENDDO ! iReg

! deallocate fixed-size temporary data arrays ---------------------------------

  DEALLOCATE( ivar,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'ivar' )
  END IF ! global%error

  DEALLOCATE( rvar,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'rvar' )
  END IF ! global%error

  DEALLOCATE( nDimPlag,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'nDimPlag' )
  END IF ! global%error

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(IF_SOLUT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

999  CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_WriteSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_WriteSolution.F90,v $
! Revision 1.6  2009/10/26 00:19:32  mtcampbe
! Updates for completion of NATIVE_MP_IO
!
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.2  2005/05/31 21:37:32  fnajjar
! Added ARV_PLAG_DISTOT for proper IO capabilities
!
! Revision 1.1  2004/12/01 20:58:22  fnajjar
! Initial revision after changing case
!
! Revision 1.13  2004/06/16 23:07:18  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.12  2004/04/09 23:15:45  fnajjar
! Added plag status to I/O
!
! Revision 1.11  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.10  2004/03/05 16:26:42  fnajjar
! Added dv(diam) and dv(spload) from tile datastructure to insure proper restart
!
! Revision 1.9  2004/02/14 21:29:01  fnajjar
! Bug fix for cvFile with incorrect index
!
! Revision 1.8  2004/02/13 23:22:07  fnajjar
! Included new cv and aiv definitions for particle burning module
!
! Revision 1.7  2003/11/21 22:43:18  fnajjar
! Removed nPclsTot and added nextIdNumber
!
! Revision 1.6  2003/05/14 00:41:21  fnajjar
! Moved pointer definitions outside IF statments
!
! Revision 1.5  2003/04/14 16:31:05  jferry
! added check that particles are used in some region
!
! Revision 1.4  2003/02/25 20:16:46  fnajjar
! Complete rewrite for proper IO capability
!
! Revision 1.3  2002/12/05 16:14:23  f-najjar
! Added dv for time factor in restart file
!
! Revision 1.2  2002/12/04 15:37:15  f-najjar
! Included restart capability for Rocpart
!
! Revision 1.1  2002/10/25 14:20:32  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







