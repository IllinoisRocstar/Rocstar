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
! Purpose: read in PLAG solution.
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
! $Id: PLAG_ReadSolution.F90,v 1.7 2009/10/26 00:19:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadSolution( regions )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal
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
             nCont, nCv, nCvTile, nDimPlag, nDimTile, nDvTile,nTiles,  &
             nextIdNumber
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: aivFile, ivar
  INTEGER, POINTER,     DIMENSION(:)   :: pCvPlagMass, pCvTileMass
  INTEGER, POINTER,     DIMENSION(:,:) :: pAiv

  REAL(RFREAL), ALLOCATABLE, DIMENSION(:,:) :: arvFile, cvFile, dvFile, rvar
  REAL(RFREAL), POINTER,     DIMENSION(:,:) :: pArv, pCv, pCvTile, pDvTile
  REAL(RFREAL) :: timediff

  TYPE(t_patch), POINTER     :: pPatch
  TYPE(t_plag),   POINTER    :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global), POINTER    :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReadSolution.F90,v $ $Revision: 1.7 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_ReadSolution',&
  'PLAG_ReadSolution.F90' )

  IF (.NOT. global%plagUsed) GOTO 999

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading PLAG solution file...'
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

! copy time to string ---------------------------------------------------------

  IF (global%flowType == FLOW_UNSTEADY) THEN
    WRITE(timeString,'(1PE11.5)') global%timeStamp
  ELSE
    WRITE(timeString,'(1PE11.5)') 0._RFREAL
  ENDIF

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    IF (global%solutFormat == FORMAT_ASCII) THEN
       IF( global%currentTime == 0.0_RFREAL) THEN
          WRITE(fname,'(A,1PE11.5)') &
               TRIM(global%inDir)//TRIM(global%casename)//'.plag_sola_',global%timeStamp
          OPEN(IF_SOLUT,file=fname,form='formatted',status='old',iostat=errorFlag)
       ELSE
          WRITE(fname,'(A,1PE11.5)') &
               TRIM(global%outDir)//TRIM(global%casename)//'.plag_sola_',global%timeStamp
          OPEN(IF_SOLUT,file=fname,form='formatted',status='old',iostat=errorFlag)
       ENDIF
    ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
      WRITE(fname,'(A,1PE11.5)') &
       TRIM(global%inDir)//TRIM(global%casename)//'.plag_solb_', global%timeStamp
      OPEN(IF_SOLUT,file=fname,form='unformatted',status='old',iostat=errorFlag)
    ENDIF ! global%solutFormat

    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )
    END IF ! global%error

  ENDIF   ! MASTERPROC

! read & check time stamp in file ----------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_ReadDataFileReal( global,IF_SOLUT,global%solutFormat,1,1,rvar )
  ENDIF

#ifdef MPI
  CALL MPI_Bcast( rvar,1,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
  IF (global%mpierr /= ERR_NONE) THEN
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  END IF ! global%mpierr
#endif

  IF (global%currentTime>0._RFREAL) THEN
    timediff = ABS(global%currentTime - rvar(1,1))
    IF (timediff.gt.1.0E-12) THEN
      WRITE(msg,1000) rvar(1,1),global%currentTime
      CALL ErrorStop( global,ERR_TIME_SOLUTION, __LINE__,msg//' File: '//TRIM(fname) )
    ENDIF
  ENDIF

! read solution data ----------------------------------------------------------

  DO iReg=1,global%nRegions

! - get dimensions and pointers for PLAG datastructure

    iLev = regions(iReg)%currLevel

    nCont = regions(iReg)%plagInput%nCont
    nAiv  = AIV_PLAG_LAST
    nArv  = ARV_PLAG_LAST
    nCv   = CV_PLAG_LAST+nCont

    pPlag => regions(iReg)%levels(iLev)%plag
    pAiv  => pPlag%aiv
    pArv  => pPlag%arv
    pCv   => pPlag%cv
    pCvPlagMass => pPlag%cvPlagMass

! - read region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      CALL RFLO_ReadDataFileInt( global,IF_SOLUT,global%solutFormat,3,1,ivar )
      iRegFile  = ivar(1,1)
      nDimPlag  = ivar(2,1)
      nextIdNumber  = ivar(3,1)

      IF (iRegFile /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__,'File: '//TRIM(fname) )
    ENDIF ! global%myProcid

! - broadcast current value of nDimPlag to other processors

#ifdef MPI
    CALL MPI_Bcast( nDimPlag,1,MPI_INTEGER,MASTERPROC,global%mpiComm,global%mpierr )
    IF (global%mpierr /= ERR_NONE) THEN
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
    END IF ! global%mpierr
#endif

! - broadcast current value of nextIdNumber to other processors

#ifdef MPI
    CALL MPI_Bcast( nextIdNumber,1,MPI_INTEGER,MASTERPROC,global%mpiComm,global%mpierr )
    IF (global%mpierr /= ERR_NONE) THEN
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
    END IF ! global%mpierr
#endif

! - activate for number of particles greater than zero

    SELECT CASE (nDimPlag)
      CASE (ONE:)

! -- master reads Aiv & sends data, others receive them

        IF (global%myProcid == MASTERPROC) THEN

          ALLOCATE( aivFile(nAiv,nDimPlag),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'aivFile' )
          END IF ! global%error

          CALL RFLO_ReadDataFileInt(  global,IF_SOLUT,global%solutFormat, &
                                      nAiv,nDimPlag,aivFile )

#ifdef MPI
          IF (regions(iReg)%procid /= MASTERPROC) THEN
            CALL MPI_Send( aivFile,nAiv*nDimPlag,MPI_INTEGER, &
                           regions(iReg)%procid,iReg, &
                           global%mpiComm,global%mpierr )
            IF (global%mpierr /= ERR_NONE ) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            END IF ! global%mpierr

          ENDIF ! regions(iReg)%procid
#endif

        ELSE ! not the master

          IF (regions(iReg)%procid == global%myProcid) THEN
            ALLOCATE( aivFile(nAiv,nDimPlag),stat=errorFlag )
            global%error = errorFlag
            IF (global%error /= ERR_NONE) THEN
              CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'aivFile' )
            ENDIF !global%error

#ifdef MPI
            CALL MPI_Recv( aivFile,nAiv*nDimPlag,MPI_INTEGER,MASTERPROC,iReg, &
                           global%mpiComm,status,global%mpierr )
            IF (global%mpierr /= ERR_NONE) THEN
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
            END IF ! global%mpierr
#endif
          END IF !regions(iReg)%procid

        END IF !global%myProcid

! -- copy solution into data structure

        IF (regions(iReg)%procid == global%myProcid) THEN

! -- copy nDimPlag to nPcls
          pPlag%nPcls = nDimPlag
          pPlag%nextIdNumber = nextIdNumber

          n=0
          DO i=1, nDimPlag
            n = n+1
            pAiv(AIV_PLAG_PIDINI,i) = aivFile(1,n)
            pAiv(AIV_PLAG_REGINI,i) = aivFile(2,n)
            pAiv(AIV_PLAG_REGCRT,i) = aivFile(3,n)
            pAiv(AIV_PLAG_ICELLS,i) = aivFile(4,n)
            pAiv(AIV_PLAG_INDEXI,i) = aivFile(5,n)
            pAiv(AIV_PLAG_INDEXJ,i) = aivFile(6,n)
            pAiv(AIV_PLAG_INDEXK,i) = aivFile(7,n)
            pAiv(AIV_PLAG_BURNSTAT,i) = aivFile(8,n)
            pAiv(AIV_PLAG_STATUS,i) = aivFile(9,n)
          ENDDO ! i
        END IF !regions(iReg)%procid

       IF (ALLOCATED(aivFile)) THEN
         DEALLOCATE( aivFile,stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'aivFile' )
         ENDIF ! global%error
       ENDIF ! aivFile

! -- master reads Arv & sends data, others receive them

       IF (global%myProcid == MASTERPROC) THEN
         ALLOCATE( arvFile(nArv,nDimPlag),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'aivFile' )
         END IF ! global%error
         CALL RFLO_ReadDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                     nArv,nDimPlag,arvFile )

#ifdef MPI
         IF (regions(iReg)%procid /= MASTERPROC) THEN
           CALL MPI_Send( arvFile,nArv*nDimPlag,MPI_RFREAL, &
                          regions(iReg)%procid,iReg, &
                          global%mpiComm,global%mpierr )
           IF (global%mpierr /= ERR_NONE ) THEN
             CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
           END IF ! global%mpierr

         ENDIF ! regions(iReg)%procid
#endif

       ELSE ! not the master

         IF (regions(iReg)%procid == global%myProcid) THEN
           ALLOCATE( arvFile(nArv,nDimPlag),stat=errorFlag )
           global%error = errorFlag
           IF (global%error /= ERR_NONE) THEN
             CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'arvFile' )
           ENDIF !global%error

#ifdef MPI
           CALL MPI_Recv( arvFile,nArv*nDimPlag,MPI_RFREAL,MASTERPROC,iReg, &
                          global%mpiComm,status,global%mpierr )
           IF (global%mpierr /= ERR_NONE) THEN
             CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
           END IF ! global%mpierr
#endif
         END IF !regions(iReg)%procid

       END IF !global%myProcid

! -- copy solution into data structure

       IF (regions(iReg)%procid == global%myProcid) THEN
         n=0
         DO i=1, nDimPlag
           n = n+1
           pArv(ARV_PLAG_SPLOAD,i) = arvFile(1,n)
           pArv(ARV_PLAG_DISTOT,i) = arvFile(2,n)
         ENDDO ! i
       END IF !regions(iReg)%procid

       IF (ALLOCATED(arvFile)) THEN
         DEALLOCATE( arvFile,stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'arvFile' )
         ENDIF ! global%error
       ENDIF ! arvFile

! -- master reads Cv & sends data, others receive them

       IF (global%myProcid == MASTERPROC) THEN
         ALLOCATE( cvFile(nCv,nDimPlag),stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
         END IF ! global%error

         CALL RFLO_ReadDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                     nCv,nDimPlag,cvFile )

#ifdef MPI
         IF (regions(iReg)%procid /= MASTERPROC) THEN
           CALL MPI_Send( cvFile,nCv*nDimPlag,MPI_RFREAL, &
                          regions(iReg)%procid,iReg, &
                          global%mpiComm,global%mpierr )
           IF (global%mpierr /= ERR_NONE ) THEN
             CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
           END IF ! global%mpierr

         ENDIF ! regions(iReg)%procid
#endif

       ELSE ! not the master

         IF (regions(iReg)%procid == global%myProcid) THEN
           ALLOCATE( cvFile(nCv,nDimPlag),stat=errorFlag )
           global%error = errorFlag
           IF (global%error /= ERR_NONE) THEN
             CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
           ENDIF !global%error

#ifdef MPI
           CALL MPI_Recv( cvFile,nCv*nDimPlag,MPI_RFREAL,MASTERPROC,iReg, &
                          global%mpiComm,status,global%mpierr )
           IF (global%mpierr /= ERR_NONE) THEN
             CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
           END IF ! global%mpierr
#endif
         END IF !regions(iReg)%procid

       END IF !global%myProcid

! -- copy solution into data structure

       IF (regions(iReg)%procid == global%myProcid) THEN
         n=0
         DO i=1, nDimPlag
           n = n+1
           pCv(CV_PLAG_XMOM,i) = cvFile(1,n)
           pCv(CV_PLAG_YMOM,i) = cvFile(2,n)
           pCv(CV_PLAG_ZMOM,i) = cvFile(3,n)
           pCv(CV_PLAG_ENER,i) = cvFile(4,n)
           pCv(CV_PLAG_XPOS,i) = cvFile(5,n)
           pCv(CV_PLAG_YPOS,i) = cvFile(6,n)
           pCv(CV_PLAG_ZPOS,i) = cvFile(7,n)
           pCv(CV_PLAG_ENERVAPOR,i) = cvFile(8,n)
           DO iCont = 1, nCont
            pCv(pCvPlagMass(iCont),i) = cvFile(CV_PLAG_LAST+iCont,n)
           ENDDO ! iCont
         ENDDO ! i
       END IF !regions(iReg)%procid

       IF (ALLOCATED(cvFile)) THEN
         DEALLOCATE( cvFile,stat=errorFlag )
         global%error = errorFlag
         IF (global%error /= ERR_NONE) THEN
           CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'cvFile' )
         ENDIF ! global%error
       ENDIF ! cvFile

    END SELECT ! nDimPlag
  ENDDO        ! iReg

! read tile datastructure -----------------------------------------------------

  DO iReg=1,global%nRegions

! - get dimensions ------------------------------------------------------------

    iLev    = regions(iReg)%currLevel
    nCont   = regions(iReg)%plagInput%nCont
    nCvTile = CV_TILE_LAST+nCont
    nDvTile = 3

    DO iPatch=1,regions(iReg)%nPatches

! - set pointers

      pPatch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType =  pPatch%bcType

      IF ( bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE ) THEN

! -- get total tile dimensions ------------------------------------------------

        n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1
        n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
        nTiles  = n1*n2

! -- set tile pointers --------------------------------------------------------

        pTilePlag   => pPatch%tilePlag
        pCvTile     => pTilePlag%cv
        pCvTileMass => pTilePlag%cvTileMass
        pDvTile     => pTilePlag%dv

! -- read dimensions (only master)

         IF (global%myProcid == MASTERPROC) THEN
           CALL RFLO_ReadDataFileInt( global,IF_SOLUT,global%solutFormat, &
                                      2,1,ivar )
           iRegFile  = ivar(1,1)
           nDimTile  = ivar(2,1)

           IF (iRegFile /= iReg) &
             CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__, &
                            'File: '//TRIM(fname) )

           IF (nDimTile > nTiles) &
             CALL ErrorStop( global,ERR_PLAG_TILESIZE,__LINE__, &
                            'File: '//TRIM(fname) )
         ENDIF ! global%myProcid

! - broadcast current value of nDimTile to other processors

#ifdef MPI
         CALL MPI_Bcast( nDimTile,1,MPI_INTEGER,MASTERPROC,global%mpiComm,global%mpierr )
         IF (global%mpierr /= ERR_NONE) THEN
           CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
         END IF ! global%mpierr
#endif

! -- activate for number of tiles greater than zero

         SELECT CASE (nDimTile)
           CASE (ONE:)

! --- master reads cv & sends data, others receive them

             IF (global%myProcid == MASTERPROC) THEN
               ALLOCATE( cvFile(nCvTile,nDimTile),stat=errorFlag )
               global%error = errorFlag
               IF (global%error /= ERR_NONE) THEN
                 CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
               END IF ! global%error

               CALL RFLO_ReadDataFileReal( global,IF_SOLUT, global%solutFormat, &
                                           nCvTile,nDimTile,cvFile )

#ifdef MPI
               IF (regions(iReg)%procid /= MASTERPROC) THEN
                 CALL MPI_Send( cvFile,nCvTile*nDimTile,MPI_RFREAL, &
                                regions(iReg)%procid,iReg, &
                                global%mpiComm,global%mpierr )
                 IF (global%mpierr /= ERR_NONE ) THEN
                   CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                 END IF ! global%mpierr

                ENDIF ! regions(iReg)%procid
#endif

             ELSE ! not the master

               IF (regions(iReg)%procid == global%myProcid) THEN
                 ALLOCATE( cvFile(nCvTile,nDimTile),stat=errorFlag )
                 global%error = errorFlag
                 IF (global%error /= ERR_NONE) THEN
                   CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
                 ENDIF !global%error

#ifdef MPI
                 CALL MPI_Recv( cvFile,nCvTile*nDimTile,MPI_RFREAL,&
                                MASTERPROC,iReg, &
                                global%mpiComm,status,global%mpierr )
                 IF (global%mpierr /= ERR_NONE) THEN
                   CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                 END IF ! global%mpierr
#endif
                END IF !regions(iReg)%procid

             END IF !global%myProcid

! --- copy solution into data structure

             IF (regions(iReg)%procid == global%myProcid) THEN
               n=0
               DO i=1, nDimTile
                 n = n+1
                 pCvTile(CV_TILE_MOMNRM,i) = cvFile(1,n)
                 pCvTile(CV_TILE_ENER,  i) = cvFile(2,n)
                 DO iCont = 1, nCont
                   pCvTile(pCvTileMass(iCont),i) = cvFile(2+iCont,n)
                 ENDDO ! iCont
               ENDDO ! i
             END IF !regions(iReg)%procid

             IF (ALLOCATED(cvFile)) THEN
               DEALLOCATE( cvFile,stat=errorFlag )
               global%error = errorFlag
               IF (global%error /= ERR_NONE) THEN
                 CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'cvFile' )
               ENDIF ! global%error
             ENDIF ! cvFile

! --- master reads dv & sends data, others receive them

             IF (global%myProcid == MASTERPROC) THEN
               ALLOCATE( dvFile(nDvTile,nDimTile),stat=errorFlag )
               global%error = errorFlag
               IF (global%error /= ERR_NONE) THEN
                 CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
               END IF ! global%error

               CALL RFLO_ReadDataFileReal( global,IF_SOLUT, global%solutFormat, &
                                           nDvTile,nDimTile,dvFile )

#ifdef MPI
               IF (regions(iReg)%procid /= MASTERPROC) THEN
                 CALL MPI_Send( dvFile,nDvTile*nDimTile,MPI_RFREAL, &
                                regions(iReg)%procid,iReg, &
                                global%mpiComm,global%mpierr )
                 IF (global%mpierr /= ERR_NONE ) THEN
                   CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                 END IF ! global%mpierr

                ENDIF ! regions(iReg)%procid
#endif

             ELSE ! not the master

               IF (regions(iReg)%procid == global%myProcid) THEN
                 ALLOCATE( dvFile(nDvTile,nDimTile),stat=errorFlag )
                 global%error = errorFlag
                 IF (global%error /= ERR_NONE) THEN
                   CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'cvFile' )
                 ENDIF !global%error

#ifdef MPI
                 CALL MPI_Recv( dvFile,nDvTile*nDimTile,MPI_RFREAL,&
                                MASTERPROC,iReg, &
                                global%mpiComm,status,global%mpierr )
                 IF (global%mpierr /= ERR_NONE) THEN
                   CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
                 END IF ! global%mpierr
#endif
                END IF !regions(iReg)%procid

             END IF !global%myProcid

! --- copy solution into data structure

             IF (regions(iReg)%procid == global%myProcid) THEN
               n=0
               DO i=1, nDimTile
                 n = n+1
                 pDvTile(DV_TILE_COUNTDOWN,i) = dvFile(1,n)
                 pDvTile(DV_TILE_DIAM    ,i) = dvFile(2,n)
                 pDvTile(DV_TILE_SPLOAD  ,i) = dvFile(3,n)
               ENDDO ! i
             END IF !regions(iReg)%procid

             IF (ALLOCATED(dvFile)) THEN
               DEALLOCATE( dvFile,stat=errorFlag )
               global%error = errorFlag
               IF (global%error /= ERR_NONE) THEN
                 CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'dvFile' )
               ENDIF ! global%error
             ENDIF ! dvFile

         END SELECT ! nDimTile

      END IF ! bcType
    ENDDO     ! iPatch
  ENDDO        ! iReg

! invoke MPI_barrier to insure all processors are insync ----------------------

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF (global%mpierr /= ERR_NONE ) THEN
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  END IF ! global%mpierr
#endif

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

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(IF_SOLUT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading PLAG solution file done...'
  END IF ! global%verbLevel

999  CONTINUE
  CALL DeregisterFunction( global )

1000 FORMAT('Time in file is= ',1PE12.5,' but it should be= ',E12.5,'.')

END SUBROUTINE PLAG_ReadSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadSolution.F90,v $
! Revision 1.7  2009/10/26 00:19:32  mtcampbe
! Updates for completion of NATIVE_MP_IO
!
! Revision 1.6  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.3  2005/05/31 21:37:32  fnajjar
! Added ARV_PLAG_DISTOT for proper IO capabilities
!
! Revision 1.2  2005/02/01 16:46:51  fnajjar
! Added IO informin that solution reading operation is complete
!
! Revision 1.1  2004/12/01 20:58:08  fnajjar
! Initial revision after changing case
!
! Revision 1.17  2004/06/16 23:07:17  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.16  2004/04/09 23:15:45  fnajjar
! Added plag status to I/O
!
! Revision 1.15  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.14  2004/03/05 16:26:42  fnajjar
! Added dv(diam) and dv(spload) from tile datastructure to insure proper restart
!
! Revision 1.13  2004/02/13 23:22:07  fnajjar
! Included new cv and aiv definitions for particle burning module
!
! Revision 1.12  2003/11/21 22:43:18  fnajjar
! Removed nPclsTot and added nextIdNumber
!
! Revision 1.11  2003/05/14 00:41:21  fnajjar
! Moved pointer definitions outside IF statments
!
! Revision 1.10  2003/04/09 15:01:29  jferry
! added check that particles are used in some region
!
! Revision 1.9  2003/02/25 23:29:29  fnajjar
! Included nTiles check
!
! Revision 1.8  2003/02/25 22:49:45  fnajjar
! Bug fix for nCont in reading Tile data
!
! Revision 1.7  2003/02/25 22:21:45  fnajjar
! Deallocate temporary arrays
!
! Revision 1.6  2003/02/25 22:06:57  fnajjar
! Bug fix for inconsistent pointers
!
! Revision 1.5  2003/02/04 19:05:27  f-najjar
! Added ifdef call around MPI_Barrier
!
! Revision 1.4  2003/01/24 16:59:35  f-najjar
! Invoke MPI_Bcast for nDimPlag and nDimTile and place MPI_Barrier at routine end
!
! Revision 1.3  2002/12/05 16:14:23  f-najjar
! Added dv for time factor in restart file
!
! Revision 1.2  2002/12/04 15:37:15  f-najjar
! Included restart capability for Rocpart
!
! Revision 1.1  2002/10/25 14:19:16  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







