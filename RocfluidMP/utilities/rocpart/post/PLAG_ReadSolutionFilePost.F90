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
! Purpose: read in PLAG solution pertaining to the main variable filed.
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
!
! Notes: only unsteady solution file format is supported.
!
!******************************************************************************
!
! $Id: PLAG_ReadSolutionFilePost.F90,v 1.4 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadSolutionFilePost( regions )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag
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
  INTEGER :: i, iCont, iReg

! ... local variables
  CHARACTER(CHRLEN+17) :: fname
  CHARACTER(CHRLEN)    :: RCSIdentString, msg, timeString

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER, PARAMETER :: ONE = 1
  INTEGER :: bcType, errorFlag, iLev, iRegFile, n, n1, n2, nAiv, nArv, &
             nCont, nCv, nCvTile, nDimPlag, nDimPlagMax
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: aivFile, ivar
  INTEGER, POINTER,     DIMENSION(:)   :: pCvPlagMass
  INTEGER, POINTER,     DIMENSION(:,:) :: pAiv

  REAL(RFREAL), ALLOCATABLE, DIMENSION(:,:) :: arvFile, cvFile, dvFile, rvar
  REAL(RFREAL), POINTER,     DIMENSION(:,:) :: pArv, pCv

  TYPE(t_patch),      POINTER :: pPatch
  TYPE(t_plag),       POINTER :: pPlag
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReadSolutionFilePost.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_ReadSolutionFilePost',&
  'PLAG_ReadSolutionFilePost.F90' )

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
      WRITE(fname,'(A,1PE11.5)') &
       TRIM(global%inDir)//TRIM(global%casename)//'.plag_sola_',global%timeStamp
      OPEN(IF_SOLUT,file=fname,form='formatted',status='old',iostat=errorFlag)
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
  IF ( global%nProcAlloc > 1 ) THEN
    CALL MPI_Bcast( rvar,1,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
    IF (global%mpierr /= ERR_NONE) THEN
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
    END IF ! global%mpierr
  ENDIF ! nProcAlloc
#endif

  IF (global%currentTime>0._RFREAL) THEN
!    IF (global%currentTime /= rvar(1,1)) THEN
    IF (ABS(global%currentTime-rvar(1,1))/global%currentTime > 1.0E-03_RFREAL) THEN
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

! - read region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      CALL RFLO_ReadDataFileInt( global,IF_SOLUT,global%solutFormat,2,1,ivar )
      iRegFile  = ivar(1,1)
      nDimPlag  = ivar(2,1)

      IF (iRegFile /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__,'File: '//TRIM(fname) )
    ENDIF ! global%myProcid

! - broadcast current value of nDimPlag to other processors

#ifdef MPI
    IF ( global%nProcAlloc > 1 ) THEN
      CALL MPI_Bcast( nDimPlag,1,MPI_INTEGER,MASTERPROC,global%mpiComm,global%mpierr )
      IF (global%mpierr /= ERR_NONE) THEN
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      END IF ! global%mpierr
    ENDIF ! nProcAlloca
#endif

    regions(iReg)%levels(iLev)%plag%nPcls = nDimPlag

    nDimPlagMax = regions(iReg)%plagInput%nPclsMax

    PRINT*,'PLAG_readSolutionFilePost: iReg, nPcls,nPclsMax = ',iReg,nDimPlag,nDimPlagMax

    IF( nDimPlag > nDimPlagMax) THEN
      WRITE(*,*)'Increase Array size to fit Solution File'
      STOP
    ENDIF ! nDimPlag

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
          pPlag => regions(iReg)%levels(iLev)%plag
          pAiv  => pPlag%aiv

! -- copy nDimPlag to nPcls
          pPlag%nPcls = nDimPlag

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
         pPlag => regions(iReg)%levels(iLev)%plag
         pArv => pPlag%arv
         n=0
         DO i=1, nDimPlag
           n = n+1
           pArv(ARV_PLAG_SPLOAD,i) = arvFile(1,n)
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
         pPlag => regions(iReg)%levels(iLev)%plag
         pCv         => pPlag%cv
         pCvPlagMass => pPlag%cvPlagMass
         nCont = regions(iReg)%plagInput%nCont
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

999  CONTINUE
  CALL DeregisterFunction( global )

1000 FORMAT('Time in file is= ',1PE12.5,' but it should be= ',E12.5,'.')

END SUBROUTINE PLAG_ReadSolutionFilePost

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadSolutionFilePost.F90,v $
! Revision 1.4  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/03/06 23:27:44  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.1  2004/12/01 22:00:46  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2004/03/05 22:09:05  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/03/02 21:51:15  jferry
! Added output of vapor energy to rplagpost output file
!
! Revision 1.1.1.1  2003/05/06 16:14:38  fnajjar
! Import of postprocessing tool for Rocpart
!
!******************************************************************************







