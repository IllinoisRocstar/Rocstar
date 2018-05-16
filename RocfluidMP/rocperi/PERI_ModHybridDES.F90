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
! ******************************************************************************
!
! Purpose: Suite of rocperi routines associated with hybrid DES.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PERI_ModHybridDES.F90,v 1.9 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE PERI_ModHybridDES

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModParameters
  USE PERI_ModParameters
#ifdef TURB
  USE TURB_ModParameters
#endif
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
#ifdef RFLO
  PUBLIC :: PERI_RFLO_ReadMean, &
            PERI_CoMeanCorrection
#endif
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PERI_ModHybridDES.F90,v $ $ $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

#ifdef RFLO
  CONTAINS

!******************************************************************************
!
! Purpose: read in mean flow solution
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions = dimensions of all regions.
!
! Output: region%levels%peri%cvMean    = mean conservative variables (current 
!                                        grid level)
!
! Notes: solution and grid speeds are read in only for the current grid level;
!        solution is also read in for all dummy cells.
!
!******************************************************************************

SUBROUTINE PERI_RFLO_ReadMean( regions )

  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal
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
  INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, iOff, ijOff, ijk
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend
  INTEGER :: nDimC, nRvar, intvar(8), errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL) :: rnik
  REAL(RFREAL), POINTER     :: cvMean(:,:)
  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), cvFile(:,:)
  LOGICAL :: doread

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'PERI_RFLO_ReadMean',&
  'PERI_ModHybridDES.F90' )


! check if reading mean flow is applicable ------------------------------------

  doread = .FALSE.
#ifdef TURB
  DO iReg=1,global%nRegions
    IF (regions(iReg)%periInput%flowKind /= PERI_FLOW_NONE) doread = .TRUE.
    IF ((doread .EQV. .TRUE.) .AND. &
        (regions(iReg)%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
      EXIT
    ELSE
      doread = .FALSE.
    ENDIF
  ENDDO
#endif
  IF (doread .EQV. .FALSE.) GOTO 999 

! allocate fixed-size temporary data arrays -----------------------------------

  nRvar = 3
  ALLOCATE( ivar(5,1),stat=errorFlag )
  ALLOCATE( rvar(nRvar,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    IF (global%solutFormat == FORMAT_ASCII) THEN
      WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.msola'
      OPEN(IF_SOLUT,file=fname,form='formatted',status='old',iostat=errorFlag)
    ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
      WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.msolb'
      OPEN(IF_SOLUT,file=fname,form='unformatted',status='old',iostat=errorFlag)
    ELSE
      CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  ENDIF   ! MASTERPROC

! read & broadcast time and initial residual in file --------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_ReadDataFileReal( global,IF_SOLUT,global%solutFormat,nRvar,1,rvar )
  ENDIF

! read solution data ----------------------------------------------------------

  DO iReg=1,global%nRegions
    iLev = regions(iReg)%currLevel

! - read region number and dimensions (only master)

    IF (global%myProcid == MASTERPROC) THEN
      CALL RFLO_ReadDataFileInt( global,IF_SOLUT,global%solutFormat,5,1,ivar )
      iRegFile  = ivar(1,1)
      ipc       = ivar(2,1)
      jpc       = ivar(3,1)
      kpc       = ivar(4,1)
      nDumCells = ivar(5,1)

! --- get dimensions and pointers

      idcbeg = 1-nDumCells
      idcend = ipc+nDumCells
      jdcbeg = 1-nDumCells
      jdcend = jpc+nDumCells
      kdcbeg = 1-nDumCells
      kdcend = kpc+nDumCells
      nDimC  = (idcend-idcbeg+1)*(jdcend-jdcbeg+1)*(kdcend-kdcbeg+1)
      
      IF (iRegFile /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__,'File: '//TRIM(fname) )
      IF ((ipc /= regions(iReg)%levels(iLev)%grid%ipc) .OR. &
          (jpc /= regions(iReg)%levels(iLev)%grid%jpc)) THEN
        WRITE(msg,1005) iReg,ipc,jpc
        CALL ErrorStop( global,ERR_GRID_DIMENSIONS,__LINE__,msg )
      ENDIF
      IF (nDumCells /= regions(iReg)%nDumCells) THEN
        WRITE(msg,1010) iReg,nDumCells,regions(iReg)%nDumCells
        CALL ErrorStop( global,ERR_GRID_DUMCELLS,__LINE__,msg )
      ENDIF
    ENDIF

! - master reads & sends data, others receive them

    IF (global%myProcid == MASTERPROC) THEN

      ALLOCATE( cvFile(CV_MIXT_NEQS,nDimC),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL RFLO_ReadDataFileReal( global,IF_SOLUT,global%solutFormat, &
                                  CV_MIXT_NEQS,nDimC,cvFile )

#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN
        intvar(1) = idcbeg
        intvar(2) = idcend
        intvar(3) = jdcbeg
        intvar(4) = jdcend
        intvar(5) = kdcbeg
        intvar(6) = kdcend
        intvar(7) = nDimC
        intvar(8) = nDumCells
        CALL MPI_Send( intvar,8,MPI_INTEGER,regions(iReg)%procid, &
                       global%nRegions+iReg,global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

        CALL MPI_Send( cvFile,CV_MIXT_NEQS*nDimC,MPI_RFREAL, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif

    ELSE   ! not the master

      IF (regions(iReg)%procid == global%myProcid) THEN
#ifdef MPI
        CALL MPI_Recv( intvar,8,MPI_INTEGER,MASTERPROC,global%nRegions+iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        idcbeg    = intvar(1)
        idcend    = intvar(2)
        jdcbeg    = intvar(3)
        jdcend    = intvar(4)
        kdcbeg    = intvar(5)
        kdcend    = intvar(6)
        nDimC     = intvar(7)
        nDumCells = intvar(8)
#endif
        ALLOCATE( cvFile(CV_MIXT_NEQS,nDimC),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#ifdef MPI
        CALL MPI_Recv( cvFile,CV_MIXT_NEQS*nDimC,MPI_RFREAL,MASTERPROC,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
      ENDIF

    ENDIF

! - copy solution into data structure

    IF (regions(iReg)%procid == global%myProcid) THEN
      cvMean => regions(iReg)%levels(iLev)%peri%cvMean
      cvMean =  0._RFREAL
      iOff   =  idcend-idcbeg+1
      ijOff  =  iOff*(jdcend-jdcbeg+1)
      n      =  0
      rnik   =  1._RFREAL/REAL((idcend-idcbeg+1-2*ndumCells)* &
                               (kdcend-kdcbeg+1-2*ndumCells))
      DO k=kdcbeg+ndumCells,kdcend-ndumCells
        DO j=jdcbeg,jdcend
          DO i=idcbeg+ndumCells,idcend-ndumCells
            n   = n + 1
            ijk = IndIJK(i,j,k,iOff,ijOff)- &
                  IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff) + 1
            cvMean(j,CV_MIXT_DENS) = cvMean(j,CV_MIXT_DENS)+cvFile(1,ijk)
            cvMean(j,CV_MIXT_XMOM) = cvMean(j,CV_MIXT_XMOM)+cvFile(2,ijk)
            cvMean(j,CV_MIXT_YMOM) = cvMean(j,CV_MIXT_YMOM)+cvFile(3,ijk)
            cvMean(j,CV_MIXT_ZMOM) = cvMean(j,CV_MIXT_ZMOM)+cvFile(4,ijk)
            cvMean(j,CV_MIXT_ENER) = cvMean(j,CV_MIXT_ENER)+cvFile(5,ijk)
          ENDDO
        ENDDO
      ENDDO
      DO j=jdcbeg,jdcend
        cvMean(j,CV_MIXT_DENS) = cvMean(j,CV_MIXT_DENS)*rnik
        cvMean(j,CV_MIXT_XMOM) = cvMean(j,CV_MIXT_XMOM)*rnik
        cvMean(j,CV_MIXT_YMOM) = cvMean(j,CV_MIXT_YMOM)*rnik
        cvMean(j,CV_MIXT_ZMOM) = cvMean(j,CV_MIXT_ZMOM)*rnik
        cvMean(j,CV_MIXT_ENER) = cvMean(j,CV_MIXT_ENER)*rnik
      ENDDO
    ENDIF      ! global%myProcid

    IF (ALLOCATED(cvFile)) THEN
      DEALLOCATE( cvFile,stat=errorFlag )
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

999  CONTINUE
  CALL DeregisterFunction( global )

1005 FORMAT('Region ',I5,', ipc= ',I6,', jpc= ',I6,', kpc= ',I6,'.')
1010 FORMAT('Region ',I5,', # dummy cells=',I2,' but should be= ',I1)

END SUBROUTINE PERI_RFLO_ReadMean

!******************************************************************************
!
! Purpose: correct mean flow by RaNS mean flow
!
! Description: none.
!
! Input: region = data current region
!
! Output: region%levels%mixt%cv = conservative variables (current grid level)
!
! Notes: tav indexes for rho,ru,rv,rw,pr,rw change if input mixtStatId in 
!        .input file change.
!
!******************************************************************************

SUBROUTINE PERI_CoMeanCorrection( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset, &
                            MixtPerf_G_CpR, MixtPerf_R_M
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, iOff, ijOff, iNOff, ijNOff, ijk, ijkN, ijkNp1
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ibc, iec
  INTEGER :: ipcbeg, jpcbeg, kpcbeg, ipcend, jpcend, kpcend
  INTEGER :: indCp, indMol, errorFlag
  INTEGER, ALLOCATABLE :: jid(:)

  REAL(RFREAL) :: mm,rgas,cpgas,g, rrho,rho,ru,rv,rw,pr,re, rAvgTim,beta
  REAL(RFREAL) :: dy1,vistau1,yn1, muel,muet,uv, upl,umi,dudy,ynd,ombeta
  REAL(RFREAL) :: restau,modtau,vistau,tottau, u,v,w, ampli, alpha,ralpha
  REAL(RFREAL) :: tavMassflux, cvmMassflux, rnik, fact
  REAL(RFREAL) :: sndMassFlux, rcvMassFlux, sndAverSize, rcvAverSize
  REAL(RFREAL), POINTER :: xyz(:,:),cv(:,:),dv(:,:),gv(:,:),tv(:,:)
  REAL(RFREAL), POINTER :: cvMean(:,:),tav(:,:),ttav(:,:)
  REAL(RFREAL), ALLOCATABLE :: rhom(:),um(:),vm(:),wm(:),rem(:),dy(:)
  REAL(RFREAL), ALLOCATABLE :: rhoa(:),ua(:),va(:),wa(:),rea(:)
  REAL(RFREAL), ALLOCATABLE :: muela(:),mueta(:),uva(:),dstress(:),yc(:),yn(:)
  LOGICAL :: docorrect
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'PERI_CoMeanCorrection',&
  'PERI_ModHybridDES.F90' )

! check if mean-flow correction is applicable ---------------------------------

  docorrect = .FALSE.
#ifdef TURB
  IF (region%periInput%flowKind /= PERI_FLOW_NONE) docorrect = .TRUE.
  IF ((docorrect .EQV. .TRUE.) .AND. &
      (region%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
  ELSE
    docorrect = .FALSE.
  ENDIF
#endif
#ifdef STATS
  IF ((docorrect .EQV. .TRUE.) .AND. &
      (global%doStat == ACTIVE)) THEN
  ELSE
    docorrect = .FALSE.
  ENDIF
#else
  docorrect = .FALSE.
#endif
  IF (docorrect .EQV. .FALSE.) GOTO 999 

! correct mean flow -----------------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iOff,ijOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz    => region%levels(iLev)%grid%xyz
  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv
  gv     => region%levels(iLev)%mixt%gv
  tv     => region%levels(iLev)%mixt%tv
  cvMean => region%levels(iLev)%peri%cvMean

#ifdef STATS
! compute mean profiles from tav or space-averaged

  ALLOCATE( rhom(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( um(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( vm(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( wm(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( rem(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( rhoa(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( ua(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( va(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( wa(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( rea(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( dy(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( yc(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( yn(jdcbeg:jdcend+1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( muela(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( mueta(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( uva(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ALLOCATE( dstress(jdcbeg:jdcend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( jid(jdcbeg-1:jdcend+1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  tav    => region%levels(iLev)%mixt%tav
  ttav   => region%levels(iLev)%turb%tav

  indCp    = region%levels(iLev)%mixt%indCp 
  indMol   = region%levels(iLev)%mixt%indMol 
  beta     = 1._RFREAL
  ombeta   = 1._RFREAL-beta
  alpha    = 1._RFREAL
  ralpha   = 1._RFREAL/alpha
  rnik     = 1._RFREAL/REAL((ipcend-ipcbeg+1)*(kpcend-kpcbeg+1))

  DO j=jdcbeg,jdcend

    rhom(j) = 0._RFREAL
    um(j)   = 0._RFREAL
    vm(j)   = 0._RFREAL
    wm(j)   = 0._RFREAL
    rem(j)  = 0._RFREAL

    DO k=kpcbeg,kpcend
      DO i=ipcbeg,ipcend
        ijk   = IndIJK(i,j,k,iOff,ijOff)
        mm    = gv(GV_MIXT_MOL ,ijk*indMol)
        rgas  = MixtPerf_R_M( mm ) 
        cpgas = gv(GV_MIXT_CP  ,ijk*indCp )
        g     = MixtPerf_G_CpR( cpgas,rgas ) 

        rho = cv(CV_MIXT_DENS,ijk)
        ru  = cv(CV_MIXT_XMOM,ijk)
        rv  = cv(CV_MIXT_YMOM,ijk)
        rw  = cv(CV_MIXT_ZMOM,ijk)
        re  = cv(CV_MIXT_ENER,ijk)

        rrho = 1._RFREAL/rho

        rhom(j) = rhom(j) + rho
        um(j)   =  um(j)  + ru*rrho
        vm(j)   =  vm(j)  + rv*rrho
        wm(j)   =  wm(j)  + rw*rrho
        rem(j)  = rem(j)  + re

      ENDDO
    ENDDO
  ENDDO

  IF (global%integrTime < 1.e-13_RFREAL) THEN
    rAvgTim = 0._RFREAL

    DO j=jdcbeg,jdcend

      rhoa(j) = 0._RFREAL
      ua(j)   = 0._RFREAL
      va(j)   = 0._RFREAL
      wa(j)   = 0._RFREAL
      rea(j)  = 0._RFREAL

      muela(j) = 0._RFREAL
      mueta(j) = 0._RFREAL
      uva(j)   = 0._RFREAL

      DO k=kpcbeg,kpcend
        DO i=ipcbeg,ipcend
          ijk   = IndIJK(i,j,k,iOff,ijOff)
          mm    = gv(GV_MIXT_MOL ,ijk*indMol)
          rgas  = MixtPerf_R_M( mm ) 
          cpgas = gv(GV_MIXT_CP  ,ijk*indCp )
          g     = MixtPerf_G_CpR( cpgas,rgas ) 

!          rho = tav(1,ijk)*rAvgTim
!          ru  = rho*tav(2,ijk)*rAvgTim
!          rv  = rho*tav(3,ijk)*rAvgTim
!          rw  = 0._RFREAL
!          pr  = rho*rgas*tav(4,ijk)*rAvgTim
!          re  = pr/(g-1._RFREAL)+0.5_RFREAL*(ru*ru+rv*rv+rw*rw)/rho

          rho = cv(CV_MIXT_DENS,ijk)
          ru  = cv(CV_MIXT_XMOM,ijk)
          rv  = cv(CV_MIXT_YMOM,ijk)
          rw  = cv(CV_MIXT_ZMOM,ijk)
          re  = cv(CV_MIXT_ENER,ijk)

          rrho = 1._RFREAL/rho

          rhoa(j) = rhoa(j) + rho
          ua(j)   =  ua(j)  + ru*rrho
          va(j)   =  va(j)  + rv*rrho
          wa(j)   =  wa(j)  + rw*rrho
          rea(j)  = rea(j)  + re

! ------- for stresses
!          muel = global%refVisc
!          muet = ttav(1,ijk)*rAvgTim
!          uv   = tav(9,ijk)*rAvgTim

          muel = tv(TV_MIXT_MUEL,ijk) 
          muet = tv(TV_MIXT_MUET,ijk) 
          uv   = cv(CV_MIXT_XMOM,ijk)*cv(CV_MIXT_YMOM,ijk)*rrho**2

          muela(j) = muela(j) + muel
          mueta(j) = mueta(j) + muet
          uva(j)   = uva(j)   + uv

        ENDDO
      ENDDO
    ENDDO
  ELSE   ! global%integrTime
    rAvgTim = 1._RFREAL/global%integrTime

    DO j=jdcbeg,jdcend

      rhoa(j) = 0._RFREAL
      ua(j)   = 0._RFREAL
      va(j)   = 0._RFREAL
      wa(j)   = 0._RFREAL
      rea(j)  = 0._RFREAL

      muela(j) = 0._RFREAL
      mueta(j) = 0._RFREAL
      uva(j)   = 0._RFREAL

      DO k=kpcbeg,kpcend
        DO i=ipcbeg,ipcend
          ijk   = IndIJK(i,j,k,iOff,ijOff)
          mm    = gv(GV_MIXT_MOL ,ijk*indMol)
          rgas  = MixtPerf_R_M( mm ) 
          cpgas = gv(GV_MIXT_CP  ,ijk*indCp )
          g     = MixtPerf_G_CpR( cpgas,rgas ) 

          rho = tav(1,ijk)*rAvgTim
          ru  = rho*tav(2,ijk)*rAvgTim
          rv  = rho*tav(3,ijk)*rAvgTim
          rw  = 0._RFREAL
          pr  = rho*rgas*tav(4,ijk)*rAvgTim
          re  = pr/(g-1._RFREAL)+0.5_RFREAL*(ru*ru+rv*rv+rw*rw)/rho

!          rho = cv(CV_MIXT_DENS,ijk)
!          ru  = cv(CV_MIXT_XMOM,ijk)
!          rv  = cv(CV_MIXT_YMOM,ijk)
!          rw  = cv(CV_MIXT_ZMOM,ijk)
!          re  = cv(CV_MIXT_ENER,ijk)

          rrho = 1._RFREAL/rho

          rhoa(j) = rhoa(j) + rho
          ua(j)   =  ua(j)  + ru*rrho
          va(j)   =  va(j)  + rv*rrho
          wa(j)   =  wa(j)  + rw*rrho
          rea(j)  = rea(j)  + re

! ------- for stresses
          muel = global%refVisc
          muet = ttav(1,ijk)*rAvgTim
          uv   = tav(9,ijk)*rAvgTim

!          muel = tv(TV_MIXT_MUEL,ijk) 
!          muet = tv(TV_MIXT_MUET,ijk) 
!          uv   = cv(CV_MIXT_XMOM,ijk)*cv(CV_MIXT_YMOM,ijk)*rrho**2

          muela(j) = muela(j) + muel
          mueta(j) = mueta(j) + muet
          uva(j)   = uva(j)   + uv

        ENDDO
      ENDDO
    ENDDO
  ENDIF  ! global%integrTime

! space averaged

  DO j=jdcbeg,jdcend
    rhom(j) = rhom(j)*rnik
    um(j)   = um(j)*rnik
    vm(j)   = vm(j)*rnik
    wm(j)   = wm(j)*rnik
    rem(j)  = rem(j)*rnik

    rhoa(j) = rhoa(j)*rnik
    ua(j)   = ua(j)*rnik
    va(j)   = va(j)*rnik
    wa(j)   = wa(j)*rnik
    rea(j)  = rea(j)*rnik

    muela(j) = muela(j)*rnik
    mueta(j) = mueta(j)*rnik
    uva(j)   = uva(j)*rnik
  ENDDO

! prepare tools for massflux ratio and stress-difference computation

  tavMassflux = 0._RFREAL
  cvmMassflux = 0._RFREAL
  DO k=1,1
    DO j=jdcbeg,jdcend
      DO i=1,1
        ijkN  = IndIJK(i ,j  ,k ,iNOff,ijNOff)
        ijkNp1= IndIJK(i ,j+1,k ,iNOff,ijNOff)
        dy(j) = (xyz(YCOORD,ijkNp1)-xyz(YCOORD,ijkN))
        yc(j) = 0.5_RFREAL*(xyz(YCOORD,ijkNp1)+xyz(YCOORD,ijkN))
        yn(j)   = xyz(YCOORD,ijkN)
        yn(j+1) = xyz(YCOORD,ijkNp1)
      ENDDO
    ENDDO
  ENDDO
  IF (region%iRegionGlobal == 1) THEN
    dy1     = dy(1)
    vistau1 = muela(1)*2._RFREAL*ua(1)/dy1
    yn1     = yn(1)
  ENDIF

#ifdef MPI
  CALL MPI_Bcast( dy1,1,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
  CALL MPI_Bcast( vistau1,1,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
  CALL MPI_Bcast( yn1,1,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
#endif

  DO j=jpcbeg,jpcend
    tavMassflux = tavMassflux + dy(j)*rhom(j)*um(j)
    cvmMassflux = cvmMassflux + dy(j)*cvMean(j,CV_MIXT_XMOM)
  ENDDO

#ifdef MPI
  sndMassFlux = tavMassflux
  CALL MPI_ALLREDUCE( sndMassFlux,rcvMassFlux,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      global%mpiComm, global%mpierr )
  tavMassflux = rcvMassFlux

  sndMassFlux = cvmMassflux
  CALL MPI_ALLREDUCE( sndMassFlux,rcvMassFlux,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      global%mpiComm, global%mpierr )
  cvmMassflux = rcvMassFlux

  IF (region%periInput%split(JCOORD) == OFF) THEN
    sndAverSize = rnik
    CALL MPI_ALLREDUCE( sndAverSize,rcvAverSize,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        global%mpiComm, global%mpierr )
    rnik = rcvAverSize
  ENDIF
#endif

! massflux ratio

  fact = tavMassflux/cvmMassflux

! stress difference

  DO j=jdcbeg,jdcend
    jid(j) = j
  ENDDO
  jid(jdcbeg-1) = jdcbeg
  jid(jdcend+1) = jdcend

  DO j=jdcbeg,jdcend
!    dy(j) = MAX( dy(j),dy1 )

    upl  = 0.5_RFREAL*(ua(jid(j+1)) + ua(j))
    umi  = 0.5_RFREAL*(ua(jid(j-1)) + ua(j))

    IF (j==jdcbeg .OR. j==jdcend) THEN
      dudy = 2._RFREAL*(upl-umi)/dy(j)
    ELSE
      dudy = (upl-umi)/dy(j)
    ENDIF

    restau = -rhoa(j)*uva(j)
    modtau = mueta(j)*dudy
    vistau = muela(j)*dudy

    tottau = MIN( restau+modtau+vistau,vistau1 )

    IF (tottau < 0._RFREAL) tottau = MAX( restau+modtau+vistau,-vistau1 )

    ynd = 1._RFREAL-(yc(j)-yn1)/global%refLength
    
!    dstress(j) = alpha*( ynd - tottau/vistau1 )/ynd

    IF (ynd >=  0._RFREAL) THEN
      dstress(j) = alpha*( ynd - tottau/vistau1 )/ &
                   MAX( restau/vistau1,0.05_RFREAL )
    ELSE
      dstress(j) = alpha*( tottau/vistau1 - ynd )/ &
                   MAX(-restau/vistau1,0.05_RFREAL )
    ENDIF

!    IF (dstress(j) < 0._RFREAL) dstress(j) = dstress(j)*ralpha
    
    IF (ABS( ynd ) < 0.1_RFREAL) THEN
      dstress(j) = 0._RFREAL
    ENDIF
write(*,100)region%iRegionGlobal,(region%iRegionGlobal-1)*8+j,dstress(j),ynd,tottau,vistau1,tottau/vistau1
!write(*,100)region%iRegionGlobal,(region%iRegionGlobal-1)*8+j,restau,modtau,vistau,vistau1,tottau/vistau1
  ENDDO

100 FORMAT( 2I5,5f20.10 )

! correct mean

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijk   = IndIJK(i,j,k,iOff,ijOff)

        ampli  = MAX( 1._RFREAL+dstress(j), 0._RFREAL )

        u = dv(DV_MIXT_UVEL,ijk)
        v = dv(DV_MIXT_VVEL,ijk)
        w = dv(DV_MIXT_WVEL,ijk)

        cv(CV_MIXT_DENS,ijk) = (cv(CV_MIXT_DENS,ijk)-rhom(j))*1._RFREAL+ &
                               beta*cvMean(j,CV_MIXT_DENS)+ ombeta*rhom(j)

        rho  = cv(CV_MIXT_DENS,ijk)
        rrho = 1._RFREAL/cvMean(j,CV_MIXT_DENS)

!        cv(CV_MIXT_XMOM,ijk) = (cv(CV_MIXT_XMOM,ijk)-rum(j))*ampli+ &
!                               beta*cvMean(j,CV_MIXT_XMOM)+ ombeta*rum(j)

!        cv(CV_MIXT_YMOM,ijk) = (cv(CV_MIXT_YMOM,ijk)-rvm(j))*ampli+ &
!                               beta*cvMean(j,CV_MIXT_YMOM)+ ombeta*rvm(j)

!        cv(CV_MIXT_ZMOM,ijk) = (cv(CV_MIXT_ZMOM,ijk)-rwm(j))*ampli+ &
!                               beta*cvMean(j,CV_MIXT_ZMOM)+ ombeta*rwm(j)

        cv(CV_MIXT_XMOM,ijk) = rho*( (u-um(j))*ampli + beta* &
                               cvMean(j,CV_MIXT_XMOM)*rrho + ombeta*um(j) )

        cv(CV_MIXT_YMOM,ijk) = rho*( (v-vm(j))*ampli + beta* &
                               cvMean(j,CV_MIXT_YMOM)*rrho + ombeta*vm(j) )

        cv(CV_MIXT_ZMOM,ijk) = rho*( (w-wm(j))*ampli + beta* &
                               cvMean(j,CV_MIXT_ZMOM)*rrho + ombeta*wm(j) )

        cv(CV_MIXT_ENER,ijk) = (cv(CV_MIXT_ENER,ijk)-rem(j))*ampli+ &
                               beta*cvMean(j,CV_MIXT_ENER)+ ombeta*rem(j)

!        mm    = gv(GV_MIXT_MOL ,ijk*indMol)
!        rgas  = MixtPerf_R_M( mm ) 
!        cpgas = gv(GV_MIXT_CP  ,ijk*indCp )
!        g     = MixtPerf_G_CpR( cpgas,rgas ) 

!        rho = cv(CV_MIXT_DENS,ijk)
!        ru  = cv(CV_MIXT_XMOM,ijk)
!        rv  = cv(CV_MIXT_YMOM,ijk)
!        rw  = cv(CV_MIXT_ZMOM,ijk)
!        pr  = dv(DV_MIXT_PRES,ijk)
!        re  = pr/(g-1._RFREAL)+0.5_RFREAL*(ru*ru+rv*rv+rw*rw)/rho

!        cv(CV_MIXT_ENER,ijk) = (re-rem(j))*1._RFREAL+ &
!                               beta*cvMean(j,CV_MIXT_ENER)+ ombeta*rem(j)

        cv(CV_MIXT_XMOM,ijk) = cv(CV_MIXT_XMOM,ijk)*fact
      ENDDO
    ENDDO
  ENDDO

! deallocate temporary arrays

  DEALLOCATE( rhom,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( um,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( vm,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( wm,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( rem,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE( rhoa,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( ua,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( va,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( wa,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( rea,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE( dy,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( yc,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( yn,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( muela,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( mueta,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( uva,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  DEALLOCATE( dstress,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE( jid,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! compute derived values for the mixture in case of dual time-stepping

  IF (global%flowType == FLOW_UNSTEADY) THEN
!    IF (global%solverType == SOLV_IMPLICIT) THEN   ! dual TST
      ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
      iec = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
      CALL MixtureProperties( region,ibc,iec,.true. )
!    ENDIF
  ENDIF
#endif

! finalize --------------------------------------------------------------------

999  CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PERI_CoMeanCorrection
#endif

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PERI_ModHybridDES

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_ModHybridDES.F90,v $
! Revision 1.9  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2005/04/21 00:46:19  wasistho
! modified PERI_RFLO_ReadMean
!
! Revision 1.6  2005/04/06 02:39:16  wasistho
! modified dstress in PERI_CoMeanCorrection
!
! Revision 1.5  2005/04/06 01:53:36  wasistho
! added amplified fluctuations
!
! Revision 1.4  2005/03/31 17:03:34  wasistho
! apply space averaged i.o. time averaged and computed dv,gv for dual TST
!
! Revision 1.3  2005/03/10 02:02:09  wasistho
! bug fixed for readMean and meanCorrection
!
! Revision 1.2  2005/03/07 18:27:56  wasistho
! insert ifdef STATS for tav
!
! Revision 1.1  2005/03/07 05:08:13  wasistho
! install hybrid DESSA turbulence model
!
!
! ******************************************************************************








