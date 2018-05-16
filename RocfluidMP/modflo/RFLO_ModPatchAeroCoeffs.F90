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
! Purpose: Collection of routines for patch aerodynamic coefficients.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLO_ModPatchAeroCoeffs.F90,v 1.5 2008/12/06 08:44:17 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModPatchAeroCoeffs

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_InitPatchAeroCoeffs, &
            RFLO_ReadPatchAeroCoeffs, &
            RFLO_WritePatchAeroCoeffs, &
            RFLO_ReadPatchAeroCoeffsReg

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModPatchAeroCoeffs.F90,v $ $Revision: 1.5 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  

! *******************************************************************************
!
! Purpose: Initialize patch aerodynamic coefficients.
!
! Description: None.
!
! Input:
!   region = current region data
!   patch  = current patch data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLO_InitPatchAeroCoeffs( region,patch )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables

! ... local variables

!********************************************************************************

  CALL RegisterFunction( region%global,&
  'RFLO_InitPatchAeroCoeffs',&
  'RFLO_ModPatchAeroCoeffs.F90' )
        
! set initial patch aero coeffs values to 0.0 -----------------------------------
        
  patch%cp           = 0.0_RFREAL
  patch%cf           = 0.0_RFREAL
  patch%ch           = 0.0_RFREAL
  patch%forceCoeffs  = 0.0_RFREAL
  patch%momentCoeffs = 0.0_RFREAL

! finalize ----------------------------------------------------------------------

  CALL DeregisterFunction(region%global)

END SUBROUTINE RFLO_InitPatchAeroCoeffs



! *******************************************************************************
!
! Purpose: Read patch coefficients.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input:
!   regions = dimensions of all regions.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLO_ReadPatchAeroCoeffs( regions )

  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, ic

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, flowModel, iFile, iRegFile, iPatchFile, n1File, n2File, iOff
  INTEGER :: n1, n2, n, nCoeffs, tag, nDimC, ijBeg, ijEnd, bcType, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), acFile(:,:)

  TYPE(t_patch),  POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,&
  'RFLO_ReadPatchAeroCoeffs',&
  'RFLO_ModPatchAeroCoeffs.F90' )
        
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Entering RFLO_ReadPatchAeroCoeffs...'
  END IF ! global%verbLevel        

! allocate fixed-size temporary data arrays -----------------------------------

  ALLOCATE( ivar(4,1),stat=errorFlag )
  ALLOCATE( rvar(1,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,&
  ERR_ALLOCATE,&
  __LINE__ )
    
  iFile = IF_PATCH_COEF   

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

! - unsteady flow

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.pcoa_', &
                                   global%timeStamp
        OPEN(iFile,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.pcob_', &
                                   global%timeStamp
        OPEN(iFile,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,&
        ERR_UNKNOWN_FORMAT,& 
        __LINE__ )
      ENDIF

! - steady flow

    ELSE
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.pcoa_', &
                                global%currentIter
        OPEN(iFile,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.pcob_', &
                                global%currentIter
        OPEN(iFile,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,&
        ERR_UNKNOWN_FORMAT,&
        __LINE__ )
      ENDIF
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,&
      ERR_FILE_OPEN,&
      __LINE__,&
      'File: '//TRIM(fname) )

  ENDIF   ! MASTERPROC

! read, broadcast and check time ----------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_ReadDataFileReal( global,iFile,global%solutFormat,1,1,rvar )
  ENDIF

#ifdef MPI
  CALL MPI_Bcast( rvar,1,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,&
  ERR_MPI_TROUBLE,&
  __LINE__ )
#endif

  IF (global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL) THEN
    IF (global%currentTime /= rvar(1,1)) THEN
      WRITE(msg,1000) rvar(1,1),global%currentTime
      CALL ErrorStop( global,ERR_TIME_SOLUTION,&
      __LINE__,&
      msg//' File: '//TRIM(fname) )
    ENDIF
  ENDIF

! read solution data ----------------------------------------------------------

  DO iReg=1,global%nRegions
    iLev      = regions(iReg)%currLevel
    flowModel = regions(iReg)%mixtInput%flowModel

    DO iPatch=1,regions(iReg)%nPatches
      patch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType = patch%bcType

      IF ((bcType < BC_REGIONCONF .OR. bcType > BC_REGIONCONF+BC_RANGE) .AND.&
          (bcType < BC_REGIONINT  .OR. bcType > BC_REGIONINT +BC_RANGE) .AND.&
          (bcType < BC_REGNONCONF .OR. bcType > BC_REGNONCONF+BC_RANGE)) THEN

! ----- get dimensions and pointers

        n1   = ABS(patch%l1end-patch%l1beg)
        n2   = ABS(patch%l2end-patch%l2beg)
        iOff = n1 + 1
        ijBeg = IndIJ( 0, 0,iOff)
        ijEnd = IndIJ(n1,n2,iOff)

        IF (flowModel==FLOW_EULER) nCoeffs = 1
        IF (flowModel==FLOW_NAVST) nCoeffs = 5
        nDimC = ijEnd - ijBeg + 1

        tag = regions(iReg)%localNumber + MPI_PATCHOFF*iPatch

! ----- read region and patch number, and dimensions (only master)

        IF (global%myProcid == MASTERPROC) THEN
          CALL RFLO_ReadDataFileInt( global,iFile,global%solutFormat,4,1,ivar )
          iRegFile   = ivar(1,1)
          iPatchFile = ivar(2,1)
          n1File     = ivar(3,1)
          n2File     = ivar(4,1)
          IF (iRegFile /= iReg) &
            CALL ErrorStop( global,ERR_REGION_NUMBER,&
            __LINE__,&
           'File: '//TRIM(fname) )
          IF (iPatchFile /= iPatch) &
            CALL ErrorStop( global,ERR_PATCH_NUMBER,&
            __LINE__, &
                            'File: '//TRIM(fname) )
          IF (n1File /= n1 .OR. n2File /= n2) THEN
            WRITE(msg,1000) iReg,iPatch, n1, n2
            CALL ErrorStop( global,ERR_PATCH_DIMENS,&
            __LINE__,&
            msg )
          ENDIF
        ENDIF

! ----- master reads & sends data, others receive them

        IF (global%myProcid == MASTERPROC) THEN

          ALLOCATE( acFile(nCoeffs,nDimC),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
          __LINE__ )

          CALL RFLO_ReadDataFileReal( global,iFile,global%solutFormat, &
                                      nCoeffs,nDimC,acFile )

#ifdef MPI
          IF (regions(iReg)%procid /= MASTERPROC) THEN
            CALL MPI_Send( acFile,nCoeffs*nDimC,MPI_RFREAL, &
                           regions(iReg)%procid,tag, &
                           global%mpiComm,global%mpierr )
            IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
            __LINE__ )
          ENDIF
#endif
        ELSE   ! not the master

          IF (regions(iReg)%procid == global%myProcid) THEN
            ALLOCATE( acFile(nCoeffs,nDimC),stat=errorFlag )
            global%error = errorFlag
            IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
            __LINE__ )

#ifdef MPI
            CALL MPI_Recv( acFile,nCoeffs*nDimC,MPI_RFREAL,MASTERPROC,tag, &
                           global%mpiComm,status,global%mpierr )
            IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
            __LINE__ )
#endif
          ENDIF ! my processor
        ENDIF   ! master or not

! ----- copy solution into data structure

        IF (regions(iReg)%procid == global%myProcid) THEN
          IF (flowModel == FLOW_EULER) THEN
            n = 0
            DO ic = ijBeg,ijEnd
              n   = n + 1
              patch%cp(ic) = acFile(1,n)
            ENDDO
          ELSEIF (flowModel == FLOW_NAVST) THEN
            n = 0
            DO ic = ijBeg,ijEnd
              n   = n + 1
              patch%cp(       ic) = acFile(1,n)
              patch%cf(XCOORD,ic) = acFile(2,n)
              patch%cf(YCOORD,ic) = acFile(3,n)
              patch%cf(ZCOORD,ic) = acFile(4,n)
              patch%ch(       ic) = acFile(5,n)
            ENDDO
          ENDIF
        ENDIF

        IF (ALLOCATED(acFile)) THEN
          DEALLOCATE( acFile,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
          __LINE__ )
        ENDIF
    
      ENDIF ! bcType
    ENDDO   ! iPatch
  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(iFile,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,&
      __LINE__,&
      'File: '//TRIM(fname) )
  ENDIF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Leaving RFLO_ReadPatchAeroCoeffs.'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', iPatch= ',I5,', n1= ',I5,', n2= ',I5,'.')

END SUBROUTINE RFLO_ReadPatchAeroCoeffs




! *******************************************************************************
!
! Purpose: Write patch coefficients.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input:
!   regions = dimensions of all regions.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLO_WritePatchAeroCoeffs( regions )

  USE ModInterfaces, ONLY : RFLO_WriteDataFileInt, RFLO_WriteDataFileReal

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, ic

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, flowModel, iFile, iRegFile, iPatchFile, n1File, n2File, iOff
  INTEGER :: n1, n2, n, nCoeffs, tag, nDimC, ijBeg, ijEnd, bcType, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), acFile(:,:)

  TYPE(t_patch),  POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_WritePatchAeroCoeffs',&
  'RFLO_ModPatchAeroCoeffs.F90' )
        
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Entering RFLO_WritePatchAeroCoeffs...'
  END IF ! global%verbLevel        

! allocate fixed-size temporary data arrays -----------------------------------

  ALLOCATE( ivar(4,1),stat=errorFlag )
  ALLOCATE( rvar(1,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
  __LINE__ )
    
  iFile = IF_PATCH_COEF   

! open solution file (only master proc.) --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

! - unsteady flow

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.pcoa_', &
                                   global%currentTime
        OPEN(iFile,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.pcob_', &
                                   global%currentTime
        OPEN(iFile,file=fname,form='unformatted',status='unknown', &
             iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,&
        __LINE__ )
      ENDIF
      rvar(1,1) = global%currentTime

! - steady flow

    ELSE
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.pcoa_', &
                                global%currentIter
        OPEN(iFile,file=fname,form='formatted',status='unknown', &
             iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.pcob_', &
                                global%currentIter
        OPEN(iFile,file=fname,form='unformatted',status='unknown', &
             iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,&
        __LINE__ )
      ENDIF
      rvar(1,1) = global%resInit
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,&
      __LINE__,&
      'File: '//TRIM(fname) )

  ENDIF   ! MASTERPROC

! read, broadcast and check time ----------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CALL RFLO_WriteDataFileReal( global,iFile,global%solutFormat,1,1,rvar )
  ENDIF

! write solution data ---------------------------------------------------------

  DO iReg=1,global%nRegions
    iLev      = regions(iReg)%currLevel
    flowModel = regions(iReg)%mixtInput%flowModel

    DO iPatch=1,regions(iReg)%nPatches
      patch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType = patch%bcType

      IF ((bcType < BC_REGIONCONF .OR. bcType > BC_REGIONCONF+BC_RANGE) .AND.&
          (bcType < BC_REGIONINT  .OR. bcType > BC_REGIONINT +BC_RANGE) .AND.&
          (bcType < BC_REGNONCONF .OR. bcType > BC_REGNONCONF+BC_RANGE)) THEN

        IF (regions(iReg)%procid==global%myProcid .OR. &
            global%myProcid==MASTERPROC) THEN

! ------- get dimensions and pointers

          n1   = ABS(patch%l1end-patch%l1beg)
          n2   = ABS(patch%l2end-patch%l2beg)
          iOff = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)

          IF (flowModel==FLOW_EULER) nCoeffs = 1
          IF (flowModel==FLOW_NAVST) nCoeffs = 5
          nDimC = ijEnd - ijBeg + 1

          tag = regions(iReg)%localNumber + MPI_PATCHOFF*iPatch

! ------- allocate temporary array

          ALLOCATE( acFile(nCoeffs,nDimC),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
          __LINE__ )

        ENDIF  ! region%procid

! ----- copy solution into data structure

        IF (regions(iReg)%procid == global%myProcid) THEN
          IF (flowModel == FLOW_EULER) THEN
            n = 0
            DO ic = ijBeg,ijEnd
              n   = n + 1
              acFile(1,n) = patch%cp(ic)
            ENDDO
          ELSEIF (flowModel == FLOW_NAVST) THEN
            n = 0
            DO ic = ijBeg,ijEnd
              n   = n + 1
              acFile(1,n) = patch%cp(       ic)
              acFile(2,n) = patch%cf(XCOORD,ic)
              acFile(3,n) = patch%cf(YCOORD,ic)
              acFile(4,n) = patch%cf(ZCOORD,ic)
              acFile(5,n) = patch%ch(       ic)
            ENDDO
          ENDIF

        ENDIF  ! region%procid

! ----- write region and patch number, and dimensions (only master)

        IF (global%myProcid == MASTERPROC) THEN
          ivar(1,1) = iReg
          ivar(2,1) = iPatch
          ivar(3,1) = n1
          ivar(4,1) = n2
          CALL RFLO_WriteDataFileInt( global,iFile,global%solutFormat,4,1,ivar )
        ENDIF

! ----- master receives and writes data, others send them

        IF (global%myProcid == MASTERPROC) THEN
#ifdef MPI
          IF (regions(iReg)%procid /= MASTERPROC) THEN
            CALL MPI_Recv( acFile,nCoeffs*nDimC,MPI_RFREAL, &
                           regions(iReg)%procid,tag, &
                           global%mpiComm,status,global%mpierr )
            IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
            __LINE__ )
          ENDIF
#endif
          CALL RFLO_WriteDataFileReal( global,iFile,global%solutFormat, &
                                       nCoeffs,nDimC,acFile )
        ELSE   ! not the master

          IF (regions(iReg)%procid == global%myProcid) THEN
#ifdef MPI
            CALL MPI_Send( acFile,nCoeffs*nDimC,MPI_RFREAL,MASTERPROC,tag, &
                           global%mpiComm,global%mpierr )
            IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
            __LINE__ )
#endif
          ENDIF ! my processor
        ENDIF   ! master or not

        IF (ALLOCATED(acFile)) THEN
          DEALLOCATE( acFile,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
          __LINE__ )
        ENDIF
    
      ENDIF ! bcType
    ENDDO   ! iPatch
  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(iFile,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,&
      __LINE__,&
      'File: '//TRIM(fname) )
  ENDIF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Leaving RFLO_WritePatchAeroCoeffs.'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', iPatch= ',I5,', n1= ',I5,', n2= ',I5,'.')

END SUBROUTINE RFLO_WritePatchAeroCoeffs



! *******************************************************************************
!
! Purpose: Read patch coefficients regionwise.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input:
!   regions = dimensions of all regions.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLO_ReadPatchAeroCoeffsReg( iReg,regions )

  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iPatch, ic

! ... local variables
  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg

  INTEGER :: iLev, flowModel, iFile, iRegFile, iPatchFile, n1File, n2File, iOff
  INTEGER :: n1, n2, n, nCoeffs, nDimC, ijBeg, ijEnd, bcType, errorFlag
  INTEGER, ALLOCATABLE :: ivar(:,:)

  REAL(RFREAL), ALLOCATABLE :: rvar(:,:), acFile(:,:)

  TYPE(t_patch),  POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_ReadPatchAeroCoeffsReg',&
  'RFLO_ModPatchAeroCoeffs.F90' )
        
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Entering RFLO_ReadPatchAeroCoeffsReg...'
  END IF ! global%verbLevel        

! allocate fixed-size temporary data arrays -----------------------------------

  ALLOCATE( ivar(4,1),stat=errorFlag )
  ALLOCATE( rvar(1,1),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
  __LINE__ )
    
  iFile = IF_PATCH_COEF   

! open solution file (only if iReg=1) -----------------------------------------

  IF (iReg == 1) THEN

! - unsteady flow

    IF (global%flowType == FLOW_UNSTEADY) THEN
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.pcoa_', &
                                   global%timeStamp
        OPEN(iFile,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.pcob_', &
                                   global%timeStamp
        OPEN(iFile,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,&
        __LINE__ )
      ENDIF

! - steady flow

    ELSE
      IF (global%solutFormat == FORMAT_ASCII) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.pcoa_', &
                                global%currentIter
        OPEN(iFile,file=fname,form='formatted',status='old',iostat=errorFlag)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.pcob_', &
                                global%currentIter
        OPEN(iFile,file=fname,form='unformatted',status='old',iostat=errorFlag)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,&
        __LINE__ )
      ENDIF
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,&
      __LINE__,&
      'File: '//TRIM(fname) )

! - read and check time -------------------------------------------------------

    IF (global%myProcid == MASTERPROC) THEN
      CALL RFLO_ReadDataFileReal( global,iFile,global%solutFormat,1,1,rvar )
    ENDIF

    IF (global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL) THEN
      IF (global%currentTime /= rvar(1,1)) THEN
        WRITE(msg,1000) rvar(1,1),global%currentTime
        CALL ErrorStop( global,ERR_TIME_SOLUTION,&
        __LINE__,&
        msg//' File: '//TRIM(fname) )
      ENDIF
    ENDIF

  ENDIF   ! 1st region

! read solution data ----------------------------------------------------------

  iLev      = regions(iReg)%currLevel
  flowModel = regions(iReg)%mixtInput%flowModel

  DO iPatch=1,regions(iReg)%nPatches
    patch => regions(iReg)%levels(iLev)%patches(iPatch)
    bcType = patch%bcType


    IF ((bcType < BC_REGIONCONF .OR. bcType > BC_REGIONCONF+BC_RANGE) .AND.&
        (bcType < BC_REGIONINT  .OR. bcType > BC_REGIONINT +BC_RANGE) .AND.&
        (bcType < BC_REGNONCONF .OR. bcType > BC_REGNONCONF+BC_RANGE)) THEN

! --- get dimensions and pointers

      n1   = ABS(patch%l1end-patch%l1beg)
      n2   = ABS(patch%l2end-patch%l2beg)
      iOff = n1 + 1
      ijBeg = IndIJ( 0, 0,iOff)
      ijEnd = IndIJ(n1,n2,iOff)

      IF (flowModel==FLOW_EULER) nCoeffs = 1
      IF (flowModel==FLOW_NAVST) nCoeffs = 5
      nDimC = ijEnd - ijBeg + 1

! --- read region and patch number, and dimensions (only master)

      CALL RFLO_ReadDataFileInt( global,iFile,global%solutFormat,4,1,ivar )
      iRegFile   = ivar(1,1)
      iPatchFile = ivar(2,1)
      n1File     = ivar(3,1)
      n2File     = ivar(4,1)
      
      IF (iRegFile /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,&
        __LINE__, &
                        'File: '//TRIM(fname) )
      IF (iPatchFile /= iPatch) &
        CALL ErrorStop( global,ERR_PATCH_NUMBER,&
        __LINE__, &
                        'File: '//TRIM(fname) )
      IF (n1File /= n1 .OR. n2File /= n2) THEN
        WRITE(msg,1000) iReg,iPatch, n1, n2
        CALL ErrorStop( global,ERR_PATCH_DIMENS,&
        __LINE__,&
        msg )
      ENDIF

! --- read data

      ALLOCATE( acFile(nCoeffs,nDimC),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL RFLO_ReadDataFileReal( global,iFile,global%solutFormat, &
                                  nCoeffs,nDimC,acFile )

! --- copy solution into data structure

      IF (flowModel == FLOW_EULER) THEN
        n = 0
        DO ic = ijBeg,ijEnd
          n   = n + 1
          patch%cp(ic) = acFile(1,n)
        ENDDO
      ELSEIF (flowModel == FLOW_NAVST) THEN
        n = 0
        DO ic = ijBeg,ijEnd
          n   = n + 1
          patch%cp(       ic) = acFile(1,n)
          patch%cf(XCOORD,ic) = acFile(2,n)
          patch%cf(YCOORD,ic) = acFile(3,n)
          patch%cf(ZCOORD,ic) = acFile(4,n)
          patch%ch(       ic) = acFile(5,n)
        ENDDO
      ENDIF

      IF (ALLOCATED(acFile)) THEN
        DEALLOCATE( acFile,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
      ENDIF
    
    ENDIF ! bcType
  ENDDO   ! iPatch

! finalize --------------------------------------------------------------------

  IF (iReg == global%nRegions) THEN
    CLOSE(iFile,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Leaving RFLO_ReadPatchAeroCoeffsReg.'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', iPatch= ',I5,', n1= ',I5,', n2= ',I5,'.')

END SUBROUTINE RFLO_ReadPatchAeroCoeffsReg



! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModPatchAeroCoeffs

! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLO_ModPatchAeroCoeffs.F90,v $
!   Revision 1.5  2008/12/06 08:44:17  mtcampbe
!   Updated license.
!
!   Revision 1.4  2008/11/19 22:17:28  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.3  2006/03/25 04:33:49  wasistho
!   incorrect placing of status in mpi send/rcv of write routine
!
!   Revision 1.2  2006/03/24 04:56:18  wasistho
!   initialized patch force moments coeffs
!
!   Revision 1.1  2006/03/22 03:06:50  wasistho
!   initial import
!
!
! ******************************************************************************










