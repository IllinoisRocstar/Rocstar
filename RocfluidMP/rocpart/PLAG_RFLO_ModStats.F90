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
! Purpose: Collection of routines for particle statistics on Eulerian grid.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_RFLO_ModStats.F90,v 1.9 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_RFLO_ModStats

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModPartLag, ONLY: t_buffer_plag
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag 
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters

#include "Indexing.h" 
  USE ModIndexing,   ONLY: IndIJKMap
  USE ModInterfaces, ONLY: RFLO_GetCellOffSet,     &
                           RFLO_GetDimensDummy,    &
                           RFLO_GetPatchIndices,   &
                           RFLO_GetPatchDirection, &
                           RFLO_GetPatchMapping,   &
                           RFLO_ReadDataFileInt,   &
                           RFLO_ReadDataFileReal,  &
                           RFLO_WriteDataFileInt,  &
                           RFLO_WriteDataFileReal

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: PLAG_RFLO_ModStats.F90,v $ $Revision: 1.9 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: PLAG_RFLO_CreateStatBuff,      &
            PLAG_RFLO_CommStatBuffWrapper, &
            PLAG_RFLO_ReadStat,            &
            PLAG_RFLO_WriteStat

! ==============================================================================
! Private functions
! ==============================================================================

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
  
  
  
  
  

  





!******************************************************************************
!
! Purpose: Clear communication requests.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_RFLO_ClearReqStatBuff( regions,iReg )

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iLev,iPatch,iRegSrc,iRequestStat,nPatches 
#ifdef MPI
    INTEGER :: status(MPI_STATUS_SIZE)
#endif
    LOGICAL :: doWait

    TYPE(t_patch),    POINTER :: pPatch
    TYPE(t_plag),     POINTER :: pPlag
    TYPE(t_region),   POINTER :: pRegion
    TYPE(t_global),   POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(iReg)%global

    CALL RegisterFunction( global,'PLAG_RFLO_ClearReqStatBuff',&
  'PLAG_RFLO_ModStats.F90' )

! ******************************************************************************
!   Set dimensions, variables, pointers
! ******************************************************************************    

#ifdef MPI
    iLev     = regions(iReg)%currLevel
    nPatches = regions(iReg)%nPatches

    pPlag => regions(iReg)%levels(iLev)%plag 

! ****************************************************************************** 
!   Wait for patch data being received by other processors
! ******************************************************************************

    DO iPatch=1,nPatches
      pPatch => regions(iReg)%levels(iLev)%patches(iPatch)

      bcType  = pPatch%bcType
      iRegSrc = pPatch%srcRegion
      iRequestStat = pPatch%bufferPlag%iRequestStat

      doWait = ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
                (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
                (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
                (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
                (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE))

      IF ( iRegSrc > 0 ) THEN
        IF ( (doWait .EQV. .TRUE.) .AND.                  &
             (regions(iRegSrc)%procid /= global%myProcid) ) THEN
          CALL MPI_Wait( pPlag%requestsStat(iRequestStat),status,global%mpierr )
          IF ( global%mpierr /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        ENDIF ! doWait 
      ENDIF ! iRegSrc

    ENDDO   ! iPatch
#endif

! ******************************************************************************
!   Finalize
! ******************************************************************************

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_ClearReqStatBuff






! ******************************************************************************
!
! Purpose: Create buffer arrays for statistics.
!
! Description: None.
!
! Input: regions = dimensions of all regions
!        iReg    = current region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLO_CreateStatBuff(regions,iReg)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg 

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,errorFlag,iLev,iPatch,iRegSrc,n1,n2,n1Src,n2Src,ndc,&
               ndcSrc,ndim,ndimSrc,nEqs,nEqsSrc,nEv,nTav

    TYPE(t_patch),       POINTER :: pPatch
    TYPE(t_buffer_plag), POINTER :: pBuffPlag
    TYPE(t_plag),        POINTER :: pPlag
    TYPE(t_global),      POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(iReg)%global

    CALL RegisterFunction( global,'PLAG_RFLO_CreateStatBuff',&
  'PLAG_RFLO_ModStats.F90' )

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,&
                               'Allocating Statistics Data Buffers for PLAG...'
    END IF ! global%verbLevel
    
    nTav = global%plagNStat

! ******************************************************************************
!   Allocate buffer data arrays for statistics
! ******************************************************************************  

    DO iLev=1,regions(iReg)%nGridLevels
      pPlag => regions(iReg)%levels(iLev)%plag
      pPlag%nRequestsStat = 0
   
      nEv = pPlag%nEv
      
      DO iPatch=1,regions(iReg)%nPatches
  
        pPatch => regions(iReg)%levels(iLev)%patches(iPatch)
        pBuffPlag => pPatch%bufferPlag
      
        bcType = pPatch%bcType

        IF ( (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
             (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
             (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc = pPatch%srcRegion
        
          IF ( regions(iRegSrc)%procid /= global%myProcid ) THEN  ! other processor
            n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 2    ! large enough
            n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 2    ! for NODES!
            n1Src   = ABS(pPatch%srcL1end-pPatch%srcL1beg) + 2
            n2Src   = ABS(pPatch%srcL2end-pPatch%srcL2beg) + 2
            nEqs    = nTav +nEv     
            nEqsSrc = nTav +nEv
            ndc     = regions(iReg   )%nDumCells
            ndcSrc  = regions(iRegSrc)%nDumCells
            ndim    = n1*n2*nEqs*ndc
            ndimSrc = n1Src*n2Src*nEqsSrc*ndcSrc

            ALLOCATE( pBuffPlag%sendBuffStat(ndimSrc),stat=errorFlag )
            global%error = errorFlag
            IF (global%error /= ERR_NONE) THEN
              CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%sendBuffStat' ) 
            END IF ! global%error

            ALLOCATE( pBuffPlag%recvBuffStat(ndim   ),stat=errorFlag )
            global%error = errorFlag
            IF (global%error /= ERR_NONE) THEN
              CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pBuffPlag%recvBuffStat' ) 
            END IF ! global%error

            pBuffPlag%nSendBuffStat = ndimSrc
            pBuffPlag%nRecvBuffStat = ndim
            pPlag%nRequestsStat     = pPlag%nRequestsStat + 1
            pBuffPlag%iRequestStat  = pPlag%nRequestsStat

          ENDIF ! regions(iRegSrc)%procid 

        ELSE IF ( (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
                  (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE)) THEN
          CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )  ! #### TEMPORARY ####

        ELSE 
          NULLIFY(pBuffPlag%sendBuffStat)
          NULLIFY(pBuffPlag%recvBuffStat)
        ENDIF  ! bcType

      ENDDO    ! iPatch

! ==============================================================================
!     Allocate array for send requests
! ==============================================================================

      ALLOCATE( pPlag%requestsStat(pPlag%nRequestsStat),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pPlag%requests' ) 
      END IF ! global%error

    ENDDO      ! iLev

! ******************************************************************************
!   Finalize
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,&
                           'Allocating Statistics Data Buffers for PLAG done...'
    END IF ! global%verbLevel

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_CreateStatBuff






  
! ******************************************************************************
!
! Purpose: Communicate buffer arrays for statistics.
!
! Description: None.
!
! Input: regions = dimensions of all regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLO_CommStatBuffWrapper(regions)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iLev,iPatch,iPatchSrc,iReg,iRegSrc,nPatches
  
    TYPE(t_patch),       POINTER :: pPatch,pPatchSrc
    TYPE(t_region),      POINTER :: pRegion,pRegionSrc
    TYPE(t_global),      POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(1)%global

    CALL RegisterFunction( global,'PLAG_RFLO_CommStatBuffWrapper',&
  'PLAG_RFLO_ModStats.F90' )

    IF (.NOT. global%plagUsed) GOTO 999

! ******************************************************************************
!   Copy statistics for regions on the same processor
! ******************************************************************************  

    DO iReg=1,global%nRegions
      pRegion => regions(iReg)
      
      IF ( pRegion%procid==global%myProcid .AND. &   ! region active and
           pRegion%active==ACTIVE ) THEN             ! on my processor

! ==============================================================================
!        Set dimensions and pointers
! ==============================================================================

        iLev     = pRegion%currLevel
        nPatches = pRegion%nPatches


! ==============================================================================         
!        Loop over patches 
! ==============================================================================
 
        DO iPatch=1,nPatches
          pPatch => pRegion%levels(iLev)%patches(iPatch)

          bcType    = pPatch%bcType
          iRegSrc   = pPatch%srcRegion
          iPatchSrc = pPatch%srcPatch

          SELECT CASE (bcType)

! ------------------------------------------------------------------------------
!           Conforming region interface
! ------------------------------------------------------------------------------

            CASE( BC_REGIONCONF:BC_REGIONCONF+BC_RANGE )
              pRegionSrc => regions(iRegSrc)
              pPatchSrc  => pRegionSrc%levels(iLev)%patches(iPatchSrc)

              IF ( regions(iRegSrc)%procid == global%myProcid ) THEN
                CALL PLAG_RFLO_CopyStatBuff( pRegion,pRegionSrc, &
                                             pPatch,pPatchSrc    )
              ENDIF ! regions(iRegSrc)%procid

! ------------------------------------------------------------------------------
!           Non-conforming region interface (integer)
! ------------------------------------------------------------------------------

            CASE( BC_REGIONINT:BC_REGIONINT+BC_RANGE )
              CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )

! ------------------------------------------------------------------------------
!           Non-conforming region interface (irregular)
! ------------------------------------------------------------------------------
             
            CASE( BC_REGNONCONF:BC_REGNONCONF+BC_RANGE)
              CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )         
          
          END SELECT! bcType 

        END DO ! iPatch
      END IF ! pRegion%procid
    END DO ! iReg      

! ******************************************************************************
!   Communicate buffer for off-processor regions
! ******************************************************************************  

! ==============================================================================
!  Send buffer data
! ==============================================================================

#ifdef MPI
    DO iReg = 1, global%nRegions
      IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
           regions(iReg)%active==ACTIVE ) THEN             ! on my processor        
        CALL PLAG_RFLO_SendStatBuffWrapper( regions, iReg )    
      ENDIF       ! regions    
    ENDDO         ! iReg

! ==============================================================================                     
!  Receive buffer data 
! ==============================================================================
 
    DO iReg = 1, global%nRegions
      IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
           regions(iReg)%active==ACTIVE ) THEN             ! on my processor        
        CALL PLAG_RFLO_RecvStatBuffWrapper( regions, iReg )
      ENDIF ! regions
    ENDDO ! iReg 

! ==============================================================================
! wait for data being received by other processors ----------------------------
! ==============================================================================

    DO iReg = 1, global%nRegions
      IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
           regions(iReg)%active==ACTIVE ) THEN             ! on my processor        
        CALL PLAG_RFLO_ClearReqStatBuff( regions, iReg )
      ENDIF ! regions  
    ENDDO ! iReg
#endif

! ******************************************************************************
!   Finalize
! ******************************************************************************

! ******************************************************************************
!   End  
! ******************************************************************************

999  CONTINUE
    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_CommStatBuffWrapper






  
! ******************************************************************************
!
! Purpose: Copy buffer statistics arrays for on-processor regions.
!
! Description: None.
!
! Input: pRegion    = current region
!        pRegionSrc = source region
!        pPatch     = current patch
!        pPatchSrc  = source patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLO_CopyStatBuff( pRegion,pRegionSrc,pPatch,pPatchSrc )

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_patch),   POINTER :: pPatch,pPatchSrc
    TYPE(t_region),  POINTER :: pRegion,pRegionSrc

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: idum,iv,i,j,k,ii,jj,kk
    INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, &
               iCOff, ijCOff, ijkD, iLev
    INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc, &
               idirSrc, jdirSrc, kdirSrc, iCOffSrc, ijCOffSrc, ijkCSrc
    INTEGER :: lb, lbs, l1SrcDir, l2SrcDir, mapMat(3,4)
    INTEGER :: ivEvBeg,ivEvEnd,ivTavBeg,ivTavEnd

    LOGICAL :: align
    
    REAL(KIND=RFREAL), DIMENSION(:,:), POINTER :: ev,evSrc,tav,tavSrc
  
    TYPE(t_plag),        POINTER :: pPlag
    TYPE(t_global),      POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction( global,'PLAG_RFLO_CopyStatBuff',&
  'PLAG_RFLO_ModStats.F90' )

! ******************************************************************************
!   Check if the source region is active
! ******************************************************************************
    
    IF ( pRegionSrc%active == OFF ) THEN
      CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
    ENDIF

! ******************************************************************************
!   Set dimensions and pointers
! ******************************************************************************

    iLev = pRegion%currLevel

    CALL RFLO_GetPatchIndices( pRegion,pPatch,iLev,ibeg,iend, &
                               jbeg,jend,kbeg,kend )
    CALL RFLO_GetPatchIndices( pRegionSrc,pPatchSrc,iLev,ibegSrc,iendSrc, &
                               jbegSrc,jendSrc,kbegSrc,kendSrc )
    CALL RFLO_GetPatchDirection( pPatch   ,idir   ,jdir,   kdir    )
    CALL RFLO_GetPatchDirection( pPatchSrc,idirSrc,jdirSrc,kdirSrc )
    CALL RFLO_GetCellOffset( pRegion   ,iLev,iCOff   ,ijCOff    )
    CALL RFLO_GetCellOffset( pRegionSrc,iLev,iCOffSrc,ijCOffSrc )
    
    ev    => pRegion%levels(iLev)%plag%ev
    evSrc => pRegionSrc%levels(iLev)%plag%ev
    
    ivEvBeg = 1
    ivEvEnd = SIZE(ev,DIM=1)
        
    tav    => pRegion%levels(iLev)%plag%tav
    tavSrc => pRegionSrc%levels(iLev)%plag%tav
    
    ivTavBeg = 1
    ivTavEnd = SIZE(tav,DIM=1)
    
! ******************************************************************************
!   Mapping between patches
! ******************************************************************************

                                           l1SrcDir =  1
    IF (pPatch%srcL1beg > pPatch%srcL1end) l1SrcDir = -1
                                           l2SrcDir =  1
    IF (pPatch%srcL2beg > pPatch%srcL2end) l2SrcDir = -1

    lb    = pPatch%lbound
    lbs   = pPatch%srcLbound
    align = pPatch%align

    CALL RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                               idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                               ibeg,iend,jbeg,jend,kbeg,kend, &
                               ibegSrc,iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc, &
                               mapMat )

! ******************************************************************************
!   Loop over dummy nodes of current patch
! ******************************************************************************
    
    DO idum=1,pRegion%nDumCells
    DO k=kbeg,kend
    DO j=jbeg,jend
    DO i=ibeg,iend
      ii      = i - idum*idir
      jj      = j - idum*jdir
      kk      = k - idum*kdir
      ijkD    = IndIJK(ii,jj,kk,iCOff,ijCOff)
      ijkCSrc = IndIJKMap(ii,jj,kk,mapMat,iCOffSrc,ijCOffSrc)

      DO iv=ivEvBeg,ivEvEnd
         ev(iv,ijkD) =  evSrc(iv,ijkCSrc)
      END DO ! iVar

      DO iv=ivTavBeg,ivTavEnd
        tav(iv,ijkD) = tavSrc(iv,ijkCSrc)
      END DO ! iVar      
    ENDDO  ! i
    ENDDO  ! j
    ENDDO  ! k
    ENDDO  ! idum

! ******************************************************************************
!   Finalize
! ******************************************************************************

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_CopyStatBuff






  

! ******************************************************************************
!
! Purpose: Read in time averaged statistics of the Lagrangian particles 
!          on Eulerian grid.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary
!
! Input: regions = dimensions of all regions
!
! Output: region%levels%plag%tav = time avg Lagrangian particle variables
!         global%integrTime      = integrated averaging time
!
! Notes: time averaged solution is read in only for the current grid level;
!        it is also read in for all dummy cells
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLO_ReadStat(regions)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: regions(:)  

! ==============================================================================
!   Locals
! ==============================================================================

  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg

  INTEGER :: iReg, i, j, k, l, n, ind

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev,iRegFile,ipc,jpc,kpc,nDumCells,nDim,iOff,ijOff,ijk
  INTEGER :: idcbeg,jdcbeg,kdcbeg,idcend,jdcend,kdcend,ijkBeg,ijkEnd
  INTEGER :: errorFlag,nTav,nTavVar

  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ivar,jvar,plagVarId

  REAL(RFREAL), POINTER, DIMENSION(:,:)     :: tav
  REAL(RFREAL), ALLOCATABLE, DIMENSION(:,:) :: rvar, tavFile

  TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(1)%global

    CALL RegisterFunction( global,'PLAG_RFLO_ReadStat',&
  'PLAG_RFLO_ModStats.F90' )

! ******************************************************************************
!   Allocate temporary data arrays 
! ******************************************************************************

    ALLOCATE( ivar(5,1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'ivar' )
    
    ALLOCATE( rvar(2,1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'rvar' )
    
    ALLOCATE( jvar(global%plagNStat+1,1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'jvar' ) 
    
    ALLOCATE( plagVarId(2,global%plagNStat+1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'plagVarId' ) 

! ******************************************************************************
!   Open statistics file (only master proc.)
! ******************************************************************************

    IF (global%myProcid == MASTERPROC) THEN

      SELECT CASE( global%solutFormat )
        CASE ( FORMAT_ASCII )
          WRITE(fname,'(A,1PE11.5)') &
            TRIM(global%inDir)//TRIM(global%casename)//'.plag_stata_', &
            global%timeStamp
          OPEN( IF_PLAG_STATS,file=fname,form='formatted',status='old', &
                iostat=errorFlag )
      
        CASE ( FORMAT_BINARY )
          WRITE(fname,'(A,1PE11.5)') &
            TRIM(global%inDir)//TRIM(global%casename)//'.plag_stat_', &
            global%timeStamp
          OPEN( IF_PLAG_STATS,file=fname,form='unformatted',status='old', &
                iostat=errorFlag )

        CASE DEFAULT
          CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

      END SELECT ! solutFormat 

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

    ENDIF ! global%myProcid

! ******************************************************************************
!   Read & broadcast current and integrated time in file, and stats ID 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC ) THEN
      CALL RFLO_ReadDataFileReal( global,IF_PLAG_STATS,global%solutFormat,2,1,rvar )
    ENDIF ! global%myProcid

#ifdef MPI
    CALL MPI_Bcast( rvar,2,MPI_RFREAL,MASTERPROC,global%mpiComm,global%mpierr )
    IF (global%mpierr /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! ==============================================================================
!   Trap error for inconsistent variables in header 
! ==============================================================================
    
    IF ( global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL ) THEN
      IF ( global%currentTime /= rvar(1,1) ) THEN
        WRITE(msg,1000) rvar(1,1),global%currentTime
        CALL ErrorStop( global,ERR_TIME_SOLUTION,__LINE__,msg//' File: '//TRIM(fname) )
      ENDIF ! currentTime
      
      IF ( global%integrTime  /= rvar(2,1) ) THEN
        WRITE(msg,2000) rvar(2,1),global%integrTime
        CALL ErrorStop( global,ERR_TIME_SOLUTION,__LINE__,msg//' File: '//TRIM(fname) )
      ENDIF ! integrTime
            
    ENDIF ! global%flowType

! ******************************************************************************
!   Read plagNStat and plagStatId from file
! ******************************************************************************
    
    IF ( global%myProcid == MASTERPROC ) THEN
      IF (global%plagNStat > 0) THEN
        CALL RFLO_ReadDataFileInt( global,IF_PLAG_STATS,global%solutFormat, &
                                   global%plagNStat+1,1,jvar )
        nTavVar  = jvar(1,1)
        IF ( nTavVar /= global%plagNStat ) THEN
         CALL ErrorStop( global,ERR_STATS_RESTART,__LINE__ )
        END IF ! nTavVar

        plagVarId(1,:) = jvar(2:global%plagNStat+1,1)
        plagVarId(2,:) = MOD(plagVarId(1,:),10)
        plagVarId(1,:) = (plagVarId(1,:)-plagVarId(2,:))/10

        DO ind=1,2
        DO l=1,global%plagNStat
          IF ( plagVarId(ind,l) /= global%plagStatId(ind,l) ) &
            CALL ErrorStop( global,ERR_STATS_RESTART,__LINE__ )
        END DO ! l
        END DO ! ind
      ENDIF ! plagNStat

    ENDIF ! myProcid 
      
! ******************************************************************************
!   Read statistics data from all regions 
! ******************************************************************************

    DO iReg=1,global%nRegions

! ==============================================================================
!     Get dimensions and pointers
! ==============================================================================

      iLev = regions(iReg)%currLevel
      CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                                jdcbeg,jdcend,kdcbeg,kdcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iOff,ijOff )
      ijkBeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
      ijkEnd = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
      nDim   = ijkEnd - ijkBeg + 1

      nTav =global%plagNStat

! ==============================================================================
!     Read region number and dimensions (only master)
! ==============================================================================

      IF (global%myProcid == MASTERPROC) THEN
        CALL RFLO_ReadDataFileInt( global,IF_PLAG_STATS,global%solutFormat,5,1,ivar )
        iRegFile  = ivar(1,1)
        ipc       = ivar(2,1)
        jpc       = ivar(3,1)
        kpc       = ivar(4,1)
        nDumCells = ivar(5,1)
      
        IF (iRegFile /= iReg) &
          CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__,'File: '//TRIM(fname) )
        
        IF ( (ipc /= regions(iReg)%levels(iLev)%grid%ipc) .OR. &
             (jpc /= regions(iReg)%levels(iLev)%grid%jpc) .OR. &
             (kpc /= regions(iReg)%levels(iLev)%grid%kpc)      ) THEN
          WRITE(msg,1005) iReg,ipc,jpc,kpc
          CALL ErrorStop( global,ERR_GRID_DIMENSIONS,__LINE__,msg )
        ENDIF ! ipc

        IF ( nDumCells /= regions(iReg)%nDumCells ) THEN
          WRITE(msg,1010) iReg,nDumCells,regions(iReg)%nDumCells
          CALL ErrorStop( global,ERR_GRID_DUMCELLS,__LINE__,msg )
        ENDIF ! nDumCells

! ==============================================================================
!      Master reads & sends data, others receive them
! ==============================================================================

        ALLOCATE( tavFile(nTav,nDim),stat=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

        CALL RFLO_ReadDataFileReal( global,IF_PLAG_STATS,global%solutFormat, &
                                    nTav,nDim,tavFile )

#ifdef MPI
        IF ( regions(iReg)%procid /= MASTERPROC ) THEN
          CALL MPI_Send( tavFile,nTav*nDim,MPI_RFREAL, &
                         regions(iReg)%procid,iReg,   &
                         global%mpiComm,global%mpierr )
          IF (global%mpierr /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        ENDIF ! regions(iReg)%procid 
#endif
      ELSE   ! not the master

        IF ( regions(iReg)%procid == global%myProcid ) THEN
          ALLOCATE( tavFile(nTav,nDim),stat=errorFlag )
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

#ifdef MPI
          CALL MPI_Recv( tavFile,nTav*nDim,MPI_RFREAL,MASTERPROC, &
                         iReg,global%mpiComm,status,global%mpierr )
          IF ( global%mpierr /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
        ENDIF ! regions(iReg)%procid

      ENDIF  ! global%myProcid

! ==============================================================================
!    Copy statistics into data structure
! ==============================================================================
     
      IF ( regions(iReg)%procid == global%myProcid ) THEN
        tav => regions(iReg)%levels(iLev)%plag%tav

        n = 0
        DO k=kdcbeg,kdcend
        DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          n   = n + 1
          ijk = IndIJK(i,j,k,iOff,ijOff)
          DO l=1,nTav
            tav(l,ijk) = tavFile(l,n)
          ENDDO ! l
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k
      ENDIF ! regions(iReg)%procid

! ==============================================================================
!     Deallocate local array
! ==============================================================================

      IF ( ALLOCATED(tavFile) ) THEN
        DEALLOCATE( tavFile,stat=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'tavFile' )
      ENDIF ! tavFile

    ENDDO     ! iReg

! ******************************************************************************
!   Deallocate temporary data arrays 
! ******************************************************************************

    DEALLOCATE( ivar,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'ivar' )
    
    DEALLOCATE( rvar,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'rvar' )
    
    DEALLOCATE( jvar,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'jvar' )
    
    DEALLOCATE( plagVarId,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'plagVarId' )

! ******************************************************************************
!   Finalize
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC ) THEN
      CLOSE(IF_PLAG_STATS,iostat=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
    ENDIF

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

! ******************************************************************************
!   Formats 
! ******************************************************************************

1000 FORMAT('Time in file is= ',1PE12.5,' but it should be= ',E12.5,'.')
1005 FORMAT('Region ',I5,', ipc= ',I6,', jpc= ',I6,', kpc= ',I6,'.')
1010 FORMAT('Region ',I5,', # dummy cells=',I2,' but should be= ',I1)
2000 FORMAT('Integration Time in file is= ',1PE12.5,' but it should be= ',E12.5,'.')

  END SUBROUTINE PLAG_RFLO_ReadStat






!******************************************************************************
!
! Purpose: Wrapper to receive data.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%plag%tav = updated plag statitics values
!                                         in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_RFLO_RecvStatBuffWrapper( regions,iReg )

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iLev,iPatch,iPatchSrc,iReg,iRegSrc,nPatches 
  
    REAL(KIND=RFREAL), DIMENSION(:,:), POINTER :: tav

    TYPE(t_patch),       POINTER :: pPatch,pPatchSrc
    TYPE(t_region),      POINTER :: pRegion,pRegionSrc
    TYPE(t_global),      POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(1)%global

    CALL RegisterFunction( global,'PLAG_RFLO_RecvStatBuffWrapper',&
  'PLAG_RFLO_ModStats.F90' )

! ******************************************************************************
!   Set dimensions, variables, pointers
! ******************************************************************************    

    pRegion => regions(iReg)
    
    iLev     = pRegion%currLevel
    nPatches = pRegion%nPatches

! ******************************************************************************
!   Receive data (regular cells) from other processors 
! ******************************************************************************

    DO iPatch=1,nPatches
      pPatch => pRegion%levels(iLev)%patches(iPatch)

      bcType    = pPatch%bcType
      iRegSrc   = pPatch%srcRegion
      iPatchSrc = pPatch%srcPatch

      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
          (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
        IF (regions(iRegSrc)%procid /= regions(iReg)%global%myProcid) THEN
          pRegionSrc => regions(iRegSrc)
          pPatchSrc  => pRegionSrc%levels(iLev)%patches(iPatchSrc)

          CALL PLAG_RFLO_RecvStatBuff( pRegion,pRegionSrc,pPatch,pPatchSrc )
        ENDIF ! procid
      
      ENDIF ! bcType 
    ENDDO   ! iPatch

! ******************************************************************************
!   Finalize
! ******************************************************************************

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_RecvStatBuffWrapper





!******************************************************************************
!
! Purpose: receive data for dummy cells from adjacent regions being
!          on another processor.
!
! Description: none.
!
! Input: pRegion     = current region
!        pRegionSrc  = source region
!        pPatch      = current patch
!        pPatchSrc   = source patch
!
! Output: pRegion%levels%plag%tav = updated statistics values
!                                   in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_RFLO_RecvStatBuff( pRegion,pRegionSrc,pPatch,pPatchSrc )

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch),       POINTER :: pPatch,pPatchSrc
    TYPE(t_region),      POINTER :: pRegion,pRegionSrc

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iLev,iPatch,iPatchSrc,iReg,iRegSrc,nEqs,nPatches
    INTEGER :: ijkBuff,ijkvBuff,ijkD,iv,ivs,lb,n1,n2,nDim
    INTEGER :: i,j,k,idum 
    INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, source, tag    
    INTEGER :: ivEvBeg,ivEvEnd,ivTavBeg,ivTavEnd,nEv,nTav
#ifdef MPI
    INTEGER :: status(MPI_STATUS_SIZE)
#endif

    REAL(KIND=RFREAL), DIMENSION(:,:), POINTER :: ev,tav   

    TYPE(t_buffer_plag), POINTER :: pBuffPlag 
    TYPE(t_global),      POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction( global,'PLAG_RFLO_RecvStatBuff',&
  'PLAG_RFLO_ModStats.F90' )

! ******************************************************************************
!   Check if the source region is active
! ******************************************************************************
    
    IF ( pRegionSrc%active == OFF ) THEN
      CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
    ENDIF

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    iLev = pRegion%currLevel
    
    ev    => pRegion%levels(iLev)%plag%ev
    tav   => pRegion%levels(iLev)%plag%tav
    pBuffPlag => pPatch%bufferPlag

! ******************************************************************************
!   Set dimensions 
! ******************************************************************************

    CALL RFLO_GetPatchIndices( pRegion,pPatch,iLev,ibeg,iend, &
                               jbeg,jend,kbeg,kend )
    CALL RFLO_GetCellOffset( pRegion,iLev,iCOff,ijCOff )

    n1   = ABS(pPatch%l1end-pPatch%l1beg) + 1   ! here, dimensions of current
    n2   = ABS(pPatch%l2end-pPatch%l2beg) + 1   ! and source patch are identical
    nDim = n1*n2*pRegion%nDumCells              ! ... but not the # of dummy cells

    nTav = global%plagNStat
    nEv  = pRegion%levels(iLev)%plag%nEv    
    nEqs = nEv +nTav       

    ivEvBeg = 1
    ivEvEnd = SIZE(ev,DIM=1)

    ivTavBeg = 1
    ivTavEnd = SIZE(tav,DIM=1)
 
! ******************************************************************************
!   Receive data
! ******************************************************************************

#ifdef MPI
    source = pRegionSrc%procid
    tag    = pRegion%localNumber + MPI_PATCHOFF*pPatchSrc%srcPatch
    IF(tag .gt. global%mpiTagMax) tag = MOD(tag,global%mpiTagMax)

    CALL MPI_Recv( pBuffPlag%recvBuffStat,nEqs*nDim,MPI_RFREAL,   &
                   source,tag,global%mpiComm,status,global%mpierr )
    IF (global%mpierr /= ERR_NONE) &
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! ******************************************************************************
!   Copy from buffer to dummy nodes 
! ******************************************************************************
    
    lb       = pPatch%lbound
    ijkBuff  = 0
    ijkvBuff = 0

    DO idum=1,pRegion%nDumCells

      SELECT CASE(lb)
      
! ==============================================================================
!       face i=const.
! ==============================================================================

        CASE(1:2)
          IF (lb == 1) i = ibeg - idum
          IF (lb == 2) i = iend + idum
          DO k=kbeg,kend
          DO j=jbeg,jend
            ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            ivs = 0
            DO iv=ivEvBeg,ivEvEnd
              ivs = ivs+1
              ijkvBuff = ijkBuff +(ivs-1)*nDim
              ev(iv,ijkD) = pBuffPlag%recvBuffStat(ijkvBuff) 
            ENDDO ! iv
            DO iv=ivTavBeg,ivTavEnd
              ivs = ivs+1
              ijkvBuff = ijkBuff +(ivs-1)*nDim
              tav(iv,ijkD) = pBuffPlag%recvBuffStat(ijkvBuff) 
            ENDDO ! iv
          ENDDO ! j
          ENDDO ! k

! ==============================================================================
!       face j=const.
! ==============================================================================

        CASE(3:4)
          IF (lb == 3) j = jbeg - idum
          IF (lb == 4) j = jend + idum
          DO i=ibeg,iend
          DO k=kbeg,kend
            ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            ivs = 0
            DO iv=ivEvBeg,ivEvEnd
              ivs = ivs+1
              ijkvBuff = ijkBuff +(ivs-1)*nDim
              ev(iv,ijkD) = pBuffPlag%recvBuffStat(ijkvBuff) 
            ENDDO ! iv
            DO iv=ivTavBeg,ivTavEnd
              ivs = ivs+1
              ijkvBuff = ijkBuff +(ivs-1)*nDim
              tav(iv,ijkD) = pBuffPlag%recvBuffStat(ijkvBuff) 
            ENDDO ! iv
          ENDDO ! k
          ENDDO ! i

! ==============================================================================
!       face k=const.
! ==============================================================================

        CASE(5:6)
          IF (lb == 5) k = kbeg - idum
          IF (lb == 6) k = kend + idum
          DO j=jbeg,jend
          DO i=ibeg,iend
            ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            ivs = 0
            DO iv=ivEvBeg,ivEvEnd
              ivs = ivs+1
              ijkvBuff = ijkBuff +(ivs-1)*nDim
              ev(iv,ijkD) = pBuffPlag%recvBuffStat(ijkvBuff) 
            ENDDO ! iv
            DO iv=ivTavBeg,ivTavEnd
              ivs = ivs+1
              ijkvBuff = ijkBuff +(ivs-1)*nDim
              tav(iv,ijkD) = pBuffPlag%recvBuffStat(ijkvBuff) 
            ENDDO ! iv
          ENDDO  ! i
          ENDDO  ! j
      END SELECT ! lb

    ENDDO ! idum

! ******************************************************************************
!   Finalize
! ******************************************************************************

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_RecvStatBuff






! ******************************************************************************
!
! Purpose: Wrapper to send buffer statistics arrays.
!
! Description: None.
!
! Input: regions = dimensions of all regions
!        iReg    = current region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLO_SendStatBuffWrapper( regions,iReg )

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: regions(:)

    INTEGER :: iReg

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iLev,iPatch,iPatchSrc,iRegSrc,nPatches

    REAL(RFREAL), POINTER, DIMENSION(:,:) :: tav
  
    TYPE(t_patch),       POINTER :: pPatch
    TYPE(t_region),      POINTER :: pRegion,pRegionSrc
    TYPE(t_global),      POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(iReg)%global

    CALL RegisterFunction( global,'PLAG_RFLO_SendStatBuffWrapper',&
  'PLAG_RFLO_ModStats.F90' )

! ******************************************************************************
!   Load communication buffers
! ******************************************************************************  

    pRegion => regions(iReg)
    
    IF ( pRegion%procid==global%myProcid .AND. &   ! region active and
         pRegion%active==ACTIVE ) THEN             ! on my processor

! ==============================================================================
!     Set dimensions and pointers
! ==============================================================================

      iLev     = pRegion%currLevel
      nPatches = pRegion%nPatches

! ==============================================================================         
!     Loop over patches 
! ==============================================================================
 
      DO iPatch=1,nPatches
        pPatch => pRegion%levels(iLev)%patches(iPatch)

        bcType    = pPatch%bcType
        iRegSrc   = pPatch%srcRegion
        iPatchSrc = pPatch%srcPatch

        SELECT CASE (bcType)

! ------------------------------------------------------------------------------
!         Conforming region interface
! ------------------------------------------------------------------------------

          CASE( BC_REGIONCONF:BC_REGIONCONF+BC_RANGE )
            pRegionSrc => regions(iRegSrc)

            IF ( pRegionSrc%procid /= global%myProcid ) THEN
              CALL PLAG_RFLO_SendStatBuffConf( pRegion,pRegionSrc,pPatch )
            ENDIF ! regions(iRegSrc)%procid

! ------------------------------------------------------------------------------
!         Non-conforming region interface (integer)
! ------------------------------------------------------------------------------

          CASE( BC_REGIONINT:BC_REGIONINT+BC_RANGE )
            CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )

! ------------------------------------------------------------------------------
!         Non-conforming region interface (irregular)
! ------------------------------------------------------------------------------
             
          CASE( BC_REGNONCONF:BC_REGNONCONF+BC_RANGE)
            CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )           
          
        END SELECT! bcType 

      END DO ! iPatch
    END IF ! pRegion%procid

! ******************************************************************************
!   Finalize
! ******************************************************************************

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_SendStatBuffWrapper






  

! ******************************************************************************
!
! Purpose: Send buffer arrays to dummy cells of the corresponding patch
!          on the adjacent region, residing on another processor.
!
! Description: None.
!
! Input: pRegion    = current region
!        pRegionSrc = source region
!        pPatch     = current patch
!
! Output: None.
!
! Notes: Pertinent for conforming regions.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLO_SendStatBuffConf( pRegion,pRegionSrc,pPatch )

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
 
    TYPE(t_patch),       POINTER :: pPatch
    TYPE(t_region),      POINTER :: pRegion,pRegionSrc

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: iv,ivs,iRequestStat,nPatches,nEqs
    INTEGER :: bcType,iPatch,iPatchSrc,iRegSrc
    INTEGER :: idum, i, j, k, ijkBuff, ijkvBuff
    INTEGER :: iLev, ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkC, &
               n1, n2, nDim, dest, tag
    INTEGER :: lb, l1SrcDir, l2SrcDir, l1Beg, l1End, l1Step, l2Beg, l2End, l2Step
    INTEGER :: ivEvBeg,ivEvEnd,ivTavBeg,ivTavEnd,nEv,nTav

    LOGICAL :: align
    
    REAL(KIND=RFREAL), DIMENSION(:,:), POINTER :: ev,tav   
    
    TYPE(t_buffer_plag), POINTER :: pBuffPlag    
    TYPE(t_plag),        POINTER :: pPlag
    TYPE(t_global),      POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction( global,'PLAG_RFLO_SendStatBuffConf',&
  'PLAG_RFLO_ModStats.F90' )

! ******************************************************************************
!   Check if the source region is active
! ******************************************************************************
    
    IF ( pRegionSrc%active == OFF ) THEN
      CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
    ENDIF

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    iLev = pRegion%currLevel
    
    pPlag     => pRegion%levels(iLev)%plag
    pBuffPlag => pPatch%bufferPlag
    ev        => pPlag%ev
    tav       => pPlag%tav    

! ******************************************************************************
!   Set dimensions and pointers
! ******************************************************************************
         
    CALL RFLO_GetPatchIndices( pRegion,pPatch,iLev,ibeg,iend, &
                               jbeg,jend,kbeg,kend )
    CALL RFLO_GetCellOffset( pRegion,iLev,iCOff,ijCOff )

    n1   = ABS(pPatch%l1end-pPatch%l1beg) + 1   ! here, dimensions of current
    n2   = ABS(pPatch%l2end-pPatch%l2beg) + 1   ! and source patch are identical
    nDim = n1*n2*pRegionSrc%nDumCells           ! ... but not the # of dummy cells

    nTav = global%plagNStat
    nEv  = pPlag%nEv    
    nEqs = nEv +nTav   
    
    ivEvBeg = 1
    ivEvEnd = SIZE(ev,DIM=1)
   
    ivTavBeg = 1
    ivTavEnd = SIZE(tav,DIM=1)
 
    iRequestStat = pBuffPlag%iRequestStat    

! ******************************************************************************
!   Mapping between patches
! ******************************************************************************

                                           l1SrcDir =  1
    IF (pPatch%srcL1beg > pPatch%srcL1end) l1SrcDir = -1
                                           l2SrcDir =  1
    IF (pPatch%srcL2beg > pPatch%srcL2end) l2SrcDir = -1

    lb    = pPatch%lbound
    align = pPatch%align

! ******************************************************************************
!   Loop over interior cells of current patch loading buffers 
! ******************************************************************************
    
    ijkBuff  = 0
    ijkvBuff = 0

    DO idum=0,pRegionSrc%nDumCells-1

      SELECT CASE(lb)
      
! ==============================================================================
!       face i=const.
! ==============================================================================

        CASE(1:2)
          IF (lb == 1) i = ibeg + idum
          IF (lb == 2) i = iend - idum
          
          IF (align) THEN
            IF (l1SrcDir > 0) THEN
              l1Beg = jbeg
              l1End = jend
            ELSE
              l1Beg = jend
              l1End = jbeg
            ENDIF ! l1SrcDir
            l1Step = l1SrcDir
            
            IF (l2SrcDir > 0) THEN
              l2Beg = kbeg
              l2End = kend
            ELSE
              l2Beg = kend
              l2End = kbeg
            ENDIF ! l2SrcDir
            l2Step = l2SrcDir

            DO k=l2Beg,l2End,l2Step
            DO j=l1Beg,l1End,l1Step
              ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
              ijkBuff = ijkBuff + 1
              ivs = 0 
              DO iv=ivEvBeg,ivEvEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = ev(iv,ijkC)
              ENDDO ! iv
              DO iv=ivTavBeg,ivTavEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = tav(iv,ijkC)
              ENDDO ! iv
            ENDDO   ! j
            ENDDO   ! k
          ELSE
            IF (l1SrcDir > 0) THEN
              l1Beg = kbeg
              l1End = kend
            ELSE
              l1Beg = kend
              l1End = kbeg
            ENDIF ! l1SrcDir
            l1Step = l1SrcDir
        
            IF (l2SrcDir > 0) THEN
              l2Beg = jbeg
              l2End = jend
            ELSE
              l2Beg = jend
              l2End = jbeg
            ENDIF ! l2SrcDir
            l2Step = l2SrcDir
       
            DO j=l2Beg,l2End,l2Step
            DO k=l1Beg,l1End,l1Step
              ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
              ijkBuff = ijkBuff + 1
              ivs = 0
              DO iv=ivEvBeg,ivEvEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = ev(iv,ijkC)
              ENDDO ! iv
              DO iv=ivTavBeg,ivTavEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = tav(iv,ijkC)
              ENDDO ! iv
            ENDDO   ! j
            ENDDO   ! k
          ENDIF ! align
      
! ==============================================================================
!       face j=const.
! ==============================================================================

        CASE(3:4)
          IF (lb == 3) j = jbeg + idum
          IF (lb == 4) j = jend - idum
          IF (align) THEN
            IF (l1SrcDir > 0) THEN
              l1Beg = kbeg
              l1End = kend
            ELSE
              l1Beg = kend
              l1End = kbeg
            ENDIF ! l1SrcDir
            l1Step = l1SrcDir
            
            IF (l2SrcDir > 0) THEN
              l2Beg = ibeg
              l2End = iend
            ELSE
              l2Beg = iend
              l2End = ibeg
            ENDIF ! l2SrcDir
            l2Step = l2SrcDir

            DO i=l2Beg,l2End,l2Step
            DO k=l1Beg,l1End,l1Step
              ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
              ijkBuff = ijkBuff + 1
              ivs = 0
              DO iv=ivEvBeg,ivEvEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = ev(iv,ijkC)
              ENDDO ! iv
              DO iv=ivTavBeg,ivTavEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = tav(iv,ijkC)
              ENDDO ! iv             
            ENDDO ! k
            ENDDO ! i
          ELSE
            IF (l1SrcDir > 0) THEN
              l1Beg = ibeg
              l1End = iend
            ELSE
              l1Beg = iend
              l1End = ibeg
            ENDIF ! l1SrcDir 
            l1Step = l1SrcDir
            
            IF (l2SrcDir > 0) THEN
              l2Beg = kbeg
              l2End = kend
            ELSE
              l2Beg = kend
              l2End = kbeg
            ENDIF ! l2SrcDir 
            l2Step = l2SrcDir

            DO k=l2Beg,l2End,l2Step
            DO i=l1Beg,l1End,l1Step
              ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
              ijkBuff = ijkBuff + 1
              ivs = 0
              DO iv=ivEvBeg,ivEvEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = ev(iv,ijkC)
              ENDDO ! iv
              DO iv=ivTavBeg,ivTavEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = tav(iv,ijkC)
              ENDDO ! iv
            ENDDO ! i
            ENDDO ! k
          ENDIF ! align

! ==============================================================================
!       face k=const.
! ==============================================================================

        CASE(5:6)
          IF (lb == 5) k = kbeg + idum
          IF (lb == 6) k = kend - idum
          IF (align) THEN
            IF (l1SrcDir > 0) THEN
              l1Beg = ibeg
              l1End = iend
            ELSE
              l1Beg = iend
              l1End = ibeg
            ENDIF ! l1SrcDir 
            l1Step = l1SrcDir
            
            IF (l2SrcDir > 0) THEN
              l2Beg = jbeg
              l2End = jend
            ELSE
              l2Beg = jend
              l2End = jbeg
            ENDIF ! l2SrcDir
            l2Step = l2SrcDir

            DO j=l2Beg,l2End,l2Step
            DO i=l1Beg,l1End,l1Step
              ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
              ijkBuff = ijkBuff + 1
              ivs = 0
              DO iv=ivEvBeg,ivEvEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = ev(iv,ijkC)
              ENDDO ! iv
              DO iv=ivTavBeg,ivTavEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = tav(iv,ijkC)
              ENDDO ! iv
            ENDDO ! i
            ENDDO ! j
          ELSE
            IF (l1SrcDir > 0) THEN
              l1Beg = jbeg
              l1End = jend
            ELSE
              l1Beg = jend
              l1End = jbeg
            ENDIF ! l1SrcDir
            l1Step = l1SrcDir
            
            IF (l2SrcDir > 0) THEN
              l2Beg = ibeg
              l2End = iend
            ELSE
              l2Beg = iend
              l2End = ibeg
            ENDIF ! l2SrcDir
            l2Step = l2SrcDir
           
            DO i=l2Beg,l2End,l2Step
            DO j=l1Beg,l1End,l1Step
              ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
              ijkBuff = ijkBuff + 1
              ivs = 0
              DO iv=ivEvBeg,ivEvEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = ev(iv,ijkC)
              ENDDO ! iv
              DO iv=ivTavBeg,ivTavEnd
                ivs = ivs +1
                ijkvBuff = ijkBuff +(ivs-1)*nDim
                pBuffPlag%sendBuffStat(ijkvBuff) = tav(iv,ijkC)
              ENDDO ! iv
            ENDDO ! j
            ENDDO ! i
          ENDIF ! align
      END SELECT ! lb

    ENDDO ! idum

! ******************************************************************************
!   Send data
! ******************************************************************************

#ifdef MPI
    dest = pRegionSrc%procid
    tag  = pRegionSrc%localNumber + MPI_PATCHOFF*pPatch%srcPatch
    IF(tag .gt. global%mpiTagMax) tag = MOD(tag,global%mpiTagMax)

    CALL MPI_Isend( pBuffPlag%sendBuffStat,nEqs*nDim,MPI_RFREAL, &
                    dest,tag,global%mpiComm, &
                    pPlag%requestsStat(iRequestStat),global%mpierr )
  IF (global%mpierr /= ERR_NONE) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! ******************************************************************************
!   Finalize
! ******************************************************************************

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_SendStatBuffConf






  



!******************************************************************************
!
! Purpose: Write time averaged statistics of the Lagrangian particles 
!          on Eulerian grid.
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary
!
! Input: regions = dimensions of all regions
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary
!
! Input: regions        = dimensions and cons. variables of all regions
!
! Output: to file
!
! Notes: solution is stored only for the current grid level; it is also
!        stored for all dummy cells; all regions are written into one file
!
!******************************************************************************

  SUBROUTINE PLAG_RFLO_WriteStat(regions)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: iReg, i, j, k, l, n

     CHARACTER(2*CHRLEN+17) :: fname

#ifdef MPI
    INTEGER :: status(MPI_STATUS_SIZE)
#endif
    INTEGER :: iLev, iRegFile, ipc, jpc, kpc, nDumCells, nDim, iOff, ijOff, ijk
    INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkBeg, ijkEnd
    INTEGER :: errorFlag, nTav
    INTEGER, ALLOCATABLE , DIMENSION(:)   :: plagVarId
    INTEGER, ALLOCATABLE , DIMENSION(:,:) :: ivar,jvar

    REAL(RFREAL), POINTER,     DIMENSION(:,:) :: tav
    REAL(RFREAL), ALLOCATABLE, DIMENSION(:,:) :: rvar, tavFile

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(1)%global

    CALL RegisterFunction( global,'PLAG_RFLO_WriteStat',&
  'PLAG_RFLO_ModStats.F90' )

    IF (.NOT. global%plagUsed) GOTO 999

! ******************************************************************************
!   Allocate temporary data arrays 
! ******************************************************************************

    ALLOCATE( ivar(5,1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'ivar' )
    
    ALLOCATE( rvar(2,1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'rvar' )
    
    ALLOCATE( jvar(global%plagNStat+1,1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'rvar' )
    
    ALLOCATE( plagVarId(global%plagNStat+1),stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'plagVarId' ) 

! ******************************************************************************
!   Write verbosity (only master proc.)
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing PLAG statistics solution file...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Open statistics file (only master proc.)
! ******************************************************************************

    IF (global%myProcid == MASTERPROC) THEN

      SELECT CASE( global%solutFormat )
        CASE ( FORMAT_ASCII )
          WRITE(fname,'(A,1PE11.5)') &
            TRIM(global%inDir)//TRIM(global%casename)//'.plag_stata_', &
            global%currentTime
          OPEN( IF_PLAG_STATS,file=fname,form='formatted',status='unknown', &
                iostat=errorFlag )
      
        CASE ( FORMAT_BINARY )
          WRITE(fname,'(A,1PE11.5)') &
            TRIM(global%inDir)//TRIM(global%casename)//'.plag_stat_', &
            global%currentTime
          OPEN( IF_PLAG_STATS,file=fname,form='unformatted',status='unknown', &
                iostat=errorFlag )
        
        CASE DEFAULT
          CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

      END SELECT ! solutFormat 

      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! ******************************************************************************
!     Write current and integrated time to file 
! ******************************************************************************

      rvar(1,1) = global%currentTime
      rvar(2,1) = global%integrTime

      CALL RFLO_WriteDataFileReal( global,IF_PLAG_STATS,global%solutFormat,&
                                   2,1,rvar )

! ******************************************************************************
!     Write plagNStat and plagStatId to file
! ******************************************************************************
      
      IF ( global%plagNStat > 0 ) THEN
        jvar(1,1)    = global%plagNStat
        plagVarId(:) = global%plagStatId(1,:)*10 + global%plagStatId(2,:)
        jvar(2:global%plagNStat+1,1) = plagVarId(1:global%plagNStat)
        CALL RFLO_WriteDataFileInt( global,IF_PLAG_STATS,global%solutFormat, &
                                    global%plagNStat+1,1,jvar )
      ENDIF ! plagNStat

    ENDIF ! global%myProcid

! ******************************************************************************
!   Write statistics data 
! ******************************************************************************

    DO iReg=1,global%nRegions

! ==============================================================================
!     Get dimensions and pointers
! ==============================================================================

      iLev = regions(iReg)%currLevel
      CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                                jdcbeg,jdcend,kdcbeg,kdcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iOff,ijOff )
      ijkBeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
      ijkEnd = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
      nDim   = ijkEnd - ijkBeg + 1

      nTav = global%plagNStat

! ==============================================================================
!    Allocate memory for data field
! ==============================================================================

      IF ( regions(iReg)%procid==global%myProcid .OR. &
           global%myProcid==MASTERPROC                ) THEN
        ALLOCATE( tavFile(nTav,nDim),stat=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'tavFile' )
      ENDIF ! regions(iReg)%procid

! ==============================================================================
!     Copy statistics into data structure
! ==============================================================================

      IF ( regions(iReg)%procid == global%myProcid ) THEN
        tav => regions(iReg)%levels(iLev)%plag%tav

        n = 0
        DO k=kdcbeg,kdcend
        DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          n   = n + 1
          ijk = IndIJK(i,j,k,iOff,ijOff)
          DO l=1,nTav
            tavFile(l,n) = tav(l,ijk)
          ENDDO ! l
        ENDDO ! i
        ENDDO ! j
        ENDDO ! k
      ENDIF ! regions(iReg)%procid 

! ==============================================================================
!     Write region number and dimensions (only master)
! ==============================================================================

      IF ( global%myProcid == MASTERPROC ) THEN
        ivar(1,1) = iReg
        ivar(2,1) = regions(iReg)%levels(iLev)%grid%ipc
        ivar(3,1) = regions(iReg)%levels(iLev)%grid%jpc
        ivar(4,1) = regions(iReg)%levels(ilev)%grid%kpc
        ivar(5,1) = regions(iReg)%nDumCells
        CALL RFLO_WriteDataFileInt( global,IF_PLAG_STATS,global%solutFormat,5,1,ivar )
      ENDIF ! global%myProcid

! ==============================================================================
!     Master receives and writes data, others send them
! ==============================================================================

      IF (global%myProcid == MASTERPROC) THEN
#ifdef MPI
        IF ( regions(iReg)%procid /= MASTERPROC ) THEN
          CALL MPI_Recv( tavFile,nTav*nDim,MPI_RFREAL, &
                         regions(iReg)%procid,iReg, &
                         global%mpiComm,status,global%mpierr )
          IF ( global%mpierr /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        ENDIF ! regions(iReg)%procid
#endif
        CALL RFLO_WriteDataFileReal( global,IF_PLAG_STATS,global%solutFormat, &
                                     nTav,nDim,tavFile )

      ELSE   ! not the master
#ifdef MPI
        IF (regions(iReg)%procid == global%myProcid) THEN
          CALL MPI_Send( tavFile,nTav*nDim,MPI_RFREAL,MASTERPROC, &
                         iReg,global%mpiComm,global%mpierr )
          IF ( global%mpierr /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
#endif
      ENDIF  ! global%myProcid

! ==============================================================================
!     Deallocate local array
! ==============================================================================

      IF (ALLOCATED(tavFile)) THEN
        DEALLOCATE( tavFile,stat=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'tavFile' )
      ENDIF

    ENDDO     ! iReg

! ******************************************************************************
!   Deallocate temporary data arrays 
! ******************************************************************************

    DEALLOCATE( ivar,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'ivar' )
    
    DEALLOCATE( rvar,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'rvar' )
   
    DEALLOCATE( jvar,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'jvar' )
   
    DEALLOCATE( plagVarId,stat=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'plagVarId' )

! ******************************************************************************
!   Write verbosity (only master proc.)
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing PLAG statistics solution file done...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Finalize
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC ) THEN
      CLOSE(IF_PLAG_STATS,iostat=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
    ENDIF

! ******************************************************************************
!   End  
! ******************************************************************************

999  CONTINUE
    CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_RFLO_WriteStat







END MODULE PLAG_RFLO_ModStats

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLO_ModStats.F90,v $
! Revision 1.9  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.8  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.5  2005/06/19 07:13:57  wasistho
! refrained from communicating and writing statistics when plag is not used
!
! Revision 1.4  2005/03/08 00:59:00  fnajjar
! Bug fix for ijkvBuff with multiple assignments
!
! Revision 1.3  2005/03/07 17:37:37  fnajjar
! Bug fix in loading iRequestStat
!
! Revision 1.2  2005/02/16 14:47:41  fnajjar
! Added infrastructure for on-processor copying and MPI-based communication
!
! Revision 1.1  2005/01/08 20:44:32  fnajjar
! Initial import for PLAG statistics
!
! ******************************************************************************
















