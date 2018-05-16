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
!
! $Id: PLAG_ReadStatPost.F90,v 1.3 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_ReadStatPost(regions,iReg)

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
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: regions(:)  
    INTEGER, INTENT(IN) :: iReg
 
! ==============================================================================
!   Locals
! ==============================================================================

  CHARACTER(2*CHRLEN+17) :: fname
  CHARACTER(CHRLEN)      :: msg

  INTEGER :: i, j, k, l, n, ind

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

    global => regions(iReg)%global

    CALL RegisterFunction( global,'PLAG_RFLO_ReadStatPost',&
  'PLAG_ReadStatPost.F90' )

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
!   Open statistics file (only if iReg=1.)
! ******************************************************************************

    IF (iReg == 1) THEN

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

    ENDIF ! iReg

! ******************************************************************************
!   Read current and integrated time in file, and stats ID (only if iReg=1)
! ******************************************************************************

    IF ( iReg == 1 ) THEN
      CALL RFLO_ReadDataFileReal( global,IF_PLAG_STATS,global%solutFormat,2,1,rvar )
      
       global%integrTime = rvar(2,1) 

! ==============================================================================
!   Trap error for inconsistent variables in header 
! ==============================================================================
    
      IF ( global%flowType==FLOW_UNSTEADY .AND. global%currentTime>0._RFREAL ) THEN
        IF (ABS(global%currentTime-rvar(1,1))/global%currentTime > 1.0E-03_RFREAL) THEN
          WRITE(msg,1000) rvar(1,1),global%currentTime
          CALL ErrorStop( global,ERR_TIME_SOLUTION,__LINE__,msg//' File: '//TRIM(fname) )
        ENDIF ! currentTime
      ENDIF ! global%flowType

    PRINT*,' PLAG_ReadStatPost: integrTime, currentime = ',global%integrTime,global%currentTime,rvar(1,1)

    ENDIF ! iReg

! ******************************************************************************
!   Read plagNStat and plagStatId from file
! ******************************************************************************
    
    IF ( iReg == 1 ) THEN
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

    ENDIF ! iReg
      
! ******************************************************************************
!   Read statistics data from all regions 
! ******************************************************************************

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
!     Read region number and dimensions 
! ==============================================================================

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
!     read data
! ==============================================================================

      ALLOCATE( tavFile(nTav,nDim),stat=errorFlag )
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL RFLO_ReadDataFileReal( global,IF_PLAG_STATS,global%solutFormat, &
                                  nTav,nDim,tavFile )

! ==============================================================================
!    Copy statistics into data structure
! ==============================================================================
     
      IF ( nTav > 0 ) tav => regions(iReg)%levels(iLev)%plag%tav

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

! ==============================================================================
!     Deallocate local array
! ==============================================================================

      IF ( ALLOCATED(tavFile) ) THEN
        DEALLOCATE( tavFile,stat=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'tavFile' )
      ENDIF ! tavFile

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

    IF ( iReg == global%nRegions ) THEN
      CLOSE(IF_PLAG_STATS,iostat=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
    ENDIF ! iReg

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

END SUBROUTINE PLAG_ReadStatPost

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadStatPost.F90,v $
! Revision 1.3  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/02/16 14:52:40  fnajjar
! Initial import
!
! ******************************************************************************







