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
! Purpose: write Tecplot file for Lagragian statistics particle field 
!          onto Eulerian-based grid.
!
! Description: none.
!
! Input: regions  = pointer to all regions
!        iReg     = region number
!
! Output: to file.
!
! Notes: 
!   1.   The output is collected in one file, but the regions are processed
!        separately to save memory.
!
! ******************************************************************************
!
! $Id: PLAG_WriteStatTecAscii.F90,v 1.3 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_WriteStatTecAscii( regions, iReg)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModGrid,       ONLY : t_grid
  USE ModPartLag,    ONLY : t_plag
  USE PLAG_ModInterfacesPost, ONLY : RFLO_CopyGeometryDummy,   &
                                     RFLO_GenerateCoarseGrids, &
				     RFLO_GetCellOffset,       &
				     RFLO_GetDimensPhysNodes,  &
                                     RFLO_GetDimensDummy,      &
                                     RFLO_GetDimensDummyNodes, &
				     RFLO_GetNodeOffset,       &
				     RFLO_GetCellOffset,       &
				     RFLO_GetDimensPhys,       &
                                     RFLO_ReadGridRegion, Aver

  USE ModParameters
  USE PLAG_ModParameters
  
  IMPLICIT NONE

#include "Indexing.h"

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN) :: iReg

! ==============================================================================  
! Locals
! ==============================================================================   

  CHARACTER(CHRLEN)   :: RCSIdentString
  CHARACTER(CHRLEN+4) :: fnameTec
  CHARACTER(256)      :: varStrTec
  CHARACTER(16)       :: compStr
  
  INTEGER, PARAMETER :: IF_PLOT_PLAGSTAT_TEC = IF_PLOT +210

  INTEGER :: iLev, ipc, jpc, kpc, ibc, iec, ibn, ien
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: errorFlag
  INTEGER :: i,j,k, ijkN
  INTEGER :: cell(8),node(8)
  INTEGER :: ic,iCont,itavCont
  INTEGER :: ijkCell,ijkDum
  INTEGER :: nCont,nEv,nTav
  INTEGER, POINTER,      DIMENSION(:)   :: cvMass

  REAL(RFREAL) :: currentTime,rTimeStat
  REAL(RFREAL) :: snumDens,sdiamL,sdiamL3,sdiamL4,smass
  REAL(RFREAL) :: xcell,ycell,zcell
  REAL(RFREAL) :: su,sv,sw,suu,svv,sww
  REAL(RFREAL), DIMENSION(:),   ALLOCATABLE :: scomp
  REAL(RFREAL), DIMENSION(:,:), POINTER :: tav
  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid)  , POINTER :: grid
  TYPE(t_plag)  , POINTER :: pPlag
  
!*******************************************************************************

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_WriteStatTecAscii.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_WriteStatTecAscii',&
  'PLAG_WriteStatTecAscii.F90' )

  IF (.NOT. global%plagUsed) GOTO 999

! ******************************************************************************
! Obtain region grid size
! ******************************************************************************

  DO iLev=2,regions(iReg)%nGridLevels
    ipc = regions(iReg)%levels(iLev-1)%grid%ipc
    jpc = regions(iReg)%levels(iLev-1)%grid%jpc
    kpc = regions(iReg)%levels(iLev-1)%grid%kpc
    regions(iReg)%levels(iLev)%grid%ipc = ipc/2
    regions(iReg)%levels(iLev)%grid%jpc = jpc/2
    regions(iReg)%levels(iLev)%grid%kpc = kpc/2
  ENDDO ! iLev

! ******************************************************************************
! Allocate memory for the grid (all grid levels)
! ******************************************************************************

  DO iLev=1,regions(iReg)%nGridLevels  
    CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
    ibn  =  IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien  =  IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

    grid => regions(iReg)%levels(iLev)%grid
    ALLOCATE( grid%xyz(3,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ENDDO ! iLev

! ******************************************************************************
! Set current grid level and get dummy cell dimensions
! ******************************************************************************

  iLev =  regions(iReg)%currLevel
  grid => regions(iReg)%levels(iLev)%grid

  CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
  ibn = IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff)
  ien = IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff)

  CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

! ******************************************************************************
! Save currentTime as global value is clobbered by RFLO_ReadGridRegion
! ******************************************************************************
  
  global%currentTime = global%timeStamp
  currentTime = global%currentTime 

! ******************************************************************************
! Read grid
! ******************************************************************************

  CALL RFLO_ReadGridRegion( iReg,regions )
  CALL RFLO_GenerateCoarseGrids( regions(iReg) )
  CALL RFLO_CopyGeometryDummy( regions(iReg) )

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
 
  pPlag => regions(iReg)%levels(iLev)%plag

  tav  => pPlag%tav

  cvMass => pPlag%cvPlagMass
  
  nCont = regions(iReg)%plagInput%nCont
  nTav  = global%plagNStat
  itavCont = 7

  rTimeStat = 1.0_RFREAL/global%integrTime

! ******************************************************************************
! Allocate local memory
! ******************************************************************************
  
  ALLOCATE( scomp(nCont),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )   

! ******************************************************************************
! Write data to Tecplot file
! ******************************************************************************

! ==============================================================================
!  Open file and write the header
! ==============================================================================

   IF (iReg == 1) THEN
     WRITE(fnameTec,'(A,ES11.5,A)') &
     TRIM(global%casename)//'.plag_stats_',currentTime,'.dat'
     OPEN(IF_PLOT_PLAGSTAT_TEC,FILE=fnameTec,STATUS='unknown',FORM='formatted',&
	  iostat=errorFlag)
     global%error = errorFlag
     IF (global%error /= 0) &
       CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fnameTec )     
     
     WRITE(IF_PLOT_PLAGSTAT_TEC,1005,err=10) &
       TRIM(global%casename),currentTime

     varStrTec = 'x y z <nDens> <diam> <u> <v> <w> <mass> '

     DO iCont=1,nCont
       SELECT CASE (iCont)
         CASE ( 0: 9)
          WRITE(compStr,'(A,I1,A)') '<comp_', iCont,'>'
         CASE (10:99)
          WRITE(compStr,'(A,I2,A)') 'comp_', iCont,'>'
         CASE DEFAULT
          WRITE(compStr,'(A)') 'comp_?>'
       END SELECT ! iCont
       varStrTec = TRIM(varStrTec)//' '//TRIM(compStr)
     ENDDO ! iCont
     
     WRITE(compStr,'(A)') '<uu>'
     varStrTec = TRIM(varStrTec)//' '//TRIM(compStr)
     
     WRITE(IF_PLOT_PLAGSTAT_TEC,1010,err=10) TRIM(varStrTec)
   ENDIF   ! iReg=1

! ==============================================================================
!  Write zone header for plotting data based on nodes
! ==============================================================================

   WRITE(IF_PLOT_PLAGSTAT_TEC,1015) iReg,ipnend-ipnbeg+1,jpnend-jpnbeg+1,&
                                    kpnend-kpnbeg+1

! ==============================================================================
! Write plotting data
! ==============================================================================

  DO k=kpnbeg,kpnend
  DO j=jpnbeg,jpnend
  DO i=ipnbeg,ipnend
    ijkN = IndIJK(i,j,k,iNOff,ijNOff)
    
    cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
    cell(2) = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
    cell(3) = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
    cell(4) = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
    cell(5) = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
    cell(6) = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
    cell(7) = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
    cell(8) = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
    
    sdiamL3  = Aver(cell,1,tav) *rTimeStat ! 01 plagStatId
    sdiamL4  = Aver(cell,2,tav) *rTimeStat ! 02 
    snumDens = Aver(cell,3,tav) *rTimeStat ! 03  
    su       = Aver(cell,4,tav) *rTimeStat ! 04
    sv       = Aver(cell,5,tav) *rTimeStat ! 05
    sw       = Aver(cell,6,tav) *rTimeStat ! 06
    smass    = Aver(cell,7,tav) *rTimeStat ! 07

    DO iCont = 1, nCont
      sComp(iCont)   = Aver(cell,itavCont+iCont,tav) *rTimeStat
      IF ( smass > 0.0_RFREAL ) THEN
        sComp(iCont) = sComp(iCont)/smass
      ELSE
        sComp(iCont) = 0.0_RFREAL
      END IF ! smass
    END DO ! iCont
    
    suu = Aver(cell,10,tav) *rTimeStat    ! 44
    suu = suu -su*su

    IF( sdiamL3 > 0.0_RFREAL) THEN
      sdiamL = sdiamL4/sdiamL3
    ELSE
      sdiamL = 0.0_RFREAL
    ENDIF ! diamL3

    WRITE(IF_PLOT_PLAGSTAT_TEC,1020,err=10) grid%xyz(XCOORD,ijkN),            &
                                            grid%xyz(YCOORD,ijkN),            &
                                            grid%xyz(ZCOORD,ijkN),            &
                                            snumDens, sdiamL,                 &
                                            su,sv,sw,smass,                   &
                                           (sComp(iCont),iCont=1,nCont),      &
                                            suu

  ENDDO ! i
  ENDDO ! j
  ENDDO ! k

  PRINT*,' PLAG_WriteStatTec: iReg, nDensSum = ',iReg,SUM(tav(3,:))*rTimeStat

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DO iLev=1,regions(iReg)%nGridLevels
    grid => regions(iReg)%levels(iLev)%grid
    DEALLOCATE( grid%xyz,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  ENDDO ! iLev

! ******************************************************************************
! Allocate local memory
! ******************************************************************************
  
  DEALLOCATE( scomp,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )   

! ******************************************************************************
! Finalize
! ******************************************************************************
 
  IF (iReg == global%nRegions) THEN
    CLOSE(IF_PLOT_PLAGSTAT_TEC,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fnameTec )
  ENDIF

  GOTO 999

10   CONTINUE

999  CONTINUE
  CALL DeregisterFunction( global )

! ******************************************************************************
! Formats
! ******************************************************************************

1005 FORMAT('TITLE="',A,'. Plag Statistics-Time: ',ES11.5,'."')
1010 FORMAT('VARIABLES= ',A)
1015 FORMAT('ZONE T="',I5.5,'", I=',I6,', J=',I6,', K=',I6,', F=POINT')
1020 FORMAT(999(1X,ES13.6))

END SUBROUTINE PLAG_WriteStatTecAscii

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_WriteStatTecAscii.F90,v $
! Revision 1.3  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/02/16 14:52:40  fnajjar
! Initial import
!
!******************************************************************************







