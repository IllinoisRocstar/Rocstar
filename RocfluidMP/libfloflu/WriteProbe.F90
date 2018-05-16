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
! Purpose: write data of a probe into a file.
!
! Description: none.
!
! Input: regions%levels%mixt = flow variables
!        iReg                = current region number
!        global%probePos     = list of probes
!
! Output: into file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: WriteProbe.F90,v 1.5 2009/01/06 21:29:46 mdbrandy Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteProbe( regions,iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetDimensPhys

#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iprobe

! ... local variables
  CHARACTER(CHRLEN+9) :: fname

  INTEGER :: errorFlag,ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, iCell, nPeul, i, j, k

  LOGICAL :: wrtProbe

  REAL(RFREAL)          :: rho, u, v, w, press, temp
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), peulCv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'WriteProbe',&
  'WriteProbe.F90' )

! determine number of smoke/species types that exist (if any)

  nPeul = 0

#ifdef RFLO
#ifdef PEUL
  IF (global%peulUsed) nPeul = regions(iReg)%peulInput%nPtypes
#endif
#endif

#ifdef RFLU
#ifdef SPEC
  IF (global%specUsed) nPeul = regions(iReg)%specInput%nSpecies
#endif
#endif

! loop over all specified probes ----------------------------------------------

  DO iprobe=1,global%nProbes

    wrtProbe = .false.

#ifdef RFLO
! - check if region number within range

    IF (global%probePos(iprobe,1)<1 .OR. &
        global%probePos(iprobe,1)>global%nRegions) &
      CALL ErrorStop( global,ERR_PROBE_LOCATION,__LINE__ )

! - prepare data

    IF (regions(global%probePos(iprobe,1))%procid==global%myProcid .AND. &
        regions(global%probePos(iprobe,1))%active==ACTIVE .AND. &
        iReg==global%probePos(iprobe,1)) THEN

! --- get dimensions

      iLev = regions(iReg)%currLevel

      CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                               jpcbeg,jpcend,kpcbeg,kpcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

! --- check if probe within index range; get cell pointer

      IF ((global%probePos(iprobe,2)<ipcbeg  .OR. &
           global%probePos(iprobe,2)>ipcend) .OR. &
          (global%probePos(iprobe,3)<jpcbeg  .OR. &
           global%probePos(iprobe,3)>jpcend) .OR. &
          (global%probePos(iprobe,4)<kpcbeg  .OR. &
           global%probePos(iprobe,4)>kpcend)) &
        CALL ErrorStop( global,ERR_PROBE_LOCATION,__LINE__ )

      i     = global%probePos(iprobe,2)
      j     = global%probePos(iprobe,3)
      k     = global%probePos(iprobe,4)
      iCell = IndIJK(i,j,k,iCOff,ijCOff)

      cv => regions(iReg)%levels(iLev)%mixt%cv
      dv => regions(iReg)%levels(iLev)%mixt%dv
#ifdef PEUL
      IF (nPeul > 0) peulCv => regions(iReg)%levels(iLev)%peul%cv
#endif
      wrtProbe = .true.
    ENDIF
#endif
#ifdef RFLU
    IF ( global%probePos(iprobe,PROBE_REGION) == &
         regions(iReg)%iRegionGlobal ) THEN

      iCell = global%probePos(iprobe,PROBE_CELL)

      cv => regions(iReg)%mixt%cv
      dv => regions(iReg)%mixt%dv
      IF (nPeul > 0) peulCv => regions(iReg)%spec%cv

      wrtProbe = .TRUE.
    END IF ! global%probePos
#endif

! - write probe data to file

    IF (wrtProbe) THEN
      rho   = cv(CV_MIXT_DENS,iCell)
      u     = cv(CV_MIXT_XMOM,iCell)/rho
      v     = cv(CV_MIXT_YMOM,iCell)/rho
      w     = cv(CV_MIXT_ZMOM,iCell)/rho
      press = dv(DV_MIXT_PRES,iCell)
      temp  = dv(DV_MIXT_TEMP,iCell)

      IF (nPeul == 0) THEN
        IF (global%flowType == FLOW_STEADY) THEN
          WRITE(IF_PROBE+iprobe-1,1000,IOSTAT=errorFlag) global%currentIter, &
                                                         rho,u,v,w,press,temp
        ELSE
          WRITE(IF_PROBE+iprobe-1,1005,IOSTAT=errorFlag) global%currentTime, &
                                                         rho,u,v,w,press,temp
        ENDIF
      ELSE
#ifdef RFLO
#ifdef PEUL
        IF (global%flowType == FLOW_STEADY) THEN
          WRITE(IF_PROBE+iprobe-1,1000,IOSTAT=errorFlag) global%currentIter,  &
                                                         rho,u,v,w,press,temp,&
                                                         peulCv(1:nPeul,iCell)
        ELSE
          WRITE(IF_PROBE+iprobe-1,1005,IOSTAT=errorFlag) global%currentTime,  &
                                                         rho,u,v,w,press,temp,&
                                                         peulCv(1:nPeul,iCell)
        ENDIF
#endif
#endif
#ifdef RFLU
        IF (global%flowType == FLOW_STEADY) THEN
          WRITE(IF_PROBE+iprobe-1,1000,IOSTAT=errorFlag) global%currentIter,  &
                                                         rho,u,v,w,press,temp,&
                                                         peulCv(1:nPeul,iCell)
        ELSE
          WRITE(IF_PROBE+iprobe-1,1005,IOSTAT=errorFlag) global%currentTime,  &
                                                         rho,u,v,w,press,temp,&
                                                         peulCv(1:nPeul,iCell)
        ENDIF
#endif
      ENDIF

      global%error = errorFlag
      IF (global%error /= 0) THEN
        CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,'Probe file' )
      ENDIF

! --- close and open probe file (instead of fflush)

      IF (global%probeOpenClose) THEN
        WRITE(fname,'(A,I4.4)') &
              TRIM(global%outDir)//TRIM(global%casename)//'.prb_',iprobe
        CLOSE(IF_PROBE+iprobe-1)
        OPEN(IF_PROBE+iprobe-1,FILE=fname,FORM='FORMATTED',STATUS='OLD', &
             POSITION='APPEND')
      ENDIF
    ENDIF   ! wrtProbe

  ENDDO     ! iprobe

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

!Modified by Mark Brandyberry to have more Decimal Places for Acoustics
1000 FORMAT(I6,1P,99E24.15)
1005 FORMAT(1PE14.7,99E24.15)

END SUBROUTINE WriteProbe

!******************************************************************************
!
! RCS Revision history:
!
! $Log: WriteProbe.F90,v $
! Revision 1.5  2009/01/06 21:29:46  mdbrandy
! Added More Decm
!
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/02/13 21:01:05  wasistho
! added ifdef PEUL
!
! Revision 1.1  2004/12/01 16:52:25  haselbac
! Initial revision after changing case
!
! Revision 1.21  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.20  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.19  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.18  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.15  2003/05/15 16:40:57  jblazek
! Changed index function call to fit into single line.
!
! Revision 1.14  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.13  2003/04/07 18:25:09  jferry
! added smoke concentrations to output
!
! Revision 1.12  2003/04/07 14:19:33  haselbac
! Removed ifdefs - now also used for RFLU
!
! Revision 1.11  2003/01/23 17:48:53  jblazek
! Changed algorithm to dump convergence, solution and probe data.
!
! Revision 1.10  2003/01/10 17:58:43  jblazek
! Added missing explicit interfaces.
!
! Revision 1.9  2002/10/07 19:24:28  haselbac
! Change use of IOSTAT, cures problem on SGIs
!
! Revision 1.8  2002/10/05 18:42:09  haselbac
! Added RFLU functionality
!
! Revision 1.7  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.4  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.3  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.2  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.1  2002/01/31 00:39:23  jblazek
! Probe output moved to common library.
!
!******************************************************************************







