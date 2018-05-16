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
! Purpose: write solution data at grid nodes to generic file.
!
! Description: none.
!
! Input: iReg   = region number
!        region = region data (dimensions, flow variables)
!
! Output: to plot file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: POST_WriteGeneric.F90,v 1.4 2008/12/06 08:44:49 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteGeneric( iReg,region )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt
  USE POST_ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, &
                                 RFLO_GetCellOffset, Aver1D, Aver, AverDiv
  USE ModMPI
  USE ModParameters
#ifdef TURB
  USE ModTurbulence, ONLY : t_turb
  USE TURB_ModParameters
#endif
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, id

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, cell(8), errorFlag

  REAL(RFREAL)          :: rho, u, v, w, press, temp, c
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), tv(:,:), gv(:,:)
#ifdef TURB
  REAL(RFREAL)          :: mut, tvort, len
  REAL(RFREAL), POINTER :: vort(:,:), lens(:)
#endif
#ifdef STATS
  REAL(RFREAL)          :: srho,su,sv,sw,spress,srr,suu,svv,sww,suv,spp
  REAL(RFREAL)          :: rtime
#endif

  TYPE(t_global), POINTER :: global
  TYPE(t_mixt)  , POINTER :: mixt
#ifdef TURB
  TYPE(t_turb)  , POINTER :: turb
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'WriteGeneric',&
  'POST_WriteGeneric.F90' )

! get dimensions --------------------------------------------------------------

  iLev = global%startLevel
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

! open file and write the header ----------------------------------------------

  IF (iReg == 1) THEN
    OPEN(IF_PLOT,file=TRIM(global%casename)//'.sol_node',status='unknown', &
         form='unformatted',iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__ )
    WRITE(IF_PLOT,err=10) global%nRegions,global%postIter,global%postTime
  ENDIF

! write zone header -----------------------------------------------------------

  WRITE(IF_PLOT,err=10) iReg,iLev, &
                        region%levels(iLev)%grid%ipc+1, &
                        region%levels(iLev)%grid%jpc+1, &
                        region%levels(iLev)%grid%kpc+1

! write data ------------------------------------------------------------------
! pointer to variables

  mixt => region%levels(iLev)%mixt
  cv   => mixt%cv
  dv   => mixt%dv
  tv   => mixt%tv
  gv   => mixt%gv
#ifdef TURB
  turb => region%levels(iLev)%turb

  IF (global%postTurbFlag) THEN 
    IF (region%turbInput%nOutField > 1) vort => turb%vort
 
    IF (region%turbInput%modelClass == MODEL_RANS) &
      lens => turb%lens
  ENDIF
#endif

#ifdef STATS
! set 1/integrTime
  rtime = 1._RFREAL/global%integrTime
#endif

! write solution

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        cell(2) = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        cell(3) = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        cell(4) = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
        cell(5) = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        cell(6) = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
        cell(7) = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
        cell(8) = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)

        rho   = Aver(cell,CV_MIXT_DENS,cv)
        u     = AverDiv(cell,CV_MIXT_XMOM,cv,CV_MIXT_DENS,cv)
        v     = AverDiv(cell,CV_MIXT_YMOM,cv,CV_MIXT_DENS,cv)
        w     = AverDiv(cell,CV_MIXT_ZMOM,cv,CV_MIXT_DENS,cv)
        press = Aver(cell,DV_MIXT_PRES,dv)
        temp  = Aver(cell,DV_MIXT_TEMP,dv)
        c     = Aver(cell,DV_MIXT_SOUN,dv)
#ifndef TURB
        WRITE(IF_PLOT,err=10) rho,u,v,w,press,temp,c
#endif
#ifdef TURB
        IF (global%postTurbFlag .EQV. .FALSE. .OR. &
            region%turbInput%modelClass == MODEL_NONE ) THEN
          WRITE(IF_PLOT,err=10) rho,u,v,w,press,temp,c
        ENDIF
        IF (global%postTurbFlag) THEN 
          mut = Aver(cell,TV_MIXT_MUET,tv)
          IF (region%turbInput%nOutField == 1) THEN
            WRITE(IF_PLOT,err=10) rho,u,v,w,press,temp,c,mut
          ELSEIF (region%turbInput%nOutField == 2) THEN
            tvort = Aver(cell,XCOORD,vort)
            WRITE(IF_PLOT,err=10) rho,u,v,w,press,temp,c,mut,tvort
          ENDIF
          IF ((region%turbInput%modelClass == MODEL_RANS ) .AND. &
              (region%turbInput%nOutField == 3)) THEN
            tvort = Aver(cell,XCOORD,vort)
            len   = Aver1D(cell,lens)
            WRITE(IF_PLOT,err=10) rho,u,v,w,press,temp,c,mut,tvort,len
          ENDIF
        ENDIF ! postTurbFlag
#endif
      ENDDO
    ENDDO
  ENDDO

#ifdef STATS
  IF (global%postStatsFlag) THEN 
    IF ((global%flowType == FLOW_UNSTEADY) .AND. &
        (global%doStat == ACTIVE)) THEN
      IF (global%mixtNStat > 0) THEN
        DO k=kpnbeg,kpnend
          DO j=jpnbeg,jpnend
            DO i=ipnbeg,ipnend
              cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
              cell(2) = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
              cell(3) = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
              cell(4) = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
              cell(5) = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
              cell(6) = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
              cell(7) = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
              cell(8) = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)

              srho   = Aver(cell,1,mixt%tav)*rtime  ! 01   mixtStatId
              su     = Aver(cell,2,mixt%tav)*rtime  ! 02 
              sv     = Aver(cell,3,mixt%tav)*rtime  ! 03
              sw     = Aver(cell,4,mixt%tav)*rtime  ! 04
              spress = Aver(cell,5,mixt%tav)*rtime  ! 06
              srr    = Aver(cell,6,mixt%tav)*rtime  ! 11
              suu    = Aver(cell,7,mixt%tav)*rtime  ! 22
              svv    = Aver(cell,8,mixt%tav)*rtime  ! 33
              sww    = Aver(cell,9,mixt%tav)*rtime  ! 44
              suv    = Aver(cell,10,mixt%tav)*rtime ! 23
              spp    = Aver(cell,11,mixt%tav)*rtime ! 66
              srr = srr - srho*srho
              suu = suu - su*su
              svv = svv - sv*sv
              sww = sww - sw*sw
              suv = suv - su*sv
              spp = spp - spress*spress

              WRITE(IF_PLOT,err=10) &
                srho,su,sv,sw,spress,srr,suu,svv,sww,suv,spp
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k
      ENDIF  ! mixtNStat >0
#ifdef TURB
      IF (global%turbNStat > 0) THEN
        DO k=kpnbeg,kpnend
          DO j=jpnbeg,jpnend
            DO i=ipnbeg,ipnend
              cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
              cell(2) = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
              cell(3) = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
              cell(4) = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
              cell(5) = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
              cell(6) = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
              cell(7) = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
              cell(8) = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)

              WRITE(IF_PLOT,err=10) &
                (Aver(cell,id,turb%tav)*rtime, id=1,global%turbNStat)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k
      ENDIF  ! turbNStat >0
#endif
    ENDIF  ! doStat
  ENDIF  ! postStatsFlag
#endif

! close file ------------------------------------------------------------------

  IF (iReg == global%nRegions) THEN
    CLOSE(IF_PLOT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__ )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__ )

999  CONTINUE

END SUBROUTINE WriteGeneric

!******************************************************************************
!
! RCS Revision history:
!
! $Log: POST_WriteGeneric.F90,v $
! Revision 1.4  2008/12/06 08:44:49  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:20:21  wasistho
! rflo_modinterfacespost to post_modinterfaces
!
! Revision 1.1  2004/12/03 02:03:16  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:32:01  wasistho
! lower to upper case
!
! Revision 1.17  2004/11/10 18:29:47  wasistho
! put rtime within ifdef STATS
!
! Revision 1.16  2004/11/10 02:20:06  wasistho
! devided accumulated tav by integrTime
!
! Revision 1.15  2004/11/09 12:16:44  wasistho
! compute <u> = <uu>-<u><u> inside routine
!
! Revision 1.14  2004/11/09 10:50:10  wasistho
! added statistics to rflopost
!
! Revision 1.13  2004/07/24 03:48:15  wasistho
! use postSection instead of command line input
!
! Revision 1.12  2004/02/11 03:26:07  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.11  2004/02/07 01:18:51  wasistho
! added turbulence related results in rocflo post processing
!
! Revision 1.10  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.9  2003/03/20 22:23:47  haselbac
! Renamed ModInterfaces
!
! Revision 1.8  2003/03/20 19:41:26  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.7  2003/03/20 19:34:37  haselbac
! Modified RegFun call to avoid probs with long 'POST_WriteGeneric.F90' names
!
! Revision 1.6  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.5  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.3  2002/02/22 20:30:39  jblazek
! Changed generic format. Enhanced Tecplot title.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/02/16 07:17:58  jblazek
! Added generic binary output format.
!
!******************************************************************************








