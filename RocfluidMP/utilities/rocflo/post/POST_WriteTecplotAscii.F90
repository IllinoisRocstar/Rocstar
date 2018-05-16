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
! Purpose: write grid (and solution) data to plot file in ASCII format.
!
! Description: none.
!
! Input: iReg    = region number
!        region  = region data (dimensions, flow variables)
!
! Output: to plot file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: POST_WriteTecplotAscii.F90,v 1.5 2008/12/06 08:44:49 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteTecplotAscii( iReg,region )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt
  USE POST_ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, &
                                 RFLO_GetNodeOffset, &
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
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN+4) :: fname

  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkN, cell(8), errorFlag

  REAL(RFREAL)          :: rho, u, v, w, press, temp, mach, c, qq
  REAL(RFREAL), POINTER :: xyz(:,:), cv(:,:), dv(:,:), tv(:,:), gv(:,:)
#ifdef TURB
  REAL(RFREAL)          :: mut, tvort, len
  REAL(RFREAL), POINTER :: vort(:,:), lens(:)
#endif
#ifdef STATS
  REAL(RFREAL) :: su, sv, sw, spress, suu, svv, sww, suv, spp
  REAL(RFREAL) :: rtime, smut, scdyn
  LOGICAL :: statsActive
#endif

  TYPE(t_global), POINTER :: global
  TYPE(t_mixt)  , POINTER :: mixt
#ifdef TURB
  TYPE(t_turb)  , POINTER :: turb
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'WriteTecplotAscii',&
  'POST_WriteTecplotAscii.F90' )

! set parameters --------------------------------------------------------------

  iLev = global%startLevel
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! statistics active?

#ifdef STATS
  statsActive = (global%postStatsFlag .AND. &
                (global%flowType == FLOW_UNSTEADY) .AND. &
                (global%doStat == ACTIVE))

! set 1/integrTime
  rtime = 1._RFREAL/global%integrTime
#endif

! open file and write the header ----------------------------------------------

  IF (iReg == 1) THEN
    fname = TRIM(global%casename)//'.dat'
    OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

    IF (global%flowType == FLOW_STEADY) THEN
      WRITE(IF_PLOT,1000,err=10) TRIM(global%casename),global%postIter
    ELSE
      WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%postTime
    ENDIF
    IF (global%postPlotType == PLOT_GRID_ONLY) THEN
      WRITE(IF_PLOT,1010,err=10) 'x y z'
    ELSE
#ifndef TURB
#ifdef STATS
      IF (statsActive) THEN
        IF (global%mixtNStat > 0) THEN
          WRITE(IF_PLOT,1010,err=10) &
            'x y z rho u v w p T M <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'
        ELSE
          WRITE(IF_PLOT,1010,err=10) &
            'x y z rho u v w p T M'
        ENDIF
      ELSE
        WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M'
      ENDIF
#else
      WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M'
#endif 
#endif 

#ifdef TURB
      IF ((global%postTurbFlag .EQV. .FALSE.) .OR. &
          region%mixtInput%turbModel==TURB_MODEL_NONE) THEN
#ifdef STATS
        IF (statsActive) THEN
          IF (global%mixtNStat > 0) THEN
            WRITE(IF_PLOT,1010,err=10) &
              'x y z rho u v w p T M <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'
          ELSE
            WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M'
          ENDIF
        ELSE
          WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M'
        ENDIF
#else
        WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M'
#endif 
      ENDIF
      IF (global%postTurbFlag) THEN 
        IF (region%turbInput%nOutField == 1) THEN
#ifdef STATS
          IF (statsActive) THEN
            IF (global%turbNStat > 0) THEN
              WRITE(IF_PLOT,1010,err=10) &
                'x y z rho u v w p T M mut <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp> <mut> <cdyn>'
            ELSE
              WRITE(IF_PLOT,1010,err=10) &
                'x y z rho u v w p T M mut <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'
            ENDIF
          ELSE
            WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M mut'
          ENDIF
#else
          WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M mut'
#endif 
        ELSEIF (region%turbInput%nOutField == 2) THEN
#ifdef STATS
          IF (statsActive) THEN
            IF (global%turbNStat > 0) THEN
              WRITE(IF_PLOT,1010,err=10) &
                'x y z rho u v w p T M mut tvort <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp> <mut> <cdyn>'
            ELSE
              WRITE(IF_PLOT,1010,err=10) &
                'x y z rho u v w p T M mut tvort <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'
            ENDIF
          ELSE
            WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M mut tvort'
          ENDIF
#else
          WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M mut tvort'
#endif 
        ENDIF
        IF ((region%turbInput%modelClass == MODEL_RANS) .AND. &
            (region%turbInput%nOutField == 3)) THEN
#ifdef STATS
          IF (statsActive) THEN
            IF (global%turbNStat > 0) THEN
              WRITE(IF_PLOT,1010,err=10) &
                'x y z rho u v w p T M mut tvort lens <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp> <mut> <cdyn>'
            ELSE
              WRITE(IF_PLOT,1010,err=10) &
                'x y z rho u v w p T M mut tvort lens <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'
            ENDIF
          ELSE
            WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M mut tvort lens'
          ENDIF
#else
          WRITE(IF_PLOT,1010,err=10) 'x y z rho u v w p T M mut tvort lens'
#endif 
        ENDIF
      ENDIF
#endif
    ENDIF
  ENDIF   ! iReg=1

! write zone header

  WRITE(IF_PLOT,1015) iReg,ipnend-ipnbeg+1,jpnend-jpnbeg+1,kpnend-kpnbeg+1

! write data ------------------------------------------------------------------
! pointer to variables

  xyz  => region%levels(iLev)%grid%xyz
  mixt => region%levels(iLev)%mixt
  cv   => mixt%cv
  dv   => mixt%dv
  tv   => mixt%tv
  gv   => mixt%gv
#ifdef TURB
  turb => region%levels(iLev)%turb

  IF (global%postTurbFlag) THEN 
    IF (region%turbInput%nOutField > 1) vort => region%levels(iLev)%turb%vort
 
    IF (region%turbInput%modelClass == MODEL_RANS) &
      lens => region%levels(iLev)%turb%lens
  ENDIF
#endif

! write grid coordinates only

  IF (global%postPlotType == PLOT_GRID_ONLY) THEN

    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend
        DO i=ipnbeg,ipnend
          ijkN = IndIJK(i,j,k,iNOff,ijNOff)
          WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                     xyz(YCOORD,ijkN), &
                                     xyz(ZCOORD,ijkN)
        ENDDO
      ENDDO
    ENDDO

! write grid and solution

  ELSE

    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend
        DO i=ipnbeg,ipnend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
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
          qq    = u*u + v*v + w*w
          mach  = SQRT(qq)/c
#ifdef STATS
          IF (statsActive) THEN
            IF (global%mixtNStat > 0) THEN
              su     = Aver(cell,2,mixt%tav)*rtime   ! 02  mixtStatId
              sv     = Aver(cell,3,mixt%tav)*rtime   ! 03
              sw     = Aver(cell,4,mixt%tav)*rtime   ! 04
              spress = Aver(cell,5,mixt%tav)*rtime   ! 06
              suu    = Aver(cell,7,mixt%tav)*rtime   ! 22
              svv    = Aver(cell,8,mixt%tav)*rtime   ! 33
              sww    = Aver(cell,9,mixt%tav)*rtime   ! 44
              suv    = Aver(cell,10,mixt%tav)*rtime  ! 23
              spp    = Aver(cell,11,mixt%tav)*rtime  ! 66
              suu = suu - su*su
              svv = svv - sv*sv
              sww = sww - sw*sw
              suv = suv - su*sv
              spp = spp - spress*spress
            ENDIF
#ifdef TURB
            IF (global%turbNStat > 0) THEN
              smut   = Aver(cell,1,turb%tav)*rtime   ! 01  turbStatId
              scdyn  = Aver(cell,2,turb%tav)*rtime   ! 03
            ENDIF
#endif
          ENDIF  ! statsActive
#endif

#ifndef TURB
#ifdef STATS
        IF (statsActive) THEN
          IF (global%mixtNStat > 0) THEN
            WRITE(IF_PLOT,1030,err=10) xyz(XCOORD,ijkN), &
                                       xyz(YCOORD,ijkN), &
                                       xyz(ZCOORD,ijkN), &
                                     rho,u,v,w,press,temp,mach, &
                                     su,sv,sw,spress,suu,svv,sww,suv,spp
          ELSE
            WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                       xyz(YCOORD,ijkN), &
                                       xyz(ZCOORD,ijkN), &
                                       rho,u,v,w,press,temp,mach
          ENDIF
        ELSE
          WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                     xyz(YCOORD,ijkN), &
                                     xyz(ZCOORD,ijkN), &
                                     rho,u,v,w,press,temp,mach
        ENDIF
#else
        WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                   xyz(YCOORD,ijkN), &
                                   xyz(ZCOORD,ijkN), &
                                   rho,u,v,w,press,temp,mach
#endif 
#endif 

#ifdef TURB
          IF ((global%postTurbFlag .EQV. .FALSE.) .OR. &
              region%turbInput%modelClass == MODEL_NONE) THEN
#ifdef STATS
            IF (statsActive) THEN
              IF (global%mixtNStat > 0) THEN
                WRITE(IF_PLOT,1030,err=10) xyz(XCOORD,ijkN), &
                                           xyz(YCOORD,ijkN), &
                                           xyz(ZCOORD,ijkN), &
                                        rho,u,v,w,press,temp,mach,&
                                        su,sv,sw,spress,suu,svv,sww,suv,spp
              ELSE
                WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                           xyz(YCOORD,ijkN), &
                                           xyz(ZCOORD,ijkN), &
                                           rho,u,v,w,press,temp,mach
              ENDIF
            ELSE
              WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                         xyz(YCOORD,ijkN), &
                                         xyz(ZCOORD,ijkN), &
                                         rho,u,v,w,press,temp,mach
            ENDIF
#else
            WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                       xyz(YCOORD,ijkN), &
                                       xyz(ZCOORD,ijkN), &
                                       rho,u,v,w,press,temp,mach
#endif 
          ENDIF
          IF (global%postTurbFlag) THEN
            mut = Aver(cell,TV_MIXT_MUET,tv)
            IF (region%turbInput%nOutField == 1) THEN
#ifdef STATS
              IF (statsActive) THEN
                IF (global%turbNStat > 0) THEN
                  WRITE(IF_PLOT,1041,err=10) xyz(XCOORD,ijkN), &
                                             xyz(YCOORD,ijkN), &
                                             xyz(ZCOORD,ijkN), &
                                     rho,u,v,w,press,temp,mach,mut, &
                                     su,sv,sw,spress,suu,svv,sww,suv,spp, &
                                     smut,scdyn
                ELSE
                  WRITE(IF_PLOT,1031,err=10) xyz(XCOORD,ijkN), &
                                             xyz(YCOORD,ijkN), &
                                             xyz(ZCOORD,ijkN), &
                                     rho,u,v,w,press,temp,mach,mut, &
                                     su,sv,sw,spress,suu,svv,sww,suv,spp
                ENDIF
              ELSE
                WRITE(IF_PLOT,1021,err=10) xyz(XCOORD,ijkN), &
                                           xyz(YCOORD,ijkN), &
                                           xyz(ZCOORD,ijkN), &
                                           rho,u,v,w,press,temp,mach,mut
              ENDIF
#else
              WRITE(IF_PLOT,1021,err=10) xyz(XCOORD,ijkN), &
                                         xyz(YCOORD,ijkN), &
                                         xyz(ZCOORD,ijkN), &
                                         rho,u,v,w,press,temp,mach,mut
#endif 
            ELSEIF (region%turbInput%nOutField == 2) THEN
              tvort = Aver(cell,XCOORD,vort)
#ifdef STATS
              IF (statsActive) THEN
                IF (global%turbNStat > 0) THEN
                  WRITE(IF_PLOT,1042,err=10) xyz(XCOORD,ijkN), &
                                             xyz(YCOORD,ijkN), &
                                             xyz(ZCOORD,ijkN), &
                                      rho,u,v,w,press,temp,mach,mut,tvort, &
                                      su,sv,sw,spress,suu,svv,sww,suv,spp, &
                                      smut,scdyn
                ELSE
                  WRITE(IF_PLOT,1032,err=10) xyz(XCOORD,ijkN), &
                                             xyz(YCOORD,ijkN), &
                                             xyz(ZCOORD,ijkN), &
                                      rho,u,v,w,press,temp,mach,mut,tvort, &
                                      su,sv,sw,spress,suu,svv,sww,suv,spp
                ENDIF
              ELSE
                WRITE(IF_PLOT,1022,err=10) xyz(XCOORD,ijkN), &
                                           xyz(YCOORD,ijkN), &
                                           xyz(ZCOORD,ijkN), &
                                           rho,u,v,w,press,temp,mach,mut,tvort
              ENDIF
#else
              WRITE(IF_PLOT,1022,err=10) xyz(XCOORD,ijkN), &
                                         xyz(YCOORD,ijkN), &
                                         xyz(ZCOORD,ijkN), &
                                         rho,u,v,w,press,temp,mach,mut,tvort
#endif 
            ENDIF
            IF ((region%turbInput%modelClass == MODEL_RANS ) .AND. &
                (region%turbInput%nOutField == 3)) THEN
              tvort = Aver(cell,XCOORD,vort)
              len   = Aver1D(cell,lens)
#ifdef STATS
              IF (statsActive) THEN
                IF (global%turbNStat > 0) THEN
                  WRITE(IF_PLOT,1043,err=10) xyz(XCOORD,ijkN), &
                                             xyz(YCOORD,ijkN), &
                                             xyz(ZCOORD,ijkN), &
                                       rho,u,v,w,press,temp,mach,mut,tvort,len, &
                                       su,sv,sw,spress,suu,svv,sww,suv,spp,smut,scdyn
                ELSE
                  WRITE(IF_PLOT,1033,err=10) xyz(XCOORD,ijkN), &
                                             xyz(YCOORD,ijkN), &
                                             xyz(ZCOORD,ijkN), &
                                       rho,u,v,w,press,temp,mach,mut,tvort,len, &
                                       su,sv,sw,spress,suu,svv,sww,suv,spp
                ENDIF
              ELSE
                WRITE(IF_PLOT,1023,err=10) xyz(XCOORD,ijkN), &
                                           xyz(YCOORD,ijkN), &
                                           xyz(ZCOORD,ijkN), &
                                           rho,u,v,w,press,temp,mach,mut,tvort,len
              ENDIF
#else
              WRITE(IF_PLOT,1023,err=10) xyz(XCOORD,ijkN), &
                                         xyz(YCOORD,ijkN), &
                                         xyz(ZCOORD,ijkN), &
                                         rho,u,v,w,press,temp,mach,mut,tvort,len
#endif 
            ENDIF
          ENDIF ! postTurbFlag
#endif
        ENDDO
      ENDDO
    ENDDO

  ENDIF   ! postPlotType

! close file, handle errors ---------------------------------------------------

  IF (iReg == global%nRegions) THEN
    CLOSE(IF_PLOT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )
  ENDIF

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,fname )

! formats ---------------------------------------------------------------------

1000 FORMAT('TITLE="',A,'. Iteration: ',I8,'."')
1005 FORMAT('TITLE="',A,'. Time: ',1PE11.5,'."')
1010 FORMAT('VARIABLES= ',A)
1015 FORMAT('ZONE T="',I5.5,'", I=',I6,', J=',I6,', K=',I6,', F=POINT')
1020 FORMAT(1P,10(1X,E13.6))
1021 FORMAT(1P,11(1X,E13.6))
1022 FORMAT(1P,12(1X,E13.6))
1023 FORMAT(1P,13(1X,E13.6))
1030 FORMAT(1P,19(1X,E13.6))
1031 FORMAT(1P,20(1X,E13.6))
1032 FORMAT(1P,21(1X,E13.6))
1033 FORMAT(1P,22(1X,E13.6))
1041 FORMAT(1P,22(1X,E13.6))
1042 FORMAT(1P,23(1X,E13.6))
1043 FORMAT(1P,24(1X,E13.6))

999  CONTINUE
END SUBROUTINE WriteTecplotAscii

!******************************************************************************
!
! RCS Revision history:
!
! $Log: POST_WriteTecplotAscii.F90,v $
! Revision 1.5  2008/12/06 08:44:49  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/12/24 04:07:59  wasistho
! put brackets when comparing logicals in if statements
!
! Revision 1.2  2004/12/03 03:20:31  wasistho
! rflo_modinterfacespost to post_modinterfaces
!
! Revision 1.1  2004/12/03 02:03:16  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:32:01  wasistho
! lower to upper case
!
! Revision 1.15  2004/11/16 04:20:44  wasistho
! replaced extension .plt to .dat
!
! Revision 1.14  2004/11/10 18:29:54  wasistho
! put rtime within ifdef STATS
!
! Revision 1.13  2004/11/10 02:19:46  wasistho
! devided accumulated tav by integrTime
!
! Revision 1.12  2004/11/09 12:17:30  wasistho
! compute <u> = <uu>-<u><u> inside routine
!
! Revision 1.11  2004/11/09 10:50:31  wasistho
! added statistics to rflopost
!
! Revision 1.10  2004/07/24 03:48:22  wasistho
! use postSection instead of command line input
!
! Revision 1.9  2004/02/11 03:26:16  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.8  2004/02/07 01:19:00  wasistho
! added turbulence related results in rocflo post processing
!
! Revision 1.7  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.6  2003/03/20 22:23:47  haselbac
! Renamed ModInterfaces
!
! Revision 1.5  2003/03/20 19:41:26  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.4  2003/03/20 19:34:37  haselbac
! Modified RegFun call to avoid probs with long 'POST_WriteTecplotAscii.F90' names
!
! Revision 1.3  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/07/20 00:42:05  jblazek
! Added ASCII Tecplot format.
!
!******************************************************************************








