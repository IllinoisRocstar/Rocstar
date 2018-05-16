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
! Purpose: write grid (and solution) data to plot file in binary format.
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
! $Id: POST_WriteTecplotBinary.F90,v 1.5 2008/12/06 08:44:49 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteTecplotBinary( iReg,region )

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
  INTEGER :: i, j, k, ii, jj, kk

! ... local variables
  CHARACTER(CHRLEN)   :: zoneName
  CHARACTER(2*CHRLEN) :: title
  CHARACTER(1)        :: nullchr

  INTEGER :: iLev, debug, vIsDouble, iret, iMax, jMax, kMax, ijkMax
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkN, cell(8), errorFlag

  REAL(RFREAL)          :: c, qq
  REAL(RFREAL), POINTER :: xyz(:,:), cv(:,:), dv(:,:), tv(:,:), gv(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: x(:,:,:), y(:,:,:), z(:,:,:), &
                                   rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:), &
                                   press(:,:,:), temp(:,:,:), mach(:,:,:)
#ifdef TURB
  REAL(RFREAL), POINTER :: vort(:,:), lens(:)
  DOUBLE PRECISION, ALLOCATABLE :: mut(:,:,:), tvort(:,:,:), len(:,:,:)
#endif
#ifdef STATS
  DOUBLE PRECISION, ALLOCATABLE :: srho(:,:,:),su(:,:,:), sv(:,:,:), sw(:,:,:), &
                                   spress(:,:,:), stemp(:,:,:), &
                                   suu(:,:,:), svv(:,:,:), sww(:,:,:), suv(:,:,:), &
                                   srr(:,:,:), spp(:,:,:), stt(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: smut(:,:,:), scdyn(:,:,:)
  REAL(RFREAL) :: rtime
  LOGICAL      :: statsActive
#endif

  TYPE(t_global), POINTER :: global
  TYPE(t_mixt)  , POINTER :: mixt
#ifdef TURB
  TYPE(t_turb)  , POINTER :: turb
#endif

! ... functions
  INTEGER :: TecIni, TecDat, TecZne, TecEnd

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'WriteTecplotBinary',&
  'POST_WriteTecplotBinary.F90' )

#ifdef NO_TECPLOT
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__ )
#else

! set parameters --------------------------------------------------------------

  nullchr   = CHAR(0)
  debug     = 0
  vIsDouble = 1
  iLev      = global%startLevel

  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  iMax   = ipnend - ipnbeg + 1
  jMax   = jpnend - jpnbeg + 1
  kMax   = kpnend - kpnbeg + 1
  ijkMax = iMax*jMax*kMax

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
    IF (global%flowType == FLOW_STEADY) THEN
      WRITE(title,1000) TRIM(global%casename),global%postIter
    ELSE
      WRITE(title,1005) TRIM(global%casename),global%postTime
    ENDIF
    IF (global%postPlotType == PLOT_GRID_ONLY) THEN
      iret = TecIni( TRIM(title)//nullchr, &
                     'x y z'//nullchr, &
                     TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                     debug,vIsDouble )
    ELSE
#ifndef TURB
#ifdef STATS
      IF (statsActive) THEN
        IF (global%mixtNStat > 0) THEN
          iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'//nullchr, &
                       'x y z rho u v w p T M <u> <v> <T> <uu> <vv> <ww> <uv> <TT>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
        ELSE
          iret = TecIni( TRIM(title)//nullchr, &
                       'x y z rho u v w p T M'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
        ENDIF
      ELSE
        iret = TecIni( TRIM(title)//nullchr, &
                     'x y z rho u v w p T M'//nullchr, &
                     TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                     debug,vIsDouble )
      ENDIF
#else
      iret = TecIni( TRIM(title)//nullchr, &
                   'x y z rho u v w p T M'//nullchr, &
                   TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                   debug,vIsDouble )
#endif ! Stats
#endif ! noTurb

#ifdef TURB
      IF ((global%postTurbFlag .EQV. .FALSE.) .OR. &
          region%mixtInput%turbModel == TURB_MODEL_NONE) THEN
#ifdef STATS
        IF (statsActive) THEN
          IF (global%mixtNStat > 0) THEN
            iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'//nullchr, &
                       'x y z rho u v w p T M <u> <v> <T> <uu> <vv> <ww> <uv> <TT>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
          ELSE
            iret = TecIni( TRIM(title)//nullchr, &
                       'x y z rho u v w p T M'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
          ENDIF
        ELSE
          iret = TecIni( TRIM(title)//nullchr, &
                       'x y z rho u v w p T M'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
        ENDIF
#else
        iret = TecIni( TRIM(title)//nullchr, &
                     'x y z rho u v w p T M'//nullchr, &
                     TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                     debug,vIsDouble )
#endif
      ELSE
        IF (region%turbInput%nOutField == 1) THEN
#ifdef STATS
          IF (statsActive) THEN
            IF (global%turbNStat > 0) THEN
              iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M mut <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp> <mut> <cdyn>'//nullchr, &
                       'x y z rho u v w p T M mut <u> <v> <T> <uu> <vv> <ww> <uv> <TT> <mut> <cdyn>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
            ELSE
              iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M mut <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'//nullchr, &
                       'x y z rho u v w p T M mut <u> <v> <T> <uu> <vv> <ww> <uv> <TT>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
            ENDIF
          ELSE
            iret = TecIni( TRIM(title)//nullchr, &
                       'x y z rho u v w p T M mut'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
          ENDIF
#else
          iret = TecIni( TRIM(title)//nullchr, &
                     'x y z rho u v w p T M mut'//nullchr, &
                     TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                     debug,vIsDouble )
#endif
        ELSEIF (region%turbInput%nOutField == 2) THEN
#ifdef STATS
          IF (statsActive) THEN
            IF (global%turbNStat > 0) THEN
              iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M mut tvort <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp> <mut> <cdyn>'//nullchr, &
                       'x y z rho u v w p T M mut tvort <u> <v> <T> <uu> <vv> <ww> <uv> <TT> <mut> <cdyn>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
            ELSE
              iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M mut tvort <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'//nullchr, &
                       'x y z rho u v w p T M mut tvort <u> <v> <T> <uu> <vv> <ww> <uv> <TT>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
            ENDIF
          ELSE
            iret = TecIni( TRIM(title)//nullchr, &
                       'x y z rho u v w p T M mut tvort'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
          ENDIF
#else
          iret = TecIni( TRIM(title)//nullchr, &
                     'x y z rho u v w p T M mut tvort'//nullchr, &
                     TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                     debug,vIsDouble )
#endif
        ENDIF
        IF ((region%turbInput%modelClass == MODEL_RANS ) .AND. &
            (region%turbInput%nOutField == 3)) THEN
#ifdef STATS
          IF (statsActive) THEN
            IF (global%turbNStat > 0) THEN
              iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M mut tvort lens <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp> <mut> <cdyn>'//nullchr, &
                       'x y z rho u v w p T M mut tvort lens <u> <v> <T> <uu> <vv> <ww> <uv> <TT> <mut> <cdyn>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
            ELSE
              iret = TecIni( TRIM(title)//nullchr, &
!                       'x y z rho u v w p T M mut tvort lens <u> <v> <w> <p> <uu> <vv> <ww> <uv> <pp>'//nullchr, &
                       'x y z rho u v w p T M mut tvort lens <u> <v> <T> <uu> <vv> <ww> <uv> <TT>'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
            ENDIF
          ELSE
            iret = TecIni( TRIM(title)//nullchr, &
                       'x y z rho u v w p T M mut tvort lens'//nullchr, &
                       TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                       debug,vIsDouble )
          ENDIF
#else
          iret = TecIni( TRIM(title)//nullchr, &
                     'x y z rho u v w p T M mut tvort lens'//nullchr, &
                     TRIM(global%casename)//'.plt'//nullchr,'.'//nullchr,  &
                     debug,vIsDouble )
#endif        ! Stats
        ENDIF ! nOutField
      ENDIF   ! postTurbFlag and turbModel
#endif    ! Turb
    ENDIF ! plotType
  ENDIF   ! iReg=1

! write zone header -----------------------------------------------------------

  WRITE(zoneName,1010) iReg

  iret = TecZne( TRIM(zoneName)//nullchr, &
                 iMax,jMax,kMax, &
                 'BLOCK'//nullchr,nullchr )

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

  IF (global%postTurbFlag .AND. &
      region%mixtInput%turbModel /= TURB_MODEL_NONE) THEN 
    IF (region%turbInput%nOutField > 1) vort => region%levels(iLev)%turb%vort
 
    IF (region%turbInput%modelClass == MODEL_RANS) &
      lens => region%levels(iLev)%turb%lens
  ENDIF
#endif

! allocate memory for temporary values

  ALLOCATE( x(iMax,jMax,kMax),stat=errorFlag )
  ALLOCATE( y(iMax,jMax,kMax),stat=errorFlag )
  ALLOCATE( z(iMax,jMax,kMax),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  IF (global%postPlotType == PLOT_GRID_FLOW) THEN
    ALLOCATE( rho  (iMax,jMax,kMax),stat=errorFlag )
    ALLOCATE( u    (iMax,jMax,kMax),stat=errorFlag )
    ALLOCATE( v    (iMax,jMax,kMax),stat=errorFlag )
    ALLOCATE( w    (iMax,jMax,kMax),stat=errorFlag )
    ALLOCATE( press(iMax,jMax,kMax),stat=errorFlag )
    ALLOCATE( temp (iMax,jMax,kMax),stat=errorFlag )
    ALLOCATE( mach (iMax,jMax,kMax),stat=errorFlag )

#ifdef TURB
    IF (global%postTurbFlag .AND. &
        region%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
      ALLOCATE( mut (iMax,jMax,kMax),stat=errorFlag )
      IF (region%turbInput%nOutField >= 2) THEN
        ALLOCATE( tvort (iMax,jMax,kMax),stat=errorFlag )
      ENDIF
      IF ((region%turbInput%modelClass == MODEL_RANS ) .AND. &
          (region%turbInput%nOutField == 3)) THEN
        ALLOCATE( len   (iMax,jMax,kMax),stat=errorFlag )
      ENDIF
    ENDIF
#endif

#ifdef STATS
    IF (statsActive) THEN
      IF (global%mixtNStat > 0) THEN
        ALLOCATE( srho  (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( su    (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( sv    (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( sw    (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( spress(iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( stemp (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( suu   (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( svv   (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( sww   (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( suv   (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( srr   (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( spp   (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( stt   (iMax,jMax,kMax),stat=errorFlag )
        sw = 0._RFREAL
      ENDIF
#ifdef TURB
      IF (global%turbNStat > 0) THEN
        ALLOCATE( smut  (iMax,jMax,kMax),stat=errorFlag )
        ALLOCATE( scdyn (iMax,jMax,kMax),stat=errorFlag )
      ENDIF
#endif
    ENDIF
#endif

  ENDIF
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! write grid coordinates

  kk = 0
  DO k=kpnbeg,kpnend
    kk = kk + 1
    jj = 0
    DO j=jpnbeg,jpnend
      jj = jj + 1
      ii = 0
      DO i=ipnbeg,ipnend
        ii   = ii + 1
        ijkN = IndIJK(i,j,k,iNOff,ijNOff)

        x(ii,jj,kk) = xyz(XCOORD,ijkN)
        y(ii,jj,kk) = xyz(YCOORD,ijkN)
        z(ii,jj,kk) = xyz(ZCOORD,ijkN)
      ENDDO
    ENDDO
  ENDDO

  iret = TecDat( ijkMax,x,vIsDouble )
  iret = TecDat( ijkMax,y,vIsDouble )
  iret = TecDat( ijkMax,z,vIsDouble )

! write solution

  IF (global%postPlotType == PLOT_GRID_FLOW) THEN
    kk = 0
    DO k=kpnbeg,kpnend
      kk = kk + 1
      jj = 0
      DO j=jpnbeg,jpnend
        jj = jj + 1
        ii = 0
        DO i=ipnbeg,ipnend
          ii      = ii + 1
          cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
          cell(2) = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
          cell(3) = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
          cell(4) = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
          cell(5) = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
          cell(6) = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
          cell(7) = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
          cell(8) = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)

          rho  (ii,jj,kk) = Aver(cell,CV_MIXT_DENS,cv)
          u    (ii,jj,kk) = AverDiv(cell,CV_MIXT_XMOM,cv,CV_MIXT_DENS,cv)
          v    (ii,jj,kk) = AverDiv(cell,CV_MIXT_YMOM,cv,CV_MIXT_DENS,cv)
          w    (ii,jj,kk) = AverDiv(cell,CV_MIXT_ZMOM,cv,CV_MIXT_DENS,cv)
          press(ii,jj,kk) = Aver(cell,DV_MIXT_PRES,dv)
          temp (ii,jj,kk) = Aver(cell,DV_MIXT_TEMP,dv)
          c               = Aver(cell,DV_MIXT_SOUN,dv)
          qq              = u(ii,jj,kk)*u(ii,jj,kk) + &
                            v(ii,jj,kk)*v(ii,jj,kk) + &
                            w(ii,jj,kk)*w(ii,jj,kk)
          mach (ii,jj,kk) = SQRT(qq)/c
#ifdef TURB
          IF (global%postTurbFlag .AND. &
              region%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
            mut  (ii,jj,kk) = Aver(cell,TV_MIXT_MUET,tv)
            IF (region%turbInput%nOutField >= 2) THEN
              tvort(ii,jj,kk) = Aver(cell,XCOORD,vort)
            ENDIF
            IF ((region%turbInput%modelClass == MODEL_RANS ) .AND. &
                (region%turbInput%nOutField == 3)) THEN
              len  (ii,jj,kk) = Aver1D(cell,lens)
            ENDIF
          ENDIF
#endif

#ifdef STATS
          IF (statsActive) THEN
            IF (global%mixtNStat > 0) THEN
!              su    (ii,jj,kk) = Aver(cell,2,mixt%tav)*rtime  ! 02  mixtStatId
!              sv    (ii,jj,kk) = Aver(cell,3,mixt%tav)*rtime  ! 03
!              sw    (ii,jj,kk) = Aver(cell,4,mixt%tav)*rtime  ! 04
!              spress(ii,jj,kk) = Aver(cell,5,mixt%tav)*rtime  ! 06
!              suu   (ii,jj,kk) = Aver(cell,7,mixt%tav)*rtime  ! 22
!              svv   (ii,jj,kk) = Aver(cell,8,mixt%tav)*rtime  ! 33
!              sww   (ii,jj,kk) = Aver(cell,9,mixt%tav)*rtime  ! 44
!              suv   (ii,jj,kk) = Aver(cell,10,mixt%tav)*rtime ! 23
!              spp   (ii,jj,kk) = Aver(cell,11,mixt%tav)*rtime ! 66
!              suu(ii,jj,kk) = suu(ii,jj,kk) - su(ii,jj,kk)*su(ii,jj,kk)
!              svv(ii,jj,kk) = svv(ii,jj,kk) - sv(ii,jj,kk)*sv(ii,jj,kk)
!              sww(ii,jj,kk) = sww(ii,jj,kk) - sw(ii,jj,kk)*sw(ii,jj,kk)
!              suv(ii,jj,kk) = suv(ii,jj,kk) - su(ii,jj,kk)*sv(ii,jj,kk)
!              spp(ii,jj,kk) = spp(ii,jj,kk) - spress(ii,jj,kk)*spress(ii,jj,kk)

              su    (ii,jj,kk) = Aver(cell,2,mixt%tav)*rtime  ! 02  mixtStatId
              sv    (ii,jj,kk) = Aver(cell,3,mixt%tav)*rtime  ! 03
              stemp (ii,jj,kk) = Aver(cell,4,mixt%tav)*rtime  ! 05
              suu   (ii,jj,kk) = Aver(cell,6,mixt%tav)*rtime  ! 22
              svv   (ii,jj,kk) = Aver(cell,7,mixt%tav)*rtime  ! 33
              sww   (ii,jj,kk) = Aver(cell,8,mixt%tav)*rtime  ! 44
              suv   (ii,jj,kk) = Aver(cell,9,mixt%tav)*rtime  ! 23
              stt   (ii,jj,kk) = Aver(cell,10,mixt%tav)*rtime ! 55
              suu(ii,jj,kk) = suu(ii,jj,kk) - su(ii,jj,kk)*su(ii,jj,kk)
              svv(ii,jj,kk) = svv(ii,jj,kk) - sv(ii,jj,kk)*sv(ii,jj,kk)
              sww(ii,jj,kk) = sww(ii,jj,kk) - sw(ii,jj,kk)*sw(ii,jj,kk)
              suv(ii,jj,kk) = suv(ii,jj,kk) - su(ii,jj,kk)*sv(ii,jj,kk)
              stt(ii,jj,kk) = stt(ii,jj,kk) - stemp(ii,jj,kk)*stemp(ii,jj,kk)
            ENDIF
#ifdef TURB
            IF (global%turbNStat > 0) THEN
              smut  (ii,jj,kk) = Aver(cell,1,turb%tav)*rtime  ! 02  turbStatId
              scdyn (ii,jj,kk) = Aver(cell,2,turb%tav)*rtime  ! 03
            ENDIF
#endif
          ENDIF  ! statsActive
#endif
        ENDDO
      ENDDO
    ENDDO

    iret = TecDat( ijkMax,rho  ,vIsDouble )
    iret = TecDat( ijkMax,u    ,vIsDouble )
    iret = TecDat( ijkMax,v    ,vIsDouble )
    iret = TecDat( ijkMax,w    ,vIsDouble )
    iret = TecDat( ijkMax,press,vIsDouble )
    iret = TecDat( ijkMax,temp ,vIsDouble )
    iret = TecDat( ijkMax,mach ,vIsDouble )

#ifdef TURB
    IF (global%postTurbFlag .AND. &
        region%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
      iret = TecDat( ijkMax,mut ,vIsDouble )
      IF (region%turbInput%nOutField >= 2) THEN
        iret = TecDat( ijkMax,tvort ,vIsDouble )
      ENDIF
      IF ((region%turbInput%modelClass == MODEL_RANS ) .AND. &
          (region%turbInput%nOutField == 3)) THEN
        iret = TecDat( ijkMax,len ,vIsDouble )
      ENDIF
    ENDIF
#endif

#ifdef STATS
    IF (statsActive) THEN
      IF (global%mixtNStat > 0) THEN
!        iret = TecDat( ijkMax,su    ,vIsDouble )
!        iret = TecDat( ijkMax,sv    ,vIsDouble )
!        iret = TecDat( ijkMax,sw    ,vIsDouble )
!        iret = TecDat( ijkMax,spress,vIsDouble )
!        iret = TecDat( ijkMax,suu   ,vIsDouble )
!        iret = TecDat( ijkMax,svv   ,vIsDouble )
!        iret = TecDat( ijkMax,sww   ,vIsDouble )
!        iret = TecDat( ijkMax,suv   ,vIsDouble )
!        iret = TecDat( ijkMax,spp   ,vIsDouble )

        iret = TecDat( ijkMax,su    ,vIsDouble )
        iret = TecDat( ijkMax,sv    ,vIsDouble )
        iret = TecDat( ijkMax,stemp ,vIsDouble )
        iret = TecDat( ijkMax,suu   ,vIsDouble )
        iret = TecDat( ijkMax,svv   ,vIsDouble )
        iret = TecDat( ijkMax,sww   ,vIsDouble )
        iret = TecDat( ijkMax,suv   ,vIsDouble )
        iret = TecDat( ijkMax,stt   ,vIsDouble )
      ENDIF
#ifdef TURB
      IF (global%turbNStat > 0) THEN  ! if (mixtNStat />0) it ll abort
        iret = TecDat( ijkMax,smut  ,vIsDouble )
        iret = TecDat( ijkMax,scdyn ,vIsDouble )
      ENDIF
#endif
    ENDIF  ! statsActive
#endif
  ENDIF  ! plotType == PLOT_GRID_FLOW

! close file ------------------------------------------------------------------

  IF (iReg == global%nRegions) THEN
    iret = TecEnd()
  ENDIF

! finalize --------------------------------------------------------------------

  DEALLOCATE( x,stat=errorFlag )
  DEALLOCATE( y,stat=errorFlag )
  DEALLOCATE( z,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  IF (global%postPlotType == PLOT_GRID_FLOW) THEN
    DEALLOCATE( rho  ,stat=errorFlag )
    DEALLOCATE( u    ,stat=errorFlag )
    DEALLOCATE( v    ,stat=errorFlag )
    DEALLOCATE( w    ,stat=errorFlag )
    DEALLOCATE( press,stat=errorFlag )
    DEALLOCATE( temp ,stat=errorFlag )
    DEALLOCATE( mach ,stat=errorFlag )
#ifdef TURB
    IF (global%postTurbFlag .AND. &
        region%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
      DEALLOCATE( mut  ,stat=errorFlag )
      IF (region%turbInput%nOutField >= 2) THEN
        DEALLOCATE( tvort ,stat=errorFlag )
      ENDIF
      IF ((region%turbInput%modelClass == MODEL_RANS ) .AND. &
          (region%turbInput%nOutField == 3)) THEN
        DEALLOCATE( len ,stat=errorFlag )
      ENDIF
    ENDIF
#endif
#ifdef STATS
    IF (statsActive) THEN
      IF (global%mixtNStat > 0) THEN
        DEALLOCATE( srho,   stat=errorFlag )
        DEALLOCATE( su,     stat=errorFlag )
        DEALLOCATE( sv,     stat=errorFlag )
        DEALLOCATE( sw,     stat=errorFlag )
        DEALLOCATE( spress, stat=errorFlag )
        DEALLOCATE( stemp,  stat=errorFlag )
        DEALLOCATE( suu,    stat=errorFlag )
        DEALLOCATE( svv,    stat=errorFlag )
        DEALLOCATE( sww,    stat=errorFlag )
        DEALLOCATE( suv,    stat=errorFlag )
        DEALLOCATE( srr,    stat=errorFlag )
        DEALLOCATE( spp,    stat=errorFlag )
        DEALLOCATE( stt,    stat=errorFlag )
      ENDIF
#ifdef TURB
      IF (global%turbNStat > 0) THEN
        DEALLOCATE( smut,   stat=errorFlag )
        DEALLOCATE( scdyn,  stat=errorFlag )
      ENDIF
#endif
    ENDIF
#endif
  ENDIF  ! plotType == PLOT_GRID_FLOW
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

#endif

  CALL DeregisterFunction( global )

! formats

1000 FORMAT(A,'. Iteration: ',I8,'.')
1005 FORMAT(A,'. Time: ',1PE11.5,'.')
1010 FORMAT(I5.5)

END SUBROUTINE WriteTecplotBinary

!******************************************************************************
!
! RCS Revision history:
!
! $Log: POST_WriteTecplotBinary.F90,v $
! Revision 1.5  2008/12/06 08:44:49  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/04/27 04:48:13  wasistho
! additinal statistics input set
!
! Revision 1.2  2004/12/03 03:20:36  wasistho
! rflo_modinterfacespost to post_modinterfaces
!
! Revision 1.1  2004/12/03 02:03:16  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:32:01  wasistho
! lower to upper case
!
! Revision 1.16  2004/11/10 18:29:59  wasistho
! put rtime within ifdef STATS
!
! Revision 1.15  2004/11/10 02:19:41  wasistho
! devided accumulated tav by integrTime
!
! Revision 1.14  2004/11/09 12:34:11  wasistho
! bug fixed, forgoten indeks ii,j,kk in getting rms of small scales
!
! Revision 1.13  2004/11/09 12:17:12  wasistho
! compute <u> = <uu>-<u><u> inside routine
!
! Revision 1.12  2004/11/09 10:50:23  wasistho
! added statistics to rflopost
!
! Revision 1.11  2004/08/26 00:44:38  wasistho
! bug fixed interconnection turbFlag and turbModel
!
! Revision 1.10  2004/07/24 03:48:28  wasistho
! use postSection instead of command line input
!
! Revision 1.9  2004/02/11 03:26:23  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.8  2004/02/07 01:19:06  wasistho
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
! Modified RegFun call to avoid probs with long 'POST_WriteTecplotBinary.F90' names
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
! Revision 1.5  2002/04/11 21:10:27  jblazek
! Set correct time when writing grid only for Tecplot.
!
! Revision 1.4  2002/02/22 20:30:39  jblazek
! Changed generic format. Enhanced Tecplot title.
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.1  2002/01/12 00:02:49  jblazek
! Added postprocessor.
!
!******************************************************************************








