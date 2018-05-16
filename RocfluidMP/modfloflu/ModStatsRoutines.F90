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
! Purpose: Suite of statistics routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModStatsRoutines.F90,v 1.11 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE ModStatsRoutines

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI

#ifdef PLAG
  USE PLAG_ModEulerian, ONLY : PLAG_CalcEulerianField
#endif
 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: StatBuildVersionString, &
            GetStatistics, &
            InitStatistics, &
            StatDataAccumulation1, &
            StatDataAccumulation2, &
            StatDataSampling, &
            StatMapping, &
            StatTimeAccumulation, &
            StatWriteMP
        
! ******************************************************************************
! Declarations and definitions
! ****************************************************************************** 
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: ModStatsRoutines.F90,v $ $Revision: 1.11 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
!******************************************************************************
!
! Purpose: Build version string for printing in header.
!
! Description: none.
!
! Input: none.
!
! Output: 
!   versionString = string containing version number and date.
!
! Notes: 
!   1. The strings are NOT to be edited by anyone except the developer of 
!      this physical module. 
!
!******************************************************************************

SUBROUTINE StatBuildVersionString( versionString )

  IMPLICIT NONE

! ... parameters
  CHARACTER(*) :: versionString

! ... local variables
  CHARACTER(LEN=2)  :: major, minor, patch
  CHARACTER(CHRLEN) :: date

!******************************************************************************
! set strings: DO NOT EDIT UNLESS YOU ARE STATS DEVELOPER

  major = '1'
  minor = '0'
  patch = '0'

  date  = '11/30/04'

! write into string

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

END SUBROUTINE StatBuildVersionString

!******************************************************************************
!
! Purpose: calling data sampling routines for time averaged statistics
!
! Description: data sampling for time averaged statistics of gas mixture 
!              and other physical module if desired; in addition
!              averaging time interval is integrated
!
! Input: regions = data of all regions
!
! Output: StatDataSampling and StatTimeAccumulation called
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE GetStatistics( regions )

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_global), POINTER :: global

#ifdef PLAG
    TYPE(t_region), POINTER :: pRegion
#endif

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'GetStatistics',&
  'ModStatsRoutines.F90' )

! perform data sampling --------------------------------------------------------

  IF (global%doStat == ACTIVE) THEN
#ifdef RFLO
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor
#endif
#ifdef RFLU
    DO iReg=1,global%nRegionsLocal
#endif
        IF (global%mixtNStat>0) CALL StatDataSampling( regions(iReg),FTYPE_MIXT )
#ifdef TURB
        IF (global%turbNStat>0) THEN
          IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
              (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) THEN
            CALL StatDataSampling( regions(iReg),FTYPE_TURB )
          ENDIF
        ENDIF
#endif
#ifdef PLAG
        IF ( global%plagUsed .EQV. .TRUE. ) THEN
          pRegion => regions(iReg)

          CALL PLAG_CalcEulerianField( pRegion )
          IF ( global%plagNStat>0 ) & 
            CALL StatDataSampling( regions(iReg),FTYPE_PLAG ) 
        ENDIF ! plagUsed
#endif
#ifdef RFLO
      ENDIF
#endif
    
    ENDDO
    CALL StatTimeAccumulation( global ) 
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE GetStatistics

!******************************************************************************
!
! Purpose: initiation of time averaged statistics
!
! Description: initiation proceeds for gas mixture and other physical 
!              module if desired; it is initiated by zero for new
!              time averaging and by old values for restarting from 
!              the previous process
!
! Input: regions = data of all regions
!
! Output: regions%levels%mixt%tav = time averaged mixture variables
!         regions%levels%turb%tav = time averaged TURB variables
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE InitStatistics( regions )

#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_ReadStat
#ifdef PLAG
  USE PLAG_RFLO_ModStats, ONLY : PLAG_RFLO_ReadStat
#endif
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY : RFLU_ReadStat
#ifdef GENX
  USE STAT_RFLU_ModRocstarAdmin, ONLY : STAT_RFLU_GenxGetData
#endif
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER :: iLev

  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'InitStatistics',&
  'ModStatsRoutines.F90' )

! check if the flow is unsteady -----------------------------------------------

  IF (global%flowType == FLOW_STEADY)  THEN
    global%doStat = OFF 
  ENDIF

! initialization of time averaging statistics ---------------------------------

  global%statBc = 0

  IF (global%doStat == ACTIVE) THEN

! - restart from previous statistics

    IF (global%reStat == ACTIVE) THEN

      IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE)  &
        WRITE(STDOUT,'(A)') SOLVER_NAME//' Restart statistics ...'
#ifndef GENX
#ifdef RFLO
      CALL RFLO_ReadStat( regions )
#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) &
        CALL PLAG_RFLO_ReadStat( regions )
#endif
#endif
#ifdef RFLU
      DO iReg=1,global%nRegionsLocal
        CALL RFLU_ReadStat( regions(iReg) )
      ENDDO
#endif

#else

#ifdef RFLU
      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL STAT_RFLU_GenxGetData( pRegion )
      ENDDO
#endif
#endif

! - initialize new statistics

    ELSE

      IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE)  &
        WRITE(STDOUT,'(A)') SOLVER_NAME//' Start new statistics ...'
      global%integrTime = 0._RFREAL      
#ifdef RFLO
      DO iReg=1,global%nRegions
        IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
            regions(iReg)%active==ACTIVE) THEN              ! on my processor
          iLev =  regions(iReg)%currLevel
          IF (global%mixtNStat > 0) THEN
            regions(iReg)%levels(iLev)%mixt%tav = 0._RFREAL
          ENDIF
#ifdef TURB
          IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
              (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) .AND. &
              (global%turbNStat > 0)) THEN
            regions(iReg)%levels(iLev)%turb%tav = 0._RFREAL
          ENDIF
#endif
#ifdef PLAG
          IF ((global%plagUsed .EQV. .TRUE.) .AND. &
              (global%plagNStat > 0)) THEN
            regions(iReg)%levels(iLev)%plag%tav = 0._RFREAL
          ENDIF
#endif
        ENDIF     ! region on this processor and active
      ENDDO       ! iReg
#endif
#ifdef RFLU
      DO iReg=1,global%nRegionsLocal
        IF (global%mixtNStat > 0) THEN
          regions(iReg)%mixt%tav =  0._RFREAL
        ENDIF
#ifdef TURB
        IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
            (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) .AND. &
            (global%turbNStat > 0)) THEN
          regions(iReg)%turb%tav =  0._RFREAL
        ENDIF
#endif
#ifdef PLAG
        IF ((global%plagUsed .EQV. .TRUE.) .AND. &
            (global%plagNStat > 0)) THEN
          regions(iReg)%plag%tav = 0._RFREAL
        ENDIF
#endif
      ENDDO       ! iReg
#endif

    ENDIF
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE InitStatistics

!******************************************************************************
!
! Purpose: accumulate selected instantaneous quantities in time
!
! Description: this file containes two routines:
!              1st routine, StatDataAccumulation1, for first moment variables
!              2nd routine, StatDataAccumulation2, for second moment variables
!
! Input 1st routine:    ijkbeg, ijkend = begin and end cell indices
!                       id             = index of variable var
!                       idtav          = index of time-accumulated quantity
!                       dTime          = time step
!                       var            = variable to be time accumulated
!                      
! Input 2nd routine:    ijkbeg, ijkend = begin and end cell indices
!                       id1            = index of first component variable
!                       id2            = index of second component variable
!                       idtav          = index of time-accumulated quantity
!                       dTime          = time step
!                       var1           = 1st component of 2nd moment variable
!                       var2           = 2nd component of 2nd moment variable 
!                      
! Output of 1st and 2nd routine: qavg = quantity accumulated in time
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE StatDataAccumulation1( ijkbeg,ijkend,id,idtav,dTime,var,qavg )

  IMPLICIT NONE

! ... parameters
  INTEGER               :: ijkbeg, ijkend, id, idtav
  REAL(RFREAL)          :: dTime
  REAL(RFREAL), POINTER :: var(:,:), qavg(:,:)

! ... loop variables
  INTEGER :: ijk

!******************************************************************************

  DO ijk=ijkbeg,ijkend
    qavg(idtav,ijk) = qavg(idtav,ijk) + dTime*var(id,ijk)
  ENDDO

END SUBROUTINE StatDataAccumulation1

! #############################################################################

SUBROUTINE StatDataAccumulation2( ijkbeg,ijkend,id1,id2,idtav, &
                                  dTime,var1,var2,qavg )

  IMPLICIT NONE

! ... parameters
  INTEGER               :: ijkbeg, ijkend, id1, id2, idtav
  REAL(RFREAL)          :: dTime
  REAL(RFREAL), POINTER :: var1(:,:), var2(:,:), qavg(:,:)

! ... loop variables
  INTEGER :: ijk

!******************************************************************************

  DO ijk=ijkbeg,ijkend
    qavg(idtav,ijk) = qavg(idtav,ijk) + dTime*var1(id1,ijk)*var2(id2,ijk)
  ENDDO

END SUBROUTINE StatDataAccumulation2

!******************************************************************************
!
! Purpose: sorting variables to be time averaged
!
! Description: sorting procedure based on selection taken by user 
!              stored in statId and mapped into statCode.
!
! Input: region              = data of current region
!        fluidType           = mixture, turb, plag, etc        
!        cv, dv, tv, gv      = variables to be time averaged
!        statCode            = mapping identifiers
!   
! Output: tav  = quantities accumulated in time
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE StatDataSampling( region,fluidType )

#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: fluidType

! ... loop variables
  INTEGER :: l

! ... local variables
  INTEGER :: iLev, iOff, ijOff
  INTEGER :: idcbeg, jdcbeg, kdcbeg, idcend, jdcend, kdcend, ijkbeg, ijkend
  INTEGER :: nStat, id1, id2
  INTEGER, POINTER      :: statCode(:,:,:)

  REAL(RFREAL), POINTER :: cv(:,:),dv(:,:),tv(:,:),gv(:,:),ev(:,:),tav(:,:)
  REAL(RFREAL), POINTER :: sv(:,:),st(:,:)
  REAL(RFREAL), POINTER :: var1(:,:),var2(:,:) 

!******************************************************************************

  CALL RegisterFunction( region%global,'StatDataSampling',&
  'ModStatsRoutines.F90' )

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel
#endif

  IF (fluidType == FTYPE_MIXT) THEN
    nStat = region%global%mixtNStat
    statCode => region%global%mixtStatCode
#ifdef RFLO
    cv    => region%levels(iLev)%mixt%cv
    dv    => region%levels(iLev)%mixt%dv
    tv    => region%levels(iLev)%mixt%tv
    gv    => region%levels(iLev)%mixt%gv
    tav   => region%levels(iLev)%mixt%tav
#endif
#ifdef RFLU
    cv    => region%mixt%cv
    dv    => region%mixt%dv
    tv    => region%mixt%tv
    gv    => region%mixt%gv
    tav   => region%mixt%tav
#endif
  ELSEIF (fluidType == FTYPE_TURB) THEN
#ifdef TURB
    nStat = region%global%turbNStat
    statCode => region%global%turbStatCode
#ifdef RFLO
    tv    => region%levels(iLev)%mixt%tv
    dv    => region%levels(iLev)%turb%dv
    sv    => region%levels(iLev)%turb%sv
    st    => region%levels(iLev)%turb%st
    tav   => region%levels(iLev)%turb%tav
#endif
#ifdef RFLU
    tv    => region%mixt%tv
    dv    => region%turb%dv
    sv    => region%turb%sv
    st    => region%turb%st
    tav   => region%turb%tav
#endif
#endif
  ELSEIF (fluidType == FTYPE_PLAG) THEN
#ifdef PLAG
    nStat = region%global%plagNStat
    statCode => region%global%plagStatCode
#ifdef RFLO
    ev    => region%levels(iLev)%plag%ev
    tav   => region%levels(iLev)%plag%tav
#endif
#ifdef RFLU
    ev    => region%plag%ev
    tav   => region%plag%tav
#endif
#endif
  ELSEIF (fluidType == FTYPE_PEUL) THEN
#ifdef PEUL
    nStat = region%global%peulNStat
    statCode => region%global%peulStatCode
#ifdef RFLO
!    cv    => region%levels(iLev)%peul%cv
!    dv    => region%levels(iLev)%peul%dv
!    tav   => region%levels(iLev)%peul%tav
#endif
#ifdef RFLU
!    cv    => region%peul%cv
!    dv    => region%peul%dv
!    tav   => region%peul%tav
#endif
#endif
  ENDIF

#ifdef RFLO
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iOff,ijOff )
  ijkbeg = IndIJK(idcbeg,jdcbeg,kdcbeg,iOff,ijOff)
  ijkend = IndIJK(idcend,jdcend,kdcend,iOff,ijOff)
#endif
#ifdef RFLU
  ijkbeg = 1
  ijkend = region%grid%nCells
#endif

! Quantities to be time-averaged is determined from the index selected by user.
! Data accumulation proceeds afterwards for each quantity.

  DO l=1,nStat

    IF ((statCode(1,1,l)==STAT_NONE).AND.(statCode(2,1,l)==STAT_NONE)) &
    GOTO 999

    IF (statCode(1,1,l)==STAT_CV) THEN
      var1 => cv
    ELSEIF (statCode(1,1,l)==STAT_DV) THEN
      var1 => dv
    ELSEIF (statCode(1,1,l)==STAT_TV) THEN
      var1 => tv
    ELSEIF (statCode(1,1,l)==STAT_GV) THEN
      var1 => gv
    ELSEIF (statCode(1,1,l)==STAT_SV) THEN
      var1 => sv
    ELSEIF (statCode(1,1,l)==STAT_ST) THEN
      var1 => st
    ELSEIF (statCode(1,1,l)==STAT_PLAGEV) THEN
      var1 => ev
    ENDIF

    IF (statCode(2,1,l)==STAT_CV) THEN
      var2 => cv
    ELSEIF (statCode(2,1,l)==STAT_DV) THEN
      var2 => dv
    ELSEIF (statCode(2,1,l)==STAT_TV) THEN
      var2 => tv
    ELSEIF (statCode(2,1,l)==STAT_GV) THEN
      var2 => gv
    ELSEIF (statCode(2,1,l)==STAT_SV) THEN
      var2 => sv
    ELSEIF (statCode(2,1,l)==STAT_ST) THEN
      var2 => st
    ELSEIF (statCode(2,1,l)==STAT_PLAGEV) THEN
      var2 => ev
    ENDIF

    id1 = statCode(1,2,l)
    id2 = statCode(2,2,l)

    IF (statCode(1,1,l)==STAT_NONE) THEN
      CALL StatDataAccumulation1( ijkbeg,ijkend,id2,l, &
                                  region%global%dtMin,var2,tav )
    ELSEIF (statCode(2,1,l)==STAT_NONE) THEN
      CALL StatDataAccumulation1( ijkbeg,ijkend,id1,l, &
                                  region%global%dtMin,var1,tav )
    ELSE
      CALL StatDataAccumulation2( ijkbeg,ijkend,id1,id2,l, &
                                  region%global%dtMin,var1,var2,tav )
    ENDIF

  ENDDO

! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( region%global )

END SUBROUTINE StatDataSampling

!******************************************************************************
!
! Purpose: mapping from mixture and physical module ID to statistics ID
!
! Description: mapping based on parameter mixtStatId, turbStatId, etc 
!              input by user
!
! Input: global%mixtStatId : mixture statistics ID from user input
!        global%turbStatId : TURB statistics ID from user input
!
! Output: global%mixtStatCode : mapped mixture statistics ID
!         global%turbStatCode : mapped TURB statistics ID
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE StatMapping( global )

#ifdef GENX
  USE ModInterfacesStatistics, ONLY : GenxStatNaming
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_StatMapping
#endif
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_StatMapping
#endif
#ifdef PEUL
!  USE ModInterfacesEulerian,   ONLY : PEUL_StatMapping
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER  :: l, n

! ... local variables
#ifdef GENX
  CHARACTER(CHRLEN), POINTER :: statName(:,:,:)
#endif
  INTEGER :: errorFlag
  INTEGER, POINTER  :: statId(:,:), statCode(:,:,:)

!******************************************************************************

  CALL RegisterFunction( global,'StatMapping',&
  'ModStatsRoutines.F90' )

! allocate mixture variables and set pointers ---------------------------------

  IF (global%mixtNStat <= 0) GOTO 111

  ALLOCATE( global%mixtStatCode(2,2,global%mixtNStat),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  statId   => global%mixtStatId
  statCode => global%mixtStatCode

#ifdef GENX
  ALLOCATE( global%mixtStatNm(2,2,global%mixtNStat),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  statName => global%mixtStatNm
#endif

! set mixture mapping ---------------------------------------------------------

  DO l=1,global%mixtNStat
    DO n=1,2
      IF (statId(n,l)==0) THEN 
        statCode(n,:,l) = STAT_NONE
      ELSE IF (statId(n,l)==1) THEN 
        statCode(n,1,l) = STAT_CV
        statCode(n,2,l) = CV_MIXT_DENS
#ifdef GENX
        statName(n,1,l) = 'rho'
        statName(n,2,l) = 'kg/m^3'
#endif
      ELSE IF (statId(n,l)==2) THEN 
#ifdef RFLO
        statCode(n,1,l) = STAT_DV
        statCode(n,2,l) = DV_MIXT_UVEL
#endif
#ifdef RFLU
        statCode(n,1,l) = STAT_CV
        statCode(n,2,l) = CV_MIXT_XVEL
#endif
#ifdef GENX
        statName(n,1,l) = 'u'
        statName(n,2,l) = 'm/s'
#endif
      ELSE IF (statId(n,l)==3) THEN 
#ifdef RFLO
        statCode(n,1,l) = STAT_DV
        statCode(n,2,l) = DV_MIXT_VVEL
#endif
#ifdef RFLU
        statCode(n,1,l) = STAT_CV 
        statCode(n,2,l) = CV_MIXT_YVEL 
#endif
#ifdef GENX
        statName(n,1,l) = 'v'
        statName(n,2,l) = 'm/s'
#endif
      ELSE IF (statId(n,l)==4) THEN 
#ifdef RFLO
        statCode(n,1,l) = STAT_DV
        statCode(n,2,l) = DV_MIXT_WVEL
#endif
#ifdef RFLU
        statCode(n,1,l) = STAT_CV 
        statCode(n,2,l) = CV_MIXT_ZVEL 
#endif
#ifdef GENX
        statName(n,1,l) = 'w'
        statName(n,2,l) = 'm/s'
#endif
      ELSE IF (statId(n,l)==5) THEN 
        statCode(n,1,l) = STAT_DV
        statCode(n,2,l) = DV_MIXT_TEMP
#ifdef GENX
        statName(n,1,l) = 'T'
        statName(n,2,l) = 'K'
#endif
      ELSE IF (statId(n,l)==6) THEN 
        statCode(n,1,l) = STAT_DV
        statCode(n,2,l) = DV_MIXT_PRES
#ifdef GENX
        statName(n,1,l) = 'p'
        statName(n,2,l) = 'N/m^2'
#endif
      ELSE IF (statId(n,l)==7) THEN 
        statCode(n,1,l) = STAT_DV
        statCode(n,2,l) = DV_MIXT_SOUN
#ifdef GENX
        statName(n,1,l) = 'c'
        statName(n,2,l) = 'm/s'
#endif
      ELSE IF (statId(n,l)==8) THEN 
        statCode(n,1,l) = STAT_TV
        statCode(n,2,l) = TV_MIXT_MUEL
#ifdef GENX
        statName(n,1,l) = 'mul'
        statName(n,2,l) = 'kg/ms'
#endif
      ELSE IF (statId(n,l)==9) THEN 
        statCode(n,1,l) = STAT_TV
        statCode(n,2,l) = TV_MIXT_TCOL
#ifdef GENX
        statName(n,1,l) = 'tcol'
        statName(n,2,l) = 'kg m/Ks^3'
#endif
      ELSE
        CALL ErrorStop( global,ERR_STATS_INDEXING,__LINE__, &
                        'mixture index out of range.' ) 
      ENDIF
    ENDDO
  ENDDO

#ifdef GENX
! defined names and units of mixture statistics

  CALL GenxStatNaming( global, FTYPE_MIXT )
#endif

111 CONTINUE

#ifdef TURB
! allocate TURB variables and set pointers ---------------------------------

  IF (global%turbNStat <= 0) GOTO 222

! set TURB mapping ---------------------------------------------------------

  CALL TURB_StatMapping( global )

222 CONTINUE
#endif

#ifdef PLAG
! allocate PLAG variables and set pointers ---------------------------------

  IF (global%plagNStat <= 0) GOTO 333

! set PLAG mapping ---------------------------------------------------------

  CALL PLAG_StatMapping( global )

333 CONTINUE
#endif

#ifdef PEUL
! allocate PEUL variables and set pointers ---------------------------------

  IF (global%peulNStat <= 0) GOTO 444

! set PEUL mapping ---------------------------------------------------------

!  CALL PEUL_StatMapping( global )

444 CONTINUE
#endif

#ifdef RFLO
! check with regard to inlet turbulence recycling -----------------------------
! stats id 01-06 must exist and lined up in the same order

  IF (global%infloNijk < NIJK_INFLOW_INIT) THEN
    DO l=1,global%mixtNStat
      IF (statId(1,l)==0) THEN 
        IF (statId(2,l) /= l) THEN
          CALL ErrorStop( global,ERR_STATS_INDEXING,__LINE__, &
              'For recycturb inflow, set stats Id 01-06 in the given order.' ) 
        ENDIF
      ENDIF
    ENDDO
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE StatMapping

!******************************************************************************
!
! Purpose: accumulation of time during time averaging process
!
! Description: time integration proceeds from STARTTIME to MAXTIME set by user
!              if statistics RESTART (global%reStat) active, the process is
!              continuation from previous process
!
! Input: global%dtMin = global minimum time advancement.
!
! Output: global%integrTime = accumulated time.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE StatTimeAccumulation( global )

  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables

!******************************************************************************

  CALL RegisterFunction( global,'StatTimeAccumulation',&
  'ModStatsRoutines.F90' )

! accumulate time steps used in time-weighted averaging process

  global%integrTime = global%integrTime + global%dtMin

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE StatTimeAccumulation

!******************************************************************************
!
! Purpose: update bc and write time averaged solution for all active physical 
!          modules
!
! Description: none.
!
! Input: regions = tav of all physical modules in all regions 
!
! Output: to file
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModStatsRoutines.F90,v 1.11 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE StatWriteMP( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_WriteStat
  USE RFLO_ModStatsBoundaryConditions, ONLY : RFLO_StatBoundaryConditionsSet
#ifdef PLAG
  USE PLAG_RFLO_ModStats, ONLY : PLAG_RFLO_WriteStat, &
                                 PLAG_RFLO_CommStatBuffWrapper
#endif
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#ifdef RFLO
#include "Indexing.h"
#endif

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables

! ... local variables

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'StatWriteMP',&
  'ModStatsRoutines.F90' )

! start -----------------------------------------------------------------------

#ifdef RFLO
#ifdef PLAG
  CALL PLAG_RFLO_CommStatBuffWrapper( regions )
#endif
  CALL RFLO_StatBoundaryConditionsSet( regions )

  CALL RFLO_WriteStat( regions )
#ifdef PLAG
  CALL PLAG_RFLO_WriteStat( regions )
#endif
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE StatWriteMP

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE ModStatsRoutines

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModStatsRoutines.F90,v $
! Revision 1.11  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/01/10 05:02:34  wasistho
! Get statistics Genx data in RFLU
!
! Revision 1.8  2005/12/20 20:42:39  wasistho
! added ifdef RFLO around kernel using global%Nijk
!
! Revision 1.7  2005/12/01 09:00:10  wasistho
! sanity check stats Id with inlet turbulence
!
! Revision 1.6  2005/06/16 03:52:18  wasistho
! activated RFLO_ModStatsBc
!
! Revision 1.5  2005/05/21 07:07:50  wasistho
! backout RFLO_ModStatsBoundaryConditions temporarily
!
! Revision 1.4  2005/05/21 01:42:47  wasistho
! added statWriteMP
!
! Revision 1.3  2005/01/11 01:28:33  wasistho
! fixed bugs, PLAG data sampling was outside regions loop
!
! Revision 1.2  2005/01/08 20:36:44  fnajjar
! Added statistics infrastructure for PLAG and activated datastructure appropriately
!
! Revision 1.1  2004/12/29 23:28:21  wasistho
! moved ModStatisticsRoutines from libfloflu to modfloflu
!
! Revision 1.1  2004/12/28 20:30:00  wasistho
! moved statistics routines into module ModStatsRoutines
!
!
! ******************************************************************************












