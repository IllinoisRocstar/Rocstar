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
! Purpose: Suite of boundary condition routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModBoundaryConditions.F90,v 1.30 2010/02/18 21:47:39 juzhang Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModBoundaryConditions

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModBndPatch, ONLY  : t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_BcondFarfield, &
            RFLO_BcondInflow, &
            RFLO_BcondInjection, &
            RFLO_BcondInjectionInit, &
            RFLO_BcondNoslipWall, &
            RFLO_BcondOutflow, &
            RFLO_BcondRotatPeriod, &
            RFLO_BcondSlipWall, &
            RFLO_BcondSymmetry, &
            RFLO_BoundaryConditionsRecv, &
            RFLO_BoundaryConditionsSend, &
            RFLO_BoundaryConditionsSet

! private : RFLO_BcondInflowVTPerf
!           RFLO_BcondInflowVPPerf
!           RFLO_BcondInflowFluc
!           RFLO_BcondInjectionAPN
 
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModBoundaryConditions.F90,v $ $Revision: 1.30 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: update values in dummy cells at far field boundary patch.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: there is no vortex correction implemented yet.
!
!******************************************************************************

SUBROUTINE RFLO_BcondFarfield( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, BcondFarfieldPerf, &
        MixtureProperties
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, n1, n2

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, lbound, iCOff, ijCOff, iNOff, ijNOff, nOff, distrib
  INTEGER :: indCp, indMol, gasModel, ijkC, ijkC1, ijkD, ijkN, i2d
  INTEGER :: inode, jnode, knode

  REAL(RFREAL)          :: sgn, ds, sxn, syn, szn
  REAL(RFREAL)          :: mach, attack, slip, press, temp, pb
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), vals(:,:), sFace(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondFarfield',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib
  indCp     = region%levels(iLev)%mixt%indCp
  indMol    = region%levels(iLev)%mixt%indMol
  gasModel = region%mixtInput%gasModel

  cv   => region%levels(iLev)%mixt%cv
  dv   => region%levels(iLev)%mixt%dv
  gv   => region%levels(iLev)%mixt%gv
  vals => patch%mixt%vals

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector

  IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
  IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
  IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)

          IF (idum == 1) THEN
            ijkC = IndIJK(i,j,k,iCOff,ijCOff)
            ijkN = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            ds   = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                        sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                        sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
            sxn  = sgn*sFace(XCOORD,ijkN)/ds
            syn  = sgn*sFace(YCOORD,ijkN)/ds
            szn  = sgn*sFace(ZCOORD,ijkN)/ds

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
            ENDIF
            i2d = distrib * IndIJ(n1,n2,nOff)

            IF (gasModel == GAS_MODEL_TCPERF) THEN
              mach   = vals(BCDAT_FARF_MACH  ,i2d)
              attack = vals(BCDAT_FARF_ATTACK,i2d)
              slip   = vals(BCDAT_FARF_SLIP  ,i2d)
              press  = vals(BCDAT_FARF_PRESS ,i2d)
              temp   = vals(BCDAT_FARF_TEMP  ,i2d)
              CALL BcondFarfieldPerf( mach,attack,slip,press,temp,sxn,syn,szn, &
                     gv(GV_MIXT_CP  ,ijkC*indCp),gv(GV_MIXT_MOL ,ijkC*indMol), &
                     cv(CV_MIXT_DENS,ijkC)      ,cv(CV_MIXT_XMOM,ijkC)       , &
                     cv(CV_MIXT_YMOM,ijkC)      ,cv(CV_MIXT_ZMOM,ijkC)       , &
                     cv(CV_MIXT_ENER,ijkC)      ,dv(DV_MIXT_PRES,ijkC)       , &
                     cv(CV_MIXT_DENS,ijkD)      ,cv(CV_MIXT_XMOM,ijkD)       , &
                     cv(CV_MIXT_YMOM,ijkD)      ,cv(CV_MIXT_ZMOM,ijkD)       , &
                     cv(CV_MIXT_ENER,ijkD)      ,pb )
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,&
                   __LINE__ )
            ENDIF

          ELSE   ! idum > 1
            ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
            ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

            cv(CV_MIXT_DENS,ijkD) = 2._RFREAL*cv(CV_MIXT_DENS,ijkC ) - &
                                              cv(CV_MIXT_DENS,ijkC1)
            cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_XMOM,ijkC ) - &
                                              cv(CV_MIXT_XMOM,ijkC1)
            cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_YMOM,ijkC ) - &
                                              cv(CV_MIXT_YMOM,ijkC1)
            cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_ZMOM,ijkC ) - &
                                              cv(CV_MIXT_ZMOM,ijkC1)
            cv(CV_MIXT_ENER,ijkD) = 2._RFREAL*cv(CV_MIXT_ENER,ijkC ) - &
                                              cv(CV_MIXT_ENER,ijkC1)
          ENDIF

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondFarfield

!******************************************************************************
!
! Purpose: update values in dummy cells at inflow boundary patch.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BcondInflow( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, BcondInflowPerf, &
        MixtureProperties
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, n1, n2

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, lbound, iCOff, ijCOff, iNOff, ijNOff, nOff, bcType, &
             distrib, bcOptFixed, bcOptType, bcOptModel
  INTEGER :: indCp, indMol, gasModel, ijkC, ijkC1, ijkD, ijkN, i2d
  INTEGER :: inode, jnode, knode
  INTEGER :: ijkF, ifluc

  REAL(RFREAL)          :: sgn, dS, sxn, syn, szn, pr, mach, amp, eps, rAvgTim
  REAL(RFREAL)          :: fluc(BCDAT_INFLOW_NELM)
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), vals(:,:), sFace(:,:)
  REAL(RFREAL), POINTER :: tav(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondInflow',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers

  ifluc  = region%global%infloNijk
  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nOff       = ABS(patch%l1end-patch%l1beg) + 1
  bcType     = patch%bcType
  distrib    = patch%mixt%distrib
  bcOptType  = patch%mixt%switches(BCSWI_INFLOW_TYPE)
  bcOptFixed = patch%mixt%switches(BCSWI_INFLOW_FIXED)
  bcOptModel = patch%mixt%switches(BCSWI_INFLOW_MODEL)
  amp        = patch%mixt%amplitude
  indCp      = region%levels(iLev)%mixt%indCp
  indMol     = region%levels(iLev)%mixt%indMol
  gasModel   = region%mixtInput%gasModel

  cv   => region%levels(iLev)%mixt%cv
  dv   => region%levels(iLev)%mixt%dv
  gv   => region%levels(iLev)%mixt%gv
  vals => patch%mixt%vals

! initial fluctuations and pointer for turbulent inflow case

  fluc    = 0._RFREAL
  eps     = 100._RFREAL*EPSILON( 1._RFREAL )

#ifdef STATS
  IF (bcOptModel == BCOPT_UNSTEADY) THEN
    rAvgTim = 1._RFREAL/(region%global%integrTime + eps)
    tav  => region%levels(iLev)%mixt%tav
!if(region%global%myProcid==MASTERPROC) write(*,*)'ok1', region%global%integrTime
  ENDIF
#endif

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector

  IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
  IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
  IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! loop over all cells of the patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)

          IF (idum == 1) THEN
            ijkC = IndIJK(i,j,k,iCOff,ijCOff)
            ijkN = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            dS   = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                        sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                        sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
            sxn  = sgn*sFace(XCOORD,ijkN)/dS
            syn  = sgn*sFace(YCOORD,ijkN)/dS
            szn  = sgn*sFace(ZCOORD,ijkN)/dS

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
            ENDIF
            i2d = distrib * IndIJ(n1,n2,nOff)

            IF (lbound==1) THEN
              ijkF = IndIJK(ibeg+ifluc-1,j,k,iCOff,ijCOff)
            ELSEIF (lbound==2) THEN
              ijkF = IndIJK(iend-ifluc+1,j,k,iCOff,ijCOff)
            ELSEIF (lbound==3) THEN
              ijkF = IndIJK(i,jbeg+ifluc-1,k,iCOff,ijCOff)
            ELSEIF (lbound==4) THEN
              ijkF = IndIJK(i,jbeg-ifluc+1,k,iCOff,ijCOff)
            ELSEIF (lbound==5) THEN
              ijkF = IndIJK(i,j,kbeg+ifluc-1,iCOff,ijCOff)
            ELSEIF (lbound==6) THEN
              ijkF = IndIJK(i,j,kbeg-ifluc+1,iCOff,ijCOff)
            ENDIF

            IF (gasModel == GAS_MODEL_TCPERF) THEN
#ifdef STATS
              IF ((bcType == BC_INFLOW_VELTEMP .OR. &
                   bcType == BC_INFLOW_VELPRESS) .AND. &
                  (region%global%integrTime > eps) .AND. &
                  (region%global%infloNijk < NIJK_INFLOW_INIT)) THEN
                CALL RFLO_BcondInflowFluc( bcOptModel, & 
                                           fluc(BCDAT_INFLOW_U), &
                                           fluc(BCDAT_INFLOW_V), &
                                           fluc(BCDAT_INFLOW_W), &
                                           fluc(BCDAT_INFLOW_T), &
                                           fluc(BCDAT_INFLOW_P), &
                                           dv(DV_MIXT_UVEL,ijkF), &
                                           dv(DV_MIXT_VVEL,ijkF), &
                                           dv(DV_MIXT_WVEL,ijkF), &
                                           dv(DV_MIXT_TEMP,ijkF), &
                                           dv(DV_MIXT_PRES,ijkF), &
                                           tav(2,ijkF)*rAvgTim, &
                                           tav(3,ijkF)*rAvgTim, &
                                           tav(4,ijkF)*rAvgTim, &
                                           tav(5,ijkF)*rAvgTim, &
                                           tav(6,ijkF)*rAvgTim )
!if (region%iRegionGlobal==1 .AND. j==17 .AND. k==17) then
!write(*,*)'tav',tav(2,ijkF)*ravgtim,tav(3,ijkF)*ravgtim,tav(4,ijkF)*ravgtim,tav(5,ijkF)*ravgtim,tav(6,ijkF)*ravgtim
!write(*,*)'dv ',dv(DV_MIXT_UVEL,ijkF),dv(DV_MIXT_VVEL,ijkF),dv(DV_MIXT_WVEL,ijkF),dv(DV_MIXT_TEMP,ijkF),dv(DV_MIXT_PRES,ijkF)
!endif
              ENDIF
#endif
              IF (bcType == BC_INFLOW_TOTANG) THEN

                IF (bcOptType == BCOPT_SUBSONIC) THEN
                  mach = 0._RFREAL
                ELSE
                  mach = vals(BCDAT_INFLOW_MACH,i2d)
                ENDIF
                IF (bcOptModel == BCOPT_UNSTEADY) THEN
                  CALL ErrorStop( region%global,ERR_VAL_BCVAL,&
                       __LINE__, &
                      'RECYCTURB > 0 not available for BC_INFLOW_TOTANG' )
                ENDIF
                CALL BcondInflowPerf( bcOptType,bcOptFixed        , & 
                                      vals(BCDAT_INFLOW_PTOT ,i2d), &
                                      vals(BCDAT_INFLOW_TTOT ,i2d), &
                                      vals(BCDAT_INFLOW_BETAH,i2d), &
                                      vals(BCDAT_INFLOW_BETAV,i2d), &
                                      mach,sxn,syn,szn            , & 
                                      gv(GV_MIXT_CP  ,ijkC*indCp ), &
                                      gv(GV_MIXT_MOL ,ijkC*indMol), &
                                      cv(CV_MIXT_DENS,ijkC)       , &
                                      cv(CV_MIXT_XMOM,ijkC)       , &
                                      cv(CV_MIXT_YMOM,ijkC)       , &
                                      cv(CV_MIXT_ZMOM,ijkC)       , &
                                      cv(CV_MIXT_DENS,ijkD)       , &
                                      cv(CV_MIXT_XMOM,ijkD)       , &
                                      cv(CV_MIXT_YMOM,ijkD)       , &
                                      cv(CV_MIXT_ZMOM,ijkD)       , &
                                      cv(CV_MIXT_ENER,ijkD),pr )

              ELSEIF (bcType == BC_INFLOW_VELTEMP) THEN

                CALL RFLO_BcondInflowVTPerf( bcOptModel, bcOptType, &
                     vals(BCDAT_INFLOW_U ,i2d)+amp*fluc(BCDAT_INFLOW_U), &
                     vals(BCDAT_INFLOW_V ,i2d)+amp*fluc(BCDAT_INFLOW_V), &
                     vals(BCDAT_INFLOW_W ,i2d)+amp*fluc(BCDAT_INFLOW_W), &
                     vals(BCDAT_INFLOW_T ,i2d)+amp*fluc(BCDAT_INFLOW_T), &
                     vals(BCDAT_INFLOW_P ,i2d)+amp*fluc(BCDAT_INFLOW_P), &
                     sxn,syn,szn                                   , & 
                     gv(GV_MIXT_CP  ,ijkC*indCp),gv(GV_MIXT_MOL ,ijkC*indMol), &
                     cv(CV_MIXT_DENS,ijkC)      ,cv(CV_MIXT_XMOM,ijkC)       , &
                     cv(CV_MIXT_YMOM,ijkC)      ,cv(CV_MIXT_ZMOM,ijkC)       , &
                     cv(CV_MIXT_ENER,ijkC)      ,dv(DV_MIXT_PRES,ijkC)       , &
                     cv(CV_MIXT_DENS,ijkD)      ,cv(CV_MIXT_XMOM,ijkD)       , &
                     cv(CV_MIXT_YMOM,ijkD)      ,cv(CV_MIXT_ZMOM,ijkD)       , &
                     cv(CV_MIXT_ENER,ijkD) )
#ifdef STATS
!                IF ((region%global%integrTime > eps) .AND. & 
!                    (region%global%infloNijk < NIJK_INFLOW_INIT)) THEN
!                    cv(CV_MIXT_XMOM,ijkD) = cv(CV_MIXT_XMOM,ijkD) + &
!                           (vals(BCDAT_INFLOW_U,i2d) - tav(2,ijkD)*rAvgTim)* &
!                            cv(CV_MIXT_DENS,ijkD)     
!                ENDIF                  
#endif
              ELSEIF (bcType == BC_INFLOW_VELPRESS) THEN

                CALL RFLO_BcondInflowVPPerf( bcOptModel, bcOptType, &
                     vals(BCDAT_INFLOW_U ,i2d)+amp*fluc(BCDAT_INFLOW_U), &
                     vals(BCDAT_INFLOW_V ,i2d)+amp*fluc(BCDAT_INFLOW_V), &
                     vals(BCDAT_INFLOW_W ,i2d)+amp*fluc(BCDAT_INFLOW_W), &
                     vals(BCDAT_INFLOW_T ,i2d)+amp*fluc(BCDAT_INFLOW_T), &
                     vals(BCDAT_INFLOW_P ,i2d)+amp*fluc(BCDAT_INFLOW_P), &
                     sxn,syn,szn                                   , & 
                     gv(GV_MIXT_CP  ,ijkC*indCp),gv(GV_MIXT_MOL ,ijkC*indMol), &
                     cv(CV_MIXT_DENS,ijkC)      ,cv(CV_MIXT_XMOM,ijkC)       , &
                     cv(CV_MIXT_YMOM,ijkC)      ,cv(CV_MIXT_ZMOM,ijkC)       , &
                     cv(CV_MIXT_ENER,ijkC)      ,dv(DV_MIXT_PRES,ijkC)       , &
                     cv(CV_MIXT_DENS,ijkD)      ,cv(CV_MIXT_XMOM,ijkD)       , &
                     cv(CV_MIXT_YMOM,ijkD)      ,cv(CV_MIXT_ZMOM,ijkD)       , &
                     cv(CV_MIXT_ENER,ijkD) )
#ifdef STATS
!                IF ((region%global%integrTime > eps) .AND. & 
!                    (region%global%infloNijk < NIJK_INFLOW_INIT)) THEN
!                    cv(CV_MIXT_XMOM,ijkD) = cv(CV_MIXT_XMOM,ijkD) + &
!                           (vals(BCDAT_INFLOW_U,i2d) - tav(2,ijkD)*rAvgTim)* &
!                            cv(CV_MIXT_DENS,ijkD)     
!                ENDIF
#endif                  
              ENDIF
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,&
                   __LINE__ )
            ENDIF

          ELSE   ! idum > 1
            ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
            ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

            cv(CV_MIXT_DENS,ijkD) = 2._RFREAL*cv(CV_MIXT_DENS,ijkC ) - &
                                              cv(CV_MIXT_DENS,ijkC1)
            cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_XMOM,ijkC ) - &
                                              cv(CV_MIXT_XMOM,ijkC1)
            cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_YMOM,ijkC ) - &
                                              cv(CV_MIXT_YMOM,ijkC1)
            cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_ZMOM,ijkC ) - &
                                              cv(CV_MIXT_ZMOM,ijkC1)
            cv(CV_MIXT_ENER,ijkD) = 2._RFREAL*cv(CV_MIXT_ENER,ijkC ) - &
                                              cv(CV_MIXT_ENER,ijkC1)
          ENDIF

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondInflow

!******************************************************************************
!
! Purpose: compute mfrate at injection boundary patch by aP^n method.
!
! Description: mfrate = solid-density*a*pressure^n
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: patch%mixt%vals(BCDAT_INJECT_MFRATE,:) = mfrate at current patch
!
!******************************************************************************

SUBROUTINE RFLO_BcondInjectionAPN( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, n1, n2

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iLev, lbound
  INTEGER :: idir, jdir, kdir, inode, jnode, knode
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, nOff, ijkC, ijkN, i2d
  INTEGER :: distrib

  REAL(RFREAL)          :: pb, sdens, coeff, power
  REAL(RFREAL)          :: sgn, dS, mrate, sxn, syn, szn
  REAL(RFREAL), POINTER :: dv(:,:), vals(:,:), sFace(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondInjectionAPN',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nOff    = ABS(patch%l1end-patch%l1beg) + 1
  distrib = patch%mixt%distrib

  dv   => region%levels(iLev)%mixt%dv
  vals => patch%mixt%vals

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector

  IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
  IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
  IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! loop over patch faces -------------------------------------------------------

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkC  = IndIJK(i,j,k,iCOff,ijCOff)
        ijkN  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
        dS    = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                     sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                     sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
        sxn   = sgn*sFace(XCOORD,ijkN)/dS
        syn   = sgn*sFace(YCOORD,ijkN)/dS
        szn   = sgn*sFace(ZCOORD,ijkN)/dS

        IF      (lbound==1 .OR. lbound==2) THEN
          n1 = j - jbeg
          n2 = k - kbeg
        ELSE IF (lbound==3 .OR. lbound==4) THEN
          n1 = k - kbeg
          n2 = i - ibeg
        ELSE IF (lbound==5 .OR. lbound==6) THEN
          n1 = i - ibeg
          n2 = j - jbeg
        ENDIF
        i2d   = distrib * IndIJ(n1,n2,nOff)
        sdens = vals(BCDAT_INJECT_SDENS ,i2d)
        coeff = vals(BCDAT_INJECT_ACOEFF,i2d)
        power = vals(BCDAT_INJECT_NPOWER,i2d)
        pb    = dv(DV_MIXT_PRES,ijkC)
        vals(BCDAT_INJECT_MFRATE,i2d) = sdens*coeff*(pb**power)

        mRate = vals(BCDAT_INJECT_MFRATE,i2d)
        vals(BCDAT_INJECT_RFVFU,i2d) = -mRate*sxn
        vals(BCDAT_INJECT_RFVFV,i2d) = -mRate*syn
        vals(BCDAT_INJECT_RFVFW,i2d) = -mRate*szn

      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondInjectionAPN

!******************************************************************************
!
! Purpose: update values in dummy cells at injection boundary patch.
!
! Description: the boundary condition uses linear extrapolation
!              of conservative variables. In addition, a strong boundary
!              condition is applied to the convective fluxes.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: if the mass flow rate is zero or negative, the corresponding
!        cell face is treated like solid wall.
!
!******************************************************************************

SUBROUTINE RFLO_BcondInjection( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, BcondInjectionPerf, &
        MixtureProperties
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, n1, n2

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, iLev, lbound
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, nOff, ijkC, ijkC1, ijkD, ijkN, i2d
  INTEGER :: distrib, gasModel, flowModel, dumExtrapol, bcType
  INTEGER :: indCp, indMol, inode, jnode, knode

  REAL(RFREAL)          :: sgn, rhoa, rhoua, rhova, rhowa, rhoea, pa
  REAL(RFREAL)          :: mRate, tBurn, rhoVrel(3), rgas, dS, sxn, syn, szn
  REAL(RFREAL)          :: uinj, vinj, winj, rhoDum, rhoEdum, maxChange
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), vals(:,:), sFace(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondInjection',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nOff        = ABS(patch%l1end-patch%l1beg) + 1
  bcType      = patch%bcType
  distrib     = patch%mixt%distrib
  dumExtrapol = patch%mixt%switches(BCSWI_INJECT_EXTRAP)
  maxChange   = patch%mixt%maxChange
  gasModel    = region%mixtInput%gasModel
  flowModel   = region%mixtInput%flowModel
  indCp       = region%levels(iLev)%mixt%indCp
  indMol      = region%levels(iLev)%mixt%indMol

  cv   => region%levels(iLev)%mixt%cv
  dv   => region%levels(iLev)%mixt%dv
  gv   => region%levels(iLev)%mixt%gv
  vals => patch%mixt%vals

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector

  IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
  IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
  IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! loop over patch and dummy cells ---------------------------------------------

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
          ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)
          ijkD  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
          ijkN  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          dS    = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                       sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                       sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
          sxn   = sgn*sFace(XCOORD,ijkN)/dS
          syn   = sgn*sFace(YCOORD,ijkN)/dS
          szn   = sgn*sFace(ZCOORD,ijkN)/dS

          IF      (lbound==1 .OR. lbound==2) THEN
            n1 = j - jbeg
            n2 = k - kbeg
          ELSE IF (lbound==3 .OR. lbound==4) THEN
            n1 = k - kbeg
            n2 = i - ibeg
          ELSE IF (lbound==5 .OR. lbound==6) THEN
            n1 = i - ibeg
            n2 = j - jbeg
          ENDIF
          i2d        = distrib * IndIJ(n1,n2,nOff)
          mRate      = vals(BCDAT_INJECT_MFRATE,i2d)
          tBurn      = vals(BCDAT_INJECT_TEMP  ,i2d)
          rhoVrel(1) = vals(BCDAT_INJECT_RFVFU ,i2d)
          rhoVrel(2) = vals(BCDAT_INJECT_RFVFV ,i2d)
          rhoVrel(3) = vals(BCDAT_INJECT_RFVFW ,i2d)

! ------- surface burning and 1st dummmy cell
! ------- or heat transfer wall

          IF ((mRate>0._RFREAL .AND. idum==1) .OR.  &
             (bcType == BC_INJECTION_HT)) THEN
            IF (gasModel == GAS_MODEL_TCPERF) THEN
              CALL BcondInjectionPerf( distrib,mRate,tBurn,rhoVrel,sxn,syn,szn,&
                                       gv(GV_MIXT_CP  ,ijkC*indCp ), &
                                       gv(GV_MIXT_MOL ,ijkC*indMol), &
                                       dv(DV_MIXT_PRES,ijkC       ), &
                                       rhoa,rhoua,rhova,rhowa,rhoea,pa, &
                                       uinj,vinj,winj )
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,&
                   __LINE__ )
            ENDIF

            rhoDum  = 2._RFREAL*rhoa  - cv(CV_MIXT_DENS,ijkC)
            rhoEdum = 2._RFREAL*rhoea - cv(CV_MIXT_ENER,ijkC)
            IF (rhoDum<=0._RFREAL .OR. rhoEdum<=0._RFREAL .OR. &
                ABS(dv(DV_MIXT_TEMP,ijkC)-tBurn)>maxChange*tBurn .OR. &
                dumExtrapol==EXTRAPOL_CONST) THEN
              cv(CV_MIXT_DENS,ijkD) = rhoa
              cv(CV_MIXT_XMOM,ijkD) = rhoua
              cv(CV_MIXT_YMOM,ijkD) = rhova
              cv(CV_MIXT_ZMOM,ijkD) = rhowa
              cv(CV_MIXT_ENER,ijkD) = rhoea
            ELSE
              cv(CV_MIXT_DENS,ijkD) = 2._RFREAL*rhoa  - cv(CV_MIXT_DENS,ijkC)
              cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*rhoua - cv(CV_MIXT_XMOM,ijkC)
              cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*rhova - cv(CV_MIXT_YMOM,ijkC)
              cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*rhowa - cv(CV_MIXT_ZMOM,ijkC)
              cv(CV_MIXT_ENER,ijkD) = 2._RFREAL*rhoea - cv(CV_MIXT_ENER,ijkC)
            ENDIF

! ------- not burning or dummy cell>1, inviscid flow

          ELSE IF (flowModel == FLOW_EULER) THEN
            rhoDum  = 2._RFREAL*cv(CV_MIXT_DENS,ijkC) - cv(CV_MIXT_DENS,ijkC1)
            rhoEdum = 2._RFREAL*cv(CV_MIXT_ENER,ijkC) - cv(CV_MIXT_ENER,ijkC1)
            IF (rhoDum<=0._RFREAL .OR. rhoEdum<=0._RFREAL .OR. &
                ABS(cv(CV_MIXT_DENS,ijkD)-rhoDum)> &
                maxChange*cv(CV_MIXT_DENS,ijkD) .OR. &
                ABS(cv(CV_MIXT_ENER,ijkD)-rhoEdum)> &
                maxChange*cv(CV_MIXT_ENER,ijkD) .OR. &
                dumExtrapol==EXTRAPOL_CONST) THEN
              cv(CV_MIXT_DENS,ijkD) = cv(CV_MIXT_DENS,ijkC)
              cv(CV_MIXT_XMOM,ijkD) = cv(CV_MIXT_XMOM,ijkC)
              cv(CV_MIXT_YMOM,ijkD) = cv(CV_MIXT_YMOM,ijkC)
              cv(CV_MIXT_ZMOM,ijkD) = cv(CV_MIXT_ZMOM,ijkC)
              cv(CV_MIXT_ENER,ijkD) = cv(CV_MIXT_ENER,ijkC)
            ELSE
              cv(CV_MIXT_DENS,ijkD) = 2._RFREAL*cv(CV_MIXT_DENS,ijkC ) - &
                                                cv(CV_MIXT_DENS,ijkC1)
              cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_XMOM,ijkC ) - &
                                                cv(CV_MIXT_XMOM,ijkC1)
              cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_YMOM,ijkC ) - &
                                                cv(CV_MIXT_YMOM,ijkC1)
              cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_ZMOM,ijkC ) - &
                                                cv(CV_MIXT_ZMOM,ijkC1)
              cv(CV_MIXT_ENER,ijkD) = 2._RFREAL*cv(CV_MIXT_ENER,ijkC ) - &
                                                cv(CV_MIXT_ENER,ijkC1)
            ENDIF

! ------- not burning or dummy cell>1, viscous flow
! ------- or heat transfer wall

          ELSE IF (flowModel == FLOW_NAVST) THEN
            IF ((mRate > 0._RFREAL) .OR.  &
             (bcType == BC_INJECTION_HT)) THEN  ! burning => extrapolate
              cv(CV_MIXT_DENS,ijkD) = 2._RFREAL*cv(CV_MIXT_DENS,ijkC ) - &
                                                cv(CV_MIXT_DENS,ijkC1)
              cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_XMOM,ijkC ) - &
                                                cv(CV_MIXT_XMOM,ijkC1)
              cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_YMOM,ijkC ) - &
                                                cv(CV_MIXT_YMOM,ijkC1)
              cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_ZMOM,ijkC ) - &
                                                cv(CV_MIXT_ZMOM,ijkC1)
              cv(CV_MIXT_ENER,ijkD) = 2._RFREAL*cv(CV_MIXT_ENER,ijkC ) - &
                                                cv(CV_MIXT_ENER,ijkC1)
            ELSE                         ! not burning => like noslip wall
              ijkC = IndIJK(i+(idum-1)*idir,j+(idum-1)*jdir,k+(idum-1)*kdir,iCOff,ijCOff)
              cv(CV_MIXT_DENS,ijkD) =  cv(CV_MIXT_DENS,ijkC)
              cv(CV_MIXT_XMOM,ijkD) = -cv(CV_MIXT_XMOM,ijkC)
              cv(CV_MIXT_YMOM,ijkD) = -cv(CV_MIXT_YMOM,ijkC)
              cv(CV_MIXT_ZMOM,ijkD) = -cv(CV_MIXT_ZMOM,ijkC)
              cv(CV_MIXT_ENER,ijkD) =  cv(CV_MIXT_ENER,ijkC)
            ENDIF
          ENDIF

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondInjection

!******************************************************************************
!
! Purpose: convert mdot into density*velocity for injection boundaries.
!
! Description: none.
!
! Input: .
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BcondInjectionInit( region )

  USE ModInterfaces, ONLY : RFLO_CopyBoundaryData, RFLO_GetPatchIndices, &
        RFLO_GetPatchDirection, RFLO_GetCellOffset, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, iPatch, i, j, k, n1, n2

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, lbound
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, nOff, ijkC, ijkN, i2d
  INTEGER :: inode, jnode, knode

  REAL(RFREAL)          :: sgn, dS, sxn, syn, szn, mRate
  REAL(RFREAL)          :: sdens, coeff, power, pb
  REAL(RFREAL), POINTER :: dv(:,:), vals(:,:), sFace(:,:)

  TYPE(t_patch), POINTER  :: patch, patchPrev
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_BcondInjectionInit',&
       'RFLO_ModBoundaryConditions.F90' )

! copy values/distribution to variables ---------------------------------------

  DO iLev=1,region%nGridLevels
    dv => region%levels(iLev)%mixt%dv

    DO iPatch=1,region%nPatches

      patch => region%levels(iLev)%patches(iPatch)
      IF (iLev > 1) patchPrev => region%levels(iLev-1)%patches(iPatch)

      IF (patch%bcType>=BC_INJECTION .AND. &
          patch%bcType<=BC_INJECTION+BC_RANGE) THEN

        lbound =  patch%lbound
        nOff   =  ABS(patch%l1end-patch%l1beg) + 1
        vals   => patch%mixt%vals

! ----- BC values as distribution

        IF (patch%mixt%distrib==BCDAT_DISTRIB) THEN
          CALL RFLO_GetPatchIndices( region,patch,iLev, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )
          CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
          CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
          CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! ------- to take the right face vector and make it point outwards

          sgn   = +1._RFREAL
          inode = 0
          jnode = 0
          knode = 0
          IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
            sgn   = -1._RFREAL
            inode = -idir
            jnode = -jdir
            knode = -kdir
          ENDIF

! ------- get the appropriate face vector

          IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
          IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
          IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! ------- mRate for injection APN

          IF (patch%bcType==BC_INJECTION_APN) THEN
            DO k=kbeg,kend
              DO j=jbeg,jend
                DO i=ibeg,iend
                  ijkC  = IndIJK(i,j,k,iCOff,ijCOff)

                  IF      (lbound==1 .OR. lbound==2) THEN
                    n1 = j - jbeg
                    n2 = k - kbeg
                  ELSE IF (lbound==3 .OR. lbound==4) THEN
                    n1 = k - kbeg
                    n2 = i - ibeg
                  ELSE IF (lbound==5 .OR. lbound==6) THEN
                    n1 = i - ibeg
                    n2 = j - jbeg
                  ENDIF
                  i2d   = IndIJ(n1,n2,nOff)
                  sdens = vals(BCDAT_INJECT_SDENS ,i2d)
                  coeff = vals(BCDAT_INJECT_ACOEFF,i2d)
                  power = vals(BCDAT_INJECT_NPOWER,i2d)
                  pb    = dv(DV_MIXT_PRES,ijkC)
                  vals(BCDAT_INJECT_MFRATE,i2d) = sdens*coeff*(pb**power)

                ENDDO  ! i
              ENDDO    ! j
            ENDDO      ! k
          ENDIF        ! bcType

! ------- loop over the cells

          DO k=kbeg,kend
            DO j=jbeg,jend
              DO i=ibeg,iend
                ijkC  = IndIJK(i,j,k,iCOff,ijCOff)
                ijkN  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
                dS    = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                             sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                             sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
                sxn   = sgn*sFace(XCOORD,ijkN)/dS
                syn   = sgn*sFace(YCOORD,ijkN)/dS
                szn   = sgn*sFace(ZCOORD,ijkN)/dS
                IF      (lbound==1 .OR. lbound==2) THEN
                  n1 = j - jbeg
                  n2 = k - kbeg
                ELSE IF (lbound==3 .OR. lbound==4) THEN
                  n1 = k - kbeg
                  n2 = i - ibeg
                ELSE IF (lbound==5 .OR. lbound==6) THEN
                  n1 = i - ibeg
                  n2 = j - jbeg
                ENDIF
                i2d   = IndIJ(n1,n2,nOff)
                mRate = vals(BCDAT_INJECT_MFRATE,i2d)
                vals(BCDAT_INJECT_RFVFU,i2d) = -mRate*sxn
                vals(BCDAT_INJECT_RFVFV,i2d) = -mRate*syn
                vals(BCDAT_INJECT_RFVFW,i2d) = -mRate*szn
              ENDDO  ! i
            ENDDO    ! j
          ENDDO      ! k

! ----- BC values constant

        ELSE
          vals(BCDAT_INJECT_RFVFU,:) = 0._RFREAL
          vals(BCDAT_INJECT_RFVFV,:) = 0._RFREAL
          vals(BCDAT_INJECT_RFVFW,:) = 0._RFREAL
        ENDIF  ! distrib

        IF (iLev > 1) THEN
          CALL RFLO_CopyBoundaryData( global,patchPrev,patch )
        ENDIF

      ENDIF    ! bcType

    ENDDO      ! iPatch
  ENDDO        ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_BcondInjectionInit

!******************************************************************************
!
! Purpose: update values in dummy cells at noslip (viscous) walls.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BcondNoslipWall( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, MixtureProperties, MixtPerf_R_M, MixtPerf_G_CpR, &
        MixtPerf_D_PRT, MixtPerf_Eo_DGPUVW, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, n1, n2

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, lbound, iCOff, ijCOff, nOff, distrib, bcOpt
  INTEGER :: indCp, indMol, gasModel, ijkC, ijkC1, ijkCi, ijkCb, ijkD, i2d
  INTEGER :: inode, jnode, knode, iNOff, ijNOff, ijkN(4)

  LOGICAL :: moveGrid

  REAL(RFREAL)          :: tWall, rgas, gamma, pb, u, v, w, temp, rho, rhoe
  REAL(RFREAL)          :: dt, dxn(4), dyn(4), dzn(4), uSurf, vSurf, wSurf
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), vals(:,:)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondNoslipWall',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib
  bcOpt     = patch%mixt%switches(BCSWI_NOSLIP_ADIABAT)
  indCp     = region%levels(iLev)%mixt%indCp
  indMol    = region%levels(iLev)%mixt%indMol
  gasModel = region%mixtInput%gasModel
  moveGrid  = region%mixtInput%moveGrid
  dt        = region%global%dtMin

  cv     => region%levels(iLev)%mixt%cv
  dv     => region%levels(iLev)%mixt%dv
  gv     => region%levels(iLev)%mixt%gv
  xyz    => region%levels(iLev)%grid%xyz
  IF (moveGrid) xyzOld => region%levels(iLev)%gridOld%xyz
  vals   => patch%mixt%vals

! to take the correct surface coordinates

  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! loop over all cells of a patch ----------------------------------------------

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkD  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
          ijkCI = IndIJK(i+(idum-1)*idir,j+(idum-1)*jdir,k+(idum-1)*kdir,iCOff,ijCOff)

! ------- wall velocity in case of moving surface

          IF (moveGrid) THEN
            ijkN(1) = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            IF      (lbound==1 .OR. lbound==2) THEN
              ijkN(2) = IndIJK(i+inode,j+jnode+1,k+knode  ,iNOff,ijNOff)
              ijkN(3) = IndIJK(i+inode,j+jnode+1,k+knode+1,iNOff,ijNOff)
              ijkN(4) = IndIJK(i+inode,j+jnode  ,k+knode+1,iNOff,ijNOff)
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              ijkN(2) = IndIJK(i+inode  ,j+jnode,k+knode+1,iNOff,ijNOff)
              ijkN(3) = IndIJK(i+inode+1,j+jnode,k+knode+1,iNOff,ijNOff)
              ijkN(4) = IndIJK(i+inode+1,j+jnode,k+knode  ,iNOff,ijNOff)
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              ijkN(2) = IndIJK(i+inode+1,j+jnode  ,k+knode,iNOff,ijNOff)
              ijkN(3) = IndIJK(i+inode+1,j+jnode+1,k+knode,iNOff,ijNOff)
              ijkN(4) = IndIJK(i+inode  ,j+jnode+1,k+knode,iNOff,ijNOff)
            ENDIF
            dxn(1) = xyz(XCOORD,ijkN(1)) - xyzOld(XCOORD,ijkN(1))
            dyn(1) = xyz(YCOORD,ijkN(1)) - xyzOld(YCOORD,ijkN(1))
            dzn(1) = xyz(ZCOORD,ijkN(1)) - xyzOld(ZCOORD,ijkN(1))
            dxn(2) = xyz(XCOORD,ijkN(2)) - xyzOld(XCOORD,ijkN(2))
            dyn(2) = xyz(YCOORD,ijkN(2)) - xyzOld(YCOORD,ijkN(2))
            dzn(2) = xyz(ZCOORD,ijkN(2)) - xyzOld(ZCOORD,ijkN(2))
            dxn(3) = xyz(XCOORD,ijkN(3)) - xyzOld(XCOORD,ijkN(3))
            dyn(3) = xyz(YCOORD,ijkN(3)) - xyzOld(YCOORD,ijkN(3))
            dzn(3) = xyz(ZCOORD,ijkN(3)) - xyzOld(ZCOORD,ijkN(3))
            dxn(4) = xyz(XCOORD,ijkN(4)) - xyzOld(XCOORD,ijkN(4))
            dyn(4) = xyz(YCOORD,ijkN(4)) - xyzOld(YCOORD,ijkN(4))
            dzn(4) = xyz(ZCOORD,ijkN(4)) - xyzOld(ZCOORD,ijkN(4))
            uSurf  = 0.25_RFREAL*(dxn(1)+dxn(2)+dxn(3)+dxn(4))/dt
            vSurf  = 0.25_RFREAL*(dyn(1)+dyn(2)+dyn(3)+dyn(4))/dt
            wSurf  = 0.25_RFREAL*(dzn(1)+dzn(2)+dzn(3)+dzn(4))/dt
          ELSE
            uSurf = 0._RFREAL
            vSurf = 0._RFREAL
            wSurf = 0._RFREAL
          ENDIF

! ------- adiabatic wall

          IF (bcOpt == BCOPT_ADIABAT) THEN
            cv(CV_MIXT_DENS,ijkD) = cv(CV_MIXT_DENS,ijkCI)
            cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*uSurf - cv(CV_MIXT_XMOM,ijkCI)
            cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*vSurf - cv(CV_MIXT_YMOM,ijkCI)
            cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*wSurf - cv(CV_MIXT_ZMOM,ijkCI)
            cv(CV_MIXT_ENER,ijkD) = cv(CV_MIXT_ENER,ijkCI)

! ------- prescribed wall temperature

          ELSE
            ijkCB = IndIJK(i,j,k,iCOff,ijCOff)
            ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
            ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
            ENDIF
            i2d   = distrib * IndIJ(n1,n2,nOff)

            tWall = vals(BCDAT_NOSLIP_TWALL,i2d)
            pb    = dv(DV_MIXT_PRES,ijkCB)
            u     = cv(CV_MIXT_XMOM,ijkCI)/cv(CV_MIXT_DENS,ijkCI)
            v     = cv(CV_MIXT_YMOM,ijkCI)/cv(CV_MIXT_DENS,ijkCI)
            w     = cv(CV_MIXT_ZMOM,ijkCI)/cv(CV_MIXT_DENS,ijkCI)
            IF (idum == 1) THEN
              temp = 2._RFREAL*tWall - dv(DV_MIXT_TEMP,ijkC)
            ELSE
              temp = 2._RFREAL*dv(DV_MIXT_TEMP,ijkC) - dv(DV_MIXT_TEMP,ijkC1)
            ENDIF

            IF (gasModel == GAS_MODEL_TCPERF) THEN
              rgas  = MixtPerf_R_M( gv(GV_MIXT_MOL,ijkCB*indMol) )
              gamma = MixtPerf_G_CpR( gv(GV_MIXT_CP,ijkCB*indCp),rgas )
              rho   = MixtPerf_D_PRT( pb,rgas,temp )
              rhoe  = rho * MixtPerf_Eo_DGPUVW( rho,gamma,pb,u,v,w )
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,&
                   __LINE__ )
            ENDIF
            cv(CV_MIXT_DENS,ijkD) = rho
            cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*uSurf - u*rho
            cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*vSurf - v*rho
            cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*wSurf - w*rho
            cv(CV_MIXT_ENER,ijkD) = rhoe
          ENDIF  ! bcOpt

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondNoslipWall

!******************************************************************************
!
! Purpose: update values in dummy cells at outflow boundary patch.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BcondOutflow( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, BcondOutflowPerf, &
        MixtureProperties
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, n1, n2

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, lbound, iCOff, ijCOff, iNOff, ijNOff, nOff, distrib, bcOpt
  INTEGER :: indCp, indMol, gasModel, ijkC, ijkC1, ijkD, ijkN, i2d
  INTEGER :: inode, jnode, knode
  INTEGER :: bcModel

  REAL(RFREAL) :: sgn, ds, sxn, syn, szn, pout
  REAL(RFREAL) :: ud, vd, wd, tempd, rOld, rrOld, uOld, vOld, wOld, &
                  rdOld, rrdOld, udOld, vdOld, wdOld, edOld, rLen, dtMin, kappa
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), vals(:,:), sFace(:,:)
  REAL(RFREAL), POINTER :: cvOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondOutflow',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib
  bcOpt     = patch%mixt%switches(BCSWI_OUTFLOW_TYPE)
  bcModel   = patch%mixt%switches(BCSWI_OUTFLOW_MODEL)
  indCp     = region%levels(iLev)%mixt%indCp
  indMol    = region%levels(iLev)%mixt%indMol
  gasModel = region%mixtInput%gasModel

  cv    => region%levels(iLev)%mixt%cv
  cvOld => region%levels(iLev)%mixt%cvOld
  dv    => region%levels(iLev)%mixt%dv
  gv    => region%levels(iLev)%mixt%gv
  vals  => patch%mixt%vals

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector

  IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
  IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
  IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)

          IF (idum == 1) THEN
            ijkC = IndIJK(i,j,k,iCOff,ijCOff)
            ijkN = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            ds   = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                        sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                        sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
            sxn  = sgn*sFace(XCOORD,ijkN)/ds
            syn  = sgn*sFace(YCOORD,ijkN)/ds
            szn  = sgn*sFace(ZCOORD,ijkN)/ds

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
            ENDIF
            i2d = distrib * IndIJ(n1,n2,nOff)

            IF (gasModel == GAS_MODEL_TCPERF) THEN
              IF (bcOpt == BCOPT_SUPERSONIC) THEN   ! cannot specify pressure
                pout = 0._RFREAL
              ELSE                                  ! take given pressure
                pout = vals(BCDAT_OUTFLOW_PRESS,i2d)
              ENDIF

              IF (bcModel == BCOPT_DEFAULT) THEN
                CALL BcondOutflowPerf( bcOpt,pout,sxn,syn,szn, &
                     gv(GV_MIXT_CP  ,ijkC*indCp),gv(GV_MIXT_MOL ,ijkC*indMol), &
                     cv(CV_MIXT_DENS,ijkC)      ,cv(CV_MIXT_XMOM,ijkC)       , &
                     cv(CV_MIXT_YMOM,ijkC)      ,cv(CV_MIXT_ZMOM,ijkC)       , &
                     cv(CV_MIXT_ENER,ijkC)      ,dv(DV_MIXT_PRES,ijkC)       , &
                     cv(CV_MIXT_DENS,ijkD)      ,cv(CV_MIXT_XMOM,ijkD)       , &
                     cv(CV_MIXT_YMOM,ijkD)      ,cv(CV_MIXT_ZMOM,ijkD)       , &
                     cv(CV_MIXT_ENER,ijkD) )
              ELSE
                ijkC1 = &
                IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

                rOld  = cvOld(CV_MIXT_DENS,ijkC )
                rrOld = 1._RFREAL/rOld
                uOld  = cvOld(CV_MIXT_XMOM,ijkC )*rrOld
                vOld  = cvOld(CV_MIXT_YMOM,ijkC )*rrOld
                wOld  = cvOld(CV_MIXT_ZMOM,ijkC )*rrOld

                rdOld = cvOld(CV_MIXT_DENS,ijkD )
                rrdOld= 1._RFREAL/rdOld
                udOld = cvOld(CV_MIXT_XMOM,ijkD )*rrdOld
                vdOld = cvOld(CV_MIXT_YMOM,ijkD )*rrdOld
                wdOld = cvOld(CV_MIXT_ZMOM,ijkD )*rrdOld
                edOld = cvOld(CV_MIXT_ENER,ijkD )*rrdOld 
                rLen  = region%global%refLength
                dtMin = region%global%dtMin
                IF (bcOpt == BCOPT_SUPERSONIC) THEN   ! cannot specify nrcoef
                  kappa = 0._RFREAL
                ELSE
                  kappa = vals(BCDAT_OUTFLOW_NRCOEF,i2d)
                ENDIF

                ud    = 2._RFREAL*dv(DV_MIXT_UVEL,ijkC ) - dv(DV_MIXT_UVEL,ijkC1)
                vd    = 2._RFREAL*dv(DV_MIXT_VVEL,ijkC ) - dv(DV_MIXT_VVEL,ijkC1)
                wd    = 2._RFREAL*dv(DV_MIXT_WVEL,ijkC ) - dv(DV_MIXT_WVEL,ijkC1)
                tempd = 2._RFREAL*dv(DV_MIXT_TEMP,ijkC ) - dv(DV_MIXT_TEMP,ijkC1)

                CALL RFLO_BcondOutflowPerf( bcOpt,pout,sxn,syn,szn, &
                     gv(GV_MIXT_CP  ,ijkC*indCp),gv(GV_MIXT_MOL ,ijkC*indMol), &
                     cv(CV_MIXT_DENS,ijkC)      ,cv(CV_MIXT_XMOM,ijkC)       , &
                     cv(CV_MIXT_YMOM,ijkC)      ,cv(CV_MIXT_ZMOM,ijkC)       , &
                     cv(CV_MIXT_ENER,ijkC)      ,dv(DV_MIXT_PRES,ijkC)       , &
                     ud, vd, wd, tempd, uOld, vOld, wOld                     , &
                     udOld, vdOld, wdOld, rdOld, edOld, rLen, dtMin, kappa   , &
                     cv(CV_MIXT_DENS,ijkD)      ,cv(CV_MIXT_XMOM,ijkD)       , &
                     cv(CV_MIXT_YMOM,ijkD)      ,cv(CV_MIXT_ZMOM,ijkD)       , &
                     cv(CV_MIXT_ENER,ijkD) )
              ENDIF  ! bcModel
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,&
                   __LINE__ )
            ENDIF  ! gasModel

          ELSE   ! idum > 1
            ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
            ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

            cv(CV_MIXT_DENS,ijkD) = 2._RFREAL*cv(CV_MIXT_DENS,ijkC ) - &
                                              cv(CV_MIXT_DENS,ijkC1)
            cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_XMOM,ijkC ) - &
                                              cv(CV_MIXT_XMOM,ijkC1)
            cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_YMOM,ijkC ) - &
                                              cv(CV_MIXT_YMOM,ijkC1)
            cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_ZMOM,ijkC ) - &
                                              cv(CV_MIXT_ZMOM,ijkC1)
            cv(CV_MIXT_ENER,ijkD) = 2._RFREAL*cv(CV_MIXT_ENER,ijkC ) - &
                                              cv(CV_MIXT_ENER,ijkC1)
          ENDIF

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondOutflow

!******************************************************************************
!
! Purpose: exchange values between interior and dummy cells of the
!          corresponding patches of a rotationally periodic boundary.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BcondRotatPeriod( region,regionSrc,patch,patchSrc )

  USE ModInterfaces, ONLY : MixtureProperties
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables


! ... local variables

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondRotatPeriod',&
       'RFLO_ModBoundaryConditions.F90' )



  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondRotatPeriod

!******************************************************************************
!
! Purpose: update values in dummy cells at slip walls.
!
! Description: the boundary condition uses extrapolation of conservative
!              variables (constant or linear). In addition, a strong
!              boundary condition is applied to the convective fluxes.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BcondSlipWall( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, MixtureProperties
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, iCOff, ijCOff, ijkC, ijkC1, ijkD
  INTEGER :: dumExtrapol, gasModel

  REAL(RFREAL)          :: rhoDum, rhoEdum, maxChange
  REAL(RFREAL)          :: uvel, vvel, wvel
  REAL(RFREAL), POINTER :: cv(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondSlipWall',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  dumExtrapol = patch%mixt%switches(BCSWI_SLIPW_EXTRAP)
  maxChange   = patch%mixt%maxChange
  gasModel   = region%mixtInput%gasModel

  cv => region%levels(iLev)%mixt%cv

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
          ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)
          ijkD  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)

          rhoDum  = 2._RFREAL*cv(CV_MIXT_DENS,ijkC) - cv(CV_MIXT_DENS,ijkC1)
          rhoEdum = 2._RFREAL*cv(CV_MIXT_ENER,ijkC) - cv(CV_MIXT_ENER,ijkC1)
          IF (rhoDum<=0._RFREAL .OR. rhoEdum<=0._RFREAL .OR. &
              ABS(cv(CV_MIXT_DENS,ijkD)-rhoDum)> &
              maxChange*cv(CV_MIXT_DENS,ijkD) .OR. &
              ABS(cv(CV_MIXT_ENER,ijkD)-rhoEdum)> &
              maxChange*cv(CV_MIXT_ENER,ijkD) .OR. &
              dumExtrapol==EXTRAPOL_CONST) THEN
            cv(CV_MIXT_DENS,ijkD) = cv(CV_MIXT_DENS,ijkC)
            cv(CV_MIXT_XMOM,ijkD) = cv(CV_MIXT_XMOM,ijkC)
            cv(CV_MIXT_YMOM,ijkD) = cv(CV_MIXT_YMOM,ijkC)
            cv(CV_MIXT_ZMOM,ijkD) = cv(CV_MIXT_ZMOM,ijkC)
            cv(CV_MIXT_ENER,ijkD) = cv(CV_MIXT_ENER,ijkC)
          ELSE
            cv(CV_MIXT_DENS,ijkD) = 2._RFREAL*cv(CV_MIXT_DENS,ijkC ) - &
                                              cv(CV_MIXT_DENS,ijkC1)
            cv(CV_MIXT_XMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_XMOM,ijkC ) - &
                                              cv(CV_MIXT_XMOM,ijkC1)
            cv(CV_MIXT_YMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_YMOM,ijkC ) - &
                                              cv(CV_MIXT_YMOM,ijkC1)
            cv(CV_MIXT_ZMOM,ijkD) = 2._RFREAL*cv(CV_MIXT_ZMOM,ijkC ) - &
                                              cv(CV_MIXT_ZMOM,ijkC1)
            cv(CV_MIXT_ENER,ijkD) = 2._RFREAL*cv(CV_MIXT_ENER,ijkC ) - &
                                              cv(CV_MIXT_ENER,ijkC1)
            IF (patch%mixt%setMotion .AND. idum==1) THEN
              uvel = 2._RFREAL*patch%mixt%bndVel(XCOORD) - & 
                        cv(CV_MIXT_XMOM,ijkC)/cv(CV_MIXT_DENS,ijkC)
              vvel = 2._RFREAL*patch%mixt%bndVel(YCOORD) - & 
                        cv(CV_MIXT_YMOM,ijkC)/cv(CV_MIXT_DENS,ijkC)
              wvel = 2._RFREAL*patch%mixt%bndVel(ZCOORD) - & 
                        cv(CV_MIXT_ZMOM,ijkC)/cv(CV_MIXT_DENS,ijkC)
              cv(CV_MIXT_XMOM,ijkD) = uvel*cv(CV_MIXT_DENS,ijkD )
              cv(CV_MIXT_YMOM,ijkD) = vvel*cv(CV_MIXT_DENS,ijkD )
              cv(CV_MIXT_ZMOM,ijkD) = wvel*cv(CV_MIXT_DENS,ijkD )
            ENDIF
          ENDIF

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondSlipWall

!******************************************************************************
!
! Purpose: update values in dummy cells at symmetry boundary patch.
!
! Description: (1) density and energy in dummy cells are set equal to values
!                  in the corresponding interior nodes
!
!              (2) velocity components in dummy cells are mirrored with respect
!                  to the boundary; it is assumed that the boundary coincides
!                  either with a x=const., y=const., or z=const. plane.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BcondSymmetry( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, MixtureProperties
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, lbound, iCOff, ijCOff, ijkC, ijkD
  INTEGER :: gasModel

  REAL(RFREAL)          :: sgn(3)
  REAL(RFREAL), POINTER :: cv(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BcondSymmetry',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions and pointers

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  gasModel = region%mixtInput%gasModel

  cv => region%levels(iLev)%mixt%cv

! to mirror the appropriate velocity components

  IF      (lbound==1 .OR. lbound==2) THEN
    sgn(1) = -1._RFREAL
    sgn(2) = +1._RFREAL
    sgn(3) = +1._RFREAL
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    sgn(1) = +1._RFREAL
    sgn(2) = -1._RFREAL
    sgn(3) = +1._RFREAL
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    sgn(1) = +1._RFREAL
    sgn(2) = +1._RFREAL
    sgn(3) = -1._RFREAL
  ENDIF

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkC = IndIJK(i+(idum-1)*idir,j+(idum-1)*jdir,k+(idum-1)*kdir,iCOff,ijCOff)
          ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)

          cv(CV_MIXT_DENS,ijkD) =        cv(CV_MIXT_DENS,ijkC)
          cv(CV_MIXT_XMOM,ijkD) = sgn(1)*cv(CV_MIXT_XMOM,ijkC)
          cv(CV_MIXT_YMOM,ijkD) = sgn(2)*cv(CV_MIXT_YMOM,ijkC)
          cv(CV_MIXT_ZMOM,ijkD) = sgn(3)*cv(CV_MIXT_ZMOM,ijkC)
          cv(CV_MIXT_ENER,ijkD) =        cv(CV_MIXT_ENER,ijkC)

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BcondSymmetry

!******************************************************************************
!
! Purpose: receive data for dummy cells from adjacent regions being
!          on another processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%mixt = updated flow values (cv,dv,tv,gv)
!                                     in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BoundaryConditionsRecv( regions,iReg )

  USE ModInterfaces, ONLY : RFLO_ReceiveDummyVals, RFLO_SetCornerEdgeCells, &
        RFLO_ExchangeCornerEdgeCells, RFLO_ReceiveCornerEdgeCells, &
        RFLO_CorrectCornerEdgeCells
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  INTEGER :: iLev, nPatches, bcType, iRegSrc, iPatchSrc

  TYPE(t_patch), POINTER :: patch, patchSrc

!******************************************************************************

  CALL RegisterFunction( regions(iReg)%global,'RFLO_BoundaryConditionsRecv', &
                         'RFLO_ModBoundaryConditions.F90' )

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! receive data (regular cells) from other processors --------------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patch%bcType
    iRegSrc   = patch%srcRegion
    iPatchSrc = patch%srcPatch

! - region interface, periodic boundary

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
        (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
        (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
      IF (regions(iRegSrc)%procid /= regions(iReg)%global%myProcid) THEN
        patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

        CALL RFLO_ReceiveDummyVals( regions(iReg),regions(iRegSrc), &
                                    patch,patchSrc )
      ENDIF
    ENDIF

  ENDDO   ! iPatch

! copy/exchange data for edge and corner cells --------------------------------

  CALL RFLO_SetCornerEdgeCells( regions(iReg) )

  DO iPatch=1,nPatches
    patch => regions(iReg)%levels(iLev)%patches(iPatch)
    bcType = patch%bcType
    IF ((bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
        (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE .AND. &
         regions(iReg)%mixtInput%flowModel==FLOW_NAVST) .OR. &
        (bcType>=BC_SYMMETRY   .AND. bcType<=BC_SYMMETRY  +BC_RANGE)) THEN
      CALL RFLO_CorrectCornerEdgeCells( regions(iReg),patch,bcType )
    ENDIF
  ENDDO   ! iPatch

  CALL RFLO_ExchangeCornerEdgeCells( regions,iReg )
  CALL RFLO_ReceiveCornerEdgeCells( regions,iReg )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( regions(iReg)%global )

END SUBROUTINE RFLO_BoundaryConditionsRecv

!******************************************************************************
!
! Purpose: send data to dummy cells of adjacent regions being located
!          on another processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: data to other processors.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BoundaryConditionsSend( regions,iReg )

  USE ModInterfaces, ONLY : RFLO_SendDummyConf, RFLO_SendDummyInt, &
        RFLO_SendDummyIreg, RFLO_SendCornerEdgeCells
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  INTEGER :: iLev, nPatches, bcType, iRegSrc, iPatchSrc

  TYPE(t_patch), POINTER  :: patch, patchSrc
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_BoundaryConditionsSend',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! send data for edge and corner cells -----------------------------------------

  CALL RFLO_SendCornerEdgeCells( regions,iReg )

! loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patch%bcType
    iRegSrc   = patch%srcRegion
    iPatchSrc = patch%srcPatch

! - conforming region interface

    IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
        CALL RFLO_SendDummyConf( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - non-conforming region interface (integer)

    ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
        CALL RFLO_SendDummyInt( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - non-conforming region interface (irregular)

    ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
        CALL RFLO_SendDummyIreg( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - translational periodicity

    ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid /= global%myProcid) THEN
        CALL RFLO_SendDummyConf( regions(iReg),regions(iRegSrc),patch )
      ENDIF

! - rotational periodicity

    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      CALL RFLO_BcondRotatPeriod( regions(iReg),regions(iRegSrc), &
                                  patch,patchSrc )

    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_BoundaryConditionsSend

!******************************************************************************
!
! Purpose: set boundary conditions or exchange data between adjacent
!          regions being on same processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%mixt = updated flow values (cv,dv,tv,gv)
!                                     in dummy cells of current region.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_BoundaryConditionsSet( regions,iReg )

  USE ModInterfaces, ONLY : RFLO_ExchangeDummyConf, RFLO_ExchangeDummyInt, &
        RFLO_ExchangeDummyIreg
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  INTEGER :: iLev, nPatches, bcType, iRegSrc, iPatchSrc

  TYPE(t_patch), POINTER  :: patch, patchSrc
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_BoundaryConditionsSet',&
       'RFLO_ModBoundaryConditions.F90' )

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patch%bcType
    iRegSrc   = patch%srcRegion
    iPatchSrc = patch%srcPatch

! - inflow

    IF (bcType>=BC_INFLOW .AND. bcType<=BC_INFLOW+BC_RANGE) THEN
      CALL RFLO_BcondInflow( regions(iReg),patch )

! - outflow

    ELSE IF (bcType>=BC_OUTFLOW .AND. bcType<=BC_OUTFLOW+BC_RANGE) THEN
      CALL RFLO_BcondOutflow( regions(iReg),patch )

! - conforming region interface

    ELSE IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid == global%myProcid) THEN
        CALL RFLO_ExchangeDummyConf( regions(iReg),regions(iRegSrc), &
                                     patch,patchSrc )
      ENDIF

! - non-conforming region interface (integer)

    ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid == global%myProcid) THEN
        CALL RFLO_ExchangeDummyInt( regions(iReg),regions(iRegSrc), &
                                    patch,patchSrc )
      ENDIF

! - non-conforming region interface (irregular)

    ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid == global%myProcid) THEN
        CALL RFLO_ExchangeDummyIreg( regions(iReg),regions(iRegSrc), &
                                     patch,patchSrc )
      ENDIF

! - slip wall

    ELSE IF (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) THEN
      CALL RFLO_BcondSlipWall( regions(iReg),patch )

! - noslip wall

    ELSE IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
      CALL RFLO_BcondNoslipWall( regions(iReg),patch )

! - far field

    ELSE IF (bcType>=BC_FARFIELD .AND. bcType<=BC_FARFIELD+BC_RANGE) THEN
      CALL RFLO_BcondFarfield( regions(iReg),patch )

! - injection

    ELSE IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
      IF (bcType==BC_INJECTION_APN) THEN 
        IF (patch%bcCoupled==BC_EXTERNAL) THEN
          CALL ErrorStop( global,ERR_EXTERNAL_FUNCT,&
               __LINE__, &
            'BC_INJECT_APN used in .bc with interaction flag in .top file' )
        ELSE
          CALL RFLO_BcondInjectionAPN( regions(iReg),patch )
        ENDIF
      ENDIF
      CALL RFLO_BcondInjection( regions(iReg),patch )

! - symmetry

    ELSE IF (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE) THEN
      CALL RFLO_BcondSymmetry( regions(iReg),patch )

! - translational periodicity

    ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      IF (regions(iRegSrc)%procid == global%myProcid) THEN
        CALL RFLO_ExchangeDummyConf( regions(iReg),regions(iRegSrc), &
                                     patch,patchSrc )
      ENDIF

! - rotational periodicity

    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

      CALL RFLO_BcondRotatPeriod( regions(iReg),regions(iRegSrc), &
                                  patch,patchSrc )

! - I dunno ...

    ELSE
      CALL ErrorStop( global,ERR_UNKNOWN_BC,&
           __LINE__ )
    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_BoundaryConditionsSet

!******************************************************************************
!
! Purpose: set outflow boundary condition for one cell.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              method of Rudy and Strikwerda: Boundary Conditions for Subsonic
!              Compressible Navier-Stokes Calculations.
!              Computers and Fluids 9, 327-338 (1981). The supersonic boundary 
!              condition consists of simple extrapolation.
!
! Input: bcOpt    = boundary treatment: subsonic, supersonic, or mixed
!        pout     = given static outlet pressure
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!
! Notes: this condition is valid only for thermally and calorically
!        perfect gas (supersonic outflow valid for all gases).
!
!******************************************************************************
!
! $Id: RFLO_ModBoundaryConditions.F90,v 1.30 2010/02/18 21:47:39 juzhang Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_BcondOutflowPerf( bcOpt,pout,sxn,syn,szn,cpgas,mol, & 
                                  rho,rhou,rhov,rhow,rhoe,press, &
                                  ud,vd,wd,tempd,uOld,vOld,wOld, udOld,vdOld, &
                                  wdOld,rdOld,edOld, refLength, dtMin, kappa, &
                                  rhob,rhoub,rhovb,rhowb,rhoeb )

  USE ModDataTypes
  USE ModParameters
  USE ModInterfaces, ONLY: MixtPerf_C_DGP, MixtPerf_Eo_DGPUVW, & 
                           MixtPerf_G_CpR, MixtPerf_R_M, MixtPerf_P_DEoGVm2

  IMPLICIT NONE

! ... parameters
  INTEGER :: bcOpt

  REAL(RFREAL) :: pout
  REAL(RFREAL) :: rho, rhou, rhov, rhow, rhoe, press
  REAL(RFREAL) :: ud, vd, wd, tempd, uOld, vOld, wOld
  REAL(RFREAL) :: udOld, vdOld, wdOld, rdOld, edOld, refLength, dtMin, kappa
  REAL(RFREAL) :: sxn, syn, szn, csound, cpgas, mol
  REAL(RFREAL) :: rhob, rhoub, rhovb, rhowb, rhoeb

! ... local variables
  REAL(RFREAL) :: rgas, gamma, gam1, u, v, w, mach, rhoc, deltp, &
                  ub, vb, wb, vnd
  REAL(RFREAL) :: pdOld, qdOld2, qnOld, qn, pb, beta

!******************************************************************************
! gas properties; velocity components; Mach number

  rgas  = MixtPerf_R_M( mol )
  gamma = MixtPerf_G_CpR( cpgas,rgas )
  gam1  = gamma - 1.0_RFREAL

  u      = rhou/rho
  v      = rhov/rho
  w      = rhow/rho
  csound = MixtPerf_C_DGP( rho,gamma,press )
  mach   = SQRT(u*u+v*v+w*w)/csound

! subsonic flow ---------------------------------------------------------------

  IF (mach < 1.0_RFREAL .AND. &
      (bcOpt == BCOPT_SUBSONIC .OR. bcOpt == BCOPT_MIXED)) THEN

    beta  = kappa*csound*(1._RFREAL - mach*mach)/refLength
    rhoc  = rho*csound
    qdOld2= udOld*udOld + vdOld*vdOld + wdOld*wdOld
    pdOld = MixtPerf_P_DEoGVm2( rdOld,edOld,gamma,qdOld2 )
    deltp = pdOld - pout
    qn    = u*sxn + v*syn + w*szn
    qnOld = uOld*sxn + vOld*syn + wOld*szn
!    qn    = ud*sxn + vd*syn + wd*szn
!    qnOld = udOld*sxn + vdOld*syn + wdOld*szn
!    qn    = 0.5_RFREAL*((ud+u)*sxn + (vd+v)*syn + (wd+w)*szn)
!    qnOld = 0.5_RFREAL*((udOld+uOld)*sxn + (vdOld+vOld)*syn + (wdOld+wOld)*szn)

    pb    = press + rhoc*(qn - qnOld) - dtMin*beta*(press - pout)
!    pb    = pdOld + rhoc*(qn - qnOld) - dtMin*beta*(pdOld - pout)
!    pb    = pout
    ub    = ud
    vb    = vd
    wb    = wd
    rhob  = pb/(rgas*tempd)

! - special treatment to prevent "deltp" from changing the sign
!   of velocity components. This may happen for very small u, v, w.

    vnd = ub*sxn + vb*syn + wb*szn
    IF ( vnd < 0.0_RFREAL ) THEN ! inflow at outflow boundary
      ub = 0.010_RFREAL*SIGN(1.0_RFREAL,u)*MIN(ABS(ub),ABS(u))
      vb = 0.010_RFREAL*SIGN(1.0_RFREAL,v)*MIN(ABS(vb),ABS(v))
      wb = 0.010_RFREAL*SIGN(1.0_RFREAL,w)*MIN(ABS(wb),ABS(w))
    END IF ! vnd

    rhoub = rhob*ub
    rhovb = rhob*vb
    rhowb = rhob*wb
    rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pb,ub,vb,wb )

! supersonic flow -------------------------------------------------------------

  ELSE
    rhob  = rho
    rhoub = rhou
    rhovb = rhov
    rhowb = rhow
    rhoeb = rhoe
  END IF ! mach

END SUBROUTINE RFLO_BcondOutflowPerf

!******************************************************************************
!
! Purpose: set inflow boundary condition for one cell based on presribed 
!          velocity and temperature.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              method of Rudy and Strikwerda: Boundary Conditions for Subsonic
!              Compressible Navier-Stokes Calculations.
!              Computers and Fluids 9, 327-338 (1981). The supersonic boundary 
!              condition consists of simple extrapolation.
!
! Input: bcOpt      = inflow bc-option: steady or unsteady
!        uInfl      = prescribed u-velocity
!        vInfl      = prescribed v-velocity
!        wInfl      = prescribed w-velocity
!        tInfl      = prescribed temperature
!        pInfl      = prescribed pressure (only used if supersonic)
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!
! Notes: this condition is valid only for a thermally and calorically
!        perfect gas (supersonic in/outflow valid for all gases).
!
!******************************************************************************

SUBROUTINE RFLO_BcondInflowVTPerf( model,bcOpt,uInfl,vInfl,wInfl,tInfl,pInfl, &
                                   sxn,syn,szn,cpgas,mol, &
                                   rho,rhou,rhov,rhow,rhoe,press, &
                                   rhob,rhoub,rhovb,rhowb,rhoeb )

  USE ModInterfaces, ONLY: MixtPerf_C_DGP, MixtPerf_C_GRT, MixtPerf_D_PRT, &
                           MixtPerf_G_CpR, MixtPerf_R_M, MixtPerf_Eo_DGPUVW

  IMPLICIT NONE

! ... parameters
  INTEGER :: model, bcOpt
  REAL(RFREAL) :: uInfl, vInfl, wInfl, tInfl, pInfl
  REAL(RFREAL) :: sxn, syn, szn, cpgas, mol
  REAL(RFREAL) :: rho, rhou, rhov, rhow, rhoe, press
  REAL(RFREAL) :: rhob, rhoub, rhovb, rhowb, rhoeb

! ... local variables
  INTEGER :: ierror
  REAL(RFREAL) :: rgas, gamma, mach
  REAL(RFREAL) :: re, ue, ve, we, pe, qn, crho0, csound
  REAL(RFREAL) :: ub, vb, wb, tb, pb

!******************************************************************************
! gas properties

  rgas  = MixtPerf_R_M( mol )
  gamma = MixtPerf_G_CpR( cpgas,rgas )

! flow values at a reference location (= interior cell)

  re  = rho
  ue  = rhou/rho
  ve  = rhov/rho
  we  = rhow/rho
  pe  = press

! flow values at inlet plane

  tb  = tInfl
  ub  = uInfl
  vb  = vInfl
  wb  = wInfl
  qn  = sxn*ub + syn*vb + szn*wb
  mach= qn/MixtPerf_C_GRT( gamma,rgas,tb )

  IF (qn > 0 .AND. model==BCOPT_STEADY) THEN
    WRITE(STDERR,*)'BC_INFLOW_VELTEMP: outflow detected at inflow plane'
#ifdef MPI
    CALL MPI_Abort( ierror )
#endif
    STOP
  ENDIF

  IF (bcOpt == BCOPT_MIXED) THEN

! - subsonic inflow (qn<0)

    IF (mach < 1._RFREAL) THEN
      csound = MixtPerf_C_DGP( re,gamma,pe )
      crho0  = csound*re

      pb  = pe-crho0*(sxn*(ub-ue)+syn*(vb-ve)+szn*(wb-we))

! - supersonic inflow (qn<0)

    ELSE

      pb    = pInfl
    ENDIF

  ELSEIF (bcOpt == BCOPT_SUBSONIC) THEN
    csound = MixtPerf_C_DGP( re,gamma,pe )
    crho0  = csound*re

    pb  = pe-crho0*(sxn*(ub-ue)+syn*(vb-ve)+szn*(wb-we))

  ELSEIF (bcOpt == BCOPT_SUPERSONIC) THEN

    pb    = pInfl
  ENDIF

  rhob  = MixtPerf_D_PRT( pb,rgas,tb )
  rhoub = rhob*ub
  rhovb = rhob*vb
  rhowb = rhob*wb
  rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pb,ub,vb,wb )

END SUBROUTINE RFLO_BcondInflowVTPerf

!******************************************************************************
!
! Purpose: set inflow boundary condition for one cell based on presribed 
!          velocity and pressure.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              method of Rudy and Strikwerda: Boundary Conditions for Subsonic
!              Compressible Navier-Stokes Calculations.
!              Computers and Fluids 9, 327-338 (1981). The supersonic boundary 
!              condition consists of simple extrapolation.
!
! Input: bcOpt      = inflow bc-option: steady or unsteady
!        uInfl      = prescribed u-velocity
!        vInfl      = prescribed v-velocity
!        wInfl      = prescribed w-velocity
!        pInfl      = prescribed pressure 
!        tInfl      = prescribed temperature (only used if supersonic)
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!
! Notes: this condition is valid only for a thermally and calorically
!        perfect gas (supersonic in/outflow valid for all gases).
!
!******************************************************************************

SUBROUTINE RFLO_BcondInflowVPPerf( model,bcOpt,uInfl,vInfl,wInfl,tInfl,pInfl, &
                                   sxn,syn,szn,cpgas,mol, &
                                   rho,rhou,rhov,rhow,rhoe,press, &
                                   rhob,rhoub,rhovb,rhowb,rhoeb )

  USE ModInterfaces, ONLY: MixtPerf_C_DGP, MixtPerf_C_GRT, MixtPerf_D_PRT, &
                           MixtPerf_G_CpR, MixtPerf_R_M, MixtPerf_Eo_DGPUVW

  IMPLICIT NONE

! ... parameters
  INTEGER :: model, bcOpt
  REAL(RFREAL) :: uInfl, vInfl, wInfl, tInfl, pInfl
  REAL(RFREAL) :: sxn, syn, szn, cpgas, mol
  REAL(RFREAL) :: rho, rhou, rhov, rhow, rhoe, press
  REAL(RFREAL) :: rhob, rhoub, rhovb, rhowb, rhoeb

! ... local variables
  INTEGER :: ierror
  REAL(RFREAL) :: rgas, gamma, mach
  REAL(RFREAL) :: re, ue, ve, we, pe, te, qn, crho0, csound
  REAL(RFREAL) :: ub, vb, wb, tb, pb, ta

!******************************************************************************
! gas properties

  rgas  = MixtPerf_R_M( mol )
  gamma = MixtPerf_G_CpR( cpgas,rgas )

! flow values at a reference location (= interior cell)

  re  = rho
  ue  = rhou/rho
  ve  = rhov/rho
  we  = rhow/rho
  pe  = press
  te  = pe/(re*rgas)

! flow values at inlet plane

  pb  = pInfl
  ub  = uInfl
  vb  = vInfl
  wb  = wInfl
  qn  = sxn*ub + syn*vb + szn*wb
  ta  = pb/(re*rgas)
  mach= qn/MixtPerf_C_GRT( gamma,rgas,ta )

  IF (qn > 0 .AND. model==BCOPT_STEADY) THEN
    WRITE(STDERR,*)'BC_INFLOW_VELPRESS: outflow detected at inflow plane'
#ifdef MPI
    CALL MPI_Abort( ierror )
#endif
    STOP
  ENDIF

  IF (bcOpt == BCOPT_MIXED) THEN

! - subsonic inflow (qn<0)

    IF (mach < 1._RFREAL) THEN
      csound = MixtPerf_C_DGP( re,gamma,pe )
      crho0  = csound*re

      tb  = te-csound/rgas*(sxn*(ub-ue)+syn*(vb-ve)+szn*(wb-we))

! - supersonic inflow (qn<0)

    ELSE

      tb    = tInfl
    ENDIF

  ELSEIF (bcOpt == BCOPT_SUBSONIC) THEN
    csound = MixtPerf_C_DGP( re,gamma,pe )
    crho0  = csound*re

    tb  = te-csound/rgas*(sxn*(ub-ue)+syn*(vb-ve)+szn*(wb-we))

  ELSEIF (bcOpt == BCOPT_SUPERSONIC) THEN

    tb    = tInfl
  ENDIF

  rhob  = MixtPerf_D_PRT( pb,rgas,tb )
  rhoub = rhob*ub
  rhovb = rhob*vb
  rhowb = rhob*wb
  rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pb,ub,vb,wb )

END SUBROUTINE RFLO_BcondInflowVPPerf

!******************************************************************************
!
! Purpose: obtain fluctuations to be recycled at inflow boundary.
!
! Description: take fluctuations = instantaneous - time-averaged-quantities,
!              then add the fluctuations to corresp. imposed inlet quantities.
!
! Input: bcOptModel = inflow-bc model: steady (laminar) or unsteady (turbulent)
!        uvel       = instantaneous u-velocity
!        vvel       = instantaneous v-velocity
!        wvel       = instantaneous w-velocity
!        temp       = instantaneous temperature (only used if supersonic)
!        pres       = instantaneous pressure 
!        muvel      = time-averaged u-velocity
!        mvvel      = time-averaged v-velocity
!        mwvel      = time-averaged w-velocity
!        mtemp      = time-averaged temperature (only used if supersonic)
!        mpres      = time-averaged pressure 
!
! Output: ufluc = u fluctuation 
!         vfluc = v fluctuation
!         wfluc = w fluctuation
!         Tfluc = temperature fluctuation
!         pfluc = pressure fluctuation
!******************************************************************************
SUBROUTINE RFLO_BcondInflowFluc( bcOptModel, ufluc,vfluc,wfluc,Tfluc,pfluc, &
                                 uvel,vvel,wvel,temp,pres, &
                                 muvel,mvvel,mwvel,mtemp,mpres )
  IMPLICIT NONE

! ... parameter
  INTEGER :: bcOptModel
  REAL(RFREAL) :: ufluc, vfluc, wfluc, Tfluc, pfluc
  REAL(RFREAL) :: uvel, vvel, wvel, temp, pres
  REAL(RFREAL) :: muvel, mvvel, mwvel, mtemp, mpres

! ... local variables

!******************************************************************************

! obtain fluctuations to be recycled

  IF (bcOptModel==BCOPT_UNSTEADY) THEN
    ufluc = uvel-muvel
    vfluc = vvel-mvvel
    wfluc = wvel-mwvel
    Tfluc = temp-mtemp
    pfluc = pres-mpres
  ELSE
    ufluc = 0._RFREAL
    vfluc = 0._RFREAL
    wfluc = 0._RFREAL
    Tfluc = 0._RFREAL
    pfluc = 0._RFREAL
  ENDIF

END SUBROUTINE RFLO_BcondInflowFluc

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLO_ModBoundaryConditions

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModBoundaryConditions.F90,v $
! Revision 1.30  2010/02/18 21:47:39  juzhang
! Heat transfer bc for non-propellant surface documented in Rocburn_PY_HT.pdf in Rocburn_PY directory is implemented within Rocburn_PY. Major changes were made to Rocburn, Rocman3, RocfluidMP/genx, RocfluidMP/modflo directories.
!
! Revision 1.29  2008/12/06 08:44:14  mtcampbe
! Updated license.
!
! Revision 1.28  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.27  2006/08/19 15:38:50  mparmar
! Renamed patch variables
!
! Revision 1.26  2006/05/05 02:19:32  wasistho
! commented redundant kernels in RFLO_BondInflow
!
! Revision 1.25  2006/03/24 05:58:41  wasistho
! activated all components inlet fluctuations
!
! Revision 1.24  2006/03/09 21:57:22  wasistho
! lower backflow coefficients in NR outflow bc
!
! Revision 1.23  2006/02/18 08:28:10  wasistho
! added bcOpt sub/supersonic/mixed in veltemp and velpress
!
! Revision 1.22  2006/01/30 20:16:20  wasistho
! added safety for using injectionAPN
!
! Revision 1.21  2006/01/25 01:36:08  wasistho
! compute rhofVf in injectionAPN and mod injectionInit
!
! Revision 1.20  2006/01/23 22:46:46  wasistho
! added bc_internal condition for injectionAPN
!
! Revision 1.19  2006/01/21 11:38:10  wasistho
! added injectionapn
!
! Revision 1.18  2005/12/07 02:25:15  wasistho
! added ifdef STATS in inflow veltemp and velpress
!
! Revision 1.17  2005/12/06 22:08:22  wasistho
! ijkF to ijkD in the last fix
!
! Revision 1.16  2005/12/06 21:54:03  wasistho
! control fluctuations for iflow veltemp and velpress
!
! Revision 1.15  2005/12/01 09:33:43  wasistho
! added condition infloNijk < NIJK_INFLOW_INIT in bcondInflow
!
! Revision 1.14  2005/11/30 23:36:06  wasistho
! increased tollerance in eps
!
! Revision 1.13  2005/11/30 20:07:14  wasistho
! compute inflow fluctuations only when stats time>0
!
! Revision 1.12  2005/11/29 22:47:52  wasistho
! bug fixed divided tav by integrTime in BcondInflow
!
! Revision 1.11  2005/11/18 07:20:11  wasistho
! multiply veltemp/pres-inflow flucts by amplitude
!
! Revision 1.10  2005/11/07 21:23:27  wasistho
! increased backflow prevention parameter
!
! Revision 1.9  2005/10/31 22:34:04  wasistho
! added bcopt_steady in detecting inlet backflow
!
! Revision 1.8  2005/10/31 21:09:35  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.7  2005/10/26 22:23:31  wasistho
! lower backflow prevention parameter
!
! Revision 1.6  2005/09/28 18:06:18  wasistho
! modified for moving slipwall
!
! Revision 1.5  2005/09/20 23:13:29  wasistho
! added fluctuation option for inflow velTemp and velPress
!
! Revision 1.4  2005/08/31 05:36:34  wasistho
! undefined nrcoef for supersonic
!
! Revision 1.3  2005/05/13 08:12:38  wasistho
! implemented new bc_outflow model based on non-reflecting condition
!
! Revision 1.2  2005/04/28 05:43:44  wasistho
! added velocity based inflow BC
!
! Revision 1.1  2004/12/28 22:51:20  wasistho
! moved RFLO_Bcond* and RFLO_BoundaryCond* routines into RFLO_ModBoundaryConditions
!
!
! ******************************************************************************



















