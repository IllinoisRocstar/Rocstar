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
! Purpose: Compute total viscous fluxes of which turbulent stress contribution
!          is based on stress model as opposed to eddy viscosity model.
!
! Description: The fluxes are computed in three consecutive sweeps (i, j and k
!              direction) at interior faces. The fluxes at the boundary patches
!              are then treated. The energy model contributions are also
!              treated on the fly.
!
! Input: region  = data of current region
!
! Output: mixt%diss = viscous fluxes added to the dissipation.
!
! Notes: This flux routine is used for the scale similarity and dyn.mixed model.
!
!******************************************************************************
!
! $Id: TURB_VFluxHybrid.F90,v 1.10 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_VFluxHybrid( region )

!$startcopy TURB_VisFluxEddy
  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModTurbulence, ONLY : t_turb
!$endcopy
  USE TURB_ModInterfaces, ONLY : TURB_lesEsgModel4, TURB_VFluxHybridPatch, &
                                 TURB_WlmTauWallMapping, TURB_WlmUpdate, &
                                 TURB_WlmFluxPatch, TURB_StatCCollector
!$startcopy TURB_VisFluxEddy
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... loop variables
  INTEGER :: i, j, k, iC, ipatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb

  INTEGER :: ijkC0,ijkC1,ijkN, indCp, ibn,ien
  REAL(RFREAL)          :: one6th,beta,engModel,rPrt,cpPrt,mueL,tCoL,tCoT
  REAL(RFREAL)          :: muf,tx,ty,tz,tgradf
  REAL(RFREAL)          :: velf(3),fd(4),sFace(3),sij(3,3)
  REAL(RFREAL)          :: velC0(3),velC1(3),tfd(4),fd4,uciTauij,taukk
  REAL(RFREAL), POINTER :: diss(:,:),tv(:,:),gv(:,:),mueT(:,:),sRate(:,:)
  REAL(RFREAL), POINTER :: trace(:),vol(:)
  LOGICAL               :: doWlm

#ifdef RFLO
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ilev,iCOff,ijCOff,iNOff,ijNOff
  REAL(RFREAL), POINTER :: avgCo(:,:), sf(:,:), dv(:,:), grad(:,:)
#endif
#ifdef RFLU
  INTEGER, POINTER      :: f2c(:,:)
  REAL(RFREAL)          :: rDens0, rDens1
  REAL(RFREAL), POINTER :: fn(:,:), cv(:,:), grad(:,:,:)
#endif
!$endcopy
!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_VFluxHybrid.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_VFluxHybrid',&
  'TURB_VFluxHybrid.F90' )

!$startcopy TURB_VisFluxEddy
! get dimensions and pointers ------------------------------------------------

#ifdef RFLO
  ilev    =  region%currLevel
  dv      => region%levels(ilev)%mixt%dv
  tv      => region%levels(ilev)%mixt%tv
  gv      => region%levels(ilev)%mixt%gv
  diss    => region%levels(ilev)%mixt%diss
  mueT    => region%levels(ilev)%turb%mueT
  trace   => region%levels(ilev)%turb%trace
  vol     => region%levels(ilev)%grid%vol
  turb    => region%levels(ilev)%turb

  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )
#endif
#ifdef RFLU
  cv      => region%mixt%cv
  tv      => region%mixt%tv
  gv      => region%mixt%gv
  diss    => region%mixt%diss  
  mueT    => region%turb%mueT
  trace   => region%turb%trace
  vol     => region%grid%vol
  turb    => region%turb
#endif

! get coefficients and parameters

  one6th = 1._RFREAL/6._RFREAL

  IF (region%turbInput%engModel==OFF) THEN
    engModel = 0._RFREAL
  ELSE
    engModel = 1._RFREAL
  ENDIF
  beta  = region%mixtInput%betrk(region%irkStep)

#ifdef RFLO
  rPrt  = 1._RFREAL/region%levels(iLev)%mixt%prTurb
  indCp = region%levels(iLev)%mixt%indCp
#endif
#ifdef RFLU
  rPrt  = 1._RFREAL/region%mixtInput%prTurb
  indCp = region%mixtInput%indCp
#endif

! compute interior fluxes through I, J and K faces

#ifdef STATS
! if desired, collect quantities of interest; for that, allocate work space
  ibn = LBOUND( trace,1 )
  ien = UBOUND( trace,1 )
  ALLOCATE( turb%stwork(2,ibn:ien) )
  turb%stwork = 0._RFREAL
#endif

  CALL ComputeFluxTot( DIRI )
#ifdef RFLO
  CALL ComputeFluxTot( DIRJ )
  CALL ComputeFluxTot( DIRK )
#endif

! get fluxes through boundaries
#ifdef RFLO  
  DO ipatch=1,region%nPatches
    patch => region%levels(ilev)%patches(ipatch)
#endif
#ifdef RFLU 
  DO ipatch=1,region%grid%nPatches
    patch => region%patches(ipatch)
#endif

    doWlm = .false.
#ifdef RFLO
    IF (patch%bcType>=BC_NOSLIPWALL .AND. &
        patch%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN   ! my boundary type
      IF (patch%valBola%switches(WLM_INPUT_MODEL) /= WLM_MODEL_NOMODEL) THEN
        doWlm = .true.
      ENDIF
    ENDIF
#endif
!$endcopy

    IF (doWlm) THEN
      CALL TURB_WlmUpdate( region,patch )
      CALL TURB_WlmTauWallMapping( region,patch )  
      CALL TURB_WlmFluxPatch( region,patch )
    ELSE
      CALL TURB_VFluxHybridPatch( region,patch )
    ENDIF
  ENDDO

!$startcopy TURB_VisFluxEddy
#ifdef STATS
! collect statistics and deallocate work space used

  IF ((region%turbInput%nSt > 0) .AND. (region%turbInput%engModel/=OFF)) THEN
    turb%stwork(:,:) = one6th*turb%stwork(:,:)
    CALL TURB_StatCCollector( region,1,2,turb%stwork )
  ENDIF
  DEALLOCATE( turb%stwork )
#endif

! if active, get contribution of energy subgrid model alpha_4
  IF (region%turbInput%engModel/=OFF) THEN
    CALL TURB_lesEsgModel4( region )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! =============================================================================
!   Flux computation subroutine
! =============================================================================

CONTAINS

  SUBROUTINE ComputeFluxTot( ijk )

! ... parameters
    INTEGER   :: ijk

! ... local variables
    INTEGER   :: ibeg,iend,jbeg,jend,kbeg,kend, iadd,jadd,kadd
    REAL(RFREAL) :: ac0, ac1
!$endcopy
!pert TURB 
    REAL(RFREAL)          :: tij(3,3),sijSij,sijTij
    REAL(RFREAL), POINTER :: tauij(:,:)
!endpert

! - Compute turbulent part tauij from scale-similarity; for Bardina model,
!   face values of cv have yet to be obtained.

#ifdef RFLO
    IF (ijk==DIRI) THEN
      tauij => region%levels(ilev)%turb%fISij
    ELSEIF (ijk==DIRJ) THEN
      tauij => region%levels(ilev)%turb%fJSij
    ELSEIF (ijk==DIRK) THEN
      tauij => region%levels(ilev)%turb%fKSij
    ENDIF
#endif
#ifdef RFLU
    tauij => region%turb%fISij
#endif

!$startcopy TURB_VisFluxEddy
! - Set limits and pointers ---------------------------------------------------

#ifdef RFLO
    IF (ijk==DIRI) THEN
      ibeg = ipcbeg+1
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = -1
      jadd = 0
      kadd = 0
      sRate => region%levels(ilev)%turb%mISij
      grad  => region%levels(ilev)%mixt%gradi
      sf    => region%levels(ilev)%grid%si
      avgCo => region%levels(iLev)%grid%c2fCoI
    ELSEIF (ijk==DIRJ) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg+1
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = 0
      jadd = -1
      kadd = 0
      sRate => region%levels(ilev)%turb%mJSij
      grad  => region%levels(ilev)%mixt%gradj
      sf    => region%levels(ilev)%grid%sj
      avgCo => region%levels(iLev)%grid%c2fCoJ
    ELSEIF (ijk==DIRK) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg+1
      kend = kpcend
      iadd = 0
      jadd = 0
      kadd = -1
      sRate => region%levels(ilev)%turb%mKSij
      grad  => region%levels(ilev)%mixt%gradk
      sf    => region%levels(ilev)%grid%sk
      avgCo => region%levels(iLev)%grid%c2fCoK
    ENDIF
#endif
#ifdef RFLU
    ibeg =  1
    iend =  region%grid%nFaces
    f2c   => region%grid%f2c
    sRate => region%turb%mISij  
    grad  => region%mixt%gradFace
    fn    => region%grid%fn
#endif

#ifdef RFLO
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          ijkC1 = IndIJK(i+iadd,j+jadd,k+kadd,iCOff,ijCOff)
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          sFace(1)= sf(XCOORD,ijkN)
          sFace(2)= sf(YCOORD,ijkN)
          sFace(3)= sf(ZCOORD,ijkN)
          ac0 = avgCo(2,ijkN)
          ac1 = avgCo(1,ijkN)
#endif
#ifdef RFLU
! - specific Rocflu, check the state of cv first
    IF (region%mixt%cvState /= CV_MIXT_STATE_CONS) &
                            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    ac0 = 0.5_RFREAL
    ac1 = 0.5_RFREAL
    DO iC=ibeg,iend
          ijkC0 = f2c(1,iC)
          ijkC1 = f2c(2,iC)
          ijkN  = iC
          sFace(1)= fn(XCOORD,ijkN)*fn(XYZMAG,ijkN)
          sFace(2)= fn(YCOORD,ijkN)*fn(XYZMAG,ijkN)
          sFace(3)= fn(ZCOORD,ijkN)*fn(XYZMAG,ijkN)
          rDens0  = 1._RFREAL/cv(CV_MIXT_DENS,ijkC0)
          rDens1  = 1._RFREAL/cv(CV_MIXT_DENS,ijkC1)
#endif
          cpPrt = (ac0*gv(GV_MIXT_CP,ijkC0*indCp) + &
                   ac1*gv(GV_MIXT_CP,ijkC1*indCp))*rPrt
          mueL  = ac0*tv(TV_MIXT_MUEL,ijkC0)+ac1*tv(TV_MIXT_MUEL,ijkC1)
          tCoL  = ac0*tv(TV_MIXT_TCOL,ijkC0)+ac1*tv(TV_MIXT_TCOL,ijkC1)
          tCoT  = cpPrt*mueT(ijk,ijkN)
#ifdef RFLO
          velC0(1)= dv(DV_MIXT_UVEL,ijkC0)
          velC0(2)= dv(DV_MIXT_VVEL,ijkC0)
          velC0(3)= dv(DV_MIXT_WVEL,ijkC0)
          velC1(1)= dv(DV_MIXT_UVEL,ijkC1)
          velC1(2)= dv(DV_MIXT_VVEL,ijkC1)
          velC1(3)= dv(DV_MIXT_WVEL,ijkC1)
          velf(1) = ac0*dv(DV_MIXT_UVEL,ijkC0)+ac1*dv(DV_MIXT_UVEL,ijkC1)
          velf(2) = ac0*dv(DV_MIXT_VVEL,ijkC0)+ac1*dv(DV_MIXT_VVEL,ijkC1)
          velf(3) = ac0*dv(DV_MIXT_WVEL,ijkC0)+ac1*dv(DV_MIXT_WVEL,ijkC1)
          tx  = grad(GR_MIXT_TX,ijkN)
          ty  = grad(GR_MIXT_TY,ijkN)
          tz  = grad(GR_MIXT_TZ,ijkN)
#endif
#ifdef RFLU
          velC0(1)= cv(CV_MIXT_XMOM,ijkC0)*rDens0
          velC0(2)= cv(CV_MIXT_YMOM,ijkC0)*rDens0
          velC0(3)= cv(CV_MIXT_ZMOM,ijkC0)*rDens0
          velC1(1)= cv(CV_MIXT_XMOM,ijkC1)*rDens1
          velC1(2)= cv(CV_MIXT_YMOM,ijkC1)*rDens1
          velC1(3)= cv(CV_MIXT_ZMOM,ijkC1)*rDens1
          velf(1) = ac0*velC0(1)+ac1*velC1(1)
          velf(2) = ac0*velC0(2)+ac1*velC1(2)
          velf(3) = ac0*velC0(3)+ac1*velC1(3)
          tx  = grad(XCOORD,GRF_MIXT_TEMP,ijkN)
          ty  = grad(YCOORD,GRF_MIXT_TEMP,ijkN)
          tz  = grad(ZCOORD,GRF_MIXT_TEMP,ijkN)
#endif

          tgradf= (tx*sFace(1)+ty*sFace(2)+tz*sFace(3))

          sij(1,1) = sRate(E11,ijkN)
          sij(1,2) = sRate(E12,ijkN)
          sij(1,3) = sRate(E13,ijkN)

          sij(2,1) = sij(1,2)
          sij(2,2) = sRate(E22,ijkN)
          sij(2,3) = sRate(E23,ijkN)

          sij(3,1) = sij(1,3)
          sij(3,2) = sij(2,3)
          sij(3,3) = sRate(E33,ijkN)

          fd(1) = mueL*(sij(1,1)*sFace(1)+sij(1,2)*sFace(2)+sij(1,3)*sFace(3)) 
          fd(2) = mueL*(sij(2,1)*sFace(1)+sij(2,2)*sFace(2)+sij(2,3)*sFace(3)) 
          fd(3) = mueL*(sij(3,1)*sFace(1)+sij(3,2)*sFace(2)+sij(3,3)*sFace(3)) 
          fd(4) = DOT_PRODUCT(fd(1:3),velf(1:3)) + tCoL*tgradf
!$endcopy
          tij(1,1) = tauij(E11,ijkN) - mueT(ijk,ijkN)*sij(1,1)
          tij(1,2) = tauij(E12,ijkN) - mueT(ijk,ijkN)*sij(1,2)
          tij(1,3) = tauij(E13,ijkN) - mueT(ijk,ijkN)*sij(1,3)

          tij(2,1) = tauij(E12,ijkN) - mueT(ijk,ijkN)*sij(2,1)
          tij(2,2) = tauij(E22,ijkN) - mueT(ijk,ijkN)*sij(2,2)
          tij(2,3) = tauij(E23,ijkN) - mueT(ijk,ijkN)*sij(2,3)

          tij(3,1) = tauij(E13,ijkN) - mueT(ijk,ijkN)*sij(3,1)
          tij(3,2) = tauij(E23,ijkN) - mueT(ijk,ijkN)*sij(3,2)
          tij(3,3) = tauij(E33,ijkN) - mueT(ijk,ijkN)*sij(3,3)

          taukk = tij(1,1)+tij(2,2)+tij(3,3)
          taukk = MAX( 0._RFREAL,taukk )
          trace(ijkC0) = trace(ijkC0) + taukk
          trace(ijkC1) = trace(ijkC1) + taukk

          tfd(1) = tij(1,1)*sFace(1)+tij(1,2)*sFace(2)+tij(1,3)*sFace(3)
          tfd(2) = tij(2,1)*sFace(1)+tij(2,2)*sFace(2)+tij(2,3)*sFace(3)
          tfd(3) = tij(3,1)*sFace(1)+tij(3,2)*sFace(2)+tij(3,3)*sFace(3)

          fd(1) = fd(1)-tfd(1)
          fd(2) = fd(2)-tfd(2)
          fd(3) = fd(3)-tfd(3)

!$startcopy TURB_VisFluxEddy
          fd4      = fd(4)
          uciTauij = DOT_PRODUCT(tfd(1:3),velC0(1:3))
          tfd(4)   = uciTauij - tCoT*tgradf
          fd(4)    = fd4 - engModel*tfd(4)

#ifdef STATS
          turb%stwork(1,ijkC0)=turb%stwork(1,ijkC0)+uciTauij    ! alp_1
          turb%stwork(2,ijkC0)=turb%stwork(2,ijkC0)+tCoT*tgradf ! alp_2+alp_3
#endif
          diss(CV_MIXT_XMOM,ijkC0) = diss(CV_MIXT_XMOM,ijkC0) + fd(1)*beta
          diss(CV_MIXT_YMOM,ijkC0) = diss(CV_MIXT_YMOM,ijkC0) + fd(2)*beta
          diss(CV_MIXT_ZMOM,ijkC0) = diss(CV_MIXT_ZMOM,ijkC0) + fd(3)*beta
          diss(CV_MIXT_ENER,ijkC0) = diss(CV_MIXT_ENER,ijkC0) + fd(4)*beta
          global%esg1Sum  = global%esg1Sum + vol(ijkC0)*uciTauij

          uciTauij = DOT_PRODUCT(tfd(1:3),velC1(1:3))
          tfd(4)   = uciTauij - tCoT*tgradf
          fd(4)    = fd4 - engModel*tfd(4)

#ifdef STATS
          turb%stwork(1,ijkC1)=turb%stwork(1,ijkC1)+uciTauij    ! alp_1
          turb%stwork(2,ijkC1)=turb%stwork(2,ijkC1)+tCoT*tgradf ! alp_2+alp_3
#endif
          diss(CV_MIXT_XMOM,ijkC1) = diss(CV_MIXT_XMOM,ijkC1) - fd(1)*beta
          diss(CV_MIXT_YMOM,ijkC1) = diss(CV_MIXT_YMOM,ijkC1) - fd(2)*beta
          diss(CV_MIXT_ZMOM,ijkC1) = diss(CV_MIXT_ZMOM,ijkC1) - fd(3)*beta
          diss(CV_MIXT_ENER,ijkC1) = diss(CV_MIXT_ENER,ijkC1) - fd(4)*beta
          global%esg1Sum  = global%esg1Sum - vol(ijkC1)*uciTauij

#ifdef RFLO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
#endif
#ifdef RFLU
    ENDDO       ! iC
#endif

!$endcopy

#ifdef RFLO
! compute equivalent eddy viscosity from scale similarity

    IF (ijk==DIRI) THEN
      ibeg = ipcbeg-1
      iend = ipcend+2
      jbeg = jpcbeg-1
      jend = jpcend+1
      kbeg = kpcbeg-1
      kend = kpcend+1
    ELSEIF (ijk==DIRJ) THEN
      ibeg = ipcbeg-1
      iend = ipcend+1
      jbeg = jpcbeg-1
      jend = jpcend+2
      kbeg = kpcbeg-1
      kend = kpcend+1
    ELSEIF (ijk==DIRK) THEN
      ibeg = ipcbeg-1
      iend = ipcend+1
      jbeg = jpcbeg-1
      jend = jpcend+1
      kbeg = kpcbeg-1
      kend = kpcend+2
    ENDIF

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
#endif
#ifdef RFLU
    ibeg =  1
    iend =  region%grid%nFaces

    DO ijkN=ibeg,iend
#endif
          sij(1,1) = sRate(E11,ijkN)
          sij(1,2) = sRate(E12,ijkN)
          sij(1,3) = sRate(E13,ijkN)

          sij(2,1) = sij(1,2)
          sij(2,2) = sRate(E22,ijkN)
          sij(2,3) = sRate(E23,ijkN)

          sij(3,1) = sij(1,3)
          sij(3,2) = sij(2,3)
          sij(3,3) = sRate(E33,ijkN)

! ------- S_{ij}*S_{ij}:
          sijSij = sij(1,1)*sij(1,1)  &
                  +sij(2,2)*sij(2,2)  &
                  +sij(3,3)*sij(3,3)  &
            +2.0d0*sij(1,2)*sij(1,2)  &
            +2.0d0*sij(1,3)*sij(1,3)  &
            +2.0d0*sij(2,3)*sij(2,3)

! ------- S_{ij}*Tau_{ij}:
          sijTij = sij(1,1)*tauij(E11,ijkN)  &
                  +sij(2,2)*tauij(E22,ijkN)  &
                  +sij(3,3)*tauij(E33,ijkN)  &
            +2.0d0*sij(1,2)*tauij(E12,ijkN)  &
            +2.0d0*sij(1,3)*tauij(E13,ijkN)  &
            +2.0d0*sij(2,3)*tauij(E23,ijkN)
        
          muf = MAX( 0._RFREAL, -(sijTij/MAX( sijSij,REAL_SMALL )) )        
          mueT(ijk,ijkN) = mueT(ijk,ijkN) + muf
#ifdef RFLO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
#endif
#ifdef RFLU
    ENDDO       ! ijkN
#endif

  END SUBROUTINE ComputeFluxTot

END SUBROUTINE TURB_VFluxHybrid

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_VFluxHybrid.F90,v $
! Revision 1.10  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2005/12/30 23:20:56  wasistho
! exclude rocflu from WLM treatment
!
! Revision 1.7  2004/10/22 23:20:33  wasistho
! collect esg1 and esg2+3 into statistics
!
! Revision 1.6  2004/08/02 19:33:58  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.5  2004/07/31 00:53:47  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.4  2004/05/29 03:38:35  wasistho
! extended computations of muet to first dummy layers
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/24 03:37:02  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.7  2004/02/26 21:22:31  wasistho
! changed turb%esg.. to global%esg..
!
! Revision 1.6  2003/08/29 01:43:00  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.5  2003/06/05 19:20:23  wasistho
! implemented heat transfer model
!
! Revision 1.4  2003/05/31 01:48:23  wasistho
! installed turb. wall layer model
!
! Revision 1.3  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.2  2003/02/04 04:00:11  wasistho
! An Esgs-bug which caused fd(4) unsymmetry fixed
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!******************************************************************************







