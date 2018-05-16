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
! Purpose: Compute viscous fluxes based on viscosity principle: mu*S_ij through 
!          current patch with mu defined at cell faces
!
! Description: this routine works in the same way as TURB_VisFluxEddy
!              but applied on region patches
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: mixt%diss = total viscous fluxes added to dissipation
!
! Notes: This flux routine is used for the fixed and dynamic Smagorinsky model.
!        This subroutine is similar to subroutine RFLO_ViscousFluxPatch for the 
!        computation of viscous-fluxes at region patches based on mu defined 
!        at cell centers 
!
!******************************************************************************
!
! $Id: TURB_VisFluxEddyPatch.F90,v 1.14 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_VisFluxEddyPatch( region,patch )

!$startcopy TURB_vFluxHybridPatch
  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModTurbulence, ONLY : t_turb
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region
  TYPE(t_patch)          :: patch

! ... loop variables
  INTEGER :: i, j, k, iC

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: ijk, bcType, indCp, ijkCB0, ijkCD, ijkNB, ijkNBG
  INTEGER :: n1, n2, nOff, j2d, aeroCoeff

  REAL(RFREAL)          :: beta,rPrt,cpPrt,mueL,tCoL,tCoT,engModel
  REAL(RFREAL)          :: tx,ty,tz,tgradf, bValFactor
  REAL(RFREAL)          :: velf(3),fd(4),sf(3),sij(3,3)
  REAL(RFREAL)          :: velCB0(3),tfd(4),uciTauij,taukk, ac0,ac1
  REAL(RFREAL)          :: rRef, vRef, rCfRef, rChRef
  REAL(RFREAL), POINTER :: diss(:,:),tv(:,:),gv(:,:),mueT(:,:),sRate(:,:)
  REAL(RFREAL), POINTER :: trace(:),vol(:)

#ifdef RFLO
  INTEGER :: inode, jnode, knode, idir, jdir, kdir
  INTEGER :: ilev, lbound, iCOff, ijCOff, iNOff, ijNOff, acId0, acId1
  REAL(RFREAL)          :: sgn
  REAL(RFREAL), POINTER :: avgCo(:,:), sFace(:,:), dv(:,:), grad(:,:)
#endif
#ifdef RFLU
  INTEGER :: ifgBeg
  REAL(RFREAL)          :: rDens0
  REAL(RFREAL), POINTER :: fn(:,:), cv(:,:), grad(:,:,:) 
#endif
!$endcopy

  REAL(RFREAL)          :: muf, modStrain

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_VisFluxEddyPatch.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_VisFluxEddyPatch',&
  'TURB_VisFluxEddyPatch.F90' )

!$startcopy TURB_vFluxHybridPatch
! get dimensions and pointers -------------------------------------------------

  bcType = patch%bcType

#ifdef RFLO
  ilev    =  region%currLevel
  lbound  =  patch%lbound
  vol     => region%levels(ilev)%grid%vol
  dv      => region%levels(ilev)%mixt%dv
  tv      => region%levels(ilev)%mixt%tv
  gv      => region%levels(ilev)%mixt%gv
  diss    => region%levels(ilev)%mixt%diss  
  mueT    => region%levels(ilev)%turb%mueT
  trace   => region%levels(ilev)%turb%trace
  turb    => region%levels(ilev)%turb
#endif  
#ifdef RFLU
  vol     => region%grid%vol
  cv      => region%mixt%cv
  tv      => region%mixt%tv
  gv      => region%mixt%gv
  diss    => region%mixt%diss  
  mueT    => region%turb%bMueT
  trace   => region%turb%trace
  turb    => region%turb
#endif  

! get coefficients -----------------------------------------------------------

  IF (region%turbInput%engModel==0) THEN
    engModel = 0._RFREAL
  ELSE
    engModel = 1._RFREAL
  ENDIF
  beta  = region%mixtInput%betrk(region%irkStep)

  IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
    bValFactor = 0._RFREAL
  ELSE
    bValFactor = 1._RFREAL
  ENDIF
  
#ifdef RFLO
  rPrt      = 1._RFREAL/region%levels(iLev)%mixt%prTurb
  indCp     = region%levels(iLev)%mixt%indCp
  aeroCoeff = global%aeroCoeffs
  nOff      = ABS(patch%l1end-patch%l1beg) + 1
#endif
#ifdef RFLU
  rPrt      = 1._RFREAL/region%mixtInput%prTurb
  indCp     = region%mixtInput%indCp
#endif

  rRef   = global%refDensity
  vRef   = global%refVelocity
  
  rCfRef = 2.0_RFREAL/(rRef*vRef*vRef)
  rChRef = 2.0_RFREAL/(rRef*vRef*vRef*vRef)

#ifdef RFLO

! take the right face vector and make it point outwards

  CALL RFLO_GetPatchIndices( region,patch,ilev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  acId0 = 2
  acId1 = 1
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
    acId0 = 1
    acId1 = 2
  ENDIF
!$endcopy

! get the appropriate face vector

  IF (lbound==1 .OR. lbound==2) THEN
    ijk   =  DIRI
    avgCo => region%levels(iLev)%grid%c2fCoI
    sFace => region%levels(ilev)%grid%si
    grad  => region%levels(ilev)%mixt%gradi
    sRate => region%levels(ilev)%turb%mISij
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    ijk   =  DIRJ
    avgCo => region%levels(iLev)%grid%c2fCoJ
    sFace => region%levels(ilev)%grid%sj
    grad  => region%levels(ilev)%mixt%gradj
    sRate => region%levels(ilev)%turb%mJSij
  ELSE
    ijk   =  DIRK
    avgCo => region%levels(iLev)%grid%c2fCoK
    sFace => region%levels(ilev)%grid%sk
    grad  => region%levels(ilev)%mixt%gradk
    sRate => region%levels(ilev)%turb%mKSij
  ENDIF

!$startcopy TURB_vFluxHybridPatch

! non-conforming region interface

  IF (bcType>=BC_regionINT .AND. bcType<=BC_regionINT+BC_RANGE) THEN

  ELSE IF (bcType>=BC_regNONCONF .AND. bcType<=BC_regNONCONF+BC_RANGE) THEN

! everything else

  ELSE

! flux in the direction normal to the patch

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! bnd cells
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)  ! bnd nodes
          ijkNBG = ijkNB

          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
          ac0 = avgCo(acId0,ijkNB)
          ac1 = avgCo(acId1,ijkNB)
#endif
#ifdef RFLU
! specific Rocflu, check the state of cv first
  IF (region%mixt%cvState /= CV_MIXT_STATE_CONS) &
                            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  ijk    = DIRI
  ibeg   = 1
  iend   = patch%nBFaces
! TEMPORARY : removing usage of bf2bg from everywhere
!  ifgBeg = patch%bf2bg(BF2BG_BEG)
  ac0    = 0.5_RFREAL
  ac1    = 0.5_RFREAL

  fn    => patch%fn
! TEMPORARY 
  grad  => patch%mixt%gradFace
  sRate => region%turb%bmISij  

  DO iC=ibeg,iend
          ijkCB0 = patch%bf2c(iC)
          ijkCD  = ijkCB0
          ijkNB  = iC
          ijkNBG = iC + ifgBeg-1

          sf(1)  = fn(XCOORD,ijkNB)*fn(XYZMAG,ijkNB)
          sf(2)  = fn(YCOORD,ijkNB)*fn(XYZMAG,ijkNB)
          sf(3)  = fn(ZCOORD,ijkNB)*fn(XYZMAG,ijkNB)           
          rDens0 = 1._RFREAL/cv(CV_MIXT_DENS,ijkCB0)
#endif

!pert TURB
          cpPrt = (ac0*gv(GV_MIXT_CP,ijkCB0*indCp) + &
                   ac1*gv(GV_MIXT_CP,ijkCD*indCp))*rPrt
          mueL  = ac0*tv(TV_MIXT_MUEL,ijkCB0)+ac1*tv(TV_MIXT_MUEL,ijkCD)
          tCoL  = ac0*tv(TV_MIXT_TCOL,ijkCB0)+ac1*tv(TV_MIXT_TCOL,ijkCD)
          tCoT  = cpPrt*mueT(ijk,ijkNBG)
!endpert

#ifdef RFLO
          velCB0(1)= dv(DV_MIXT_UVEL,ijkCB0)
          velCB0(2)= dv(DV_MIXT_VVEL,ijkCB0)
          velCB0(3)= dv(DV_MIXT_WVEL,ijkCB0)
          velf(1)= ac0*dv(DV_MIXT_UVEL,ijkCB0)+ac1*dv(DV_MIXT_UVEL,ijkCD)
          velf(2)= ac0*dv(DV_MIXT_VVEL,ijkCB0)+ac1*dv(DV_MIXT_VVEL,ijkCD)
          velf(3)= ac0*dv(DV_MIXT_WVEL,ijkCB0)+ac1*dv(DV_MIXT_WVEL,ijkCD)
          tx  = grad(GR_MIXT_TX,ijkNB)
          ty  = grad(GR_MIXT_TY,ijkNB)
          tz  = grad(GR_MIXT_TZ,ijkNB)
#endif
#ifdef RFLU
          velCB0(1)= cv(CV_MIXT_XMOM,ijkCB0)*rDens0
          velCB0(2)= cv(CV_MIXT_YMOM,ijkCB0)*rDens0
          velCB0(3)= cv(CV_MIXT_ZMOM,ijkCB0)*rDens0
          velf(1)  = velCB0(1)
          velf(2)  = velCB0(2)
          velf(3)  = velCB0(3)
          velf(1:3)= bValFactor*velf(1:3)
          tx  = grad(XCOORD,GRF_MIXT_TEMP,ijkNBG)
          ty  = grad(YCOORD,GRF_MIXT_TEMP,ijkNBG)
          tz  = grad(ZCOORD,GRF_MIXT_TEMP,ijkNBG)
#endif

          tgradf = tx*sf(1)+ty*sf(2)+tz*sf(3)

          sij(1,1) = sRate(E11,ijkNBG)
          sij(1,2) = sRate(E12,ijkNBG)
          sij(1,3) = sRate(E13,ijkNBG)

          sij(2,1) = sij(1,2)
          sij(2,2) = sRate(E22,ijkNBG)
          sij(2,3) = sRate(E23,ijkNBG)

          sij(3,1) = sij(1,3)
          sij(3,2) = sij(2,3)
          sij(3,3) = sRate(E33,ijkNBG)

          fd(1) = mueL*(sij(1,1)*sf(1)+sij(1,2)*sf(2)+sij(1,3)*sf(3))
          fd(2) = mueL*(sij(2,1)*sf(1)+sij(2,2)*sf(2)+sij(2,3)*sf(3))
          fd(3) = mueL*(sij(3,1)*sf(1)+sij(3,2)*sf(2)+sij(3,3)*sf(3))
          fd(4) = DOT_PRODUCT(fd(1:3),velf(1:3)) + tCoL*tgradf
!wlmCheckprobe--------------------------------------------------------------
!  IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
!    write(*,*) region%procId, i, j, k, mueL*sij(1,2), mueL*sij(3,2), &
!               mueL*sij(2,1), mueL*sij(2,3), tCoL*tgradf
!  ENDIF
!---------------------------------------------------------------------------
!$endcopy
! Yoshizawa k-model for eddy viscosity type LES
          modStrain = sqrt(sij(1,1)*sij(1,1) &
                          +sij(2,2)*sij(2,2) &
                          +sij(1,1)*sij(2,2) &
                          +sij(1,2)*sij(1,2) &
                          +sij(1,3)*sij(1,3) &
                          +sij(2,3)*sij(2,3))
          taukk = 2._RFREAL*mueT(ijk,ijkNBG)*modStrain
          trace(ijkCB0) = trace(ijkCB0) + taukk

          muf      = -mueT(ijk,ijkNBG)/mueL
          tfd(1)   = muf*fd(1)
          tfd(2)   = muf*fd(2)
          tfd(3)   = muf*fd(3)
          uciTauij = DOT_PRODUCT(tfd(1:3),velCB0(1:3))
          tfd(4)   = uciTauij - tCoT*tgradf

!$startcopy TURB_vFluxHybridPatch
          fd(1) = fd(1)-tfd(1)
          fd(2) = fd(2)-tfd(2)
          fd(3) = fd(3)-tfd(3)
          fd(4) = fd(4)-engModel*tfd(4)

#ifdef STATS
          turb%stwork(1,ijkCB0)=turb%stwork(1,ijkCB0)+uciTauij    ! alp_1
          turb%stwork(2,ijkCB0)=turb%stwork(2,ijkCB0)+tCoT*tgradf ! alp_2+alp_3
#endif
          diss(CV_MIXT_XMOM,ijkCB0) = diss(CV_MIXT_XMOM,ijkCB0)+fd(1)*beta
          diss(CV_MIXT_YMOM,ijkCB0) = diss(CV_MIXT_YMOM,ijkCB0)+fd(2)*beta
          diss(CV_MIXT_ZMOM,ijkCB0) = diss(CV_MIXT_ZMOM,ijkCB0)+fd(3)*beta
          diss(CV_MIXT_ENER,ijkCB0) = diss(CV_MIXT_ENER,ijkCB0)+fd(4)*beta
          global%esg1Sum  = global%esg1Sum + vol(ijkCB0)*uciTauij

! ------- Set friction and heat-transfer coefficients

#ifdef RFLO
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
          j2d  = aeroCoeff * IndIJ(n1,n2,nOff)
          patch%cf(XCOORD,j2d) = rCfRef*fd(1)
          patch%cf(YCOORD,j2d) = rCfRef*fd(2)
          patch%cf(ZCOORD,j2d) = rCfRef*fd(3)
          patch%ch(       j2d) = rChRef*fd(4)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
  ENDIF         ! bcType
#endif
#ifdef RFLU
          patch%cf(XCOORD,iC) = rCfRef*fd(1)
          patch%cf(YCOORD,iC) = rCfRef*fd(2)
          patch%cf(ZCOORD,iC) = rCfRef*fd(3)
          patch%ch(       iC) = rChRef*fd(4)               
  ENDDO         ! iC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )
!$endcopy
END SUBROUTINE TURB_VisFluxEddyPatch

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_VisFluxEddyPatch.F90,v $
! Revision 1.14  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2006/08/19 15:40:40  mparmar
! Removed bf2bg and moved region%mixt%bGradFace to patch%mixt%gradFace
!
! Revision 1.11  2006/03/13 03:43:23  wasistho
! defined aero coeffs
!
! Revision 1.10  2004/10/22 23:20:39  wasistho
! collect esg1 and esg2+3 into statistics
!
! Revision 1.9  2004/08/04 22:07:19  wasistho
! bugfixed: ac0 and ac1 are common to flo and flu
!
! Revision 1.8  2004/08/02 23:06:36  wasistho
! mv libfloflu/viscousFluxEddy(Patch) to rocflo/RFLO_viscousFlux(Patch)
!
! Revision 1.7  2004/08/02 21:06:17  wasistho
! moved averaging coeff local variables to within ifdef RFLO
!
! Revision 1.6  2004/08/02 19:34:18  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.5  2004/07/31 00:53:38  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.4  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.3  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/24 03:37:03  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.7  2004/02/26 21:22:23  wasistho
! changed turb%esg.. to global%esg..
!
! Revision 1.6  2004/01/21 03:42:18  wasistho
! get rid of grad(T)*n = 0 for adiabatic wall
!
! Revision 1.5  2003/08/29 01:43:17  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.4  2003/06/05 19:21:16  wasistho
! adiabatic switch
!
! Revision 1.3  2003/05/31 01:48:23  wasistho
! installed turb. wall layer model
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!******************************************************************************







