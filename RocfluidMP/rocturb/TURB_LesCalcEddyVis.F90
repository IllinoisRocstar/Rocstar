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
! Purpose: Compute turbulent viscosity mu_t.
!
! Description: The eddy viscosity is computed at faces using the model 
!              coefficient depending on the model selected by user.
!
! Input: region  = data of current region
!        ibn,ien = begin and end nodes index
!        ijk     = integer identifying wich face is being treated
!
! Output: mueT(ijk,:) = eddy viscosity at ijk face
!
! Notes: The face values contribution are collected here for model coeffs. 
!        at cell centers.
!
!******************************************************************************
!
! $Id: TURB_LesCalcEddyVis.F90,v 1.17 2009/08/12 04:15:59 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesCalcEddyVis( region,ibn,ien,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModTurbulence, ONLY      : t_turb
  USE ModGlobal, ONLY          : t_global
#ifdef RFLU
  USE ModBndPatch, ONLY        : t_patch
  USE TURB_ModInterfaces, ONLY : TURB_FluLesC2F
#endif
#ifdef RFLO
  USE ModInterfaces, ONLY      : RFLO_GetDimensPhys, &
                                 RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesGenC2F
                                 
#include "Indexing.h"
#endif
  USE TURB_ModInterfaces, ONLY : TURB_LesCoefDynSmag, TURB_LesCoefDynMixd
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: region
#endif
  INTEGER                 :: ibn,ien,ijk

! ... loop variables
  INTEGER :: i, j, k, m, ijkN, iPatch, ifl

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb
#ifdef RFLU
  TYPE(t_patch), POINTER  :: patch
#endif

  INTEGER :: turbModel, zofid, XCO, YCO, ZCO, errorFlag
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ijkC0, ijkC1
  REAL(RFREAL)          :: two3rd,cSmag,delFac2,delta2, zofFac,bzofFac
  REAL(RFREAL)          :: vorMag,volFace,rhoFace,muTurb,rhoMin,muTMin
  REAL(RFREAL)          :: sxx,sxy,sxz,syy,syz,szz
  REAL(RFREAL), POINTER :: fvol(:),fVar(:,:),cModel(:,:),mueT(:,:)
  REAL(RFREAL), POINTER :: tdv(:,:),sij(:,:),zof(:,:,:)
#ifdef RFLU
  REAL(RFREAL), POINTER :: bfVol(:),bfVar(:,:),bCModel(:,:),bMueT(:,:),bSij(:,:)
  REAL(RFREAL), POINTER :: bzof(:,:,:)
#endif

#ifdef RFLO
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend, iadd, jadd, kadd
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff, filterType
#endif
#ifdef RFLU
  INTEGER :: bcType
  INTEGER :: ifg,ifgt,ifgBeg,ifgtBeg
  INTEGER :: nPatches,nFaces,nFacesTot,nBFaces,nBFacesTot
  REAL(RFREAL) :: bValFactor
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesCalcEddyVis.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'Turb_LesCalcEddyVis',&
  'TURB_LesCalcEddyVis.F90' )

! get parameters and coefficients -------------------------------------------

  XCO        = XCOORD
  YCO        = YCOORD
  ZCO        = ZCOORD
  zofid      = ZOF_LES_EDDYVIS
  turbModel  = region%mixtInput%turbModel
  cSmag      = region%turbInput%cSmag
  delFac2    = region%turbInput%delFac2
#ifdef RFLO
  filterType = region%turbInput%filterType
  iLev       = region%currLevel
#endif

  rhoMin = 10000._RFREAL
  muTMin = 10000._RFREAL
  two3rd = 2._RFREAL/3._RFREAL

! allocate required arrays within this scope

#ifdef RFLO
  turb   => region%levels(ilev)%turb

  ALLOCATE( turb%coef(1,ibn:ien),stat=errorFlag )
  ALLOCATE( turb%fVar(CV_TURB_NELM,ibn:ien),stat=errorFlag )    
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#endif
#ifdef RFLU
  nPatches   = region%grid%nPatches
  nFaces     = region%grid%nFaces
  nFacesTot  = region%grid%nFacesTot
  nBFaces    = 0
  nBFacesTot = 0

  DO iPatch = 1,nPatches
    patch => region%patches(iPatch)

    nBFaces    = nBFaces    + patch%nBTris    + patch%nBQuads
    nBFacesTot = nBFacesTot + patch%nBTrisTot + patch%nBQuadsTot
  END DO ! iPatch

  IF (ijk /= DIRI) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'RFLU but ijk/=DIRI' )
  ENDIF
  turb   => region%turb

  ALLOCATE( turb%coef( 1,ibn:ien),stat=errorFlag )
  ALLOCATE( turb%bCoef(1,nBFaces),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%fVar( CV_TURB_NELM, nFaces),stat=errorFlag )    
  ALLOCATE( turb%bfVar(CV_TURB_NELM,nBFaces),stat=errorFlag )    
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#endif

! get model coefficients -----------------------------------------------

  IF (turbModel==TURB_MODEL_FIXSMAG) THEN
    
#ifdef RFLO
    CALL TURB_FloLesGenC2F( region,ijk )
    DO ijkN = ibn,ien
      turb%fVar(CV_TURB_DENS,ijkN) = 1._RFREAL/turb%fVar(CV_TURB_DENS,ijkN)
      turb%coef(1,ijkN) = cSmag*cSmag
    ENDDO
#endif
#ifdef RFLU
    CALL TURB_FluLesC2F( region )
    DO ijkN = ibn,ien
      turb%coef(1,ijkN)             = cSmag*cSmag
      turb%fVar(CV_TURB_DENS,ijkN)  = 1._RFREAL/turb%fVar(CV_TURB_DENS,ijkN)
    ENDDO
    DO ijkN = 1,nBFaces
      turb%bCoef(1,ijkN)            = cSmag*cSmag
      turb%bfVar(CV_TURB_DENS,ijkN) = 1._RFREAL/turb%bfVar(CV_TURB_DENS,ijkN)
    ENDDO
#endif
  ELSEIF (turbModel==TURB_MODEL_DYNSMAG) THEN   
    CALL TURB_LesCoefDynSmag( region,ibn,ien,ijk )
  ELSEIF (turbModel==TURB_MODEL_DYNMIXD) THEN   
    CALL TURB_LesCoefDynMixd( region,ibn,ien,ijk )
  ENDIF

! apply spatial tripping to eddy viscosity
    
#ifdef RFLO
  IF (ijk==DIRI) THEN
    zof => turb%zofi
  ELSEIF (ijk==DIRJ) THEN
    zof => turb%zofj
  ELSEIF (ijk==DIRK) THEN
    zof => turb%zofk
  ENDIF    
  DO ijkN = ibn,ien
    zofFac = zof(XCO,zofid,ijkN)*zof(YCO,zofid,ijkN)*zof(ZCO,zofid,ijkN)
    turb%coef(1,ijkN) = zofFac*turb%coef(1,ijkN)
  ENDDO
#endif
#ifdef RFLU
  zof  => turb%zofi
  bzof => turb%bZofi
  DO ijkN = ibn,ien
    zofFac = zof(XCO,zofid,ijkN)*zof(YCO,zofid,ijkN)*zof(ZCO,zofid,ijkN)
    turb%coef(1,ijkN)  = zofFac*turb%coef(1,ijkN)
  ENDDO
  DO ijkN = 1,nBFaces
    bzofFac = bzof(XCO,zofid,ijkN)*bzof(YCO,zofid,ijkN)*bzof(ZCO,zofid,ijkN)
    turb%bCoef(1,ijkN) = bzofFac*turb%bCoef(1,ijkN)
  ENDDO
#endif

! get dimensions and pointers

  tdv    => turb%dv
  fVar   => turb%fVar
  cModel => turb%coef
  mueT   => turb%mueT

#ifdef RFLO

  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

  IF (ijk==DIRI) THEN
! - pointers and indices for i-direction
    ibeg = ipcbeg-1
    iend = ipcend+2
    jbeg = jpcbeg-1
    jend = jpcend+1
    kbeg = kpcbeg-1
    kend = kpcend+1
    iadd = -1
    jadd = 0
    kadd = 0
    fvol => turb%fvolI
    sij  => turb%mISij
  ELSEIF (ijk==DIRJ) THEN
! - pointers and indices for j-direction
    ibeg = ipcbeg-1
    iend = ipcend+1
    jbeg = jpcbeg-1
    jend = jpcend+2
    kbeg = kpcbeg-1
    kend = kpcend+1
    iadd = 0
    jadd = -1
    kadd = 0
    fvol => turb%fvolJ
    sij  => turb%mJSij
  ELSEIF (ijk==DIRK) THEN
! - pointers and indices for k-direction
    ibeg = ipcbeg-1
    iend = ipcend+1
    jbeg = jpcbeg-1
    jend = jpcend+1
    kbeg = kpcbeg-1
    kend = kpcend+2
    iadd = 0
    jadd = 0
    kadd = -1
    fvol => turb%fvolK
    sij  => turb%mKSij
  END IF

! compute eddy viscosity; note that face density is stored in fVar as 1/rho
! which is leftover from lesLij; so we devide by fVar(1,:) instead of multiply

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend

        ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkC1 = ijkC0 + iadd + jadd*iCOff + kadd*ijCOff
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
#endif
#ifdef RFLU

  fvol => turb%fvolI
  sij  => turb%mISij

  DO ijkN =ibn,ien
#endif

        sxx = sij(E11,ijkN)
        sxy = sij(E12,ijkN)
        sxz = sij(E13,ijkN)

        syy = sij(E22,ijkN)
        syz = sij(E23,ijkN)

        szz = sij(E33,ijkN)

        vorMag = 0.5_RFREAL*(sxx*sxx+syy*syy+szz*szz)+sxy*sxy+sxz*sxz+syz*syz
        vorMag = SQRT(vorMag)

        volFace= fvol(ijkN)
        delta2 = volFace**two3rd  ! another option: delta2=3.*volFace**two3rd

!        rhoFace= 1._RFREAL/fVar(CV_TURB_DENS,ijkN)
        muTurb = cModel(1,ijkN)/fVar(CV_TURB_DENS,ijkN)*delFac2*delta2*vorMag

! Temporary clipping fix 
        muTurb = MAX(muTurb,REAL_SMALL)
!

!        rhoMin = MIN(rhoFace,rhoMin)
!        muTMin = MIN(muTurb ,muTMin)

        mueT(ijk,ijkN) = muTurb

#ifdef RFLO
        tdv(DV_TURB_CDYN,ijkC0) = tdv(DV_TURB_CDYN,ijkC0)+cModel(1,ijkN)
        tdv(DV_TURB_CDYN,ijkC1) = tdv(DV_TURB_CDYN,ijkC1)+cModel(1,ijkN)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
#endif
#ifdef RFLU
  ENDDO       ! ijkN

! sum model coefficient from faces; the division by nfaces is in GetTvCell
  DO ijkN = 1,region%grid%nFaces
    ijkC0 = region%grid%f2c(1,ijkN)
    ijkC1 = region%grid%f2c(2,ijkN)

    tdv(DV_TURB_CDYN,ijkC0) = tdv(DV_TURB_CDYN,ijkC0)+cModel(1,ijkN)
    tdv(DV_TURB_CDYN,ijkC1) = tdv(DV_TURB_CDYN,ijkC1)+cModel(1,ijkN)
  ENDDO  ! ijkN
#endif

! RFLU boundary treatment =====================================================

#ifdef RFLU
  bfVol   => turb%bfVolI
  bfVar   => turb%bfVar
  bCModel => turb%bCoef
  bMueT   => turb%bMueT
  bSij    => turb%bmISij

  DO iPatch = 1,region%grid%nPatches
    patch  => region%patches(iPatch)
! TEMPORARY : removing usage of bf2bg from everywhere
!    ifgBeg =  patch%bf2bg(BF2BG_BEG)
!    ifgtBeg=  patch%bf2bgTot(BF2BG_BEG)

    bcType = patch%bcType
    IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE .OR. &
        bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
      IF (turbModel==TURB_MODEL_FIXSMAG) THEN
        bValFactor = 1._RFREAL
      ELSE
        bValFactor = 0._RFREAL
      ENDIF ! turbModel
    ELSE
      bValFactor = 1._RFREAL
    ENDIF

! - sum model coefficient from bFaces; the division by nfaces is in GetTvCell
!   Note, bValFactor is needed for Rocflu to set eddy viscosity explicitly
!   to zero at injection and noslip-walls if desired.

    DO ifl = 1,patch%nBFaces
      ijkC0 = patch%bf2c(ifl)
      ifg   = ifl + ifgBeg-1

      sxx = bSij(E11,ifg)
      sxy = bSij(E12,ifg)
      sxz = bSij(E13,ifg)

      syy = bSij(E22,ifg)
      syz = bSij(E23,ifg)

      szz = bSij(E33,ifg)

      vorMag = 0.5_RFREAL*(sxx*sxx+syy*syy+szz*szz)+sxy*sxy+sxz*sxz+syz*syz
      vorMag = SQRT(vorMag)

      volFace= bfVol(ifg)
      delta2 = volFace**two3rd  ! another option: delta2=3.*volFace**two3rd

!      rhoFace= 1._RFREAL/bfVar(CV_TURB_DENS,ifg)
      muTurb = cModel(1,ifg)/bfVar(CV_TURB_DENS,ifg)*delFac2*delta2*vorMag
! Temporary clipping fix
      muTurb = MAX(muTurb,REAL_SMALL)
!

!      rhoMin = MIN(rhoFace,rhoMin)
!      muTMin = MIN(muTurb ,muTMin)

      bMueT(ijk,ifg) = bValFactor*muTurb
      tdv(DV_TURB_CDYN,ijkC0) = tdv(DV_TURB_CDYN,ijkC0)+bCModel(1,ifg)
    ENDDO    ! ifl
  ENDDO      ! iPatch
#endif

! check minimum density

  IF (rhoMin < 0._RFREAL) THEN
    write(*,*) 'TURB_LesCalcEddyVis: rhoMin < 0'
  ENDIF
  IF (muTMin < 0._RFREAL) THEN
    write(*,*) 'TURB_LesCalcEddyVis: muTMin < 0'
  ENDIF

! deallocate temporary arrays

  DEALLOCATE( turb%coef, turb%fVar )
#ifdef RFLU
  DEALLOCATE( turb%bCoef, turb%bfVar )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesCalcEddyVis

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_LesCalcEddyVis.F90,v $
! Revision 1.17  2009/08/12 04:15:59  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.16  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2006/08/19 15:40:59  mparmar
! Removed bf2bg, bf2bgTot
!
! Revision 1.13  2006/01/17 17:51:55  wasistho
! applied tripping to all eddy viscosity models
!
! Revision 1.12  2006/01/13 07:13:49  wasistho
! set rflu bValFactor=1 for FixSM and 0 else
!
! Revision 1.11  2006/01/13 06:53:42  wasistho
! set rflu bValFactor=1 for more stable sim
!
! Revision 1.10  2006/01/12 09:48:31  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.9  2006/01/01 00:11:03  wasistho
! multiplied muturb by bValFactor for rocflu
!
! Revision 1.8  2005/12/29 19:51:00  wasistho
! modified face indexing for sij
!
! Revision 1.7  2004/07/30 22:33:54  wasistho
! replaced floLesUniC2F by floLesGenC2F
!
! Revision 1.6  2004/05/28 02:00:17  wasistho
! update unstructured grid LES
!
! Revision 1.5  2004/05/17 20:47:39  wasistho
! compute first layer dummy mu_t instead of extrapolated
!
! Revision 1.4  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.3  2004/03/19 02:48:36  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2003/10/09 23:10:45  wasistho
! renamed CalcEddyVis to LesCalcEddyVis
!
! Revision 1.5  2003/08/29 01:42:51  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.4  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.3  2002/10/16 07:47:41  wasistho
! Enable Fix Smagorinsky
!
! Revision 1.2  2002/10/16 01:59:03  wasistho
! Changed global%error flag
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







