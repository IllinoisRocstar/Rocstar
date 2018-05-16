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
! Purpose: Collection of AUSM flux routines.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModAUSMFlux.F90,v 1.11 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModAUSMFlux

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_AUSM_ComputeFlux
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModAUSMFlux.F90,v $ $Revision: 1.11 $' 
              
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Wrapper function for AUSM flux functions.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux(pRegion)
                       
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================
  
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux',&
  'RFLU_ModAUSMFlux.F90')

! ******************************************************************************
! Call flux functions
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%gasModel ) 

! ==============================================================================
!   Thermally and calorically perfect gas
! ==============================================================================

    CASE ( GAS_MODEL_TCPERF )
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
        CASE ( DISCR_ORDER_1 ) 
          CALL RFLU_AUSM_ComputeFlux1(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_AUSM_ComputeFlux2_TCP(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder    

! ==============================================================================
!   Mixture of thermally and calorically perfect gases
! ==============================================================================

    CASE ( GAS_MODEL_MIXT_TCPERF ) 
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
        CASE ( DISCR_ORDER_1 ) 
          CALL RFLU_AUSM_ComputeFlux1(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_AUSM_ComputeFlux2_MTCP(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder

! ==============================================================================
!   Pseudo-gas
! ==============================================================================

    CASE ( GAS_MODEL_MIXT_PSEUDO ) 
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
! TO DO
!        CASE ( DISCR_ORDER_1 ) 
!          CALL RFLU_AUSM_ComputeFlux1_MPSD(pRegion)
! END TO DO 
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_AUSM_ComputeFlux2_MPSD(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder                           

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)             
  END SELECT ! pRegion%mixtInput%gasModel       

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux







! ******************************************************************************
!
! Purpose: Compute convective fluxes using AUSM+ scheme.
!
! Description: None.
!
! Input: 
!   nx          x-component of face normal
!   ny          y-component of face normal
!   nz          z-component of face normal
!   nm          Magnitude of face normal
!   fs          Face speed
!   rl          Density of left state
!   ul          x-component of velocity of left state
!   vl          y-component of velocity of left state
!   wl          z-component of velocity of left state   
!   Hl		Total enthalpy of left state
!   al		Speed of sound of left state
!   pl          Pressure of left state
!   rr          Density of right state
!   ur          x-component of velocity of right state
!   vr          y-component of velocity of right state
!   wr          z-component of velocity of right state  
!   pr          Pressure of right state
!   Hr		Total enthalpy of right state
!   ar		Speed of sound of right state
!
! Output: 
!   flx         Fluxes
!   vf          Face velocities
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSM+, AIAA Paper
!      95-1701, 1995
!   2. Do not use computation of face speed of sound which leads to exact 
!      capturing of isolated normal shock waves because of robustness problems
!      for unsteady flows and because that formulation is not applicable to 
!      anything but calorically and thermally perfect gases.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur, &
                                  vr,wr,pr,Hr,ar,flx,vf)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), INTENT(IN) :: al,ar,fs,Hl,Hr,nm,nx,ny,nz,pl,pr,rl,rr,ul,ur,vl, &
                              vr,wl,wr
  REAL(RFREAL), INTENT(OUT) :: flx(5),vf(3)

! ==============================================================================
! Locals
! ==============================================================================

  REAL(RFREAL) :: af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,pf,ql,qr,vml,vmr, &
                  wtl,wtr

! ******************************************************************************
! Start, compute face state
! ******************************************************************************

  ql = ul*nx + vl*ny + wl*nz - fs
  qr = ur*nx + vr*ny + wr*nz - fs
    
  af = 0.5_RFREAL*(al+ar) ! NOTE not using original formulation, see note

  ml  = ql/af
  mla = ABS(ml)

  mr  = qr/af
  mra = ABS(mr)    

  IF ( mla <= 1.0_RFREAL ) THEN 
    mlp = 0.25_RFREAL*(ml+1.0_RFREAL)*(ml+1.0_RFREAL) & 
        + 0.125_RFREAL*(ml*ml-1.0_RFREAL)*(ml*ml-1.0_RFREAL)
    wtl = 0.25_RFREAL*(ml+1.0_RFREAL)*(ml+1.0_RFREAL)*(2.0_RFREAL-ml) & 
        + 0.1875_RFREAL*ml*(ml*ml-1.0_RFREAL)*(ml*ml-1.0_RFREAL)
  ELSE
    mlp = 0.5_RFREAL*(ml+mla)
    wtl = 0.5_RFREAL*(1.0_RFREAL+ml/mla)
  END IF ! mla

  IF ( mra <= 1.0_RFREAL ) THEN 
    mrm = -0.25_RFREAL*(mr-1.0_RFREAL)*(mr-1.0_RFREAL) & 
          -0.125_RFREAL*(mr*mr-1.0_RFREAL)*(mr*mr-1.0_RFREAL)
    wtr = 0.25_RFREAL*(mr-1.0_RFREAL)*(mr-1.0_RFREAL)*(2.0_RFREAL+mr) & 
        - 0.1875_RFREAL*mr*(mr*mr-1.0_RFREAL)*(mr*mr-1.0_RFREAL)
  ELSE
    mrm = 0.5_RFREAL*(mr-mra)
    wtr = 0.5_RFREAL*(1.0_RFREAL-mr/mra)
  END IF ! mla

  mf  = mlp + mrm
  mfa = ABS(mf)
  mfp = 0.5_RFREAL*(mf+mfa)
  mfm = 0.5_RFREAL*(mf-mfa)

  pf = wtl*pl + wtr*pr 

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

  vf(1) = mfp*ul + mfm*ur
  vf(2) = mfp*vl + mfm*vr
  vf(3) = mfp*wl + mfm*wr    

  flx(1) = (af*(mfp*rl    + mfm*rr   )        )*nm
  flx(2) = (af*(mfp*rl*ul + mfm*rr*ur) + pf*nx)*nm
  flx(3) = (af*(mfp*rl*vl + mfm*rr*vr) + pf*ny)*nm
  flx(4) = (af*(mfp*rl*wl + mfm*rr*wr) + pf*nz)*nm
  flx(5) = (af*(mfp*rl*Hl + mfm*rr*Hr) + pf*fs)*nm

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_AUSM_FluxFunction









! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate AUSM+ scheme.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine can also be used for a mixture of thermally and calorically 
!      perfect gases because, the face states being identical to the cell 
!      states, the mixture properties are constant.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux1(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_Ho_CpTUVW, & 
                           RFLU_CentralFirstPatch

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cpl,cpr,fs,Hl,Hr,irl,irr,nm,nx,ny,nz,pl,pr,rl,rr, &
                  tl,tr,ul,ur,vl,vr,wl,wr
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux1',&
  'RFLU_ModAUSMFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux1")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pRegion%grid%indGs

  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol
  indSd  = pRegion%mixtInput%indSd

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    cpl = pGv(GV_MIXT_CP,indCp*c1)
    
    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = pCv(CV_MIXT_XMOM,c1)*irl
    vl  = pCv(CV_MIXT_YMOM,c1)*irl
    wl  = pCv(CV_MIXT_ZMOM,c1)*irl
    
    pl  = pDv(DV_MIXT_PRES,c1)
    tl  = pDv(DV_MIXT_TEMP,c1)
    al  = pDv(DV_MIXT_SOUN,c1)

    Hl  = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    cpr = pGv(GV_MIXT_CP,indCp*c2)

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur  = pCv(CV_MIXT_XMOM,c2)*irr
    vr  = pCv(CV_MIXT_YMOM,c2)*irr
    wr  = pCv(CV_MIXT_ZMOM,c2)*irr
    
    pr  = pDv(DV_MIXT_PRES,c2)
    tr  = pDv(DV_MIXT_TEMP,c2)    
    ar  = pDv(DV_MIXT_SOUN,c2)
    
    Hr  = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf)

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flx(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    CALL RFLU_CentralFirstPatch(pRegion,pPatch)
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux1")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux1







! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSM+ scheme 
!   for mixture of thermally and calorically perfect gaseous species and 
!   particulate species.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux2_MPSD(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_G_CpR, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_M_R, & 
                           MixtPerf_R_M, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, &
                           RFLU_CentralSecondPatch

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cpl,cpr,dx,dy,dz,fs,gcl,gcr,gl,gr,Hl,Hr,irl,irr,mml, &
                  mmr,nm,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef SPEC
  INTEGER :: iSpec
  REAL(RFREAL) :: gcg,immg,mmg,phip,Yg,Y1,Y2
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec 
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGcSpec 
  TYPE(t_spec_type), POINTER :: pSpecType
#endif  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux2_MPSD',&
  'RFLU_ModAUSMFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux2_MPSD")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pGrid%indGs
  
  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol 
  indSd  = pRegion%mixtInput%indSd
  
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

#ifdef SPEC
  pCvSpec => pRegion%spec%cv  
  pGcSpec => pRegion%spec%gradCell    
#endif

  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_PSEUDO ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  IF ( (indCp /= 1) .OR. (indMol /= 1) ) THEN 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END IF ! indCp

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)

#ifdef SPEC
    immg = 0.0_RFREAL
    cpl  = 0.0_RFREAL
    Yg   = 0.0_RFREAL
    phip = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y1 = irl*pCvSpec(iSpec,c1) + pGcSpec(XCOORD,iSpec,c1)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c1)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c1)*dz

      IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
        cpl  = cpl  +    Y1*pSpecType%pMaterial%spht
        phip = phip + rl*Y1/pSpecType%pMaterial%dens
      ELSE 
        immg = immg + Y1/pSpecType%pMaterial%molw
        cpl  = cpl  + Y1*pSpecType%pMaterial%spht
        Yg   = Yg   + Y1                                                              
      END IF ! pSpecType%discreteFlag
    END DO ! iSpec

    mmg = Yg/immg
    gcg = MixtPerf_R_M(mmg)
    gcl = Yg/(1.0_RFREAL-phip)*gcg
    mml = MixtPerf_M_R(gcl)
#endif                    

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)
        
    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz    

    tl = MixtPerf_T_DPR(rl,pl,gcl)
    al = MixtPerf_C_GRT(gl,gcl,tl)
    Hl = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)
    
! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)

#ifdef SPEC
    immg = 0.0_RFREAL
    cpr  = 0.0_RFREAL
    Yg   = 0.0_RFREAL
    phip = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y2 = irr*pCvSpec(iSpec,c2) + pGcSpec(XCOORD,iSpec,c2)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c2)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c2)*dz

      IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
        cpr  = cpr  +    Y2*pSpecType%pMaterial%spht
        phip = phip + rr*Y2/pSpecType%pMaterial%dens
      ELSE 
        immg = immg + Y2/pSpecType%pMaterial%molw
        cpr  = cpr  + Y2*pSpecType%pMaterial%spht
        Yg   = Yg   + Y2                                                              
      END IF ! pSpecType%discreteFlag
    END DO ! iSpec

    mmg = Yg/immg
    gcg = MixtPerf_R_M(mmg)
    gcr = Yg/(1.0_RFREAL-phip)*gcg
    mmr = MixtPerf_M_R(gcr)
#endif  
    
    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)
        
    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz

    tr = MixtPerf_T_DPR(rr,pr,gcr)
    ar = MixtPerf_C_GRT(gr,gcr,tr)
    Hr = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)
    
! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    gl = cpl/(cpl-gcl)
    gr = cpr/(cpr-gcr)
    
    IF ( ABS(gr-gl) < 0.05_RFREAL ) THEN 
      CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                                  rr,ur,vr,wr,pr,Hr,ar,flx,vf)

      pMf(indMf*ifg) = flx(1)

      pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
      pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
      pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
      pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
      pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

      pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
      pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
      pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
      pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
      pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)

      pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
      pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
      pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)

      pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
      pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
      pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)
    ELSE 
      CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                                  rr,ur,vr,wr,pr,Hr,ar,flx,vf)

      pMf(indMf*ifg) = 0.5_RFREAL*flx(1)

      pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
      pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
      pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
      pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
      pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

      pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
      pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
      pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
      
      CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                                  rr,ur,vr,wr,pr,Hr,ar,flx,vf)

      pMf(indMf*ifg) = pMf(indMf*ifg) + 0.5_RFREAL*flx(1)

      pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
      pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
      pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
      pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
      pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)            
      
      pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
      pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
      pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)      
    END IF ! ABS(gr-gl)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux2_MPSD")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux2_MPSD









! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSM+ scheme 
!   for mixture of thermally and calorically perfect gases.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux2_MTCP(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_G_CpR, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_R_M, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, & 
                           RFLU_CentralSecondPatch

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cpl,cpr,dx,dy,dz,fs,gcl,gcr,gl,gr,Hl,Hr,irl,irr, &
                  nm,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef SPEC
  INTEGER :: iSpec
  REAL(RFREAL) :: mml,mmr,Y1,Y2
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec 
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGcSpec 
  TYPE(t_spec_type), POINTER :: pSpecType
#endif  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux2_MTCP',&
  'RFLU_ModAUSMFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux2_MTCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pGrid%indGs
  
  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol 
  indSd  = pRegion%mixtInput%indSd
  
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

#ifdef SPEC
  pCvSpec => pRegion%spec%cv  
  pGcSpec => pRegion%spec%gradCell    
#endif

  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_TCPERF ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  IF ( (indCp /= 1) .OR. (indMol /= 1) ) THEN 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END IF ! indCp

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)

#ifdef SPEC
    mml = 0.0_RFREAL
    cpl = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y1 = irl*pCvSpec(iSpec,c1) + pGcSpec(XCOORD,iSpec,c1)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c1)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c1)*dz

      mml = mml + Y1/pSpecType%pMaterial%molw
      cpl = cpl + Y1*pSpecType%pMaterial%spht
    END DO ! iSpec
    
    mml = 1.0_RFREAL/mml
    gcl = MixtPerf_R_M(mml)
    gl  = MixtPerf_G_CpR(cpl,gcl)    
#endif                   

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)
        
    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz    

    tl = MixtPerf_T_DPR(rl,pl,gcl)
    al = MixtPerf_C_GRT(gl,gcl,tl)
    Hl = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)

#ifdef SPEC
    mmr = 0.0_RFREAL
    cpr = 0.0_RFREAL


    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y2 = irr*pCvSpec(iSpec,c2) + pGcSpec(XCOORD,iSpec,c2)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c2)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c2)*dz


      mmr = mmr + Y2/pSpecType%pMaterial%molw
      cpr = cpr + Y2*pSpecType%pMaterial%spht
    END DO ! iSpec
    
    mmr = 1.0_RFREAL/mmr 
    gcr = MixtPerf_R_M(mmr)
    gr  = MixtPerf_G_CpR(cpr,gcr)           
#endif
              
    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)
        
    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz

    tr = MixtPerf_T_DPR(rr,pr,gcr)
    ar = MixtPerf_C_GRT(gr,gcr,tr)
    Hr = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf)

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flx(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux2_MTCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux2_MTCP








! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSM+ scheme
!   for thermally and calorically perfect gas.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!  1. Only applicable to thermally and calorically perfect gas because would
!     need to recompute mixture properties at faces depending on face states.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux2_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_R_CpG, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, & 
                           RFLU_CentralSecondPatch

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indGs,indMf,indSd,iPatch
  REAL(RFREAL) :: al,ar,cp,dx,dy,dz,fs,gc,g,Hl,Hr,irl,irr,nm,nx,ny,nz,pl,pr, &
                  rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux2_TCP',&
  'RFLU_ModAUSMFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux2_TCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs = pGrid%indGs
  
  indMf = pRegion%mixtInput%indMfMixt  
  indSd = pRegion%mixtInput%indSd
  
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_TCPERF ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  cp = global%refCp
  g  = global%refGamma  
  gc = MixtPerf_R_CpG(cp,g)    

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = pCv(CV_MIXT_XMOM,c1)*irl
    vl  = pCv(CV_MIXT_YMOM,c1)*irl
    wl  = pCv(CV_MIXT_ZMOM,c1)*irl
    pl  = pDv(DV_MIXT_PRES,c1)
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)
        
    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz

    tl = MixtPerf_T_DPR(rl,pl,gc)
    al = MixtPerf_C_GRT(g,gc,tl)
    Hl = MixtPerf_Ho_CpTUVW(cp,tl,ul,vl,wl)
    
! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur  = pCv(CV_MIXT_XMOM,c2)*irr
    vr  = pCv(CV_MIXT_YMOM,c2)*irr
    wr  = pCv(CV_MIXT_ZMOM,c2)*irr
    pr  = pDv(DV_MIXT_PRES,c2)

    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)
        
    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz
    
    tr = MixtPerf_T_DPR(rr,pr,gc)
    ar = MixtPerf_C_GRT(g,gc,tr)
    Hr = MixtPerf_Ho_CpTUVW(cp,tr,ur,vr,wr)
    
! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, & 
                                wr,pr,Hr,ar,flx,vf)

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flx(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux2_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux2_TCP






  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModAUSMFlux


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModAUSMFlux.F90,v $
! Revision 1.11  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:31  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/05/01 22:20:28  haselbac
! Removed debug statements
!
! Revision 1.8  2006/05/01 21:00:47  haselbac
! Rewrite for consistency and cleanliness
!
! Revision 1.7  2006/04/15 16:58:42  haselbac
! Added capability of running 1st order boundary fluxes with 2nd order volume 
! fluxes
!
! Revision 1.6  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.5  2005/11/15 17:25:07  haselbac
! Bug fix: Missing declarations lead to compilation failure without SPEC=1
!
! Revision 1.4  2005/11/14 16:58:05  haselbac
! Added flux for pseudo-gas, reordered routines
!
! Revision 1.3  2005/11/10 02:26:06  haselbac
! Complete rewrite to allow for computations with variable properties
!
! Revision 1.2  2005/07/19 20:06:29  haselbac
! Modified af to prevent crashes for Skews problem
!
! Revision 1.1  2005/07/14 21:39:02  haselbac
! Initial revision
!
! ******************************************************************************
  
  











