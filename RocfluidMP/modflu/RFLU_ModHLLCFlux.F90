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
! Purpose: Collection of HLLC flux routines.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModHLLCFlux.F90,v 1.10 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModHLLCFlux

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
  PUBLIC :: RFLU_HLLC_ComputeFlux
   
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModHLLCFlux.F90,v $ $Revision: 1.10 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Wrapper function for HLLC flux functions.
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

SUBROUTINE RFLU_HLLC_ComputeFlux(pRegion)
                       
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

  CALL RegisterFunction(global,'RFLU_HLLC_ComputeFlux',&
  'RFLU_ModHLLCFlux.F90')

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
          CALL RFLU_HLLC_ComputeFlux1_TCP(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_HLLC_ComputeFlux2_TCP(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder    

! ==============================================================================
!   Mixture of thermally and calorically perfect gases
! ==============================================================================

    CASE ( GAS_MODEL_MIXT_TCPERF ) 
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
        CASE ( DISCR_ORDER_1 ) 
          CALL RFLU_HLLC_ComputeFlux1_MTCP(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_HLLC_ComputeFlux2_MTCP(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder

! ==============================================================================
!   Mixture of gas, vapor, and liquid
! ==============================================================================

    CASE ( GAS_MODEL_MIXT_GASLIQ ) 
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
        CASE ( DISCR_ORDER_1 ) 
          CALL RFLU_HLLC_ComputeFlux1_GL(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_HLLC_ComputeFlux2_GL(pRegion)               
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

END SUBROUTINE RFLU_HLLC_ComputeFlux





! ******************************************************************************
!
! Purpose: Compute central convective fluxes for multiphase mixture using 
!   first-order accurate version of HLLC scheme.
!
! Description: None.
!
! Input: 
!   pRegion     Data of current region.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_HLLC_ComputeFlux1_GL(pRegion)

  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtLiq_C2_Bp, &
                           MixtPerf_C2_GRT, &
                           MixtPerf_Cv_CpR, &
                           MixtPerf_D_PRT, &
                           MixtPerf_R_CpG, &
                           MixtPerf_R_M, &                
                           RFLU_CentralFirstPatch_GL 
                       
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
  INTEGER, DIMENSION(:,:), POINTER :: f2c
  REAL(RFREAL) :: Bgh2,Blh2,Bp,Bt,Bvh2,cml,cmr,cms,cvg,cvl,cvv,Cgh2,Clh2, &
                  Cvh2,Cvm,el,er,fs,irl,irr,nm,nx,ny,nz,ph,pl,po,pr,ps,ql, &
                  qr,qs,res,rgh,rgl,rgr,rh,rl,rlh,ro,rr,rs,rus,rvh,rvl,rvr, &
                  rvs,rws,rYgh,rYgl,rYgr,rYvh,rYvl,rYvr,Rg,Rv,sl,sm,sr, &
                  term,th,tl,to,tr,ul,ur,us,vfgh,vfgl,vfgr,vflh,vfvh,vfvl, &
                  vfvr,vl,vr,vs,Vms,wl,wr,ws,wt1,wt2
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pDvMixt,pRhs,pSd
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_patch), POINTER :: pPatch  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_HLLC_ComputeFlux1_GL',&
  'RFLU_ModHLLCFlux.F90')

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pGrid%indGs
  
  indMf  = pRegion%mixtInput%indMfMixt  
  indSd  = pRegion%mixtInput%indSd
  
  pCvMixt => pRegion%mixt%cv
  pDvMixt => pRegion%mixt%dv
  pCvSpec => pRegion%spec%cv
    
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd 

  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel
  
! ******************************************************************************
! Define constants
! ******************************************************************************

  ro  = global%refDensityLiq
  po  = global%refPressLiq
  to  = global%refTempLiq
  Bp  = global%refBetaPLiq
  Bt  = global%refBetaTLiq
  cvl = global%refCvLiq

  Rg  = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
  cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht,Rg)

  Rv  = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
  cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht,Rv)

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
 
! =============================================================================   
!   Compute left and right states     
! =============================================================================

! -----------------------------------------------------------------------------
!   Left state
! -----------------------------------------------------------------------------    
  
    rl   = pCvMixt(CV_MIXT_DENS,c1)
    irl  = 1.0_RFREAL/rl
    
    ul   = pCvMixt(CV_MIXT_XMOM,c1)*irl
    vl   = pCvMixt(CV_MIXT_YMOM,c1)*irl
    wl   = pCvMixt(CV_MIXT_ZMOM,c1)*irl
    el   = pCvMixt(CV_MIXT_ENER,c1)*irl
    pl   = pDvMixt(DV_MIXT_PRES,c1)
    tl   = pDvMixt(DV_MIXT_TEMP,c1)
    cml  = pDvMixt(DV_MIXT_SOUN,c1)
    
    rYgl = pCvSpec(1,c1)
    rYvl = pCvSpec(2,c1)   
     
    rvl  = MixtPerf_D_PRT(pl,Rv,tl)
    rgl  = MixtPerf_D_PRT(pl,Rg,tl)  
      
    vfvl = rYvl/rvl
    vfgl = rYgl/rgl

    ql  = ul*nx + vl*ny + wl*nz - fs

! -----------------------------------------------------------------------------
!   Right state
! -----------------------------------------------------------------------------    
    
    rr   = pCvMixt(CV_MIXT_DENS,c2)
    irr  = 1.0_RFREAL/rr
    
    ur   = pCvMixt(CV_MIXT_XMOM,c2)*irr
    vr   = pCvMixt(CV_MIXT_YMOM,c2)*irr
    wr   = pCvMixt(CV_MIXT_ZMOM,c2)*irr
    er   = pCvMixt(CV_MIXT_ENER,c2)*irr
    pr   = pDvMixt(DV_MIXT_PRES,c2)
    tr   = pDvMixt(DV_MIXT_TEMP,c2)
    cmr  = pDvMixt(DV_MIXT_SOUN,c2)
 
    rYgr = pCvSpec(1,c2)
    rYvr = pCvSpec(2,c2)   
    
    rvr  = MixtPerf_D_PRT(pr,Rv,tr)
    rgr  = MixtPerf_D_PRT(pr,Rg,tr) 
               
    vfvr = rYvr/rvr
    vfgr = rYgr/rgr
 
    qr  = ur*nx + vr*ny + wr*nz - fs
   
! ==============================================================================
!   Compute weights
! ==============================================================================

    wt1 = SQRT(rr/rl)
    wt2 = 1.0_RFREAL/(1.0_RFREAL + wt1)   
   
! ==============================================================================
!   Compute averages
! ==============================================================================

    us  = wt2*(ul + wt1*ur)
    vs  = wt2*(vl + wt1*vr)
    ws  = wt2*(wl + wt1*wr)
    qs  = wt2*(ql + wt1*qr)

    Vms = us*us + vs*vs + ws*ws

    rh   = wt1*rl
    ph   = wt2*(pl   + wt1*pr)  
    th   = wt2*(tl   + wt1*tr)  
    rYgh = wt2*(rYgl + wt1*rYgr)
    rYvh = wt2*(rYvl + wt1*rYvr)    

    rlh = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,ph,po,th,to)
    rvh = MixtPerf_D_PRT(ph,Rv,th)
    rgh = MixtPerf_D_PRT(ph,Rg,th)

    vfgh = rYgh/rgh
    vfvh = rYvh/rvh
  
    IF ( vfgh > 1.0_RFREAL) THEN
      vfgh = 1.0_RFREAL
      vfvh = 0.0_RFREAL
    ELSE IF ( vfgh < 0.0_RFREAL) THEN
      vfgh = 0.0_RFREAL
    ELSE IF ( vfvh > 1.0_RFREAL) THEN
      vfvh = 1.0_RFREAL
      vfgh = 0.0_RFREAL
    ELSE IF ( vfvh < 0.0_RFREAL) THEN
      vfvh = 0.0_RFREAL
    END IF ! vfgh

    vflh = 1.0_RFREAL - vfvh - vfgh

    IF ( vflh > 1.0_RFREAL ) THEN
      vflh = 1.0_RFREAL
    ELSE IF ( vflh < 0.0_RFREAL) THEN
      vflh = 0.0_RFREAL
    END IF ! vflh  

    Clh2 = MixtLiq_C2_Bp(Bp) 
    Cvh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,th)
    Cgh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,th)

    Blh2 = -Bt/Bp
    Bvh2 = rvh*Rv
    Bgh2 = rgh*Rg

    Cvm = (rlh*vflh*cvl + rvh*vfvh*cvv + rgh*vfgh*cvg)/rh
 
    cms = MixtGasLiq_C(Cvm,rh,ph,rlh,rvh,rgh,vflh,vfvh,vfgh,Clh2,Cvh2,Cgh2, &
                       Blh2,Bvh2,Bgh2)

    sl = MIN(ql-cml,qs-cms)
    sr = MAX(qr+cmr,qs+cms) 

! ==============================================================================
!   Compute fluxes
! ==============================================================================
      
    IF ( sl > 0.0_RFREAL ) THEN 
      flx(1) =  ql* rl                    *nm
      flx(2) = (ql* rl*ul + pl*nx        )*nm
      flx(3) = (ql* rl*vl + pl*ny        )*nm
      flx(4) = (ql* rl*wl + pl*nz        )*nm
      flx(5) = (ql*(rl*el + pl)   + pl*fs)*nm      
    ELSE IF ( sr < 0.0_RFREAL ) THEN  
      flx(1) =  qr* rr                    *nm
      flx(2) = (qr* rr*ur + pr*nx        )*nm
      flx(3) = (qr* rr*vr + pr*ny        )*nm
      flx(4) = (qr* rr*wr + pr*nz        )*nm
      flx(5) = (qr*(rr*er + pr)   + pr*fs)*nm          
    ELSE 
      sm = (rr*qr*(sr-qr) - rl*ql*(sl-ql) - pr + pl)/(rr*(sr-qr) - rl*(sl-ql))
      ps = pl + rl*(ql-sl)*(ql-sm) 
    
      IF ( sm > 0.0_RFREAL ) THEN 
        term = 1.0_RFREAL/(sl-sm)
        
        rs  = term* (sl-ql)*rl
        rus = term*((sl-ql)*rl*ul + (ps    - pl   )*nx)
        rvs = term*((sl-ql)*rl*vl + (ps    - pl   )*ny)
        rws = term*((sl-ql)*rl*wl + (ps    - pl   )*nz) 
        res = term*((sl-ql)*rl*el + (ps*sm - pl*ql)   )        
      ELSE 
        term = 1.0_RFREAL/(sr-sm)
        
        rs  = term* (sr-qr)*rr
        rus = term*((sr-qr)*rr*ur + (ps    - pr   )*nx)
        rvs = term*((sr-qr)*rr*vr + (ps    - pr   )*ny)
        rws = term*((sr-qr)*rr*wr + (ps    - pr   )*nz) 
        res = term*((sr-qr)*rr*er + (ps*sm - pr*qr)   )        
      END IF ! sm
      
      flx(1) =  sm* rs                  *nm
      flx(2) = (sm* rus + ps*nx        )*nm
      flx(3) = (sm* rvs + ps*ny        )*nm
      flx(4) = (sm* rws + ps*nz        )*nm 
      flx(5) = (sm*(res + ps)   + ps*fs)*nm      
    END IF ! sl
   
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

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + us*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + ws*pMf(indMf*ifg)
   
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - us*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - ws*pMf(indMf*ifg)
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    CALL RFLU_CentralFirstPatch_GL(pRegion,pPatch)
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_HLLC_ComputeFlux1_GL








! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate version of
!   HLLC scheme for mixture of thermally and calorically perfect gases.
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

SUBROUTINE RFLU_HLLC_ComputeFlux1_MTCP(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtPerf_C_DGP, &
                           MixtPerf_C_GHoVm2, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_R_M, & 
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
  REAL(RFREAL) :: al,ar,as,cpl,cpr,el,er,fs,gcl,gcr,gl,gr,gs,Hl,Hr,Hs,irl, &
                  irr,mml,mmr,nm,nx,ny,nz,ql,qr,qs,pl,pr,ps,rl,rr,res,rs, &
                  rus,rvs,rws,sl,sm,sr,term,ul,ur,us,vl,Vml,Vmr,Vms,vr,vs, &
                  wl,wr,ws,wt1,wt2
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef SPEC
  INTEGER :: iSpec
  REAL(RFREAL) :: cps,gcs,mms,Y1,Y2,Ys
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec 
  TYPE(t_spec_type), POINTER :: pSpecType
#endif  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_HLLC_ComputeFlux1_MTCP',&
  'RFLU_ModHLLCFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::HLLC_ComputeFlux1_MTCP")
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
  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd 

#ifdef SPEC
  pCvSpec => pRegion%spec%cv     
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

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    mml = pGv(GV_MIXT_MOL,c1)
    cpl = pGv(GV_MIXT_CP ,c1)
    
    gcl = MixtPerf_R_M(mml)
    gl  = MixtPerf_G_CpR(cpl,gcl)    

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)
    al = pDv(DV_MIXT_SOUN,c1)    

    el  = MixtPerf_Eo_DGPUVW(rl,gl,pl,ul,vl,wl)

    Vml = ul*ul + vl*vl + wl*wl
    ql  = ul*nx + vl*ny + wl*nz - fs
    Hl  = el + irl*pl

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    mmr = pGv(GV_MIXT_MOL,c2)
    cpr = pGv(GV_MIXT_CP ,c2)
    
    gcr = MixtPerf_R_M(mmr)
    gr  = MixtPerf_G_CpR(cpr,gcr)           
    
    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)
    ar = pDv(DV_MIXT_SOUN,c2)

    er  = MixtPerf_Eo_DGPUVW(rr,gr,pr,ur,vr,wr)

    Vmr = ur*ur + vr*vr + wr*wr
    qr  = ur*nx + vr*ny + wr*nz - fs
    Hr  = er + irr*pr
    
! ==============================================================================
!   Compute weights
! ==============================================================================

    wt1 = SQRT(rr/rl)
    wt2 = 1.0_RFREAL/(1.0_RFREAL + wt1)

! ==============================================================================
!   Compute averages. NOTE this average differs from that of the TCP case in 
!   that need to recompute gas properties to compute speed of sound. Choose to
!   apply Roe-average also to mass fractions, although there is no analysis to
!   back this up. It is not clear what is the best way to do this.
! ==============================================================================

    us = wt2*(ul + wt1*ur)
    vs = wt2*(vl + wt1*vr)
    ws = wt2*(wl + wt1*wr)
    qs = wt2*(ql + wt1*qr)
    Hs = wt2*(Hl + wt1*Hr)    

#ifdef SPEC
    mms = 0.0_RFREAL
    cps = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y1 = irl*pCvSpec(iSpec,c1) 
      Y2 = irr*pCvSpec(iSpec,c2)

      Ys = wt2*(Y1 + wt1*Y2)

      mms = mms + Ys/pSpecType%pMaterial%molw
      cps = cps + Ys*pSpecType%pMaterial%spht
    END DO ! iSpec
    
    mms = 1.0_RFREAL/mms 
    gcs = MixtPerf_R_M(mms)
    gs  = MixtPerf_G_CpR(cps,gcs)              
#endif

    Vms = us*us + vs*vs + ws*ws
    as  = MixtPerf_C_GHoVm2(gs,Hs,Vms)    
    sl  = MIN(ql-al,qs-as)
    sr  = MAX(qr+ar,qs+as)

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    IF ( sl > 0.0_RFREAL ) THEN 
      flx(1) =  ql* rl                    *nm
      flx(2) = (ql* rl*ul + pl*nx        )*nm
      flx(3) = (ql* rl*vl + pl*ny        )*nm
      flx(4) = (ql* rl*wl + pl*nz        )*nm
      flx(5) = (ql*(rl*el + pl)   + pl*fs)*nm
    ELSE IF ( sr < 0.0_RFREAL ) THEN  
      flx(1) =  qr* rr                    *nm
      flx(2) = (qr* rr*ur + pr*nx        )*nm
      flx(3) = (qr* rr*vr + pr*ny        )*nm
      flx(4) = (qr* rr*wr + pr*nz        )*nm
      flx(5) = (qr*(rr*er + pr)   + pr*fs)*nm    
    ELSE 
      sm = (rr*qr*(sr-qr) - rl*ql*(sl-ql) - pr + pl)/(rr*(sr-qr) - rl*(sl-ql))
      ps = pl + rl*(ql-sl)*(ql-sm) 
    
      IF ( sm > 0.0_RFREAL ) THEN 
        term = 1.0_RFREAL/(sl-sm)
        
        rs  = term* (sl-ql)*rl
        rus = term*((sl-ql)*rl*ul + (ps    - pl   )*nx)
        rvs = term*((sl-ql)*rl*vl + (ps    - pl   )*ny)
        rws = term*((sl-ql)*rl*wl + (ps    - pl   )*nz) 
        res = term*((sl-ql)*rl*el + (ps*sm - pl*ql)   )               
      ELSE 
        term = 1.0_RFREAL/(sr-sm)
        
        rs  = term* (sr-qr)*rr
        rus = term*((sr-qr)*rr*ur + (ps    - pr   )*nx)
        rvs = term*((sr-qr)*rr*vr + (ps    - pr   )*ny)
        rws = term*((sr-qr)*rr*wr + (ps    - pr   )*nz) 
        res = term*((sr-qr)*rr*er + (ps*sm - pr*qr)   )      
      END IF ! sm
      
      flx(1) =  sm* rs                  *nm
      flx(2) = (sm* rus + ps*nx        )*nm
      flx(3) = (sm* rvs + ps*ny        )*nm
      flx(4) = (sm* rws + ps*nz        )*nm 
      flx(5) = (sm*(res + ps)   - ps*fs)*nm           
    END IF ! sl

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

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + us*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + ws*pMf(indMf*ifg)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - us*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - ws*pMf(indMf*ifg)
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
  CALL FPROFILER_ENDS("RFLU::HLLC_ComputeFlux1_MTCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_HLLC_ComputeFlux1_MTCP








! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate version of 
!   HLLC scheme for thermally and calorically perfect gas.
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

SUBROUTINE RFLU_HLLC_ComputeFlux1_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GHoVm2, & 
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

  INTEGER :: c1,c2,ifg,indGs,indMf,indSd,iPatch
  INTEGER, DIMENSION(:,:), POINTER :: f2c
  REAL(RFREAL) :: al,ar,as,el,er,fs,g,Hl,Hr,Hs,irl,irr,nm,nx,ny,nz,ql,qr, &
                  qs,pl,pr,ps,rl,rr,res,rs,rus,rvs,rws,sl,sm,sr,term,ul,ur, &
                  us,vl,Vml,Vmr,Vms,vr,vs,wl,wr,ws,wt1,wt2
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pRhs,pSd
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_HLLC_ComputeFlux1_TCP',&
  'RFLU_ModHLLCFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::HLLC_ComputeFlux1_TCP")
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
  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd 

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_TCPERF ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  g = global%refGamma

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

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = pCv(CV_MIXT_XMOM,c1)*irl
    vl  = pCv(CV_MIXT_YMOM,c1)*irl
    wl  = pCv(CV_MIXT_ZMOM,c1)*irl
    el  = pCv(CV_MIXT_ENER,c1)*irl
    pl  = pDv(DV_MIXT_PRES,c1)
    al  = pDv(DV_MIXT_SOUN,c1)

    Vml = ul*ul + vl*vl + wl*wl
    ql  = ul*nx + vl*ny + wl*nz - fs
    Hl  = el + irl*pl

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur  = pCv(CV_MIXT_XMOM,c2)*irr
    vr  = pCv(CV_MIXT_YMOM,c2)*irr
    wr  = pCv(CV_MIXT_ZMOM,c2)*irr
    er  = pCv(CV_MIXT_ENER,c2)*irr
    pr  = pDv(DV_MIXT_PRES,c2)
    ar  = pDv(DV_MIXT_SOUN,c2)

    Vmr = ur*ur + vr*vr + wr*wr
    qr  = ur*nx + vr*ny + wr*nz - fs
    Hr  = er + irr*pr

! ==============================================================================
!   Compute weights
! ==============================================================================

    wt1 = SQRT(rr/rl)
    wt2 = 1.0_RFREAL/(1.0_RFREAL + wt1)

! ==============================================================================
!   Compute averages
! ==============================================================================

    us  = wt2*(ul + wt1*ur)
    vs  = wt2*(vl + wt1*vr)
    ws  = wt2*(wl + wt1*wr)
    qs  = wt2*(ql + wt1*qr)
    Hs  = wt2*(Hl + wt1*Hr)

    Vms = us*us + vs*vs + ws*ws
    as  = MixtPerf_C_GHoVm2(g,Hs,Vms)
    sl  = MIN(ql-al,qs-as)
    sr  = MAX(qr+ar,qs+as)

! ==============================================================================
!   Compute fluxes
! ==============================================================================

    IF ( sl > 0.0_RFREAL ) THEN 
      flx(1) =  ql* rl                    *nm
      flx(2) = (ql* rl*ul + pl*nx        )*nm
      flx(3) = (ql* rl*vl + pl*ny        )*nm
      flx(4) = (ql* rl*wl + pl*nz        )*nm
      flx(5) = (ql*(rl*el + pl)   + pl*fs)*nm
    ELSE IF ( sr < 0.0_RFREAL ) THEN  
      flx(1) =  qr* rr                    *nm
      flx(2) = (qr* rr*ur + pr*nx        )*nm
      flx(3) = (qr* rr*vr + pr*ny        )*nm
      flx(4) = (qr* rr*wr + pr*nz        )*nm
      flx(5) = (qr*(rr*er + pr)   + pr*fs)*nm    
    ELSE 
      sm = (rr*qr*(sr-qr) - rl*ql*(sl-ql) - pr + pl)/(rr*(sr-qr) - rl*(sl-ql))
      ps = pl + rl*(ql-sl)*(ql-sm) 
    
      IF ( sm > 0.0_RFREAL ) THEN 
        term = 1.0_RFREAL/(sl-sm)
        
        rs  = term* (sl-ql)*rl
        rus = term*((sl-ql)*rl*ul + (ps    - pl   )*nx)
        rvs = term*((sl-ql)*rl*vl + (ps    - pl   )*ny)
        rws = term*((sl-ql)*rl*wl + (ps    - pl   )*nz) 
        res = term*((sl-ql)*rl*el + (ps*sm - pl*ql)   )               
      ELSE 
        term = 1.0_RFREAL/(sr-sm)
        
        rs  = term* (sr-qr)*rr
        rus = term*((sr-qr)*rr*ur + (ps    - pr   )*nx)
        rvs = term*((sr-qr)*rr*vr + (ps    - pr   )*ny)
        rws = term*((sr-qr)*rr*wr + (ps    - pr   )*nz) 
        res = term*((sr-qr)*rr*er + (ps*sm - pr*qr)   )      
      END IF ! sm
      
      flx(1) =  sm* rs                  *nm
      flx(2) = (sm* rus + ps*nx        )*nm
      flx(3) = (sm* rvs + ps*ny        )*nm
      flx(4) = (sm* rws + ps*nz        )*nm 
      flx(5) = (sm*(res + ps)   - ps*fs)*nm           
    END IF ! sl

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

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + us*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + ws*pMf(indMf*ifg)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - us*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - ws*pMf(indMf*ifg)
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
  CALL FPROFILER_ENDS("RFLU::HLLC_ComputeFlux1_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_HLLC_ComputeFlux1_TCP








! ******************************************************************************
!
! Purpose: Compute central convective fluxes for multiphase mixture using 
!   second-order accurate version of HLLC scheme.
!
! Description: None.
!
! Input: 
!   pRegion     Data of current region.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_HLLC_ComputeFlux2_GL(pRegion)

  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtGasLiq_P, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtLiq_C2_Bp, &
                           MixtPerf_C2_GRT, &
                           MixtPerf_Cv_CpR, & 
                           MixtPerf_D_PRT, &
                           MixtPerf_R_CpG, &
                           MixtPerf_R_M, &     
                           RFLU_CentralFirstPatch_GL, &           
                           RFLU_CentralSecondPatch_GL                    
  
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

  INTEGER ::  c1,c2,ifg,indGs,indMf,indSd,iPatch
  REAL(RFREAL) :: Bgh2,Bgl2,Bgr2,Blh2,Bll2,Blr2,Bp,Bt,Bvh2,Bvl2,Bvr2,cml, &
                  cmr,cms,cvg,cvl,cvm,cvv,Cgh2,Cgl2,Cgr2,Clh2,Cll2,Clr2, &
                  Cvh2,Cvl2,Cvr2,dx,dy,dz,el,er,fs,irl,irr,nm,nx,ny,nz,ph, &
                  pl,po,pr,ps,ql,qr,qs,rel,rer,res,rgh,rgl,rgr,rh,rl,rlh, &
                  rll,rlr,ro,rr,rs,rus,rvh,rvl,rvr,rvs,rws,rYgl,rYgr,rYll, &
                  rYlr,rYvl,rYvr,Rg,Rv,sl,sm,sr,term,th,tl,to,tr,ul,ur,us, &
                  vfgh,vfgl,vfgr,vflh,vfll,vflr,vfvh,vfvl,vfvr,vl,vr,vs,Vl2, &
                  Vms,Vr2,wl,wr,ws,wt1,wt2,xc,yc,zc,Ygl,Ygr,Yll,Ylr,Yvl,Yvr
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pDvMixt,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradMixt,pGradSpec
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_patch), POINTER :: pPatch  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_HLLC_ComputeFlux2_GL',&
  'RFLU_ModHLLCFlux.F90')

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pGrid%indGs
  
  indMf  = pRegion%mixtInput%indMfMixt  
  indSd  = pRegion%mixtInput%indSd
  
  pCvMixt => pRegion%mixt%cv
  pDvMixt => pRegion%mixt%dv
  pCvSpec => pRegion%spec%cv
  
  pGradMixt => pRegion%mixt%gradCell
  pGradSpec => pRegion%spec%gradCell    
  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd 

  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel
  
! ******************************************************************************
! Define constants
! ******************************************************************************

  ro  = global%refDensityLiq
  po  = global%refPressLiq
  to  = global%refTempLiq 
  Bp  = global%refBetaPLiq
  Bt  = global%refBetaTLiq
  cvl = global%refCvLiq

  Rg  = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw) 
  cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht,Rg)

  Rv  = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
  cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht,Rv)

! *****************************************************************************
! Compute fluxes through interior faces
! *****************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! =============================================================================  
!   Get face geometry and grid speed
! =============================================================================  
      
    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)
  
    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)
  
    fs = pGrid%gs(indGs*ifg)  
  
! =============================================================================   
!   Compute left and right states     
! =============================================================================

! -----------------------------------------------------------------------------
!   Left state
! -----------------------------------------------------------------------------    
    
    rl  = pCvMixt(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    ul  = pCvMixt(CV_MIXT_XMOM,c1)*irl
    vl  = pCvMixt(CV_MIXT_YMOM,c1)*irl
    wl  = pCvMixt(CV_MIXT_ZMOM,c1)*irl
    tl  = pDvMixt(DV_MIXT_TEMP,c1)
    
    Ygl = pCvSpec(1,c1)*irl
    Yvl = pCvSpec(2,c1)*irl
         
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)    

    rl  = rl  + pGradMixt(XCOORD,GRC_MIXT_DENS,c1)*dx & 
              + pGradMixt(YCOORD,GRC_MIXT_DENS,c1)*dy & 
              + pGradMixt(ZCOORD,GRC_MIXT_DENS,c1)*dz 
    ul  = ul  + pGradMixt(XCOORD,GRC_MIXT_XVEL,c1)*dx & 
              + pGradMixt(YCOORD,GRC_MIXT_XVEL,c1)*dy & 
              + pGradMixt(ZCOORD,GRC_MIXT_XVEL,c1)*dz     
    vl  = vl  + pGradMixt(XCOORD,GRC_MIXT_YVEL,c1)*dx & 
              + pGradMixt(YCOORD,GRC_MIXT_YVEL,c1)*dy & 
              + pGradMixt(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl  = wl  + pGradMixt(XCOORD,GRC_MIXT_ZVEL,c1)*dx & 
              + pGradMixt(YCOORD,GRC_MIXT_ZVEL,c1)*dy & 
              + pGradMixt(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    tl  = tl  + pGradMixt(XCOORD,GRC_MIXT_TEMP,c1)*dx &
              + pGradMixt(YCOORD,GRC_MIXT_TEMP,c1)*dy &
              + pGradMixt(ZCOORD,GRC_MIXT_TEMP,c1)*dz    

    Ygl = Ygl + pGradSpec(XCOORD,1,c1)*dx &
              + pGradSpec(YCOORD,1,c1)*dy &
              + pGradSpec(ZCOORD,1,c1)*dz
    Yvl = Yvl + pGradSpec(XCOORD,2,c1)*dx &
              + pGradSpec(YCOORD,2,c1)*dy &
              + pGradSpec(ZCOORD,2,c1)*dz

    IF ( Ygl > 1.0_RFREAL ) THEN
      Ygl = 1.0_RFREAL
      Yvl = 0.0_RFREAL
    ELSE IF ( Ygl < 0.0_RFREAL ) THEN
      Ygl = 0.0_RFREAL
    ELSE IF ( Yvl > 1.0_RFREAL ) THEN
      Yvl = 1.0_RFREAL
      Ygl = 0.0_RFREAL
    ELSE IF ( Yvl < 0.0_RFREAL ) THEN
      Yvl = 0.0_RFREAL
    END IF ! Ygl

    Yll = 1.0_RFREAL - Ygl - Yvl

    IF ( Yll > 1.0_RFREAL) THEN
      Yll = 1.0_RFREAL
    ELSE IF ( Yll < 0.0_RFREAL) THEN
      Yll = 0.0_RFREAL
    END IF ! Yll
    
    Vl2 = ul*ul + vl*vl + wl*wl
    ql  = ul*nx + vl*ny + wl*nz - fs

    Cll2 = MixtLiq_C2_Bp(Bp)
    Cvl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tl)
    Cgl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tl)

    rYgl = rl*Ygl
    rYvl = rl*Yvl 
    rYll = rl*Yll       

    pl = MixtGasLiq_P(rYll,rYvl,rYgl,Cll2,Cvl2,Cgl2,rl,ro,po,to,Bp,Bt,tl)

    rll = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pl,po,tl,to)
    rvl = MixtPerf_D_PRT(pl,Rv,tl)
    rgl = MixtPerf_D_PRT(pl,Rg,tl)

    vfgl = rYgl/rgl
    vfvl = rYvl/rvl
    vfll = rYll/rll      
   
    cvm = (rYll*cvl + rYvl*cvv + rYgl*cvg)/rl   
    rel  = rl*cvm*tl + 0.5*rl*Vl2 
    el   = rel/rl
      
    Bll2 = -Bt/Bp
    Bvl2 = rvl*Rv
    Bgl2 = rgl*Rg
    
    cml = MixtGasLiq_C(cvm,rl,pl,rll,rvl,rgl,vfll,vfvl,vfgl,Cll2,Cvl2,Cgl2, &
                       Bll2,Bvl2,Bgl2)

! -----------------------------------------------------------------------------
!   Right state
! -----------------------------------------------------------------------------    
    
    rr  = pCvMixt(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    ur  = pCvMixt(CV_MIXT_XMOM,c2)*irr
    vr  = pCvMixt(CV_MIXT_YMOM,c2)*irr
    wr  = pCvMixt(CV_MIXT_ZMOM,c2)*irr
    tr  = pDvMixt(DV_MIXT_TEMP,c2)  
 
    Ygr = pCvSpec(1,c2)*irr
    Yvr = pCvSpec(2,c2)*irr 
       
    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)    

    rr  = rr  + pGradMixt(XCOORD,GRC_MIXT_DENS,c2)*dx &
              + pGradMixt(YCOORD,GRC_MIXT_DENS,c2)*dy &
              + pGradMixt(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur  = ur  + pGradMixt(XCOORD,GRC_MIXT_XVEL,c2)*dx & 
              + pGradMixt(YCOORD,GRC_MIXT_XVEL,c2)*dy & 
              + pGradMixt(ZCOORD,GRC_MIXT_XVEL,c2)*dz     
    vr  = vr  + pGradMixt(XCOORD,GRC_MIXT_YVEL,c2)*dx & 
              + pGradMixt(YCOORD,GRC_MIXT_YVEL,c2)*dy & 
              + pGradMixt(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr  = wr  + pGradMixt(XCOORD,GRC_MIXT_ZVEL,c2)*dx & 
              + pGradMixt(YCOORD,GRC_MIXT_ZVEL,c2)*dy & 
              + pGradMixt(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    tr  = tr  + pGradMixt(XCOORD,GRC_MIXT_TEMP,c2)*dx &
              + pGradMixt(YCOORD,GRC_MIXT_TEMP,c2)*dy &
              + pGradMixt(ZCOORD,GRC_MIXT_TEMP,c2)*dz   
 
    Ygr = Ygr + pGradSpec(XCOORD,1,c2)*dx &
              + pGradSpec(YCOORD,1,c2)*dy &
              + pGradSpec(ZCOORD,1,c2)*dz
    Yvr = Yvr + pGradSpec(XCOORD,2,c2)*dx &
              + pGradSpec(YCOORD,2,c2)*dy &
              + pGradSpec(ZCOORD,2,c2)*dz

    IF ( Ygr > 1.0_RFREAL ) THEN
      Ygr = 1.0_RFREAL
      Yvr = 0.0_RFREAL
    ELSE IF ( Ygr < 0.0_RFREAL ) THEN
      Ygr = 0.0_RFREAL
    ELSE IF ( Yvr > 1.0_RFREAL ) THEN
      Yvr = 1.0_RFREAL
      Ygr = 0.0_RFREAL
    ELSE IF ( Yvr < 0.0_RFREAL ) THEN
      Yvr = 0.0_RFREAL
    END IF ! Ygr
 
    Ylr = 1.0_RFREAL - Ygr - Yvr

    IF ( Ylr > 1.0_RFREAL ) THEN
      Ylr = 1.0_RFREAL
    ELSE IF ( Ylr < 0.0_RFREAL ) THEN 
      Ylr = 0.0_RFREAL
    END IF ! Ylr

    Vr2 = ur*ur + vr*vr + wr*wr
    qr  = ur*nx + vr*ny + wr*nz - fs

    Clr2 = MixtLiq_C2_Bp(Bp)
    Cvr2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tr)
    Cgr2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tr)

    rYgr = rr*Ygr
    rYvr = rr*Yvr
    rYlr = rr*Ylr      

    pr = MixtGasLiq_P(rYlr,rYvr,rYgr,Clr2,Cvr2,Cgr2,rr,ro,po,to,Bp,Bt,tr)

    rlr = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pr,po,tr,to)
    rvr = MixtPerf_D_PRT(pr,Rv,tr)
    rgr = MixtPerf_D_PRT(pr,Rg,tr)

    vfgr = rYgr/rgr
    vfvr = rYvr/rvr
    vflr = rYlr/rlr    

    cvm = (rYlr*cvl + rYvr*cvv + rYgr*cvg)/rr
    rer  = rr*cvm*tr + 0.5*rr*Vr2
    er   = rer/rr
 
    Blr2 = -Bt/Bp
    Bvr2 = rvr*Rv
    Bgr2 = rgr*Rg
    
    cmr = MixtGasLiq_C(cvm,rr,pr,rlr,rvr,rgr,vflr,vfvr,vfgr,Clr2,Cvr2,Cgr2, &
                       Blr2,Bvr2,Bgr2)

! ==============================================================================
!   Compute weights
! ==============================================================================

    wt1 = SQRT(rr/rl)
    wt2 = 1.0_RFREAL/(1.0_RFREAL + wt1)   
   
! ==============================================================================
!   Compute averages
! ==============================================================================

    us  = wt2*(ul + wt1*ur)
    vs  = wt2*(vl + wt1*vr)
    ws  = wt2*(wl + wt1*wr)
    qs  = wt2*(ql + wt1*qr)

    Vms = us*us + vs*vs + ws*ws
   
    rh  = wt1*rl
    ph  = wt2*(pl + wt1*pr)  
    th  = wt2*(tl + wt1*tr)  
    
    vfgh = wt2*(vfgl + wt1*vfgr)
    vfvh = wt2*(vfvl + wt1*vfvr)
 
    IF ( vfgh > 1.0_RFREAL ) THEN
      vfgh = 1.0_RFREAL
      vfvh = 0.0_RFREAL
    ELSE IF ( vfgh < 0.0_RFREAL ) THEN
      vfgh = 0.0_RFREAL
    ELSE IF ( vfvh > 1.0_RFREAL ) THEN
      vfvh = 1.0_RFREAL
      vfgh = 0.0_RFREAL
    ELSE IF ( vfvh < 0.0_RFREAL ) THEN
      vfvh = 0.0_RFREAL
    END IF ! vgfh

    vflh = 1.0_RFREAL - vfvh - vfgh

    IF ( vflh > 1.0_RFREAL ) THEN
      vflh = 1.0_RFREAL
    ELSE IF ( vflh < 0.0_RFREAL ) THEN 
      vflh = 0.0_RFREAL
    END IF ! vflh
 
    Clh2 = MixtLiq_C2_Bp(Bp) 
    Cvh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,th)
    Cgh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,th)

    rlh = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,ph,po,th,to)
    rvh = MixtPerf_D_PRT(ph,Rv,th)
    rgh = MixtPerf_D_PRT(ph,Rg,th)

    Blh2 = -Bt/Bp
    Bvh2 = rvh*Rv
    Bgh2 = rgh*Rg

    cvm = (rlh*vflh*cvl + rvh*vfvh*cvv + rgh*vfgh*cvg)/rh
 
    cms = MixtGasLiq_C(cvm,rh,ph,rlh,rvh,rgh,vflh,vfvh,vfgh,Clh2,Cvh2,Cgh2, &
                       Blh2,Bvh2,Bgh2)

    sl = MIN(ql-cml,qs-cms)
    sr = MAX(qr+cmr,qs+cms) 
   
! ==============================================================================
!   Compute fluxes
! ==============================================================================
      
    IF ( sl > 0.0_RFREAL ) THEN 
      flx(1) =  ql* rl                    *nm
      flx(2) = (ql* rl*ul + pl*nx        )*nm
      flx(3) = (ql* rl*vl + pl*ny        )*nm
      flx(4) = (ql* rl*wl + pl*nz        )*nm
      flx(5) = (ql*(rl*el + pl)   + pl*fs)*nm
    ELSE IF ( sr < 0.0_RFREAL ) THEN  
      flx(1) =  qr* rr                    *nm
      flx(2) = (qr* rr*ur + pr*nx        )*nm
      flx(3) = (qr* rr*vr + pr*ny        )*nm
      flx(4) = (qr* rr*wr + pr*nz        )*nm
      flx(5) = (qr*(rr*er + pr)   + pr*fs)*nm 
    ELSE 
      sm = (rr*qr*(sr-qr) - rl*ql*(sl-ql) - pr + pl)/(rr*(sr-qr) - rl*(sl-ql))
      ps = pl + rl*(ql-sl)*(ql-sm) 
    
      IF ( sm > 0.0_RFREAL ) THEN 
        term = 1.0_RFREAL/(sl-sm)
        
        rs    = term* (sl-ql)*rl
        rus   = term*((sl-ql)*rl*ul + (ps    - pl   )*nx)
        rvs   = term*((sl-ql)*rl*vl + (ps    - pl   )*ny)
        rws   = term*((sl-ql)*rl*wl + (ps    - pl   )*nz) 
        res   = term*((sl-ql)*rl*el + (ps*sm - pl*ql)   )
      ELSE 
        term = 1.0_RFREAL/(sr-sm)
        
        rs    = term* (sr-qr)*rr
        rus   = term*((sr-qr)*rr*ur + (ps    - pr   )*nx)
        rvs   = term*((sr-qr)*rr*vr + (ps    - pr   )*ny)
        rws   = term*((sr-qr)*rr*wr + (ps    - pr   )*nz) 
        res   = term*((sr-qr)*rr*er + (ps*sm - pr*qr)   ) 
      END IF ! sm
            
      flx(1) =  sm* rs                  *nm
      flx(2) = (sm* rus + ps*nx        )*nm
      flx(3) = (sm* rvs + ps*ny        )*nm
      flx(4) = (sm* rws + ps*nz        )*nm 
      flx(5) = (sm*(res + ps)   - ps*fs)*nm  
    END IF ! sl

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

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + us*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + ws*pMf(indMf*ifg)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - us*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - ws*pMf(indMf*ifg)
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************
 
  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch_GL(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch_GL(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************    

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_HLLC_ComputeFlux2_GL
  
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate version of
!   HLLC scheme for mixture of thermally and calorically perfect gases.
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

SUBROUTINE RFLU_HLLC_ComputeFlux2_MTCP(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtPerf_C_DGP, &
                           MixtPerf_C_GHoVm2, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_R_M, & 
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
  REAL(RFREAL) :: al,ar,as,dx1,dx2,dy1,dy2,dz1,dz2,el,er,fs,gcl,gcr,gl,gr,gs, &
                  Hl,Hr,Hs,irl,irr,nm,nx,ny,nz,ql,qr,qs,pl,pr,ps,rl,rr,res, &
                  rs,rus,rvs,rws,sl,sm,sr,term,ul,ur,us,vl,Vml,Vmr,Vms,vr,vs, &
                  wl,wr,ws,wt1,wt2,xc,yc,zc
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef SPEC
  INTEGER :: iSpec
  REAL(RFREAL) :: cpl,cpr,cps,gcs,mml,mmr,mms,Y1,Y2,Ys
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec 
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGcSpec 
  TYPE(t_spec_type), POINTER :: pSpecType
#endif  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_HLLC_ComputeFlux2_MTCP',&
  'RFLU_ModHLLCFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::HLLC_ComputeFlux2_MTCP")
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

    dx1 = xc - pGrid%cofg(XCOORD,c1)
    dy1 = yc - pGrid%cofg(YCOORD,c1)
    dz1 = zc - pGrid%cofg(ZCOORD,c1)

#ifdef SPEC
    mml = 0.0_RFREAL
    cpl = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y1 = irl*pCvSpec(iSpec,c1) + pGcSpec(XCOORD,iSpec,c1)*dx1 & 
                                 + pGcSpec(YCOORD,iSpec,c1)*dy1 & 
                                 + pGcSpec(ZCOORD,iSpec,c1)*dz1

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

    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx1 &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy1 &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz1
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx1 &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy1 &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz1
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx1 &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy1 &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz1
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx1 &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy1 &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz1
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx1 &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy1 &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz1

    el  = MixtPerf_Eo_DGPUVW(rl,gl,pl,ul,vl,wl)
    al  = MixtPerf_C_DGP(rl,gl,pl)

    Vml = ul*ul + vl*vl + wl*wl
    ql  = ul*nx + vl*ny + wl*nz - fs
    Hl  = el + irl*pl

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    dx2 = xc - pGrid%cofg(XCOORD,c2)
    dy2 = yc - pGrid%cofg(YCOORD,c2)
    dz2 = zc - pGrid%cofg(ZCOORD,c2)
    
#ifdef SPEC
    mmr = 0.0_RFREAL
    cpr = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y2 = irr*pCvSpec(iSpec,c2) + pGcSpec(XCOORD,iSpec,c2)*dx2 & 
                                 + pGcSpec(YCOORD,iSpec,c2)*dy2 & 
                                 + pGcSpec(ZCOORD,iSpec,c2)*dz2

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

    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx2 &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy2 &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz2
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx2 &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy2 &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz2
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx2 &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy2 &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz2
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx2 &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy2 &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz2
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx2 &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy2 &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz2

    er  = MixtPerf_Eo_DGPUVW(rr,gr,pr,ur,vr,wr)
    ar  = MixtPerf_C_DGP(rr,gr,pr)

    Vmr = ur*ur + vr*vr + wr*wr
    qr  = ur*nx + vr*ny + wr*nz - fs
    Hr  = er + irr*pr
    
! ==============================================================================
!   Compute weights
! ==============================================================================

    wt1 = SQRT(rr/rl)
    wt2 = 1.0_RFREAL/(1.0_RFREAL + wt1)

! ==============================================================================
!   Compute averages. NOTE this average differs from that of the TCP case in 
!   that need to recompute gas properties to compute speed of sound. Choose to
!   apply Roe-average also to mass fractions, although there is no analysis to
!   back this up. It is not clear what is the best way to do this.
! ==============================================================================

    us = wt2*(ul + wt1*ur)
    vs = wt2*(vl + wt1*vr)
    ws = wt2*(wl + wt1*wr)
    qs = wt2*(ql + wt1*qr)
    Hs = wt2*(Hl + wt1*Hr)    

#ifdef SPEC
    mms = 0.0_RFREAL
    cps = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y1 = irl*pCvSpec(iSpec,c1) + pGcSpec(XCOORD,iSpec,c1)*dx1 & 
                                 + pGcSpec(YCOORD,iSpec,c1)*dy1 & 
                                 + pGcSpec(ZCOORD,iSpec,c1)*dz1
                                 
      Y2 = irr*pCvSpec(iSpec,c2) + pGcSpec(XCOORD,iSpec,c2)*dx2 & 
                                 + pGcSpec(YCOORD,iSpec,c2)*dy2 & 
                                 + pGcSpec(ZCOORD,iSpec,c2)*dz2

      Ys = wt2*(Y1 + wt1*Y2)

      mms = mms + Ys/pSpecType%pMaterial%molw
      cps = cps + Ys*pSpecType%pMaterial%spht
    END DO ! iSpec
    
    mms = 1.0_RFREAL/mms 
    gcs = MixtPerf_R_M(mms)
    gs  = MixtPerf_G_CpR(cps,gcs)              
#endif

    Vms = us*us + vs*vs + ws*ws
    as  = MixtPerf_C_GHoVm2(gs,Hs,Vms)    
    sl  = MIN(ql-al,qs-as)
    sr  = MAX(qr+ar,qs+as)

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    IF ( sl > 0.0_RFREAL ) THEN 
      flx(1) =  ql* rl                    *nm
      flx(2) = (ql* rl*ul + pl*nx        )*nm
      flx(3) = (ql* rl*vl + pl*ny        )*nm
      flx(4) = (ql* rl*wl + pl*nz        )*nm
      flx(5) = (ql*(rl*el + pl)   + pl*fs)*nm
    ELSE IF ( sr < 0.0_RFREAL ) THEN  
      flx(1) =  qr* rr                    *nm
      flx(2) = (qr* rr*ur + pr*nx        )*nm
      flx(3) = (qr* rr*vr + pr*ny        )*nm
      flx(4) = (qr* rr*wr + pr*nz        )*nm
      flx(5) = (qr*(rr*er + pr)   + pr*fs)*nm    
    ELSE 
      sm = (rr*qr*(sr-qr) - rl*ql*(sl-ql) - pr + pl)/(rr*(sr-qr) - rl*(sl-ql))
      ps = pl + rl*(ql-sl)*(ql-sm) 
    
      IF ( sm > 0.0_RFREAL ) THEN 
        term = 1.0_RFREAL/(sl-sm)
        
        rs  = term* (sl-ql)*rl
        rus = term*((sl-ql)*rl*ul + (ps    - pl   )*nx)
        rvs = term*((sl-ql)*rl*vl + (ps    - pl   )*ny)
        rws = term*((sl-ql)*rl*wl + (ps    - pl   )*nz) 
        res = term*((sl-ql)*rl*el + (ps*sm - pl*ql)   )               
      ELSE 
        term = 1.0_RFREAL/(sr-sm)
        
        rs  = term* (sr-qr)*rr
        rus = term*((sr-qr)*rr*ur + (ps    - pr   )*nx)
        rvs = term*((sr-qr)*rr*vr + (ps    - pr   )*ny)
        rws = term*((sr-qr)*rr*wr + (ps    - pr   )*nz) 
        res = term*((sr-qr)*rr*er + (ps*sm - pr*qr)   )      
      END IF ! sm
      
      flx(1) =  sm* rs                  *nm
      flx(2) = (sm* rus + ps*nx        )*nm
      flx(3) = (sm* rvs + ps*ny        )*nm
      flx(4) = (sm* rws + ps*nz        )*nm 
      flx(5) = (sm*(res + ps)   - ps*fs)*nm           
    END IF ! sl

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

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + us*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + ws*pMf(indMf*ifg)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - us*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - ws*pMf(indMf*ifg)
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
  CALL FPROFILER_ENDS("RFLU::HLLC_ComputeFlux2_MTCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_HLLC_ComputeFlux2_MTCP


  
  
  
  
  

! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate version of
!   HLLC scheme for thermally and calorically perfect gas.
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

SUBROUTINE RFLU_HLLC_ComputeFlux2_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_DGP, &
                           MixtPerf_C_GHoVm2, &
                           MixtPerf_Eo_DGPUVW, &
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
  REAL(RFREAL) :: al,ar,as,cp,dx,dy,dz,el,er,fs,g,Hl,Hr,Hs,irl,irr,nm, &
                  nx,ny,nz,ql,qr,qs,pl,pr,ps,rl,rr,res,rs,rus,rvs,rws, &
                  sl,sm,sr,term,ul,ur,us,vl,Vml,Vmr,Vms,vr,vs,wl,wr,ws, &
                  wt1,wt2,xc,yc,zc
  REAL(RFREAL) :: flx(5)
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

  CALL RegisterFunction(global,'RFLU_HLLC_ComputeFlux2_TCP',&
  'RFLU_ModHLLCFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::HLLC_ComputeFlux2_TCP")
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

    el  = MixtPerf_Eo_DGPUVW(rl,g,pl,ul,vl,wl)
    al  = MixtPerf_C_DGP(rl,g,pl)

    Vml = ul*ul + vl*vl + wl*wl
    ql  = ul*nx + vl*ny + wl*nz - fs
    Hl  = el + irl*pl

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

    er  = MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
    ar  = MixtPerf_C_DGP(rr,g,pr)

    Vmr = ur*ur + vr*vr + wr*wr
    qr  = ur*nx + vr*ny + wr*nz - fs
    Hr  = er + irr*pr

! ==============================================================================
!   Compute weights
! ==============================================================================

    wt1 = SQRT(rr/rl)
    wt2 = 1.0_RFREAL/(1.0_RFREAL + wt1)

! ==============================================================================
!   Compute averages
! ==============================================================================

    us  = wt2*(ul + wt1*ur)
    vs  = wt2*(vl + wt1*vr)
    ws  = wt2*(wl + wt1*wr)
    qs  = wt2*(ql + wt1*qr)
    Hs  = wt2*(Hl + wt1*Hr)

    Vms = us*us + vs*vs + ws*ws
    as  = MixtPerf_C_GHoVm2(g,Hs,Vms)
    sl  = MIN(ql-al,qs-as)
    sr  = MAX(qr+ar,qs+as)

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    IF ( sl > 0.0_RFREAL ) THEN 
      flx(1) =  ql* rl                    *nm
      flx(2) = (ql* rl*ul + pl*nx        )*nm
      flx(3) = (ql* rl*vl + pl*ny        )*nm
      flx(4) = (ql* rl*wl + pl*nz        )*nm
      flx(5) = (ql*(rl*el + pl)   + pl*fs)*nm
    ELSE IF ( sr < 0.0_RFREAL ) THEN  
      flx(1) =  qr* rr                    *nm
      flx(2) = (qr* rr*ur + pr*nx        )*nm
      flx(3) = (qr* rr*vr + pr*ny        )*nm
      flx(4) = (qr* rr*wr + pr*nz        )*nm
      flx(5) = (qr*(rr*er + pr)   + pr*fs)*nm    
    ELSE 
      sm = (rr*qr*(sr-qr) - rl*ql*(sl-ql) - pr + pl)/(rr*(sr-qr) - rl*(sl-ql))
      ps = pl + rl*(ql-sl)*(ql-sm) 
    
      IF ( sm > 0.0_RFREAL ) THEN 
        term = 1.0_RFREAL/(sl-sm)
        
        rs  = term* (sl-ql)*rl
        rus = term*((sl-ql)*rl*ul + (ps    - pl   )*nx)
        rvs = term*((sl-ql)*rl*vl + (ps    - pl   )*ny)
        rws = term*((sl-ql)*rl*wl + (ps    - pl   )*nz) 
        res = term*((sl-ql)*rl*el + (ps*sm - pl*ql)   )               
      ELSE 
        term = 1.0_RFREAL/(sr-sm)
        
        rs  = term* (sr-qr)*rr
        rus = term*((sr-qr)*rr*ur + (ps    - pr   )*nx)
        rvs = term*((sr-qr)*rr*vr + (ps    - pr   )*ny)
        rws = term*((sr-qr)*rr*wr + (ps    - pr   )*nz) 
        res = term*((sr-qr)*rr*er + (ps*sm - pr*qr)   )      
      END IF ! sm
      
      flx(1) =  sm* rs                  *nm
      flx(2) = (sm* rus + ps*nx        )*nm
      flx(3) = (sm* rvs + ps*ny        )*nm
      flx(4) = (sm* rws + ps*nz        )*nm 
      flx(5) = (sm*(res + ps)   - ps*fs)*nm           
    END IF ! sl

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

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + us*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + ws*pMf(indMf*ifg)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - us*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vs*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - ws*pMf(indMf*ifg)
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
  CALL FPROFILER_ENDS("RFLU::HLLC_ComputeFlux2_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_HLLC_ComputeFlux2_TCP

  
  
  
  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModHLLCFlux

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModHLLCFlux.F90,v $
! Revision 1.10  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/05/01 22:15:29  haselbac
! Bug fix: Moved gs out of SPEC declaration, removed debug statements
!
! Revision 1.7  2006/05/01 21:00:47  haselbac
! Rewrite for consistency and cleanliness
!
! Revision 1.6  2006/04/15 16:58:42  haselbac
! Added capability of running 1st order boundary fluxes with 2nd order volume 
! fluxes
!
! Revision 1.5  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.4  2006/03/30 20:49:40  haselbac
! Bug fixes for perf gas routines, cosmetics elsewhere
!
! Revision 1.3  2006/03/26 20:22:02  haselbac
! Added support for GL model
!
! Revision 1.2  2005/07/07 22:45:01  haselbac
! Added profiling calls
!
! Revision 1.1  2005/05/16 20:36:29  haselbac
! Initial revision
! ******************************************************************************
  
  













