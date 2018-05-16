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
! Purpose: Collection of Roe flux routines.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModRoeFlux.F90,v 1.10 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRoeFlux

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE RFLU_ModEntropyFixes

  USE RFLU_ModNSCBC

  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: RFLU_ROE_ComputeFlux, &
            RFLU_ROE_ComputeFluxC1, & 
            RFLU_ROE_ComputeFluxC2_TCP, & 
            RFLU_ROE_ComputeFluxD1_TCP, & 
            RFLU_ROE_ComputeFluxD2_TCP, &         
            RFLU_ROE_ComputeFlux1_GL, & 
            RFLU_ROE_ComputeFlux1_TCP, &             
            RFLU_ROE_ComputeFlux2_GL, &
            RFLU_ROE_ComputeFlux2_TCP
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModRoeFlux.F90,v $ $Revision: 1.10 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS








! ******************************************************************************
!
! Purpose: Wrapper function for Roe flux functions.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!   fluxPart	Part of flux (central or dissipation)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ROE_ComputeFlux(pRegion,fluxPart)
                       
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: fluxPart
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================
  
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFlux',&
  'RFLU_ModRoeFlux.F90')

! ******************************************************************************
! Call flux functions
! ******************************************************************************

  SELECT CASE ( global%flowType ) 

! ==============================================================================
!   Unsteady flow
! ==============================================================================

    CASE ( FLOW_UNSTEADY )   
      SELECT CASE ( pRegion%mixtInput%gasModel ) 

! ------------------------------------------------------------------------------    
!       Thermally and calorically perfect gas         
! ------------------------------------------------------------------------------    

        CASE ( GAS_MODEL_TCPERF )
          SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
            CASE ( DISCR_ORDER_1 ) 
              CALL RFLU_ROE_ComputeFlux1_TCP(pRegion)
            CASE ( DISCR_ORDER_2 )
              CALL RFLU_ROE_ComputeFlux2_TCP(pRegion)   
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%spaceOrder          

! ------------------------------------------------------------------------------    
!       Mixture of gas, vapor, and liquid     
! ------------------------------------------------------------------------------    

        CASE ( GAS_MODEL_MIXT_GASLIQ )
          SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
            CASE ( DISCR_ORDER_1 ) 
              CALL RFLU_ROE_ComputeFlux1_GL(pRegion)
            CASE ( DISCR_ORDER_2 )
              CALL RFLU_ROE_ComputeFlux2_GL(pRegion)               
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%spaceOrder           

! ------------------------------------------------------------------------------    
!       Default           
! ------------------------------------------------------------------------------    

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)             
      END SELECT ! pRegion%mixtInput%gasModel  

! ==============================================================================
!   Steady flow. NOTE for explicit solver, distinguish between central and
!   dissipative parts of flux.
! ==============================================================================

    CASE ( FLOW_STEADY )
      SELECT CASE ( global%solverType ) 

! ------------------------------------------------------------------------------    
!       Explicit solver
! ------------------------------------------------------------------------------    

        CASE ( SOLV_EXPLICIT ) 
          SELECT CASE ( fluxPart ) 

! --------- Central part of flux -----------------------------------------------

            CASE ( FLUX_PART_CENTRAL ) 
              SELECT CASE ( pRegion%mixtInput%gasModel ) 
                CASE ( GAS_MODEL_TCPERF )
                  SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
                    CASE ( DISCR_ORDER_1 ) 
                      CALL RFLU_ROE_ComputeFluxC1(pRegion)
                    CASE ( DISCR_ORDER_2 )
                      CALL RFLU_ROE_ComputeFluxC2_TCP(pRegion)              
                    CASE DEFAULT
                      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                  END SELECT ! pRegion%mixtInput%spaceOrder          
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)             
              END SELECT ! pRegion%mixtInput%gasModel              

! --------- Dissipative part of flux -------------------------------------------

            CASE ( FLUX_PART_DISSIP ) 
              SELECT CASE ( pRegion%mixtInput%gasModel ) 
                CASE ( GAS_MODEL_TCPERF )
                  SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
                    CASE ( DISCR_ORDER_1 ) 
                      CALL RFLU_ROE_ComputeFluxD1_TCP(pRegion)
                    CASE ( DISCR_ORDER_2 )
                      CALL RFLU_ROE_ComputeFluxD2_TCP(pRegion)           
                    CASE DEFAULT
                      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                  END SELECT ! pRegion%mixtInput%spaceOrder          
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)             
              END SELECT ! pRegion%mixtInput%gasModel            

! --------- Default ------------------------------------------------------------

            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
          END SELECT ! fluxPart

! ------------------------------------------------------------------------------    
!       Implicit solver. NOTE compute entire flux.             
! ------------------------------------------------------------------------------    

        CASE ( SOLV_IMPLICIT_NK ) 
          SELECT CASE ( pRegion%mixtInput%gasModel ) 
            CASE ( GAS_MODEL_TCPERF )
              SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
                CASE ( DISCR_ORDER_1 ) 
                  CALL RFLU_ROE_ComputeFlux1_TCP(pRegion)
                CASE ( DISCR_ORDER_2 )
                  CALL RFLU_ROE_ComputeFlux2_TCP(pRegion)               
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
              END SELECT ! pRegion%mixtInput%spaceOrder          
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)             
          END SELECT ! pRegion%mixtInput%gasModel              

! ------------------------------------------------------------------------------    
!       Default 
! ------------------------------------------------------------------------------    

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)   
      END SELECT ! global%solverType                               
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
  END SELECT ! global%flowType      

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFlux








! ******************************************************************************
!
! Purpose: Compute central convective fluxes using first-order accurate
!   accurate approximation by average of variables.
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

SUBROUTINE RFLU_ROE_ComputeFluxC1(pRegion)

  USE ModInterfaces, ONLY: RFLU_CentralFirstPatch

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

  INTEGER :: c1,c2,ifg,indGs,iPatch
  REAL(RFREAL) :: el,er,fs,irl,irr,nm,nx,ny,nz,ql,qr,rl,rr,ul,ur,vl,vr,wl,wr, &
                  pl,pr
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pRhs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFluxC1',&
  'RFLU_ModRoeFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::ROE_ComputeFluxC1")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs = pGrid%indGs

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  
  pRhs => pRegion%mixt%rhs

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

    ql  = ul*nx + vl*ny + wl*nz - fs

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

    qr  = ur*nx + vr*ny + wr*nz - fs

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    flx(1) = 0.5_RFREAL*(ql* rl                  + qr* rr                 )*nm
    flx(2) = 0.5_RFREAL*(ql* rl*ul       + pl*nx + qr* rr*ur       + pr*nx)*nm
    flx(3) = 0.5_RFREAL*(ql* rl*vl       + pl*ny + qr* rr*vr       + pr*ny)*nm
    flx(4) = 0.5_RFREAL*(ql* rl*wl       + pl*nz + qr* rr*wr       + pr*nz)*nm
    flx(5) = 0.5_RFREAL*(ql*(rl*el + pl) + pl*fs + qr*(rr*er + pr) + pr*fs)*nm

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
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      CALL RFLU_NSCBC_CompFirstPatchFlux(pRegion,pPatch)
    ELSE
      CALL RFLU_CentralFirstPatch(pRegion,pPatch)
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::ROE_ComputeFluxC1")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFluxC1










! ******************************************************************************
!
! Purpose: Compute central convective fluxes using second-order accurate
!   accurate approximation by average of variables for thermally and 
!   calorically perfect gas.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assumes ideal gas and constant Cp and R! This assumption is necessary
!      because need to compute face states from reconstructed variables.
!
! ******************************************************************************

SUBROUTINE RFLU_ROE_ComputeFluxC2_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_Eo_DGPUVW, & 
                           RFLU_CentralFirstPatch, &
                           RFLU_CentralSecondPatch

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indGs,iPatch
  REAL(RFREAL) :: dx,dy,dz,el,er,fs,g,irl,irr,nm,nx,ny,nz,ql,qr,rl,rr,ul,ur, &
                  vl,vr,wl,wr,pl,pr,xc,yc,zc
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFluxC2_TCP',&
  'RFLU_ModRoeFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::ROE_ComputeFluxC2_TCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs = pGrid%indGs

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  
  pGc  => pRegion%mixt%gradCell    
  pRhs => pRegion%mixt%rhs  

! ******************************************************************************
! Define constants
! ******************************************************************************

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

    el = MixtPerf_Eo_DGPUVW(rl,g,pl,ul,vl,wl)
    ql = ul*nx + vl*ny + wl*nz - fs

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

    er = MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
    qr = ur*nx + vr*ny + wr*nz - fs

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    flx(1) = 0.5_RFREAL*(ql* rl                  + qr* rr                 )*nm
    flx(2) = 0.5_RFREAL*(ql* rl*ul       + pl*nx + qr* rr*ur       + pr*nx)*nm
    flx(3) = 0.5_RFREAL*(ql* rl*vl       + pl*ny + qr* rr*vr       + pr*ny)*nm
    flx(4) = 0.5_RFREAL*(ql* rl*wl       + pl*nz + qr* rr*wr       + pr*nz)*nm
    flx(5) = 0.5_RFREAL*(ql*(rl*el + pl) + pl*fs + qr*(rr*er + pr) + pr*fs)*nm

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
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      SELECT CASE ( pPatch%spaceOrder )
        CASE ( 1 ) 
          CALL RFLU_NSCBC_CompFirstPatchFlux(pRegion,pPatch)
        CASE ( 2 ) 
          CALL RFLU_NSCBC_CompSecondPatchFlux(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%spaceOrder
    ELSE
      SELECT CASE ( pPatch%spaceOrder )
        CASE ( 1 ) 
          CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
        CASE ( 2 ) 
          CALL RFLU_CentralSecondPatch(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%spaceOrder
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::ROE_ComputeFluxC2_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFluxC2_TCP







! ******************************************************************************
!
! Purpose: Compute dissipative fluxes using first-order accurate Roe scheme
!   for thermally and calorically perfect gas.
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

SUBROUTINE RFLU_ROE_ComputeFluxD1_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GHoVm2, &
                           MixtPerf_Ho_CpTUVW

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

  INTEGER :: c1,c2,ifg,indGs
  REAL(RFREAL) :: ah,betrk,cp,dissFact,dp,dq,dr,du,de,dv1,dv2,dv5,dw,efc, &
                  epsentr,fs,g,Hh,Hl,Hr,irl,irr,l1,l2,l5,nm,nx,ny,nz,pl, &
                  pr,qh,ql,qr,rl,rh,rr,sh,term,tl,tr,t1,t2,t3,t5,ul,uh,ur, &
                  vl,vh,vr,wl,wh,wr,wt
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDiss,pDv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid    

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFluxD1_TCP',&
  'RFLU_ModRoeFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::ROE_ComputeFluxD1_TCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs = pGrid%indGs

  pCv   => pRegion%mixt%cv
  pDv   => pRegion%mixt%dv
  
  pDiss => pRegion%mixt%diss

! ******************************************************************************
! Define constants
! ******************************************************************************

  cp = global%refCp
  g  = global%refGamma

  epsentr  = pRegion%mixtInput%epsentr
  betrk    = pRegion%mixtInput%betrk(pRegion%irkStep)
  dissFact = pRegion%mixtInput%dissFact
  term     = 0.5_RFREAL*betrk*dissFact

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
    pl  = pDv(DV_MIXT_PRES,c1)
    tl  = pDv(DV_MIXT_TEMP,c1)

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
    tr  = pDv(DV_MIXT_TEMP,c2)

    Hr = MixtPerf_Ho_CpTUVW(cp,tr,ur,vr,wr)

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)
    wt = rl/(rl + rh)

    uh = wt*ul + (1.0_RFREAL-wt)*ur
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    Hh = wt*Hl + (1.0_RFREAL-wt)*Hr

    qh = uh*nx + vh*ny + wh*nz - fs
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    ah = MixtPerf_C_GHoVm2(g,Hh,2.0_RFREAL*sh) ! NOTE factor of 2

! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh - ah)
    l2 = ABS(qh     )
    l5 = ABS(qh + ah)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr = rr - rl
    du = ur - ul
    de = vr - vl ! use de to denote dv because dv is already used
    dw = wr - wl
    dp = pr - pl
    dq = du*nx + de*ny + dw*nz ! Note: no gs contribution

    dv1 = (dp - rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv5 = (dp + rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv2 = dr - dp/(ah*ah)

    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l2*rh
    t5 = l5*dv5

! ==============================================================================
!   Compute fluxes
! ==============================================================================

    flx(1) = term*(t1            + t2                                   + t5           )*nm
    flx(2) = term*(t1*(uh-nx*ah) + t2*uh + t3*(du-nx*dq)                + t5*(uh+nx*ah))*nm
    flx(3) = term*(t1*(vh-ny*ah) + t2*vh + t3*(de-ny*dq)                + t5*(vh+ny*ah))*nm
    flx(4) = term*(t1*(wh-nz*ah) + t2*wh + t3*(dw-nz*dq)                + t5*(wh+nz*ah))*nm
    flx(5) = term*(t1*(Hh-qh*ah) + t2*sh + t3*(uh*du+vh*de+wh*dw-qh*dq) + t5*(Hh+qh*ah))*nm

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pDiss(CV_MIXT_DENS,c1) = pDiss(CV_MIXT_DENS,c1) + flx(1)
    pDiss(CV_MIXT_XMOM,c1) = pDiss(CV_MIXT_XMOM,c1) + flx(2)
    pDiss(CV_MIXT_YMOM,c1) = pDiss(CV_MIXT_YMOM,c1) + flx(3)
    pDiss(CV_MIXT_ZMOM,c1) = pDiss(CV_MIXT_ZMOM,c1) + flx(4)
    pDiss(CV_MIXT_ENER,c1) = pDiss(CV_MIXT_ENER,c1) + flx(5)

    pDiss(CV_MIXT_DENS,c2) = pDiss(CV_MIXT_DENS,c2) - flx(1)
    pDiss(CV_MIXT_XMOM,c2) = pDiss(CV_MIXT_XMOM,c2) - flx(2)
    pDiss(CV_MIXT_YMOM,c2) = pDiss(CV_MIXT_YMOM,c2) - flx(3)
    pDiss(CV_MIXT_ZMOM,c2) = pDiss(CV_MIXT_ZMOM,c2) - flx(4)
    pDiss(CV_MIXT_ENER,c2) = pDiss(CV_MIXT_ENER,c2) - flx(5)
  END DO ! ifg

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::ROE_ComputeFluxD1_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFluxD1_TCP








! ******************************************************************************
!
! Purpose: Compute dissipative fluxes using second-order accurate Roe scheme
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
!   1. Assumes ideal gas and constant Cp and R!
!
! ******************************************************************************

SUBROUTINE RFLU_ROE_ComputeFluxD2_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GHoVm2, & 
                           MixtPerf_Ho_CpTUVW, &
                           MixtPerf_R_CpG, &
                           MixtPerf_T_DPR

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indGs
  REAL(RFREAL) :: ah,betrk,cp,de,dissFact,dp,dq,dr,du,dv1,dv2,dv5,dw,dx, &
                  dy,dz,efc,epsentr,fs,g,gc,Hh,Hl,Hr,irl,irr,l1,l2,l5,nm,nx, &
                  ny,nz,pl,pr,qh,ql,qr,rl,rh,rr,sh,term,tl,tr,t1,t2,t3, &
                  t5,ul,uh,ur,vl,vh,vr,wl,wh,wr,wt,xc,yc,zc
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDiss,pDv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid    

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFluxD2_TCP',&
  'RFLU_ModRoeFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::ROE_ComputeFluxD2_TCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs = pGrid%indGs

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  
  pGc   => pRegion%mixt%gradCell  
  pDiss => pRegion%mixt%diss
  
! ******************************************************************************
! Define constants
! ******************************************************************************

  cp = global%refCp
  g  = global%refGamma
  gc = MixtPerf_R_CpG(cp,g)

  epsentr  = pRegion%mixtInput%epsentr
  betrk    = pRegion%mixtInput%betrk(pRegion%irkStep)
  dissFact = pRegion%mixtInput%dissFact
  term     = 0.5_RFREAL*betrk*dissFact

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
    Hr = MixtPerf_Ho_CpTUVW(cp,tr,ur,vr,wr)

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)
    wt = rl/(rl + rh)

    uh = wt*ul + (1.0_RFREAL-wt)*ur
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    Hh = wt*Hl + (1.0_RFREAL-wt)*Hr

    qh = uh*nx + vh*ny + wh*nz - fs
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    ah = MixtPerf_C_GHoVm2(g,Hh,2.0_RFREAL*sh) ! NOTE factor of 2

! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh - ah)
    l2 = ABS(qh     )
    l5 = ABS(qh + ah)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr = rr - rl
    du = ur - ul
    de = vr - vl ! use de to denote dv because dv is already used
    dw = wr - wl
    dp = pr - pl
    dq = du*nx + de*ny + dw*nz ! Note: no gs contribution

    dv1 = (dp - rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv5 = (dp + rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv2 = dr - dp/(ah*ah)

    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l2*rh
    t5 = l5*dv5

! ==============================================================================
!   Compute fluxes
! ==============================================================================

    flx(1) = term*(t1            + t2                                   + t5           )*nm
    flx(2) = term*(t1*(uh-nx*ah) + t2*uh + t3*(du-nx*dq)                + t5*(uh+nx*ah))*nm
    flx(3) = term*(t1*(vh-ny*ah) + t2*vh + t3*(de-ny*dq)                + t5*(vh+ny*ah))*nm
    flx(4) = term*(t1*(wh-nz*ah) + t2*wh + t3*(dw-nz*dq)                + t5*(wh+nz*ah))*nm
    flx(5) = term*(t1*(Hh-qh*ah) + t2*sh + t3*(uh*du+vh*de+wh*dw-qh*dq) + t5*(Hh+qh*ah))*nm

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pDiss(CV_MIXT_DENS,c1) = pDiss(CV_MIXT_DENS,c1) + flx(1)
    pDiss(CV_MIXT_XMOM,c1) = pDiss(CV_MIXT_XMOM,c1) + flx(2)
    pDiss(CV_MIXT_YMOM,c1) = pDiss(CV_MIXT_YMOM,c1) + flx(3)
    pDiss(CV_MIXT_ZMOM,c1) = pDiss(CV_MIXT_ZMOM,c1) + flx(4)
    pDiss(CV_MIXT_ENER,c1) = pDiss(CV_MIXT_ENER,c1) + flx(5)

    pDiss(CV_MIXT_DENS,c2) = pDiss(CV_MIXT_DENS,c2) - flx(1)
    pDiss(CV_MIXT_XMOM,c2) = pDiss(CV_MIXT_XMOM,c2) - flx(2)
    pDiss(CV_MIXT_YMOM,c2) = pDiss(CV_MIXT_YMOM,c2) - flx(3)
    pDiss(CV_MIXT_ZMOM,c2) = pDiss(CV_MIXT_ZMOM,c2) - flx(4)
    pDiss(CV_MIXT_ENER,c2) = pDiss(CV_MIXT_ENER,c2) - flx(5)
  END DO ! ifg

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::ROE_ComputeFluxD2_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFluxD2_TCP










! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate Roe scheme for
!   multiphase flow.
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

SUBROUTINE RFLU_ROE_ComputeFlux1_GL(pRegion)

  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtLiq_C2_Bp, &
                           MixtLiq_D_DoBpPPoBtTTo, &
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

  INTEGER :: c1,c2,ifg,indGs,indMf,iPatch
  REAL(RFREAL) :: Bgh2,Blh2,Bp,Bt,Bvh2,cmh,cml,cmr,cvg,cvl,cvv,Cgh2,Clh2, &
                  Cvh2,Cvm,dissFact,de,dp,dq,dr,dtemp,du,dv1,dv2,dv3,dv4, &
                  dv5,dvfg,dvfv,dw,efc,egh,el,elh,epsentr,er,evh,factor1, &
                  factor2,factor4,factor5,fdv1,fdv2,fdv3,fs,fxn1,fxn2, &
                  fxn3,fxn4,fxn5,fxn6,invc2p,invc2pb,irl,irr,l1,l2,l3,l4, &  
                  l5,nm,nx,ny,nz,ph,pl,po,pr,qh,ql,qr,ra11,ra12x,ra12y, &
                  ra12z,ra13,ra14,ra15,ra21,ra22x,ra22y,ra22z,ra23,ra24, &
                  ra25,ra31,ra32x,ra32y,ra32z,ra33,ra34,ra35,ra41,ra42x, &
                  ra42y,ra42z,ra43,ra44,ra45,ra51,ra52x,ra52y,ra52z,ra53, &
                  ra54,ra55,rgh,rgl,rgr,rh,rl,rlh,ro,rr,rvh,rvl,rvr, &
                  rYgl,rYgr,rYvl,rYvr,Rg,Rv,sh,t1,t2,t3,t4,t5,t6,th,tl, &
                  to,tr,uh,ul,ur,vfgh,vfgl,vfgr,vflh,vfvh,vfvl,vfvr,vh, & 
                  vl,vr,wh,wl,wr,wt              
  REAL(RFREAL) :: flxConv(5),flxDiss(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: fn,pCvMixt,pCvSpec,pDvMixt,pRhs
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_grid), POINTER :: pGrid   
 
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFlux1_GL',&
  'RFLU_ModRoeFlux.F90')

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************
 
  pGrid => pRegion%grid 
 
  indGs    = pRegion%grid%indGs
  indMf    = pRegion%mixtInput%indMfMixt
  epsentr  = pRegion%mixtInput%epsentr
  dissFact = pRegion%mixtInput%dissFact

  pCvMixt => pRegion%mixt%cv
  pDvMixt => pRegion%mixt%dv
  pCvSpec => pRegion%spec%cv
  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs

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

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------
    
    rl  = pCvMixt(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    ul  = pCvMixt(CV_MIXT_XMOM,c1)*irl
    vl  = pCvMixt(CV_MIXT_YMOM,c1)*irl
    wl  = pCvMixt(CV_MIXT_ZMOM,c1)*irl
    el  = pCvMixt(CV_MIXT_ENER,c1)*irl
    pl  = pDvMixt(DV_MIXT_PRES,c1)
    tl  = pDvMixt(DV_MIXT_TEMP,c1)
    cml = pDvMixt(DV_MIXT_SOUN,c1)
    
    rYgl = pCvSpec(1,c1)
    rYvl = pCvSpec(2,c1)    
    rvl  = MixtPerf_D_PRT(pl,Rv,tl)
    rgl  = MixtPerf_D_PRT(pl,Rg,tl)    
    vfvl = rYvl/rvl
    vfgl = rYgl/rgl
 
    ql = ul*nx + vl*ny + wl*nz - fs

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------
   
    rr  = pCvMixt(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    ur  = pCvMixt(CV_MIXT_XMOM,c2)*irr
    vr  = pCvMixt(CV_MIXT_YMOM,c2)*irr
    wr  = pCvMixt(CV_MIXT_ZMOM,c2)*irr
    er  = pCvMixt(CV_MIXT_ENER,c2)*irr
    pr  = pDvMixt(DV_MIXT_PRES,c2)
    tr  = pDvMixt(DV_MIXT_TEMP,c2)
    cmr = pDvMixt(DV_MIXT_SOUN,c2)
 
    rYgr = pCvSpec(1,c2)
    rYvr = pCvSpec(2,c2)   
    rvr  = MixtPerf_D_PRT(pr,Rv,tr)
    rgr  = MixtPerf_D_PRT(pr,Rg,tr)            
    vfvr = rYvr/rvr
    vfgr = rYgr/rgr
 
    qr = ur*nx + vr*ny + wr*nz - fs
   
! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)                     
    wt = rl/(rl + rh)                     
    uh = wt*ul + (1.0_RFREAL-wt)*ur     
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    ph = wt*pl + (1.0_RFREAL-wt)*pr      
    th = wt*tl + (1.0_RFREAL-wt)*tr       

    vfgh = wt*vfgl + (1.0_RFREAL-wt)*vfgr
    vfvh = wt*vfvl + (1.0_RFREAL-wt)*vfvr
    vflh = 1.0_RFREAL - vfvh - vfgh

    qh = uh*nx + vh*ny + wh*nz - fs 
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    elh = cvl*th + sh
    evh = cvv*th + sh
    egh = cvg*th + sh

    Clh2 = MixtLiq_C2_Bp(Bp) 
    Cvh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,th)
    Cgh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,th)

    rlh = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,ph,po,th,to)
    rvh = MixtPerf_D_PRT(ph,Rv,th)
    rgh = MixtPerf_D_PRT(ph,Rg,th)

    Blh2 = -Bt/Bp
    Bvh2 = rvh*Rv
    Bgh2 = rgh*Rg

    Cvm = (rlh*vflh*cvl + rvh*vfvh*cvv + rgh*vfgh*cvg)/rh
 
    cmh = MixtGasLiq_C(Cvm,rh,ph,rlh,rvh,rgh,vflh,vfvh,vfgh,Clh2,Cvh2,Cgh2, &
                       Blh2,Bvh2,Bgh2)
    
! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh      )
    l2 = ABS(qh      )
    l3 = ABS(qh      )
    l4 = ABS(qh - cmh)
    l5 = ABS(qh + cmh)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5 
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l3  = EntropyFixHartenHyman(l3,efc)
    l4  = EntropyFixHartenHyman(l4,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr    = rr - rl     
    du    = ur - ul
    de    = vr - vl
    dw    = wr - wl  
    dq    = du*nx + de*ny + dw*nz ! Note: no gs contribution
    dp    = pr - pl 
    dtemp = tr - tl
    dvfv  = vfvr - vfvl
    dvfg  = vfgr - vfgl

    fdv1 = ph/(rh*cmh*cmh*rh*Cvm*th)
    fdv2 = (vfgh*(rgh*Cgh2 - rh*cmh*cmh  &
         + (Bgh2*ph)/(rh*Cvm)))/(rgh*Cgh2*rh*cmh*cmh)
    fdv3 = (vfvh*(rvh*Cvh2 - rh*cmh*cmh  &
         + (Bvh2*ph)/(rh*Cvm)))/(rvh*Cvh2*rh*cmh*cmh)

    dv1  = (dtemp/th) - (dp*fdv1)
    dv2  = dvfg - dp*fdv2
    dv3  = dvfv - dp*fdv3
    dv4  = -dq/cmh + dp/(rh*cmh*cmh) 
    dv5  =  dq/cmh + dp/(rh*cmh*cmh) 

    invc2p  = vflh/Clh2 + vfgh/Cgh2 + vfvh/Cvh2 
    invc2pb = (Blh2*vflh)/Clh2 + (Bgh2*vfgh)/Cgh2 & 
            + (Bvh2*vfvh)/Cvh2

    ra11  = -th*invc2pb
    ra12x = -th*invc2pb*uh
    ra12y = -th*invc2pb*vh
    ra12z = -th*invc2pb*wh
    ra13  = rlh*vflh*cvl*th + Bt*cvl*th*th*vflh - th*invc2pb*sh
    ra14  = -(Bgh2*vfgh*th)/Cgh2
    ra15  = -(Bvh2*vfvh*th)/Cvh2

    ra21  =  rgh - rlh
    ra22x = (rgh - rlh)*uh
    ra22y = (rgh - rlh)*vh
    ra22z = (rgh - rlh)*wh
    ra23  =  rgh*egh - rlh*elh
    ra24  =  rgh
    ra25  =  0.0_RFREAL    

    ra31  =  rvh - rlh
    ra32x = (rvh - rlh)*uh
    ra32y = (rvh - rlh)*vh
    ra32z = (rvh - rlh)*wh
    ra33  =  rvh*evh - rlh*elh
    ra34  =  0.0_RFREAL
    ra35  =  rvh    
    
    fxn1  = ((rgh - rlh)*vfgh*(rgh*Cgh2 - rh*cmh*cmh &
          + (Bgh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rgh*Cgh2)
    fxn2  = ((rvh - rlh)*vfvh*(rvh*Cvh2 - rh*cmh*cmh &
          + (Bvh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rvh*Cvh2)
    factor1 = 0.5_RFREAL*rh*cmh*cmh*invc2p & 
            - (ph*invc2pb)/(2.0_RFREAL*rh*Cvm) + fxn1 + fxn2

    fxn3  = 0.5_RFREAL*rh*cmh*cmh*((elh*vflh)/Clh2 &
          + (egh*vfgh)/Cgh2 + (evh*vfvh)/Cvh2)
    fxn4  = (ph*(rlh*vflh*cvl + Bt*cvl*th*vflh &
          -   (sh*invc2pb)))/(2.0_RFREAL*rh*Cvm)
    fxn5  = ((rgh*egh - rlh*elh)*vfgh*(rgh*Cgh2 - rh*cmh*cmh &
          + (Bgh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rgh*Cgh2)
    fxn6  = ((rvh*evh - rlh*elh)*vfvh*(rvh*Cvh2 - rh*cmh*cmh &
          + (Bvh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rvh*Cvh2)
    factor2 = fxn3 + fxn4 + fxn5 + fxn6

    factor4 = 0.5_RFREAL*rgh*vfgh
    factor5 = 0.5_RFREAL*rvh*vfvh

    ra41  = factor1
    ra42x = factor1*uh - 0.5_RFREAL*rh*cmh*nx 
    ra42y = factor1*vh - 0.5_RFREAL*rh*cmh*ny
    ra42z = factor1*wh - 0.5_RFREAL*rh*cmh*nz
    ra43  = factor2    - 0.5_RFREAL*rh*cmh*qh
    ra44  = factor4
    ra45  = factor5

    ra51  = factor1
    ra52x = factor1*uh + 0.5_RFREAL*rh*cmh*nx 
    ra52y = factor1*vh + 0.5_RFREAL*rh*cmh*ny
    ra52z = factor1*wh + 0.5_RFREAL*rh*cmh*nz
    ra53  = factor2    + 0.5_RFREAL*rh*cmh*qh
    ra54  = factor4
    ra55  = factor5   
 
    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l3*dv3
    t4 = l4*dv4
    t5 = l5*dv5
    t6 = l1*rh 

! ==============================================================================
!   Compute fluxes
! ==============================================================================

! ------------------------------------------------------------------------------
!   Central part
! ------------------------------------------------------------------------------
  
    flxConv(1) = 0.5_RFREAL*(ql*rl            + qr*rr           )*nm
    flxConv(2) = 0.5_RFREAL*(ql*rl*ul + pl*nx + qr*rr*ur + pr*nx)*nm
    flxConv(3) = 0.5_RFREAL*(ql*rl*vl + pl*ny + qr*rr*vr + pr*ny)*nm
    flxConv(4) = 0.5_RFREAL*(ql*rl*wl + pl*nz + qr*rr*wr + pr*nz)*nm
    flxConv(5) = 0.5_RFREAL*(ql*(rl*el + pl) + pl*fs  + qr*(rr*er + pr) & 
                           + pr*fs)*nm

! ------------------------------------------------------------------------------
!   Dissipative part
! ------------------------------------------------------------------------------

    flxDiss(1) = 0.5_RFREAL*dissFact*(ra11*t1  + ra21*t2  + ra31*t3 &
                                    + ra41*t4  + ra51*t5)*nm 
    flxDiss(2) = 0.5_RFREAL*dissFact*(ra12x*t1 + ra22x*t2 + ra32x*t3 &
                                    + ra42x*t4 + ra52x*t5 + t6*(du-nx*dq))*nm
    flxDiss(3) = 0.5_RFREAL*dissFact*(ra12y*t1 + ra22y*t2 + ra32y*t3 &
                                    + ra42y*t4 + ra52y*t5 + t6*(de-ny*dq))*nm
    flxDiss(4) = 0.5_RFREAL*dissFact*(ra12z*t1 + ra22z*t2 + ra32z*t3 &
                                    + ra42z*t4 + ra52z*t5 + t6*(dw-nz*dq))*nm
    flxDiss(5) = 0.5_RFREAL*dissFact*(ra13*t1  + ra23*t2  + ra33*t3 &
                                    + ra43*t4  + ra53*t5 &
                                    + t6*(uh*du+vh*de+wh*dw-qh*dq))*nm

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flxConv(1) - flxDiss(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + (flxConv(5) - flxDiss(5))

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - (flxConv(5) - flxDiss(5))
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

END SUBROUTINE RFLU_ROE_ComputeFlux1_GL 








! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate Roe scheme
!  for calorically and thermally perfect gas.
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

SUBROUTINE RFLU_ROE_ComputeFlux1_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GHoVm2, &
                           MixtPerf_Ho_CpTUVW, &
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
  REAL(RFREAL) :: ah,cp,dissFact,dp,dq,dr,du,de,dv1,dv2,dv5,dw,efc, &
                  epsentr,fs,g,Hh,Hl,Hr,irl,irr,l1,l2,l5,nm,nx,ny,nz,pl, &
                  pr,qh,ql,qr,rl,rh,rr,sh,term,tl,tr,t1,t2,t3,t5,ul,uh,ur,vl, &
                  vh,vr,wl,wh,wr,wt
  REAL(RFREAL) :: flxConv(5),flxDiss(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid    
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFlux1_TCP',&
  'RFLU_ModRoeFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::ROE_ComputeFlux1_TCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs    = pRegion%grid%indGs
  indMf    = pRegion%mixtInput%indMfMixt
  indSd    = pRegion%mixtInput%indSd
  epsentr  = pRegion%mixtInput%epsentr
  dissFact = pRegion%mixtInput%dissFact

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv

  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd   

! ******************************************************************************
! Define constants
! ******************************************************************************

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

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)
    tl = pDv(DV_MIXT_TEMP,c1)

    Hl = MixtPerf_Ho_CpTUVW(cp,tl,ul,vl,wl)

    ql = ul*nx + vl*ny + wl*nz - fs

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)
    tr = pDv(DV_MIXT_TEMP,c2)

    Hr = MixtPerf_Ho_CpTUVW(cp,tr,ur,vr,wr)

    qr = ur*nx + vr*ny + wr*nz - fs

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)
    wt = rl/(rl + rh)

    uh = wt*ul + (1.0_RFREAL-wt)*ur
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    Hh = wt*Hl + (1.0_RFREAL-wt)*Hr

    qh = uh*nx + vh*ny + wh*nz - fs
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    ah = MixtPerf_C_GHoVm2(g,Hh,2.0_RFREAL*sh) ! NOTE factor of 2

! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh - ah)
    l2 = ABS(qh     )
    l5 = ABS(qh + ah)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr = rr - rl
    du = ur - ul
    de = vr - vl ! use de to denote dv because dv is already used
    dw = wr - wl
    dp = pr - pl
    dq = du*nx + de*ny + dw*nz ! Note: no gs contribution

    dv1 = (dp - rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv5 = (dp + rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv2 = dr - dp/(ah*ah)

    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l2*rh
    t5 = l5*dv5

! ==============================================================================
!   Compute fluxes
! ==============================================================================

! ------------------------------------------------------------------------------
!   Central part
! ------------------------------------------------------------------------------

    flxConv(1) = 0.5_RFREAL*(ql*rl            + qr*rr           )*nm
    flxConv(2) = 0.5_RFREAL*(ql*rl*ul + pl*nx + qr*rr*ur + pr*nx)*nm
    flxConv(3) = 0.5_RFREAL*(ql*rl*vl + pl*ny + qr*rr*vr + pr*ny)*nm
    flxConv(4) = 0.5_RFREAL*(ql*rl*wl + pl*nz + qr*rr*wr + pr*nz)*nm
    flxConv(5) = 0.5_RFREAL*(ql*rl*Hl + pl*fs + qr*rr*Hr + pr*fs)*nm

! ------------------------------------------------------------------------------
!   Dissipative part
! ------------------------------------------------------------------------------

    flxDiss(1) = 0.5_RFREAL*dissFact*(t1            + t2                                   + t5           )*nm
    flxDiss(2) = 0.5_RFREAL*dissFact*(t1*(uh-nx*ah) + t2*uh + t3*(du-nx*dq)                + t5*(uh+nx*ah))*nm
    flxDiss(3) = 0.5_RFREAL*dissFact*(t1*(vh-ny*ah) + t2*vh + t3*(de-ny*dq)                + t5*(vh+ny*ah))*nm
    flxDiss(4) = 0.5_RFREAL*dissFact*(t1*(wh-nz*ah) + t2*wh + t3*(dw-nz*dq)                + t5*(wh+nz*ah))*nm
    flxDiss(5) = 0.5_RFREAL*dissFact*(t1*(Hh-qh*ah) + t2*sh + t3*(uh*du+vh*de+wh*dw-qh*dq) + t5*(Hh+qh*ah))*nm

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flxConv(1) - flxDiss(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + (flxConv(5) - flxDiss(5))

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - (flxConv(5) - flxDiss(5))

! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + uh*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vh*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + wh*pMf(indMf*ifg)

    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - uh*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vh*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - wh*pMf(indMf*ifg)
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      CALL RFLU_NSCBC_CompFirstPatchFlux(pRegion,pPatch)
    ELSE
      CALL RFLU_CentralFirstPatch(pRegion,pPatch)
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::ROE_ComputeFlux1_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFlux1_TCP









! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate Roe scheme for
!   multiphase flow.
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

SUBROUTINE RFLU_ROE_ComputeFlux2_GL(pRegion)

  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtLiq_C2_Bp, &
                           MixtLiq_D_DoBpPPoBtTTo, &
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

  INTEGER :: c1,c2,ifg,indGs,indMf,iPatch
  REAL(RFREAL) :: Bgh2,Bgl2,Bgr2,Blh2,Bll2,Blr2,Bp,Bt,Bvh2,Bvl2,Bvr2,cmh, &
                  cml,cmr,cvg,cvl,cvv,Cgh2,Cgl2,Cgr2,Clh2,Cll2,Clr2,Cvh2, &
                  Cvl2,Cvr2,Cvm,dissFact,de,dp,dq,dr,dtemp,du,dv1,dv2,dv3, &
                  dv4,dv5,dvfg,dvfv,dw,dx,dy,dz,efc,egh,el,elh,epsentr,er, &
                  evh,factor1,factor2,factor4,factor5,fdv1,fdv2,fdv3,fs, &
                  fxn1,fxn2,fxn3,fxn4,fxn5,fxn6,invc2p,invc2pb,irl,irr,l1, &
                  l2,l3,l4,l5,nm,nx,ny,nz,ph,pl,po,pr,qh,ql,qr,ra11,ra12x, &
                  ra12y,ra12z,ra13,ra14,ra15,ra21,ra22x,ra22y,ra22z,ra23, &
                  ra24,ra25,ra31,ra32x,ra32y,ra32z,ra33,ra34,ra35,ra41, &
                  ra42x,ra42y,ra42z,ra43,ra44,ra45,ra51,ra52x,ra52y,ra52z, &
                  ra53,ra54,ra55,rel,rer,rgh,rgl,rgr,rh,rl,rlh,rll,rlr,ro, &
                  rr,rvh,rvl,rvr,rYgl,rYgr,rYvl,rYvr,Rg,Rv,sh,t1,t2,t3,t4, &
                  t5,t6,th,tl,to,tr,uh,ul,ur,vfgh,vfgl,vfgr,vflh,vfll,vflr, &
                  vfvh,vfvl,vfvr,vh,vl,vr,Vl2,Vr2,Ygl,Ygr,Yll,Ylr,Yvl,Yvr, & 
                  wh,wl,wr,wt,xc,yc,zc  
  REAL(RFREAL) :: flxConv(5),flxDiss(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pDvMixt,pRhs
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradMixt,pGradSpec
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFlux2_GL',&
  'RFLU_ModRoeFlux.F90')

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs    = pRegion%grid%indGs
  indMf    = pRegion%mixtInput%indMfMixt
  epsentr  = pRegion%mixtInput%epsentr
  dissFact = pRegion%mixtInput%dissFact

  pCvMixt => pRegion%mixt%cv
  pDvMixt => pRegion%mixt%dv
  pCvSpec => pRegion%spec%cv

  pGradMixt => pRegion%mixt%gradCell
  pGradSpec => pRegion%spec%gradCell

  pMf   => pRegion%mixt%mfMixt
  pRhs  => pRegion%mixt%rhs

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

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! -----------------------------------------------------------------------------
!   Left state
! -----------------------------------------------------------------------------    
  
    rl  = pCvMixt(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    ul = pCvMixt(CV_MIXT_XMOM,c1)*irl
    vl = pCvMixt(CV_MIXT_YMOM,c1)*irl
    wl = pCvMixt(CV_MIXT_ZMOM,c1)*irl
    pl = pDvMixt(DV_MIXT_PRES,c1)
    tl = pDvMixt(DV_MIXT_TEMP,c1)
    
    Ygl = pCvSpec(1,c1)*irl
    Yvl = pCvSpec(2,c1)*irl
         
    dx = xc - pGrid%cofg(XCOORD,c1)
    dy = yc - pGrid%cofg(YCOORD,c1)
    dz = zc - pGrid%cofg(ZCOORD,c1)    

    pl = pl + pGradMixt(XCOORD,GRC_MIXT_PRES,c1)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_PRES,c1)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_PRES,c1)*dz 
    ul = ul + pGradMixt(XCOORD,GRC_MIXT_XVEL,c1)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_XVEL,c1)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_XVEL,c1)*dz     
    vl = vl + pGradMixt(XCOORD,GRC_MIXT_YVEL,c1)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_YVEL,c1)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGradMixt(XCOORD,GRC_MIXT_ZVEL,c1)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_ZVEL,c1)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    tl = tl + pGradMixt(XCOORD,GRC_MIXT_TEMP,c1)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_TEMP,c1)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_TEMP,c1)*dz 
                    
    Ygl = Ygl + pGradSpec(XCOORD,1,c1)*dx &
              + pGradSpec(YCOORD,1,c1)*dy &
              + pGradSpec(ZCOORD,1,c1)*dz
    Yvl = Yvl + pGradSpec(XCOORD,2,c1)*dx &
              + pGradSpec(YCOORD,2,c1)*dy &
              + pGradSpec(ZCOORD,2,c1)*dz

    Yll = 1.0_RFREAL - Ygl - Yvl

    Vl2 = ul*ul + vl*vl + wl*wl
    ql  = ul*nx + vl*ny + wl*nz - fs

    rll = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pl,po,tl,to)
    rvl = MixtPerf_D_PRT(pl,Rv,tl)
    rgl = MixtPerf_D_PRT(pl,Rg,tl) 
    rl  = 1.0_RFREAL/(Yll/rll + Yvl/rvl + Ygl/rgl)

    vfgl = rl*Ygl/rgl
    vfvl = rl*Yvl/rvl
    vfll = 1.0_RFREAL - vfgl - vfvl

    rel = (rll*vfll*cvl + rgl*vfgl*cvg + rvl*vfvl*cvv)*tl + 0.5_RFREAL*rl*Vl2   
    el = rel/rl
      
    Cll2 = MixtLiq_C2_Bp(Bp) 
    Cvl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tl)
    Cgl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tl)

    Bll2 = -Bt/Bp
    Bvl2 = rvl*Rv
    Bgl2 = rgl*Rg
    
    Cvm = (rll*vfll*cvl + rvl*vfvl*cvv + rgl*vfgl*cvg)/rl
 
    cml = MixtGasLiq_C(Cvm,rl,pl,rll,rvl,rgl,vfll,vfvl,vfgl,Cll2,Cvl2,Cgl2, &
                       Bll2,Bvl2,Bgl2)
    
!----------------------------------------------------------------------------
!   Right state
! -----------------------------------------------------------------------------    
    
    rr  = pCvMixt(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    ur = pCvMixt(CV_MIXT_XMOM,c2)*irr
    vr = pCvMixt(CV_MIXT_YMOM,c2)*irr
    wr = pCvMixt(CV_MIXT_ZMOM,c2)*irr
    pr = pDvMixt(DV_MIXT_PRES,c2)
    tr = pDvMixt(DV_MIXT_TEMP,c2)
    
    Ygr = pCvSpec(1,c2)*irr
    Yvr = pCvSpec(2,c2)*irr 
       
    dx = xc - pGrid%cofg(XCOORD,c2)
    dy = yc - pGrid%cofg(YCOORD,c2)
    dz = zc - pGrid%cofg(ZCOORD,c2)    

    pr = pr + pGradMixt(XCOORD,GRC_MIXT_PRES,c2)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_PRES,c2)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_PRES,c2)*dz 
    ur = ur + pGradMixt(XCOORD,GRC_MIXT_XVEL,c2)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_XVEL,c2)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_XVEL,c2)*dz     
    vr = vr + pGradMixt(XCOORD,GRC_MIXT_YVEL,c2)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_YVEL,c2)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGradMixt(XCOORD,GRC_MIXT_ZVEL,c2)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_ZVEL,c2)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    tr = tr + pGradMixt(XCOORD,GRC_MIXT_TEMP,c2)*dx & 
            + pGradMixt(YCOORD,GRC_MIXT_TEMP,c2)*dy & 
            + pGradMixt(ZCOORD,GRC_MIXT_TEMP,c2)*dz 
                    
    Ygr = Ygr + pGradSpec(XCOORD,1,c2)*dx &
              + pGradSpec(YCOORD,1,c2)*dy &
              + pGradSpec(ZCOORD,1,c2)*dz
    Yvr = Yvr + pGradSpec(XCOORD,2,c2)*dx &
              + pGradSpec(YCOORD,2,c2)*dy &
              + pGradSpec(ZCOORD,2,c2)*dz

    Ylr = 1.0_RFREAL - Ygr - Yvr
 
    Vr2 = ur*ur + vr*vr + wr*wr
    qr  = ur*nx + vr*ny + wr*nz - fs
   
    rlr = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pr,po,tr,to)
    rvr = MixtPerf_D_PRT(pr,Rv,tr)
    rgr = MixtPerf_D_PRT(pr,Rg,tr)
    rr  = 1.0_RFREAL/(Ylr/rlr + Yvr/rvr + Ygr/rgr)  
   
    vfgr = rr*Ygr/rgr
    vfvr = rr*Yvr/rvr
    vflr = 1.0_RFREAL - vfgr - vfvr

    rer = (rlr*vflr*cvl + rgr*vfgr*cvg + rvr*vfvr*cvv)*tr + 0.5_RFREAL*rr*Vr2
    er  = rer/rr
    
    Clr2 = MixtLiq_C2_Bp(Bp) 
    Cvr2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tr)
    Cgr2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tr)

    Blr2 = -Bt/Bp
    Bvr2 = rvr*Rv
    Bgr2 = rgr*Rg
    
    Cvm = (rlr*vflr*cvl + rvr*vfvr*cvv + rgr*vfgr*cvg)/rr
 
    cmr = MixtGasLiq_C(Cvm,rr,pr,rlr,rvr,rgr,vflr,vfvr,vfgr,Clr2,Cvr2,Cgr2, &
                       Blr2,Bvr2,Bgr2)

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)                     
    wt = rl/(rl + rh)                     
    uh = wt*ul + (1.0_RFREAL-wt)*ur     
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    ph = wt*pl + (1.0_RFREAL-wt)*pr      
    th = wt*tl + (1.0_RFREAL-wt)*tr       

    vfgh = wt*vfgl + (1.0_RFREAL-wt)*vfgr
    vfvh = wt*vfvl + (1.0_RFREAL-wt)*vfvr
    vflh = 1.0_RFREAL - vfvh - vfgh

    qh = uh*nx + vh*ny + wh*nz - fs 
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    elh = cvl*th + sh
    evh = cvv*th + sh
    egh = cvg*th + sh

    Clh2 = MixtLiq_C2_Bp(Bp) 
    Cvh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,th)
    Cgh2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,th)

    rlh = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,ph,po,th,to)
    rvh = MixtPerf_D_PRT(ph,Rv,th)
    rgh = MixtPerf_D_PRT(ph,Rg,th)

    Blh2 = -Bt/Bp
    Bvh2 = rvh*Rv
    Bgh2 = rgh*Rg

    Cvm = (rlh*vflh*cvl + rvh*vfvh*cvv + rgh*vfgh*cvg)/rh
 
    cmh = MixtGasLiq_C(Cvm,rh,ph,rlh,rvh,rgh,vflh,vfvh,vfgh,Clh2,Cvh2,Cgh2, &
                       Blh2,Bvh2,Bgh2) 

! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh      )
    l2 = ABS(qh      )
    l3 = ABS(qh      )
    l4 = ABS(qh - cmh)
    l5 = ABS(qh + cmh)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5 
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l3  = EntropyFixHartenHyman(l3,efc)
    l4  = EntropyFixHartenHyman(l4,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr    = rr - rl     
    du    = ur - ul
    de    = vr - vl
    dw    = wr - wl  
    dq    = du*nx + de*ny + dw*nz   ! Note: no gs contribution
    dp    = pr - pl 
    dtemp = tr - tl
    dvfv  = vfvr - vfvl
    dvfg  = vfgr - vfgl

    fdv1 = ph/(rh*cmh*cmh*rh*Cvm*th)
    fdv2 = (vfgh*(rgh*Cgh2 - rh*cmh*cmh &
         + (Bgh2*ph)/(rh*Cvm)))/(rgh*Cgh2*rh*cmh*cmh)
    fdv3 = (vfvh*(rvh*Cvh2 - rh*cmh*cmh &
         + (Bvh2*ph)/(rh*Cvm)))/(rvh*Cvh2*rh*cmh*cmh)

    dv1  = (dtemp/th) - (dp*fdv1)
    dv2  = dvfg - dp*fdv2
    dv3  = dvfv - dp*fdv3
    dv4  = -dq/cmh + dp/(rh*cmh*cmh) 
    dv5  =  dq/cmh + dp/(rh*cmh*cmh) 

    invc2p  = vflh/Clh2 + vfgh/Cgh2 + vfvh/Cvh2 
    invc2pb = (Blh2*vflh)/Clh2 + (Bgh2*vfgh)/Cgh2 &
            + (Bvh2*vfvh)/Cvh2

    ra11  = -th*invc2pb
    ra12x = -th*invc2pb*uh
    ra12y = -th*invc2pb*vh
    ra12z = -th*invc2pb*wh
    ra13  = rlh*vflh*cvl*th + Bt*cvl*th*th*vflh - th*invc2pb*sh
    ra14  = -(Bgh2*vfgh*th)/Cgh2
    ra15  = -(Bvh2*vfvh*th)/Cvh2

    ra21  =  rgh - rlh
    ra22x = (rgh - rlh)*uh
    ra22y = (rgh - rlh)*vh
    ra22z = (rgh - rlh)*wh
    ra23  =  rgh*egh - rlh*elh
    ra24  =  rgh
    ra25  =  0.0_RFREAL    

    ra31  =  rvh - rlh
    ra32x = (rvh - rlh)*uh
    ra32y = (rvh - rlh)*vh
    ra32z = (rvh - rlh)*wh
    ra33  =  rvh*evh - rlh*elh
    ra34  =  0.0_RFREAL
    ra35  =  rvh    
    
    fxn1 = ((rgh - rlh)*vfgh*(rgh*Cgh2 - rh*cmh*cmh &
         + (Bgh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rgh*Cgh2)
    fxn2 = ((rvh - rlh)*vfvh*(rvh*Cvh2 - rh*cmh*cmh &
         + (Bvh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rvh*Cvh2)
    factor1 = 0.5_RFREAL*rh*cmh*cmh*invc2p & 
            - (ph*invc2pb)/(2.0_RFREAL*rh*Cvm) + fxn1 + fxn2

    fxn3 = 0.5_RFREAL*rh*cmh*cmh*((elh*vflh)/Clh2 &
         + (egh*vfgh)/Cgh2 + (evh*vfvh)/Cvh2)
    fxn4 = (ph*(rlh*vflh*cvl + Bt*cvl*th*vflh &
         - (sh*invc2pb)))/(2.0_RFREAL*rh*Cvm)
    fxn5 = ((rgh*egh - rlh*elh)*vfgh*(rgh*Cgh2 - rh*cmh*cmh &
         + (Bgh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rgh*Cgh2)
    fxn6 = ((rvh*evh - rlh*elh)*vfvh*(rvh*Cvh2 - rh*cmh*cmh &
         +  (Bvh2*ph)/(rh*Cvm)))/(2.0_RFREAL*rvh*Cvh2)
    factor2 = fxn3 + fxn4 + fxn5 + fxn6

    factor4 = 0.5_RFREAL*rgh*vfgh
    factor5 = 0.5_RFREAL*rvh*vfvh

    ra41  = factor1
    ra42x = factor1*uh - 0.5_RFREAL*rh*cmh*nx 
    ra42y = factor1*vh - 0.5_RFREAL*rh*cmh*ny
    ra42z = factor1*wh - 0.5_RFREAL*rh*cmh*nz
    ra43  = factor2    - 0.5_RFREAL*rh*cmh*qh
    ra44  = factor4
    ra45  = factor5

    ra51  = factor1
    ra52x = factor1*uh + 0.5_RFREAL*rh*cmh*nx 
    ra52y = factor1*vh + 0.5_RFREAL*rh*cmh*ny
    ra52z = factor1*wh + 0.5_RFREAL*rh*cmh*nz
    ra53  = factor2    + 0.5_RFREAL*rh*cmh*qh
    ra54  = factor4
    ra55  = factor5   
 
    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l3*dv3
    t4 = l4*dv4
    t5 = l5*dv5
    t6 = l1*rh 

! ==============================================================================
!   Compute fluxes
! ==============================================================================

! ------------------------------------------------------------------------------
!   Central part
! ------------------------------------------------------------------------------
  
    flxConv(1) = 0.5_RFREAL*(ql*rl            + qr*rr           )*nm
    flxConv(2) = 0.5_RFREAL*(ql*rl*ul + pl*nx + qr*rr*ur + pr*nx)*nm
    flxConv(3) = 0.5_RFREAL*(ql*rl*vl + pl*ny + qr*rr*vr + pr*ny)*nm
    flxConv(4) = 0.5_RFREAL*(ql*rl*wl + pl*nz + qr*rr*wr + pr*nz)*nm
    flxConv(5) = 0.5_RFREAL*(ql*(rl*el + pl) + pl*fs  + qr*(rr*er + pr) &
                           + pr*fs)*nm

! ------------------------------------------------------------------------------
!   Dissipative part
! ------------------------------------------------------------------------------

    flxDiss(1) = 0.5_RFREAL*dissFact*(ra11*t1  + ra21*t2  + ra31*t3 &
                                    + ra41*t4  + ra51*t5)*nm 
    flxDiss(2) = 0.5_RFREAL*dissFact*(ra12x*t1 + ra22x*t2 + ra32x*t3 &
                                    + ra42x*t4 + ra52x*t5 + t6*(du-nx*dq))*nm
    flxDiss(3) = 0.5_RFREAL*dissFact*(ra12y*t1 + ra22y*t2 + ra32y*t3 &
                                    + ra42y*t4 + ra52y*t5 + t6*(de-ny*dq))*nm
    flxDiss(4) = 0.5_RFREAL*dissFact*(ra12z*t1 + ra22z*t2 + ra32z*t3 &
                                    + ra42z*t4 + ra52z*t5 + t6*(dw-nz*dq))*nm
    flxDiss(5) = 0.5_RFREAL*dissFact*(ra13*t1  + ra23*t2  + ra33*t3 &
                                    + ra43*t4  + ra53*t5 &
                                    + t6*(uh*du+vh*de+wh*dw-qh*dq))*nm

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flxConv(1) - flxDiss(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + (flxConv(5) - flxDiss(5))

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - (flxConv(5) - flxDiss(5))
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      SELECT CASE ( pPatch%spaceOrder )
        CASE ( 1 ) 
          CALL RFLU_NSCBC_CompFirstPatchFlux(pRegion,pPatch)
        CASE ( 2 ) 
          CALL RFLU_NSCBC_CompSecondPatchFlux(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%spaceOrder
    ELSE
      SELECT CASE ( pPatch%spaceOrder )
        CASE ( 1 ) 
          CALL RFLU_CentralFirstPatch_GL(pRegion,pPatch)      
        CASE ( 2 ) 
          CALL RFLU_CentralSecondPatch_GL(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%spaceOrder
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFlux2_GL 










! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate Roe scheme
!  for calorically and thermally perfect gas.
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

SUBROUTINE RFLU_ROE_ComputeFlux2_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GHoVm2, &
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
  REAL(RFREAL) :: ah,cp,de,dissFact,dp,dq,dr,du,dv1,dv2,dv5,dw,dx,dy,dz, &
                  efc,epsentr,fs,g,gc,Hh,Hl,Hr,irl,irr,l1,l2,l5,nm,nx,ny,nz, &
                  pl,pr,qh,ql,qr,rl,rh,rr,sh,term,tl,tr,t1,t2,t3,t5,ul, &
                  uh,ur,vl,vh,vr,wl,wh,wr,wt,xc,yc,zc
  REAL(RFREAL) :: flxConv(5),flxDiss(5)
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

  CALL RegisterFunction(global,'RFLU_ROE_ComputeFlux2_TCP',&
  'RFLU_ModRoeFlux.F90')

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::ROE_ComputeFlux2_TCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs    = pRegion%grid%indGs
  indMf    = pRegion%mixtInput%indMfMixt
  indSd    = pRegion%mixtInput%indSd
  epsentr  = pRegion%mixtInput%epsentr
  dissFact = pRegion%mixtInput%dissFact

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd   
    
! ******************************************************************************
! Define constants
! ******************************************************************************

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

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)

    dx = xc - pGrid%cofg(XCOORD,c1)
    dy = yc - pGrid%cofg(YCOORD,c1)
    dz = zc - pGrid%cofg(ZCOORD,c1)

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
    Hl = MixtPerf_Ho_CpTUVW(cp,tl,ul,vl,wl)

    ql = ul*nx + vl*ny + wl*nz - fs

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)

    dx = xc - pGrid%cofg(XCOORD,c2)
    dy = yc - pGrid%cofg(YCOORD,c2)
    dz = zc - pGrid%cofg(ZCOORD,c2)

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
    Hr = MixtPerf_Ho_CpTUVW(cp,tr,ur,vr,wr)

    qr = ur*nx + vr*ny + wr*nz - fs

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)
    wt = rl/(rl + rh)

    uh = wt*ul + (1.0_RFREAL-wt)*ur
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    Hh = wt*Hl + (1.0_RFREAL-wt)*Hr

    qh = uh*nx + vh*ny + wh*nz - fs
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    ah = MixtPerf_C_GHoVm2(g,Hh,2.0_RFREAL*sh) ! NOTE factor of 2

! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh - ah)
    l2 = ABS(qh     )
    l5 = ABS(qh + ah)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr = rr - rl
    du = ur - ul
    de = vr - vl ! use de to denote dv because dv is already used
    dw = wr - wl
    dp = pr - pl
    dq = du*nx + de*ny + dw*nz ! Note: no gs contribution

    dv1 = (dp - rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv5 = (dp + rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv2 = dr - dp/(ah*ah)

    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l2*rh
    t5 = l5*dv5

! =============================================================================
!   Compute fluxes
! =============================================================================

! ------------------------------------------------------------------------------
!   Central part
! ------------------------------------------------------------------------------

    flxConv(1) = 0.5_RFREAL*(ql*rl            + qr*rr           )*nm
    flxConv(2) = 0.5_RFREAL*(ql*rl*ul + pl*nx + qr*rr*ur + pr*nx)*nm
    flxConv(3) = 0.5_RFREAL*(ql*rl*vl + pl*ny + qr*rr*vr + pr*ny)*nm
    flxConv(4) = 0.5_RFREAL*(ql*rl*wl + pl*nz + qr*rr*wr + pr*nz)*nm
    flxConv(5) = 0.5_RFREAL*(ql*rl*Hl + pl*fs + qr*rr*Hr + pr*fs)*nm

! ------------------------------------------------------------------------------
!   Dissipative part
! ------------------------------------------------------------------------------

    flxDiss(1) = 0.5_RFREAL*dissFact*(t1            + t2                                   + t5           )*nm
    flxDiss(2) = 0.5_RFREAL*dissFact*(t1*(uh-nx*ah) + t2*uh + t3*(du-nx*dq)                + t5*(uh+nx*ah))*nm
    flxDiss(3) = 0.5_RFREAL*dissFact*(t1*(vh-ny*ah) + t2*vh + t3*(de-ny*dq)                + t5*(vh+ny*ah))*nm
    flxDiss(4) = 0.5_RFREAL*dissFact*(t1*(wh-nz*ah) + t2*wh + t3*(dw-nz*dq)                + t5*(wh+nz*ah))*nm
    flxDiss(5) = 0.5_RFREAL*dissFact*(t1*(Hh-qh*ah) + t2*sh + t3*(uh*du+vh*de+wh*dw-qh*dq) + t5*(Hh+qh*ah))*nm

! =============================================================================
!   Store mass flux
! =============================================================================

    pMf(indMf*ifg) = flxConv(1) - flxDiss(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + (flxConv(5) - flxDiss(5))

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - (flxConv(1) - flxDiss(1))
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - (flxConv(2) - flxDiss(2))
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - (flxConv(3) - flxDiss(3))
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - (flxConv(4) - flxDiss(4))
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - (flxConv(5) - flxDiss(5))

! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + uh*pMf(indMf*ifg)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vh*pMf(indMf*ifg)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + wh*pMf(indMf*ifg)
   
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - uh*pMf(indMf*ifg)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vh*pMf(indMf*ifg)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - wh*pMf(indMf*ifg)
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      SELECT CASE ( pPatch%spaceOrder )
        CASE ( 1 ) 
          CALL RFLU_NSCBC_CompFirstPatchFlux(pRegion,pPatch)
        CASE ( 2 ) 
          CALL RFLU_NSCBC_CompSecondPatchFlux(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%spaceOrder
    ELSE
      SELECT CASE ( pPatch%spaceOrder )
        CASE ( 1 ) 
          CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
        CASE ( 2 ) 
          CALL RFLU_CentralSecondPatch(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%spaceOrder
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::ROE_ComputeFlux2_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ROE_ComputeFlux2_TCP








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModRoeFlux

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRoeFlux.F90,v $
! Revision 1.10  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/08/19 15:39:18  mparmar
! Added computation of boundary flux using boundary variables
!
! Revision 1.7  2006/05/03 17:55:42  haselbac
! Bug fix: Missing pGrid def in RFLU_ROE_ComputeFlux2_TCP
!
! Revision 1.6  2006/05/01 21:00:47  haselbac
! Rewrite for consistency and cleanliness
!
! Revision 1.5  2006/04/15 16:58:42  haselbac
! Added capability of running 1st order boundary fluxes with 2nd order volume fluxes
!
! Revision 1.4  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.3  2006/03/26 20:22:09  haselbac
! Added support for GL model
!
! Revision 1.2  2005/07/07 22:45:01  haselbac
! Added profiling calls
!
! Revision 1.1  2005/05/16 20:36:29  haselbac
! Initial revision
!
! ******************************************************************************
  
  















