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
! $Id: RFLU_ModViscousFlux.F90,v 1.10 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModViscousFlux

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_EnforceHeatFlux, & 
            RFLU_ViscousFluxes, &  
            RFLU_ViscousFluxesPatches
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModViscousFlux.F90,v $ $Revision: 1.10 $' 
              
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Enforce heat flux.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  tv           Transport variables
!  tvIndxCond   Index to conductity entry in tv
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_EnforceHeatFlux(pRegion,tv,tvIndxCond)
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: tvIndxCond
  REAL(RFREAL), DIMENSION(:,:) :: tv
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: c1,distrib,ifg,ifgBeg,ifgEnd,ifl,iPatch
  REAL(RFREAL) :: cond,dtdn,dtdx,dtdy,dtdz,nx,ny,nz
  REAL(RFREAL) :: tvf(tvIndxCond)
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EnforceHeatFlux',&
  'RFLU_ModViscousFlux.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

    distrib = pPatch%mixt%distrib

! ==============================================================================   
!   Select boundary type
! ==============================================================================   

    SELECT CASE ( pPatch%bcType )    

! ------------------------------------------------------------------------------
!     No-slip wall with imposed heat flux
! ------------------------------------------------------------------------------

      CASE ( BC_NOSLIPWALL_HFLUX )

! ----- Loop over faces and compute gradients ----------------------------------

        DO ifl = 1,pPatch%nBFaces
          c1 = pPatch%bf2c(ifl)

! ------- Get face geometry 
    
          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)       

! ------- Get face state 
   
          cond = tv(tvIndxCond,c1)

! TEMPORARY
!          CALL RFLU_InterpCells2FacePatch(pRegion,pPatch,ifl, &
!               pRegion%mixt%tv(tvIndxCond:tvIndxCond,:), &
!               tvf(tvIndxCond:tvIndxCond))
! 
!          cond = tvf(tvIndxCond)
! END TEMPORARY
    
! ------- Compute and set gradient

          dtdn = -pPatch%mixt%vals(BCDAT_NOSLIP_Q,distrib*ifl)/cond

          dtdx = dtdn*nx
          dtdy = dtdn*ny
          dtdz = dtdn*nz                  
    
! ------- Set gradient     

          pPatch%mixt%gradFace(XCOORD,GRBF_MIXT_TEMP,ifl) = dtdx
          pPatch%mixt%gradFace(YCOORD,GRBF_MIXT_TEMP,ifl) = dtdy
          pPatch%mixt%gradFace(ZCOORD,GRBF_MIXT_TEMP,ifl) = dtdz
        
        END DO ! ifl
        
! ------------------------------------------------------------------------------
!     Any other patches
! ------------------------------------------------------------------------------

      CASE DEFAULT
        
    END SELECT ! pPatch%bcType
  END DO ! iPatches

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EnforceHeatFlux







! ******************************************************************************
!
! Purpose: Compute viscous fluxes for actual faces.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  tv           Transport variables
!  tvIndxVisc   Index to viscosity entry in tv
!  tvIndxCond   Index to conductity entry in tv
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ViscousFluxes(pRegion,tv,tvIndxVisc,tvIndxCond)
   
  USE RFLU_ModInterpolation, ONLY: RFLU_InterpCells2Face 
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: tvIndxCond,tvIndxVisc
  REAL(RFREAL), DIMENSION(:,:) :: tv
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: c1,c2,ifg
  REAL(RFREAL), PARAMETER :: TWO_THIRDS = 2.0_RFREAL/3.0_RFREAL
  REAL(RFREAL) :: beta,cond,divTerm,dtdx,dtdy,dtdz,dudx,dudy,dudz,dvdx,dvdy, &
                  dvdz,dwdx,dwdy,dwdz,nm,nx,ny,nz,s11,s12,s13,s21,s22,s23, &
                  s31,s32,s33,u,v,visc,w
  REAL(RFREAL) :: fd(4)
  REAL(RFREAL) :: cvf(CV_MIXT_XVEL:CV_MIXT_ZVEL),tvf(tvIndxVisc:tvIndxCond)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDiss
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ViscousFluxes',&
  'RFLU_ModViscousFlux.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv
  pDiss => pRegion%mixt%diss
    
  beta = pRegion%mixtInput%betrk(pRegion%irkStep)  

! ******************************************************************************
! Check state of conserved state variable vector (defensive coding)
! ******************************************************************************

  IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_DUVWT ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%mixt%cvState
  
! ******************************************************************************  
! Loop over faces and compute viscous fluxes
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================   
!   Get face geometry
! ==============================================================================    
    
    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg) 
    nm = pGrid%fn(XYZMAG,ifg)        

! ==============================================================================   
!   Get face state
! ==============================================================================

    u = 0.5_RFREAL*(pCv(CV_MIXT_XVEL,c1) + pCv(CV_MIXT_XVEL,c2))
    v = 0.5_RFREAL*(pCv(CV_MIXT_YVEL,c1) + pCv(CV_MIXT_YVEL,c2))
    w = 0.5_RFREAL*(pCv(CV_MIXT_ZVEL,c1) + pCv(CV_MIXT_ZVEL,c2))

    visc = 0.5_RFREAL*(tv(tvIndxVisc,c1) + tv(tvIndxVisc,c2))    
    cond = 0.5_RFREAL*(tv(tvIndxCond,c1) + tv(tvIndxCond,c2))

! TEMPORARY
!    CALL RFLU_InterpCells2Face(pRegion,ifg, &
!                               pRegion%mixt%cv(CV_MIXT_XVEL:CV_MIXT_ZVEL,:), &
!                               cvf(CV_MIXT_XVEL:CV_MIXT_ZVEL))
! 
!    CALL RFLU_InterpCells2Face(pRegion,ifg, &
!                               pRegion%mixt%tv(tvIndxVisc:tvIndxCond,:), &
!                               tvf(tvIndxVisc:tvIndxCond))
! 
!    u = cvf(CV_MIXT_XVEL)
!    v = cvf(CV_MIXT_YVEL)
!    w = cvf(CV_MIXT_ZVEL)
!    
!    visc = tvf(tvIndxVisc)
!    cond = tvf(tvIndxCond)
! END TEMPORARY

! ==============================================================================   
!   Get gradients
! ==============================================================================

    dudx = pRegion%mixt%gradFace(XCOORD,GRF_MIXT_XVEL,ifg)
    dudy = pRegion%mixt%gradFace(YCOORD,GRF_MIXT_XVEL,ifg)
    dudz = pRegion%mixt%gradFace(ZCOORD,GRF_MIXT_XVEL,ifg)
    
    dvdx = pRegion%mixt%gradFace(XCOORD,GRF_MIXT_YVEL,ifg)
    dvdy = pRegion%mixt%gradFace(YCOORD,GRF_MIXT_YVEL,ifg)
    dvdz = pRegion%mixt%gradFace(ZCOORD,GRF_MIXT_YVEL,ifg)

    dwdx = pRegion%mixt%gradFace(XCOORD,GRF_MIXT_ZVEL,ifg)
    dwdy = pRegion%mixt%gradFace(YCOORD,GRF_MIXT_ZVEL,ifg)
    dwdz = pRegion%mixt%gradFace(ZCOORD,GRF_MIXT_ZVEL,ifg)

    dtdx = pRegion%mixt%gradFace(XCOORD,GRF_MIXT_TEMP,ifg)
    dtdy = pRegion%mixt%gradFace(YCOORD,GRF_MIXT_TEMP,ifg)
    dtdz = pRegion%mixt%gradFace(ZCOORD,GRF_MIXT_TEMP,ifg)
  
! ==============================================================================   
!   Compute fluxes
! ==============================================================================
  
    divTerm = TWO_THIRDS*(dudx + dvdy + dwdz)
  
    s11 = 2.0_RFREAL*dudx - divTerm
    s12 = dudy + dvdx
    s13 = dudz + dwdx 
  
    s21 = s12
    s22 = 2.0_RFREAL*dvdy - divTerm
    s23 = dvdz + dwdy  
  
    s31 = s13
    s32 = s23
    s33 = 2.0_RFREAL*dwdz - divTerm
  
    fd(1) = visc*(s11*nx + s12*ny + s13*nz)*nm
    fd(2) = visc*(s21*nx + s22*ny + s23*nz)*nm
    fd(3) = visc*(s31*nx + s32*ny + s33*nz)*nm
    
    fd(4) = u*fd(1) + v*fd(2) + w*fd(3) + cond*(dtdx*nx + dtdy*ny + dtdz*nz)*nm   

! ==============================================================================   
!   Accumulate into residual     
! ==============================================================================

    pDiss(CV_MIXT_XMOM,c1) = pDiss(CV_MIXT_XMOM,c1) + beta*fd(1)
    pDiss(CV_MIXT_YMOM,c1) = pDiss(CV_MIXT_YMOM,c1) + beta*fd(2)
    pDiss(CV_MIXT_ZMOM,c1) = pDiss(CV_MIXT_ZMOM,c1) + beta*fd(3)
    pDiss(CV_MIXT_ENER,c1) = pDiss(CV_MIXT_ENER,c1) + beta*fd(4)

    pDiss(CV_MIXT_XMOM,c2) = pDiss(CV_MIXT_XMOM,c2) - beta*fd(1)
    pDiss(CV_MIXT_YMOM,c2) = pDiss(CV_MIXT_YMOM,c2) - beta*fd(2)
    pDiss(CV_MIXT_ZMOM,c2) = pDiss(CV_MIXT_ZMOM,c2) - beta*fd(3)
    pDiss(CV_MIXT_ENER,c2) = pDiss(CV_MIXT_ENER,c2) - beta*fd(4)
  END DO ! ifg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ViscousFluxes









! ******************************************************************************
!
! Purpose: Compute viscous fluxes for actual faces on boundary patches.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  tv           Transport variables
!  tvIndxVisc   Index to viscosity entry in tv
!  tvIndxCond   Index to conductity entry in tv
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ViscousFluxesPatches(pRegion,tv,tvIndxVisc,tvIndxCond)
    
  USE RFLU_ModInterpolation, ONLY: RFLU_InterpCells2FacePatch  
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: tvIndxCond,tvIndxVisc
  REAL(RFREAL), DIMENSION(:,:) :: tv
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: c1,ifg,ifgBeg,ifgEnd,ifl,iPatch
  REAL(RFREAL), PARAMETER :: TWO_THIRDS = 2.0_RFREAL/3.0_RFREAL
  REAL(RFREAL) :: beta,cond,divTerm,dtdx,dtdy,dtdz,dudx,dudy,dudz,dvdx,dvdy, &
                  dvdz,dwdx,dwdy,dwdz,iCfRef,iChRef,nm,nx,ny,nz,rRef,s11, &
                  s12,s13,s21,s22,s23,s31,s32,s33,visc,vRef
  REAL(RFREAL) :: fd(4)
  REAL(RFREAL) :: tvf(tvIndxVisc:tvIndxCond)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDiss
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ViscousFluxesPatches',&
  'RFLU_ModViscousFlux.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv
  pDiss => pRegion%mixt%diss
  
  beta = pRegion%mixtInput%betrk(pRegion%irkStep)  

  rRef = global%refDensity
  vRef = global%refVelocity
  
  iCfRef = 2.0_RFREAL/(rRef*vRef*vRef)
  iChRef = 2.0_RFREAL/(rRef*vRef*vRef*vRef)

! ******************************************************************************
! Check state of conserved state variable vector (defensive coding)
! ******************************************************************************

  IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_DUVWT ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%mixt%cvState
  
! ******************************************************************************  
! Loop over patches 
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

! ==============================================================================   
!   Select boundary type
! ==============================================================================   

    SELECT CASE ( pPatch%bcType )    

! ------------------------------------------------------------------------------
!     No-slip wall
! ------------------------------------------------------------------------------

      CASE ( BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP )

! ----- Loop over faces and compute gradients ----------------------------------

        DO ifl = 1,pPatch%nBFaces
          c1 = pPatch%bf2c(ifl)

! ------- Get face geometry 
    
          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl) 
          nm = pPatch%fn(XYZMAG,ifl)        

! ------- Get face state 

          visc = tv(tvIndxVisc,c1)   
          cond = tv(tvIndxCond,c1)

! TEMPORARY
!          CALL RFLU_InterpCells2FacePatch(pRegion,pPatch,ifl, &
!               pRegion%mixt%tv(tvIndxVisc:tvIndxCond,:), &
!               tvf(tvIndxVisc:tvIndxCond))
! 
!          visc = tvf(tvIndxVisc)   
!          cond = tvf(tvIndxCond)
! END TEMPORARY

! ------- Get gradients 

          dudx = pPatch%mixt%gradFace(XCOORD,GRBF_MIXT_XVEL,ifl)
          dudy = pPatch%mixt%gradFace(YCOORD,GRBF_MIXT_XVEL,ifl)
          dudz = pPatch%mixt%gradFace(ZCOORD,GRBF_MIXT_XVEL,ifl)

          dvdx = pPatch%mixt%gradFace(XCOORD,GRBF_MIXT_YVEL,ifl)
          dvdy = pPatch%mixt%gradFace(YCOORD,GRBF_MIXT_YVEL,ifl)
          dvdz = pPatch%mixt%gradFace(ZCOORD,GRBF_MIXT_YVEL,ifl)

          dwdx = pPatch%mixt%gradFace(XCOORD,GRBF_MIXT_ZVEL,ifl)
          dwdy = pPatch%mixt%gradFace(YCOORD,GRBF_MIXT_ZVEL,ifl)
          dwdz = pPatch%mixt%gradFace(ZCOORD,GRBF_MIXT_ZVEL,ifl)

          dtdx = pPatch%mixt%gradFace(XCOORD,GRBF_MIXT_TEMP,ifl)
          dtdy = pPatch%mixt%gradFace(YCOORD,GRBF_MIXT_TEMP,ifl)
          dtdz = pPatch%mixt%gradFace(ZCOORD,GRBF_MIXT_TEMP,ifl)
   
! ------- Compute fluxes 
  
          divTerm = TWO_THIRDS*(dudx + dvdy + dwdz)

          s11 = 2.0_RFREAL*dudx - divTerm
          s12 = dudy + dvdx
          s13 = dudz + dwdx 

          s21 = s12
          s22 = 2.0_RFREAL*dvdy - divTerm
          s23 = dvdz + dwdy  

          s31 = s13
          s32 = s23
          s33 = 2.0_RFREAL*dwdz - divTerm

          fd(1) = visc*(s11 *nx + s12 *ny + s13 *nz)
          fd(2) = visc*(s21 *nx + s22 *ny + s23 *nz)
          fd(3) = visc*(s31 *nx + s32 *ny + s33 *nz)
          fd(4) = cond*(dtdx*nx + dtdy*ny + dtdz*nz)
    
! ------- Set friction and heat-transfer coefficients

          pPatch%cf(XCOORD,ifl) = -iCfRef*fd(1)
          pPatch%cf(YCOORD,ifl) = -iCfRef*fd(2)
          pPatch%cf(ZCOORD,ifl) = -iCfRef*fd(3)       
          
          pPatch%ch(ifl) = iChRef*fd(4)               
    
! ------- Accumulate into residual     

          pDiss(CV_MIXT_XMOM,c1) = pDiss(CV_MIXT_XMOM,c1) + beta*fd(1)*nm
          pDiss(CV_MIXT_YMOM,c1) = pDiss(CV_MIXT_YMOM,c1) + beta*fd(2)*nm
          pDiss(CV_MIXT_ZMOM,c1) = pDiss(CV_MIXT_ZMOM,c1) + beta*fd(3)*nm
          pDiss(CV_MIXT_ENER,c1) = pDiss(CV_MIXT_ENER,c1) + beta*fd(4)*nm      
        
        END DO ! ifl
        
! ------------------------------------------------------------------------------
!     Any other patches
! ------------------------------------------------------------------------------

      CASE DEFAULT
        
    END SELECT ! pPatch%bcType
  END DO ! iPatches

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ViscousFluxesPatches








  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModViscousFlux


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModViscousFlux.F90,v $
! Revision 1.10  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2007/08/07 15:18:30  rfiedler
! Bug fix: Wrong loop limit in RFLU_EnforceHeatFlux.
!
! Revision 1.7  2006/08/19 15:39:19  mparmar
! Renamed bGradFace, removed bf2bg, used GRBF_ for boundary grad arrays
!
! Revision 1.6  2006/04/07 15:19:21  haselbac
! Removed tabs
!
! Revision 1.5  2005/10/16 18:03:37  haselbac
! Bug fix: Missing definition of pGrid
!
! Revision 1.4  2005/10/16 17:16:42  haselbac
! Bug fix: distrib set in wrong place
!
! Revision 1.3  2005/10/14 14:08:41  haselbac
! Added RFLU_EnforceHeatFlux - temporary until have proper constr reconstr
!
! Revision 1.2  2005/10/05 14:14:04  haselbac
! No longer distinguish between isothermal and adiabatic walls
!
! Revision 1.1  2005/05/16 20:36:29  haselbac
! Initial revision
!
! ******************************************************************************
  
  









