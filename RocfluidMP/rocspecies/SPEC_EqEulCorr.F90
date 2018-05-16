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
! Purpose: Compute Equilibrium Eulerian correction to species mass fluxes and 
!   overall mass flux for interior faces.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to data of current region
!   iSpec          Index of species
!
! Output: None.
!
! Notes:
!  1. By passing iSpec, rather than looping over all iSpec inside the face
!     loop, this routine is optimized for the case of a few Eq Eul species 
!     among many non Eq Eul (e.g., gas species).
!
! ******************************************************************************
!
! $Id: SPEC_EqEulCorr.F90,v 1.8 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_EqEulCorr(pRegion,iSpec)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid  
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModSpecies, ONLY: t_spec_type  

  USE SPEC_ModParameters

  USE SPEC_ModInterfaces, ONLY: SPEC_EqEulCorrPatch

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iSpec
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,ifl,iPatch,iSpecEEv
  REAL(RFREAL) :: Eo1,Eo2,flx,iDens,ir1,ir2,nx,ny,nz,nm,phi1,phi2,p1,p2,rY, &
                  rY1,rY2,r1,r2,tauCoef,term,Ts1,Ts2,ug,ug1,ug2,us,us1,us2, &
                  vg,vg1,vg2,vs,vs1,vs2,wg,wg1,wg2,ws,ws1,ws2             
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pDvMixt,pSd, &
                                           rhsMixt,rhsSpec
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_spec_type), POINTER :: pSpecType  

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_EqEulCorr.F90,v $ $Revision: 1.8 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_EqEulCorr',&
  'SPEC_EqEulCorr.F90')

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid   => pRegion%grid
  pSd     => pRegion%mixt%sd
  pCvMixt => pRegion%mixt%cv
  pCvSpec => pRegion%spec%cv
  pDvMixt => pRegion%mixt%cv  
  rhsMixt => pRegion%mixt%rhs
  rhsSpec => pRegion%spec%rhs
  pEEv    => pRegion%spec%eev  

  pSpecType => pRegion%specInput%specType(iSpec)  

  tauCoef  = pSpecType%tauCoefficient
  iSpecEEv = pSpecType%iSpec2iSpecEEv
  iDens    = 1.0_RFREAL/pSpecType%pMaterial%dens

! ******************************************************************************
! Checks: Defensive coding, should never occur
! ******************************************************************************

  IF ( pSpecType%discreteFlag .EQV. .FALSE. ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pSpecType%discreteFlag

  IF ( pSpecType%velocityMethod /= SPEC_METHV_EQEUL ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pSpecType%velocityMethod

  IF ( tauCoef <= 0.0_RFREAL ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! tauCoef

  IF ( pRegion%mixtInput%indSd /= 1 ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pRegion%mixtInput%indSd

  IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%mixt%cvState
  
  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

  DO ifl = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifl)
    c2 = pGrid%f2c(2,ifl)

    nx = pGrid%fn(XCOORD,ifl)
    ny = pGrid%fn(YCOORD,ifl)
    nz = pGrid%fn(ZCOORD,ifl)
    nm = pGrid%fn(XYZMAG,ifl)

! ==============================================================================
!   Get data
! ==============================================================================

    ir1 = 1.0_RFREAL/pCvMixt(CV_MIXT_DENS,c1)
    ug1 = ir1*pCvMixt(CV_MIXT_XMOM,c1)
    vg1 = ir1*pCvMixt(CV_MIXT_YMOM,c1)
    wg1 = ir1*pCvMixt(CV_MIXT_ZMOM,c1)
    Eo1 = ir1*pCvMixt(CV_MIXT_ENER,c1)
    p1  =     pDvMixt(DV_MIXT_PRES,c1)
    
    ir2 = 1.0_RFREAL/pCvMixt(CV_MIXT_DENS,c2)
    ug2 = ir2*pCvMixt(CV_MIXT_XMOM,c2)
    vg2 = ir2*pCvMixt(CV_MIXT_YMOM,c2)
    wg2 = ir2*pCvMixt(CV_MIXT_ZMOM,c2) 
    Eo2 = ir2*pCvMixt(CV_MIXT_ENER,c2) 
    p2  =     pDvMixt(DV_MIXT_PRES,c2)                   

    rY1 = pCvSpec(iSpec,c1)
    rY2 = pCvSpec(iSpec,c2)
    
    phi1 = rY1*iDens
    phi2 = rY2*iDens

    us1 = pEEv(EEV_SPEC_XVEL,iSpecEEv,c1)
    vs1 = pEEv(EEV_SPEC_YVEL,iSpecEEv,c1)
    ws1 = pEEv(EEV_SPEC_ZVEL,iSpecEEv,c1) 
    Ts1 = pEEv(EEV_SPEC_TEMP,iSpecEEv,c1) 
    
    us2 = pEEv(EEV_SPEC_XVEL,iSpecEEv,c2)
    vs2 = pEEv(EEV_SPEC_YVEL,iSpecEEv,c2)
    ws2 = pEEv(EEV_SPEC_ZVEL,iSpecEEv,c2)
    Ts2 = pEEv(EEV_SPEC_TEMP,iSpecEEv,c2)    

! ==============================================================================
!   Compute fluxes
! ==============================================================================

    rY = 0.5_RFREAL*(rY1 + rY2)

    ug = 0.5_RFREAL*(ug1 + ug2)
    vg = 0.5_RFREAL*(vg1 + vg2)
    wg = 0.5_RFREAL*(wg1 + wg2)

    us = 0.5_RFREAL*(us1 + us2)
    vs = 0.5_RFREAL*(vs1 + vs2)
    ws = 0.5_RFREAL*(ws1 + ws2)
            
! TEMPORARY
    flx = ((us-ug)*nx + (vs-vg)*ny + (ws-wg)*nz)*nm    
!    flx = 0.0_RFREAL
! END TEMPORARY

    term = rY*flx

! ==============================================================================
!   Accumulate into residuals
! ==============================================================================

    rhsSpec(iSpec       ,c1) = rhsSpec(iSpec       ,c1) + term
    rhsSpec(iSpec       ,c2) = rhsSpec(iSpec       ,c2) - term

    rhsMixt(CV_MIXT_DENS,c1) = rhsMixt(CV_MIXT_DENS,c1) + term
    rhsMixt(CV_MIXT_DENS,c2) = rhsMixt(CV_MIXT_DENS,c2) - term    
    
    rhsMixt(CV_MIXT_XMOM,c1) = rhsMixt(CV_MIXT_XMOM,c1) + term*ug1
    rhsMixt(CV_MIXT_XMOM,c2) = rhsMixt(CV_MIXT_XMOM,c2) - term*ug2  
    
    rhsMixt(CV_MIXT_YMOM,c1) = rhsMixt(CV_MIXT_YMOM,c1) + term*vg1
    rhsMixt(CV_MIXT_YMOM,c2) = rhsMixt(CV_MIXT_YMOM,c2) - term*vg2  
    
    rhsMixt(CV_MIXT_ZMOM,c1) = rhsMixt(CV_MIXT_ZMOM,c1) + term*wg1
    rhsMixt(CV_MIXT_ZMOM,c2) = rhsMixt(CV_MIXT_ZMOM,c2) - term*wg2
    
    rhsMixt(CV_MIXT_ENER,c1) = rhsMixt(CV_MIXT_ENER,c1) + term*Eo1 + flx*p1*phi1
    rhsMixt(CV_MIXT_ENER,c2) = rhsMixt(CV_MIXT_ENER,c2) - term*Eo2 - flx*p2*phi2                
  END DO ! ifl

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch=1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    CALL SPEC_EqEulCorrPatch(pRegion,pPatch,iSpec)
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_EqEulCorr

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_EqEulCorr.F90,v $
! Revision 1.8  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.5  2005/11/27 02:08:47  haselbac
! Rewrote to take advantage of EEv
!
! Revision 1.4  2005/11/17 14:40:39  haselbac
! Added corrections to momentum and energy equations
!
! Revision 1.3  2005/11/10 02:34:15  haselbac
! Added support for gravity
!
! Revision 1.2  2005/03/31 17:17:15  haselbac
! Changed computation of correction, cosmetics
!
! Revision 1.1  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! ******************************************************************************







