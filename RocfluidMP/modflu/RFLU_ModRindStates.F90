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
! Purpose: Collection of routines for setting rind states.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModRindStates.F90,v 1.6 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRindStates

  USE ModDataTypes
  USE ModParameters
  USE ModGlobal, ONLY: t_global

  USE ModInterfaces, ONLY: MixtPerf_C_DGP, &
                           MixtPerf_C_GRT, &  
                           MixtPerf_D_PRT, & 
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_Eo_DGPVm, &  
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_G_CpR, &                           
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_R_M 

  IMPLICIT NONE
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================  
! Procedures
! ==============================================================================  

  PUBLIC :: RFLU_SetRindStateFarfieldPerf, & 
            RFLU_SetRindStateSlipWallPerf, & 
            RFLU_SetRindStateInjectPerf

! ==============================================================================  
! Data
! ==============================================================================  

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModRindStates.F90,v $ $Revision: 1.6 $'        
       
! ******************************************************************************
! Procedures
! ******************************************************************************
                
  CONTAINS
  







! ******************************************************************************
!
! Purpose: Set rind state for farfield boundaries and perfect gas.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   cpGas       Specific heat at constant pressure
!   mmGas       Molecular mass
!   nx,ny,nz    Components of unit normal vector
!   machInf     Mach number at infinity
!   pInf        Pressure at infinity
!   tInf        Temperature at infinity
!   alphaInf    Angle of attack
!   betaInf     Sideslip angle
!   corrFlag    Flag for vortex-correction
!   liftCoef    Lift coefficient
!   xc          x-coordinate 
!   yc          y-coordinate
!   zc          z-coordinate
!   rl          Density at boundary
!   rul         x-momentum component at boundary
!   rvl         y-momentum component at boundary
!   rwl         z-momentum component at boundary
!   rel         Total internal energy at boundary
!
! Output: 
!   rr          Density 
!   rur         x-momentum component 
!   rvr         y-momentum component 
!   rwr         z-momentum component 
!   rer         Total internal energy 
!   pr          Pressure 
!
! Notes: 
!   1. Valid only for thermally and calorically perfect gas.
!   2. Valid only for two-dimensional flows.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetRindStateFarfieldPerf(global,cpGas,mmGas,nx,ny,nz, &
                                           machInf,pInf,tInf,alphaInf, &
                                           betaInf,corrFlag,liftCoef,xc,yc, &
                                           zc,rl,rul,rvl,rwl,rel,rr,rur,rvr, &
                                           rwr,rer,pr)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    LOGICAL, INTENT(IN) :: corrFlag
    REAL(RFREAL), INTENT(IN) :: alphaInf,betaInf,cpGas,liftCoef,machInf, &
                                mmGas,nx,ny,nz,pInf,rl,rel,rul,rvl,rwl, & 
                                tInf,xc,yc,zc
    REAL(RFREAL), INTENT(OUT) :: pr,rer,rr,rur,rvr,rwr 
    TYPE(t_global), POINTER :: global

! ==============================================================================  
!   Locals 
! ==============================================================================  

    REAL(RFREAL) :: al,corr,corrTerm,denom,dist,dq2,dx,dy,el,gGas,gm1og, & 
                    gogm1,gGasTerm3,numer,pb,pi,pl,qi,ql,rb,rGas,ri,sl2, &
                    term,theta,ub,ui,ul,vb,vi,vl,wb,wi,wl    

! ******************************************************************************
!   Compute gas properties
! ******************************************************************************

    rGas = MixtPerf_R_M(mmGas)
    gGas = MixtPerf_G_CpR(cpGas,rGas)
          
    gm1og = (gGas-1.0_RFREAL)/gGas
    gogm1 = 1.0_RFREAL/gm1og  
          
    corrTerm = global%forceRefLength/(4.0_RFREAL*global%pi)
          
! ******************************************************************************
!   Interior state at boundary
! ******************************************************************************
        
    ul = rul/rl
    vl = rvl/rl
    wl = rwl/rl
    ql = ul*nx + vl*ny + wl*nz
 
    el  = rel/rl   
    sl2 = ul*ul + vl*vl + wl*wl
    pl  = MixtPerf_P_DEoGVm2(rl,el,gGas,sl2)
          
! ******************************************************************************
!   Compute state at infinity without vortex correction
! ******************************************************************************
 
    ri = MixtPerf_D_PRT(pInf,rGas,tInf)
    pi = pInf    
    
    qi = machInf*MixtPerf_C_GRT(gGas,rGas,tInf)
    ui = qi*COS(alphaInf)*COS(betaInf)
    vi = qi*SIN(alphaInf)*COS(betaInf)
    wi = qi*              SIN(betaInf)    

! ******************************************************************************
!   Compute state at infinity with vortex correction. NOTE the correction is 
!   assumed to be two-dimensional and the aerofoil center of pressure is assumed
!   to be located at (x,y) = (0.25,0.0).
! ******************************************************************************
 
    IF ( corrFlag .EQV. .TRUE. ) THEN 
      dx = xc - 0.25_RFREAL
      dy = yc
      
      dist  = SQRT(dx**2 + dy**2)
      theta = ATAN2(dy,dx)
      
      numer = liftCoef*qi*SQRT(1.0_RFREAL-machInf**2)
      denom = dist*(1.0_RFREAL - (machInf*SIN(theta-alphaInf))**2)
      corr  = corrTerm*numer/denom
      
      ui = ui + corr*SIN(theta)
      vi = vi - corr*COS(theta)
            
      dq2 = qi*qi - (ui*ui + vi*vi)            
      pi  = (pi**gm1og + 0.5_RFREAL*gm1og*ri/pi**(1.0_RFREAL/gGas)*dq2)**gogm1         
      ri  = ri*(pi/pInf)**gGas
    END IF ! corrFlag

! ******************************************************************************
!   Compute right state at boundary
! ******************************************************************************

! ==============================================================================  
!   Subsonic flow
! ==============================================================================  
   
    IF ( machInf < 1.0_RFREAL ) THEN
      al = MixtPerf_C_DGP(rl,gGas,pl)

! ------------------------------------------------------------------------------
!     Subsonic inflow
! ------------------------------------------------------------------------------

      IF ( ql < 0.0_RFREAL ) THEN
        pb = 0.5_RFREAL*(pi+pl-rl*al*((ui-ul)*nx+(vi-vl)*ny+(wi-wl)*nz))

        rb = ri -    (pi - pb)/(al*al)        
        ub = ui - nx*(pi - pb)/(rl*al)
        vb = vi - ny*(pi - pb)/(rl*al)
        wb = wi - nz*(pi - pb)/(rl*al)

! ------------------------------------------------------------------------------
!     Subsonic outflow
! ------------------------------------------------------------------------------

      ELSE
        pb = pi
      
        rb = rl -    (pl-pi)/(al*al)
        ub = ul + nx*(pl-pi)/(rl*al)
        vb = vl + ny*(pl-pi)/(rl*al)
        wb = wl + nz*(pl-pi)/(rl*al)
      END IF ! ql

      rr  = rb
      rur = rb*ub
      rvr = rb*vb
      rwr = rb*wb
      rer = rb*MixtPerf_Eo_DGPUVW(rb,gGas,pb,ub,vb,wb)
      pr  = pb

! ==============================================================================  
!   Supersonic flow
! ==============================================================================  

    ELSE

! ------------------------------------------------------------------------------
!     Supersonic inflow
! ------------------------------------------------------------------------------

      IF ( ql < 0.0_RFREAL ) THEN
        rr  = ri
        rur = ri*ui
        rvr = ri*vi
        rwr = ri*wi
        rer = ri*MixtPerf_Eo_DGPVm(ri,gGas,pi,qi)
        pr  = pi

! ------------------------------------------------------------------------------
!     Supersonic outflow
! ------------------------------------------------------------------------------

      ELSE
        rr  = rl
        rur = rul
        rvr = rvl
        rwr = rwl
        rer = rel
        pr  = pl
      END IF ! ql
    END IF ! machInf
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetRindStateFarfieldPerf











! ******************************************************************************
!
! Purpose: Set rind state for injection boundaries and perfect gas.
!
! Description: None.
!
! Input:
!   cpGas       Specific heat at constant pressure
!   mmGas       Molecular mass
!   nx,ny,nz    Components of unit normal vector
!   mInj        Injection mass flux
!   tInj        Injection temperature
!   pl          Pressure
!   fs          Grid speed
!
! Output: 
!   rl          Density
!   ul          x-velocity component
!   vl          y-velocity component
!   wl          z-velocity component
!   Hl          Stagnation enthalpy per unit mass
!
! Notes: 
!   1. Valid only for thermally and calorically perfect gas.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetRindStateInjectPerf(cpGas,mmGas,nx,ny,nz,mInj,tInj,pl, &
                                         fs,rl,ul,vl,wl,Hl)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    REAL(RFREAL), INTENT(IN) :: cpGas,fs,mInj,mmGas,nx,ny,nz,pl,tInj
    REAL(RFREAL), INTENT(OUT) :: Hl,rl,ul,vl,wl 

! ==============================================================================  
!   Locals 
! ==============================================================================  

    REAL(RFREAL) :: gGas,ql,rGas    
          
! ******************************************************************************
!   Compute wall pressure
! ******************************************************************************

    rGas = MixtPerf_R_M(mmGas)
    gGas = MixtPerf_G_CpR(cpGas,rGas)
 
    rl = MixtPerf_D_PRT(pl,rGas,tInj)
 
    ql = -mInj/rl + fs
    ul = ql*nx
    vl = ql*ny
    wl = ql*nz
     
    Hl = MixtPerf_Ho_CpTUVW(cpGas,tInj,ul,vl,wl)
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetRindStateInjectPerf








! ******************************************************************************
!
! Purpose: Set rind state for slip-wall boundaries and perfect gas.
!
! Description: None.
!
! Input:
!   cpGas       Specific heat at constant pressure
!   mmGas       Molecular mass
!   nx,ny,nz    Components of unit normal vector
!   rl          Density
!   rul         x-momentum component
!   rvl         y-momentum component
!   rwl         z-momentum component
!   fs          Grid speed
!   pl          Pressure
!
! Output: 
!   pl          Pressure
!
! Notes: 
!   1. Valid only for thermally and calorically perfect gas.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetRindStateSlipWallPerf(cpGas,mmGas,nx,ny,nz,rl,rul,rvl, &
                                           rwl,fs,pl)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    REAL(RFREAL), INTENT(IN) :: cpGas,fs,mmGas,nx,ny,nz,rl,rul,rvl,rwl
    REAL(RFREAL), INTENT(INOUT) :: pl

! ==============================================================================  
!   Locals 
! ==============================================================================  

    REAL(RFREAL) :: al,gGas,irl,ql,rGas,term,ul,vl,wl    
          
! ******************************************************************************
!   Compute wall pressure
! ******************************************************************************

    rGas = MixtPerf_R_M(mmGas)
    gGas = MixtPerf_G_CpR(cpGas,rGas)
 
    irl = 1.0_RFREAL/rl          
    ul  = irl*rul
    vl  = irl*rvl
    wl  = irl*rwl 
    ql  = ul*nx + vl*ny + wl*nz - fs
 
    al  = MixtPerf_C_DGP(rl,gGas,pl)

    IF ( ql < 0.0_RFREAL ) THEN           
      term = 1.0_RFREAL + 0.5_RFREAL*(gGas-1.0_RFREAL)*ql/al          
      pl   = pl*term**(2.0_RFREAL*gGas/(gGas-1.0_RFREAL))
    ELSE 
      term = (gGas+1.0_RFREAL)/4.0_RFREAL
      pl   = pl + term*rl*ql*(ql + SQRT(al*al + term*term*ql*ql)/term)
    END IF ! ql
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetRindStateSlipWallPerf





END MODULE RFLU_ModRindStates

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRindStates.F90,v $
! Revision 1.6  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.3  2004/12/27 23:29:50  haselbac
! Added setting of rind state for farf bc
!
! Revision 1.2  2004/10/19 19:28:31  haselbac
! Added procedure to set rind state for injecting boundaries
!
! Revision 1.1  2004/04/14 02:05:11  haselbac
! Initial revision
!
! ******************************************************************************






