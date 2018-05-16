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
! Purpose: Initialize flow field in a region using hard-coded values.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowHardCode.F90,v 1.38 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowHardCode(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
    
  USE RFLU_ModBessel  
  USE RFLU_ModExactFlow  
  USE RFLU_ModFlowHardCode        
    
  USE ModInterfaces, ONLY: MixtPerf_Cv_CpR, &
                           MixtPerf_D_CGP, &
                           MixtPerf_D_PRT, &
                           MixtGasLiq_Eo_CvmTVm2, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_P_DRT, &
                           MixtPerf_R_CpG, &  
                           MixtPerf_R_M
         
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

  CHARACTER(CHRLEN) :: errorString,RCSIdentString
  INTEGER :: iBc,icg,im,in,indCp,indMol,iq
  REAL(RFREAL) :: aTot,const,cp,cpRef,Cvm,cvg,cvl,cvv,d,dInc,dMin,dOffs,dTot, &
                  dummyReal,etaqm,g,gc,gcRef,gRef,gx,gy,gz,H,height,L,Mi,mInj, &
                  mw,nx,ny,nz,omega,p,pg,pl,pMin,pOffs,pTot,pv,r,rg,rGas,ri, &
                  rl,ro,rv,rVap,t,tTot,theta,u,ur,ut,uz,v,vInj,Vm2,w,x,xMin, &
                  y,yp,yMin,z,zMin
  REAL(RFREAL) :: xx,uo_c,A_c,L_l,c,A,po,uo,xo,Rc_l,Rc,A_cl,radius,psi,r1,r2,um
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pGv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowHardCode.F90,v $ $Revision: 1.38 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowHardCode', &
                        'RFLU_InitFlowHardCode.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from hard code...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pCv        => pRegion%mixt%cv
  pGv        => pRegion%mixt%gv    
  pMixtInput => pRegion%mixtInput

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol  

! ******************************************************************************
! Set constant gas properties. NOTE used only in some cases
! ******************************************************************************

  cpRef = global%refCp
  gRef  = global%refGamma  
  gcRef = MixtPerf_R_CpG(cpRef,gRef)

! ******************************************************************************
! Initialize flow field based on user input and fluid model
! ******************************************************************************

  SELECT CASE ( pMixtInput%fluidModel )
  
! ==============================================================================
!   Incompressible fluid model
! ==============================================================================  
  
    CASE ( FLUID_MODEL_INCOMP ) 
      pRegion%mixt%cvState = CV_MIXT_STATE_PRIM
      
! TEMPORARY, to be replaced by proper code
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
! END TEMPORARY

! ==============================================================================
!   Compressible fluid model
! ==============================================================================  
    
    CASE ( FLUID_MODEL_COMP )
      pRegion%mixt%cvState = CV_MIXT_STATE_CONS

        SELECT CASE ( global%casename )

! ------------------------------------------------------------------------------
!       Cylinder shock diffraction
! ------------------------------------------------------------------------------

        CASE ( "cylds" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! x         

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Diffracting shock issuing from open-ended shock tube
! ------------------------------------------------------------------------------

! ----- Without backchamber ----------------------------------------------------

        CASE ( "ds_7p25"   , "ds_20p0"   , "ds_50p0"   , "ds_125p0"   , & 
	       "ds_7p25_v2", "ds_20p0_v2", "ds_50p0_v2", "ds_125p0_v2", &
	       "ds_7p25_v3", "ds_20p0_v3", "ds_50p0_v3", "ds_125p0_v3", &
	       "ds_7p25_v4", "ds_20p0_v4", "ds_50p0_v4", "ds_125p0_v4", &
	       "ds_7p25_v5", "ds_20p0_v5", "ds_50p0_v5", "ds_125p0_v5" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! x

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg
	  
! ----- With backchamber -----------------------------------------------------	  
	  	  
        CASE ( "ds_7p25_v6", "ds_20p0_v6", "ds_50p0_v6", "ds_125p0_v6" )
          DO icg = 1,pGrid%nCellsTot
            x =     pGrid%cofg(XCOORD,icg)
	    y = ABS(pGrid%cofg(YCOORD,icg))

            IF ( (x < pMixtInput%prepRealVal1) .AND. & 
	         (y < pMixtInput%prepRealVal2) ) THEN
              d = pMixtInput%prepRealVal3
              u = pMixtInput%prepRealVal4
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal5
            ELSE
              d = pMixtInput%prepRealVal6
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            END IF ! x

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg	  

! ------------------------------------------------------------------------------
!       hexahedral channel mesh
! ------------------------------------------------------------------------------

        CASE ( "vortexNSCBC" ) 

!         a = speed of sound
!         C = -0.0005_RFREAL*a
!         Rc = 0.15_RFREAL
!         u0 = 1.1_RFREAL*a
!         pinf = pMixtInput%prepRealVal3

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

!            radius = (x*x+y*y)**(0.5_RFREAL)
!         psi = C*exp(-(x*x+y*y)/(2.0_RFREAL*Rc*Rc))

!            IF ( radius < Rc ) THEN
!              d = pMixtInput%prepRealVal2
!              u = u0 - x*psi/(d*Rc*Rc)
!              v = y*psi/(d*Rc*Rc)
!              w = 0.0_RFREAL
!              p = pinf + d*C*psi/(Rc*Rc)
!            ELSE
!              d = pMixtInput%prepRealVal2
!              u = u0
!              v = 0.0_RFREAL 
!              w = 0.0_RFREAL
!              p = pinf
!            END IF ! radius

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Generic compressible gravity current
! ------------------------------------------------------------------------------

        CASE ( "gcgc" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( (x < pMixtInput%prepRealVal1) .AND. & 
                 (y < pMixtInput%prepRealVal2) ) THEN 
              d = pMixtInput%prepRealVal3
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            END IF ! x
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg     
	  
! ------------------------------------------------------------------------------
!       Generic multiphase wall jet
! ------------------------------------------------------------------------------

        CASE ( "gmpjet" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
	    z = pGrid%cofg(ZCOORD,icg)

            IF ( (    x  < 0.00_RFREAL) .AND. & 
                 (ABS(y) < 0.55_RFREAL) .AND. & 
		 (ABS(z) < 0.55_RFREAL) ) THEN 
              d = 34.7418238572_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 5.0E+6_RFREAL
            ELSE 
              d = 1.20903522664_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! x
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg     	  

! ------------------------------------------------------------------------------
!       Gradient test 
! ------------------------------------------------------------------------------

! ----- Linear function --------------------------------------------------------

        CASE ( "gtlin" ) 
          xMin = MINVAL(pGrid%cofg(XCOORD,1:pGrid%nCellsTot))        
          yMin = MINVAL(pGrid%cofg(YCOORD,1:pGrid%nCellsTot))
          zMin = MINVAL(pGrid%cofg(ZCOORD,1:pGrid%nCellsTot))
        
          CALL RFLU_SetExactFlowLinear(xMin,yMin,zMin,1,dMin,gx,gy,gz)
          CALL RFLU_SetExactFlowLinear(xMin,yMin,zMin,5,pMin,gx,gy,gz)
                        
          IF ( dMin < 0.0_RFREAL ) THEN 
            dOffs = -2.0_RFREAL*dMin
          ELSE 
            dOffs = 0.0_RFREAL
          END IF ! dMin  
          
          IF ( pMin < 0.0_RFREAL ) THEN 
            pOffs = -2.0_RFREAL*pMin
          ELSE 
            pOffs = 0.0_RFREAL
          END IF ! pMin                        
                                                                      
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                        
     
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            CALL RFLU_SetExactFlowLinear(x,y,z,1,d,gx,gy,gz)
            CALL RFLU_SetExactFlowLinear(x,y,z,2,u,gx,gy,gz) 
            CALL RFLU_SetExactFlowLinear(x,y,z,3,v,gx,gy,gz)
            CALL RFLU_SetExactFlowLinear(x,y,z,4,w,gx,gy,gz)
            CALL RFLU_SetExactFlowLinear(x,y,z,5,p,gx,gy,gz) 
                                                                           
            d = d + dOffs 
            p = p + pOffs                  
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ----- Trigonometric function -------------------------------------------------

        CASE ( "gttri" ) 
          nx = pRegion%mixtInput%prepRealVal1
          ny = pRegion%mixtInput%prepRealVal2
          nz = pRegion%mixtInput%prepRealVal3            
        
          dOffs = 2.0_RFREAL ! Assume amplitude is unity
          pOffs = 2.0_RFREAL ! Assume amplitude is unity
                                                                      
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                        
     
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,1,d,gx,gy,gz)
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,2,u,gx,gy,gz) 
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,3,v,gx,gy,gz)
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,4,w,gx,gy,gz)
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,5,p,gx,gy,gz)
                                                                            
            d = d + dOffs 
            p = p + pOffs                  
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Jet Flow
! ------------------------------------------------------------------------------

        CASE ( "jet" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            yp = 8.9E-05_RFREAL + 0.0175_RFREAL*(x - 2.5E-04_RFREAL)

            IF ( (y <= 8.9E-05_RFREAL) .AND. (y >= -8.9E-05_RFREAL) .AND. &
                 (x >= 0.0E+00_RFREAL) .AND. (x <=  2.5E-04_RFREAL) ) THEN

              rl = 880.0_RFREAL
              pl = 1.0_RFREAL
              rg = 11.69_RFREAL
              pg = 0.0_RFREAL
              rv = 0.722_RFREAL
              pv = 0.0_RFREAL
  
              d  = rl*pl + rg*pg + rv*pv
              H  = 8.9E-05_RFREAL
              u  = 300.0*(1.0_RFREAL - (y*y)/(H*H))

              IF ( ABS(y) < 5.5E-05_RFREAL ) THEN
                u = 300.0_RFREAL
              END IF ! ABS(y)

              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 288.0_RFREAL
            ELSE
              rl = 880.0_RFREAL
              pl = 0.0_RFREAL
              rg = 11.69_RFREAL
              pg = 1.0_RFREAL
              rv = 0.722_RFREAL
              pv = 0.0_RFREAL
          
              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 288.0_RFREAL
            END IF ! y
         
            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w
 
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Kieffer jet. NOTE configuration 1 has jet axis along z, configuration 2
!       has jet axis along x.
! ------------------------------------------------------------------------------

        CASE ( "kjet" ) 
          DO icg = 1,pGrid%nCellsTot
            z = pGrid%cofg(ZCOORD,icg)
            
            IF ( z < -0.0254_RFREAL ) THEN 
              d = 7.96_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 7.25E+5_RFREAL
            ELSE 
              d = 1.16871466381_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! z        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)             
          END DO ! icg

        CASE ( "kjet2","kjet2v3","kjet2v4","kjet2v5" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            
            IF ( x < -0.0127_RFREAL ) THEN 
              d = 7.96_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 7.25E+5_RFREAL
            ELSE 
              d = 1.1220002356_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)             
          END DO ! icg

        CASE ( "kjet2v3mp","kjet2v4mp","kjet2v5mp" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE
              d = pMixtInput%prepRealVal5
              u = pMixtInput%prepRealVal6
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            END IF ! x

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Multiphase Shocktube: Water-Air
! ------------------------------------------------------------------------------

        CASE ( "MShock_H2O_Air001" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq
          
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.7_RFREAL ) THEN
              rl = 1500.0_RFREAL
              pl = 1.0_RFREAL
              rg = 0.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 515.0_RFREAL
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 150.0_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 2.25_RFREAL 
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d 
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

        CASE ( "MShock_H2O_Air002" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq        
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.7_RFREAL ) THEN
              rl = 1500.0_RFREAL
              pl = 1.0_RFREAL
              rg = 0.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 515.0_RFREAL
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 50.0_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 6.73_RFREAL
            END IF ! x 

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Multiphase Shocktube: Air-Air-He
! ------------------------------------------------------------------------------

        CASE ( "MShock_Air_Air_He001" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.25_RFREAL ) THEN
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.7017_RFREAL 
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 98.956_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 307.13_RFREAL
            ELSE IF ( (x > 0.25_RFREAL) .AND. (x < 0.5_RFREAL) ) THEN
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.2763       
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.0_RFREAL 
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 0.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.1760_RFREAL
              pv = 1.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.5588_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Moving shock
! ------------------------------------------------------------------------------

        CASE ( "mvsh" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! z
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg     

! ------------------------------------------------------------------------------
!       NSCBC farfield test
! ------------------------------------------------------------------------------

        CASE ( "farf" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL
          uo_c = 0.0_RFREAL

          r1 = 0.006125_RFREAL
          r2 = 0.091875_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            radius = SQRT(x*x+y*y)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            um= uo + A*(SIN(global%pi*(radius-r1)/(r2-r1)))**9.0_RFREAL
            u = um*COS(ATAN2(y,x))
            p = po + (um - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = um*SIN(ATAN2(y,x))
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       NSCBC tests
! ------------------------------------------------------------------------------

        CASE ( "nscbc1" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            u = uo + A*(SIN(global%pi*x/L))**9.0_RFREAL
            p = po + (u - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = 0.0_RFREAL
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc2" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            u = uo 
            p = po + ro*c*A*((SIN(global%pi*x/L))**9.0_RFREAL)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = 0.0_RFREAL
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc3" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            y = pGrid%cofg(YCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            v = uo + A*(SIN(global%pi*y/L))**9.0_RFREAL
            p = po + (v - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            u = 0.0_RFREAL
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc4" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            y = pGrid%cofg(YCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            u = 0.0_RFREAL 
            p = po + ro*c*A*((SIN(global%pi*y/L))**9.0_RFREAL)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = uo
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc5" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            xx = ((x*x+y*y)**0.5)*COS(ATAN2(y,x)-global%pi/4.0_RFREAL)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
            
            um = uo + A*(SIN(global%pi*xx/L))**9.0_RFREAL
                   
            u = um*COS(global%pi/4.0_RFREAL)
            p = po + (um - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = um*SIN(global%pi/4.0_RFREAL)
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc6" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            xx = ((x*x+y*y)**0.5)*COS(ATAN2(y,x)-global%pi/4.0_RFREAL)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
            
            um = uo + A*(SIN(global%pi*xx/L))**9.0_RFREAL
                   
            u = uo*COS(global%pi/4.0_RFREAL)
            p = po + ro*c*A*((SIN(global%pi*xx/L))**9.0_RFREAL)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = uo*SIN(global%pi/4.0_RFREAL)
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc7" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4
          xo = pMixtInput%prepRealVal5

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
      
            x = x + xo
                     
            IF ( x < L ) THEN          
              u = uo 
              p = po + ro*c*A*((SIN(global%pi*x/L))**9.0_RFREAL)
              d = ro*(p/po)**(1.0_RFREAL/g)
              v = 0.0_RFREAL
              w = 0.0_RFREAL
            ELSE
              u = uo 
              p = po
              d = ro
              v = 0.0_RFREAL
              w = 0.0_RFREAL
            END IF ! x

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Shock-tube testing on nscbc
! ------------------------------------------------------------------------------

        CASE ( "nscbc8" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = 0.0_RFREAL 
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal3
            ELSE 
              d = pMixtInput%prepRealVal4
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal5
            END IF ! x         

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------          
!       Nozzle cavity
! ------------------------------------------------------------------------------

        CASE ( "ncavity" )   
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
             
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            
            IF ( x <= 3.0E-03_RFREAL ) THEN
              rl = 882.655_RFREAL 
              pl = 1.0_RFREAL
              rg = 24.14836_RFREAL
              pg = 0.0_RFREAL
              rv = 15.547_RFREAL 
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 293.0_RFREAL
            ELSE
              rl = 880.0_RFREAL  
              pl = 0.0_RFREAL
              rg = 24.14836_RFREAL 
              pg = 1.0_RFREAL
              rv = 15.547_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 293_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Proudman-Culick flow. NOTE this problem is two-dimensional and assumed 
!       to lie in the x-y plane, and that the injection boundary is located at 
!       y = -height.
! ------------------------------------------------------------------------------

        CASE ( "onera_c0", "onera_c0_2d_100x50" )
          CALL RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

          height = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj, &
                                               pTot,d,u,v,w,p) 

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg 
          
! ------------------------------------------------------------------------------
!       Culick flow 
! ------------------------------------------------------------------------------

        CASE ( "onera_c0_3d" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)

            CALL RFLU_ComputeExactFlowCulick(global,x,y,z,d,u,v,w,p) 

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg            

! ------------------------------------------------------------------------------
!       Pipe acoustics. NOTE the pipe is assumed to have the x-coordinate 
!       running down the axis. 
! ------------------------------------------------------------------------------

        CASE ( "pipeacoust" )
          CALL RFLU_GetParamsHardCodePAcoust(pTot,aTot)
        
          dTot = MixtPerf_D_CGP(aTot,gRef,pTot)        
                                          
          L  = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert)) 
          ro = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))        
                                                        
          im  = MAX(pMixtInput%prepIntVal1,1)
          in  = MAX(pMixtInput%prepIntVal2,1)
          iq  = MAX(pMixtInput%prepIntVal3,1)
          iBc = MAX(MIN(pMixtInput%prepIntVal4,1),0)
                  
          const = MAX(pMixtInput%prepRealVal1,0.0_RFREAL)        
                                                                 
          CALL RFLU_JYZOM(im,iq,dummyReal,etaqm,dummyReal,dummyReal)       
                                                                
          omega = aTot*SQRT((in*global%pi/L)**2 + (etaqm/ro)**2)         
        
          IF ( global%verbLevel > VERBOSE_LOW ) THEN           
            WRITE(STDOUT,'(A,5X,A,1X,I2)'   ) SOLVER_NAME, &
                  'Boundary condition:',iBc          
            WRITE(STDOUT,'(A,5X,A,3(1X,I2))') SOLVER_NAME, &
                  'Mode:',im,in,iq                   
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Total density (kg/m^3):   ',dTot
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Total pressure (N/m^2):   ',pTot            
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Angular frequency (rad/s):',omega
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Constant (-):             ',const                  
          END IF ! global%verbLevel                         
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)

            CALL RFLU_ComputeExactFlowPAcoust(global,z,y,x,global%currentTime, &
                                              L,ro,iBc,im,in,iq,etaqm,omega, &
                                              dTot,pTot,aTot,const,d,u,v,w,p) 
                                    
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,gRef,p,u,v,w)               
          END DO ! icg    
          
! ------------------------------------------------------------------------------
!       Ringleb flow. NOTE this problem is two-dimensional and assumed to lie in 
!       the x-y plane and that the exact solution is restricted to gamma = 1.4.
! ------------------------------------------------------------------------------

        CASE ( "ringleb" )
          CALL RFLU_GetParamsHardCodeRingleb(pTot,tTot)

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            CALL RFLU_ComputeExactFlowRingleb(x,y,gcRef,pTot,tTot,d,u,v,w,p)  

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,1.4_RFREAL,p,u,v,w)               
          END DO ! icg    

! ------------------------------------------------------------------------------
!       Shock-bubble interaction (Quirk and Karni 1996)
! ------------------------------------------------------------------------------

        CASE ( "ShockBubble" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( ((x-0.085_RFREAL)**2 + y**2) <= 6.25E-04_RFREAL )THEN
              rl = 1000.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.232127_RFREAL
              pv = 1.0_RFREAL

              d = rl*pl + rg*pg + rv*pv
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              t = 273.0_RFREAL
            ELSE IF ( x < 0.050_RFREAL ) THEN
              rl = 1000.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.75666_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d = rl*pl + rg*pg + rv*pv
              u = 110.49_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              t = 311.366_RFREAL
            ELSE
              rl = 1000_RFREAL
              pl = 0.0_RFREAL
              rg = 1.2763_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d = rl*pl + rg*pg + rv*pv
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              t = 273.0_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u + v*v + w*w
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCV(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg


! ------------------------------------------------------------------------------
!       Skews diffracting shock
! ------------------------------------------------------------------------------

        CASE ( "skews_ms2p0","skews_ms3p0","skews_ms4p0" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! z
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Sommerfeld shock-particle-interaction 
! ------------------------------------------------------------------------------

        CASE ( "somm_spi" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! z
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Sphere shock diffraction
! ------------------------------------------------------------------------------

        CASE ( "sphds" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! x         
                              
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)                              
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Shocktubes
! ------------------------------------------------------------------------------

! ----- Generic 1d/2d ----------------------------------------------------------  
  
        CASE ( "stg1d","stg2d" )
          DO icg = 1,pGrid%nCellsTot           
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE IF ( (x >= pMixtInput%prepRealVal1 ) .AND. &
                      (x  < pMixtInput%prepRealVal11) ) THEN  
              d = pMixtInput%prepRealVal5
              u = pMixtInput%prepRealVal6
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            ELSE
              d = pMixtInput%prepRealVal12
              u = pMixtInput%prepRealVal13
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal14
            END IF ! x         

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)	    
          END DO ! icg 

! ----- Sod case 1 -------------------------------------------------------------  
  
        CASE ( "st_sod1","st_sod1_mp2" )
          DO icg = 1,pGrid%nCellsTot           
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.5_RFREAL ) THEN 
              d = 1.0_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            ELSE 
              d = 0.125_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+4_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)	    
          END DO ! icg                       
  
! ----- Sod case 2 -------------------------------------------------------------  
  
        CASE ( "st_sod2","st_sod2_mp2" )
          DO icg = 1,pGrid%nCellsTot           
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.0_RFREAL ) THEN 
              d = 1.0_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            ELSE 
              d = 0.01_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+3_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg     
  
! ------------------------------------------------------------------------------
!       Supersonic vortex flow. NOTE this problem is two-dimensional and assumed 
!       to lie in the x-y plane. 
! ------------------------------------------------------------------------------
  
        CASE ( "ssvorth20x5l1"   ,"ssvortp20x5l1"   , &
               "ssvorth20x5l3"   ,"ssvortp20x5l3"   , &
               "ssvorth40x10l1"  ,"ssvortp40x10l1"  , & 
               "ssvorth40x10l3"  ,"ssvortp40x10l3"  , & 
               "ssvorth80x20l1"  ,"ssvortp80x20l1"  , & 
               "ssvorth80x20l3"  ,"ssvortp80x20l3"  , & 
               "ssvorth160x40l1" ,"ssvortp160x40l1" , &
               "ssvorth160x40l3" ,"ssvortp160x40l3" , &
               "ssvorth320x80l1" ,"ssvortp320x80l1" , &
               "ssvorth320x80l3" ,"ssvortp320x80l3" , &
               "ssvorth640x160l1","ssvortp640x160l1", &
               "ssvorth640x160l3","ssvortp640x160l3" )
          CALL RFLU_GetParamsHardCodeSsVortex(ri,Mi,pTot,tTot)     

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            CALL RFLU_ComputeExactFlowSsVortex(x,y,gRef,gcRef,ri,Mi,pTot, & 
                                               tTot,d,u,v,w,p)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,gRef,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Multiphase Riemann problem: Two rarefactions (Toro case 1) 
! ------------------------------------------------------------------------------

        CASE ( "Two_Rarefaction" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.5_RFREAL ) THEN
              rl = 1000_RFREAL
              pl = 0.1_RFREAL
              rg = 1.2342_RFREAL
              pg = 0.9_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = -200.0_RFREAL             
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.0_RFREAL
            ELSE
              rl = 1000_RFREAL
              pl = 0.1_RFREAL
              rg = 1.2342_RFREAL
              pg = 0.9_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 200.0_RFREAL             
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.0_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCV(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Simple volcano model
! ------------------------------------------------------------------------------

        CASE ( "volcmod2dv3" )
          DO icg = 1,pGrid%nCellsTot
            y = pGrid%cofg(YCOORD,icg)

            IF ( y < -500.0_RFREAL ) THEN
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 12.5E+6_RFREAL
              t = 600.0_RFREAL
            ELSE
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 0.87E+5_RFREAL
              t = 288.15_RFREAL
            END IF ! y

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            d = MixtPerf_D_PRT(p,gc,t)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Vortex test
! ------------------------------------------------------------------------------

        CASE ( "vort" ) 

          Rc_l = 0.15_RFREAL
          A_cl = -0.0005_RFREAL
          uo_c = 1.1_RFREAL
!          uo_c = 0.3_RFREAL

          Rc = Rc_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_cl*c*pMixtInput%prepRealVal1
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            radius = SQRT(x*x+y*y)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            IF ( radius <= Rc ) THEN
              psi = A*EXP(-(x*x+y*y)/(2.0_RFREAL*Rc*Rc))                           
            ELSE
              psi = 0.0_RFREAL
            END IF ! radius   

            d = ro 
            u = uo - y*psi/(d*Rc*Rc)
            v = x*psi/(d*Rc*Rc)
            w = 0.0_RFREAL
            p = po + (d*A*psi)/(Rc*Rc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Solid-body vortex  
! ------------------------------------------------------------------------------

        CASE ( "vortex_part" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            d = 1.0_RFREAL
            u = y
            v = -x
            w = 0.0_RFREAL
            p = 1.0E+5_RFREAL  

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Woodward-Colella ramp. NOTE the initial condition is one-dimensional 
!       and assumed to lie along the x-axis.
! ------------------------------------------------------------------------------

        CASE ( "wcramp","wcrampsc","wcrampsm" )        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < -0.1_RFREAL ) THEN 
              d = 8.0_RFREAL
              u = 2608.87906904_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 116.5E+5_RFREAL
            ELSE 
              d = 1.4_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)                   
          END DO ! icg
          
! ------------------------------------------------------------------------------
!       Gas-liquid mixture: Sod shocktube
! ------------------------------------------------------------------------------

        CASE ( "2DShock001" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.5_RFREAL ) THEN
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.0_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 348.43_RFREAL          
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 0.125_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 278.7_RFREAL 
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Default - must be due to input error
! ------------------------------------------------------------------------------

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
      END SELECT ! global%casename

! ==============================================================================
!   Default
! ==============================================================================  
    
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
  END SELECT ! pMixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from hard code done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowHardCode

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowHardCode.F90,v $
! Revision 1.38  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.37  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.36  2007/04/05 01:13:16  haselbac
! Added stg1d, modified code to allow 2nd interface
!
! Revision 1.35  2007/02/27 13:15:00  haselbac
! Added init for stg1d
!
! Revision 1.34  2007/02/16 20:01:05  haselbac
! Added init for somm_spi case
!
! Revision 1.33  2006/08/19 19:44:09  haselbac
! Significant clean-up and cosmetic changes
!
! Revision 1.32  2006/08/19 15:41:11  mparmar
! Added initializations for vortexNSCBC, farf, nscbc[1-8], vort
!
! Revision 1.31  2006/05/06 16:50:56  haselbac
! Added MP versions of kjet2 cases
!
! Revision 1.30  2006/04/20 20:57:35  haselbac
! Added kjet2v5 case
!
! Revision 1.29  2006/04/19 19:26:23  haselbac
! Fixed order of CASE statements
!
! Revision 1.28  2006/04/13 18:10:12  haselbac
! Added treatment of Culick flow
!
! Revision 1.27  2006/03/30 20:53:04  haselbac
! Changed ShockBubble hard-code, cosmetics
!
! Revision 1.26  2006/03/28 00:36:36  haselbac
! Added kjet2v4 case
!
! Revision 1.25  2006/03/26 20:22:21  haselbac
! Added cases for GL model
!
! Revision 1.24  2006/03/08 23:39:19  haselbac
! Added kjet2v3 case
!
! Revision 1.23  2006/01/06 22:16:49  haselbac
! Added routines to init linear and trig fn for grad testing
!
! Revision 1.22  2005/12/07 17:58:19  fnajjar
! Bug fix for gas constant defs in vortex_part
!
! Revision 1.21  2005/12/07 16:58:37  haselbac
! Added init for vortex_part
!
! Revision 1.20  2005/11/11 16:58:08  haselbac
! Added init for generic 2d shocktube
!
! Revision 1.19  2005/11/10 16:58:33  haselbac
! Bug fix for gmpjet
!
! Revision 1.18  2005/11/10 02:42:20  haselbac
! Added support for variable props, new cases
!
! Revision 1.17  2005/10/09 15:12:17  haselbac
! Added 2d C0 case
!
! Revision 1.16  2005/09/28 22:54:26  mparmar
! Changed cylds init to use prepRealValx
!
! Revision 1.15  2005/09/13 21:38:15  haselbac
! Changed sphds init to use prepRealValx
!
! Revision 1.14  2005/08/25 16:22:36  haselbac
! Bug fix in ds v6 init
!
! Revision 1.13  2005/08/25 03:42:48  haselbac
! Added v6 to ds cases
!
! Revision 1.12  2005/08/24 01:37:04  haselbac
! Added v5 for ds cases
!
! Revision 1.11  2005/08/21 16:03:06  haselbac
! Added v4 to ds cases
!
! Revision 1.10  2005/08/17 23:00:00  haselbac
! Added v3 versions of diffracting shock cases
!
! Revision 1.9  2005/08/08 00:12:46  haselbac
! Modified ds case initialization, now use prepRealValx
!
! Revision 1.8  2005/08/06 00:56:42  haselbac
! Added v2 to ds case names
!
! Revision 1.7  2005/07/05 19:47:58  mparmar
! Added init for diffracting shock over sphere
!
! Revision 1.6  2005/07/04 14:58:37  haselbac
! Added init for mvsh case
!
! Revision 1.5  2005/07/01 16:21:36  haselbac
! Added init for kjet2
!
! Revision 1.4  2005/06/14 00:57:45  haselbac
! Changed init for Skews case
!
! Revision 1.3  2005/04/29 12:51:24  haselbac
! Adapted to changes in interface of RFLU_ComputeExactFlowPAcoust
!
! Revision 1.2  2005/04/20 14:45:37  haselbac
! Extended init of pipeacoust case
!
! Revision 1.1  2005/04/15 15:08:16  haselbac
! Initial revision
!
! ******************************************************************************







