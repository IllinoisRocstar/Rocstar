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
! Purpose: Set dependent variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  icgBeg       Beginning cell index
!  icgEnd       Ending cell index
!
! Output: None.
!
! Notes: 
!   1. Restricted to perfect gas.
!
! ******************************************************************************
!
! $Id: RFLU_SetDependentVars.F90,v 1.11 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetDependentVars(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtGasLiq_P, &
                           MixtLiq_C2_Bp, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtPerf_C_GRT, &
                           MixtPerf_C2_GRT, &
                           MixtPerf_Cv_CpR, &                       
                           MixtPerf_D_PRT, &
                           MixtPerf_G_CpR, & 
                           MixtPerf_P_DEoGVm2, & 
                           MixtPerf_R_M, &
                           MixtPerf_T_CvEoVm2, & 
                           MixtPerf_T_DPR
                                 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: icgBeg,icgEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indCp,indMol
  REAL(RFREAL) :: Eo,g,gc,ir,r,rg,term,Vm2
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  INTEGER :: iSpec
  REAL(RFREAL) :: Bg2,Bl2,Bp,Bt,Bv2,Cg2,Cl2,cpg,cpm,Cv2,cvg,cvl,Cvm,cvv, &
                  gcg,gcm,gcv,gm,immg,mmg,mmm,phip,po,rl,ro,rv,rYg,rYi,rYl, &
		  rYv,to,vfg,vfl,vfv,Yg
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_type), POINTER :: pSpecType
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetDependentVars.F90,v $ $Revision: 1.11 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetDependentVars',&
  'RFLU_SetDependentVars.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pCv => pRegion%mixt%cv
  pDv => pRegion%mixt%dv 
  pGv => pRegion%mixt%gv 
  
#ifdef SPEC  
  pCvSpec => pRegion%spec%cv
#endif

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol
  
! ******************************************************************************  
! Compute dependent variables
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )
  
! ==============================================================================  
!   Incompressible fluid model
! ==============================================================================  
   
    CASE ( FLUID_MODEL_INCOMP ) 

! ==============================================================================  
!   Compressible fluid model. NOTE check state of solution vector.
! ==============================================================================  

    CASE ( FLUID_MODEL_COMP )
      SELECT CASE ( pRegion%mixtInput%gasModel ) 
      
! ------------------------------------------------------------------------------
!       Thermally and calorically perfect gas (single species)
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_TCPERF ) 
          IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%mixt%cvState    

          DO icg = icgBeg,icgEnd
            gc = MixtPerf_R_M(pGv(GV_MIXT_MOL,icg*indMol))
            g  = MixtPerf_G_CpR(pGv(GV_MIXT_CP,icg*indCp),gc)

            r  = pCv(CV_MIXT_DENS,icg)
            ir = 1.0_RFREAL/r
            Eo = ir*pCv(CV_MIXT_ENER,icg)

            Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                         pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                         pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

            pDv(DV_MIXT_PRES,icg) = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2)
            pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gc)
            pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_GRT(g,gc, &
                                                   pDv(DV_MIXT_TEMP,icg))    
          END DO ! icg

! ------------------------------------------------------------------------------
!       Mixture of thermally and calorically perfect gases
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_MIXT_TCPERF ) 
          IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%mixt%cvState    

          DO icg = icgBeg,icgEnd
            gc = MixtPerf_R_M(pGv(GV_MIXT_MOL,icg))
            g  = MixtPerf_G_CpR(pGv(GV_MIXT_CP,icg),gc)

            r  = pCv(CV_MIXT_DENS,icg)
            ir = 1.0_RFREAL/r
            Eo = ir*pCv(CV_MIXT_ENER,icg)

            Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                         pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                         pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

            pDv(DV_MIXT_PRES,icg) = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2)
            pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gc)
            pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_GRT(g,gc, &
                                                   pDv(DV_MIXT_TEMP,icg))    
          END DO ! icg

! ------------------------------------------------------------------------------
!       Pseudogas
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_MIXT_PSEUDO ) 
#ifdef SPEC        
          IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%mixt%cvState    

          DO icg = icgBeg,icgEnd
            mmm = pGv(GV_MIXT_MOL,icg)
            cpm = pGv(GV_MIXT_CP ,icg)            
            gcm = MixtPerf_R_M(mmm)
            gm  = MixtPerf_G_CpR(cpm,gcm)

            r  = pCv(CV_MIXT_DENS,icg)
            ir = 1.0_RFREAL/r
            Eo = ir*pCv(CV_MIXT_ENER,icg)

            cpg  = 0.0_RFREAL
            Yg   = 0.0_RFREAL
            phip = 0.0_RFREAL

            DO iSpec = 1,pRegion%specInput%nSpecies
              pSpecType => pRegion%specInput%specType(iSpec)

              rYi = pCvSpec(iSpec,icg)

              IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
                phip = phip + rYi/pSpecType%pMaterial%dens
              ELSE 
                cpg  = cpg  + ir*rYi*pSpecType%pMaterial%spht 
                Yg   = Yg   + ir*rYi                                                              
              END IF ! pSpecType%discreteFlag
            END DO ! iSpec

            cpg = cpg/Yg
            
            Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                         pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                         pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

            term = MixtPerf_P_DEoGVm2(r,Eo,gm,Vm2)
            pDv(DV_MIXT_PRES,icg) = cpm/cpg*term/(1.0_RFREAL-phip)
            
            term = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gcm)            
            pDv(DV_MIXT_TEMP,icg) = term*(1.0_RFREAL-phip)
            
            term = MixtPerf_C_GRT(gm,gcm,pDv(DV_MIXT_TEMP,icg))
            pDv(DV_MIXT_SOUN,icg) = term/(1.0_RFREAL-phip)
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif

! ------------------------------------------------------------------------------
!       Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_MIXT_GASLIQ ) 
#ifdef SPEC
          IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%mixt%cvState  

          ro  = global%refDensityLiq
          po  = global%refPressLiq
          to  = global%refTempLiq
          Bp  = global%refBetaPLiq
          Bt  = global%refBetaTLiq
          cvl = global%refCvLiq

          gcg = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 gcg)

          gcv = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 gcv)

          DO icg = icgBeg,icgEnd
            r   = pCv(CV_MIXT_DENS,icg)
            ir  = 1.0_RFREAL/r
            Eo  = ir*pCv(CV_MIXT_ENER,icg)
            
            rYg = r*pCvSpec(1,icg)
            rYv = r*pCvSpec(2,icg) 
            rYl = r - rYg - rYv

            Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                         pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                         pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))
            Cvm = (rYl*cvl + rYv*cvv + rYg*cvg)/r

            pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_CvEoVm2(Cvm,Eo,Vm2) 

            Cl2 = MixtLiq_C2_Bp(Bp) 
            Cv2 = MixtPerf_C2_GRT(1.0_RFREAL,gcv,pDv(DV_MIXT_TEMP,icg))
            Cg2 = MixtPerf_C2_GRT(1.0_RFREAL,gcg,pDv(DV_MIXT_TEMP,icg))
            pDv(DV_MIXT_PRES,icg) = MixtGasLiq_P(rYl,rYv,rYg,Cl2,Cv2,Cg2,r, &
                                                 ro,po,to,Bp,Bt, &
                                                 pDv(DV_MIXT_TEMP,icg))

            rl = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pDv(DV_MIXT_PRES,icg),po, &
                                        pDv(DV_MIXT_TEMP,icg),to)
            rv = MixtPerf_D_PRT(pDv(DV_MIXT_PRES,icg),gcv, &
                                pDv(DV_MIXT_TEMP,icg))
            rg = MixtPerf_D_PRT(pDv(DV_MIXT_PRES,icg),gcg, &
                                pDv(DV_MIXT_TEMP,icg))

            vfg = rYg/rg
            vfv = rYv/rv
            vfl = rYl/rl    

            Bl2 = -Bt/Bp
            Bv2 = rv*gcv
            Bg2 = rg*gcg

            pDv(DV_MIXT_SOUN,icg) = MixtGasLiq_C(Cvm,r,pDv(DV_MIXT_PRES,icg), &
                                                 rl,rv,rg,vfl,vfv,vfg,Cl2, &
                                                 Cv2,Cg2,Bl2,Bv2,Bg2)
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------     

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
      END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================  
!   Default
! ==============================================================================  

    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetDependentVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetDependentVars.F90,v $
! Revision 1.11  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/08/16 19:16:39  haselbac
! Bug fix: bad check-in, declarations missing
!
! Revision 1.8  2006/08/15 15:23:19  haselbac
! Bug fix: dv was computed incorrectly for pseudo-gas
!
! Revision 1.7  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.6  2006/03/26 20:21:40  haselbac
! Added support for GL model
!
! Revision 1.5  2005/11/14 16:55:17  haselbac
! Added support for pseudo-gas model
!
! Revision 1.4  2005/11/10 02:14:49  haselbac
! Added SELECT CASE on gasModel
!
! Revision 1.3  2005/04/15 15:06:20  haselbac
! Added range arguments, removed pGrid declaration
!
! Revision 1.2  2004/11/14 19:41:20  haselbac
! Added setting of dv for incompressible fluid model
!
! Revision 1.1  2004/10/19 19:23:53  haselbac
! Initial revision
!
! ******************************************************************************







