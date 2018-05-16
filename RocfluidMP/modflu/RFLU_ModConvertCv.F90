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
! Purpose: Suite of routines to convert state vector.
!
! Description: None.
!
! Notes: 
!   1. The routines which convert scalar conserved vectors assume that that 
!      form contains variables per unit volume (i.e., density times variable
!      per unit mass). 
!
! ******************************************************************************
!
! $Id: RFLU_ModConvertCv.F90,v 1.10 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModConvertCv

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
    
  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_ConvertCvCons2Prim, & 
            RFLU_ConvertCvPrim2Cons, & 
            RFLU_ScalarConvertCvCons2Prim, &
            RFLU_ScalarConvertCvPrim2Cons
   
  SAVE    
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


  
! ******************************************************************************
!
! Purpose: Convert conserved state vector to primitive variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ConvertCvCons2Prim(pRegion,cvStateFuture)

    USE ModInterfaces, ONLY: MixtPerf_R_M, &
                             MixtPerf_T_DPR

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: cvStateFuture
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,indMol
    REAL(RFREAL) :: gc,ir,mw,p,r
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvertCvCons2Prim',&
  'RFLU_ModConvertCv.F90')

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::ConvertCvCons2Prim")
#endif

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv
    pDv   => pRegion%mixt%dv
    pGv   => pRegion%mixt%gv

    indMol = pRegion%mixtInput%indMol

! ******************************************************************************
!   Actual conversion
! ******************************************************************************

    SELECT CASE ( pRegion%mixt%cvState )

! ==============================================================================
!     Convert from conservative to primitive form 
! ==============================================================================

      CASE ( CV_MIXT_STATE_CONS ) 
        SELECT CASE ( cvStateFuture )

! ------------------------------------------------------------------------------      
!         Convert to duvwp form        
! ------------------------------------------------------------------------------

          CASE ( CV_MIXT_STATE_DUVWP )
            pRegion%mixt%cvState = CV_MIXT_STATE_DUVWP          
          
            DO icg = 1,pGrid%nCellsTot
              ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,icg)

              pCv(CV_MIXT_XVEL,icg) = ir*pCv(CV_MIXT_XMOM,icg)
              pCv(CV_MIXT_YVEL,icg) = ir*pCv(CV_MIXT_YMOM,icg)
              pCv(CV_MIXT_ZVEL,icg) = ir*pCv(CV_MIXT_ZMOM,icg)
              pCv(CV_MIXT_PRES,icg) = pDv(DV_MIXT_PRES,icg)
            END DO ! icg

! ------------------------------------------------------------------------------          
!         Convert to duvwt form             
! ------------------------------------------------------------------------------

          CASE (CV_MIXT_STATE_DUVWT)
            pRegion%mixt%cvState = CV_MIXT_STATE_DUVWT

            SELECT CASE ( pRegion%mixtInput%fluidModel ) 

! ----------- Compressible fluid model -----------------------------------------

              CASE ( FLUID_MODEL_COMP )              
                SELECT CASE ( pRegion%mixtInput%gasModel )
                
! --------------- TC perfect gas or mixture thereof, pseudo-gas                 
                
                  CASE ( GAS_MODEL_TCPERF, &
                         GAS_MODEL_MIXT_TCPERF, & 
                         GAS_MODEL_MIXT_PSEUDO )
                    DO icg = 1,pGrid%nCellsTot                
                      r  = pCv(CV_MIXT_DENS,icg)
                      p  = pDv(DV_MIXT_PRES,icg)                
                      ir = 1.0_RFREAL/r

                      pCv(CV_MIXT_XVEL,icg) = ir*pCv(CV_MIXT_XMOM,icg)
                      pCv(CV_MIXT_YVEL,icg) = ir*pCv(CV_MIXT_YMOM,icg)
                      pCv(CV_MIXT_ZVEL,icg) = ir*pCv(CV_MIXT_ZMOM,icg)

                      mw = pGv(GV_MIXT_MOL,indMol*icg)
                      gc = MixtPerf_R_M(mw)

                      pCv(CV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,p,gc)
                    END DO ! icg 
                    
! --------------- Gas-liquid mixture                    
                    
                  CASE ( GAS_MODEL_MIXT_GASLIQ ) 
                    DO icg = 1,pGrid%nCellsTot
                      ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,icg)

                      pCv(CV_MIXT_XVEL,icg) = ir*pCv(CV_MIXT_XMOM,icg)
                      pCv(CV_MIXT_YVEL,icg) = ir*pCv(CV_MIXT_YMOM,icg)
                      pCv(CV_MIXT_ZVEL,icg) = ir*pCv(CV_MIXT_ZMOM,icg)

                      pCv(CV_MIXT_TEMP,icg) = pDv(DV_MIXT_TEMP,icg) 
                    END DO ! icg                                                   
                  CASE DEFAULT
                    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                END SELECT ! pRegion%mixtInput%gasModel              

! ----------- Default ----------------------------------------------------------

              CASE DEFAULT
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)     
            END SELECT ! pRegion%mixtInput%fluidModel 
             
! ------------------------------------------------------------------------------          
!         Default            
! ------------------------------------------------------------------------------
  
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)   
        END SELECT ! cvStateFuture

! ==============================================================================
!     Error - invalid input
! ==============================================================================

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixt%cvStateFuture

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ConvertCvCons2Prim")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ConvertCvCons2Prim




! ******************************************************************************
!
! Purpose: Convert primitive state vector to consverved variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes: 
!   1. Strictly speaking, cvStateFuture is not needed (there is only one
!      state for conserved variables), but kept for consistency with 
!      RFLU_ConvertCvCons2Prim.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ConvertCvPrim2Cons(pRegion,cvStateFuture)

    USE ModInterfaces, ONLY: MixtPerf_Cv_CpR, &
                             MixtGasLiq_Eo_CvmTVm2, &
                             MixtPerf_Eo_DGPUVW, &
                             MixtPerf_G_CpR, &
                             MixtPerf_R_M

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: cvStateFuture
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,indCp,indMol
    REAL(RFREAL) :: cp,Cvm,cvg,cvl,cvv,g,gc,mw,p,r,Rg,Rv,rYg,rYl,rYv,t,u,v, &
                    Vm2,w
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
#ifdef SPEC
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
#endif    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvertCvPrim2Cons',&
  'RFLU_ModConvertCv.F90')

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::ConvertCvPrim2Cons")
#endif

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv
    pDv   => pRegion%mixt%dv
    pGv   => pRegion%mixt%gv

#ifdef SPEC
    pCvSpec => pRegion%spec%cv
#endif

    indCp  = pRegion%mixtInput%indCp
    indMol = pRegion%mixtInput%indMol

! ******************************************************************************
!   Actual conversion
! ******************************************************************************
  
    IF ( pRegion%mixt%cvState == CV_MIXT_STATE_DUVWP .OR. & 
         pRegion%mixt%cvState == CV_MIXT_STATE_DUVWT ) THEN 

! ==============================================================================
!     Convert from primitive to conservative form 
! ==============================================================================

      SELECT CASE ( cvStateFuture )
        CASE ( CV_MIXT_STATE_CONS )
          pRegion%mixt%cvState = CV_MIXT_STATE_CONS          

          SELECT CASE ( pRegion%mixtInput%gasModel )

! ------------------------------------------------------------------------------
!           Thermally and calorically perfect gas or pseudo-gas
! ------------------------------------------------------------------------------

            CASE ( GAS_MODEL_TCPERF, & 
                   GAS_MODEL_MIXT_TCPERF, & 
                   GAS_MODEL_MIXT_PSEUDO )
              DO icg = 1,pGrid%nCellsTot
                r = pCv(CV_MIXT_DENS,icg)
                u = pCv(CV_MIXT_XVEL,icg)
                v = pCv(CV_MIXT_YVEL,icg)
                w = pCv(CV_MIXT_ZVEL,icg)
                p = pDv(DV_MIXT_PRES,icg)

                pCv(CV_MIXT_XMOM,icg) = r*u
                pCv(CV_MIXT_YMOM,icg) = r*v
                pCv(CV_MIXT_ZMOM,icg) = r*w

                cp = pGv(GV_MIXT_CP,indCp*icg)
                mw = pGv(GV_MIXT_MOL,indMol*icg)                  
                gc = MixtPerf_R_M(mw)
                g  = MixtPerf_G_CpR(cp,gc)

                pCv(CV_MIXT_ENER,icg) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)
              END DO ! icg 

! ------------------------------------------------------------------------------
!           Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------

            CASE ( GAS_MODEL_MIXT_GASLIQ )  
#ifdef SPEC                          
              Rg  = MixtPerf_R_M(pRegion%specInput%specType(1)% &
                                 pMaterial%molw)
              cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)% &
                                    pMaterial%spht,Rg)
              Rv  = MixtPerf_R_M(pRegion%specInput%specType(2)% &
                                 pMaterial%molw)
              cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)% &
                                    pMaterial%spht,Rv)
              cvl = global%refCvLiq

              DO icg = 1,pGrid%nCellsTot
                r = pCv(CV_MIXT_DENS,icg)
                u = pCv(CV_MIXT_XVEL,icg)
                v = pCv(CV_MIXT_YVEL,icg)
                w = pCv(CV_MIXT_ZVEL,icg)
                t = pDv(DV_MIXT_TEMP,icg)

                pCv(CV_MIXT_XMOM,icg) = r*u
                pCv(CV_MIXT_YMOM,icg) = r*v
                pCv(CV_MIXT_ZMOM,icg) = r*w

                rYg = pCvSpec(1,icg)
                rYv = pCvSpec(2,icg)
                rYl = r - rYg - rYv 

                Cvm = (rYl*cvl + rYg*cvg + rYv*cvv)/r
                Vm2 = u*u + v*v + w*w

                pCv(CV_MIXT_ENER,icg) = r*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
              END DO ! icg
#else
              CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                             'Can only be used with species module.')              
#endif

! ------------------------------------------------------------------------------
!           Other or invalid gas models             
! ------------------------------------------------------------------------------

            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput
          
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! cvStateFuture

! ==============================================================================
!   Error - invalid input
! ==============================================================================

    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! pRegion%mixt%cvState                   
                         
! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ConvertCvPrim2Cons")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ConvertCvPrim2Cons




! ******************************************************************************
!
! Purpose: Convert conserved state vector to primitive variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvScal              State vector of conserved variables
!   cvScalStateCurrent  Current state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ScalarConvertCvCons2Prim(pRegion,cvScal,cvScalStateCurrent)

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(INOUT) :: cvScalStateCurrent
    REAL(RFREAL), DIMENSION(:,:) :: cvScal
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,iScal,nScal
    REAL(RFREAL) :: ir
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ScalarConvertCvCons2Prim',&
  'RFLU_ModConvertCv.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv

    nScal = SIZE(cvScal,1)

! ******************************************************************************
!   Actual conversion
! ******************************************************************************

    SELECT CASE ( cvScalStateCurrent )

! ==============================================================================
!     Conserved, so convert to primitive variables     
! ==============================================================================

      CASE ( CV_MIXT_STATE_CONS )
        cvScalStateCurrent = CV_MIXT_STATE_PRIM         

        DO icg = 1,pGrid%nCellsTot
          ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,icg)

          DO iScal = 1,nScal
            cvScal(iScal,icg) = ir*cvScal(iScal,icg)
          END DO ! iScal
        END DO ! icg

! ==============================================================================
!     Error - invalid input
! ==============================================================================

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! cvScalStateCurrent

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ScalarConvertCvCons2Prim




! ******************************************************************************
!
! Purpose: Convert primitive state vector to conserved variables.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   cvScal              State vector of conserved variables
!   cvScalStateCurrent  Current state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ScalarConvertCvPrim2Cons(pRegion,cvScal,cvScalStateCurrent)

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(INOUT) :: cvScalStateCurrent
    REAL(RFREAL), DIMENSION(:,:) :: cvScal
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,iScal,nScal
    REAL(RFREAL) :: r
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ScalarConvertCvPrim2Cons',&
  'RFLU_ModConvertCv.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    pCv   => pRegion%mixt%cv

    nScal = SIZE(cvScal,1)

! ******************************************************************************
!   Actual conversion
! ******************************************************************************

    SELECT CASE ( cvScalStateCurrent )

! ==============================================================================
!     Primitive, so convert to conserved variables     
! ==============================================================================

      CASE ( CV_MIXT_STATE_PRIM )
        cvScalStateCurrent = CV_MIXT_STATE_CONS         

        DO icg = 1,pGrid%nCellsTot
          r = pCv(CV_MIXT_DENS,icg)

          DO iScal = 1,nScal
            cvScal(iScal,icg) = r*cvScal(iScal,icg)
          END DO ! iScal
        END DO ! icg

! ==============================================================================
!     Error - invalid input
! ==============================================================================

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! cvScalStateCurrent

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ScalarConvertCvPrim2Cons





! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModConvertCv


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModConvertCv.F90,v $
! Revision 1.10  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.7  2006/03/26 20:22:00  haselbac
! Added support for GL model, cosmetics
!
! Revision 1.6  2005/11/14 16:59:08  haselbac
! Added support for pseudo-gas model
!
! Revision 1.5  2005/11/10 02:26:54  haselbac
! Cosmetics only
!
! Revision 1.4  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.3  2005/07/07 22:45:01  haselbac
! Added profiling calls
!              
! Revision 1.2  2004/01/29 22:57:35  haselbac  
! Added routines for scalars                   
!
! Revision 1.1  2002/09/09 15:17:58  haselbac  
! Initial revision                             
!
! ******************************************************************************










