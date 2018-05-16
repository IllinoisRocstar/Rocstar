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
! Purpose: Collection of routines related to explicit multistage schemes.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModExplicitMultiStage.F90,v 1.19 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModExplicitMultiStage

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt,t_mixt_input
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_EMS_ComputeResidual, & 
            RFLU_EMS_SetCvOld, & 
            RFLU_EMS_SetDiss, & 
            RFLU_EMS_SetRhs, & 
            RFLU_EMS_UpdateConservedVars
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModExplicitMultiStage.F90,v $ $Revision: 1.19 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  




! ******************************************************************************
!
! Purpose: Compute residual.
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

SUBROUTINE RFLU_EMS_ComputeResidual(pRegion)

  USE RFLU_ModConvertCv
  USE RFLU_ModDifferentiationBFaces
  USE RFLU_ModDifferentiationCells
  USE RFLU_ModDifferentiationFaces  
  USE RFLU_ModLimiters
  USE RFLU_ModViscousFlux
  USE RFLU_ModWENO     

  USE ModInterfaces, ONLY: RFLU_ComputeFluxInv

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

  INTEGER :: iPatch
  INTEGER :: varInfoGradCells(CV_MIXT_DENS:CV_MIXT_PRES), & 
             varInfoGradFaces(CV_MIXT_XVEL:CV_MIXT_TEMP)
  TYPE(t_global), POINTER :: global
  TYPE(t_mixt), POINTER :: pMixt
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EMS_ComputeResidual',&
  'RFLU_ModExplicitMultiStage.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pMixt      => pRegion%mixt
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Set dissipation array
! ******************************************************************************

  CALL RFLU_EMS_SetDiss(pRegion)

! ******************************************************************************
! Compute cell gradients
! ******************************************************************************

  IF ( pMixtInput%spaceOrder > DISCR_ORDER_1 ) THEN
    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

    varInfoGradCells(CV_MIXT_DENS) = V_MIXT_DENS
    varInfoGradCells(CV_MIXT_XVEL) = V_MIXT_XVEL
    varInfoGradCells(CV_MIXT_YVEL) = V_MIXT_YVEL
    varInfoGradCells(CV_MIXT_ZVEL) = V_MIXT_ZVEL
    varInfoGradCells(CV_MIXT_PRES) = V_MIXT_PRES
                          
    CALL RFLU_ComputeGradCellsWrapper(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                      GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                      varInfoGradCells,pMixt%cv,pMixt%gradCell)

    SELECT CASE ( pRegion%mixtInput%reconst ) 
      CASE ( RECONST_NONE )
      CASE ( RECONST_WENO_SIMPLE ) 
        CALL RFLU_WENOGradCellsWrapper(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                       pRegion%mixt%gradCell)
        CALL RFLU_LimitGradCellsSimple(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                      GRC_MIXT_DENS,GRC_MIXT_PRES, & 
                                      pRegion%mixt%cv,pRegion%mixt%cvInfo, &
                                      pRegion%mixt%gradCell)
      CASE ( RECONST_WENO_XYZ ) 
        CALL RFLU_WENOGradCellsXYZWrapper(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                          pRegion%mixt%gradCell)
        CALL RFLU_LimitGradCellsSimple(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                      GRC_MIXT_DENS,GRC_MIXT_PRES, & 
                                      pRegion%mixt%cv,pRegion%mixt%cvInfo, &
                                      pRegion%mixt%gradCell)
      CASE ( RECONST_LIM_BARTHJESP )  
        CALL RFLU_CreateLimiter(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                pRegion%mixt%lim)
        CALL RFLU_ComputeLimiterBarthJesp(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                          GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                          pRegion%mixt%cv, &
                                          pRegion%mixt%gradCell, &
                                          pRegion%mixt%lim)
        CALL RFLU_LimitGradCells(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                 pRegion%mixt%gradCell,pRegion%mixt%lim)
        CALL RFLU_DestroyLimiter(pRegion,pRegion%mixt%lim)
      CASE ( RECONST_LIM_VENKAT )  
        CALL RFLU_CreateLimiter(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                pRegion%mixt%lim)
        CALL RFLU_ComputeLimiterVenkat(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                       GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                       pRegion%mixt%cv, &
                                       pRegion%mixt%gradCell, &
                                       pRegion%mixt%lim)
        CALL RFLU_LimitGradCells(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                 pRegion%mixt%gradCell,pRegion%mixt%lim)
        CALL RFLU_DestroyLimiter(pRegion,pRegion%mixt%lim)
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%reconst

    CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)
  END IF ! pMixtInput%spaceOrder

! ******************************************************************************
! Convective fluxes (dissipative part)
! ******************************************************************************

  IF ( pMixtInput%ldiss(pRegion%irkStep) /=0 ) THEN 
    CALL RFLU_ComputeFluxInv(pRegion,FLUX_PART_DISSIP) 
  END IF ! pMixtInput%ldiss

! ******************************************************************************
! Viscous fluxes
! ******************************************************************************

  IF ( (pMixtInput%flowModel == FLOW_NAVST) .AND. & 
       (pMixtInput%ldiss(pRegion%irkStep) /=0) ) THEN 
    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWT)  
    
    varInfoGradFaces(CV_MIXT_XVEL) = V_MIXT_XVEL
    varInfoGradFaces(CV_MIXT_YVEL) = V_MIXT_YVEL
    varInfoGradFaces(CV_MIXT_ZVEL) = V_MIXT_ZVEL
    varInfoGradFaces(CV_MIXT_TEMP) = V_MIXT_TEMP      
    
    CALL RFLU_ComputeGradFacesWrapper(pRegion,CV_MIXT_XVEL,CV_MIXT_TEMP, & 
                                      GRF_MIXT_XVEL,GRF_MIXT_TEMP, &
                                      pRegion%mixt%cv, &
                                      pRegion%mixt%gradFace)                              

    IF ( pRegion%grid%nFacesConstr > 0 ) THEN                                  
      CALL RFLU_ComputeGradFacesConstr(pRegion,CV_MIXT_XVEL,CV_MIXT_TEMP, & 
                                       GRF_MIXT_XVEL,GRF_MIXT_TEMP, &
                                       varInfoGradFaces,pRegion%mixt%cv, &
                                       pRegion%mixt%gradFace) 
    END IF ! pRegion%grid%nFacesConstr                                
                                                                                             
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        CALL RFLU_ComputeGradBFacesWrapper(pRegion,pPatch,CV_MIXT_XVEL, &
                                           CV_MIXT_TEMP,GRBF_MIXT_XVEL, &
                                           GRBF_MIXT_TEMP,pRegion%mixt%cv, &
                                           pPatch%mixt%gradFace)
        IF ( pPatch%cReconst /= CONSTR_NONE ) THEN                                          
          CALL RFLU_ComputeBFGradConstrWrapper(pRegion,pPatch,CV_MIXT_XVEL, &
                                               CV_MIXT_TEMP,GRBF_MIXT_XVEL, & 
                                               GRBF_MIXT_TEMP,varInfoGradFaces, &
                                               pRegion%mixt%cv, &
                                               pPatch%mixt%gradFace) 
        END IF ! pPatch%cReconst                                             
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
                                      
    IF ( pMixtInput%turbModel == TURB_MODEL_NONE ) THEN ! Only laminar                                    
      CALL RFLU_EnforceHeatFlux(pRegion,pRegion%mixt%tv,TV_MIXT_TCOL)
      CALL RFLU_ViscousFluxes(pRegion,pRegion%mixt%tv,TV_MIXT_MUEL,TV_MIXT_TCOL)
      CALL RFLU_ViscousFluxesPatches(pRegion,pRegion%mixt%tv,TV_MIXT_MUEL, &
                                     TV_MIXT_TCOL)
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! pMixtInput%turbModel                                                          

    CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)
  END IF ! pMixtInput%flowModel

! ******************************************************************************
! Set right-hand side array
! ******************************************************************************

  CALL RFLU_EMS_SetRhs(pRegion)

! ******************************************************************************
! Convective fluxes (central part)
! ******************************************************************************

  CALL RFLU_ComputeFluxInv(pRegion,FLUX_PART_CENTRAL)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EMS_ComputeResidual

  









! ******************************************************************************
!
! Purpose: Set old values of conserved variables.
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

SUBROUTINE RFLU_EMS_SetCvOld(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE (t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt), POINTER :: pMixt

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EMS_SetCvOld',&
  'RFLU_ModExplicitMultiStage.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pMixt => pRegion%mixt

! ******************************************************************************
! Set old values
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    pMixt%cvOld(CV_MIXT_DENS,icg) = pMixt%cv(CV_MIXT_DENS,icg)
    pMixt%cvOld(CV_MIXT_XMOM,icg) = pMixt%cv(CV_MIXT_XMOM,icg)
    pMixt%cvOld(CV_MIXT_YMOM,icg) = pMixt%cv(CV_MIXT_YMOM,icg)
    pMixt%cvOld(CV_MIXT_ZMOM,icg) = pMixt%cv(CV_MIXT_ZMOM,icg)
    pMixt%cvOld(CV_MIXT_ENER,icg) = pMixt%cv(CV_MIXT_ENER,icg)
  END DO ! icg  

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EMS_SetCvOld
  
  
  






! ******************************************************************************
!
! Purpose: Set dissipation array.
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

SUBROUTINE RFLU_EMS_SetDiss(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE (t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg
  REAL(RFREAL) :: term
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt), POINTER :: pMixt
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EMS_SetDiss',&
  'RFLU_ModExplicitMultiStage.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pMixt      => pRegion%mixt
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Set old values
! ******************************************************************************

  IF ( pRegion%irkStep == 1 ) THEN 
    DO icg = 1,pGrid%nCellsTot
      pMixt%diss(CV_MIXT_DENS,icg) = 0.0_RFREAL
      pMixt%diss(CV_MIXT_XMOM,icg) = 0.0_RFREAL
      pMixt%diss(CV_MIXT_YMOM,icg) = 0.0_RFREAL
      pMixt%diss(CV_MIXT_ZMOM,icg) = 0.0_RFREAL
      pMixt%diss(CV_MIXT_ENER,icg) = 0.0_RFREAL
    END DO ! icg 
  ELSE 
    IF ( pMixtInput%ldiss(pRegion%irkStep) /=0 ) THEN
      term = 1.0_RFREAL - pMixtInput%betrk(pRegion%irkStep)

      DO icg = 1,pGrid%nCellsTot
        pMixt%diss(CV_MIXT_DENS,icg) = term*pMixt%diss(CV_MIXT_DENS,icg)
        pMixt%diss(CV_MIXT_XMOM,icg) = term*pMixt%diss(CV_MIXT_XMOM,icg)
        pMixt%diss(CV_MIXT_YMOM,icg) = term*pMixt%diss(CV_MIXT_YMOM,icg)
        pMixt%diss(CV_MIXT_ZMOM,icg) = term*pMixt%diss(CV_MIXT_ZMOM,icg)
        pMixt%diss(CV_MIXT_ENER,icg) = term*pMixt%diss(CV_MIXT_ENER,icg)
      END DO ! icg
    END IF ! pMixtInput%ldiss    
  END IF ! pRegion%irkStep 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EMS_SetDiss

  
  
 
  
  



! ******************************************************************************
!
! Purpose: Set right-hand side.
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

SUBROUTINE RFLU_EMS_SetRhs(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE (t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt), POINTER :: pMixt

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EMS_SetRhs',&
  'RFLU_ModExplicitMultiStage.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pMixt => pRegion%mixt

! ******************************************************************************
! Set right-hand side
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    pMixt%rhs(CV_MIXT_DENS,icg) = -pMixt%diss(CV_MIXT_DENS,icg)
    pMixt%rhs(CV_MIXT_XMOM,icg) = -pMixt%diss(CV_MIXT_XMOM,icg)
    pMixt%rhs(CV_MIXT_YMOM,icg) = -pMixt%diss(CV_MIXT_YMOM,icg)
    pMixt%rhs(CV_MIXT_ZMOM,icg) = -pMixt%diss(CV_MIXT_ZMOM,icg)
    pMixt%rhs(CV_MIXT_ENER,icg) = -pMixt%diss(CV_MIXT_ENER,icg)
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EMS_SetRhs








! ******************************************************************************
!
! Purpose: Update conserved variables.
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

SUBROUTINE RFLU_EMS_UpdateConservedVars(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE (t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg
  REAL(RFREAL) :: term
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt), POINTER :: pMixt
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EMS_UpdateConservedVars',&
  'RFLU_ModExplicitMultiStage.F90')

! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pGrid      => pRegion%grid
  pMixt      => pRegion%mixt
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Set old values
! ******************************************************************************

  DO icg = 1,pGrid%nCells
    term = pMixtInput%ark(pRegion%irkStep)*pMixtInput%cfl*pRegion%dt(icg)/pGrid%vol(icg)        

    pMixt%cv(CV_MIXT_DENS,icg) = pMixt%cvOld(CV_MIXT_DENS,icg) - term*pMixt%rhs(CV_MIXT_DENS,icg)
    pMixt%cv(CV_MIXT_XMOM,icg) = pMixt%cvOld(CV_MIXT_XMOM,icg) - term*pMixt%rhs(CV_MIXT_XMOM,icg)
    pMixt%cv(CV_MIXT_YMOM,icg) = pMixt%cvOld(CV_MIXT_YMOM,icg) - term*pMixt%rhs(CV_MIXT_YMOM,icg)
    pMixt%cv(CV_MIXT_ZMOM,icg) = pMixt%cvOld(CV_MIXT_ZMOM,icg) - term*pMixt%rhs(CV_MIXT_ZMOM,icg)
    pMixt%cv(CV_MIXT_ENER,icg) = pMixt%cvOld(CV_MIXT_ENER,icg) - term*pMixt%rhs(CV_MIXT_ENER,icg)
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EMS_UpdateConservedVars




  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModExplicitMultiStage


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModExplicitMultiStage.F90,v $
! Revision 1.19  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.18  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.17  2006/08/19 15:39:07  mparmar
! Moved region%mixt%bGradFace to patch%mixt%gradFace
!
! Revision 1.16  2006/05/01 21:01:19  haselbac
! Converted to new flux routine
!
! Revision 1.15  2006/04/27 15:10:32  haselbac
! Added VENKAT limiter option
!
! Revision 1.14  2006/04/15 16:59:24  haselbac
! Added RECONST_NONE and cReconst IF for bf grads
!
! Revision 1.13  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.12  2006/04/07 14:49:41  haselbac
! Adapted to changes in f and bf diff modules and new WENO module
!
! Revision 1.11  2006/01/06 22:12:44  haselbac
! Added call to cell grad wrapper routine
!
! Revision 1.10  2005/11/14 16:59:29  haselbac
! Added support for pseudo-gas model
!
! Revision 1.9  2005/11/10 02:28:49  haselbac
! Adapted to changes in AUSM flux module
!
! Revision 1.8  2005/10/27 19:18:35  haselbac
! Adapted to changes in gradient routines
!
! Revision 1.7  2005/10/14 14:06:07  haselbac
! Added call to RFLU_EnforceHeatFlux
!
! Revision 1.6  2005/10/05 13:57:06  haselbac
! Adapted to changes in gradient computation
!
! Revision 1.5  2005/07/14 21:43:57  haselbac
! Added limiters, new ENO scheme, and AUSM flux function
!
! Revision 1.4  2005/07/11 19:32:26  mparmar
! Adapted call to RFLU_LimitGradCells to changes in modules
!
! Revision 1.3  2005/07/07 22:45:01  haselbac
! Added profiling calls
!
! Revision 1.2  2005/05/19 18:18:44  haselbac
! Cosmetics only
!
! Revision 1.1  2005/05/16 20:36:29  haselbac
! Initial revision
!
! ******************************************************************************
  











