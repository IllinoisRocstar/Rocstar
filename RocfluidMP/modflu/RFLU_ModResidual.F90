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
! Purpose: Collection of routines related to residual computation.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModResidual.F90,v 1.18 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModResidual

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModMixture, ONLY: t_mixt,t_mixt_input
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_RES_ComputeResidual, &
            RFLU_GetResidualSupport1, & 
            RFLU_GetResidualSupport2, &  
            RFLU_RES_SetCvOld, & 
            RFLU_RES_SetDiss, & 
            RFLU_RES_SetRhs
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModResidual.F90,v $ $Revision: 1.18 $' 
              
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

SUBROUTINE RFLU_RES_ComputeResidual(pRegion)

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

  CALL RegisterFunction(global,'RFLU_RES_ComputeResidual',&
  'RFLU_ModResidual.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pMixt      => pRegion%mixt
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Set dissipation array
! ******************************************************************************

  CALL RFLU_RES_SetDiss(pRegion)

! ******************************************************************************
! Compute cell gradients. NOTE need to be done here so that back-up viscous 
! flux computation can make use of them.
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
                                       pRegion%mixt%cv,pRegion%mixt%gradCell, &
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
! Viscous fluxes
! ******************************************************************************

  IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN 
    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWT)                                

    varInfoGradFaces(CV_MIXT_XVEL) = V_MIXT_XVEL
    varInfoGradFaces(CV_MIXT_YVEL) = V_MIXT_YVEL
    varInfoGradFaces(CV_MIXT_ZVEL) = V_MIXT_ZVEL
    varInfoGradFaces(CV_MIXT_TEMP) = V_MIXT_TEMP      
    
    CALL RFLU_ComputeGradFacesWrapper(pRegion,CV_MIXT_XVEL,CV_MIXT_TEMP, & 
                                      GRF_MIXT_XVEL,GRF_MIXT_TEMP, &
                                      pRegion%mixt%cv,pRegion%mixt%gradFace)                              

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

  CALL RFLU_RES_SetRhs(pRegion)

! ******************************************************************************
! Inviscid fluxes
! ******************************************************************************

  CALL RFLU_ComputeFluxInv(pRegion,FLUX_PART_BOTH)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_RES_ComputeResidual

  









! ******************************************************************************
!
! Purpose: Determine support of residual for first-order scheme
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   icg                 Cell index 
!   rsSizeMax           Maximum allowed size of residual support
!
! Output: 
!   rs                  Residual support
!   rsSize              Actual size of residual support
!
! Notes: 
!   1. At present, only take into account cell-to-cell stencil.
!
! ******************************************************************************

SUBROUTINE RFLU_GetResidualSupport1(pRegion,icg,rs,rsSizeMax,rsSize)
  
  USE ModDataTypes
  USE ModParameters
  USE ModDataStruct, ONLY: t_region 
  USE ModGlobal, ONLY: t_global  
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icg
  INTEGER, INTENT(INOUT) :: rsSizeMax
  INTEGER, INTENT(OUT) :: rsSize
  INTEGER, DIMENSION(:) :: rs
  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
! Locals
! ==============================================================================
    
  INTEGER :: c1,c2,errorFlag,icl,ict,ifg,ifl,iPatch,nFaces   
  INTEGER, DIMENSION(:,:), POINTER :: pC2f
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GetResidualSupport1',&
  'RFLU_ModResidual.F90')

  pGrid => pRegion%grid

! ******************************************************************************
! Set up cell-to-face pointer
! ******************************************************************************

  ict = pGrid%cellGlob2Loc(1,icg)
  icl = pGrid%cellGlob2Loc(2,icg)

  SELECT CASE ( ict ) 
    CASE ( CELL_TYPE_TET )
      pC2f => pGrid%tet2f(:,:,icl)
    CASE ( CELL_TYPE_HEX ) 
      pC2f => pGrid%hex2f(:,:,icl)
    CASE ( CELL_TYPE_PRI ) 
      pC2f => pGrid%pri2f(:,:,icl)      
    CASE ( CELL_TYPE_PYR ) 
      pC2f => pGrid%pyr2f(:,:,icl)          
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
  END SELECT ! ict
  
  nFaces = SIZE(pC2f,2)

! ******************************************************************************
! Initialize
! ******************************************************************************

  rsSize = 1
  
  rs(rsSize) = icg

! ******************************************************************************
! Add face neighbors
! ******************************************************************************

  DO ifl = 1,nFaces
    iPatch = pC2f(1,ifl)
    ifg    = pC2f(2,ifl)

    IF ( iPatch == 0 ) THEN ! Interior face
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

      IF ( c1 == icg ) THEN
        IF ( rsSize < rsSizeMax ) THEN       
          rsSize = rsSize + 1
        ELSE 
! TEMPORARY
          WRITE(*,*) 'ERROR - About to exceed size of array rs!'
          STOP
! END TEMPORARY          
        END IF ! rsSize
      
        rs(rsSize) = c2
      ELSE IF ( c2 == icg ) THEN
        IF ( rsSize < rsSizeMax ) THEN       
          rsSize = rsSize + 1
        ELSE 
! TEMPORARY
          WRITE(*,*) 'ERROR - About to exceed size of array rs!'
          STOP
! END TEMPORARY          
        END IF ! rsSize

        rs(rsSize) = c1  
      ELSE ! defensive programming
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! c1                          
    END IF ! iPatch
  END DO ! ifl

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetResidualSupport1









! ******************************************************************************
!
! Purpose: Determine support of residual for second-order scheme
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   icg                 Cell index 
!   rsSizeMax           Maximum allowed size of residual support
!
! Output: 
!   rs                  Residual support
!   rsSize              Actual size of residual support
!
! Notes: 
!   1. At present, only take into account cell-to-cell stencil, i.e., only
!      inviscid flux contribution. This is because cell-to-cell stencil is 
!      typically larger than face-to-cell stenci, at least for WENO scheme.
!
! ******************************************************************************

SUBROUTINE RFLU_GetResidualSupport2(pRegion,icg,rs,rsSizeMax,rsSize)
  
  USE ModDataTypes
  USE ModParameters
  USE ModDataStruct, ONLY: t_region 
  USE ModGlobal, ONLY: t_global  
  USE ModGrid, ONLY: t_grid
  USE ModError
  
  USE ModSortSearch

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icg
  INTEGER, INTENT(INOUT) :: rsSizeMax
  INTEGER, INTENT(OUT) :: rsSize
  INTEGER, DIMENSION(:) :: rs
  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
! Locals
! ==============================================================================
    
  INTEGER :: errorFlag,icg2,icg3,iLevel,iLoc,irs,isl,nLevels,rsSaveSize, & 
             rsSaveSizeMax,rsTempSize,rsTempSizeMax  
  INTEGER, DIMENSION(:), ALLOCATABLE :: rsSave,rsTemp  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GetResidualSupport2',&
  'RFLU_ModResidual.F90')

  pGrid => pRegion%grid

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  rsSaveSizeMax = 1000
  rsTempSizeMax = 1000 ! Must be equal to or larger than <rsSaveSizeMax>

  ALLOCATE(rsSave(rsSaveSizeMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'rsSave')
  END IF ! global%errorFlag

  ALLOCATE(rsTemp(rsTempSizeMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'rsTemp')
  END IF ! global%errorFlag

! ******************************************************************************
! Initialize
! ******************************************************************************

! TEMPORARY
  nLevels = 3
! END TEMPORARY

  rsSize = 1
  
  rs(rsSize) = icg

! ******************************************************************************
! Cell-to-cell stencil for cell icg itself
! ******************************************************************************

  DO isl = 1,pGrid%c2cs(icg)%nCellMembs
    IF ( rsSize < rsSizeMax ) THEN 
      rsSize = rsSize + 1
       
      rs(rsSize) = pGrid%c2cs(icg)%cellMembs(isl)  
    ELSE 
! TEMPORARY
      WRITE(*,*) 'ERROR - About to exceed size of array rs!'
      STOP
! END TEMPORARY    
    END IF ! rsSize
  END DO ! isl

! ==============================================================================
! Save initial stencil and sort for searching 
! ==============================================================================

  rsSaveSize = rsSize  
  
  DO irs = 1,rsSaveSize
    rsSave(irs) = rs(irs)
  END DO ! irs
  
  CALL QuickSortInteger(rs(1:rsSize),rsSize)

! ******************************************************************************
! Loop over levels
! ******************************************************************************

  DO iLevel = 2,nLevels
    rsTempSize = 0  
  
! ==============================================================================
!   Loop over members in saved stencil 
! ==============================================================================  
  
    DO irs = 1,rsSaveSize
      icg2 = rsSave(irs)
      
! ------------------------------------------------------------------------------
!     Loop over stencil members
! ------------------------------------------------------------------------------
      
      DO isl = 1,pGrid%c2cs(icg2)%nCellMembs
        icg3 = pGrid%c2cs(icg2)%cellMembs(isl) 
        
! ----- Add to temporary stencil if not already in original or saved stencil ---     
        
        CALL BinarySearchInteger(rs(1:rsSize),rsSize,icg3,iLoc)
                
        IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
          IF ( rsTempSize > 1 ) THEN 
            CALL BinarySearchInteger(rsTemp(1:rsTempSize),rsTempSize,icg3,iLoc)
          ELSE 
            iLoc = ELEMENT_NOT_FOUND
          END IF ! rsTempSize
                
          IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
            IF ( rsTempSize < rsTempSizeMax ) THEN
              rsTempSize = rsTempSize + 1

              rsTemp(rsTempSize) = icg3

              IF ( rsTempSize > 1 ) THEN 
                CALL QuickSortInteger(rsTemp(1:rsTempSize),rsTempSize)
              END IF ! rsTempSize
            ELSE 
! TEMPORARY
              WRITE(*,*) 'ERROR - About to exceed size of array rsTemp!'
              STOP
! END TEMPORARY            
            END IF ! rsTempSize
          END IF ! iLoc
        END IF ! iLoc 
      END DO ! isl      
    END DO ! irs

! ==============================================================================
!   Loop over members in saved stencil, add to original stencil and sort 
! ==============================================================================  
  
    DO irs = 1,rsTempSize
      rsSave(irs) = rsTemp(irs)
    
      IF ( rsSize < rsSizeMax ) THEN 
        rsSize = rsSize + 1

        rs(rsSize) = rsTemp(irs)
      ELSE 
! TEMPORARY
        WRITE(*,*) 'ERROR - About to exceed size of array rs!'
        STOP
! END TEMPORARY       
      END IF ! rsSize
    END DO ! irs
   
    rsSaveSize = rsTempSize
  
    CALL QuickSortInteger(rs(1:rsSize),rsSize)
  END DO ! iLevel

! DEBUG
!  WRITE(*,*) rs(1:rsSize)
! END DEBUG

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(rsSave,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'rsSave')
  END IF ! global%errorFlag

  DEALLOCATE(rsTemp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'rsTemp')
  END IF ! global%errorFlag

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetResidualSupport2













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

SUBROUTINE RFLU_RES_SetCvOld(pRegion)

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

  CALL RegisterFunction(global,'RFLU_RES_SetCvOld',&
  'RFLU_ModResidual.F90')

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

END SUBROUTINE RFLU_RES_SetCvOld
  
  
  






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

SUBROUTINE RFLU_RES_SetDiss(pRegion)

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

  CALL RegisterFunction(global,'RFLU_RES_SetDiss',&
  'RFLU_ModResidual.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pMixt => pRegion%mixt

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
  END IF ! pRegion%irkStep 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_RES_SetDiss

  
  
 
  
  



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

SUBROUTINE RFLU_RES_SetRhs(pRegion)

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

  CALL RegisterFunction(global,'RFLU_RES_SetRhs',&
  'RFLU_ModResidual.F90')

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

END SUBROUTINE RFLU_RES_SetRhs



  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModResidual


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModResidual.F90,v $
! Revision 1.18  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.17  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.16  2006/08/19 15:39:16  mparmar
! Moved region%mixt%bGradFace to patch%mixt%gradFace
!
! Revision 1.15  2006/05/01 21:02:07  haselbac
! Converted to new flux routine
!
! Revision 1.14  2006/04/27 15:10:32  haselbac
! Added VENKAT limiter option
!
! Revision 1.13  2006/04/15 17:00:38  haselbac
! Added RECONST_NONE and cReconst IF for bf grads
!
! Revision 1.12  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.11  2006/04/07 14:49:42  haselbac
! Adapted to changes in f and bf diff modules and new WENO module
!
! Revision 1.10  2006/01/06 22:13:11  haselbac
! Added call to cell grad wrapper routine
!
! Revision 1.9  2005/11/14 16:59:29  haselbac
! Added support for pseudo-gas model
!
! Revision 1.8  2005/11/10 02:28:34  haselbac
! Adapted to changes in AUSM flux module
!
! Revision 1.7  2005/10/27 19:18:46  haselbac
! Adapted to changes in gradient routines
!
! Revision 1.6  2005/10/14 14:08:09  haselbac
! Added call to RFLU_EnforceHeatFlux
!
! Revision 1.5  2005/10/05 14:05:51  haselbac
! Adapted to changes in gradient computation
!
! Revision 1.4  2005/08/17 20:25:22  hdewey2
! Moved RFLU_GetResidualSupport2 from RFLU_ModPETScNewtonKrylov, added RFLU_GetResidualSupport1
!
! Revision 1.3  2005/07/14 21:44:28  haselbac
! Added limiters, new ENO scheme, and AUSM flux function
!
! Revision 1.2  2005/07/11 19:32:40  mparmar
! Adapted call to RFLU_LimitGradCells to changes in modules
!
! Revision 1.1  2005/05/19 18:17:52  haselbac
! Initial revision
!
! ******************************************************************************
  












