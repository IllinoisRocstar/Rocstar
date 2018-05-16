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
!******************************************************************************
!
! Purpose: Compute total viscous fluxes and add them to dissipation.
!
! Description: total viscous flux is superposition of laminar and turbulent:
!              fv_total = sigma_ij + tau_ij
!                       = mu_l.S_ij + tau_ij
!              for eddy viscosity type turbulence model:
!              fv_total = (mu_l+mu_t).S_ij, otherwise
!              fv_total = mu_l.S_ij + m_ij, where m_ij is a model for tau_ij
!
! Input: region = data of current region.
!
! Output: mixt%diss = total viscous fluxes + num. dissipation
!
! Notes: 
!   1. This routine is similar to routine ViscousFluxes for laminar flow.
!   2. THIS ROUTINE IS BROKEN FOR ROCFLU - CANNOT ADAPT EASILY TO CHANGES
!      IN ROCFLU, WILL NOT BE USED IN ROCFLU.
!
!******************************************************************************
!
! $Id: TURB_coViscousFluxesFlu.F90,v 1.9 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_coViscousFluxes( region ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModTurbulence, ONLY : t_turb
  USE ModGlobal, ONLY     : t_global
  USE RFLU_ModViscousFlux, ONLY: RFLU_ViscousFluxes, RFLU_ViscousFluxesPatches  
  USE RFLU_ModConvertCv, ONLY  : RFLU_ConvertCvCons2Prim, & 
                                 RFLU_ConvertCvPrim2Cons  
  USE TURB_ModInterfaces, ONLY : TURB_CalcStrainRate, TURB_GetModelStressCell, &
                                 TURB_GetTvCell,      TURB_LesFluxFixSmag, &
                                 TURB_LesFluxScalSim, TURB_LesGetEddyVis, &
                                 TURB_RansSAVisFlux,  TURB_RansSAGetEddyVis, &
                                 TURB_RansTotalTv,    TURB_VisFluxEddy
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: region

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb
  TYPE(t_patch), POINTER  :: patch
 
  INTEGER :: turbModel, modelClass, errorFlag, prevCvState
  INTEGER :: nPatches, nCellsTot, nFaces, nFacesTot, nBFaces, nBFacesTot
  INTEGER :: ibc, iec, ibn, ien, gradIndx(3)

  REAL(RFREAL), POINTER :: gradi(:,:,:), bGradi(:,:,:)
 
!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_coViscousFluxesFlu.F90,v $ $Revision: 1.9 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_CoViscousFluxes',&
  'TURB_coViscousFluxesFlu.F90' )

! Specific Rocflu ------------------------------------------------------------ 
! check the state of cv first and convert to conservative if not yet

  IF (region%mixt%cvState /= CV_MIXT_STATE_CONS) THEN
    prevCvState = region%mixt%cvState
    CALL RFLU_ConvertCvPrim2Cons(region,CV_MIXT_STATE_CONS)
  ENDIF

! get cell and node dimensions -----------------------------------------------
  nCellsTot = region%grid%nCellsTot
  nFaces    = region%grid%nFaces
  nFacesTot = region%grid%nFacesTot
  ibc       = 1
  iec       = nCellsTot
  ibn       = 1
  ien       = nFaces
  nPatches  = region%grid%nPatches

  nBFaces    = 0
  nBFacesTot = 0

  DO iPatch = 1,nPatches
    patch => region%patches(iPatch)

    nBFaces    = nBFaces    + patch%nBTris    + patch%nBQuads
    nBFacesTot = nBFacesTot + patch%nBTrisTot + patch%nBQuadsTot
  END DO ! iPatch

! get pointers and parameters
 
  turbModel =  region%mixtInput%turbModel
  modelClass=  region%turbInput%modelClass
  turb      => region%turb

! get mixture strain rate tensor and store in mISij

  gradi  => region%mixt%gradFace
! TEMPORARY 
!  bGradi => region%mixt%bGradFace

  gradIndx(1) = GRF_MIXT_XVEL
  gradIndx(2) = GRF_MIXT_YVEL
  gradIndx(3) = GRF_MIXT_ZVEL

  ALLOCATE( turb%mISij(TENSOR_SYMM_NELM, nFaces ),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( turb%bmISij(TENSOR_SYMM_NELM,nBFaces),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  CALL TURB_CalcStrainRate( region,1, nFaces,gradIndx, gradi, turb%mISij )
  IF (nPatches > 0) &
  CALL TURB_CalcStrainRate( region,1,nBFaces,gradIndx,bGradi,turb%bmISij )

! get new non-uniform filter and averaging coefficients if the grid moves

!  IF (region%mixtInput%moveGrid) THEN  ! better performed in TURB_CalcMetric
!    IF (region%irkStep == 1) THEN
!      CALL TURB_FluLesMoveGrid
!    ENDIF
!  ENDIF

! allocate and initiate LES arrays required within this scope

  IF (modelClass == MODEL_LES) THEN
    ALLOCATE( turb%mueT( DIRI, nFaces),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( turb%bMueT(DIRI,nBFaces),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( turb%trace(ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    turb%dv      = 0._RFREAL  ! cell dynamic coefficient is stored in dv
    turb%mueT    = 0._RFREAL
    turb%bMueT   = 0._RFREAL
    turb%trace   = 0._RFREAL
  ENDIF

! get total viscous flux based on selected turbulence model

  IF (turbModel==TURB_MODEL_FIXSMAG) THEN
    CALL TURB_LesFluxFixSmag( region,ibn,ien )

  ELSEIF ((turbModel==TURB_MODEL_DYNSMAG) .OR. &
          (turbModel==TURB_MODEL_DYNMIXD)) THEN

! - allocate arrays for strain rate of filtered velocities in LesGetEddyVis,
!   and for space of tauij (in Dynamic Mixed model)

    ALLOCATE( turb%fISij( TENSOR_SYMM_NELM, nFaces),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( turb%bfISij(TENSOR_SYMM_NELM,nBFaces),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - get eddy viscosity of dynamic Smagorinsky and dynamic Mixed models 
    CALL TURB_LesGetEddyVis( region,ibc,iec,ibn,ien )

! - get viscous fluxes
    IF (turbModel==TURB_MODEL_DYNSMAG) THEN
      CALL TURB_VisFluxEddy( region )
    ELSEIF (turbModel==TURB_MODEL_DYNMIXD) THEN
      CALL TURB_VFluxHybrid( region )
    ENDIF

! - deallocate retired arrays
    DEALLOCATE( turb%fISij, turb%bfISij )

  ELSEIF (turbModel==TURB_MODEL_SCALSIM) THEN
    CALL TURB_LesFluxScalSim( region,ibn,ien )


  ELSEIF ((turbModel==TURB_MODEL_SA).OR. &
          (turbModel==TURB_MODEL_DESSA).OR. &
          (turbModel==TURB_MODEL_HDESSA)) THEN

! - viscous fluxes of SA equation
    CALL TURB_RansSAVisFlux( region )

! - allocate total tv for NS viscous fluxes
    ALLOCATE( turb%tv(TVT_RANS_NELM,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - compute total RaNS tv and use it to obtain NS viscous fluxes 
    CALL TURB_RansTotalTv( region,TVT_RANS_MUE,TVT_RANS_TCO,turb%tv )
    CALL RFLU_ViscousFluxes( region,turb%tv,TVT_RANS_MUE,TVT_RANS_TCO )
    CALL RFLU_ViscousFluxesPatches( region,turb%tv,TVT_RANS_MUE,TVT_RANS_TCO )

! - deallocate retired arrays
    DEALLOCATE( turb%tv )

  ENDIF

! finalize viscous flux treatment
! if desired, interpolate model stresses to cell centers for statistics

  IF (region%turbInput%nSv > 0) CALL TURB_GetModelStressCell( region )
  DEALLOCATE( turb%mISij, turb%bmISij )

  IF (modelClass == MODEL_LES) THEN
! - interpolate transport variables at cell centers
    CALL TURB_GetTvCell( region )

! - deallocate LES arrays
    DEALLOCATE( turb%mueT, turb%bMueT, turb%trace )
  ENDIF

! convert cv back to the previous state before entering this routine

  IF (region%mixt%cvState /= prevCvState) &
      CALL RFLU_ConvertCvCons2Prim( region,prevCvState )  

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CoViscousFluxes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_coViscousFluxesFlu.F90,v $
! Revision 1.9  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/08/19 15:40:43  mparmar
! Commented use of region%mixt%bGradFace
!
! Revision 1.6  2005/12/29 19:52:06  wasistho
! modified allocation for bMISij
!
! Revision 1.5  2005/12/20 20:44:19  wasistho
! adapted to changing in Rocflu on viscous fluxes routines
!
! Revision 1.4  2005/03/07 05:03:58  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.3  2004/05/28 01:58:35  wasistho
! update unstructured grid LES
!
! Revision 1.2  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.1  2004/03/25 04:42:58  wasistho
! prepared for RFLU
!
!
!
!******************************************************************************







