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
! Notes: This routine is similar to routine ViscousFluxes for laminar flow.
!
!******************************************************************************
!
! $Id: TURB_coViscousFluxesFlo.F90,v 1.11 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_CoViscousFluxes( region ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModTurbulence, ONLY : t_turb
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensDummyNodes, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset, &
                            RFLO_CalcGradVector
  USE ModInterfaces, ONLY : RFLO_ViscousFlux
  USE TURB_ModInterfaces, ONLY : TURB_CalcStrainRate, TURB_GetModelStressCell, &
                                 TURB_GetTvCell,      TURB_LesFluxFixSmag, &
                                 TURB_LesFluxScalSim, TURB_LesGetEddyVis, &
                                 TURB_RansSAVisFlux,  TURB_RansSAGetEddyVis, &
                                 TURB_RansTotalTv,    TURB_VisFluxEddy, & 
                                 TURB_FloFaceVolume,  TURB_FloFaceWidth, &
                                 TURB_FloLesGenCoCC,  TURB_FloLesGenCoFF, &
                                 TURB_FloWlmMetric
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), TARGET :: region

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb
  
  INTEGER :: iLev, turbModel, modelClass
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, errorFlag
  INTEGER :: iCOff, ijCOff,iNOff, ijNOff, ibc, iec, ibn, ien, gradIndx(9)

  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:),  gradk(:,:)
 
!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_coViscousFluxesFlo.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_CoViscousFluxes',&
  'TURB_coViscousFluxesFlo.F90' )

! get node dimensions ------------------------------------------------------

  iLev      =  region%currLevel
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

! get pointers and parameters
 
  turbModel =  region%mixtInput%turbModel
  modelClass=  region%turbInput%modelClass
  turb      => region%levels(iLev)%turb

! get mixture strain rate tensor and store in mISij, mJSij and mKSij

  gradi  => region%levels(iLev)%mixt%gradi
  gradj  => region%levels(iLev)%mixt%gradj
  gradk  => region%levels(iLev)%mixt%gradk

  gradIndx(1) = GR_MIXT_UX
  gradIndx(2) = GR_MIXT_VX
  gradIndx(3) = GR_MIXT_WX
  gradIndx(4) = GR_MIXT_UY
  gradIndx(5) = GR_MIXT_VY
  gradIndx(6) = GR_MIXT_WY
  gradIndx(7) = GR_MIXT_UZ
  gradIndx(8) = GR_MIXT_VZ
  gradIndx(9) = GR_MIXT_WZ

  ALLOCATE( turb%mISij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
  ALLOCATE( turb%mJSij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
  ALLOCATE( turb%mKSij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  CALL TURB_CalcStrainRate( region,ibn,ien,gradIndx,gradi,gradj,gradk, &
                            turb%mISij,turb%mJSij,turb%mKSij )

! get new non-uniform filter and averaging coefficients if the grid move

  IF (region%mixtInput%moveGrid) THEN
    IF (region%irkStep == 1) THEN
      CALL TURB_LesMoveGrid
    ENDIF
  ENDIF

! allocate and initiate LES arrays required within this scope

  IF (modelClass == MODEL_LES) THEN
    ALLOCATE( turb%mueT(DIRI:DIRK,ibn:ien),stat=errorFlag )
    ALLOCATE( turb%trace(ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    turb%mueT    = 0._RFREAL
    turb%trace   = 0._RFREAL
    turb%dv      = 0._RFREAL
  ENDIF

! get total viscous flux based on selected turbulence model

  IF (turbModel==TURB_MODEL_FIXSMAG) THEN
    CALL TURB_LesFluxFixSmag( region,ibn,ien )

  ELSEIF ((turbModel==TURB_MODEL_DYNSMAG) .OR. &
          (turbModel==TURB_MODEL_DYNMIXD)) THEN

! - allocate arrays for strain rate of filtered velocities in LesGetEddyVis
    ALLOCATE( turb%fISij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
    ALLOCATE( turb%fJSij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
    ALLOCATE( turb%fKSij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )
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
    DEALLOCATE( turb%fISij,turb%fJSij,turb%fKSij )

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
    CALL RFLO_ViscousFlux( region,TVT_RANS_MUE,TVT_RANS_TCO,turb%tv )

! - deallocate retired arrays
    DEALLOCATE( turb%tv )

  ENDIF

! finalize viscous flux treatment
! if desired, interpolate model stresses to cell centers for statistics

  IF (region%turbInput%nSv > 0) CALL TURB_GetModelStressCell( region )
  DEALLOCATE( turb%mISij,turb%mJSij,turb%mKSij )

  IF (modelClass == MODEL_LES) THEN
! - interpolate transport variables at cell centers
    CALL TURB_GetTvCell( region )

! - deallocate LES arrays
    DEALLOCATE( turb%mueT,turb%trace )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! =============================================================================
!   Calls to filtering and averaging coefficients routines
! =============================================================================

CONTAINS

  SUBROUTINE TURB_LesMoveGrid

! - local variables
    TYPE(t_patch), POINTER :: patch
    INTEGER :: iPatch
    LOGICAL :: doWlm

! - compute face volumes needed for lesMij and lesCalcEddyVis
    IF ((turbModel==TURB_MODEL_FIXSMAG) .OR. &
        (turbModel==TURB_MODEL_DYNSMAG) .OR. &
        (turbModel==TURB_MODEL_DYNMIXD)) THEN
      CALL TURB_FloFaceVolume( region,DIRI )
      CALL TURB_FloFaceVolume( region,DIRJ )
      CALL TURB_FloFaceVolume( region,DIRK )
    ENDIF

! - compute filter coefficients
    IF (((turbModel==TURB_MODEL_SCALSIM)  .OR. &
         (turbModel==TURB_MODEL_DYNSMAG)  .OR. &
         (turbModel==TURB_MODEL_DYNMIXD)) .AND. &
         (region%turbInput%filterType == FILTYPE_NONUNIF)) THEN

      ALLOCATE( turb%workI(2,ibn:ien),stat=errorFlag )
      ALLOCATE( turb%workJ(2,ibn:ien),stat=errorFlag )
      ALLOCATE( turb%workK(2,ibn:ien),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL TURB_FloFaceWidth( region )
      CALL TURB_FloLesGenCoCC( region )
      CALL TURB_FloLesGenCoFF( region )
      DEALLOCATE( turb%workI,turb%workJ,turb%workK )
    ENDIF

! - if wlm is active, recompute mapping coefficients

    DO iPatch=1,region%nPatches
      patch  => region%levels(iLev)%patches(iPatch)

      doWlm = .false.
      IF (patch%bcType>=BC_NOSLIPWALL .AND. &
          patch%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN ! my boundary type
        IF (patch%valBola%switches(WLM_INPUT_MODEL) /= WLM_MODEL_NOMODEL) THEN
          doWlm = .true.
        ENDIF
      ENDIF

      IF (doWlm) THEN
! ----- compute mapping coeffs. from body fitted to cartesian and wlm-metric
        CALL TURB_FloWlmMetric( region,patch )
      ENDIF    ! doWlm
    ENDDO      ! iPatch

  END SUBROUTINE TURB_LesMoveGrid

END SUBROUTINE TURB_CoViscousFluxes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_coViscousFluxesFlo.F90,v $
! Revision 1.11  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2005/03/07 05:03:51  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.8  2004/08/04 02:47:30  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.7  2004/08/02 23:09:21  wasistho
! mv libfloflu/viscousFluxEddy(Patch) to rocflo/RFLO_viscousFlux(Patch)
!
! Revision 1.6  2004/06/03 02:10:39  wasistho
! enabled non-uniform fix-Smagorinsky
!
! Revision 1.5  2004/05/17 20:33:04  wasistho
! reordering gradIndx for more efficient cache memory addressing
!
! Revision 1.4  2004/03/19 02:54:58  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/13 03:12:50  wasistho
! get rid of flo/flu identifier in TURB_Co.. routines
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:33:31  wasistho
! changed turb nomenclature
!
! Revision 1.14  2004/02/26 21:25:07  wasistho
! delete esg1Sum
!
! Revision 1.13  2004/01/23 00:26:44  wasistho
! insert condition istage=1 for calling TURB_LesMoveGrid
!
! Revision 1.12  2004/01/22 03:53:27  wasistho
! move TURB_ransSAGetEddyVis from TURB_viscousFluxes to TURB_ransClearSendRequests
!
! Revision 1.11  2003/10/09 23:07:50  wasistho
! renamed CalcEddyVis to LesCalcEddyVis
!
! Revision 1.10  2003/10/07 02:08:02  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.9  2003/09/12 20:08:03  wasistho
!  == should be /= in TURB_LesMoveGrid
!
! Revision 1.8  2003/08/29 01:42:36  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.7  2003/05/31 01:48:23  wasistho
! installed turb. wall layer model
!
! Revision 1.6  2003/05/24 02:05:36  wasistho
! turbulence statistics expanded
!
! Revision 1.5  2003/05/16 05:43:05  wasistho
! modified array range of CC-filtered
!
! Revision 1.4  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.3  2002/10/16 07:47:58  wasistho
! Enable Fix Smagorinsky
!
! Revision 1.2  2002/10/16 02:01:28  wasistho
! Changed global%error flag
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************







