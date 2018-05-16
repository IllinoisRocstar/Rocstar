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
! Purpose: computes the RHS for the multiphase injection algorithm
!          specific to RFLU.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_InjcTileCalcRhs.F90,v 1.8 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_InjcTileCalcRhs( pRegion )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
 
  INTEGER           :: bcType,c1,distrib,errorFlag,i2dVals,iCont,iCvMass,ifc, &
                       iPatch,iTile,nCont,nPatches
                       
  INTEGER, POINTER, DIMENSION(:) :: bf2c,pCvTileMass

  REAL(RFREAL)  :: area,cp,heatCapSum,heatCapSumR,injcVelRatio,&
                   massFluxLimit,massFluxSum,massFluxSumR,minj, &
                   nm,nx,ny,nz,rhoMixtbCond,tinj,tileTemp,tileVelNrm               
  REAL(RFREAL), POINTER, DIMENSION(:)   :: injcMassFluxRatio, injcTemp, specHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv, fn, pRhs, vals
    
  TYPE(t_global), POINTER    :: global
  TYPE(t_patch), POINTER     :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag  

! TEMPORARY
  REAL(RFREAL)  :: heatCapRatio,massFluxRatio,rhoPlag
! END TEMPORARY

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_InjcTileCalcRhs.F90,v $ $Revision: 1.8 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InjcTileCalcRhs',&
  'PLAG_RFLU_InjcTileCalcRhs.F90')

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  nCont    = pRegion%plagInput%nCont
  nPatches = pRegion%grid%nPatches

  injcVelRatio = pRegion%plagInput%injcVelRatio  
  injcMassFluxRatio => pRegion%plagInput%injcMassFluxRatio
  specHeat          => pRegion%plagInput%spht
  injcTemp          => pRegion%plagInput%injcTemp

  cv => pRegion%mixt%cv

  massFluxLimit = 1.0E-10_RFREAL

! ******************************************************************************  
! Loop over patches 
! ******************************************************************************

  DO iPatch=1,nPatches

    pPatch => pRegion%patches(iPatch)
    vals   => pPatch%mixt%vals
    
    bf2c   => pPatch%bf2c
    fn     => pPatch%fn
    
    bcType = pPatch%bcType
    distrib = pPatch%mixt%distrib

!===============================================================================
!   Select injection boundary condition 
!===============================================================================

    IF ( (bcType >= BC_INJECTION .AND. bcType <= BC_INJECTION + BC_RANGE) .OR. &
         (bcType >= BC_INFLOW    .AND. bcType <= BC_INFLOW    + BC_RANGE)      ) THEN 

!-------------------------------------------------------------------------------
!     Loop over face patches
!-------------------------------------------------------------------------------

      DO ifc = 1, pPatch%nBFaces
        iTile = ifc
        c1 = bf2c(ifc)

        nx = fn(XCOORD,ifc)
        ny = fn(YCOORD,ifc)
        nz = fn(ZCOORD,ifc)
        nm = fn(XYZMAG,ifc)

        i2dVals = distrib * ifc
        mInj       = vals(BCDAT_INJECT_MFRATE,i2dVals)
        tInj       = vals(BCDAT_INJECT_TEMP  ,i2dVals)

!-------------------------------------------------------------------------------
!       Compute mixture density on boundary 
!        Using a low-order projection based on cell value
!-------------------------------------------------------------------------------
        
        rhoMixtbCond = cv(CV_MIXT_DENS,c1)

!-------------------------------------------------------------------------------
!       Update tile Rhs datastructure 
!-------------------------------------------------------------------------------

        pTilePlag => pPatch%tilePlag
 
        pRhs        => pTilePlag%rhs
        pCvTileMass => pTilePlag%cvTileMass

        SELECT CASE(bcType)

!-------------------------------------------------------------------------------
!        Injection 
!-------------------------------------------------------------------------------

         CASE(BC_INJECTION) 
          massFluxSum = SUM ( mInj * injcMassFluxRatio(:) )
          heatCapSum  = DOT_PRODUCT( specHeat, mInj * injcMassFluxRatio(:) )

          tileTemp    = tInj
          tileVelNrm  = injcVelRatio * mInj / rhoMixtbCond
          area        = nm

          DO iCont = 1, nCont
            iCvMass = pCvTileMass(iCont)
            pRhs(iCvMass,iTile) = -area *mInj *injcMassFluxRatio(iCont) 
          END DO ! iCont            
            
          pRhs(CV_TILE_MOMNRM,iTile) = -area *massFluxSum *tileVelNrm
            
          pRhs(CV_TILE_ENER  ,iTile) = -area *                             &
                                    ( 0.5_RFREAL*massFluxSum*tileVelNrm**2 &
                                    +            heatCapSum *tileTemp      )

!-------------------------------------------------------------------------------
!        Inflow based on velocity and temperature 
!-------------------------------------------------------------------------------

         CASE(BC_INFLOW_VELTEMP) 
          area        = nm

          massFluxSum = 0.0_RFREAL
          heatCapSum  = 0.0_RFREAL
          tileVelNrm  = 0.0_RFREAL
          tileTemp    = 0.0_RFREAL

          DO iCont = 1, nCont
            rhoPlag =  pRegion%plagInput%dens(iCont)
            massFluxSum = massFluxSum + injcMassFluxRatio(iCont) *rhoPlag
            heatCapSum  = heatCapSum  + specHeat(iCont) *injcMassFluxRatio(iCont) *rhoPlag 
          END DO ! iCont

! ==============================================================================
!         Set inverse of massFluxSum to avoid division by zero 
! ==============================================================================
  
          IF ( massFluxSum > massFluxLimit ) THEN
            massFluxSumR = 1.0_RFREAL/massFluxSum
            heatCapSumR       = 1.0_RFREAL/heatCapSum
          ELSE
            massFluxSumR = 1.0_RFREAL
            heatCapSumR       = 1.0_RFREAL
          END IF ! massFluxSum

          DO iCont = 1, nCont
            iCvMass = pCvTileMass(iCont)
            rhoPlag =  pRegion%plagInput%dens(iCont)

            massFluxRatio =  injcMassFluxRatio(iCont) *rhoPlag *massFluxSumR
            heatCapRatio  =  injcMassFluxRatio(iCont) *rhoPlag &
                          *  specHeat(iCont) *heatCapSumR

            tileVelNrm    = tileVelNrm +massFluxRatio *injcMassFluxRatio(iCont)
            tileTemp      = tileTemp   +heatCapRatio  *injcTemp(iCont)

            pRhs(iCvMass,iTile) = -area * injcMassFluxRatio(iCont) *rhoPlag
          END DO ! iCont

          pRhs(CV_TILE_MOMNRM,iTile) = -area *massFluxSum *tileVelNrm

          pRhs(CV_TILE_ENER  ,iTile) = -area *                              &
                                     ( 0.5_RFREAL*massFluxSum*tileVelNrm**2 &
                                     +            heatCapSum *tileTemp      )

!-------------------------------------------------------------------------------
!        Trap error for default 
!-------------------------------------------------------------------------------

         CASE DEFAULT 
           CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! bcType

      END DO ! ifc                                                        
    ENDIF    ! bcType 
  ENDDO      ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_InjcTileCalcRhs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InjcTileCalcRhs.F90,v $
! Revision 1.8  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2007/03/09 15:16:31  fnajjar
! Fixed bug for massFluxSum being zero and avoiding division by zero
!
! Revision 1.5  2006/09/18 20:35:30  fnajjar
! Activated tile datastructure for inflow bc and created proper particle bc
!
! Revision 1.4  2006/08/19 15:40:07  mparmar
! Renamed patch variables
!
! Revision 1.3  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.2  2004/12/07 22:56:24  fnajjar
! Removed rhoVrel being obsolete
!
! Revision 1.1  2004/03/08 22:36:08  fnajjar
! Initial import of RFLU-specific injection routines
!
!******************************************************************************







