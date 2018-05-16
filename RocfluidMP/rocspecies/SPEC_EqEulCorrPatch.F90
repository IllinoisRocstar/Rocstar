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
!   overall mass flux for boundary faces.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to data of current region
!   pPatch         Pointer to data of current patch
!   iSpec          Index of species
!
! Output: None.
!
! Notes:
!   1. Outflow boundary conditions have a modified flux due to the Equilibrium
!      Eulerian velocity correction
!   2. Inflow  boundary conditions have unmodified fluxes: the total mass flux 
!      is specified, and has already been accounted for
!
! ******************************************************************************
!
! $Id: SPEC_EqEulCorrPatch.F90,v 1.6 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_EqEulCorrPatch(pRegion,pPatch,iSpec)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_bcvalues,t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iSpec
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_patch), POINTER :: pPatch

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: bcType,c1,ifl
  REAL(RFREAL) :: flx,gx,gy,gz,nx,ny,nz,nm,sd1x,sd1y,sd1z,sdf,tauCoef
  REAL(RFREAL), DIMENSION(:), POINTER :: pMfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pSd,pCvSpec,pTvMixt,rhsMixt,rhsSpec
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_EqEulCorrPatch.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_EqEulCorrPatch',&
  'SPEC_EqEulCorrPatch.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid   => pRegion%grid
  pMfMixt => pPatch%mfMixt

  pSd     => pRegion%mixt%sd
  pCvSpec => pRegion%spec%cv
  pTvMixt => pRegion%mixt%tv  
  rhsMixt => pRegion%mixt%rhs  
  rhsSpec => pRegion%spec%rhs

  bcType  =  pPatch%bcType

  tauCoef =  pRegion%specInput%specType(iSpec)%tauCoefficient

! ******************************************************************************
! Checks: Defensive coding, should never occur
! ******************************************************************************

  IF ( pRegion%specInput%specType(iSpec)%discreteFlag .EQV. .FALSE. ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pRegion%specInput

  IF ( pRegion%specInput%specType(iSpec)%velocityMethod /= &
       SPEC_METHV_EQEUL ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pRegion%specInput

  IF ( tauCoef <= 0.0_RFREAL ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! tauCoef

  IF ( pRegion%mixtInput%indSd /= 1 ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pRegion%mixtInput%indSd

  IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%mixt%cvState

! ******************************************************************************
! Set gravity vector for settling velocity
! ******************************************************************************

  IF ( (global%accelOn .EQV. .TRUE.) .AND. &  
       (pRegion%specInput%specType(iSpec)%settlingFlag .EQV. .TRUE.) ) THEN 
    gx = global%accelX
    gy = global%accelY
    gz = global%accelZ
  ELSE 
    gx = 0.0_RFREAL
    gy = 0.0_RFREAL
    gz = 0.0_RFREAL            
  END IF ! global%accelOn
  
! ******************************************************************************
! Select boundary type
! ******************************************************************************

  SELECT CASE ( bcType )

! ==============================================================================
!   Inflow: No additional terms
! ==============================================================================

    CASE ( BC_INFLOW_TOTANG,BC_INFLOW_VELTEMP )

! ==============================================================================
!   Outflow: No test on pMfMixt: for consistency with RFLU_ScalarFirstPatch
! ==============================================================================

    CASE ( BC_OUTFLOW )   
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        sd1x = pSd(SD_XMOM,c1)
        sd1y = pSd(SD_YMOM,c1)
        sd1z = pSd(SD_ZMOM,c1)

        sdf = ((sd1x-gx)*nx + (sd1y-gy)*ny + (sd1z-gz)*nz)/pTvMixt(TV_MIXT_MUEL,c1)
        
        flx = -tauCoef*sdf*pCvSpec(iSpec,c1)*nm

        rhsSpec(iSpec,c1) = rhsSpec(iSpec,c1) + flx
        
        rhsMixt(CV_MIXT_DENS,c1) = rhsMixt(CV_MIXT_DENS,c1) + flx          
      END DO ! ifl

! ==============================================================================
!   Slip wall: No additional terms
! ==============================================================================

    CASE ( BC_SLIPWALL )

! ==============================================================================
!   No-slip wall: No additional terms
! ==============================================================================

    CASE ( BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP )

! ==============================================================================
!   Farfield
! ==============================================================================

    CASE ( BC_FARFIELD )
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        IF ( pMfMixt(ifl) > 0.0_RFREAL ) THEN ! Outflow
          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)
          nm = pPatch%fn(XYZMAG,ifl)

          sd1x = pSd(SD_XMOM,c1)
          sd1y = pSd(SD_YMOM,c1)
          sd1z = pSd(SD_ZMOM,c1)

          sdf = ((sd1x-gx)*nx + (sd1y-gy)*ny + (sd1z-gz)*nz)/pTvMixt(TV_MIXT_MUEL,c1)

          flx = -tauCoef*sdf*pCvSpec(iSpec,c1)*nm

          rhsSpec(iSpec,c1) = rhsSpec(iSpec,c1) + flx
          
          rhsMixt(CV_MIXT_DENS,c1) = rhsMixt(CV_MIXT_DENS,c1) + flx            
        END IF ! pMfMixt
      END DO ! ifl

! ==============================================================================
!   Injection: No additional terms
! ==============================================================================

    CASE ( BC_INJECTION )
    
! ==============================================================================
!   Virtual 
! ==============================================================================

    CASE ( BC_VIRTUAL )    

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! bcType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_EqEulCorrPatch

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_EqEulCorrPatch.F90,v $
! Revision 1.6  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/11/10 02:34:41  haselbac
! Added support for gravity, cleaned up patch CASE statement
!
! Revision 1.3  2005/10/05 14:24:10  haselbac
! Bug fix: Missing definition of pTvMixt
!
! Revision 1.2  2005/03/31 17:17:41  haselbac
! Changed computation of correction, cosmetics
!
! Revision 1.1  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! ******************************************************************************







