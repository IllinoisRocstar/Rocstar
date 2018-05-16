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
! Purpose: Compute viscous fluxes of scalar for actual faces.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  nVarScal     Number of scalar variables
!  tvScal       Scalar transport variables
!  gradScal     Scalar face gradients
!  resScal      Scalar residuals
!
! Output: None.
!
! Notes: 
!   1. The viscosity must be a dynamic viscosity!
!
!******************************************************************************
!
! $Id: RFLU_ScalarViscousFluxes.F90,v 1.4 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ScalarViscousFluxes(pRegion,nVarScal,tvScal,gradScal,resScal)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
   
  IMPLICIT NONE

! *****************************************************************************
! Declarations
! *****************************************************************************

! =============================================================================  
! Arguments 
! =============================================================================  

  INTEGER, INTENT(IN) :: nVarScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: tvScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal 
  REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: gradScal   
  TYPE(t_region), POINTER :: pRegion

! =============================================================================  
! Locals
! =============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,ifg,iVarScal
  REAL(RFREAL) :: beta,flx,mu,nm,nx,ny,nz
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ScalarViscousFluxes.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ScalarViscousFluxes',&
  'RFLU_ScalarViscousFluxes.F90')

! *****************************************************************************
! Set variables and pointers
! *****************************************************************************

  pGrid => pRegion%grid
    
  beta = pRegion%mixtInput%betrk(pRegion%irkStep)  
  
! *****************************************************************************  
! Loop over faces and compute viscous fluxes
! *****************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! =============================================================================   
!   Get face geometry
! =============================================================================    
    
    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg) 
    nm = pGrid%fn(XYZMAG,ifg)        

! =============================================================================   
!   Get face state. NOTE this simple average is not accurate for non-uniform 
!   grids, so this will have to be replaced by a proper interpolation.
! =============================================================================

    DO iVarScal = 1,nVarScal
      mu = 0.5_RFREAL*(tvScal(iVarScal,c1) + tvScal(iVarScal,c2)) 
     
! -----------------------------------------------------------------------------
!     Compute fluxes
! -----------------------------------------------------------------------------
    
      flx = mu*(  gradScal(XCOORD,iVarScal,ifg)*nx     & 
                + gradScal(YCOORD,iVarScal,ifg)*ny     & 
                + gradScal(ZCOORD,iVarScal,ifg)*nz)*nm
  
! -----------------------------------------------------------------------------
!     Accumulate into residual     
! -----------------------------------------------------------------------------

      resScal(iVarScal,c1) = resScal(iVarScal,c1) + beta*flx    
      resScal(iVarScal,c2) = resScal(iVarScal,c2) - beta*flx
    END DO ! iVarScal
  END DO ! ifg

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ScalarViscousFluxes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ScalarViscousFluxes.F90,v $
! Revision 1.4  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2004/01/29 22:56:17  haselbac
! Initial revision
!
!******************************************************************************







