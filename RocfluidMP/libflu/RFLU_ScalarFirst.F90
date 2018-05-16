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
! Purpose: Compute first-order accurate discretization of scalar inviscid flux.
!
! Description: None.
!
! Input: 
!  pRegion              Pointer to region data
!  nVarScal             Number of scalar variables
!  cvScal               Conserved scalar variables
!  resScal              Residual due to central scalar fluxes
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ScalarFirst.F90,v 1.4 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ScalarFirst(pRegion,nVarScal,cvScal,resScal)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters
    
  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: nVarScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal
  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,ifc,iVarScal
  REAL(RFREAL) :: flx,mf,mfn,mfp,sl,sr
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ScalarFirst.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ScalarFirst',&
  'RFLU_ScalarFirst.F90')
  
! *****************************************************************************
! Checks: Defensive coding, should never occur
! *****************************************************************************
 
  IF ( pRegion%mixtInput%indMfMixt /= 1 ) THEN 
    CALL ErrorStop(global,ERR_INDMFMIXT_INVALID,__LINE__)
  END IF ! pRegion%mixtInput%indMfMixt

! *****************************************************************************
! Set dimensions and pointers 
! *****************************************************************************

  pGrid => pRegion%grid
  pMf   => pRegion%mixt%mfMixt

! *****************************************************************************
! Compute fluxes
! *****************************************************************************

  DO ifc = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifc)
    c2 = pGrid%f2c(2,ifc) 

! =============================================================================  
!   Get mass flux 
! =============================================================================  

    mf  = pMf(ifc)
    mfp = MAX(mf,0.0_RFREAL)
    mfn = MIN(mf,0.0_RFREAL)
         
! =============================================================================  
!   Compute flux and accumulate into residual
! =============================================================================  
     
    DO iVarScal = 1,nVarScal
      sl = cvScal(iVarScal,c1) 
      sr = cvScal(iVarScal,c2) 
     
      flx = mfp*sl + mfn*sr

      resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
      resScal(iVarScal,c2) = resScal(iVarScal,c2) - flx
    END DO ! iVarScal
  END DO  ! ifc
  
! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ScalarFirst

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ScalarFirst.F90,v $
! Revision 1.4  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2004/01/29 22:56:07  haselbac
! Initial revision
!
!******************************************************************************







