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
! Purpose: Compute convective fluxes using centered scheme.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InvFlux_I.F90,v 1.3 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InvFlux_I(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

  USE ModInterfaces, ONLY: RFLU_CentralFirstPatch, & 
                           RFLU_GetCvLoc

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
  INTEGER :: cvMixtXVel,cvMixtYVel,cvMixtZVel,c1,c2,ifg,iPatch
  REAL(RFREAL) :: term
  REAL(RFREAL) :: flx(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pVf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pRhs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InvFlux_I.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InvFlux_I',&
  'RFLU_InvFlux_I.F90')

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  pCv  => pRegion%mixt%cv
  pRhs => pRegion%mixt%rhs
  pVf  => pRegion%mixt%vfMixt    

  cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)
  cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)
  cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

    term = 0.5_RFREAL*pVf(ifg)*pGrid%fn(XYZMAG,ifg)        

    flx(1) = pCv(cvMixtXVel,c2)*term
    flx(2) = pCv(cvMixtYVel,c2)*term
    flx(3) = pCv(cvMixtZVel,c2)*term

    pRhs(cvMixtXVel,c1) = pRhs(cvMixtXVel,c1) - flx(1)
    pRhs(cvMixtYVel,c1) = pRhs(cvMixtYVel,c1) - flx(2)
    pRhs(cvMixtZVel,c1) = pRhs(cvMixtZVel,c1) - flx(3)

    flx(1) = -pCv(cvMixtXVel,c1)*term
    flx(2) = -pCv(cvMixtYVel,c1)*term
    flx(3) = -pCv(cvMixtZVel,c1)*term

    pRhs(cvMixtXVel,c2) = pRhs(cvMixtXVel,c2) - flx(1)
    pRhs(cvMixtYVel,c2) = pRhs(cvMixtYVel,c2) - flx(2)
    pRhs(cvMixtZVel,c2) = pRhs(cvMixtZVel,c2) - flx(3)
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

! TO DO 
!  DO iPatch = 1,region%grid%nPatches
!    CALL RFLU_CentralFirstPatch(region,region%patches(iPatch))
!  END DO ! iPatch
! END TO DO 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InvFlux_I

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InvFlux_I.F90,v $
! Revision 1.3  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/19 15:41:11  haselbac
! Initial revision
!
! ******************************************************************************







