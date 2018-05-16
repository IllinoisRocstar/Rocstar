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
! Purpose: Collect utility routines for extraction of data from flow solution.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModExtractFlowDataUtils.F90,v 1.4 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModExtractFlowDataUtils

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region

#ifdef PLAG
  USE PLAG_ModParameters
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModExtractFlowDataUtils.F90,v $ $Revision: 1.4 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ExtractShockLocation1D

! ==============================================================================
! Private functions
! ==============================================================================



! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS




! ******************************************************************************
!
! Purpose: Extract shock location for 1D cases.
!
! Description: Look for location of largest density gradient and then, 
!   starting from this position, determine position of zero second derivative, 
!   which is then taken to be shock location.
!
! Input:
!   pRegion		Pointer to region
!   icgBeg		Beginning cell index
!   icgEnd		Ending cell index
!   nCellsX		Number of cells in x-direction
!
! Output: 
!   xs			Shock location
!
! Notes: 
!   1. This routine can be used for 2d cases, but then the array numbering must
!      be such that the indices follow the longest direction.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractShockLocation1D(pRegion,icgBeg,icgEnd,nCellsX,xs)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icgBeg,icgEnd,nCellsX
  REAL(RFREAL), INTENT(OUT) :: xs  
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,icg,icl,iclOffs,iclShock
  INTEGER :: dummy(1)
  REAL(RFREAL) :: idx,r,rm1,rp1,rxx,rxxp1,x,xp1
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: gradx,gradxx
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractShockLocation1D', &
                        'RFLU_ModExtractFlowDataUtils.F90')

! ******************************************************************************
! Set 
! ******************************************************************************

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv

  iclOffs = 4

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  ALLOCATE(gradx(nCellsX-2),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradx')
  END IF ! global%error

  ALLOCATE(gradxx(nCellsX-2),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradxx')
  END IF ! global%error

! ******************************************************************************
! Compute first and second density derivatives and find location of maximum 
! (absolute) value of first derivative
! ******************************************************************************

  DO icl = 1,nCellsX-2
    icg = icgBeg + icl 

    rp1 = pCv(CV_MIXT_DENS,icg+1)       
    r   = pCv(CV_MIXT_DENS,icg  )
    rm1 = pCv(CV_MIXT_DENS,icg-1)       

    gradx(icl)  = 0.5_RFREAL*(rp1-rm1)
    gradxx(icl) = rp1-2.0_RFREAL*r+rm1
  END DO ! icl    

  dummy = MAXLOC(ABS(gradx(1:nCellsX-2)))
  iclShock = dummy(1) 

! ******************************************************************************
! Starting from this location, search for zero crossing of second derivative 
! and then use linear approximation to compute location of zero crossing. NOTE 
! this search is carried out for iclOffs points on either side of location of 
! maximum (absolute) value of first derivative, so need to make sure do not step 
! out of bounds. If that is the case, simply pick shock location to be that cell
! with maximum (absolute) value of first derivative.
! ******************************************************************************

  IF ( (iclShock <= (nCellsX-iclOffs-2)) .AND. & 
       (iclShock >= (iclOffs+1)) ) THEN  
    DO icl = iclShock-iclOffs,iclShock+iclOffs-1
      icg = icgBeg + icl

      IF ( SIGN(1.0_RFREAL,gradxx(icl)) /= SIGN(1.0_RFREAL,gradxx(icl+1)) ) THEN 
        rxxp1 = gradxx(icl+1)
        rxx   = gradxx(icl)
  
        xp1 = pGrid%cofg(XCOORD,icg+1)
        x   = pGrid%cofg(XCOORD,icg)
  
        xs = (x*rxxp1-xp1*rxx)/(rxxp1-rxx)
      END IF ! SIGN  
    END DO ! icl
  ELSE 
    xs = pGrid%cofg(XCOORD,icgBeg+iclShock)
  END IF ! icl

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(gradx,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradx')
  END IF ! global%error

  DEALLOCATE(gradxx,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradxx')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractShockLocation1D





END MODULE RFLU_ModExtractFlowDataUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModExtractFlowDataUtils.F90,v $
! Revision 1.4  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/12 00:27:49  haselbac
! Changed to use zero of 2nd derivative to find shock location
!
! Revision 1.1  2007/04/05 01:38:07  haselbac
! Initial revision
!
! ******************************************************************************







