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
! Purpose: Calculate the face centroid of a triangular or quadrilateral face.
!
! Description: File contains the following subroutines:
!  - FaceCentroidTria = Centroid of a triangle
!  - FaceCentroidQuad = Centroid of a quadrilateral
!
! Input: xyzNodes = coordinates (1st index) of the nodes (2nd index) of the
!                   polygon (clockwise ordered).
!
! Output: fCenX, fCenY, fCenZ = x-,y-,z-component of the face centroid.
!
! Notes:
!   1. Quadrilateral face centroid is only approximate for distorted faces.
!
!******************************************************************************
!
! $Id: FaceCentroid.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE FaceCentroidTria( xyz,fCenX,fCenY,fCenZ )

  USE ModDataTypes

  IMPLICIT NONE

! ... parameters

  REAL(RFREAL), INTENT(IN) :: xyz(3,3)
  REAL(RFREAL), INTENT(OUT) :: fCenX, fCenY, fCenZ

! ... local variables

  REAL(RFREAL), PARAMETER :: thrd = 1.0_RFREAL/3.0_RFREAL 
  
!******************************************************************************

  fCenX = thrd*(xyz(1,1) + xyz(1,2) + xyz(1,3))
  fCenY = thrd*(xyz(2,1) + xyz(2,2) + xyz(2,3))
  fCenZ = thrd*(xyz(3,1) + xyz(3,2) + xyz(3,3))

END SUBROUTINE FaceCentroidTria

! #############################################################################
! #############################################################################

SUBROUTINE FaceCentroidQuad( xyz,fCenX,fCenY,fCenZ )

  USE ModDataTypes
  
  IMPLICIT NONE

! ... parameters

  REAL(RFREAL), INTENT(IN) :: xyz(3,4)
  REAL(RFREAL), INTENT(OUT) :: fCenX, fCenY, fCenZ

! ... local variables

!******************************************************************************

  fCenX = 0.25_RFREAL*(xyz(1,1) + xyz(1,2) + xyz(1,3) + xyz(1,4))
  fCenY = 0.25_RFREAL*(xyz(2,1) + xyz(2,2) + xyz(2,3) + xyz(2,4))
  fCenZ = 0.25_RFREAL*(xyz(3,1) + xyz(3,2) + xyz(3,3) + xyz(3,4))

END SUBROUTINE FaceCentroidQuad

!******************************************************************************
!
! RCS Revision history:
!
! $Log: FaceCentroid.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:32  haselbac
! Initial revision after changing case
!
! Revision 1.1  2002/03/14 19:01:10  haselbac
! Initial revision
!
!******************************************************************************






