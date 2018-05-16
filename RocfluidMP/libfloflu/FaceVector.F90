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
! Purpose: calculate the face vector of a polygon.
!
! Description: file contains the following subroutines:
!
!  - faceVectorTria = triangle
!  - FaceVectorQuad = quadrilateral
!
! Input: xyzNodes = coordinates (1st index) of the nodes (2nd index) of the
!                   polygon (clockwise ordered).
!
! Output: fVecX, fVecY, fVecZ = x-,y-,z-component of the face vector.
!
! Notes: Gauss formula is used to calculate the face vector. The coordinates
!        should be ordered clockwise (when looking from inside the volume) to
!        obtain an outward facing vector.
!
!******************************************************************************
!
! $Id: FaceVector.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE FaceVectorTria( xyzNodes,fVecX,fVecY,fVecZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,3)
  REAL(RFREAL) :: fVecX, fVecY, fVecZ

! ... local variables
  REAL(RFREAL) :: dxa,dxb,dya,dyb,dza,dzb

!******************************************************************************

  dxa = xyzNodes(1,1) - xyzNodes(1,2)
  dya = xyzNodes(2,1) - xyzNodes(2,2) 
  dza = xyzNodes(3,1) - xyzNodes(3,2)  
  
  dxb = xyzNodes(1,1) - xyzNodes(1,3)
  dyb = xyzNodes(2,1) - xyzNodes(2,3) 
  dzb = xyzNodes(3,1) - xyzNodes(3,3)
    
  fVecX =  0.5_RFREAL*(dya*dzb - dyb*dza)
  fVecY = -0.5_RFREAL*(dxa*dzb - dxb*dza)
  fVecZ =  0.5_RFREAL*(dxa*dyb - dxb*dya)         

END SUBROUTINE FaceVectorTria

! #############################################################################
! #############################################################################

SUBROUTINE FaceVectorQuad( xyzNodes,fVecX,fVecY,fVecZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,4)
  REAL(RFREAL) :: fVecX, fVecY, fVecZ

! ... local variables
  REAL(RFREAL) :: dxa, dya, dza, dxb, dyb, dzb

!******************************************************************************

  dxa = xyzNodes(1,3) - xyzNodes(1,1)
  dya = xyzNodes(2,3) - xyzNodes(2,1)
  dza = xyzNodes(3,3) - xyzNodes(3,1)

  dxb = xyzNodes(1,2) - xyzNodes(1,4)
  dyb = xyzNodes(2,2) - xyzNodes(2,4)
  dzb = xyzNodes(3,2) - xyzNodes(3,4)

  fVecX = 0.5_RFREAL*(dza*dyb-dya*dzb)
  fVecY = 0.5_RFREAL*(dxa*dzb-dza*dxb)
  fVecZ = 0.5_RFREAL*(dya*dxb-dxa*dyb)

END SUBROUTINE FaceVectorQuad

!******************************************************************************
!
! RCS Revision history:
!
! $Log: FaceVector.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:34  haselbac
! Initial revision after changing case
!
! Revision 1.2  2003/03/27 14:28:47  haselbac
! Substantially simplified faceVectorTria
!
! Revision 1.1  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
!******************************************************************************






