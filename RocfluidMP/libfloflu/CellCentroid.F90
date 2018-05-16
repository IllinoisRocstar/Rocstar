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
! Purpose: calculate the centroid of a cell.
!
! Description: file contains the following subroutines:
!
!  - CentroidTetra   = tetrahedron
!  - CentroidPyramid = pyramid
!  - CentroidPrism   = prism
!  - CentroidHexa    = hexahedra
!
! Input: xyzNodes = coordinates (1st index) of the nodes (2nd index) of the
!                   control volume (must be ordered!).
!
! Output: cofgX,cofgY,cofgZ = x-,y-,z-coordinate of the centroid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: CellCentroid.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE CentroidTetra( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,4)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

! ... local variables


!******************************************************************************

END SUBROUTINE CentroidTetra

! #############################################################################
! #############################################################################

SUBROUTINE CentroidPyramid( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,5)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

! ... local variables


!******************************************************************************

END SUBROUTINE CentroidPyramid

! #############################################################################
! #############################################################################

SUBROUTINE CentroidPrism( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,6)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

! ... local variables


!******************************************************************************

END SUBROUTINE CentroidPrism

! #############################################################################
! #############################################################################

SUBROUTINE CentroidHexa( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,8)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

!******************************************************************************

  cofgX = 0.125_RFREAL*(xyzNodes(1,1)+xyzNodes(1,2)+xyzNodes(1,3)+ &
                        xyzNodes(1,4)+xyzNodes(1,5)+xyzNodes(1,6)+ &
                        xyzNodes(1,7)+xyzNodes(1,8))
  cofgY = 0.125_RFREAL*(xyzNodes(2,1)+xyzNodes(2,2)+xyzNodes(2,3)+ &
                        xyzNodes(2,4)+xyzNodes(2,5)+xyzNodes(2,6)+ &
                        xyzNodes(2,7)+xyzNodes(2,8))
  cofgZ = 0.125_RFREAL*(xyzNodes(3,1)+xyzNodes(3,2)+xyzNodes(3,3)+ &
                        xyzNodes(3,4)+xyzNodes(3,5)+xyzNodes(3,6)+ &
                        xyzNodes(3,7)+xyzNodes(3,8))

END SUBROUTINE CentroidHexa

!******************************************************************************
!
! RCS Revision history:
!
! $Log: CellCentroid.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:18  haselbac
! Initial revision after changing case
!
! Revision 1.1  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
!******************************************************************************






