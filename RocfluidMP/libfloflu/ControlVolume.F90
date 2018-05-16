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
! Purpose: calculate the volume of a cell.
!
! Description: file contains the following subroutines:
!
!  - VolumeTetra   = tetrahedron
!  - VolumePyramid = pyramid
!  - VolumePrism   = prism
!  - VolumeHexa    = hexahedron
!
! Input: xyzNodes = coordinates (1st index) of the nodes (2nd index) of the
!                   control volume (must be ordered!);
!        faceVecs = face vectors (outward facing).
!
! Output: volume = volume of the cell.
!
! Notes: volume is calculated using the divergence theorem as described in:
!        Bruner, C.W.S.: Geometric Properties of Arbitrary Polyhedra in Terms
!        of Face Geometry. AIAA Journal, 33 (1995), p. 1350. See developer`s
!        guide for numbering of nodes.
!
!******************************************************************************
!
! $Id: ControlVolume.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE VolumeTetra( xyzNodes,faceVecs,volume )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,4), faceVecs(3,4)
  REAL(RFREAL) :: volume

! ... local variables


!******************************************************************************

END SUBROUTINE VolumeTetra

! #############################################################################
! #############################################################################

SUBROUTINE VolumePyramid( xyzNodes,faceVecs,volume )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,5), faceVecs(3,5)
  REAL(RFREAL) :: volume

! ... local variables


!******************************************************************************

END SUBROUTINE VolumePyramid

! #############################################################################
! #############################################################################

SUBROUTINE VolumePrism( xyzNodes,faceVecs,volume )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,6), faceVecs(3,5)
  REAL(RFREAL) :: volume

! ... local variables


!******************************************************************************

END SUBROUTINE VolumePrism

! #############################################################################
! #############################################################################

SUBROUTINE VolumeHexa( xyzNodes,faceVecs,volume )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,8), faceVecs(3,6)
  REAL(RFREAL) :: volume

! ... local variables
  REAL(RFREAL) :: rm1(3), rm2(3), rm3(3), rm4(3), rm5(3), rm6(3)
  REAL(RFREAL) :: rms1, rms2, rms3, rms4, rms5, rms6

!******************************************************************************

  rm1(1) = 0.25_RFREAL*(xyzNodes(1,1)+xyzNodes(1,2)+xyzNodes(1,3)+xyzNodes(1,4))
  rm1(2) = 0.25_RFREAL*(xyzNodes(2,1)+xyzNodes(2,2)+xyzNodes(2,3)+xyzNodes(2,4))
  rm1(3) = 0.25_RFREAL*(xyzNodes(3,1)+xyzNodes(3,2)+xyzNodes(3,3)+xyzNodes(3,4))

  rm2(1) = 0.25_RFREAL*(xyzNodes(1,5)+xyzNodes(1,8)+xyzNodes(1,7)+xyzNodes(1,6))
  rm2(2) = 0.25_RFREAL*(xyzNodes(2,5)+xyzNodes(2,8)+xyzNodes(2,7)+xyzNodes(2,6))
  rm2(3) = 0.25_RFREAL*(xyzNodes(3,5)+xyzNodes(3,8)+xyzNodes(3,7)+xyzNodes(3,6))

  rm3(1) = 0.25_RFREAL*(xyzNodes(1,1)+xyzNodes(1,5)+xyzNodes(1,6)+xyzNodes(1,2))
  rm3(2) = 0.25_RFREAL*(xyzNodes(2,1)+xyzNodes(2,5)+xyzNodes(2,6)+xyzNodes(2,2))
  rm3(3) = 0.25_RFREAL*(xyzNodes(3,1)+xyzNodes(3,5)+xyzNodes(3,6)+xyzNodes(3,2))

  rm4(1) = 0.25_RFREAL*(xyzNodes(1,8)+xyzNodes(1,4)+xyzNodes(1,3)+xyzNodes(1,7))
  rm4(2) = 0.25_RFREAL*(xyzNodes(2,8)+xyzNodes(2,4)+xyzNodes(2,3)+xyzNodes(2,7))
  rm4(3) = 0.25_RFREAL*(xyzNodes(3,8)+xyzNodes(3,4)+xyzNodes(3,3)+xyzNodes(3,7))

  rm5(1) = 0.25_RFREAL*(xyzNodes(1,1)+xyzNodes(1,4)+xyzNodes(1,8)+xyzNodes(1,5))
  rm5(2) = 0.25_RFREAL*(xyzNodes(2,1)+xyzNodes(2,4)+xyzNodes(2,8)+xyzNodes(2,5))
  rm5(3) = 0.25_RFREAL*(xyzNodes(3,1)+xyzNodes(3,4)+xyzNodes(3,8)+xyzNodes(3,5))

  rm6(1) = 0.25_RFREAL*(xyzNodes(1,2)+xyzNodes(1,6)+xyzNodes(1,7)+xyzNodes(1,3))
  rm6(2) = 0.25_RFREAL*(xyzNodes(2,2)+xyzNodes(2,6)+xyzNodes(2,7)+xyzNodes(2,3))
  rm6(3) = 0.25_RFREAL*(xyzNodes(3,2)+xyzNodes(3,6)+xyzNodes(3,7)+xyzNodes(3,3))

  rms1 = rm1(1)*faceVecs(1,1) + rm1(2)*faceVecs(2,1) + rm1(3)*faceVecs(3,1)
  rms2 = rm2(1)*faceVecs(1,2) + rm2(2)*faceVecs(2,2) + rm2(3)*faceVecs(3,2)
  rms3 = rm3(1)*faceVecs(1,3) + rm3(2)*faceVecs(2,3) + rm3(3)*faceVecs(3,3)
  rms4 = rm4(1)*faceVecs(1,4) + rm4(2)*faceVecs(2,4) + rm4(3)*faceVecs(3,4)
  rms5 = rm5(1)*faceVecs(1,5) + rm5(2)*faceVecs(2,5) + rm5(3)*faceVecs(3,5)
  rms6 = rm6(1)*faceVecs(1,6) + rm6(2)*faceVecs(2,6) + rm6(3)*faceVecs(3,6)

  volume = (rms1+rms2+rms3+rms4+rms5+rms6)/3._RFREAL

END SUBROUTINE VolumeHexa

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ControlVolume.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:22  haselbac
! Initial revision after changing case
!
! Revision 1.2  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.1  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
!******************************************************************************






