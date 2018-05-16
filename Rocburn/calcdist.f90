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
MODULE M_CALCDIST


  USE  M_ROCBURN_INTERFACE_DATA, ONLY : list_block, DBL
  IMPLICIT NONE

CONTAINS

  SUBROUTINE CALCDIST_2D( G_b, xyz_2d, dist_2d) !  For computing distance.
    TYPE( list_block), POINTER  :: G_b
    REAL(DBL), INTENT (IN)      :: xyz_2d(:,:)
    REAL(DBL), INTENT (OUT)     :: dist_2d(:)

    dist_2d(:) = xyz_2d(1,:)


  END SUBROUTINE CALCDIST_2D

END MODULE M_CALCDIST






