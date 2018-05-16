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
! Purpose: Define the derived data types related to grids in ROCFLU.
!
! Description: None
!
! Notes: None
!
! ******************************************************************************
!
! $Id: RFLU_ModGrid.F90,v 1.15 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2001-2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGrid

  USE ModDataTypes
  USE ModParameters
  
  IMPLICIT NONE

! ******************************************************************************
! Boundary type
! ******************************************************************************  
  
  TYPE bound_t
    INTEGER :: bType,nBTris,nBQuads
    INTEGER, DIMENSION(:,:), POINTER :: bTri2v,bQuad2v
    CHARACTER(CHRLEN) :: bName
  END TYPE bound_t

  TYPE(bound_t), DIMENSION(:), ALLOCATABLE :: bound

! ******************************************************************************
! Mapping from faces to vertices for ROCFLU grid definition
! ******************************************************************************

  INTEGER, DIMENSION(4,4), PARAMETER :: f2vTet = &
    RESHAPE((/1,2,3,VERT_NONE,2,4,3,VERT_NONE,1,3,4,VERT_NONE,1,4,2, & 
              VERT_NONE/), (/4,4/))
  INTEGER, DIMENSION(4,6), PARAMETER :: f2vHex = &
    RESHAPE((/1,4,3,2,1,2,6,5,2,3,7,6,3,4,8,7,1,5,8,4,5,6,7,8/), (/4,6/))
  INTEGER, DIMENSION(4,5), PARAMETER :: f2vPri = &
    RESHAPE((/1,3,2,VERT_NONE,1,2,5,4,2,3,6,5,1,4,6,3,4,5,6, & 
              VERT_NONE/), (/4,5/))
  INTEGER, DIMENSION(4,5), PARAMETER :: f2vPyr = &
    RESHAPE((/1,4,3,2,1,2,5,VERT_NONE,2,3,5,VERT_NONE,3,4,5,VERT_NONE,1,5,4, & 
              VERT_NONE/), (/4,5/))

! ******************************************************************************
! Mapping from faces to opposite faces for ROCFLU grid definition
! ******************************************************************************

  INTEGER, DIMENSION(6), PARAMETER :: f2fOppHex = (/6,4,5,2,3,1/) 
  INTEGER, DIMENSION(5), PARAMETER :: f2fOppPri = & 
    (/5,FACE_NONE,FACE_NONE,FACE_NONE,1/)
              
! ******************************************************************************
! Mapping from cell edges to vertices for ROCFLU grid definition
! ******************************************************************************

  INTEGER, DIMENSION(2,6), PARAMETER  :: ce2vTet = & 
    RESHAPE((/1,2,2,3,1,3,1,4,2,4,3,4/), (/2,6/))
  INTEGER, DIMENSION(2,12), PARAMETER :: ce2vHex = & 
    RESHAPE((/1,2,2,3,3,4,1,4,5,6,6,7,7,8,5,8,1,5,2,6,3,7,4,8/), (/2,12/))
  INTEGER, DIMENSION(2,9), PARAMETER  :: ce2vPri = & 
    RESHAPE((/1,2,2,3,1,3,1,4,2,5,3,6,4,5,5,6,4,6/), (/2,9/))
  INTEGER, DIMENSION(2,8), PARAMETER  :: ce2vPyr = & 
    RESHAPE((/1,2,2,3,3,4,1,4,1,5,2,5,3,5,4,5/), (/2,8/))

END MODULE RFLU_ModGrid

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGrid.F90,v $
! Revision 1.15  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2004/09/27 01:39:26  haselbac
! Added arrays for opposing faces
!
! Revision 1.12  2004/07/06 15:14:41  haselbac
! Removed data types for grid generators and patchDimens
!
! Revision 1.11  2004/02/16 01:10:27  haselbac
! Added PARAMETER qualifyer to avoid problems on Apple OS X
!
! Revision 1.10  2003/08/19 22:49:07  haselbac
! Added COBALT grid type
!
! Revision 1.9  2003/03/19 16:47:09  haselbac
! Added gridTETMESH type
!
! Revision 1.8  2003/03/15 18:10:40  haselbac
! Added VERT_NONE, *2vTEC, changed patchDimens
!
! Revision 1.7  2003/01/28 16:32:05  haselbac
! Added patchDimens and ModParameters
!
! Revision 1.6  2002/10/27 19:07:39  haselbac
! Added c2v arrays for edge list construction
!
! Revision 1.5  2002/03/01 16:22:32  haselbac
! Added face to vertex mapping arrays
!
! Revision 1.4  2002/01/15 15:52:07  haselbac
! Added grid types for grid generators
!
! ******************************************************************************






