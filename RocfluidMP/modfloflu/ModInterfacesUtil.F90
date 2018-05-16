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
! Purpose: set explicit interfaces to utility subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesUtil.F90,v 1.3 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesUtil

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvar )
    USE ModDataTypes
    INTEGER               :: ijkD, ijkN1, ijkN2, iFBeg, iFEnd
    REAL(RFREAL), POINTER :: fvar(:,:)
  END SUBROUTINE AverageVecVar

  SUBROUTINE CentroidHexa( xyzNodes,cofgX,cofgY,cofgZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,8)
    REAL(RFREAL) :: cofgX, cofgY, cofgZ
  END SUBROUTINE CentroidHexa

  SUBROUTINE DescaleGridSpeeds(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region    
  END SUBROUTINE DescaleGridSpeeds

  SUBROUTINE FaceCentroidQuad( xyzNodes,fCenX,fCenY,fCenZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,4)
    REAL(RFREAL) :: fCenX, fCenY, fCenZ
  END SUBROUTINE FaceCentroidQuad

  SUBROUTINE FaceCentroidTria( xyzNodes,fCenX,fCenY,fCenZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,3)
    REAL(RFREAL) :: fCenX, fCenY, fCenZ
  END SUBROUTINE FaceCentroidTria

  SUBROUTINE FaceVectorQuad( xyzNodes,fVecX,fVecY,fVecZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,4)
    REAL(RFREAL) :: fVecX, fVecY, fVecZ
  END SUBROUTINE FaceVectorQuad

  SUBROUTINE FaceVectorTria( xyzNodes,fVecX,fVecY,fVecZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,3)
    REAL(RFREAL) :: fVecX, fVecY, fVecZ
  END SUBROUTINE FaceVectorTria

  SUBROUTINE ScaleGridSpeeds(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region    
  END SUBROUTINE ScaleGridSpeeds

  SUBROUTINE ScaleRotateVector( global,vect )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    REAL(RFREAL), POINTER   :: vect(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ScaleRotateVector

  SUBROUTINE SplitQuadFace( global,xyz1,xyz2,xyz3,xyz4,splitFlag )
    USE ModDataTypes
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
    INTEGER, INTENT(OUT)    :: splitFlag
    REAL(RFREAL), DIMENSION(3), INTENT(IN) :: xyz1,xyz2,xyz3,xyz4
  END SUBROUTINE SplitQuadFace

  SUBROUTINE VolumeHexa( xyzNodes,faceVecs,volume )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,8), faceVecs(3,6)
    REAL(RFREAL) :: volume
  END SUBROUTINE VolumeHexa

  END INTERFACE

END MODULE ModInterfacesUtil

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesUtil.F90,v $
! Revision 1.3  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
!******************************************************************************






