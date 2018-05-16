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
! Purpose: calculate face vectors.
!
! Description: face vectors are computed using Gauss` formula. This results
!              in an average vector in the case of non-planar faces.
!
! Input: region%levels%grid = dimensions, coordinates (current region)
!
! Output: region%levels%grid%si/j/k = face vectors
!
! Notes: face vectors are calculated for all dummy cells, including
!        corner and edge cells.
!
!******************************************************************************
!
! $Id: RFLO_CalcFaceVectors.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcFaceVectors( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, &
                            FaceVectorQuad
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, i, j, k

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: iNOff, ijNOff, corner(4)

  REAL(RFREAL)          :: xyzQuad(3,4)
  REAL(RFREAL), POINTER :: xyz(:,:), si(:,:), sj(:,:), sk(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcFaceVectors',&
  'RFLO_CalcFaceVectors.F90' )

! loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    xyz => region%levels(iLev)%grid%xyz
    si  => region%levels(iLev)%grid%si
    sj  => region%levels(iLev)%grid%sj
    sk  => region%levels(iLev)%grid%sk

! - i-direction

    DO k=kdnbeg,kdnend-1
      DO j=jdnbeg,jdnend-1
        DO i=idnbeg,idnend
          corner(1) = IndIJK(i,j  ,k  ,iNOff,ijNOff)
          corner(2) = IndIJK(i,j  ,k+1,iNOff,ijNOff)
          corner(3) = IndIJK(i,j+1,k+1,iNOff,ijNOff)
          corner(4) = IndIJK(i,j+1,k  ,iNOff,ijNOff)
          xyzQuad(1:3,1) = xyz(1:3,corner(1))
          xyzQuad(1:3,2) = xyz(1:3,corner(2))
          xyzQuad(1:3,3) = xyz(1:3,corner(3))
          xyzQuad(1:3,4) = xyz(1:3,corner(4))
          CALL FaceVectorQuad( xyzQuad, &
                               si(XCOORD,corner(1)), &
                               si(YCOORD,corner(1)), &
                               si(ZCOORD,corner(1)) )
          si(XYZMAG,corner(1)) = 1._RFREAL/ &
                         SQRT( si(XCOORD,corner(1))*si(XCOORD,corner(1)) + &
                               si(YCOORD,corner(1))*si(YCOORD,corner(1)) + &
                               si(ZCOORD,corner(1))*si(ZCOORD,corner(1)) )
        ENDDO
      ENDDO
    ENDDO

! - j-direction

    DO k=kdnbeg,kdnend-1
      DO j=jdnbeg,jdnend
        DO i=idnbeg,idnend-1
          corner(1) = IndIJK(i  ,j,k  ,iNOff,ijNOff)
          corner(2) = IndIJK(i+1,j,k  ,iNOff,ijNOff)
          corner(3) = IndIJK(i+1,j,k+1,iNOff,ijNOff)
          corner(4) = IndIJK(i  ,j,k+1,iNOff,ijNOff)
          xyzQuad(1:3,1) = xyz(1:3,corner(1))
          xyzQuad(1:3,2) = xyz(1:3,corner(2))
          xyzQuad(1:3,3) = xyz(1:3,corner(3))
          xyzQuad(1:3,4) = xyz(1:3,corner(4))
          CALL FaceVectorQuad( xyzQuad, &
                               sj(XCOORD,corner(1)), &
                               sj(YCOORD,corner(1)), &
                               sj(ZCOORD,corner(1)) )
          sj(XYZMAG,corner(1)) = 1._RFREAL/ &
                         SQRT( sj(XCOORD,corner(1))*sj(XCOORD,corner(1)) + &
                               sj(YCOORD,corner(1))*sj(YCOORD,corner(1)) + &
                               sj(ZCOORD,corner(1))*sj(ZCOORD,corner(1)) )
        ENDDO
      ENDDO
    ENDDO

! - k-direction

    DO k=kdnbeg,kdnend
      DO j=jdnbeg,jdnend-1
        DO i=idnbeg,idnend-1
          corner(1) = IndIJK(i  ,j  ,k,iNOff,ijNOff)
          corner(2) = IndIJK(i  ,j+1,k,iNOff,ijNOff)
          corner(3) = IndIJK(i+1,j+1,k,iNOff,ijNOff)
          corner(4) = IndIJK(i+1,j  ,k,iNOff,ijNOff)
          xyzQuad(1:3,1) = xyz(1:3,corner(1))
          xyzQuad(1:3,2) = xyz(1:3,corner(2))
          xyzQuad(1:3,3) = xyz(1:3,corner(3))
          xyzQuad(1:3,4) = xyz(1:3,corner(4))
          CALL FaceVectorQuad( xyzQuad, &
                               sk(XCOORD,corner(1)), &
                               sk(YCOORD,corner(1)), &
                               sk(ZCOORD,corner(1)) )
          sk(XYZMAG,corner(1)) = 1._RFREAL/ &
                         SQRT( sk(XCOORD,corner(1))*sk(XCOORD,corner(1)) + &
                               sk(YCOORD,corner(1))*sk(YCOORD,corner(1)) + &
                               sk(ZCOORD,corner(1))*sk(ZCOORD,corner(1)) )
        ENDDO
      ENDDO
    ENDDO

  ENDDO   ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcFaceVectors

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcFaceVectors.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/24 05:04:50  wasistho
! computed inverse magnitude of face vectors
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.9  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
!******************************************************************************







