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
! Purpose: calculate cell volumes.
!
! Description: cell volumes are calculated based on the divergence theorem.
!              See the corresponding subroutine for the reference.
!
! Input: region%levels%grid = dimensions, coordinates, face vectors
!                             (current region)
!
! Output: region%levels%grid%vol = cell volumes
!
! Notes: volumes are first calculated for all cells. In a second step,
!        volumes are copied from interior to dummy cells at all but
!        inter-region or periodic boundaries. Face vectors must be already
!        defined before this routine is called.
!
!******************************************************************************
!
! $Id: RFLO_CalcControlVolumes.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcControlVolumes( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, VolumeHexa, &
        RFLO_copyVectorPatches, RFLO_copyVectorEdges, RFLO_copyVectorCorners
  USE ModError
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, i, j, k

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkCell, corner(8)

  REAL(RFREAL)          :: xyzHexa(3,8), faceVecs(3,6)
  REAL(RFREAL), POINTER :: xyz(:,:), si(:,:), sj(:,:), sk(:,:), vol(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcControlVolumes',&
  'RFLO_CalcControlVolumes.F90' )

! volumes of all cells

  DO iLev=1,region%nGridLevels               ! all grid levels

    CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    xyz => region%levels(iLev)%grid%xyz
    si  => region%levels(iLev)%grid%si
    sj  => region%levels(iLev)%grid%sj
    sk  => region%levels(iLev)%grid%sk
    vol => region%levels(iLev)%grid%vol

    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkCell   = IndIJK(i,j,k,iCOff,ijCOff)
          corner(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          corner(2) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
          corner(3) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
          corner(4) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          corner(5) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          corner(6) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
          corner(7) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
          corner(8) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)

          xyzHexa(1:3,1) = xyz(1:3,corner(1))
          xyzHexa(1:3,2) = xyz(1:3,corner(2))
          xyzHexa(1:3,3) = xyz(1:3,corner(3))
          xyzHexa(1:3,4) = xyz(1:3,corner(4))
          xyzHexa(1:3,5) = xyz(1:3,corner(5))
          xyzHexa(1:3,6) = xyz(1:3,corner(6))
          xyzHexa(1:3,7) = xyz(1:3,corner(7))
          xyzHexa(1:3,8) = xyz(1:3,corner(8))

          faceVecs(1:3,1) =  si(1:3,corner(1))
          faceVecs(1:3,2) = -si(1:3,corner(5))
          faceVecs(1:3,3) =  sj(1:3,corner(1))
          faceVecs(1:3,4) = -sj(1:3,corner(4))
          faceVecs(1:3,5) =  sk(1:3,corner(1))
          faceVecs(1:3,6) = -sk(1:3,corner(2))

          CALL VolumeHexa( xyzHexa,faceVecs,vol(ijkCell) )
        ENDDO
      ENDDO
    ENDDO

! - copy to dummy cells at patches

    CALL RFLO_copyVectorPatches( iLev,region,vol )

! - copy to dummy cells at edges

    CALL RFLO_copyVectorEdges( iLev,region,vol )

! - copy to dummy cells at corners

    CALL RFLO_copyVectorCorners( iLev,region,vol )

  ENDDO   ! iLev

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcControlVolumes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcControlVolumes.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
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







