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
! Purpose: calculate grid speeds at the faces of the control volume.
!
! Description: none.
!
! Input: region%levels%grid(Old) = dimensions, coordinates, face vectors.
!
! Output: region%levels%grid%si/j/kVel = grid speeds.
!
! Notes: grid speeds are calculated only for the faces of the physical cells.
!
!******************************************************************************
!
! $Id: RFLO_CalcGridSpeeds.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGridSpeeds( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset, &
                            FaceVectorQuad, VolumeHexa
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, i, j, k

! ... local variables
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iNOff, ijNOff, corner(4)

  REAL(RFREAL)          :: xyzQuad(3,4), xyzHexa(3,8), faceVecs(3,6), dVol, dt
  REAL(RFREAL), POINTER :: xyzOld(:,:), xyz(:,:), siOld(:,:), sjOld(:,:), &
                           skOld(:,:), si(:,:), sj(:,:), sk(:,:), &
                           siVel(:), sjVel(:), skVel(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGridSpeeds',&
  'RFLO_CalcGridSpeeds.F90' )

! loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                  jpnbeg,jpnend,kpnbeg,kpnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    xyzOld => region%levels(iLev)%gridOld%xyz
    siOld  => region%levels(iLev)%gridOld%si
    sjOld  => region%levels(iLev)%gridOld%sj
    skOld  => region%levels(iLev)%gridOld%sk
    xyz    => region%levels(iLev)%grid%xyz
    si     => region%levels(iLev)%grid%si
    sj     => region%levels(iLev)%grid%sj
    sk     => region%levels(iLev)%grid%sk
    siVel  => region%levels(iLev)%grid%siVel
    sjVel  => region%levels(iLev)%grid%sjVel
    skVel  => region%levels(iLev)%grid%skVel

    dt = region%global%dtMin

! - i-direction

    DO k=kpnbeg,kpnend-1
      DO j=jpnbeg,jpnend-1
        DO i=ipnbeg,ipnend
          corner(1) = IndIJK(i,j  ,k  ,iNOff,ijNOff)
          corner(2) = IndIJK(i,j  ,k+1,iNOff,ijNOff)
          corner(3) = IndIJK(i,j+1,k+1,iNOff,ijNOff)
          corner(4) = IndIJK(i,j+1,k  ,iNOff,ijNOff)

          xyzHexa(1:3,1) = xyz   (1:3,corner(1))
          xyzHexa(1:3,2) = xyz   (1:3,corner(2))
          xyzHexa(1:3,3) = xyz   (1:3,corner(3))
          xyzHexa(1:3,4) = xyz   (1:3,corner(4))
          xyzHexa(1:3,5) = xyzOld(1:3,corner(1))
          xyzHexa(1:3,6) = xyzOld(1:3,corner(2))
          xyzHexa(1:3,7) = xyzOld(1:3,corner(3))
          xyzHexa(1:3,8) = xyzOld(1:3,corner(4))

          faceVecs(1:3,1) =  si   (1:3,corner(1))
          faceVecs(1:3,2) = -siOld(1:3,corner(1))

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(1))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(2))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(2))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(1))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,3), &
                               faceVecs(YCOORD,3), &
                               faceVecs(ZCOORD,3) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(3))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(4))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(4))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(3))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,4), &
                               faceVecs(YCOORD,4), &
                               faceVecs(ZCOORD,4) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(4))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(1))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(1))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(4))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,5), &
                               faceVecs(YCOORD,5), &
                               faceVecs(ZCOORD,5) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(2))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(3))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(3))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(2))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,6), &
                               faceVecs(YCOORD,6), &
                               faceVecs(ZCOORD,6) )

          CALL VolumeHexa( xyzHexa,faceVecs,dVol )

          siVel(corner(1)) = dVol/dt
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k

! - j-direction

    DO k=kpnbeg,kpnend-1
      DO j=jpnbeg,jpnend
        DO i=ipnbeg,ipnend-1
          corner(1) = IndIJK(i+1,j,k  ,iNOff,ijNOff)
          corner(2) = IndIJK(i+1,j,k+1,iNOff,ijNOff)
          corner(3) = IndIJK(i  ,j,k+1,iNOff,ijNOff)
          corner(4) = IndIJK(i  ,j,k  ,iNOff,ijNOff)

          xyzHexa(1:3,1) = xyz   (1:3,corner(1))
          xyzHexa(1:3,2) = xyz   (1:3,corner(2))
          xyzHexa(1:3,3) = xyz   (1:3,corner(3))
          xyzHexa(1:3,4) = xyz   (1:3,corner(4))
          xyzHexa(1:3,5) = xyzOld(1:3,corner(1))
          xyzHexa(1:3,6) = xyzOld(1:3,corner(2))
          xyzHexa(1:3,7) = xyzOld(1:3,corner(3))
          xyzHexa(1:3,8) = xyzOld(1:3,corner(4))

          faceVecs(1:3,1) =  sj   (1:3,corner(4))
          faceVecs(1:3,2) = -sjOld(1:3,corner(4))

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(1))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(2))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(2))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(1))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,3), &
                               faceVecs(YCOORD,3), &
                               faceVecs(ZCOORD,3) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(3))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(4))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(4))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(3))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,4), &
                               faceVecs(YCOORD,4), &
                               faceVecs(ZCOORD,4) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(4))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(1))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(1))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(4))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,5), &
                               faceVecs(YCOORD,5), &
                               faceVecs(ZCOORD,5) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(2))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(3))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(3))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(2))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,6), &
                               faceVecs(YCOORD,6), &
                               faceVecs(ZCOORD,6) )

          CALL VolumeHexa( xyzHexa,faceVecs,dVol )

          sjVel(corner(4)) = dVol/dt
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k

! - k-direction

    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend-1
        DO i=ipnbeg,ipnend-1
          corner(1) = IndIJK(i+1,j  ,k,iNOff,ijNOff)
          corner(2) = IndIJK(i  ,j  ,k,iNOff,ijNOff)
          corner(3) = IndIJK(i  ,j+1,k,iNOff,ijNOff)
          corner(4) = IndIJK(i+1,j+1,k,iNOff,ijNOff)

          xyzHexa(1:3,1) = xyz   (1:3,corner(1))
          xyzHexa(1:3,2) = xyz   (1:3,corner(2))
          xyzHexa(1:3,3) = xyz   (1:3,corner(3))
          xyzHexa(1:3,4) = xyz   (1:3,corner(4))
          xyzHexa(1:3,5) = xyzOld(1:3,corner(1))
          xyzHexa(1:3,6) = xyzOld(1:3,corner(2))
          xyzHexa(1:3,7) = xyzOld(1:3,corner(3))
          xyzHexa(1:3,8) = xyzOld(1:3,corner(4))

          faceVecs(1:3,1) =  sk   (1:3,corner(2))
          faceVecs(1:3,2) = -skOld(1:3,corner(2))

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(1))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(2))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(2))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(1))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,3), &
                               faceVecs(YCOORD,3), &
                               faceVecs(ZCOORD,3) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(3))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(4))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(4))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(3))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,4), &
                               faceVecs(YCOORD,4), &
                               faceVecs(ZCOORD,4) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(4))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(1))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(1))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(4))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,5), &
                               faceVecs(YCOORD,5), &
                               faceVecs(ZCOORD,5) )

          xyzQuad(1:3,1)  = xyzOld(1:3,corner(2))
          xyzQuad(1:3,2)  = xyzOld(1:3,corner(3))
          xyzQuad(1:3,3)  = xyz   (1:3,corner(3))
          xyzQuad(1:3,4)  = xyz   (1:3,corner(2))
          CALL FaceVectorQuad( xyzQuad, &
                               faceVecs(XCOORD,6), &
                               faceVecs(YCOORD,6), &
                               faceVecs(ZCOORD,6) )

          CALL VolumeHexa( xyzHexa,faceVecs,dVol )

          skVel(corner(2)) = dVol/dt
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k

  ENDDO   ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcGridSpeeds

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcGridSpeeds.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.9  2003/11/20 16:40:37  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.4  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.1  2002/08/27 20:44:31  jblazek
! Implemented calculation of grid speeds.
!
!******************************************************************************







