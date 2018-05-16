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
! Purpose: redistribute the interior grid based on new discretization
!          of the boundaries.
!
! Description: the method used is linear transfinite interpolation (TFI).
!
! Input: region     = grid dimensions
!        boundMoved = flag for boundaries of a region which have moved
!        edgeMoved  = flag for edges whose nodes have moved
!        arcLen12   = arclength between i=const. boundaries for each j, k
!        arcLen34   = arclength between j=const. boundaries for each k, i
!        arcLen56   = arclength between k=const. boundaries for each i, j
!        xyzOld     = grid from previous time step
!        xyz        = deformations at the boundaries of region.
!
! Output: xyz = new grid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ChangeInteriorGrid.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ChangeInteriorGrid( region,boundMoved,edgeMoved, &
                                    arcLen12,arcLen34,arcLen56,xyzOld,xyz )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6), edgeMoved(12)

  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijkN, imjkN, ijmkN, ijkmN, iNOff, ijNOff, errorFlag

  REAL(RFREAL) :: phii, phii1, phij, phij1, phik, phik1, dsi
  REAL(RFREAL) :: v1(3), v2(3), v3(3), v12(3), v13(3), v23(3), v123(3)
  REAL(RFREAL), ALLOCATABLE :: dsj(:), dsk(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ChangeInteriorGrid',&
  'RFLO_ChangeInteriorGrid.F90' )

! get dimensions, allocate temporary storage ----------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  ALLOCATE( dsj(ipnbeg:ipnend)              ,stat=errorFlag )
  ALLOCATE( dsk(ipnbeg:ipnend,jpnbeg:jpnend),stat=errorFlag )
  region%global%error = errorFlag
  IF (region%global%error /= 0) &
    CALL ErrorStop( region%global,ERR_ALLOCATE,&
    __LINE__ )

! interpolate displacements inside region -------------------------------------

  dsk(:,:) = 0._RFREAL
  DO k=kpnbeg+1,kpnend-1
    dsj(:) = 0._RFREAL
    DO j=jpnbeg+1,jpnend-1
      dsi  = 0._RFREAL
      DO i=ipnbeg+1,ipnend-1
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        imjkN   = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ijmkN   = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijkmN   = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

        dsi      = dsi      + &
                     SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,imjkN))**2 + &
                          (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,imjkN))**2 + &
                          (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,imjkN))**2)
        dsj(i)   = dsj(i)   + &
                     SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijmkN))**2 + &
                          (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijmkN))**2 + &
                          (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijmkN))**2)
        dsk(i,j) = dsk(i,j) + &
                     SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkmN))**2 + &
                          (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkmN))**2 + &
                          (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkmN))**2)

        phii  = dsi/arcLen12(j,k)
        phii1 = 1._RFREAL - phii
        phij  = dsj(i)/arcLen34(k,i)
        phij1 = 1._RFREAL - phij
        phik  = dsk(i,j)/arcLen56(i,j)
        phik1 = 1._RFREAL - phik

        v1(:)   = phii1*xyz(:,IndIJK(ipnbeg,j,k,iNOff,ijNOff)) + &
                  phii *xyz(:,IndIJK(ipnend,j,k,iNOff,ijNOff))
        v2(:)   = phij1*xyz(:,IndIJK(i,jpnbeg,k,iNOff,ijNOff)) + &
                  phij *xyz(:,IndIJK(i,jpnend,k,iNOff,ijNOff))
        v3(:)   = phik1*xyz(:,IndIJK(i,j,kpnbeg,iNOff,ijNOff)) + &
                  phik *xyz(:,IndIJK(i,j,kpnend,iNOff,ijNOff))

        v12(:)  = phii1*phij1*xyz(:,IndIJK(ipnbeg,jpnbeg,k,iNOff,ijNOff)) + &
                  phii1*phij *xyz(:,IndIJK(ipnbeg,jpnend,k,iNOff,ijNOff)) + &
                  phii *phij1*xyz(:,IndIJK(ipnend,jpnbeg,k,iNOff,ijNOff)) + &
                  phii *phij *xyz(:,IndIJK(ipnend,jpnend,k,iNOff,ijNOff))
        v13(:)  = phii1*phik1*xyz(:,IndIJK(ipnbeg,j,kpnbeg,iNOff,ijNOff)) + &
                  phii1*phik *xyz(:,IndIJK(ipnbeg,j,kpnend,iNOff,ijNOff)) + &
                  phii *phik1*xyz(:,IndIJK(ipnend,j,kpnbeg,iNOff,ijNOff)) + &
                  phii *phik *xyz(:,IndIJK(ipnend,j,kpnend,iNOff,ijNOff))
        v23(:)  = phij1*phik1*xyz(:,IndIJK(i,jpnbeg,kpnbeg,iNOff,ijNOff)) + &
                  phij1*phik *xyz(:,IndIJK(i,jpnbeg,kpnend,iNOff,ijNOff)) + &
                  phij *phik1*xyz(:,IndIJK(i,jpnend,kpnbeg,iNOff,ijNOff)) + &
                  phij *phik *xyz(:,IndIJK(i,jpnend,kpnend,iNOff,ijNOff))

        v123(:) = phii1*phij1*phik1* &
                    xyz(:,IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff)) + &
                  phii1*phij1*phik * &
                    xyz(:,IndIJK(ipnbeg,jpnbeg,kpnend,iNOff,ijNOff)) + &
                  phii1*phij *phik1* &
                    xyz(:,IndIJK(ipnbeg,jpnend,kpnbeg,iNOff,ijNOff)) + &
                  phii *phij1*phik1* &
                    xyz(:,IndIJK(ipnend,jpnbeg,kpnbeg,iNOff,ijNOff)) + &
                  phii1*phij *phik * &
                    xyz(:,IndIJK(ipnbeg,jpnend,kpnend,iNOff,ijNOff)) + &
                  phii *phij1*phik * &
                    xyz(:,IndIJK(ipnend,jpnbeg,kpnend,iNOff,ijNOff)) + &
                  phii *phij *phik1* &
                    xyz(:,IndIJK(ipnend,jpnend,kpnbeg,iNOff,ijNOff)) + &
                  phii *phij *phik * &
                    xyz(:,IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff))

        xyz(:,ijkN) = v1(:) + v2(:) + v3(:) - v12(:) - v13(:) - v23(:) + v123(:)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! move grid -------------------------------------------------------------------

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijkN = IndIJK(i,j,k,iNOff,ijNOff)
        xyz(XCOORD,ijkN) = region%levels(1)%gridOld%xyz(XCOORD,ijkN) + &
                           xyz(XCOORD,ijkN)
        xyz(YCOORD,ijkN) = region%levels(1)%gridOld%xyz(YCOORD,ijkN) + &
                           xyz(YCOORD,ijkN)
        xyz(ZCOORD,ijkN) = region%levels(1)%gridOld%xyz(ZCOORD,ijkN) + &
                           xyz(ZCOORD,ijkN)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

  DEALLOCATE( dsj,stat=errorFlag )
  DEALLOCATE( dsk,stat=errorFlag )
  region%global%error = errorFlag
  IF (region%global%error /= 0) &
    CALL ErrorStop( region%global,ERR_DEALLOCATE,&
    __LINE__ )

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ChangeInteriorGrid

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ChangeInteriorGrid.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/05 19:04:15  wasistho
! set computational space coordinates from initial grid
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.8  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/14 22:05:10  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.3  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
!******************************************************************************







