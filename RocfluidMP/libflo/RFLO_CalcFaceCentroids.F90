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
! Purpose: calculate face centroids (optionally).
!
! Description: none.
!
! Input: region%levels%grid = dimensions, coordinates (current region)
!
! Output: region%levels%grid%cfcI,J,K = i,j,k-face centroids
!
!
!******************************************************************************
!
! $Id: RFLO_CalcFaceCentroids.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcFaceCentroids( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, i, j, k, l

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: iNOff, ijNOff
  INTEGER, PARAMETER :: IDIR=1, JDIR=2, KDIR=3

  REAL(RFREAL), POINTER :: xyz(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcFaceCentroids',&
       'RFLO_CalcFaceCentroids.F90' )

! loop over all grid levels

  DO iLev=1,region%nGridLevels

    IF (region%global%calcFaceCtr) THEN
      CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                     jdnbeg,jdnend,kdnbeg,kdnend )
      CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

      xyz  => region%levels(iLev)%grid%xyz

! --- compute i, j, k-face centroids

      CALL ComputeFaceCtr( IDIR )
      CALL ComputeFaceCtr( JDIR )
      CALL ComputeFaceCtr( KDIR )

    ENDIF  ! calcFaceCtr
  ENDDO    ! iLev

! finalize

  CALL DeregisterFunction( region%global )

! =============================================================================
!   Face centroids computation subroutine
! =============================================================================

CONTAINS

  SUBROUTINE ComputeFaceCtr( ijk )

! ... parameters
    INTEGER   :: ijk 

! ... local variables
    INTEGER   :: ibeg,iend,jbeg,jend,kbeg,kend
    INTEGER   :: i2,i3,i4,j2,j3,j4,k2,k3,k4
    INTEGER   :: corner(4)
    REAL(RFREAL), POINTER :: faceCtr(:,:)

! - Set limits and pointers ---------------------------------------------------
    IF (ijk==IDIR) THEN
      ibeg = idnbeg
      iend = idnend
      jbeg = jdnbeg
      jend = jdnend-1
      kbeg = kdnbeg
      kend = kdnend-1
      i2   = 0
      i3   = 0
      i4   = 0
      j2   = 0
      j3   = 1
      j4   = 1
      k2   = 1
      k3   = 1
      k4   = 0
      faceCtr => region%levels(iLev)%grid%cfcI
    ELSEIF (ijk==JDIR) THEN
      ibeg = idnbeg
      iend = idnend-1
      jbeg = jdnbeg
      jend = jdnend
      kbeg = kdnbeg
      kend = kdnend-1
      i2   = 0
      i3   = 1
      i4   = 1
      j2   = 0
      j3   = 0
      j4   = 0
      k2   = 1
      k3   = 1
      k4   = 0
      faceCtr => region%levels(iLev)%grid%cfcJ
    ELSEIF (ijk==KDIR) THEN
      ibeg = idnbeg
      iend = idnend-1
      jbeg = jdnbeg
      jend = jdnend-1
      kbeg = kdnbeg
      kend = kdnend
      i2   = 0
      i3   = 1
      i4   = 1
      j2   = 1
      j3   = 1
      j4   = 0
      k2   = 0
      k3   = 0
      k4   = 0
      faceCtr => region%levels(iLev)%grid%cfcK
    ENDIF

! - define face centroids

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          corner(1) = IndIJK(i    ,j    ,k    ,iNOff,ijNOff)
          corner(2) = IndIJK(i+i2 ,j+j2 ,k+k2 ,iNOff,ijNOff)
          corner(3) = IndIJK(i+i3 ,j+j3 ,k+k3 ,iNOff,ijNOff)
          corner(4) = IndIJK(i+i4 ,j+j4 ,k+k4 ,iNOff,ijNOff)
         
          DO l=1,3
            faceCtr(l,corner(1)) = 0.25_RFREAL*(xyz(l,corner(1))+ &
                                                xyz(l,corner(2))+ &
                                                xyz(l,corner(3))+ &
                                                xyz(l,corner(4)))
          ENDDO         
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

    IF (ijk==IDIR) THEN

      DO j=jdnbeg,jdnend-1
        DO i=idnbeg,idnend
          corner(1) = IndIJK(i    ,j    ,kdnend  ,iNOff,ijNOff)
          corner(2) = IndIJK(i    ,j    ,kdnend-1,iNOff,ijNOff)
          DO l=1,3
            faceCtr(l,corner(1)) = faceCtr(l,corner(2))
          ENDDO
        ENDDO
      ENDDO
      DO k=kdnbeg,kdnend
        DO i=idnbeg,idnend
          corner(1) = IndIJK(i    ,jdnend  ,k    ,iNOff,ijNOff)
          corner(2) = IndIJK(i    ,jdnend-1,k    ,iNOff,ijNOff)
          DO l=1,3
            faceCtr(l,corner(1)) = faceCtr(l,corner(2))
          ENDDO
        ENDDO
      ENDDO

    ELSEIF (ijk==JDIR) THEN

      DO k=kdnbeg,kdnend-1
        DO j=jdnbeg,jdnend
          corner(1) = IndIJK(idnend  ,j    ,k    ,iNOff,ijNOff)
          corner(2) = IndIJK(idnend-1,j    ,k    ,iNOff,ijNOff)
          DO l=1,3
            faceCtr(l,corner(1)) = faceCtr(l,corner(2))
          ENDDO
        ENDDO
      ENDDO
      DO j=jdnbeg,jdnend
        DO i=idnbeg,idnend
          corner(1) = IndIJK(i    ,j     ,kdnend  ,iNOff,ijNOff)
          corner(2) = IndIJK(i    ,j     ,kdnend-1,iNOff,ijNOff)
          DO l=1,3
            faceCtr(l,corner(1)) = faceCtr(l,corner(2))
          ENDDO
        ENDDO
      ENDDO

    ELSEIF (ijk==KDIR) THEN

      DO k=kdnbeg,kdnend
        DO j=jdnbeg,jdnend-1
          corner(1) = IndIJK(idnend  ,j    ,k    ,iNOff,ijNOff)
          corner(2) = IndIJK(idnend-1,j    ,k    ,iNOff,ijNOff)
          DO l=1,3
            faceCtr(l,corner(1)) = faceCtr(l,corner(2))
          ENDDO
        ENDDO
      ENDDO
      DO k=kdnbeg,kdnend
        DO i=idnbeg,idnend
          corner(1) = IndIJK(i    ,jdnend  ,k    ,iNOff,ijNOff)
          corner(2) = IndIJK(i    ,jdnend-1,k    ,iNOff,ijNOff)
          DO l=1,3
            faceCtr(l,corner(1)) = faceCtr(l,corner(2))
          ENDDO
        ENDDO
      ENDDO

    ENDIF   ! ijk

  END SUBROUTINE ComputeFaceCtr

END SUBROUTINE RFLO_CalcFaceCentroids

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcFaceCentroids.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/01/12 09:47:02  wasistho
! extrapolate to idnend, jdnend, kdnend
!
! Revision 1.1  2005/10/20 06:57:44  wasistho
! initial import RFLO_CalcFaceCentroids
!
!
!******************************************************************************







