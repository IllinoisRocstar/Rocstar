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
! Purpose: copy geometry from a patch of the source region to dummy nodes
!          of the current region`s patch (both regions are on the same
!          processor).
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: region%levels(1)%grid%xyz = grid coordinates.
!
! Notes: operation carried out for the finest grid only. Intended for
!        conforming grid boundaries (also periodic) only. In the case
!        of rotationally periodic boundary, it is assumed that the
!        axis of rotation coincides with the x-axis.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeGeometryCopy.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeGeometryCopy( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModIndexing, ONLY   : IndIJKMap
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetPatchDirection, &
                            RFLO_GetNodeOffset, RFLO_GetPatchMapping
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: idum, i, j, k, ii, jj, kk

! ... local variables
  INTEGER :: bcType, ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, &
             iNOff, ijNOff, ijkN, ijkNB
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc, &
             idirSrc, jdirSrc, kdirSrc, iNOffSrc, ijNOffSrc, ijkNSrc, ijkNBSrc
  INTEGER :: lb, lbs, l1SrcDir, l2SrcDir, mapMat(3,4)

  LOGICAL :: align

  REAL(RFREAL)          :: dx, dy, dz, sina, cosa, v1(3), v2(3)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzSrc(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ExchangeGeometryCopy',&
  'RFLO_ExchangeGeometryCopy.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( region%global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers

  CALL RFLO_GetPatchIndicesNodes( region,patch,1,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchIndicesNodes( regionSrc,patchSrc,1,ibegSrc,iendSrc, &
                                  jbegSrc,jendSrc,kbegSrc,kendSrc )
  CALL RFLO_GetPatchDirection( patch   ,idir   ,jdir,   kdir    )
  CALL RFLO_GetPatchDirection( patchSrc,idirSrc,jdirSrc,kdirSrc )
  CALL RFLO_GetNodeOffset( region   ,1,iNOff   ,ijNOff    )
  CALL RFLO_GetNodeOffset( regionSrc,1,iNOffSrc,ijNOffSrc )

  xyz    => region%levels(1)%grid%xyz
  xyzSrc => regionSrc%levels(1)%grid%xyz

  lb     = patch%lbound
  lbs    = patch%srcLbound
  align  = patch%align
  bcType = patch%bcType

! mapping between patches

  l1SrcDir = 1
  IF (patch%srcL1beg > patch%srcL1end) THEN
    l1SrcDir = -1
  ELSE IF (patch%srcL1beg == patch%srcL1end) THEN
    IF (lbs==1 .OR. lbs==2) THEN
      v1(:) = xyzSrc(:,IndIJK(ibegSrc,jendSrc,kbegSrc,iNOffSrc,ijNOffSrc)) - &
              xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc))
    ELSE IF (lbs==3 .OR. lbs==4) THEN
      v1(:) = xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kendSrc,iNOffSrc,ijNOffSrc)) - &
              xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc))
    ELSE IF (lbs==5 .OR. lbs==6) THEN
      v1(:) = xyzSrc(:,IndIJK(iendSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc)) - &
              xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc))
    ENDIF
    IF (lb==1 .OR. lb==2) THEN
      IF (align) THEN
        v2(:) = xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ELSE
        v2(:) = xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ENDIF
    ELSE IF (lb==3 .OR. lb==4) THEN
      IF (align) THEN
        v2(:) = xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ELSE
        v2(:) = xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ENDIF
    ELSE IF (lb==5 .OR. lb==6) THEN
      IF (align) THEN
        v2(:) = xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ELSE
        v2(:) = xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ENDIF
    ENDIF
    IF (DOT_PRODUCT(v1,v2) < 0._RFREAL) l1SrcDir = -1
  ENDIF

  l2SrcDir =  1
  IF (patch%srcL2beg > patch%srcL2end) THEN
    l2SrcDir = -1
  ELSE IF (patch%srcL2beg == patch%srcL2end) THEN
    IF (lbs==1 .OR. lbs==2) THEN
      v1(:) = xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kendSrc,iNOffSrc,ijNOffSrc)) - &
              xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc))
    ELSE IF (lbs==3 .OR. lbs==4) THEN
      v1(:) = xyzSrc(:,IndIJK(iendSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc)) - &
              xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc))
    ELSE IF (lbs==5 .OR. lbs==6) THEN
      v1(:) = xyzSrc(:,IndIJK(ibegSrc,jendSrc,kbegSrc,iNOffSrc,ijNOffSrc)) - &
              xyzSrc(:,IndIJK(ibegSrc,jbegSrc,kbegSrc,iNOffSrc,ijNOffSrc))
    ENDIF
    IF (lb==1 .OR. lb==2) THEN
      IF (align) THEN
        v2(:) = xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ELSE
        v2(:) = xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ENDIF
    ELSE IF (lb==3 .OR. lb==4) THEN
      IF (align) THEN
        v2(:) = xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ELSE
        v2(:) = xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ENDIF
    ELSE IF (lb==5 .OR. lb==6) THEN
      IF (align) THEN
        v2(:) = xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ELSE
        v2(:) = xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)) - &
                xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
      ENDIF
    ENDIF
    IF (DOT_PRODUCT(v1,v2) < 0._RFREAL) l2SrcDir = -1
  ENDIF

  ibegSrc = ibegSrc + idirSrc    ! shift by one node into the interior
  iendSrc = iendSrc + idirSrc    ! (otherwise would start at the boundary)
  jbegSrc = jbegSrc + jdirSrc
  jendSrc = jendSrc + jdirSrc
  kbegSrc = kbegSrc + kdirSrc
  kendSrc = kendSrc + kdirSrc

  CALL RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                             idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                             ibeg,iend,jbeg,jend,kbeg,kend, &
                             ibegSrc,iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc, &
                             mapMat )

! loop over dummy nodes of current patch

  DO idum=1,region%nDumCells

! - inter-region boundary

    IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ii      = i - idum*idir
            jj      = j - idum*jdir
            kk      = k - idum*kdir
            ijkN    = IndIJK(ii,jj,kk,iNOff,ijNOff)
            ijkNSrc = IndIJKMap(ii,jj,kk,mapMat,iNOffSrc,ijNOffSrc)
            xyz(XCOORD,ijkN) = xyzSrc(XCOORD,ijkNSrc)
            xyz(YCOORD,ijkN) = xyzSrc(YCOORD,ijkNSrc)
            xyz(ZCOORD,ijkN) = xyzSrc(ZCOORD,ijkNSrc)
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k

! - translational periodicity

    ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ii       = i - idum*idir
            jj       = j - idum*jdir
            kk       = k - idum*kdir
            ijkN     = IndIJK(ii,jj,kk,iNOff,ijNOff)
            ijkNB    = IndIJK(i ,j ,k ,iNOff,ijNOff)
            ijkNSrc  = IndIJKMap(ii,jj,kk,mapMat,iNOffSrc,ijNOffSrc)
            ijkNBSrc = IndIJKMap(i ,j ,k ,mapMat,iNOffSrc,ijNOffSrc)
            dx       = xyzSrc(XCOORD,ijkNSrc) - xyzSrc(XCOORD,ijkNBSrc)
            dy       = xyzSrc(YCOORD,ijkNSrc) - xyzSrc(YCOORD,ijkNBSrc)
            dz       = xyzSrc(ZCOORD,ijkNSrc) - xyzSrc(ZCOORD,ijkNBSrc)
            xyz(XCOORD,ijkN) = xyz(XCOORD,ijkNB) + dx
            xyz(YCOORD,ijkN) = xyz(YCOORD,ijkNB) + dy
            xyz(ZCOORD,ijkN) = xyz(ZCOORD,ijkNB) + dz
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k

! - rotational periodicity (around x-axis)

    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      cosa = COS(patch%periodAngle)
      sina = SIN(patch%periodAngle)
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ii      = i - idum*idir
            jj      = j - idum*jdir
            kk      = k - idum*kdir
            ijkN    = IndIJK(ii,jj,kk,iNOff,ijNOff)
            ijkNSrc = IndIJKMap(ii,jj,kk,mapMat,iNOffSrc,ijNOffSrc)
            xyz(XCOORD,ijkN) = xyzSrc(XCOORD,ijkNSrc)
            xyz(YCOORD,ijkN) = cosa*xyzSrc(YCOORD,ijkNSrc) - &
                               sina*xyzSrc(ZCOORD,ijkNSrc)
            xyz(ZCOORD,ijkN) = sina*xyzSrc(YCOORD,ijkNSrc) + &
                               cosa*xyzSrc(ZCOORD,ijkNSrc)
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDIF

  ENDDO   ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ExchangeGeometryCopy

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeGeometryCopy.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.11  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.6  2002/12/06 22:29:26  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
! Revision 1.5  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.2  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
!******************************************************************************







