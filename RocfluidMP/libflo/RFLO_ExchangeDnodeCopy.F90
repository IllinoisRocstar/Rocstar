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
! Purpose: copy node movement from a patch of the source region to the nodes
!          of the current region`s patch (both regions are on the same
!          processor).
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc
!        average   = average source and current node movement (true/false)
!        dNode     = node movement of current region
!        dNodeSrc  = node movement of source region.
!
! Output: dNode = resulting node movement of current region.
!
! Notes: operation is carried out for the finest grid only. Intended for
!        conforming grid boundaries (also periodic) only. In the case
!        of rotationally periodic boundary, it is assumed that the
!        axis of rotation coincides with the x-axis.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeDnodeCopy.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeDnodeCopy( region,regionSrc,patch,patchSrc, &
                                   average,dNode,dNodeSrc )

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
  LOGICAL :: average

  REAL(RFREAL), POINTER :: dNode(:,:), dNodeSrc(:,:)

  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: bcType, ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, &
             iNOff, ijNOff, ijkN, ijkNB
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc, &
             idirSrc, jdirSrc, kdirSrc, iNOffSrc, ijNOffSrc, ijkNSrc, ijkNBSrc
  INTEGER :: lb, lbs, l1SrcDir, l2SrcDir, mapMat(3,4)

  LOGICAL :: align

  REAL(RFREAL) :: dy, dz, sina, cosa, v1(3), v2(3)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzSrc(:,:), dOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ExchangeDnodeCopy',&
  'RFLO_ExchangeDnodeCopy.F90' )

! check if the source region is active

  IF (regionSrc%active==OFF .AND. (.NOT.regionSrc%mixtInput%moveGrid)) THEN
    CALL ErrorStop( region%global,ERR_SRCREGION_OFF,&
    __LINE__ )
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

  xyz    => region%levels(1)%gridOld%xyz
  dOld   => region%levels(1)%grid%xyzOld
  xyzSrc => regionSrc%levels(1)%gridOld%xyz

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

! inter-region boundary or translational periodicity

  IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
      (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkNSrc = IndIJKMap(i,j,k,mapMat,iNOffSrc,ijNOffSrc)
          IF (average) THEN
            IF (ABS(dNodeSrc(XCOORD,ijkNSrc)-dOld(XCOORD,ijkN))<1.E-30_RFREAL &
              .OR. ABS(dNode(XCOORD,ijkN)-dOld(XCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)
            ELSE
              dNode(XCOORD,ijkN) = 0.5_RFREAL*(dNode(XCOORD,ijkN)+ &
                                               dNodeSrc(XCOORD,ijkNSrc))
              dNodeSrc(XCOORD,ijkNSrc) = dNode(XCOORD,ijkN)
            ENDIF
            IF (ABS(dNodeSrc(YCOORD,ijkNSrc)-dOld(YCOORD,ijkN))<1.E-30_RFREAL &
              .OR. ABS(dNode(YCOORD,ijkN)-dOld(YCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            ELSE
              dNode(YCOORD,ijkN) = 0.5_RFREAL*(dNode(YCOORD,ijkN)+ &
                                               dNodeSrc(YCOORD,ijkNSrc))
              dNodeSrc(YCOORD,ijkNSrc) = dNode(YCOORD,ijkN)
            ENDIF
            IF (ABS(dNodeSrc(ZCOORD,ijkNSrc)-dOld(ZCOORD,ijkN))<1.E-30_RFREAL &
              .OR. ABS(dNode(ZCOORD,ijkN)-dOld(ZCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
            ELSE
              dNode(ZCOORD,ijkN) = 0.5_RFREAL*(dNode(ZCOORD,ijkN)+ &
                                               dNodeSrc(ZCOORD,ijkNSrc))
              dNodeSrc(ZCOORD,ijkNSrc) = dNode(ZCOORD,ijkN)
            ENDIF
          ELSE    ! no average
            IF (ABS(dNodeSrc(XCOORD,ijkNSrc)) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dNodeSrc(XCOORD,ijkNSrc)
            ENDIF
            IF (ABS(dNodeSrc(YCOORD,ijkNSrc)) > ABS(dNode(YCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(YCOORD,ijkN) = dNodeSrc(YCOORD,ijkNSrc)
            ENDIF
            IF (ABS(dNodeSrc(ZCOORD,ijkNSrc)) > ABS(dNode(ZCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(ZCOORD,ijkN) = dNodeSrc(ZCOORD,ijkNSrc)
            ENDIF
          ENDIF
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k

! rotational periodicity (around x-axis)

  ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
    cosa = COS(patch%periodAngle)
    sina = SIN(patch%periodAngle)
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkNSrc = IndIJKMap(i,j,k,mapMat,iNOffSrc,ijNOffSrc)
          dy      = cosa*dNodeSrc(YCOORD,ijkNSrc) - &
                    sina*dNodeSrc(ZCOORD,ijkNSrc)
          dz      = sina*dNodeSrc(YCOORD,ijkNSrc) + &
                    cosa*dNodeSrc(ZCOORD,ijkNSrc)
          IF (average) THEN
            dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)  ! no extra movement here
            dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
          ELSE
            IF (ABS(dNodeSrc(XCOORD,ijkNSrc)) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dNodeSrc(XCOORD,ijkNSrc)
            ENDIF
            IF (ABS(dy) > ABS(dNode(YCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(YCOORD,ijkN) = dy
            ENDIF
            IF (ABS(dz) > ABS(dNode(ZCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(ZCOORD,ijkN) = dz
            ENDIF
          ENDIF
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k

  ENDIF

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ExchangeDnodeCopy

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeDnodeCopy.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.9  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/08/25 21:51:23  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.5  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.4  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.3  2003/05/06 20:05:38  jblazek
! Corrected bug in grid motion (corner "averaging").
!
! Revision 1.2  2003/03/17 23:48:47  jblazek
! Fixed bug in exchange of node displacements.
!
! Revision 1.1  2003/03/14 22:05:10  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
!******************************************************************************







