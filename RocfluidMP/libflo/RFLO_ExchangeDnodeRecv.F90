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
! Purpose: receive movement of interface nodes from adjacent region on
!          another processor.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc
!        average   = average source and current node movement (true/false).
!
! Output: dNode = grid movement of the interface.
!
! Notes: operation carried out for the finest grid only. Intended for
!        conforming grid boundaries (also periodic) only. In the case
!        of rotationally periodic boundary, it is assumed that the
!        axis of rotation coincides with the x-axis.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeDnodeRecv.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeDnodeRecv( region,regionSrc,patch,patchSrc, &
                                   average,dNode )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: average

  REAL(RFREAL), POINTER :: dNode(:,:)

  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: i, j, k, ijkBuff

! ... local variables
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: bcType, ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, &
             ijkN, ijkNB, iBnd, jBnd, kBnd, n1, n2, nDim, lb, source, tag

  REAL(RFREAL)          :: dx, dy, dz, dyr, dzr, cosa, sina
  REAL(RFREAL), POINTER :: xyz(:,:), dOld(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ExchangeDnodeRecv',&
  'RFLO_ExchangeDnodeRecv.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,&
    __LINE__ )
  ENDIF

! get dimensions and pointers

  CALL RFLO_GetPatchIndicesNodes( region,patch,1,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,1,iNOff,ijNOff )

  xyz  => region%levels(1)%gridOld%xyz
  dOld => region%levels(1)%grid%xyzOld

  n1   = ABS(patch%l1end-patch%l1beg) + 2   ! here, dimensions of current
  n2   = ABS(patch%l2end-patch%l2beg) + 2   ! and source patch are identical
  nDim = n1*n2

! receive data

#ifdef MPI
  source = regionSrc%procid
  tag    = region%localNumber + MPI_PATCHOFF*patchSrc%srcPatch
  CALL MPI_Recv( patch%mixt%recvBuff,3*nDim,MPI_RFREAL,source,tag, &
                 global%mpiComm,status,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
  __LINE__ )
#endif

! copy from buffer to dummy nodes ---------------------------------------------

  lb      = patch%lbound
  bcType  = patch%bcType
  ijkBuff = 0

! face i=const.

  IF (lb==1 .OR. lb==2) THEN

    IF (lb == 1) THEN
      i    = ibeg
      iBnd = ibeg
    ELSE IF (lb == 2) THEN
      i    = iend
      iBnd = iend
    ENDIF

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
      DO k=kbeg,kend
        DO j=jbeg,jend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkBuff = ijkBuff + 1
          dx      = patch%mixt%recvBuff(ijkBuff       )
          dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
          dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
          IF (average) THEN
            IF (ABS(dx-dOld(XCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(XCOORD,ijkN)-dOld(XCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)
            ELSE
              dNode(XCOORD,ijkN) = 0.5_RFREAL*(dNode(XCOORD,ijkN)+dx)
            ENDIF
            IF (ABS(dy-dOld(YCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(YCOORD,ijkN)-dOld(YCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            ELSE
              dNode(YCOORD,ijkN) = 0.5_RFREAL*(dNode(YCOORD,ijkN)+dy)
            ENDIF
            IF (ABS(dz-dOld(ZCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(ZCOORD,ijkN)-dOld(ZCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
            ELSE
              dNode(ZCOORD,ijkN) = 0.5_RFREAL*(dNode(ZCOORD,ijkN)+dz)
            ENDIF
          ELSE    ! no average
            IF (ABS(dx) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dx
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
        ENDDO
      ENDDO
    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      cosa = COS(patch%periodAngle)
      sina = SIN(patch%periodAngle)
      DO k=kbeg,kend
        DO j=jbeg,jend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkBuff = ijkBuff + 1
          dx      = patch%mixt%recvBuff(ijkBuff       )
          dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
          dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
          dyr     = cosa*dy - sina*dz
          dzr     = sina*dy + cosa*dz
          IF (average) THEN
            dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)
            dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
          ELSE    ! no average
            IF (ABS(dx ) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dx
            ENDIF
            IF (ABS(dyr) > ABS(dNode(YCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(YCOORD,ijkN) = dyr
            ENDIF
            IF (ABS(dzr) > ABS(dNode(ZCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(ZCOORD,ijkN) = dzr
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

! face j=const.

  ELSE IF (lb==3 .OR. lb==4) THEN

    IF (lb == 3) THEN
      j    = jbeg
      jBnd = jbeg
    ELSE IF (lb == 4) THEN
      j    = jend
      jBnd = jend
    ENDIF

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
      DO i=ibeg,iend
        DO k=kbeg,kend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkBuff = ijkBuff + 1
          dx      = patch%mixt%recvBuff(ijkBuff       )
          dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
          dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
          IF (average) THEN
            IF (ABS(dx-dOld(XCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(XCOORD,ijkN)-dOld(XCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)
            ELSE
              dNode(XCOORD,ijkN) = 0.5_RFREAL*(dNode(XCOORD,ijkN)+dx)
            ENDIF
            IF (ABS(dy-dOld(YCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(YCOORD,ijkN)-dOld(YCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            ELSE
              dNode(YCOORD,ijkN) = 0.5_RFREAL*(dNode(YCOORD,ijkN)+dy)
            ENDIF
            IF (ABS(dz-dOld(ZCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(ZCOORD,ijkN)-dOld(ZCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
            ELSE
              dNode(ZCOORD,ijkN) = 0.5_RFREAL*(dNode(ZCOORD,ijkN)+dz)
            ENDIF
          ELSE    ! no average
            IF (ABS(dx) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dx
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
        ENDDO
      ENDDO
    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      cosa = COS(patch%periodAngle)
      sina = SIN(patch%periodAngle)
      DO i=ibeg,iend
        DO k=kbeg,kend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkBuff = ijkBuff + 1
          dx      = patch%mixt%recvBuff(ijkBuff       )
          dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
          dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
          dyr     = cosa*dy - sina*dz
          dzr     = sina*dy + cosa*dz
          IF (average) THEN
            dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)
            dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
          ELSE    ! no average
            IF (ABS(dx ) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dx
            ENDIF
            IF (ABS(dyr) > ABS(dNode(YCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(YCOORD,ijkN) = dyr
            ENDIF
            IF (ABS(dzr) > ABS(dNode(ZCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(ZCOORD,ijkN) = dzr
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

! face k=const.

  ELSE IF (lb==5 .OR. lb==6) THEN

    IF (lb == 5) THEN
      k    = kbeg
      kBnd = kbeg
    ELSE IF (lb == 6) THEN
      k    = kend
      kBnd = kend
    ENDIF

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkBuff = ijkBuff + 1
          dx      = patch%mixt%recvBuff(ijkBuff       )
          dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
          dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
          IF (average) THEN
            IF (ABS(dx-dOld(XCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(XCOORD,ijkN)-dOld(XCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)
            ELSE
              dNode(XCOORD,ijkN) = 0.5_RFREAL*(dNode(XCOORD,ijkN)+dx)
            ENDIF
            IF (ABS(dy-dOld(YCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(YCOORD,ijkN)-dOld(YCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            ELSE
              dNode(YCOORD,ijkN) = 0.5_RFREAL*(dNode(YCOORD,ijkN)+dy)
            ENDIF
            IF (ABS(dz-dOld(ZCOORD,ijkN))<1.E-30_RFREAL .OR. &
                ABS(dNode(ZCOORD,ijkN)-dOld(ZCOORD,ijkN))<1.E-30_RFREAL) THEN
              dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
            ELSE
              dNode(ZCOORD,ijkN) = 0.5_RFREAL*(dNode(ZCOORD,ijkN)+dz)
            ENDIF
          ELSE    ! no average
            IF (ABS(dx) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dx
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
        ENDDO
      ENDDO
    ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
      cosa = COS(patch%periodAngle)
      sina = SIN(patch%periodAngle)
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          ijkBuff = ijkBuff + 1
          dx      = patch%mixt%recvBuff(ijkBuff       )
          dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
          dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
          dyr     = cosa*dy - sina*dz
          dzr     = sina*dy + cosa*dz
          IF (average) THEN
            dNode(XCOORD,ijkN) = dOld(XCOORD,ijkN)
            dNode(YCOORD,ijkN) = dOld(YCOORD,ijkN)
            dNode(ZCOORD,ijkN) = dOld(ZCOORD,ijkN)
          ELSE    ! no average
            IF (ABS(dx ) > ABS(dNode(XCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(XCOORD,ijkN) = dx
            ENDIF
            IF (ABS(dyr) > ABS(dNode(YCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(YCOORD,ijkN) = dyr
            ENDIF
            IF (ABS(dzr) > ABS(dNode(ZCOORD,ijkN))) THEN
              region%levels(1)%grid%boundMoved(lb) = .true.
              dNode(ZCOORD,ijkN) = dzr
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  ENDIF   ! lb

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExchangeDnodeRecv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeDnodeRecv.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:38:02  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.8  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/08/25 21:51:23  jblazek
! Full version of global grid motion scheme.
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







