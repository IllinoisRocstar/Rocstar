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
! Purpose: send node movement of an interface to the adjacent region
!          which is on another processor.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = region to send data to
!        patch     = current patch of region
!        dNode     = node movement of current region.
!
! Output: patch%mixt%sendBuff = send buffer.
!
! Notes: operation carried out for the finest grid only. Intended for
!        conforming grid boundaries (also periodic) only.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeDnodeSend.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeDnodeSend( region,regionSrc,patch,dNode )

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
  REAL(RFREAL), POINTER :: dNode(:,:)

  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, ijkBuff

! ... local variables
  INTEGER :: bcType, ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, ijkN, &
             ijkNB, n1, n2, nDim, dest, tag, request, iBnd, jBnd, kBnd
  INTEGER :: lb, l1SrcDir, l2SrcDir, l1Beg, l1End, l1Step, l2Beg, l2End, l2Step

  LOGICAL :: align

  REAL(RFREAL) :: v1(3), v2(3)
  REAL(RFREAL), POINTER :: xyz(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ExchangeDnodeSend',&
  'RFLO_ExchangeDnodeSend.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,&
    __LINE__ )
  ENDIF

! get dimensions and pointers

  CALL RFLO_GetPatchIndicesNodes( region,patch,1,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,1,iNOff,ijNOff )

  xyz => region%levels(1)%gridOld%xyz

  n1   = ABS(patch%l1end-patch%l1beg) + 2   ! here, dimensions of current
  n2   = ABS(patch%l2end-patch%l2beg) + 2   ! and source patch are identical
  nDim = n1*n2

  lb     = patch%lbound
  align  = patch%align
  bcType = patch%bcType

! mapping between patches

  l1SrcDir = 1
  IF (patch%srcL1beg > patch%srcL1end) THEN
    l1SrcDir = -1
  ELSE IF (patch%srcL1beg == patch%srcL1end) THEN
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
    v1(:) = patch%l1VecSrc(:)
    IF (DOT_PRODUCT(v1,v2) < 0._RFREAL) l1SrcDir = -1
  ENDIF

  l2SrcDir =  1
  IF (patch%srcL2beg > patch%srcL2end) THEN
    l2SrcDir = -1
  ELSE IF (patch%srcL2beg == patch%srcL2end) THEN
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
    v1(:) = patch%l2VecSrc(:)
    IF (DOT_PRODUCT(v1,v2) < 0._RFREAL) l2SrcDir = -1
  ENDIF

! loop over interior nodes of current patch -----------------------------------

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

    IF (align) THEN
      IF (l1SrcDir > 0) THEN
        l1Beg = jbeg
        l1End = jend
      ELSE
        l1Beg = jend
        l1End = jbeg
      ENDIF
      l1Step = l1SrcDir
      IF (l2SrcDir > 0) THEN
        l2Beg = kbeg
        l2End = kend
      ELSE
        l2Beg = kend
        l2End = kbeg
      ENDIF
      l2Step = l2SrcDir
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
        DO k=l2Beg,l2End,l2Step
          DO j=l1Beg,l1End,l1Step
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            patch%mixt%sendBuff(ijkBuff       ) = dNode(XCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+  nDim) = dNode(YCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+2*nDim) = dNode(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDIF
    ELSE   ! not aligned
      IF (l1SrcDir > 0) THEN
        l1Beg = kbeg
        l1End = kend
      ELSE
        l1Beg = kend
        l1End = kbeg
      ENDIF
      l1Step = l1SrcDir
      IF (l2SrcDir > 0) THEN
        l2Beg = jbeg
        l2End = jend
      ELSE
        l2Beg = jend
        l2End = jbeg
      ENDIF
      l2Step = l2SrcDir
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
        DO j=l2Beg,l2End,l2Step
          DO k=l1Beg,l1End,l1Step
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            patch%mixt%sendBuff(ijkBuff       ) = dNode(XCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+  nDim) = dNode(YCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+2*nDim) = dNode(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDIF
    ENDIF   ! align

! face j=const.

  ELSE IF (lb==3 .OR. lb==4) THEN

    IF (lb == 3) THEN
      j    = jbeg
      jBnd = jbeg
    ELSE IF (lb == 4) THEN
      j    = jend
      jBnd = jend
    ENDIF

    IF (align) THEN
      IF (l1SrcDir > 0) THEN
        l1Beg = kbeg
        l1End = kend
      ELSE
        l1Beg = kend
        l1End = kbeg
      ENDIF
      l1Step = l1SrcDir
      IF (l2SrcDir > 0) THEN
        l2Beg = ibeg
        l2End = iend
      ELSE
        l2Beg = iend
        l2End = ibeg
      ENDIF
      l2Step = l2SrcDir
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
        DO i=l2Beg,l2End,l2Step
          DO k=l1Beg,l1End,l1Step
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            patch%mixt%sendBuff(ijkBuff       ) = dNode(XCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+  nDim) = dNode(YCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+2*nDim) = dNode(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDIF
    ELSE   ! not aligned
      IF (l1SrcDir > 0) THEN
        l1Beg = ibeg
        l1End = iend
      ELSE
        l1Beg = iend
        l1End = ibeg
      ENDIF
      l1Step = l1SrcDir
      IF (l2SrcDir > 0) THEN
        l2Beg = kbeg
        l2End = kend
      ELSE
        l2Beg = kend
        l2End = kbeg
      ENDIF
      l2Step = l2SrcDir
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
        DO k=l2Beg,l2End,l2Step
          DO i=l1Beg,l1End,l1Step
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            patch%mixt%sendBuff(ijkBuff       ) = dNode(XCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+  nDim) = dNode(YCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+2*nDim) = dNode(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDIF
    ENDIF   ! align

! face k=const.

  ELSE IF (lb==5 .OR. lb==6) THEN

    IF (lb == 5) THEN
      k    = kbeg
      kBnd = kbeg
    ELSE IF (lb == 6) THEN
      k    = kend
      kBnd = kend
    ENDIF

    IF (align) THEN
      IF (l1SrcDir > 0) THEN
        l1Beg = ibeg
        l1End = iend
      ELSE
        l1Beg = iend
        l1End = ibeg
      ENDIF
      l1Step = l1SrcDir
      IF (l2SrcDir > 0) THEN
        l2Beg = jbeg
        l2End = jend
      ELSE
        l2Beg = jend
        l2End = jbeg
      ENDIF
      l2Step = l2SrcDir
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
        DO j=l2Beg,l2End,l2Step
          DO i=l1Beg,l1End,l1Step
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            patch%mixt%sendBuff(ijkBuff       ) = dNode(XCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+  nDim) = dNode(YCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+2*nDim) = dNode(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDIF
    ELSE   ! not aligned
      IF (l1SrcDir > 0) THEN
        l1Beg = jbeg
        l1End = jend
      ELSE
        l1Beg = jend
        l1End = jbeg
      ENDIF
      l1Step = l1SrcDir
      IF (l2SrcDir > 0) THEN
        l2Beg = ibeg
        l2End = iend
      ELSE
        l2Beg = iend
        l2End = ibeg
      ENDIF
      l2Step = l2SrcDir
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE)) THEN
        DO i=l2Beg,l2End,l2Step
          DO j=l1Beg,l1End,l1Step
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            patch%mixt%sendBuff(ijkBuff       ) = dNode(XCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+  nDim) = dNode(YCOORD,ijkN)
            patch%mixt%sendBuff(ijkBuff+2*nDim) = dNode(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDIF
    ENDIF   ! align

  ENDIF     ! lb

! send data -------------------------------------------------------------------

#ifdef MPI
  dest = regionSrc%procid
  tag  = regionSrc%localNumber + MPI_PATCHOFF*patch%srcPatch
  CALL MPI_Isend( patch%mixt%sendBuff,3*nDim,MPI_RFREAL,dest,tag, &
                  global%mpiComm,global%requests(patch%mixt%iRequest), &
                  global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
  __LINE__ )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExchangeDnodeSend

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeDnodeSend.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:38:04  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.5  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.1  2003/03/14 22:05:10  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
!******************************************************************************







