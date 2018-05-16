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
! Purpose: send geometry of interior nodes to adjacent region
!          which is on another processor.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = region to send data to
!        patch     = current patch of region.
!
! Output: patch%mixt%sendBuff = send buffer.
!
! Notes: operation carried out for the finest grid only. Intended for
!        conforming grid boundaries (also periodic) only.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeGeometrySend.F90,v 1.4 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeGeometrySend( region,regionSrc,patch )

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
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, ijkBuff

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

  CALL RegisterFunction( global,'RFLO_ExchangeGeometrySend',&
  'RFLO_ExchangeGeometrySend.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers

  CALL RFLO_GetPatchIndicesNodes( region,patch,1,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,1,iNOff,ijNOff )

  xyz => region%levels(1)%grid%xyz

  n1   = ABS(patch%l1end-patch%l1beg) + 2   ! here, dimensions of current
  n2   = ABS(patch%l2end-patch%l2beg) + 2   ! and source patch are identical
  nDim = n1*n2*regionSrc%nDumCells          ! ... but not the # of dummy cells

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

  DO idum=1,regionSrc%nDumCells

! - face i=const.

    IF (lb==1 .OR. lb==2) THEN

      IF (lb == 1) THEN
        i    = ibeg + idum
        iBnd = ibeg
      ELSE IF (lb == 2) THEN
        i    = iend - idum
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
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          DO k=l2Beg,l2End,l2Step
            DO j=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN)
            ENDDO
          ENDDO
        ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
          DO k=l2Beg,l2End,l2Step
            DO j=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i   ,j,k,iNOff,ijNOff)
              ijkNB   = IndIJK(iBnd,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN ) - &
                                                       xyz(XCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN ) - &
                                                       xyz(YCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN ) - &
                                                       xyz(ZCOORD,ijkNB)
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
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          DO j=l2Beg,l2End,l2Step
            DO k=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN)
            ENDDO
          ENDDO
        ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
          DO j=l2Beg,l2End,l2Step
            DO k=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i   ,j,k,iNOff,ijNOff)
              ijkNB   = IndIJK(iBnd,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN ) - &
                                                       xyz(XCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN ) - &
                                                       xyz(YCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN ) - &
                                                       xyz(ZCOORD,ijkNB)
            ENDDO
          ENDDO
        ENDIF
      ENDIF   ! align

! - face j=const.

    ELSE IF (lb==3 .OR. lb==4) THEN

      IF (lb == 3) THEN
        j    = jbeg + idum
        jBnd = jbeg
      ELSE IF (lb == 4) THEN
        j    = jend - idum
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
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          DO i=l2Beg,l2End,l2Step
            DO k=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN)
            ENDDO
          ENDDO
        ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
          DO i=l2Beg,l2End,l2Step
            DO k=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j   ,k,iNOff,ijNOff)
              ijkNB   = IndIJK(i,jBnd,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN ) - &
                                                       xyz(XCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN ) - &
                                                       xyz(YCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN ) - &
                                                       xyz(ZCOORD,ijkNB)
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
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          DO k=l2Beg,l2End,l2Step
            DO i=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN)
            ENDDO
          ENDDO
        ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
          DO k=l2Beg,l2End,l2Step
            DO i=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j   ,k,iNOff,ijNOff)
              ijkNB   = IndIJK(i,jBnd,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN ) - &
                                                       xyz(XCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN ) - &
                                                       xyz(YCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN ) - &
                                                       xyz(ZCOORD,ijkNB)
            ENDDO
          ENDDO
        ENDIF
      ENDIF   ! align

! - face k=const.

    ELSE IF (lb==5 .OR. lb==6) THEN

      IF (lb == 5) THEN
        k    = kbeg + idum
        kBnd = kbeg
      ELSE IF (lb == 6) THEN
        k    = kend - idum
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
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          DO j=l2Beg,l2End,l2Step
            DO i=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN)
            ENDDO
          ENDDO
        ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
          DO j=l2Beg,l2End,l2Step
            DO i=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k   ,iNOff,ijNOff)
              ijkNB   = IndIJK(i,j,kBnd,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN ) - &
                                                       xyz(XCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN ) - &
                                                       xyz(YCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN ) - &
                                                       xyz(ZCOORD,ijkNB)
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
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          DO i=l2Beg,l2End,l2Step
            DO j=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN)
            ENDDO
          ENDDO
        ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
          DO i=l2Beg,l2End,l2Step
            DO j=l1Beg,l1End,l1Step
              ijkN    = IndIJK(i,j,k   ,iNOff,ijNOff)
              ijkNB   = IndIJK(i,j,kBnd,iNOff,ijNOff)
              ijkBuff = ijkBuff + 1
              patch%mixt%sendBuff(ijkBuff       ) = xyz(XCOORD,ijkN ) - &
                                                       xyz(XCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+  nDim) = xyz(YCOORD,ijkN ) - &
                                                       xyz(YCOORD,ijkNB)
              patch%mixt%sendBuff(ijkBuff+2*nDim) = xyz(ZCOORD,ijkN ) - &
                                                       xyz(ZCOORD,ijkNB)
            ENDDO
          ENDDO
        ENDIF
      ENDIF   ! align

    ENDIF     ! lb

  ENDDO       ! idum

! send data -------------------------------------------------------------------

#ifdef MPI
  dest = regionSrc%procid
  tag  = regionSrc%localNumber + MPI_PATCHOFF*patch%srcPatch
  CALL MPI_Isend( patch%mixt%sendBuff,3*nDim,MPI_RFREAL,dest,tag, &
                  global%mpiComm,global%requests(patch%mixt%iRequest), &
                  global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExchangeGeometrySend

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeGeometrySend.F90,v $
! Revision 1.4  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:39:39  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.14  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.10  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.9  2003/03/14 22:05:11  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.8  2002/12/06 22:29:26  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
! Revision 1.7  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.6  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.5  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.3  2002/03/30 00:50:49  jblazek
! Cleaned up with flint.
!
! Revision 1.2  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.1  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
!******************************************************************************







