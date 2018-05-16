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
! Purpose: send values to dummy cells of the corresponding patch
!          of the adjacent region which is on another processor.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = region to send data to
!        patch     = current patch of region.
!
! Output: patch%peul%SendBuff = send buffer.
!
! Notes: intended for conforming grid boundaries only.
!
!******************************************************************************
!
! $Id: PEUL_SendDummyConf.F90,v 1.4 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_SendDummyConf( region,regionSrc,patch )

  USE ModDataTypes
  USE ModPartEul, ONLY    : t_buffer_peul
  USE ModBndPatch, ONLY   : t_patch
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset
  USE ModPartEul, ONLY    : t_peul
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region,regionSrc
  TYPE(t_patch),  INTENT(INOUT) :: patch

! ... loop variables
  INTEGER :: iCv, idum, i, j, k, ijkBuff

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkC, &
             n1, n2, nDim, dest, tagPeul
  INTEGER :: lb, l1SrcDir, l2SrcDir, l1Beg, l1End, l1Step, l2Beg, l2End, l2Step

  INTEGER :: iRequestPeul, nCv, nDimSendBuff

  LOGICAL :: align

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pSendBuffEul
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv

  TYPE(t_peul),        POINTER :: pPeul
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_SendDummyConf.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction( global,'RFLO_SendDummyConf',&
  'PEUL_SendDummyConf.F90' )

! check if the source region is active ----------------------------------------

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  pPeul        => region%levels(iLev)%peul
  pCv          => pPeul%cv
  pSendBuffEul => patch%bufferPeul%sendBuff

  n1   = ABS(patch%l1end-patch%l1beg) + 1   ! here, dimensions of current
  n2   = ABS(patch%l2end-patch%l2beg) + 1   ! and source patch are identical
  nDim = n1*n2*regionSrc%nDumCells          ! ... but not the # of dummy cells

  nCv  =  region%levels(iLev)%peul%nCv

  IF (nCv /= region%peulInput%nPtypes) &
    CALL ErrorStop( region%global,ERR_PEUL_NPMISMATCH,__LINE__ )

  nDimSendBuff = nCv*nDim

! mapping between patches -----------------------------------------------------

                                       l1SrcDir =  1
  IF (patch%srcL1beg > patch%srcL1end) l1SrcDir = -1
                                       l2SrcDir =  1
  IF (patch%srcL2beg > patch%srcL2end) l2SrcDir = -1

  lb    = patch%lbound
  align = patch%align

! loop over interior cells of current patch -----------------------------------

  ijkBuff = 0

  DO iCv = 1, nCv
    DO idum=0,regionSrc%nDumCells-1

! - face i=const. -------------------------------------------------------------

      IF (lb==1 .OR. lb==2) THEN

        SELECT CASE(lb)
          CASE(1)
            i = ibeg + idum
          CASE(2)
            i = iend - idum
        END SELECT ! lb

        IF (align) THEN
          IF (l1SrcDir > 0) THEN
            l1Beg = jbeg
            l1End = jend
          ELSE
            l1Beg = jend
            l1End = jbeg
          ENDIF ! l1SrcDir
          l1Step = l1SrcDir

          IF (l2SrcDir > 0) THEN
            l2Beg = kbeg
            l2End = kend
          ELSE
            l2Beg = kend
            l2End = kbeg
          ENDIF ! l2SrcDir
          l2Step = l2SrcDir

          DO k=l2Beg,l2End,l2Step
          DO j=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            pSendBuffEul(ijkBuff) = pCv(iCv,ijkC)
          ENDDO ! j
          ENDDO ! k

        ELSE ! align
          IF (l1SrcDir > 0) THEN
            l1Beg = kbeg
            l1End = kend
          ELSE
            l1Beg = kend
            l1End = kbeg
          ENDIF ! l1SrcDir
          l1Step = l1SrcDir

          IF (l2SrcDir > 0) THEN
            l2Beg = jbeg
            l2End = jend
          ELSE
            l2Beg = jend
            l2End = jbeg
          ENDIF ! l2SrcDir
          l2Step = l2SrcDir

          DO j=l2Beg,l2End,l2Step
          DO k=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            pSendBuffEul(ijkBuff) = pCv(iCv,ijkC)
          ENDDO ! k
          ENDDO ! j

         ENDIF ! align

! - face j=const. -------------------------------------------------------------

      ELSE IF (lb==3 .OR. lb==4) THEN

        SELECT CASE(lb)
          CASE(3)
            j = jbeg + idum
          CASE(4)
            j = jend - idum
        END SELECT ! lb

        IF (align) THEN
          IF (l1SrcDir > 0) THEN
            l1Beg = kbeg
            l1End = kend
          ELSE
            l1Beg = kend
            l1End = kbeg
          ENDIF ! l1SrcDir
          l1Step = l1SrcDir

          IF (l2SrcDir > 0) THEN
            l2Beg = ibeg
            l2End = iend
          ELSE
            l2Beg = iend
            l2End = ibeg
          ENDIF ! l2SrcDir
          l2Step = l2SrcDir

          DO i=l2Beg,l2End,l2Step
          DO k=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            pSendBuffEul(ijkBuff) = pCv(iCv,ijkC)
          ENDDO ! k
          ENDDO ! i

        ELSE ! align
          IF (l1SrcDir > 0) THEN
            l1Beg = ibeg
            l1End = iend
          ELSE
            l1Beg = iend
            l1End = ibeg
          ENDIF ! l1SrcDir
          l1Step = l1SrcDir

          IF (l2SrcDir > 0) THEN
            l2Beg = kbeg
            l2End = kend
          ELSE
            l2Beg = kend
            l2End = kbeg
          ENDIF ! l2SrcDir
          l2Step = l2SrcDir

          DO k=l2Beg,l2End,l2Step
          DO i=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            pSendBuffEul(ijkBuff) = pCv(iCv,ijkC)
          ENDDO ! i
          ENDDO ! k

        ENDIF ! align

! - face k=const. -------------------------------------------------------------

      ELSE IF (lb==5 .OR. lb==6) THEN

        SELECT CASE(lb)
          CASE(5)
            k = kbeg + idum
          CASE(6)
            k = kend - idum
        END SELECT ! lb

        IF (align) THEN
          IF (l1SrcDir > 0) THEN
            l1Beg = ibeg
            l1End = iend
          ELSE
            l1Beg = iend
            l1End = ibeg
          ENDIF ! l1SrcDir
          l1Step = l1SrcDir

          IF (l2SrcDir > 0) THEN
            l2Beg = jbeg
            l2End = jend
          ELSE
            l2Beg = jend
            l2End = jbeg
          ENDIF ! l2SrcDir
          l2Step = l2SrcDir

          DO j=l2Beg,l2End,l2Step
          DO i=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            pSendBuffEul(ijkBuff) = pCv(iCv,ijkC)
          ENDDO ! i
          ENDDO ! j

        ELSE ! align
          IF (l1SrcDir > 0) THEN
            l1Beg = jbeg
            l1End = jend
          ELSE
            l1Beg = jend
            l1End = jbeg
          ENDIF ! l1SrcDir
          l1Step = l1SrcDir

          IF (l2SrcDir > 0) THEN
            l2Beg = ibeg
            l2End = iend
          ELSE
            l2Beg = iend
            l2End = ibeg
          ENDIF ! l2SrcDir
          l2Step = l2SrcDir

          DO i=l2Beg,l2End,l2Step
          DO j=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            pSendBuffEul(ijkBuff) = pCv(iCv,ijkC)
          ENDDO ! j
          ENDDO ! i
        ENDIF ! align

      ENDIF    ! lb

    ENDDO      ! idum
  ENDDO        ! iCv

! send data

#ifdef MPI
  dest    = regionSrc%procid
  tagPeul = regionSrc%localNumber + MPI_PATCHOFF*patch%srcPatch &
          * PEUL_TAG_SHIFT
  iRequestPeul = patch%bufferPeul%iRequest
  CALL MPI_Isend( pSendBuffEul,nDimSendBuff,MPI_RFREAL,      &
                  dest,tagPeul,global%mpiComm,               &
                  pPeul%requests(iRequestPeul),global%mpierr )
  IF (global%mpierr /= ERR_NONE) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_SendDummyConf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_SendDummyConf.F90,v $
! Revision 1.4  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:26  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/12/01 21:09:56  haselbac
! Initial revision after changing case
!
! Revision 1.4  2004/04/15 16:04:03  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.3  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1  2003/04/09 14:34:25  fnajjar
! Initial Import of MPI-based rocsmoke
!
!******************************************************************************







