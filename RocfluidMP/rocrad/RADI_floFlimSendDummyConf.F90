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
! Purpose: send FLD radiation data to dummy cells of the corresponding patch
!          of the adjacent region located on another processor.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = region to send data to
!        patch     = current patch of region.
!
! Output: patch%valRadi%nSendBuff = send buffer.
!
! Notes: intended for conforming grid boundaries only.
!
!******************************************************************************
!
! $Id: RADI_floFlimSendDummyConf.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimSendDummyConf( region,regionSrc,patch )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModBndPatch, ONLY   : t_patch
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, l, ijkBuff

! ... local variables
  INTEGER :: iLev, ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkC, &
             n1, n2, nDim, dest, tag
  INTEGER :: lb, l1SrcDir, l2SrcDir, l1Beg, l1End, l1Step, l2Beg, l2End, l2Step
  INTEGER :: nCv

  LOGICAL :: align

  REAL(RFREAL), POINTER :: rcv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_FloFlimSendDummyConf',&
  'RADI_floFlimSendDummyConf.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  rcv  => region%levels(iLev)%radi%cv
  nCv  =  region%radiInput%nCv

  n1   = ABS(patch%l1end-patch%l1beg) + 1   ! here, dimensions of current
  n2   = ABS(patch%l2end-patch%l2beg) + 1   ! and source patch are identical
  nDim = n1*n2*regionSrc%nDumCells          ! ... but not the # of dummy cells

! mapping between patches

                                       l1SrcDir =  1
  IF (patch%srcL1beg > patch%srcL1end) l1SrcDir = -1
                                       l2SrcDir =  1
  IF (patch%srcL2beg > patch%srcL2end) l2SrcDir = -1

  lb    = patch%lbound
  align = patch%align

! loop over interior cells of current patch -----------------------------------

  ijkBuff = 0

  DO idum=0,regionSrc%nDumCells-1

! - face i=const.

    IF (lb==1 .OR. lb==2) THEN

      IF (lb == 1) i = ibeg + idum
      IF (lb == 2) i = iend - idum
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
        DO k=l2Beg,l2End,l2Step
          DO j=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            DO l=1,nCv
              patch%valRadi%sendBuff(ijkBuff+(l-1)*nDim) = rcv(l,ijkC)
            ENDDO
          ENDDO
        ENDDO
      ELSE
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
        DO j=l2Beg,l2End,l2Step
          DO k=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            DO l=1,nCv
              patch%valRadi%sendBuff(ijkBuff+(l-1)*nDim) = rcv(l,ijkC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! - face j=const.

    ELSE IF (lb==3 .OR. lb==4) THEN

      IF (lb == 3) j = jbeg + idum
      IF (lb == 4) j = jend - idum
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
        DO i=l2Beg,l2End,l2Step
          DO k=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            DO l=1,nCv
              patch%valRadi%sendBuff(ijkBuff+(l-1)*nDim) = rcv(l,ijkC)
            ENDDO
          ENDDO
        ENDDO
      ELSE
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
        DO k=l2Beg,l2End,l2Step
          DO i=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            DO l=1,nCv
              patch%valRadi%sendBuff(ijkBuff+(l-1)*nDim) = rcv(l,ijkC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! - face k=const.

    ELSE IF (lb==5 .OR. lb==6) THEN

      IF (lb == 5) k = kbeg + idum
      IF (lb == 6) k = kend - idum
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
        DO j=l2Beg,l2End,l2Step
          DO i=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            DO l=1,nCv
              patch%valRadi%sendBuff(ijkBuff+(l-1)*nDim) = rcv(l,ijkC)
            ENDDO
          ENDDO
        ENDDO
      ELSE
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
        DO i=l2Beg,l2End,l2Step
          DO j=l1Beg,l1End,l1Step
            ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
            ijkBuff = ijkBuff + 1
            DO l=1,nCv
              patch%valRadi%sendBuff(ijkBuff+(l-1)*nDim) = rcv(l,ijkC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

    ENDIF    ! lb

  ENDDO      ! idum

! send data

#ifdef MPI
  dest = regionSrc%procid
  tag  = regionSrc%localNumber + MPI_PATCHOFF*patch%srcPatch* &
         RADI_TAG_SHIFT
  CALL MPI_Isend( patch%valRadi%sendBuff,nCv*nDim,MPI_RFREAL, &
                  dest,tag,global%mpiComm, &
                  global%requests(patch%valRadi%iRequest),global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FloFlimSendDummyConf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimSendDummyConf.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







