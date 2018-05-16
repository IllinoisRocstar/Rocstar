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
! Purpose: receive FLD radiation data to patch dummy cells from neighboring
!          region / periodic boundary located on a different processor.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: region%levels%radi = FLD radiation variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_floFlimRecvDummyVals.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimRecvDummyVals( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: idum, i, j, k, l, ijkBuff

! ... local variables
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: lb, ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkD, &
             n1, n2, nDim, source, tag, iLev, nCv

  REAL(RFREAL), POINTER :: rcv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_FloFlimRecvDummyVals',&
  'RADI_floFlimRecvDummyVals.F90' )

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

  n1   =  ABS(patch%l1end-patch%l1beg) + 1   ! here, dimensions of current
  n2   =  ABS(patch%l2end-patch%l2beg) + 1   ! and source patch are identical
  nDim =  n1*n2*region%nDumCells             ! ... but not the # of dummy cells

! receive data

#ifdef MPI
  source = regionSrc%procid
  tag    = region%localNumber + MPI_PATCHOFF*patchSrc%srcPatch* &
           RADI_TAG_SHIFT
  CALL MPI_Recv( patch%valRadi%recvBuff,nCv*nDim,MPI_RFREAL, &
                 source,tag,global%mpiComm,status,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! copy from buffer to dummy nodes ---------------------------------------------

  lb      = patch%lbound
  ijkBuff = 0

  DO idum=1,region%nDumCells

! - face i=const.

    IF (lb==1 .OR. lb==2) THEN
      IF (lb == 1) i = ibeg - idum
      IF (lb == 2) i = iend + idum
      DO k=kbeg,kend
        DO j=jbeg,jend
          ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
          ijkBuff = ijkBuff + 1
          DO l=1,nCv
            rcv(l,ijkD) = patch%valRadi%recvBuff(ijkBuff+(l-1)*nDim)
          ENDDO
        ENDDO
      ENDDO

! - face j=const.

    ELSE IF (lb==3 .OR. lb==4) THEN
      IF (lb == 3) j = jbeg - idum
      IF (lb == 4) j = jend + idum
      DO i=ibeg,iend
        DO k=kbeg,kend
          ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
          ijkBuff = ijkBuff + 1
          DO l=1,nCv
            rcv(l,ijkD) = patch%valRadi%recvBuff(ijkBuff+(l-1)*nDim)
          ENDDO
        ENDDO
      ENDDO

! - face k=const.

    ELSE IF (lb==5 .OR. lb==6) THEN
      IF (lb == 5) k = kbeg - idum
      IF (lb == 6) k = kend + idum
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
          ijkBuff = ijkBuff + 1
          DO l=1,nCv
            rcv(l,ijkD) = patch%valRadi%recvBuff(ijkBuff+(l-1)*nDim)
          ENDDO
        ENDDO
      ENDDO
    ENDIF   ! lb

  ENDDO     ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FloFlimRecvDummyVals

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimRecvDummyVals.F90,v $
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







