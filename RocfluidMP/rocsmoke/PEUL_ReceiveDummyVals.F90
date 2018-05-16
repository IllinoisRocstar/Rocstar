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
! Purpose: receive values to the dummy cells of the patch from neighboring
!          region / periodic boundary, which is located on a different
!          processor specific to PEUL Module.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: region%levels%peul = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ReceiveDummyVals.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReceiveDummyVals( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModPartEul, ONLY    : t_buffer_peul
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  TYPE(t_region), INTENT(IN)    :: regionSrc
  TYPE(t_patch),  INTENT(INOUT) :: patch
  TYPE(t_patch),  INTENT(IN)    :: patchSrc

! ... loop variables
  INTEGER :: iCv, idum, i, j, k, ijkBuff

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: lb, ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkD, &
             n1, n2, nCv, nDim, nDimRecvBuff, source, tagPeul, iLev

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pRecvBuffEul
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_ReceiveDummyVals.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_ReceiveDummyVals',&
  'PEUL_ReceiveDummyVals.F90' )

! check if the source region is active ----------------------------------------

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  pCv          => region%levels(iLev)%peul%cv
  pRecvBuffEul => patch%bufferPeul%recvBuff

  n1   = ABS(patch%l1end-patch%l1beg) + 1   ! here, dimensions of current
  n2   = ABS(patch%l2end-patch%l2beg) + 1   ! and source patch are identical
  nDim = n1*n2*region%nDumCells             ! ... but not the # of dummy cells

  nCv  =  region%levels(iLev)%peul%nCv

  IF (nCv /= region%peulInput%nPtypes) &
    CALL ErrorStop( region%global,ERR_PEUL_NPMISMATCH,__LINE__ )

  nDimRecvBuff = nCv*nDim

! receive data ----------------------------------------------------------------

#ifdef MPI
  source  = regionSrc%procid
  tagPeul = region%localNumber + MPI_PATCHOFF*patchSrc%srcPatch &
          * PEUL_TAG_SHIFT
  CALL MPI_Recv( pRecvBuffEul,nDimRecvBuff,MPI_RFREAL, &
                 source,tagPeul,global%mpiComm,status,global%mpierr )
  IF (global%mpierr /= ERR_NONE) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! copy from buffer to dummy nodes ---------------------------------------------

  lb      = patch%lbound
  ijkBuff = 0

  DO iCv = 1, nCv
    DO idum=1,region%nDumCells

! - face i=const. -------------------------------------------------------------

      IF (lb==1 .OR. lb==2) THEN

        SELECT CASE(lb)
          CASE(1)
            i = ibeg - idum
          CASE(2)
            i = iend + idum
        END SELECT ! lb

        DO k=kbeg,kend
        DO j=jbeg,jend
          ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
          ijkBuff = ijkBuff + 1
          pCv(iCv,ijkD) = pRecvBuffEul(ijkBuff)
        ENDDO ! j
        ENDDO ! k

! - face j=const. -------------------------------------------------------------

      ELSE IF (lb==3 .OR. lb==4) THEN

        SELECT CASE(lb)
          CASE(3)
            j = jbeg - idum
          CASE(4)
            j = jend + idum
        END SELECT ! lb

        DO i=ibeg,iend
        DO k=kbeg,kend
          ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
          ijkBuff = ijkBuff + 1
          pCv(iCv,ijkD) = pRecvBuffEul(ijkBuff)
        ENDDO ! k
        ENDDO ! i

! - face k=const.

      ELSE IF (lb==5 .OR. lb==6) THEN

        SELECT CASE(lb)
          CASE(5)
            k = kbeg - idum
          CASE(6)
            k = kend + idum
        END SELECT ! lb

        DO j=jbeg,jend
        DO i=ibeg,iend
          ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
          ijkBuff = ijkBuff + 1
          pCv(iCv,ijkD) = pRecvBuffEul(ijkBuff)
        ENDDO ! j
        ENDDO ! i
      ENDIF   ! lb

    ENDDO     ! idum
  ENDDO     ! iCv

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ReceiveDummyVals

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReceiveDummyVals.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:52  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/04/15 16:04:03  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.4  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/09 20:57:24  fnajjar
! Fixed sign error for ibeg, iend
!
! Revision 1.1  2003/04/09 14:34:25  fnajjar
! Initial Import of MPI-based rocsmoke
!
!******************************************************************************







