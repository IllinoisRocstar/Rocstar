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
!          processor.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: region%levels%mixt = flow variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReceiveDummyVals.F90,v 1.5 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReceiveDummyVals( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset, &
                            MixtureProperties
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: idum, i, j, k, ijkBuff

! ... local variables
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: lb, ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkD, &
             n1, n2, nDim, source, tag, iLev, gasModel

  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ReceiveDummyVals',&
  'RFLO_ReceiveDummyVals.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  gasModel = region%mixtInput%gasModel

  cv => region%levels(iLev)%mixt%cv

  n1   = ABS(patch%l1end-patch%l1beg) + 1   ! here, dimensions of current
  n2   = ABS(patch%l2end-patch%l2beg) + 1   ! and source patch are identical
  nDim = n1*n2*region%nDumCells             ! ... but not the # of dummy cells

! receive data

#ifdef MPI
  source = regionSrc%procid
  tag    = region%localNumber + MPI_PATCHOFF*patchSrc%srcPatch
  CALL MPI_Recv( patch%mixt%recvBuff,CV_MIXT_NEQS*nDim,MPI_RFREAL, &
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
          cv(CV_MIXT_DENS,ijkD) = patch%mixt%recvBuff(ijkBuff       )
          cv(CV_MIXT_XMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+  nDim)
          cv(CV_MIXT_YMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+2*nDim)
          cv(CV_MIXT_ZMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+3*nDim)
          cv(CV_MIXT_ENER,ijkD) = patch%mixt%recvBuff(ijkBuff+4*nDim)
          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF
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
          cv(CV_MIXT_DENS,ijkD) = patch%mixt%recvBuff(ijkBuff       )
          cv(CV_MIXT_XMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+  nDim)
          cv(CV_MIXT_YMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+2*nDim)
          cv(CV_MIXT_ZMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+3*nDim)
          cv(CV_MIXT_ENER,ijkD) = patch%mixt%recvBuff(ijkBuff+4*nDim)
          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF
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
          cv(CV_MIXT_DENS,ijkD) = patch%mixt%recvBuff(ijkBuff       )
          cv(CV_MIXT_XMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+  nDim)
          cv(CV_MIXT_YMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+2*nDim)
          cv(CV_MIXT_ZMOM,ijkD) = patch%mixt%recvBuff(ijkBuff+3*nDim)
          cv(CV_MIXT_ENER,ijkD) = patch%mixt%recvBuff(ijkBuff+4*nDim)
          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF
        ENDDO
      ENDDO
    ENDIF   ! lb

  ENDDO     ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReceiveDummyVals

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReceiveDummyVals.F90,v $
! Revision 1.5  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:39:40  mparmar
! Renamed patch variables
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.16  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.12  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.11  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.10  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.9  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.8  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.7  2002/04/01 19:36:08  jblazek
! Added routine to clear send requests.
!
! Revision 1.6  2002/03/30 00:50:49  jblazek
! Cleaned up with flint.
!
! Revision 1.5  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.4  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.1  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
!******************************************************************************







