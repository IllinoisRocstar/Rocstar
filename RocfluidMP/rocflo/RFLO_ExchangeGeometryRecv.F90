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
! Purpose: receive geometry of dummy nodes from adjacent region on
!          another processor.
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
! $Id: RFLO_ExchangeGeometryRecv.F90,v 1.4 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeGeometryRecv( region,regionSrc,patch,patchSrc )

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
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: idum, i, j, k, ijkBuff

! ... local variables
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: bcType, ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, &
             ijkN, ijkNB, iBnd, jBnd, kBnd, n1, n2, nDim, lb, source, tag

  REAL(RFREAL)          :: dx, dy, dz, cosa, sina
  REAL(RFREAL), POINTER :: xyz(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ExchangeGeometryRecv',&
  'RFLO_ExchangeGeometryRecv.F90' )

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
  nDim = n1*n2*region%nDumCells             ! ... but not the # of dummy cells

! receive data

#ifdef MPI
  source = regionSrc%procid
  tag    = region%localNumber + MPI_PATCHOFF*patchSrc%srcPatch
  CALL MPI_Recv( patch%mixt%recvBuff,3*nDim,MPI_RFREAL,source,tag, &
                 global%mpiComm,status,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! copy from buffer to dummy nodes ---------------------------------------------

  lb      = patch%lbound
  bcType  = patch%bcType
  ijkBuff = 0

  DO idum=1,region%nDumCells

! - face i=const.

    IF (lb==1 .OR. lb==2) THEN

      IF (lb == 1) THEN
        i    = ibeg - idum
        iBnd = ibeg
      ELSE IF (lb == 2) THEN
        i    = iend + idum
        iBnd = iend
      ENDIF

      IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
        DO k=kbeg,kend
          DO j=jbeg,jend
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            xyz(XCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff       )
            xyz(YCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff+  nDim)
            xyz(ZCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff+2*nDim)
          ENDDO
        ENDDO
      ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
        DO k=kbeg,kend
          DO j=jbeg,jend
            ijkN    = IndIJK(i   ,j,k,iNOff,ijNOff)
            ijkNB   = IndIJK(iBnd,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            dx      = patch%mixt%recvBuff(ijkBuff       )
            dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
            dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
            xyz(XCOORD,ijkN) = xyz(XCOORD,ijkNB) + dx
            xyz(YCOORD,ijkN) = xyz(YCOORD,ijkNB) + dy
            xyz(ZCOORD,ijkN) = xyz(ZCOORD,ijkNB) + dz
          ENDDO
        ENDDO
      ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
        cosa = COS(patch%periodAngle)
        sina = SIN(patch%periodAngle)
        DO k=kbeg,kend
          DO j=jbeg,jend
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
            dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
            xyz(XCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff)
            xyz(YCOORD,ijkN) = cosa*dy - sina*dz
            xyz(ZCOORD,ijkN) = sina*dy + cosa*dz
          ENDDO
        ENDDO
      ENDIF

! - face j=const.

    ELSE IF (lb==3 .OR. lb==4) THEN

      IF (lb == 3) THEN
        j    = jbeg - idum
        jBnd = jbeg
      ELSE IF (lb == 4) THEN
        j    = jend + idum
        jBnd = jend
      ENDIF

      IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
        DO i=ibeg,iend
          DO k=kbeg,kend
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            xyz(XCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff       )
            xyz(YCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff+  nDim)
            xyz(ZCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff+2*nDim)
          ENDDO
        ENDDO
      ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
        DO i=ibeg,iend
          DO k=kbeg,kend
            ijkN    = IndIJK(i,j   ,k,iNOff,ijNOff)
            ijkNB   = IndIJK(i,jBnd,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            dx      = patch%mixt%recvBuff(ijkBuff       )
            dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
            dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
            xyz(XCOORD,ijkN) = xyz(XCOORD,ijkNB) + dx
            xyz(YCOORD,ijkN) = xyz(YCOORD,ijkNB) + dy
            xyz(ZCOORD,ijkN) = xyz(ZCOORD,ijkNB) + dz
          ENDDO
        ENDDO
      ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
        cosa = COS(patch%periodAngle)
        sina = SIN(patch%periodAngle)
        DO i=ibeg,iend
          DO k=kbeg,kend
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
            dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
            xyz(XCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff)
            xyz(YCOORD,ijkN) = cosa*dy - sina*dz
            xyz(ZCOORD,ijkN) = sina*dy + cosa*dz
          ENDDO
        ENDDO
      ENDIF

! - face k=const.

    ELSE IF (lb==5 .OR. lb==6) THEN

      IF (lb == 5) THEN
        k    = kbeg - idum
        kBnd = kbeg
      ELSE IF (lb == 6) THEN
        k    = kend + idum
        kBnd = kend
      ENDIF

      IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            xyz(XCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff       )
            xyz(YCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff+  nDim)
            xyz(ZCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff+2*nDim)
          ENDDO
        ENDDO
      ELSE IF (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) THEN
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN    = IndIJK(i,j,k   ,iNOff,ijNOff)
            ijkNB   = IndIJK(i,j,kBnd,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            dx      = patch%mixt%recvBuff(ijkBuff       )
            dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
            dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
            xyz(XCOORD,ijkN) = xyz(XCOORD,ijkNB) + dx
            xyz(YCOORD,ijkN) = xyz(YCOORD,ijkNB) + dy
            xyz(ZCOORD,ijkN) = xyz(ZCOORD,ijkNB) + dz
          ENDDO
        ENDDO
      ELSE IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
        cosa = COS(patch%periodAngle)
        sina = SIN(patch%periodAngle)
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
            ijkBuff = ijkBuff + 1
            dy      = patch%mixt%recvBuff(ijkBuff+  nDim)
            dz      = patch%mixt%recvBuff(ijkBuff+2*nDim)
            xyz(XCOORD,ijkN) = patch%mixt%recvBuff(ijkBuff)
            xyz(YCOORD,ijkN) = cosa*dy - sina*dz
            xyz(ZCOORD,ijkN) = sina*dy + cosa*dz
          ENDDO
        ENDDO
      ENDIF

    ENDIF   ! lb

  ENDDO     ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExchangeGeometryRecv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeGeometryRecv.F90,v $
! Revision 1.4  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:39:37  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.13  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.8  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.7  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.6  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.4  2002/04/01 19:36:08  jblazek
! Added routine to clear send requests.
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







