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
! Purpose: receive RaNS data to the dummy cells of the patch from neighboring
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
! Output: region%levels%turb = RaNS variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_floRansRecvDummyVals.F90,v 1.7 2009/04/07 15:00:28 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansRecvDummyVals( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
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

  REAL(RFREAL), POINTER :: tcv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_FloRansRecvDummyVals',&
  'TURB_floRansRecvDummyVals.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  tcv  => region%levels(iLev)%turb%cv
  nCv  =  region%turbInput%nCv

  n1   =  ABS(patch%l1end-patch%l1beg) + 1   ! here, dimensions of current
  n2   =  ABS(patch%l2end-patch%l2beg) + 1   ! and source patch are identical
  nDim =  n1*n2*region%nDumCells             ! ... but not the # of dummy cells

! receive data

#ifdef MPI
  source = regionSrc%procid
  tag    = region%localNumber + MPI_PATCHOFF*patchSrc%srcPatch + &
           TURB_TAG_SHIFT
  IF(tag .GT. global%mpiTagMax) tag = MOD(tag,global%mpiTagMax)
  CALL MPI_Recv( patch%turb%recvBuff,nCv*nDim,MPI_RFREAL, &
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
            tcv(l,ijkD) = patch%turb%recvBuff(ijkBuff+(l-1)*nDim)
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
            tcv(l,ijkD) = patch%turb%recvBuff(ijkBuff+(l-1)*nDim)
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
            tcv(l,ijkD) = patch%turb%recvBuff(ijkBuff+(l-1)*nDim)
          ENDDO
        ENDDO
      ENDDO
    ENDIF   ! lb

  ENDDO     ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloRansRecvDummyVals

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansRecvDummyVals.F90,v $
! Revision 1.7  2009/04/07 15:00:28  mtcampbe
! Fixed possible tag errors.
!
! Revision 1.6  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.3  2006/08/19 15:40:47  mparmar
! Renamed patch variables
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







