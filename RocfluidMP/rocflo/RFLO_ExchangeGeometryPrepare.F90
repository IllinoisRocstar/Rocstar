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
! Purpose: exchange orientation of l1- and l2-directions between regions
!          on different processors.
!
! Description: none.
!
! Input: regions%levels%grid = coordinates (physical cells) and dimensions
!
! Output: regions%levels%patches%l1/2Vec = vectors in l1/2-direction.
!
! Notes: only region interfaces with continuous grid are treated here.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeGeometryPrepare.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeGeometryPrepare( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset, &
                            RFLO_ClearSendRequests
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER :: bcType, iRegSrc, iPatchSrc, lb, dest, source, tag
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif

  REAL(RFREAL) :: vec(6)
  REAL(RFREAL), POINTER :: xyz(:,:)

  TYPE(t_patch), POINTER  :: patch, patchSrc
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ExchangeGeometryPrepare',&
  'RFLO_ExchangeGeometryPrepare.F90' )

! send direction vectors ------------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

! --- loop over patches (finest grid level)

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType

! ----- conforming inter-region boundary, translational
!       and rotational periodicity

        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc = patch%srcRegion
          IF (regions(iRegSrc)%active == OFF) &
            CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )

          IF (regions(iRegSrc)%procid /= global%myProcid) THEN
            CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,1,ibeg,iend, &
                                            jbeg,jend,kbeg,kend )
            CALL RFLO_GetNodeOffset( regions(iReg),1,iNOff,ijNOff )
            lb  =  patch%lbound
            xyz => regions(iReg)%levels(1)%grid%xyz
            IF (lb==1 .OR. lb==2) THEN
              vec(1:3) = xyz(1:3,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)) - &
                         xyz(1:3,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
              vec(4:6) = xyz(1:3,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)) - &
                         xyz(1:3,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
            ELSE IF (lb==3 .OR. lb==4) THEN
              vec(1:3) = xyz(1:3,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)) - &
                         xyz(1:3,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
              vec(4:6) = xyz(1:3,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)) - &
                         xyz(1:3,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
            ELSE IF (lb==5 .OR. lb==6) THEN
              vec(1:3) = xyz(1:3,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)) - &
                         xyz(1:3,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
              vec(4:6) = xyz(1:3,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)) - &
                         xyz(1:3,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
            ENDIF
#ifdef MPI
            dest = regions(iRegSrc)%procid
            tag  = regions(iRegSrc)%localNumber + MPI_PATCHOFF*patch%srcPatch
            CALL MPI_Isend( vec,6,MPI_RFREAL,dest,tag,global%mpiComm, &
                            global%requests(patch%mixt%iRequest), &
                            global%mpierr )
            IF (global%mpierr /= 0) &
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! receive direction vectors ---------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

! --- loop over patches (finest grid level)

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType

! ----- conforming inter-region boundary, translational
!       or rotational periodicity

        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)

          IF (regions(iRegSrc)%procid /= global%myProcid) THEN
#ifdef MPI
            source = regions(iRegSrc)%procid
            tag    = regions(iReg)%localNumber + MPI_PATCHOFF*patchSrc%srcPatch
            CALL MPI_Recv( vec,6,MPI_RFREAL,source,tag,global%mpiComm, &
                           status,global%mpierr )
            IF (global%mpierr /= 0) &
              CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

            patch%l1VecSrc(1:3) = vec(1:3)
            patch%l2VecSrc(1:3) = vec(4:6)
#endif
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! clear send requests ---------------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      CALL RFLO_ClearSendRequests( regions,iReg,.true. )
    ENDIF
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExchangeGeometryPrepare

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeGeometryPrepare.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:39:35  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.1  2002/12/06 22:32:12  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
!******************************************************************************







