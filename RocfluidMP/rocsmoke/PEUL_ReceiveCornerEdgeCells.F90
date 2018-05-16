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
! Purpose: receives values for edge and corner cells from an adjacent
!          region.
!
! Description: this is for the case if the other region is located
!              on a different processor.
!
! Input: regions = data of all regions
!        iReg    = current region.
!
! Output: region%levels%peul%cv = cv variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ReceiveCornerEdgeCells.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReceiveCornerEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER    :: regions(:)
  INTEGER,        INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: ir, iedge, icorner, i, j, k, ijk, iCv

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegSrc, ibuff, nCv, nDim, source, tag
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkC

  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: level
  TYPE(t_dCellTransf), POINTER :: rcvPeulEcCell

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_ReceiveCornerEdgeCells.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PEUL_ReceiveCornerEdgeCells',&
  'PEUL_ReceiveCornerEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev  =  regions(iReg)%currLevel
  level => regions(iReg)%levels(iLev)
  cv    => level%peul%cv
  nCv   =  level%peul%nCv

  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

! copy data from buffer to dummy cells

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (level%rcvPeulEcCells(ir)%nCells > 0) THEN

      rcvPeulEcCell => level%rcvPeulEcCells(ir)
      nDim       =  rcvPeulEcCell%nCells
      ibuff      =  0

#ifdef MPI
      source = regions(ir)%procid
      tag    = regions(iReg)%localNumber + PEUL_TAG_SHIFT
      CALL MPI_Recv( rcvPeulEcCell%buff,nCv*nDim,MPI_RFREAL, &
                     source,tag,global%mpiComm,status,global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! --- edges

      DO iedge=1,12
        IF (level%edgeCells(iedge)%interact) THEN
          CALL RFLO_GetEdgeCellsIndices( regions(iReg),iLev,iedge, &
                                         ibeg,iend,jbeg,jend,kbeg,kend )

          ijk = 0
          DO k=kbeg,kend
            DO j=jbeg,jend
              DO i=ibeg,iend
                ijk     = ijk + 1
                ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
                iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
                IF (iRegSrc == ir) THEN
                  ibuff = ibuff + 1
                  IF (level%edgeCells(iedge)%cells(ijk)%rotate) THEN
                    ! rotational periodicity
                  ELSE
                    DO iCv=1,nCv
                      cv(iCv,ijkC) = rcvPeulEcCell%buff(ibuff+(iCv-1)*nDim)
                    ENDDO ! iCv
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF  ! interact
      ENDDO    ! iedge

! --- corners

      DO icorner=1,8
        IF (level%cornerCells(icorner)%interact) THEN
          CALL RFLO_GetCornerCellsIndices( regions(iReg),iLev,icorner, &
                                           ibeg,iend,jbeg,jend,kbeg,kend )

          ijk = 0
          DO k=kbeg,kend
            DO j=jbeg,jend
              DO i=ibeg,iend
                ijk     = ijk + 1
                ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
                iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
                IF (iRegSrc == ir) THEN
                  ibuff = ibuff + 1
                  IF (level%cornerCells(icorner)%cells(ijk)%rotate) THEN
                    ! rotational periodicity
                  ELSE
                    DO iCv=1,nCv
                      cv(iCv,ijkC) = rcvPeulEcCell%buff(ibuff+(iCv-1)*nDim)
                    ENDDO ! iCv
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF  ! interact
      ENDDO    ! icorner

    ENDIF      ! some cells to receive
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ReceiveCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReceiveCornerEdgeCells.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:50  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
!******************************************************************************







