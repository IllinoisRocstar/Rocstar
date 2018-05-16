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
! Output: new values of FLD radiation variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_floFlimRecvCornEdgeCells.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimRecvCornEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices
  USE ModError
  USE ModMPI
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: ir, iedge, icorner, i, j, k, l, ijk

! ... local variables
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, nCv, iRegSrc, ibuff, nDim, source, tag
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkC

  REAL(RFREAL), POINTER :: rcv(:,:)

  TYPE(t_global), POINTER      :: global
  TYPE(t_level), POINTER       :: level
  TYPE(t_dCellTransf), POINTER :: rcvRadiEcCell

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RADI_FloFlimRecvCornEdgeCells',&
  'RADI_floFlimRecvCornEdgeCells.F90' )

  iLev  =  regions(iReg)%currLevel
  nCv   =  regions(iReg)%radiInput%nCv

  level => regions(iReg)%levels(iLev)
  rcv   => level%radi%cv

  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

! copy data from buffer to dummy cells

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (level%rcvRadiEcCells(ir)%nCells > 0) THEN

      rcvRadiEcCell => level%rcvRadiEcCells(ir)
      nDim          =  rcvRadiEcCell%nCells
      ibuff         =  0

#ifdef MPI
      source = regions(ir)%procid
      tag    = regions(iReg)%localNumber + RADI_TAG_SHIFT
      CALL MPI_Recv( rcvRadiEcCell%buff,nCv*nDim,MPI_RFREAL, &
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
                    DO l=1,nCv
                      rcv(l,ijkC) = rcvRadiEcCell%buff(ibuff+(l-1)*nDim)
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF  ! interact
      ENDDO    ! iedge

! --- corners (currently not participating)

!      DO icorner=1,8
!        IF (level%cornerCells(icorner)%interact) THEN
!          CALL RFLO_GetCornerCellsIndices( regions(iReg),iLev,icorner, &
!                                           ibeg,iend,jbeg,jend,kbeg,kend )

!          ijk = 0
!          DO k=kbeg,kend
!            DO j=jbeg,jend
!              DO i=ibeg,iend
!                ijk     = ijk + 1
!                ijkC    = IndIJK(i,j,k,iCOff,ijCOff)
!                iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
!                IF (iRegSrc == ir) THEN
!                  ibuff = ibuff + 1
!                  IF (level%cornerCells(icorner)%cells(ijk)%rotate) THEN
!                    ! rotational periodicity
!                  ELSE
!                    DO l=1,nCv
!                      rcv(l,ijkC) = rcvRadiEcCell%buff(ibuff+(l-1)*nDim)
!                    ENDDO
!                  ENDIF
!                ENDIF
!              ENDDO
!            ENDDO
!          ENDDO

!        ENDIF  ! interact
!      ENDDO    ! icorner

    ENDIF      ! some cells to receive
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FloFlimRecvCornEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimRecvCornEdgeCells.F90,v $
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







