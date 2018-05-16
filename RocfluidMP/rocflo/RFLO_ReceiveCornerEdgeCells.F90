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
! Output: new values of conservative variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReceiveCornerEdgeCells.F90,v 1.4 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReceiveCornerEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices, MixtureProperties
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: ir, iedge, icorner, i, j, k, ijk

! ... local variables
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, iRegSrc, ibuff, nDim, source, tag, gasModel
   INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff, ijkC

  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global), POINTER      :: global
  TYPE(t_level), POINTER       :: level
  TYPE(t_dCellTransf), POINTER :: recvEcCell

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_ReceiveCornerEdgeCells',&
  'RFLO_ReceiveCornerEdgeCells.F90' )

  iLev      =  regions(iReg)%currLevel
  gasModel =  regions(iReg)%mixtInput%gasModel
  level     => regions(iReg)%levels(iLev)
  cv        => level%mixt%cv

  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

! copy data from buffer to dummy cells

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (level%recvEcCells(ir)%nCells > 0) THEN

      recvEcCell => level%recvEcCells(ir)
      nDim       =  recvEcCell%nCells
      ibuff      =  0

#ifdef MPI
      source = regions(ir)%procid
      tag    = regions(iReg)%localNumber
      CALL MPI_Recv( recvEcCell%buff,CV_MIXT_NEQS*nDim,MPI_RFREAL, &
                     source,tag,global%mpiComm,status,global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! --- edges

      DO iedge=1,12
        IF (level%edgeCells(iedge)%interact .AND. &
            level%edgeCells(iedge)%degenrt==DEGENERAT_NONE) THEN
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
                    cv(CV_MIXT_DENS,ijkC) = recvEcCell%buff(ibuff       )
                    cv(CV_MIXT_XMOM,ijkC) = recvEcCell%buff(ibuff+  nDim)
                    cv(CV_MIXT_YMOM,ijkC) = recvEcCell%buff(ibuff+2*nDim)
                    cv(CV_MIXT_ZMOM,ijkC) = recvEcCell%buff(ibuff+3*nDim)
                    cv(CV_MIXT_ENER,ijkC) = recvEcCell%buff(ibuff+4*nDim)
                  ENDIF

                  IF (gasModel == GAS_MODEL_TCPERF) THEN
                    CALL MixtureProperties( regions(iReg),ijkC,ijkC,.false. )
                  ELSE
                    CALL MixtureProperties( regions(iReg),ijkC,ijkC,.true.  )
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF  ! interact
      ENDDO    ! iedge

! --- corners

      DO icorner=1,8
        IF (level%cornerCells(icorner)%interact .AND. &
            level%cornerCells(icorner)%degenrt==DEGENERAT_NONE) THEN
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
                    cv(CV_MIXT_DENS,ijkC) = recvEcCell%buff(ibuff       )
                    cv(CV_MIXT_XMOM,ijkC) = recvEcCell%buff(ibuff+  nDim)
                    cv(CV_MIXT_YMOM,ijkC) = recvEcCell%buff(ibuff+2*nDim)
                    cv(CV_MIXT_ZMOM,ijkC) = recvEcCell%buff(ibuff+3*nDim)
                    cv(CV_MIXT_ENER,ijkC) = recvEcCell%buff(ibuff+4*nDim)
                  ENDIF

                  IF (gasModel == GAS_MODEL_TCPERF) THEN
                    CALL MixtureProperties( regions(iReg),ijkC,ijkC,.false. )
                  ELSE
                    CALL MixtureProperties( regions(iReg),ijkC,ijkC,.true.  )
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

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReceiveCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReceiveCornerEdgeCells.F90,v $
! Revision 1.4  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.13  2004/09/02 01:35:10  wasistho
! replaced 0 with DEGENERAT_NONE
!
! Revision 1.12  2004/09/02 01:23:13  wasistho
! do not exchange and send/recv for degenerated e/c
!
! Revision 1.11  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.6  2003/02/28 21:04:27  jblazek
! Corrected bug in send/recv of CE cells on single processor.
!
! Revision 1.5  2003/02/14 22:32:37  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.4  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/31 20:20:13  jblazek
! Added treatment of edge & corner cells.
!
!******************************************************************************







