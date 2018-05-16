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
! Purpose: send values to edge and corner cells of an adjacent region.
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
! $Id: RFLO_SendCornerEdgeCells.F90,v 1.3 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_SendCornerEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: ir, iedge, icorner, ijk

! ... local variables
  INTEGER :: iLev, iRegSrc, icell, ibuff, nDim, dest, tag

  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global), POINTER      :: global
  TYPE(t_level), POINTER       :: level
  TYPE(t_dCellTransf), POINTER :: sendEcCell

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_SendCornerEdgeCells',&
  'RFLO_SendCornerEdgeCells.F90' )

  iLev =  regions(iReg)%currLevel
  cv   => regions(iReg)%levels(iLev)%mixt%cv

! fill send buffers

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (regions(iReg)%levels(iLev)%sendEcCells(ir)%nCells > 0) THEN

      sendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
      level      => regions(ir)%levels(iLev)
      nDim       =  sendEcCell%nCells
      ibuff      =  0

      DO iedge=1,12
        IF (level%edgeCells(iedge)%interact .AND. &
            level%edgeCells(iedge)%degenrt==DEGENERAT_NONE) THEN
          DO ijk=1,UBOUND(level%edgeCells(iedge)%cells,1)
            iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
            icell   = level%edgeCells(iedge)%cells(ijk)%srcCell
            IF (iRegSrc == iReg) THEN
              ibuff = ibuff + 1
              sendEcCell%buff(ibuff       ) = cv(CV_MIXT_DENS,icell)
              sendEcCell%buff(ibuff+  nDim) = cv(CV_MIXT_XMOM,icell)
              sendEcCell%buff(ibuff+2*nDim) = cv(CV_MIXT_YMOM,icell)
              sendEcCell%buff(ibuff+3*nDim) = cv(CV_MIXT_ZMOM,icell)
              sendEcCell%buff(ibuff+4*nDim) = cv(CV_MIXT_ENER,icell)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO icorner=1,8
        IF (level%cornerCells(icorner)%interact .AND. &
            level%cornerCells(icorner)%degenrt==DEGENERAT_NONE) THEN
          DO ijk=1,UBOUND(level%cornerCells(icorner)%cells,1)
            iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
            icell   = level%cornerCells(icorner)%cells(ijk)%srcCell
            IF (iRegSrc == iReg) THEN
              ibuff = ibuff + 1
              sendEcCell%buff(ibuff       ) = cv(CV_MIXT_DENS,icell)
              sendEcCell%buff(ibuff+  nDim) = cv(CV_MIXT_XMOM,icell)
              sendEcCell%buff(ibuff+2*nDim) = cv(CV_MIXT_YMOM,icell)
              sendEcCell%buff(ibuff+3*nDim) = cv(CV_MIXT_ZMOM,icell)
              sendEcCell%buff(ibuff+4*nDim) = cv(CV_MIXT_ENER,icell)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

#ifdef MPI
      dest = regions(ir)%procid
      tag  = regions(ir)%localNumber
      CALL MPI_Isend( sendEcCell%buff,CV_MIXT_NEQS*nDim,MPI_RFREAL, &
                      dest,tag,global%mpiComm, &
                      global%requests(sendEcCell%iRequest),global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

    ENDIF      ! some cells to send
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_SendCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_SendCornerEdgeCells.F90,v $
! Revision 1.3  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.9  2004/09/02 01:35:17  wasistho
! replaced 0 with DEGENERAT_NONE
!
! Revision 1.8  2004/09/02 01:23:32  wasistho
! do not exchange and send/recv for degenerated e/c
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







