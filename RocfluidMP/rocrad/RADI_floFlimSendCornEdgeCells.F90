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
! Output: new values of FLD radiation variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_floFlimSendCornEdgeCells.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimSendCornEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModError
  USE ModMPI
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: ir, iedge, icorner, ijk, l

! ... local variables
  INTEGER :: iLev, nCv, iRegSrc, icell, ibuff, nDim, dest, tag

  REAL(RFREAL), POINTER :: rcv(:,:)

  TYPE(t_global), POINTER      :: global
  TYPE(t_level), POINTER       :: level
  TYPE(t_dCellTransf), POINTER :: sndRadiEcCell

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RADI_FloFlimSendCornEdgeCells',&
  'RADI_floFlimSendCornEdgeCells.F90' )

  iLev =  regions(iReg)%currLevel
  nCv  =  regions(iReg)%radiInput%nCv
  rcv  => regions(iReg)%levels(iLev)%radi%cv

! fill send buffers

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (regions(iReg)%levels(iLev)%sndRadiEcCells(ir)%nCells > 0) THEN

      sndRadiEcCell => regions(iReg)%levels(iLev)%sndRadiEcCells(ir)
      level         => regions(ir)%levels(iLev)
      nDim          =  sndRadiEcCell%nCells
      ibuff         =  0

      DO iedge=1,12
        IF (level%edgeCells(iedge)%interact) THEN
          DO ijk=1,UBOUND(level%edgeCells(iedge)%cells,1)
            iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
            icell   = level%edgeCells(iedge)%cells(ijk)%srcCell
            IF (iRegSrc == iReg) THEN
              ibuff = ibuff + 1
              DO l=1,nCv
                sndRadiEcCell%buff(ibuff+(l-1)*nDim) = rcv(l,icell)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO

! --- corners are not participating currently

!      DO icorner=1,8
!        IF (level%cornerCells(icorner)%interact) THEN
!          DO ijk=1,UBOUND(level%cornerCells(icorner)%cells,1)
!            iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
!            icell   = level%cornerCells(icorner)%cells(ijk)%srcCell
!            IF (iRegSrc == iReg) THEN
!              ibuff = ibuff + 1
!              DO l=1,nCv
!                sndRadiEcCell%buff(ibuff+(l-1)*nDim) = rcv(l,icell)
!              ENDDO
!            ENDIF
!          ENDDO
!        ENDIF
!      ENDDO

#ifdef MPI
      dest = regions(ir)%procid
      tag  = regions(ir)%localNumber + RADI_TAG_SHIFT
      CALL MPI_Isend( sndRadiEcCell%buff,nCv*nDim,MPI_RFREAL, &
                      dest,tag,global%mpiComm, &
                      global%requests(sndRadiEcCell%iRequest),global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

    ENDIF      ! some cells to send
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FloFlimSendCornEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimSendCornEdgeCells.F90,v $
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







