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
! Purpose: send smoke values to edge and corner cells of an adjacent region.
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
! $Id: PEUL_SendCornerEdgeCells.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_SendCornerEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER    :: regions(:)
  INTEGER,        INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: ir, iedge, icorner, ijk, iCv

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, iRegSrc, icell, ibuff, nCv, nDim, dest, tag

  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: level
  TYPE(t_dCellTransf), POINTER :: sndPeulEcCell

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_SendCornerEdgeCells.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PEUL_SendCornerEdgeCells',&
  'PEUL_SendCornerEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev =  regions(iReg)%currLevel
  cv   => regions(iReg)%levels(iLev)%peul%cv
  nCv  =  regions(iReg)%levels(iLev)%peul%nCv

! fill send buffers

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (regions(iReg)%levels(iLev)%sndPeulEcCells(ir)%nCells > 0) THEN

      sndPeulEcCell => regions(iReg)%levels(iLev)%sndPeulEcCells(ir)
      level         => regions(ir)%levels(iLev)
      nDim          =  sndPeulEcCell%nCells
      ibuff         =  0

      DO iedge=1,12
        IF (level%edgeCells(iedge)%interact) THEN
          DO ijk=1,UBOUND(level%edgeCells(iedge)%cells,1)
            iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
            icell   = level%edgeCells(iedge)%cells(ijk)%srcCell
            IF (iRegSrc == iReg) THEN
              ibuff = ibuff + 1
              DO iCv=1,nCv
                sndPeulEcCell%buff(ibuff+(iCv-1)*nDim) = cv(iCv,icell)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO icorner=1,8
        IF (level%cornerCells(icorner)%interact) THEN
          DO ijk=1,UBOUND(level%cornerCells(icorner)%cells,1)
            iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
            icell   = level%cornerCells(icorner)%cells(ijk)%srcCell
            IF (iRegSrc == iReg) THEN
              ibuff = ibuff + 1
              DO iCv=1,nCv
                sndPeulEcCell%buff(ibuff+(iCv-1)*nDim) = cv(iCv,icell)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO

#ifdef MPI
      dest = regions(ir)%procid
      tag  = regions(ir)%localNumber + PEUL_TAG_SHIFT
      CALL MPI_Isend( sndPeulEcCell%buff,nCv*nDim,MPI_RFREAL, &
                      dest,tag,global%mpiComm, &
                      global%requests(sndPeulEcCell%iRequest),global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

    ENDIF      ! some cells to send
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_SendCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_SendCornerEdgeCells.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:55  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
!******************************************************************************







