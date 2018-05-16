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
! Output: new values of RaNS variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_floRansSendCornEdgeCells.F90,v 1.6 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansSendCornEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: ir, iedge, icorner, ijk, l

! ... local variables
  INTEGER :: iLev, nCv, iRegSrc, icell, ibuff, nDim, dest, tag

  REAL(RFREAL), POINTER :: tcv(:,:)

  TYPE(t_global), POINTER      :: global
  TYPE(t_level), POINTER       :: level
  TYPE(t_dCellTransf), POINTER :: sndTurbEcCell

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'TURB_FloRansSendCornEdgeCells',&
  'TURB_floRansSendCornEdgeCells.F90' )

  iLev =  regions(iReg)%currLevel
  nCv  =  regions(iReg)%turbInput%nCv
  tcv  => regions(iReg)%levels(iLev)%turb%cv

! fill send buffers

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (regions(iReg)%levels(iLev)%sndTurbEcCells(ir)%nCells > 0) THEN

      sndTurbEcCell => regions(iReg)%levels(iLev)%sndTurbEcCells(ir)
      level         => regions(ir)%levels(iLev)
      nDim          =  sndTurbEcCell%nCells
      ibuff         =  0

      DO iedge=1,12
        IF (level%edgeCells(iedge)%interact.eqv..true.) THEN
          DO ijk=1,UBOUND(level%edgeCells(iedge)%cells,1)
            iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
            icell   = level%edgeCells(iedge)%cells(ijk)%srcCell
            IF (iRegSrc == iReg) THEN
              ibuff = ibuff + 1
              DO l=1,nCv
                sndTurbEcCell%buff(ibuff+(l-1)*nDim) = tcv(l,icell)
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
!                sndTurbEcCell%buff(ibuff+(l-1)*nDim) = tcv(l,icell)
!              ENDDO
!            ENDIF
!          ENDDO
!        ENDIF
!      ENDDO

#ifdef MPI
      dest = regions(ir)%procid
      tag  = regions(ir)%localNumber + TURB_TAG_SHIFT
      IF(tag .GT. global%mpiTagMax) tag = MOD(tag,global%mpiTagMax)
      CALL MPI_Isend( sndTurbEcCell%buff,nCv*nDim,MPI_RFREAL, &
                      dest,tag,global%mpiComm, &
                      global%requests(sndTurbEcCell%iRequest),global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

    ENDIF      ! some cells to send
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloRansSendCornEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansSendCornEdgeCells.F90,v $
! Revision 1.6  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.5  2009/04/07 15:00:29  mtcampbe
! Fixed possible tag errors.
!
! Revision 1.4  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.2  2004/02/18 22:50:30  wasistho
! added TURB_TAG_SHIFT in mpi_send tag
!
! Revision 1.1  2004/01/23 00:39:06  wasistho
! added communication routines for RaNS edge/corners
!
!
!******************************************************************************







