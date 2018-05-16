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
! Purpose: copy values to edge and corner cells of an adjacent region.
!
! Description: this is for the case if the other region is located
!              on the same processor.
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
! $Id: TURB_floRansExchCornEdgeCells.F90,v 1.5 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansExchCornEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iedge, icorner, ijk, i, j, k, l

! ... local variables
  INTEGER :: icell, iRegSrc, iLev, nCv, iCOff, ijCOff, ijkD
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend

  REAL(RFREAL), POINTER :: tcv(:,:), tcvSrc(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER  :: level

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'TURB_FloRansExchCornEdgeCells',&
  'TURB_floRansExchCornEdgeCells.F90' )

! get dimensions and pointers

  iLev = regions(iReg)%currLevel
  nCv  = regions(iReg)%turbInput%nCv

  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

  level => regions(iReg)%levels(iLev)
  tcv   => level%turb%cv

! edge cells

  DO iedge=1,12
    IF (level%edgeCells(iedge)%interact.eqv..true.) THEN
      CALL RFLO_GetEdgeCellsIndices( regions(iReg),iLev,iedge, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )
      ijk = 0
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijk     =  ijk + 1
            ijkD    =  IndIJK(i,j,k,iCOff,ijCOff)
            icell   =  level%edgeCells(iedge)%cells(ijk)%srcCell
            iRegSrc =  level%edgeCells(iedge)%cells(ijk)%srcRegion
            IF (iRegSrc > 0) THEN
              IF (regions(iRegSrc)%procid == global%myProcid) THEN
                tcvSrc => regions(iRegSrc)%levels(iLev)%turb%cv

                IF (level%edgeCells(iedge)%cells(ijk)%rotate .eqv. .true.) THEN
                  ! rotational periodicity
                ELSE
                  DO l=1,nCv
                    tcv(l,ijkD) = tcvSrc(l,icell)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDIF        ! interact
  ENDDO          ! iedge

! corner cells (currently not participating)

!  DO icorner=1,8
!    IF (level%cornerCells(icorner)%interact) THEN
!      CALL RFLO_GetCornerCellsIndices( regions(iReg),iLev,icorner, &
!                                       ibeg,iend,jbeg,jend,kbeg,kend )
!      ijk = 0
!      DO k=kbeg,kend
!        DO j=jbeg,jend
!          DO i=ibeg,iend
!            ijk     =  ijk + 1
!            ijkD    =  IndIJK(i,j,k,iCOff,ijCOff)
!            icell   =  level%cornerCells(icorner)%cells(ijk)%srcCell
!            iRegSrc =  level%cornerCells(icorner)%cells(ijk)%srcRegion
!            IF (iRegSrc > 0) THEN
!              IF (regions(iRegSrc)%procid == global%myProcid) THEN
!                tcvSrc => regions(iRegSrc)%levels(iLev)%turb%cv

!                IF (level%cornerCells(icorner)%cells(ijk)%rotate) THEN
!                  ! rotational periodicity
!                ELSE
!                  DO l=1,nCv
!                    tcv(l,ijkD) = tcvSrc(l,icell)
!                  ENDDO
!                ENDIF
!              ENDIF
!            ENDIF
!          ENDDO  ! i
!        ENDDO    ! j
!      ENDDO      ! k
!    ENDIF        ! interact
!  ENDDO          ! icorner

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloRansExchCornEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansExchCornEdgeCells.F90,v $
! Revision 1.5  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
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
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.1  2004/01/23 00:39:06  wasistho
! added communication routines for RaNS edge/corners
!
!
!******************************************************************************







