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
! Purpose: copy smoke values to edge and corner cells of an adjacent region.
!
! Description: this is for the case if the other region is located
!              on the same processor.
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
! $Id: PEUL_ExchangeCornerEdgeCells.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ExchangeCornerEdgeCells( regions,iReg )

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
  TYPE(t_region), POINTER    :: regions(:)
  INTEGER,        INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: iedge, icorner, ijk, i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: icell, iRegSrc, iLev, iCOff, ijCOff, ijkD
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend

  REAL(RFREAL), POINTER :: cv(:,:), cvSrc(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_level),  POINTER :: level

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_ExchangeCornerEdgeCells.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PEUL_ExchangeCornerEdgeCells',&
  'PEUL_ExchangeCornerEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = regions(iReg)%currLevel

  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

  level => regions(iReg)%levels(iLev)
  cv    => level%peul%cv

! edge cells

  DO iedge=1,12
    IF (level%edgeCells(iedge)%interact) THEN
      CALL RFLO_GetEdgeCellsIndices( regions(iReg),iLev,iedge, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )
      ijk = 0
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijk     = ijk + 1
            ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
            icell   = level%edgeCells(iedge)%cells(ijk)%srcCell
            iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
            IF (iRegSrc > 0) THEN
              IF (regions(iRegSrc)%procid == global%myProcid) THEN
                cvSrc => regions(iRegSrc)%levels(iLev)%peul%cv

                IF (level%edgeCells(iedge)%cells(ijk)%rotate) THEN
                  ! rotational periodicity
                ELSE
                  cv(:,ijkD) = cvSrc(:,icell)
                ENDIF

              ENDIF
            ENDIF
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDIF        ! interact
  ENDDO          ! iedge

! corner cells

  DO icorner=1,8
    IF (level%cornerCells(icorner)%interact) THEN
      CALL RFLO_GetCornerCellsIndices( regions(iReg),iLev,icorner, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )
      ijk = 0
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijk     = ijk + 1
            ijkD    = IndIJK(i,j,k,iCOff,ijCOff)
            icell   = level%cornerCells(icorner)%cells(ijk)%srcCell
            iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
            IF (iRegSrc > 0) THEN
              IF (regions(iRegSrc)%procid == global%myProcid) THEN
                cvSrc => regions(iRegSrc)%levels(iLev)%peul%cv

                IF (level%cornerCells(icorner)%cells(ijk)%rotate) THEN
                  ! rotational periodicity
                ELSE
                  cv(:,ijkD) = cvSrc(:,icell)
                ENDIF

              ENDIF
            ENDIF
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDIF        ! interact
  ENDDO          ! icorner

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ExchangeCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ExchangeCornerEdgeCells.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:34  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
!******************************************************************************







