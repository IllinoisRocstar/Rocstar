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
! Output: new values of conservative variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeCornerEdgeCells.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeCornerEdgeCells( regions,iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices, MixtureProperties
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iedge, icorner, ijk, i, j, k

! ... local variables
  INTEGER :: icell, iRegSrc, iLev, gasModel, iCOff, ijCOff, ijkD
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend

  REAL(RFREAL), POINTER :: cv(:,:), cvSrc(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER  :: level

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_ExchangeCornerEdgeCells',&
  'RFLO_ExchangeCornerEdgeCells.F90' )

! get dimensions and pointers

  iLev = regions(iReg)%currLevel

  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

  gasModel =  regions(iReg)%mixtInput%gasModel
  level     => regions(iReg)%levels(iLev)
  cv        => level%mixt%cv

! edge cells

  DO iedge=1,12
    IF (level%edgeCells(iedge)%interact .AND. &
        level%edgeCells(iedge)%degenrt==DEGENERAT_NONE) THEN
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
                cvSrc => regions(iRegSrc)%levels(iLev)%mixt%cv

                IF (level%edgeCells(iedge)%cells(ijk)%rotate) THEN
                  ! rotational periodicity
                ELSE
                  cv(CV_MIXT_DENS,ijkD) = cvSrc(CV_MIXT_DENS,icell)
                  cv(CV_MIXT_XMOM,ijkD) = cvSrc(CV_MIXT_XMOM,icell)
                  cv(CV_MIXT_YMOM,ijkD) = cvSrc(CV_MIXT_YMOM,icell)
                  cv(CV_MIXT_ZMOM,ijkD) = cvSrc(CV_MIXT_ZMOM,icell)
                  cv(CV_MIXT_ENER,ijkD) = cvSrc(CV_MIXT_ENER,icell)
                ENDIF

                IF (gasModel == GAS_MODEL_TCPERF) THEN
                  CALL MixtureProperties( regions(iReg),ijkD,ijkD,.false. )
                ELSE
                  CALL MixtureProperties( regions(iReg),ijkD,ijkD,.true.  )
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
    IF (level%cornerCells(icorner)%interact .AND. &
        level%cornerCells(icorner)%degenrt==DEGENERAT_NONE) THEN
      CALL RFLO_GetCornerCellsIndices( regions(iReg),iLev,icorner, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )
      ijk = 0
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijk     =  ijk + 1
            ijkD    =  IndIJK(i,j,k,iCOff,ijCOff)
            icell   =  level%cornerCells(icorner)%cells(ijk)%srcCell
            iRegSrc =  level%cornerCells(icorner)%cells(ijk)%srcRegion
            IF (iRegSrc > 0) THEN
              IF (regions(iRegSrc)%procid == global%myProcid) THEN
                cvSrc => regions(iRegSrc)%levels(iLev)%mixt%cv

                IF (level%cornerCells(icorner)%cells(ijk)%rotate) THEN
                  ! rotational periodicity
                ELSE
                  cv(CV_MIXT_DENS,ijkD) = cvSrc(CV_MIXT_DENS,icell)
                  cv(CV_MIXT_XMOM,ijkD) = cvSrc(CV_MIXT_XMOM,icell)
                  cv(CV_MIXT_YMOM,ijkD) = cvSrc(CV_MIXT_YMOM,icell)
                  cv(CV_MIXT_ZMOM,ijkD) = cvSrc(CV_MIXT_ZMOM,icell)
                  cv(CV_MIXT_ENER,ijkD) = cvSrc(CV_MIXT_ENER,icell)
                ENDIF

                IF (gasModel == GAS_MODEL_TCPERF) THEN
                  CALL MixtureProperties( regions(iReg),ijkD,ijkD,.false. )
                ELSE
                  CALL MixtureProperties( regions(iReg),ijkD,ijkD,.true.  )
                ENDIF
              ENDIF
            ENDIF
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDIF        ! interact
  ENDDO          ! icorner

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExchangeCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeCornerEdgeCells.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.13  2004/09/02 01:35:02  wasistho
! replaced 0 with DEGENERAT_NONE
!
! Revision 1.12  2004/09/02 01:22:52  wasistho
! do not exchange and send/recv for degenerated e/c
!
! Revision 1.11  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.6  2003/02/14 22:32:37  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.5  2003/02/03 19:20:47  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
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







