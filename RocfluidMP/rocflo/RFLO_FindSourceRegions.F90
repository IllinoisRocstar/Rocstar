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
! Purpose: find source regions for edge and corner cells.
!
! Description: none.
!
! Input: regions = dimensions and topology of all regions.
!
! Output: regions = pointers to edge and corner cells for each region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_FindSourceRegions.F90,v 1.5 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_FindSourceRegions( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices, RFLO_FindSourceCell, &
                            RFLO_FindSourceCellInvert
  USE RFLO_ModDegenerateCornEdge, ONLY: RFLO_FindDegeneratCell
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev, iedge, icorner, i, j, k, ijk, nt

! ... local variables
  LOGICAL :: found, rotate

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, errorFlag
  INTEGER :: ic, jc, kc, icell, iRegTemp, iRegSrc

  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER  :: level

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_FindSourceRegions',&
  'RFLO_FindSourceRegions.F90' )

! search for source regions ===================================================

  global%degenrtEc = .FALSE.

  DO iReg=1,global%nRegions
    DO iLev=1,regions(iReg)%nGridLevels

      level => regions(iReg)%levels(iLev)

      CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                               jpcbeg,jpcend,kpcbeg,kpcend )

! --- edge cells --------------------------------------------------------------

      DO iedge=1,12
        CALL RFLO_GetEdgeCellsIndices( regions(iReg),iLev,iedge, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )

! ----- allocate memory for cells, reset interaction and degeneration

        ijk = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)
        ALLOCATE( level%edgeCells(iedge)%cells(ijk),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

        level%edgeCells(iedge)%interact = .false.
        level%edgeCells(iedge)%degenrt  = DEGENERAT_NONE
        level%edgeCells(iedge)%interType= EDGE_INTERACT_FULL

! ----- loop over all cells of iedge

        ijk = 0
        DO k=kbeg,kend
          DO j=jbeg,jend
            DO i=ibeg,iend
              ijk      = ijk + 1
              ic       = i
              jc       = j
              kc       = k
              iRegTemp = iReg
              DO nt=1,2
                CALL RFLO_FindSourceCell( regions,iRegTemp,iLev,ic,jc,kc, &
                                          icell,found,rotate,iRegSrc )
                IF (found) EXIT
                iRegTemp = iRegSrc
              ENDDO
              IF (found) THEN
                level%edgeCells(iedge)%interact             = .true.
                level%edgeCells(iedge)%cells(ijk)%srcCell   = icell
                level%edgeCells(iedge)%cells(ijk)%srcRegion = iRegSrc
                level%edgeCells(iedge)%cells(ijk)%rotate    = rotate
              ELSE
                level%edgeCells(iedge)%interType            = EDGE_INTERACT_PART
                level%edgeCells(iedge)%cells(ijk)%srcCell   = -999999
                level%edgeCells(iedge)%cells(ijk)%srcRegion = -999999
                level%edgeCells(iedge)%cells(ijk)%rotate    = .false.
              ENDIF  ! primary found

! ----------- this part is added to detect degenerative edges

              IF (found .AND. &
                  level%edgeCells(iedge)%degenrt == DEGENERAT_NONE) THEN
                CALL RFLO_FindDegeneratCell( regions,0,iedge,iReg,iLev, &
                                             iRegSrc,icell,ic,jc,kc )

                IF (level%edgeCells(iedge)%degenrt == DEGENERAT_NONE) THEN
                  ic       = i
                  jc       = j
                  kc       = k
                  iRegTemp = iReg
                  DO nt=1,2
                    CALL RFLO_FindSourceCellInvert( regions,iRegTemp,iLev, &
                                          ic,jc,kc,icell,found,rotate,iRegSrc )
                    IF (found) EXIT
                    iRegTemp = iRegSrc
                  ENDDO
                  IF (found .AND. (iRegSrc/= &
                      level%edgeCells(iedge)%cells(ijk)%srcRegion .OR. &
                      icell/= &
                      level%edgeCells(iedge)%cells(ijk)%srcCell)) THEN
                    level%edgeCells(iedge)%degenrt = DEGENERAT_DETACH
                  ENDIF  ! assign negative degenrt
                ENDIF  ! search for negative degeneration

              ENDIF  ! search for positive degeneration

! ----------- end degeneration kernel 

            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k
        IF (.NOT. level%edgeCells(iedge)%interact) THEN
          DEALLOCATE( level%edgeCells(iedge)%cells )
        ENDIF
      ENDDO    ! iedge

! --- corner cells ------------------------------------------------------------

      DO icorner=1,8
        CALL RFLO_GetCornerCellsIndices( regions(iReg),iLev,icorner, &
                                         ibeg,iend,jbeg,jend,kbeg,kend )

! ----- allocate memory for cells, reset interaction

        ijk = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)
        ALLOCATE( level%cornerCells(icorner)%cells(ijk),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

        level%cornerCells(icorner)%interact = .false.
        level%cornerCells(icorner)%degenrt  = DEGENERAT_NONE

! ----- loop over all cells of icorner

        ijk = 0
        DO k=kbeg,kend
          DO j=jbeg,jend
            DO i=ibeg,iend
              ijk      = ijk + 1
              ic       = i
              jc       = j
              kc       = k
              iRegTemp = iReg
              DO nt=1,3
                CALL RFLO_FindSourceCell( regions,iRegTemp,iLev,ic,jc,kc, &
                                          icell,found,rotate,iRegSrc )
                IF (found) EXIT
                iRegTemp = iRegSrc
              ENDDO
              IF (found) THEN
                level%cornerCells(icorner)%interact             = .true.
                level%cornerCells(icorner)%cells(ijk)%srcCell   = icell
                level%cornerCells(icorner)%cells(ijk)%srcRegion = iRegSrc
                level%cornerCells(icorner)%cells(ijk)%rotate    = rotate
              ELSE
                level%cornerCells(icorner)%cells(ijk)%srcCell   = -999999
                level%cornerCells(icorner)%cells(ijk)%srcRegion = -999999
                level%cornerCells(icorner)%cells(ijk)%rotate    = .false.
              ENDIF  ! primary found

! ----------- this part is added to detect degenerative corners

              IF (found .AND. &
                  level%cornerCells(icorner)%degenrt == DEGENERAT_NONE) THEN
                CALL RFLO_FindDegeneratCell( regions,1,icorner,iReg,iLev, &
                                             iRegSrc,icell,ic,jc,kc )

                IF (level%cornerCells(icorner)%degenrt == DEGENERAT_NONE) THEN
                  ic       = i
                  jc       = j
                  kc       = k
                  iRegTemp = iReg
                  DO nt=1,3
                    CALL RFLO_FindSourceCellInvert( regions,iRegTemp,iLev, &
                                          ic,jc,kc,icell,found,rotate,iRegSrc )
                    IF (found) EXIT
                    iRegTemp = iRegSrc
                  ENDDO
                  IF (found .AND. (iRegSrc/= &
                      level%cornerCells(icorner)%cells(ijk)%srcRegion .OR. &
                      icell/= &
                      level%cornerCells(icorner)%cells(ijk)%srcCell)) THEN
                    level%cornerCells(icorner)%degenrt = DEGENERAT_DETACH
                  ENDIF  ! assign negative degenrt
                ENDIF  ! search for negative degeneration

              ENDIF  ! search for positive degeneration

! ----------- end degeneration kernel 

            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k
        IF (.NOT. level%cornerCells(icorner)%interact) THEN
          DEALLOCATE( level%cornerCells(icorner)%cells )
        ENDIF
      ENDDO    ! icorner

    ENDDO      ! iLev
  ENDDO        ! iReg

! # cells sent to/received from another processor =============================

  IF (global%nProcAlloc > 1) THEN    ! only if multiple processors

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid == global%myProcid) THEN
        DO iLev=1,regions(iReg)%nGridLevels

          level => regions(iReg)%levels(iLev)

          ALLOCATE( level%sendEcCells(global%nRegions),stat=errorFlag )
          ALLOCATE( level%recvEcCells(global%nRegions),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          level%sendEcCells(:)%nCells = 0
          level%recvEcCells(:)%nCells = 0

          DO iedge=1,12
            IF (level%edgeCells(iedge)%interact) THEN
              DO ijk=1,UBOUND(level%edgeCells(iedge)%cells,1)
                iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
                IF (iRegSrc > 0) THEN
                  IF (regions(iRegSrc)%procid /= global%myProcid) &
                    level%recvEcCells(iRegSrc)%nCells = &
                      level%recvEcCells(iRegSrc)%nCells + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          DO icorner=1,8
            IF (level%cornerCells(icorner)%interact) THEN
              DO ijk=1,UBOUND(level%cornerCells(icorner)%cells,1)
                iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
                IF (iRegSrc > 0) THEN
                  IF (regions(iRegSrc)%procid /= global%myProcid) &
                    level%recvEcCells(iRegSrc)%nCells = &
                      level%recvEcCells(iRegSrc)%nCells + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO

        ENDDO    ! iLev
      ENDIF      ! my processor
    ENDDO        ! iReg

! - cells to send

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid /= global%myProcid) THEN
        DO iLev=1,regions(iReg)%nGridLevels

          level => regions(iReg)%levels(iLev)

          DO iedge=1,12
            IF (level%edgeCells(iedge)%interact) THEN
              DO ijk=1,UBOUND(level%edgeCells(iedge)%cells,1)
                iRegSrc = level%edgeCells(iedge)%cells(ijk)%srcRegion
                IF (iRegSrc > 0) THEN
                  IF (regions(iRegSrc)%procid == global%myProcid) &
                    regions(iRegSrc)%levels(iLev)%sendEcCells(iReg)%nCells = &
                      regions(iRegSrc)%levels(iLev)%sendEcCells(iReg)%nCells + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          DO icorner=1,8
            IF (level%cornerCells(icorner)%interact) THEN
              DO ijk=1,UBOUND(level%cornerCells(icorner)%cells,1)
                iRegSrc = level%cornerCells(icorner)%cells(ijk)%srcRegion
                IF (iRegSrc > 0) THEN
                  IF (regions(iRegSrc)%procid == global%myProcid) &
                    regions(iRegSrc)%levels(iLev)%sendEcCells(iReg)%nCells = &
                      regions(iRegSrc)%levels(iLev)%sendEcCells(iReg)%nCells + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO

        ENDDO    ! iLev
      ENDIF      ! not my processor
    ENDDO        ! iReg

  ENDIF

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_FindSourceRegions

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_FindSourceRegions.F90,v $
! Revision 1.5  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/11/11 07:19:09  wasistho
! moved RFLO_FindDegeneratCell to modflo module
!
! Revision 1.2  2005/06/29 22:52:43  wasistho
! find value of interType
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.13  2004/08/21 00:32:22  wasistho
! insert kernel to search degenerated edge/corners
!
! Revision 1.12  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.7  2003/02/14 22:32:37  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.6  2003/02/06 23:56:11  jblazek
! Corrected bug for iRegSrc<0.
!
! Revision 1.5  2003/02/03 19:20:47  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.4  2002/10/25 18:36:47  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.3  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.2  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.6  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.5  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.4  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.3  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.2  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************







