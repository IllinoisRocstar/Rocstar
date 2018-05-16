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
! Purpose: copy a matrix from the boundary cells to the dummy cells.
!
! Description: file contains the following subroutines:
!
!  - RFLO_copyMatrixPatches = copy for patches
!  - RFLO_copyMatrixEdges   = copy for edge cells
!  - RFLO_copyMatrixCorners = copy for corner cells
!
! Input: iLev   = grid level
!        region = dimensions of current region
!        mat    = data to be copied
!
! Output: mat = data copied to dummy cells.
!
! Notes: it is assumed that the first dimension of the matrix is to be
!        copied completely. No averaging is done. Thus, these routines are
!        primarily intended to copy geometrical quantities rather than
!        flow variables. All edge and corner cells are copied, regardless
!        of their connection to other regions. Routines are used now to copy
!        the spectral radii.
!
!******************************************************************************
!
! $Id: RFLO_CopyMatrixDummy.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_copyMatrixPatches( iLev,region,mat )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetPatchIndices, &
                            RFLO_GetPatchDirection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  INTEGER :: iLev

  REAL(RFREAL), POINTER :: mat(:,:)

! ... loop variables
  INTEGER :: iPatch, i, j, k, idum

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iCOff, ijCOff, ijkCell, ijkDum, bcType

  TYPE(t_patch), POINTER :: patch

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_copyMatrixPatches',&
  'RFLO_CopyMatrixDummy.F90' )

! loop over patches (all but inter-region and periodic boundaries)

  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  DO iPatch=1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    bcType = patch%bcType
    IF ((bcType<BC_REGIONCONF .OR. bcType>BC_REGIONCONF+BC_RANGE) .AND.&
        (bcType<BC_REGIONINT  .OR. bcType>BC_REGIONINT +BC_RANGE) .AND.&
        (bcType<BC_REGNONCONF .OR. bcType>BC_REGNONCONF+BC_RANGE) .AND.&
        (bcType<BC_TRA_PERI   .OR. bcType>BC_TRA_PERI  +BC_RANGE) .AND.&
        (bcType<BC_ROT_PERI   .OR. bcType>BC_ROT_PERI  +BC_RANGE)) THEN
      CALL RFLO_GetPatchIndices( region,patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
      DO idum=1,region%nDumCells
        DO k=kbeg,kend
          DO j=jbeg,jend
            DO i=ibeg,iend
              ijkCell = IndIJK(i,j,k,iCOff,ijCOff)
              ijkDum  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
              mat(:,ijkDum) = mat(:,ijkCell)
            ENDDO  ! i
          ENDDO    ! j
        ENDDO      ! k
      ENDDO        ! idum
    ENDIF          ! bcType
  ENDDO            ! iPatch

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_copyMatrixPatches

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_copyMatrixEdges( iLev,region,mat )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetEdgeCellsIndices
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  INTEGER :: iLev

  REAL(RFREAL), POINTER :: mat(:,:)

! ... loop variables
  INTEGER :: i, j, k, iedge

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: iCOff, ijCOff, ijkCell, ijkDum

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_copyMatrixEdges',&
  'RFLO_CopyMatrixDummy.F90' )

! loop over edge cells without source region
! or with multiple source regions

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  DO iedge=1,12
    CALL RFLO_GetEdgeCellsIndices( region,iLev,iedge, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          IF (iedge ==  1) ijkCell = IndIJK(ipcbeg,jpcbeg,k,iCOff,ijCOff)
          IF (iedge ==  2) ijkCell = IndIJK(ipcbeg,j,kpcend,iCOff,ijCOff)
          IF (iedge ==  3) ijkCell = IndIJK(ipcbeg,jpcend,k,iCOff,ijCOff)
          IF (iedge ==  4) ijkCell = IndIJK(ipcbeg,j,kpcbeg,iCOff,ijCOff)
          IF (iedge ==  5) ijkCell = IndIJK(ipcend,jpcbeg,k,iCOff,ijCOff)
          IF (iedge ==  6) ijkCell = IndIJK(ipcend,j,kpcend,iCOff,ijCOff)
          IF (iedge ==  7) ijkCell = IndIJK(ipcend,jpcend,k,iCOff,ijCOff)
          IF (iedge ==  8) ijkCell = IndIJK(ipcend,j,kpcbeg,iCOff,ijCOff)
          IF (iedge ==  9) ijkCell = IndIJK(i,jpcbeg,kpcbeg,iCOff,ijCOff)
          IF (iedge == 10) ijkCell = IndIJK(i,jpcbeg,kpcend,iCOff,ijCOff)
          IF (iedge == 11) ijkCell = IndIJK(i,jpcend,kpcbeg,iCOff,ijCOff)
          IF (iedge == 12) ijkCell = IndIJK(i,jpcend,kpcend,iCOff,ijCOff)
          ijkDum        = IndIJK(i,j,k,iCOff,ijCOff)
          mat(:,ijkDum) = mat(:,ijkCell)
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO          ! iedge

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_copyMatrixEdges

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_copyMatrixCorners( iLev,region,mat )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, VolumeHexa, &
        RFLO_GetPatchIndices, RFLO_GetPatchDirection, RFLO_GetEdgeCellsIndices,&
        RFLO_GetCornerCellsIndices
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  INTEGER :: iLev

  REAL(RFREAL), POINTER :: mat(:,:)

! ... loop variables
  INTEGER :: i, j, k, icorn

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: iCOff, ijCOff, ijkCell, ijkDum

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_copyMatrixCorners',&
  'RFLO_CopyMatrixDummy.F90' )

! loop over corner cells without source region
! or with multiple source regions

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  DO icorn=1,8
    CALL RFLO_GetCornerCellsIndices( region,iLev,icorn, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          IF (icorn == 1) ijkCell = IndIJK(ipcbeg,jpcbeg,kpcbeg,iCOff,ijCOff)
          IF (icorn == 2) ijkCell = IndIJK(ipcbeg,jpcbeg,kpcend,iCOff,ijCOff)
          IF (icorn == 3) ijkCell = IndIJK(ipcbeg,jpcend,kpcend,iCOff,ijCOff)
          IF (icorn == 4) ijkCell = IndIJK(ipcbeg,jpcend,kpcbeg,iCOff,ijCOff)
          IF (icorn == 5) ijkCell = IndIJK(ipcend,jpcbeg,kpcbeg,iCOff,ijCOff)
          IF (icorn == 6) ijkCell = IndIJK(ipcend,jpcbeg,kpcend,iCOff,ijCOff)
          IF (icorn == 7) ijkCell = IndIJK(ipcend,jpcend,kpcend,iCOff,ijCOff)
          IF (icorn == 8) ijkCell = IndIJK(ipcend,jpcend,kpcbeg,iCOff,ijCOff)
          ijkDum        = IndIJK(i,j,k,iCOff,ijCOff)
          mat(:,ijkDum) = mat(:,ijkCell)
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO          ! icorn

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_copyMatrixCorners

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CopyMatrixDummy.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.16  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.12  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.11  2003/02/03 19:20:47  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.10  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.9  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.8  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.7  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.6  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.5  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.3  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.2  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.1  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
!******************************************************************************









