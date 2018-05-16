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
! ******************************************************************************
!
! Purpose: Suite for Laplacian-smoothing routines.
!
! Description: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModLaplaceSmoothing.F90,v 1.13 2009/08/27 14:04:50 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModLaplaceSmoothing

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_LaplaceGridSmoo, &
            RFLO_LaplaceGridSolve, &
            RFLO_LaplaceGridPatch, &
            RFLO_LaplaceGridJump

! private :
!           RFLO_LaplaceGridOrtho
!           RFLO_ProjectQuadCorner
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModLaplaceSmoothing.F90,v $ $Revision: 1.13 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: smooth the distribution of grid points by solving simplified
!          Laplace equation in physical space.
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions%levels%grid%xyz = new grid coordinates
!         resid = convergence of the Jacobi iteration.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridSmoo( regions,resid )

  USE ModInterfaces, ONLY : RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_ExchangeDnodeCopy, &
        RFLO_ExchangeDnodeSend, RFLO_ExchangeDnodeRecv, RFLO_ClearSendRequests,&
        RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset, &
        RFLO_CalcCellCentroids, RFLO_CalcFaceCentroids, RFLO_ChangeInteriorGrid
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  REAL(RFREAL) :: resid

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, ijk, i, j, k, iter

! ... local variables
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff
  INTEGER :: bcType, iRegSrc, iPatchSrc

  REAL(RFREAL) :: dx, dy, dz
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:), xyzTemp(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER   :: grid, gridOld, gridSrc
  TYPE(t_patch), POINTER  :: patch, patchSrc

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_LaplaceGridSmoo',&
       'RFLO_ModLaplaceSmoothing.F90' )

! smooth grid region-wise -----------------------------------------------------

  resid = 0._RFREAL

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      IF (global%moveGridScheme==MOVEGRID_FOMS) THEN

! ----- enforce grid orthogonality
        CALL RFLO_LaplaceGridOrtho( regions(iReg) )

      ELSE
! ----- compute movements in the interior and along the boundaries
        CALL RFLO_LaplaceGridSolve( regions(iReg) )

      ENDIF

! --- zero out movements along certain boundaries

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_INFLOW     .AND. bcType<=BC_INFLOW    +BC_RANGE) .OR. &
            (bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW   +BC_RANGE) .OR. &
            (bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL  +BC_RANGE) .OR. &
            (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
            (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE) .OR. &
            (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_SYMMETRY_FREE .AND. bcType<=BC_SYMMETRY +BC_RANGE)) THEN
          CALL RFLO_LaplaceGridPatch( regions(iReg),patch )
        ENDIF  ! bcType
      ENDDO    ! iPatch

      CALL RFLO_GetDimensPhysNodes( regions(iReg),1,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
      CALL RFLO_GetNodeOffset( regions(iReg),1,iNOff,ijNOff )

      xyz    => regions(iReg)%levels(1)%grid%xyz
      xyzOld => regions(iReg)%levels(1)%grid%xyzOld

      DO k=kpnbeg,kpnend
        DO j=jpnbeg,jpnend
          DO i=ipnbeg,ipnend
            ijk   = IndIJK(i,j,k,iNOff,ijNOff)
            dx    = xyz(XCOORD,ijk) - xyzOld(XCOORD,ijk)
            dy    = xyz(YCOORD,ijk) - xyzOld(YCOORD,ijk)
            dz    = xyz(ZCOORD,ijk) - xyzOld(ZCOORD,ijk)
            resid = resid + dx*dx + dy*dy +dz*dz
          ENDDO
        ENDDO
      ENDDO

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! fix interfaces between regions ----------------------------------------------
! copy / send deformations

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)
          gridSrc   => regions(iRegSrc)%levels(1)%grid

          IF (regions(iRegSrc)%procid == global%myProcid) THEN
            CALL RFLO_ExchangeDnodeCopy( regions(iReg),regions(iRegSrc), &
                                         patch,patchSrc,.true., &
                                         grid%xyz,gridSrc%xyz )
          ELSE
            CALL RFLO_ExchangeDnodeSend( regions(iReg),regions(iRegSrc), &
                                         patch,grid%xyz )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! receive deformations

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)
          gridSrc   => regions(iRegSrc)%levels(1)%grid

          IF (regions(iRegSrc)%procid /= global%myProcid) THEN
            CALL RFLO_ExchangeDnodeRecv( regions(iReg),regions(iRegSrc), &
                                         patch,patchSrc,.true.,grid%xyz )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! clear send requests

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving
      CALL RFLO_ClearSendRequests( regions,iReg,.true. )
    ENDIF
  ENDDO

! update grid, dummy, corner and edge cells -----------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! --- change xyz from deformations to coordinates

      xyz    => regions(iReg)%levels(1)%grid%xyz
      xyzOld => regions(iReg)%levels(1)%gridOld%xyz

      DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
        xyz(XCOORD,ijk) = xyz(XCOORD,ijk) + xyzOld(XCOORD,ijk)
        xyz(YCOORD,ijk) = xyz(YCOORD,ijk) + xyzOld(YCOORD,ijk)
        xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) + xyzOld(ZCOORD,ijk)
      ENDDO

! --- update coarse grids and dummy cells

      CALL RFLO_GenerateCoarseGrids( regions(iReg) )   ! coarsen finest grid
      CALL RFLO_CopyGeometryDummy( regions(iReg) )     ! copy to dummy nodes
      CALL RFLO_ExtrapolateGeometry( regions(iReg) )   ! extrapolate
    ENDIF     ! region on this processor and active, grid moving
  ENDDO       ! iReg

  CALL RFLO_ExchangeGeometry( regions )                ! exchange geometry

! recompute cell and face centroids for FOMS grid scheme ----------------------

!  IF (global%moveGridScheme/=MOVEGRID_FOMS) GOTO 777
  GOTO 777

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE .AND. &             ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN           ! and moving
      CALL RFLO_CalcCellCentroids( regions(iReg) )       ! cell centroids
      CALL RFLO_CalcFaceCentroids( regions(iReg) )       ! face centroids
      CALL RFLO_CalcFaceVectors( regions(iReg) )
    ENDIF   ! region on this processor and active, grid moving
  ENDDO     ! iReg

  DO iter=1,1

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! --- enforce grid orthogonality

      CALL RFLO_LaplaceGridOrtho( regions(iReg) )

! --- zero out movements and recored jumps along certain boundaries

!      regions(iReg)%levels(1)%grid%xyzTemp = 0._RFREAL

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_INFLOW  .AND. bcType<=BC_INFLOW    +BC_RANGE).OR. &
            (bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW   +BC_RANGE).OR. &
            (bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL  +BC_RANGE).OR. &
            (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE).OR. &
            (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE).OR. &
            (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE).OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE).OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE).OR. &
            (bcType>=BC_SYMMETRY_FREE .AND. bcType<=BC_SYMMETRY +BC_RANGE)) THEN
          CALL RFLO_LaplaceGridPatch( regions(iReg),patch )
        ENDIF
!        CALL RFLO_LaplaceGridJump( regions(iReg),patch )
      ENDDO    ! iPatch

! --- change the interior movements

!      grid    => regions(iReg)%levels(1)%grid
!      gridOld => regions(iReg)%levels(1)%gridOld

!      CALL RFLO_ChangeInteriorGrid( regions(iReg),grid%boundMoved, &
!                                    grid%edgeMoved,grid%arcLen12, &
!                                    grid%arcLen34,grid%arcLen56, &
!                                    gridOld%xyzOld,grid%xyzTemp )

!      CALL RFLO_GetDimensPhysNodes( regions(iReg),1,ipnbeg,ipnend, &
!                                    jpnbeg,jpnend,kpnbeg,kpnend )
!      CALL RFLO_GetNodeOffset( regions(iReg),1,iNOff,ijNOff )
!IF (iReg==13) write(*,*)'ORTHO1',grid%xyz(:,IndIJK(10,10,10,iNOff,ijNOff))+ gridOld%xyz(:,IndIJK(10,10,10,iNOff,ijNOff)), grid%xyzOrth(:,IndIJK(10,10,10,iNOff,ijNOff))+ gridOld%xyz(:,IndIJK(10,10,10,iNOff,ijNOff))

!      DO k=kpnbeg,kpnend
!        DO j=jpnbeg,jpnend
!          DO i=ipnbeg,ipnend
!            ijk = IndIJK(i,j,k,iNOff,ijNOff)
!            grid%xyz(:,ijk) = grid%xyzTemp(:,ijk) - &
!                              gridOld%xyz(:,ijk) + grid%xyzOrth(:,ijk)
!          ENDDO
!        ENDDO
!      ENDDO
!IF (iReg==13) write(*,*)'ORTHO2',grid%xyz(:,IndIJK(10,10,10,iNOff,ijNOff)),grid%xyzOrth(:,IndIJK(10,10,10,iNOff,ijNOff))

! --- compute residuals

      CALL RFLO_GetDimensPhysNodes( regions(iReg),1,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
      CALL RFLO_GetNodeOffset( regions(iReg),1,iNOff,ijNOff )

      xyz    => regions(iReg)%levels(1)%grid%xyz
      xyzOld => regions(iReg)%levels(1)%grid%xyzOld

      DO k=kpnbeg,kpnend
        DO j=jpnbeg,jpnend
          DO i=ipnbeg,ipnend
            ijk   = IndIJK(i,j,k,iNOff,ijNOff)
            dx    = xyz(XCOORD,ijk) - xyzOld(XCOORD,ijk)
            dy    = xyz(YCOORD,ijk) - xyzOld(YCOORD,ijk)
            dz    = xyz(ZCOORD,ijk) - xyzOld(ZCOORD,ijk)
            resid = resid + dx*dx + dy*dy +dz*dz
          ENDDO
        ENDDO
      ENDDO

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! fix interfaces between regions ----------------------------------------------
! copy / send deformations

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)
          gridSrc   => regions(iRegSrc)%levels(1)%grid

          IF (regions(iRegSrc)%procid == global%myProcid) THEN
            CALL RFLO_ExchangeDnodeCopy( regions(iReg),regions(iRegSrc), &
                                         patch,patchSrc,.true., &
                                         grid%xyz,gridSrc%xyz )
          ELSE
            CALL RFLO_ExchangeDnodeSend( regions(iReg),regions(iRegSrc), &
                                         patch,grid%xyz )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! receive deformations

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)
          gridSrc   => regions(iRegSrc)%levels(1)%grid

          IF (regions(iRegSrc)%procid /= global%myProcid) THEN
            CALL RFLO_ExchangeDnodeRecv( regions(iReg),regions(iRegSrc), &
                                         patch,patchSrc,.true.,grid%xyz )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! clear send requests

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving
      CALL RFLO_ClearSendRequests( regions,iReg,.true. )
    ENDIF
  ENDDO

! update grid, dummy, corner and edge cells -----------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! --- change xyz from deformations to coordinates

      xyz    => regions(iReg)%levels(1)%grid%xyz
      xyzOld => regions(iReg)%levels(1)%gridOld%xyz

      DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
        xyz(XCOORD,ijk) = xyz(XCOORD,ijk) + xyzOld(XCOORD,ijk)
        xyz(YCOORD,ijk) = xyz(YCOORD,ijk) + xyzOld(YCOORD,ijk)
        xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) + xyzOld(ZCOORD,ijk)
      ENDDO

! --- update coarse grids and dummy cells

      CALL RFLO_GenerateCoarseGrids( regions(iReg) )   ! coarsen finest grid
      CALL RFLO_CopyGeometryDummy( regions(iReg) )     ! copy to dummy nodes
      CALL RFLO_ExtrapolateGeometry( regions(iReg) )   ! extrapolate
    ENDIF     ! region on this processor and active, grid moving
  ENDDO       ! iReg

  CALL RFLO_ExchangeGeometry( regions )                ! exchange geometry

  ENDDO  ! iter

! finalize --------------------------------------------------------------------

777 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_LaplaceGridSmoo


!******************************************************************************
!
! Purpose: zero out movements obtained by LaplaceGridSolve along a boundary.
!
! Description: none.
!
! Input: region = data of current region, grid movements
!        patch  = current patch.
!
! Output: region%levels%grid%xyz = new grid movements.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridPatch( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, ijkNB

  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_LaplaceGridPatch',&
       'RFLO_ModLaplaceSmoothing.F90' )

! get dimensions and pointers

  iLev = 1

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz    => region%levels(iLev)%grid%xyz
  xyzOld => region%levels(iLev)%grid%xyzOld

! new = old

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkNB  = IndIJK(i,j,k,iNOff,ijNOff)
        xyz(XCOORD,ijkNB) = xyzOld(XCOORD,ijkNB)
        xyz(YCOORD,ijkNB) = xyzOld(YCOORD,ijkNB)
        xyz(ZCOORD,ijkNB) = xyzOld(ZCOORD,ijkNB)
      ENDDO
    ENDDO
  ENDDO

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_LaplaceGridPatch

!******************************************************************************
!
! Purpose: obtain jump of surface grid motion due to enforcing back boundary
!          surface motion.
!
! Description: none.
!
! Input: region = data of current region, grid movements
!        patch  = current patch.
!
! Output: region%levels%grid%xyz = new grid movements.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridJump( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, ijkNB

  REAL(RFREAL), POINTER :: xyz(:,:), xyzOrth(:,:), xyzTemp(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_LaplaceGridJump',&
       'RFLO_ModLaplaceSmoothing.F90' )

! get dimensions and pointers

  iLev = 1

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz     => region%levels(iLev)%grid%xyz
  xyzOrth => region%levels(iLev)%grid%xyzOrth
  xyzTemp => region%levels(iLev)%grid%xyzTemp

! jump = new-(orth-orig)

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkNB  = IndIJK(i,j,k,iNOff,ijNOff)
        xyzTemp(XCOORD,ijkNB) = xyz(XCOORD,ijkNB) - xyzOrth(XCOORD,ijkNB)
        xyzTemp(YCOORD,ijkNB) = xyz(YCOORD,ijkNB) - xyzOrth(YCOORD,ijkNB)
        xyzTemp(ZCOORD,ijkNB) = xyz(ZCOORD,ijkNB) - xyzOrth(ZCOORD,ijkNB)
      ENDDO
    ENDDO
  ENDDO

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_LaplaceGridJump


!******************************************************************************
!
! Purpose: conduct one Jacobi iteration to obtain new grid movements
!          in the flow domain (boundaries included).
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates.
!
! Output: region%levels%grid%xyz = grid movements.
!
! Notes: on entry, xyz holds node coordinates from a previous smoothing
!        step. On exit however, xyz contains only the grid motion. 
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridSolve( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetDimensPhysNodes, &
                            RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, iLev, iNOff, ijNOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, ibn, ien
  INTEGER :: ijk, ip1, ip2, im1, im2, jp1, jp2, jm1, jm2, kp1, kp2, km1, km2
  INTEGER :: moveBlock, method, in, ico, ic(8)
  INTEGER, POINTER :: idgen(:)

  LOGICAL, POINTER :: bndMoved(:)

  REAL(RFREAL) :: rxi, ryi, rzi, rxj, ryj, rzj, rxk, ryk, rzk, wi, wj, wk, &
                  si, sj, sk, sim, sjm, skm, sip, sjp, skp, d, p, eps
  REAL(RFREAL) :: denom(8), dist(8,8)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:), xyzOrig(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_LaplaceGridSolve',&
       'RFLO_ModLaplaceSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  xyz      => region%levels(iLev)%grid%xyz
  xyzOld   => region%levels(iLev)%grid%xyzOld
  xyzOrig  => region%levels(iLev)%gridOld%xyz
  bndMoved => region%levels(iLev)%grid%boundMoved
  idgen    => region%levels(iLev)%grid%ijkDgen

  p        =  region%global%moveGridPower
  eps      =  EPSILON( 1._RFREAL )

! reset motion vectors --------------------------------------------------------

  DO ijk=ibn,ien
    xyzOld(XCOORD,ijk) = xyz(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
    xyzOld(YCOORD,ijk) = xyz(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
    xyzOld(ZCOORD,ijk) = xyz(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
  ENDDO

! move block corners locally

  moveBlock = 0
  method    = 1

  IF (moveBlock==1) THEN
    ic(1) = IndIJK(ipnbeg ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
    ic(2) = IndIJK(ipnbeg ,jpnbeg ,kpnend ,iNOff,ijNOff)
    ic(3) = IndIJK(ipnbeg ,jpnend ,kpnend ,iNOff,ijNOff)
    ic(4) = IndIJK(ipnbeg ,jpnend ,kpnbeg ,iNOff,ijNOff)
    ic(5) = IndIJK(ipnend ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
    ic(6) = IndIJK(ipnend ,jpnbeg ,kpnend ,iNOff,ijNOff)
    ic(7) = IndIJK(ipnend ,jpnend ,kpnend ,iNOff,ijNOff)
    ic(8) = IndIJK(ipnend ,jpnend ,kpnbeg ,iNOff,ijNOff)

    DO ico = 1,8
      denom(ico) = 0._RFREAL
      DO in = 1,8
        IF (ico/=in) THEN
          dist(in,ico) = SQRT( (xyz(XCOORD,ic(in))-xyz(XCOORD,ic(ico)))**2 + &
                               (xyz(YCOORD,ic(in))-xyz(YCOORD,ic(ico)))**2 + &
                               (xyz(ZCOORD,ic(in))-xyz(ZCOORD,ic(ico)))**2 )
          dist(in,ico) = 1._RFREAL/dist(in,ico)
          denom(ico)   = denom(ico) + dist(in,ico)
        ENDIF
      ENDDO
    ENDDO

    IF (.NOT. bndMoved(1)) THEN
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(1)) = 0._RFREAL
        xyz(YCOORD,ic(1)) = 0._RFREAL
        xyz(ZCOORD,ic(1)) = 0._RFREAL
        DO in=1,8
          IF (in/=1) THEN
            xyz(XCOORD,ic(1))= xyz(XCOORD,ic(1))+dist(in,1)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(1))= xyz(YCOORD,ic(1))+dist(in,1)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(1))= xyz(ZCOORD,ic(1))+dist(in,1)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(1)) = xyz(XCOORD,ic(1))/denom(1)
        xyzOld(YCOORD,ic(1)) = xyz(YCOORD,ic(1))/denom(1)
        xyzOld(ZCOORD,ic(1)) = xyz(ZCOORD,ic(1))/denom(1)
      ENDIF
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(2)) = 0._RFREAL
        xyz(YCOORD,ic(2)) = 0._RFREAL
        xyz(ZCOORD,ic(2)) = 0._RFREAL
        DO in=1,8
          IF (in/=2) THEN
            xyz(XCOORD,ic(2))= xyz(XCOORD,ic(2))+dist(in,2)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(2))= xyz(YCOORD,ic(2))+dist(in,2)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(2))= xyz(ZCOORD,ic(2))+dist(in,2)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(2)) = xyz(XCOORD,ic(2))/denom(2)
        xyzOld(YCOORD,ic(2)) = xyz(YCOORD,ic(2))/denom(2)
        xyzOld(ZCOORD,ic(2)) = xyz(ZCOORD,ic(2))/denom(2)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(3)) = 0._RFREAL
        xyz(YCOORD,ic(3)) = 0._RFREAL
        xyz(ZCOORD,ic(3)) = 0._RFREAL
        DO in=1,8
          IF (in/=3) THEN
            xyz(XCOORD,ic(3))= xyz(XCOORD,ic(3))+dist(in,3)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(3))= xyz(YCOORD,ic(3))+dist(in,3)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(3))= xyz(ZCOORD,ic(3))+dist(in,3)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(3)) = xyz(XCOORD,ic(3))/denom(3)
        xyzOld(YCOORD,ic(3)) = xyz(YCOORD,ic(3))/denom(3)
        xyzOld(ZCOORD,ic(3)) = xyz(ZCOORD,ic(3))/denom(3)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(4)) = 0._RFREAL
        xyz(YCOORD,ic(4)) = 0._RFREAL
        xyz(ZCOORD,ic(4)) = 0._RFREAL
        DO in=1,8
          IF (in/=4) THEN
            xyz(XCOORD,ic(4))= xyz(XCOORD,ic(4))+dist(in,4)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(4))= xyz(YCOORD,ic(4))+dist(in,4)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(4))= xyz(ZCOORD,ic(4))+dist(in,4)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(4)) = xyz(XCOORD,ic(4))/denom(4)
        xyzOld(YCOORD,ic(4)) = xyz(YCOORD,ic(4))/denom(4)
        xyzOld(ZCOORD,ic(4)) = xyz(ZCOORD,ic(4))/denom(4)
      ENDIF
    ENDIF

    IF (.NOT. bndMoved(2)) THEN
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(5)) = 0._RFREAL
        xyz(YCOORD,ic(5)) = 0._RFREAL
        xyz(ZCOORD,ic(5)) = 0._RFREAL
        DO in=1,8
          IF (in/=5) THEN
            xyz(XCOORD,ic(5))= xyz(XCOORD,ic(5))+dist(in,5)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(5))= xyz(YCOORD,ic(5))+dist(in,5)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(5))= xyz(ZCOORD,ic(5))+dist(in,5)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(5)) = xyz(XCOORD,ic(5))/denom(5)
        xyzOld(YCOORD,ic(5)) = xyz(YCOORD,ic(5))/denom(5)
        xyzOld(ZCOORD,ic(5)) = xyz(ZCOORD,ic(5))/denom(5)
      ENDIF
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(6)) = 0._RFREAL
        xyz(YCOORD,ic(6)) = 0._RFREAL
        xyz(ZCOORD,ic(6)) = 0._RFREAL
        DO in=1,8
          IF (in/=6) THEN
            xyz(XCOORD,ic(6))= xyz(XCOORD,ic(6))+dist(in,6)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(6))= xyz(YCOORD,ic(6))+dist(in,6)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(6))= xyz(ZCOORD,ic(6))+dist(in,6)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(6)) = xyz(XCOORD,ic(6))/denom(6)
        xyzOld(YCOORD,ic(6)) = xyz(YCOORD,ic(6))/denom(6)
        xyzOld(ZCOORD,ic(6)) = xyz(ZCOORD,ic(6))/denom(6)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(7)) = 0._RFREAL
        xyz(YCOORD,ic(7)) = 0._RFREAL
        xyz(ZCOORD,ic(7)) = 0._RFREAL
        DO in=1,8
          IF (in/=7) THEN
            xyz(XCOORD,ic(7))= xyz(XCOORD,ic(7))+dist(in,7)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(7))= xyz(YCOORD,ic(7))+dist(in,7)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(7))= xyz(ZCOORD,ic(7))+dist(in,7)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(7)) = xyz(XCOORD,ic(7))/denom(7)
        xyzOld(YCOORD,ic(7)) = xyz(YCOORD,ic(7))/denom(7)
        xyzOld(ZCOORD,ic(7)) = xyz(ZCOORD,ic(7))/denom(7)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(8)) = 0._RFREAL
        xyz(YCOORD,ic(8)) = 0._RFREAL
        xyz(ZCOORD,ic(8)) = 0._RFREAL
        DO in=1,8
          IF (in/=8) THEN
            xyz(XCOORD,ic(8))= xyz(XCOORD,ic(8))+dist(in,8)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(8))= xyz(YCOORD,ic(8))+dist(in,8)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(8))= xyz(ZCOORD,ic(8))+dist(in,8)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(8)) = xyz(XCOORD,ic(8))/denom(8)
        xyzOld(YCOORD,ic(8)) = xyz(YCOORD,ic(8))/denom(8)
        xyzOld(ZCOORD,ic(8)) = xyz(ZCOORD,ic(8))/denom(8)
      ENDIF
    ENDIF   ! bndMoved
  ENDIF     ! moveBlock corners

! compute new coordinates -----------------------------------------------------

  IF (method==1) THEN

! - New weighted Laplacian Smoothing of Bono Wasistho:

    DO k=kpnbeg,kpnend
     DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ip1  = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ip2  = IndIJK(i+2,j  ,k  ,iNOff,ijNOff)
        im1  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        im2  = IndIJK(i-2,j  ,k  ,iNOff,ijNOff)
        jp1  = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        jp2  = IndIJK(i  ,j+2,k  ,iNOff,ijNOff)
        jm1  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        jm2  = IndIJK(i  ,j-2,k  ,iNOff,ijNOff)
        kp1  = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        kp2  = IndIJK(i  ,j  ,k+2,iNOff,ijNOff)
        km1  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        km2  = IndIJK(i  ,j  ,k-2,iNOff,ijNOff)

        sim = SQRT( (xyzOrig(XCOORD,im1)+xyzOld(XCOORD,im1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,im1)+xyzOld(YCOORD,im1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,im1)+xyzOld(ZCOORD,im1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p       
        sip = SQRT( (xyzOrig(XCOORD,ip1)+xyzOld(XCOORD,ip1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,ip1)+xyzOld(YCOORD,ip1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,ip1)+xyzOld(ZCOORD,ip1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p

        sjm = SQRT( (xyzOrig(XCOORD,jm1)+xyzOld(XCOORD,jm1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,jm1)+xyzOld(YCOORD,jm1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,jm1)+xyzOld(ZCOORD,jm1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p       
        sjp = SQRT( (xyzOrig(XCOORD,jp1)+xyzOld(XCOORD,jp1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,jp1)+xyzOld(YCOORD,jp1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,jp1)+xyzOld(ZCOORD,jp1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p

        skm = SQRT( (xyzOrig(XCOORD,km1)+xyzOld(XCOORD,km1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,km1)+xyzOld(YCOORD,km1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,km1)+xyzOld(ZCOORD,km1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p       
        skp = SQRT( (xyzOrig(XCOORD,kp1)+xyzOld(XCOORD,kp1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,kp1)+xyzOld(YCOORD,kp1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,kp1)+xyzOld(ZCOORD,kp1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p

        sim = 1._RFREAL/(sim+eps)       
        sip = 1._RFREAL/(sip+eps)
        sjm = 1._RFREAL/(sjm+eps)       
        sjp = 1._RFREAL/(sjp+eps)
        skm = 1._RFREAL/(skm+eps)       
        skp = 1._RFREAL/(skp+eps)

        rxi = sim*xyzOld(XCOORD,im1) + sip*xyzOld(XCOORD,ip1)
        ryi = sim*xyzOld(YCOORD,im1) + sip*xyzOld(YCOORD,ip1)
        rzi = sim*xyzOld(ZCOORD,im1) + sip*xyzOld(ZCOORD,ip1)

        rxj = sjm*xyzOld(XCOORD,jm1) + sjp*xyzOld(XCOORD,jp1)
        ryj = sjm*xyzOld(YCOORD,jm1) + sjp*xyzOld(YCOORD,jp1)
        rzj = sjm*xyzOld(ZCOORD,jm1) + sjp*xyzOld(ZCOORD,jp1)

        rxk = skm*xyzOld(XCOORD,km1) + skp*xyzOld(XCOORD,kp1)
        ryk = skm*xyzOld(YCOORD,km1) + skp*xyzOld(YCOORD,kp1)
        rzk = skm*xyzOld(ZCOORD,km1) + skp*xyzOld(ZCOORD,kp1)

        d = 1._RFREAL/(sim+sip+sjm+sjp+skm+skp)

!        xyz(XCOORD,ijk) = (rxi+rxj+rxk)*d
!        xyz(YCOORD,ijk) = (ryi+ryj+ryk)*d
!        xyz(ZCOORD,ijk) = (rzi+rzj+rzk)*d

        xyz(XCOORD,ijk) = (1._RFREAL-REAL( idgen(ijk) ))*(rxi+rxj+rxk)*d + &
                          REAL( idgen(ijk) )*(xyz(XCOORD,ijk)-xyzOrig(XCOORD,ijk))
        xyz(YCOORD,ijk) = (1._RFREAL-REAL( idgen(ijk) ))*(ryi+ryj+ryk)*d + &
                          REAL( idgen(ijk) )*(xyz(YCOORD,ijk)-xyzOrig(YCOORD,ijk))
        xyz(ZCOORD,ijk) = (1._RFREAL-REAL( idgen(ijk) ))*(rzi+rzj+rzk)*d + &
                          REAL( idgen(ijk) )*(xyz(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))

      ENDDO   ! i
     ENDDO    ! j
    ENDDO     ! k

  ELSEIF (method==2) THEN

! - Old method by Jiri Blazek

    DO k=kpnbeg,kpnend
     DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ip1  = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ip2  = IndIJK(i+2,j  ,k  ,iNOff,ijNOff)
        im1  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        im2  = IndIJK(i-2,j  ,k  ,iNOff,ijNOff)
        jp1  = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        jp2  = IndIJK(i  ,j+2,k  ,iNOff,ijNOff)
        jm1  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        jm2  = IndIJK(i  ,j-2,k  ,iNOff,ijNOff)
        kp1  = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        kp2  = IndIJK(i  ,j  ,k+2,iNOff,ijNOff)
        km1  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        km2  = IndIJK(i  ,j  ,k-2,iNOff,ijNOff)

        rxi = xyzOld(XCOORD,im1) + xyzOld(XCOORD,ip1)
        ryi = xyzOld(YCOORD,im1) + xyzOld(YCOORD,ip1)
        rzi = xyzOld(ZCOORD,im1) + xyzOld(ZCOORD,ip1)

        rxj = xyzOld(XCOORD,jm1) + xyzOld(XCOORD,jp1)
        ryj = xyzOld(YCOORD,jm1) + xyzOld(YCOORD,jp1)
        rzj = xyzOld(ZCOORD,jm1) + xyzOld(ZCOORD,jp1)

        rxk = xyzOld(XCOORD,km1) + xyzOld(XCOORD,kp1)
        ryk = xyzOld(YCOORD,km1) + xyzOld(YCOORD,kp1)
        rzk = xyzOld(ZCOORD,km1) + xyzOld(ZCOORD,kp1)

        si = (xyzOld(XCOORD,ip1)+xyzOrig(XCOORD,ip1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,ip1)+xyzOrig(YCOORD,ip1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,ip1)+xyzOrig(ZCOORD,ip1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,im1)+xyzOrig(XCOORD,im1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,im1)+xyzOrig(YCOORD,im1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,im1)+xyzOrig(ZCOORD,im1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,ip2)+xyzOrig(XCOORD,ip2)- &
              xyzOld(XCOORD,ip1)-xyzOrig(XCOORD,ip1))**2 + &
             (xyzOld(YCOORD,ip2)+xyzOrig(YCOORD,ip2)- &
              xyzOld(YCOORD,ip1)-xyzOrig(YCOORD,ip1))**2 + &
             (xyzOld(ZCOORD,ip2)+xyzOrig(ZCOORD,ip2)- &
              xyzOld(ZCOORD,ip1)-xyzOrig(ZCOORD,ip1))**2 + &
             (xyzOld(XCOORD,im2)+xyzOrig(XCOORD,im2)- &
              xyzOld(XCOORD,im1)-xyzOrig(XCOORD,im1))**2 + &
             (xyzOld(YCOORD,im2)+xyzOrig(YCOORD,im2)- &
              xyzOld(YCOORD,im1)-xyzOrig(YCOORD,im1))**2 + &
             (xyzOld(ZCOORD,im2)+xyzOrig(ZCOORD,im2)- &
              xyzOld(ZCOORD,im1)-xyzOrig(ZCOORD,im1))**2

        sj = (xyzOld(XCOORD,jp1)+xyzOrig(XCOORD,jp1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,jp1)+xyzOrig(YCOORD,jp1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,jp1)+xyzOrig(ZCOORD,jp1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,jm1)+xyzOrig(XCOORD,jm1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,jm1)+xyzOrig(YCOORD,jm1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,jm1)+xyzOrig(ZCOORD,jm1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,jp2)+xyzOrig(XCOORD,jp2)- &
              xyzOld(XCOORD,jp1)-xyzOrig(XCOORD,jp1))**2 + &
             (xyzOld(YCOORD,jp2)+xyzOrig(YCOORD,jp2)- &
              xyzOld(YCOORD,jp1)-xyzOrig(YCOORD,jp1))**2 + &
             (xyzOld(ZCOORD,jp2)+xyzOrig(ZCOORD,jp2)- &
              xyzOld(ZCOORD,jp1)-xyzOrig(ZCOORD,jp1))**2 + &
             (xyzOld(XCOORD,jm2)+xyzOrig(XCOORD,jm2)- &
              xyzOld(XCOORD,jm1)-xyzOrig(XCOORD,jm1))**2 + &
             (xyzOld(YCOORD,jm2)+xyzOrig(YCOORD,jm2)- &
              xyzOld(YCOORD,jm1)-xyzOrig(YCOORD,jm1))**2 + &
             (xyzOld(ZCOORD,jm2)+xyzOrig(ZCOORD,jm2)- &
              xyzOld(ZCOORD,jm1)-xyzOrig(ZCOORD,jm1))**2

        sk = (xyzOld(XCOORD,kp1)+xyzOrig(XCOORD,kp1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,kp1)+xyzOrig(YCOORD,kp1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,kp1)+xyzOrig(ZCOORD,kp1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,km1)+xyzOrig(XCOORD,km1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,km1)+xyzOrig(YCOORD,km1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,km1)+xyzOrig(ZCOORD,km1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,kp2)+xyzOrig(XCOORD,kp2)- &
              xyzOld(XCOORD,kp1)-xyzOrig(XCOORD,kp1))**2 + &
             (xyzOld(YCOORD,kp2)+xyzOrig(YCOORD,kp2)- &
              xyzOld(YCOORD,kp1)-xyzOrig(YCOORD,kp1))**2 + &
             (xyzOld(ZCOORD,kp2)+xyzOrig(ZCOORD,kp2)- &
              xyzOld(ZCOORD,kp1)-xyzOrig(ZCOORD,kp1))**2 + &
             (xyzOld(XCOORD,km2)+xyzOrig(XCOORD,km2)- &
              xyzOld(XCOORD,km1)-xyzOrig(XCOORD,km1))**2 + &
             (xyzOld(YCOORD,km2)+xyzOrig(YCOORD,km2)- &
              xyzOld(YCOORD,km1)-xyzOrig(YCOORD,km1))**2 + &
             (xyzOld(ZCOORD,km2)+xyzOrig(ZCOORD,km2)- &
              xyzOld(ZCOORD,km1)-xyzOrig(ZCOORD,km1))**2

        wi  = 1._RFREAL/(si)**p
        wj  = 1._RFREAL/(sj)**p
        wk  = 1._RFREAL/(sk)**p
        d   = 2._RFREAL*(wi+wj+wk)

!        xyz(XCOORD,ijk) = (wi*rxi+wj*rxj+wk*rxk)/d
!        xyz(YCOORD,ijk) = (wi*ryi+wj*ryj+wk*ryk)/d
!        xyz(ZCOORD,ijk) = (wi*rzi+wj*rzj+wk*rzk)/d

        xyz(XCOORD,ijk) = (1._RFREAL-REAL( idgen(ijk) ))* &
                          (wi*rxi+wj*rxj+wk*rxk)/d + &
                          REAL( idgen(ijk) )*(xyz(XCOORD,ijk)-xyzOrig(XCOORD,ijk))
        xyz(YCOORD,ijk) = (1._RFREAL-REAL( idgen(ijk) ))* &
                          (wi*ryi+wj*ryj+wk*ryk)/d + &
                          REAL( idgen(ijk) )*(xyz(YCOORD,ijk)-xyzOrig(YCOORD,ijk))
        xyz(ZCOORD,ijk) = (1._RFREAL-REAL( idgen(ijk) ))* &
                          (wi*rzi+wj*rzj+wk*rzk)/d + &
                          REAL( idgen(ijk) )*(xyz(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))
      ENDDO   ! i
     ENDDO    ! j
    ENDDO     ! k
  ENDIF       ! method
!if (region%iRegionGlobal==13) write(*,*)'xyzOrth0',ipnbeg,ipnend,jpnbeg,jpnend,kpnbeg,kpnend,xyz(XCOORD,IndIJK(1 ,1 ,23 ,iNOff,ijNOff))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_LaplaceGridSolve


!******************************************************************************
!
! Purpose: Enforce orthogonal Laplacian smoothing
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates, and 
!                 previous grid movements
!
! Output: region%levels%grid%xyz = new grid movements.
!
! Notes: on entry xyz holds the actual grid, on exit xyz holds the grid motion
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridOrtho( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensDummyNodes, &
                            RFLO_GetDimensPhysNodes, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset, &
                            RFLO_CalcFaceVectors
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, ibn, ien, ibc, iec, iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijk, ijkC, ijkN, ijkm, ijkp, ijka(3), errorFlag
  INTEGER :: ijkN1, ijkN2, ijkN3, ijkN4, ijkN5, ijkN6, ijkN7, ijkN8
  INTEGER :: ijkp1, ijkm1, ijkp2, ijkm2, ijkp3, ijkm3, ijkp4, ijkm4
  INTEGER, POINTER :: idgen(:)

  REAL(RFREAL) :: rnfac, ocell, eps, oo6, oo8, smod, fac1, fac2
  REAL(RFREAL) :: usf(XCOORD:ZCOORD), fc(XCOORD:ZCOORD)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOrig(:,:), xyzOld(:,:), xyzOrth(:,:)
  REAL(RFREAL), POINTER :: xyzTemp(:,:), si(:,:), sj(:,:), sk(:,:)
!  REAL(RFREAL), POINTER :: cofg(:,:)   ! *
  REAL(RFREAL), ALLOCATABLE :: cofg(:,:)
  REAL(RFREAL), ALLOCATABLE :: xyzWork(:,:), xyzBuf(:,:)
  REAL(RFREAL), ALLOCATABLE :: fci(:,:), fcj(:,:), fck(:,:)

  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_LaplaceGridOrtho',&
       'RFLO_ModLaplaceSmoothing.F90' )

! get local parameters --------------------------------------------------------

  rnfac = 1._RFREAL/12._RFREAL
  oo6   = 1._RFREAL/6._RFREAL
  oo8   = 1._RFREAL/8._RFREAL
  ocell = global%moveGridOrthCell
  eps   = 100._RFREAL*EPSILON(1._RFREAL)
  fac1  = region%global%moveGridWeight
  fac2  = 1._RFREAL-fac1

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  xyz     => region%levels(iLev)%grid%xyz
  xyzOrig => region%levels(iLev)%gridOld%xyz
  xyzOld  => region%levels(iLev)%grid%xyzOld
  xyzOrth => region%levels(iLev)%grid%xyzOrth
  xyzTemp => region%levels(iLev)%grid%xyzTemp
  si      => region%levels(iLev)%grid%si
  sj      => region%levels(iLev)%grid%sj
  sk      => region%levels(iLev)%grid%sk
!  cofg    => region%levels(iLev)%grid%cofg            ! *
  idgen   => region%levels(iLev)%grid%ijkDgen

  ALLOCATE( xyzWork(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( xyzBuf(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( fci(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( fcj(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( fck(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( cofg(3,ibc:iec),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  xyzWork = 0._RFREAL

! reset motion vectors --------------------------------------------------------

  DO ijk=ibn,ien
    xyzOld(XCOORD,ijk) = xyz(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
    xyzOld(YCOORD,ijk) = xyz(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
    xyzOld(ZCOORD,ijk) = xyz(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
  ENDDO

!  xyzTemp = xyz
  xyzTemp = xyzOrig + xyzOld
  xyzBuf  = xyz

! get cell center -------------------------------------------------------------

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC     = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkN1    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        ijkN3    = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
        ijkN4    = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkN5    = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijkN6    = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
        ijkN7    = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
        ijkN8    = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
        cofg(:,ijkC) = xyzTemp(:,ijkN1)+xyzTemp(:,ijkN2)+xyzTemp(:,ijkN3)+ &
                       xyzTemp(:,ijkN4)+xyzTemp(:,ijkN5)+xyzTemp(:,ijkN6)+ &
                       xyzTemp(:,ijkN7)+xyzTemp(:,ijkN8)
        cofg(:,ijkC) = oo8*cofg(:,ijkC)
      ENDDO
    ENDDO
  ENDDO

! perturb grid ----------------------------------------------------------------
        
  DO k=kdnbeg+1,kdnend-1
    DO j=jdnbeg+1,jdnend-1
      DO i=idnbeg+1,idnend-1
        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN1    = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkN3    = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        ijkN4    = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ijkN5    = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijkN6    = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        xyzOrth(:,ijkN) = xyzTemp(:,ijkN1)+xyzTemp(:,ijkN2)+xyzTemp(:,ijkN3)+ &
                          xyzTemp(:,ijkN4)+xyzTemp(:,ijkN5)+xyzTemp(:,ijkN6)
        xyz(:,ijkN)     = oo6*xyzOrth(:,ijkN)
      ENDDO
    ENDDO
  ENDDO

! get face vectors ------------------------------------------------------------

  CALL RFLO_CalcFaceVectors( region )

! weighted face-centers -------------------------------------------------------

  DO k=kdnbeg+1,kdnend-1
    DO j=jdnbeg+1,jdnend-1
      DO i=idnbeg+1,idnbeg+1+(idnend-1-idnbeg-1)/2
        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkp     = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm     = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        fci(:,ijkN) = fac1*cofg(:,ijkm) + fac2*cofg(:,ijkp)
      ENDDO
      DO i=idnbeg+1+(idnend-1-idnbeg-1)/2+1,idnend-1
        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkp     = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm     = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        fci(:,ijkN) = fac2*cofg(:,ijkm) + fac1*cofg(:,ijkp)
      ENDDO
    ENDDO
  ENDDO
  DO k=kdnbeg+1,kdnend-1
    DO j=jdnbeg+1,jdnbeg+1+(jdnend-1-jdnbeg-1)/2
      DO i=idnbeg+1,idnend-1
        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkp     = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm     = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        fcj(:,ijkN) = fac1*cofg(:,ijkm) + fac2*cofg(:,ijkp)
      ENDDO
    ENDDO
    DO j=jdnbeg+1+(jdnend-1-jdnbeg-1)/2+1,jdnend-1
      DO i=idnbeg+1,idnend-1
        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkp     = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm     = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        fcj(:,ijkN) = fac2*cofg(:,ijkm) + fac1*cofg(:,ijkp)
      ENDDO
    ENDDO
  ENDDO
  DO k=kdnbeg+1,kdnbeg+1+(kdnend-1-kdnbeg-1)/2
    DO j=jdnbeg+1,jdnend-1
      DO i=idnbeg+1,idnend-1
        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkp     = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm     = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        fck(:,ijkN) = fac1*cofg(:,ijkm) + fac2*cofg(:,ijkp)
      ENDDO
    ENDDO
  ENDDO
  DO k=kdnbeg+1+(kdnend-1-kdnbeg-1)/2+1,kdnend-1
    DO j=jdnbeg+1,jdnend-1
      DO i=idnbeg+1,idnend-1
        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkp     = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm     = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        fck(:,ijkN) = fac2*cofg(:,ijkm) + fac1*cofg(:,ijkp)
      ENDDO
    ENDDO
  ENDDO

! orthogonal Laplacian --------------------------------------------------------

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend

! ----- I-Face

        ijkN1    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijkN3    = IndIJK(i  ,j-1,k-1,iNOff,ijNOff)
        ijkN4    = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        ijkp1    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm1    = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        ijkp2    = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        ijkm2    = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
        ijkp3    = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
        ijkm3    = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
        ijkp4    = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        ijkm4    = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
!IF (i==2.AND.j==2.AND.k==2) write(*,*)'orthoI'
        CALL AccumulateXyzWork( si, &
             fci(:,ijkN1),fci(:,ijkN2),fci(:,ijkN3),fci(:,ijkN4) )

! ----- J-Face

        ijkN1    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        ijkN3    = IndIJK(i-1,j  ,k-1,iNOff,ijNOff)
        ijkN4    = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ijkp1    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm1    = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        ijkp2    = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        ijkm2    = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
        ijkp3    = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
        ijkm3    = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
        ijkp4    = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        ijkm4    = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
!IF (i==2.AND.j==2.AND.k==2) write(*,*)'orthoJ'
        CALL AccumulateXyzWork( sj, &
             fcj(:,ijkN1),fcj(:,ijkN2),fcj(:,ijkN3),fcj(:,ijkN4) )

! ----- K-Face

        ijkN1    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ijkN3    = IndIJK(i-1,j-1,k  ,iNOff,ijNOff)
        ijkN4    = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijkp1    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkm1    = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        ijkp2    = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        ijkm2    = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
        ijkp3    = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
        ijkm3    = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
        ijkp4    = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        ijkm4    = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
!IF (i==2.AND.j==2.AND.k==2) write(*,*)'orthoK'
        CALL AccumulateXyzWork( sk, &
             fck(:,ijkN1),fck(:,ijkN2),fck(:,ijkN3),fck(:,ijkN4) )

      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! average accumulated grid coordinates

  xyz     = xyzBuf-xyzOrig
  ijka(:) = 1

  DO k=kpnbeg+ijka(3),kpnend-ijka(3)
    DO j=jpnbeg+ijka(2),jpnend-ijka(2)
      DO i=ipnbeg+ijka(1),ipnend-ijka(1)
        ijkN = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        xyzOrth(:,ijkN) = ocell*rnfac*xyzWork(:,ijkN)
!                          +(1._RFREAL-ocell)*xyz(:,ijkN)
        xyz(:,ijkN)     = xyzOrth(:,ijkN) 
!        xyz(:,ijkN)     = (1._RFREAL-REAL( idgen(ijkN) ))*xyzOrth(:,ijkN) + &
!                          REAL( idgen(ijkN) )*xyz(:,ijkN)
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! deallocate temporary memory -------------------------------------------------

  DEALLOCATE( xyzWork,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( xyzBuf,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( fci,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( fcj,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( fck,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( cofg,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

CONTAINS

! ==============================================================================
! subroutine for accumulation of projected cell corners coordinates
! ==============================================================================

  SUBROUTINE AccumulateXyzWork( sf,fc1,fc2,fc3,fc4 )

    IMPLICIT NONE
    REAL(RFREAL), POINTER :: sf(:,:)
    REAL(RFREAL) :: fc1(XCOORD:ZCOORD), fc2(XCOORD:ZCOORD)
    REAL(RFREAL) :: fc3(XCOORD:ZCOORD), fc4(XCOORD:ZCOORD)

    smod    = SQRT( sf(XCOORD,ijkN1)*sf(XCOORD,ijkN1) + &
                    sf(YCOORD,ijkN1)*sf(YCOORD,ijkN1) + &
                    sf(ZCOORD,ijkN1)*sf(ZCOORD,ijkN1) )
    usf(:)  = sf(:,ijkN1)/(smod+eps)
!    fc1(:)   = 0.5_RFREAL*(cofg(:,ijkp1)+cofg(:,ijkm1))
    xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + &
                      DOT_PRODUCT( fc1(:)-xyzTemp(:,ijkN1), usf(:) )*usf(:)
!IF (i==2.AND.j==2.AND.k==2) write(*,*)'ortho1',xyzWork(:,ijkN1),usf(:)

    smod    = SQRT( sf(XCOORD,ijkN2)*sf(XCOORD,ijkN2) + &
                    sf(YCOORD,ijkN2)*sf(YCOORD,ijkN2) + &
                    sf(ZCOORD,ijkN2)*sf(ZCOORD,ijkN2) )
    usf(:)  = sf(:,ijkN2)/(smod+eps)
!    fc2(:)   = 0.5_RFREAL*(cofg(:,ijkp2)+cofg(:,ijkm2))
    xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + &
                      DOT_PRODUCT( fc2(:)-xyzTemp(:,ijkN1), usf(:) )*usf(:)
!IF (i==2.AND.j==2.AND.k==2) write(*,*)'ortho2',xyzWork(:,ijkN1),usf(:)

    smod    = SQRT( sf(XCOORD,ijkN3)*sf(XCOORD,ijkN3) + &
                    sf(YCOORD,ijkN3)*sf(YCOORD,ijkN3) + &
                    sf(ZCOORD,ijkN3)*sf(ZCOORD,ijkN3) )
    usf(:)  = sf(:,ijkN3)/(smod+eps)
!    fc3(:)   = 0.5_RFREAL*(cofg(:,ijkp3)+cofg(:,ijkm3))
    xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + &
                      DOT_PRODUCT( fc3(:)-xyzTemp(:,ijkN1), usf(:) )*usf(:)
!IF (i==2.AND.j==2.AND.k==2) write(*,*)'ortho3',xyzWork(:,ijkN1),usf(:)

    smod    = SQRT( sf(XCOORD,ijkN4)*sf(XCOORD,ijkN4) + &
                    sf(YCOORD,ijkN4)*sf(YCOORD,ijkN4) + &
                    sf(ZCOORD,ijkN4)*sf(ZCOORD,ijkN4) )
    usf(:)  = sf(:,ijkN4)/(smod+eps)
!    fc4(:)   = 0.5_RFREAL*(cofg(:,ijkp4)+cofg(:,ijkm4))
    xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + &
                      DOT_PRODUCT( fc4(:)-xyzTemp(:,ijkN1), usf(:) )*usf(:)
!IF (i==2.AND.j==2.AND.k==2) write(*,*)'ortho4',xyzWork(:,ijkN1),usf(:)

  END SUBROUTINE AccumulateXyzWork

END SUBROUTINE RFLO_LaplaceGridOrtho

!******************************************************************************
!
! Purpose: Enforce orthogonality to Laplacian smoothing
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates, and 
!                 previous grid movements
!
! Output: region%levels%grid%xyz = new grid movements.
!
! Notes: on entry xyz holds the actual grid, on exit xyz holds the grid motion
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridOrtho1( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
                            RFLO_GetDimensPhysNodes, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, ibn, ien, iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijk, ijkC1, ijkC2, ijkC3, ijkC4, ijkN1, ijkN2, ijkN3, ijkN4
  INTEGER :: ijkC1m, ijkC2m, ijkC3m, ijkC4m, ia, ja, ka, errorFlag

  REAL(RFREAL) :: rnfac, ocell, eps
  REAL(RFREAL) :: cofgp(XCOORD:ZCOORD), cofgm(XCOORD:ZCOORD), cfc(XCOORD:ZCOORD)
  REAL(RFREAL) :: xyzc(XCOORD:ZCOORD), xyzp(XCOORD:ZCOORD)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOrig(:,:), xyzOld(:,:), xyzOrth(:,:)
  REAL(RFREAL), POINTER :: cofg(:,:), cfcI(:,:), cfcJ(:,:), cfcK(:,:)
  REAL(RFREAL), ALLOCATABLE :: xyzWork(:,:)

  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_LaplaceGridOrtho1',&
       'RFLO_ModLaplaceSmoothing.F90' )

! get local parameters --------------------------------------------------------

  rnfac = 1._RFREAL/12._RFREAL
  ocell = global%moveGridOrthCell
  eps   = 100._RFREAL*EPSILON(1._RFREAL)

! get dimensions and pointers -------------------------------------------------

  iLev = 1
  ia   = 1
  ja   = 1
  ka   = 1

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  xyz     => region%levels(iLev)%grid%xyz
  xyzOrig => region%levels(iLev)%gridOld%xyz
  xyzOld  => region%levels(iLev)%grid%xyzOld
  xyzOrth => region%levels(iLev)%grid%xyzOrth
  cofg    => region%levels(iLev)%grid%cofg
  cfcI    => region%levels(iLev)%grid%cfcI
  cfcJ    => region%levels(iLev)%grid%cfcJ
  cfcK    => region%levels(iLev)%grid%cfcK

  ALLOCATE( xyzWork(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  xyzWork = 0._RFREAL

! reset motion vectors --------------------------------------------------------

!  xyzOrth = xyz
  xyzOrth = xyz + xyzOrig

!  DO ijk=ibn,ien
!    xyzOld(XCOORD,ijk)  = xyz(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
!    xyzOld(YCOORD,ijk)  = xyz(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
!    xyzOld(ZCOORD,ijk)  = xyz(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
!    xyz(XCOORD,ijk)  = xyz(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
!    xyz(YCOORD,ijk)  = xyz(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
!    xyz(ZCOORD,ijk)  = xyz(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
!  ENDDO

!if (region%iRegionGlobal==1) write(*,*)'xyzOrth0',xyzOrth(:,IndIJK(2,4,2 ,iNOff,ijNOff))-xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))

! reset motion vectors --------------------------------------------------------

  DO k=kpnbeg+ka,kpnend-ka
    DO j=jpnbeg+ja,jpnend-ja
      DO i=ipnbeg+ia,ipnend-ia

! ----- I-Face

        ijkC1    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkC1m   = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        ijkC2    = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        ijkC2m   = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
        ijkC3    = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
        ijkC3m   = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
        ijkC4    = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        ijkC4m   = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
        ijkN1    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        ijkN3    = IndIJK(i  ,j-1,k-1,iNOff,ijNOff)
        ijkN4    = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)

        CALL AccumulateXyzWork( cfcI )
!if (i==2 .AND. j==4 .AND. k==2) write(*,*)'xyzOrthI',xyzWork(:,IndIJK(2,4,2 ,iNOff,ijNOff))/4-xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))

! ----- J-Face

        ijkC1    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkC1m   = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        ijkC2    = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        ijkC2m   = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
        ijkC3    = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
        ijkC3m   = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
        ijkC4    = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        ijkC4m   = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
        ijkN1    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ijkN3    = IndIJK(i-1,j  ,k-1,iNOff,ijNOff)
        ijkN4    = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

        CALL AccumulateXyzWork( cfcJ )
!if (i==2 .AND. j==4 .AND. k==2) write(*,*)'xyzOrthJ',xyzWork(:,IndIJK(2,4,2 ,iNOff,ijNOff))/8-xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))

! ----- K-Face

        ijkC1    = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        ijkC1m   = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        ijkC2    = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        ijkC2m   = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
        ijkC3    = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
        ijkC3m   = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
        ijkC4    = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        ijkC4m   = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
        ijkN1    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijkN3    = IndIJK(i-1,j-1,k  ,iNOff,ijNOff)
        ijkN4    = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)

        CALL AccumulateXyzWork( cfcK )
!if (i==2 .AND. j==4 .AND. k==2) write(*,*)'xyzOrthK',xyzWork(:,IndIJK(2,4,2 ,iNOff,ijNOff))/12-xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))

! ----- average accumulated grid coordinates

!if (region%iRegionGlobal==13 .AND. i==10 .AND. j==10 .AND. k==10) THEN
!        xyzc(:)  = rnfac*xyzWork(:,ijkN1)
!        xyz(XCOORD,ijkN1) = xyzc(XCOORD) - xyzOrig(XCOORD,ijkN1)
!        xyz(YCOORD,ijkN1) = xyzc(YCOORD) - xyzOrig(YCOORD,ijkN1)
!        xyz(ZCOORD,ijkN1) = xyzc(ZCOORD) - xyzOrig(ZCOORD,ijkN1)
!endif
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! average accumulated grid coordinates

  DO k=kpnbeg+ka,kpnend-ka
    DO j=jpnbeg+ja,jpnend-ja
      DO i=ipnbeg+ia,ipnend-ia
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        xyzOrth(:,ijk)  = rnfac*xyzWork(:,ijk) - xyzOrig(:,ijk)
        xyz(XCOORD,ijk) = xyzOrth(XCOORD,ijk)
        xyz(YCOORD,ijk) = xyzOrth(YCOORD,ijk)
        xyz(ZCOORD,ijk) = xyzOrth(ZCOORD,ijk)
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

!if (region%iRegionGlobal==1) write(*,*)'xyzOrth5',xyzOrth(:,IndIJK(2,4,2 ,iNOff,ijNOff))+xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))
! deallocate temporary memory -------------------------------------------------

  DEALLOCATE( xyzWork,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

CONTAINS

! ==============================================================================
! subroutine for accumulation of projected cell corners coordinates
! ==============================================================================

  SUBROUTINE AccumulateXyzWork( cfac )

    IMPLICIT NONE
    REAL(RFREAL), POINTER :: cfac(:,:)

        cofgp(:) = cofg(:,ijkC1)
        cofgm(:) = cofg(:,ijkC1m)
        cfc(:)   = cfac(:,ijkN1)
!        cfc(:)   = 0.5_RFREAL*(cofgp(:)+cofgm(:))
        xyzc(:)  = xyzOrth(:,ijkN1)
        CALL RFLO_ProjectQuadCorner( ocell,eps,cofgp,cofgm,cfc,xyzc,xyzp )
        xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + xyzp(:)
!write(*,*)'xyzOrth1',xyzp(:),xyzc(:)
!write(*,*)'xyzOrth1',xyzWork(:,ijkN1)
        cofgp(:) = cofg(:,ijkC2)
        cofgm(:) = cofg(:,ijkC2m)
        cfc(:)   = cfac(:,ijkN2)
!        cfc(:)   = 0.5_RFREAL*(cofgp(:)+cofgm(:))
        xyzc(:)  = xyzOrth(:,ijkN1)
        CALL RFLO_ProjectQuadCorner( ocell,eps,cofgp,cofgm,cfc,xyzc,xyzp )
        xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + xyzp(:)
!write(*,*)'xyzOrth2',xyzp(:),xyzc(:)
!write(*,*)'xyzOrth2',xyzWork(:,ijkN1)
        cofgp(:) = cofg(:,ijkC3)
        cofgm(:) = cofg(:,ijkC3m)
        cfc(:)   = cfac(:,ijkN3)
!        cfc(:)   = 0.5_RFREAL*(cofgp(:)+cofgm(:))
        xyzc(:)  = xyzOrth(:,ijkN1)
        CALL RFLO_ProjectQuadCorner( ocell,eps,cofgp,cofgm,cfc,xyzc,xyzp )
        xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + xyzp(:)
!write(*,*)'xyzOrth3',xyzp(:),xyzc(:)
!write(*,*)'xyzOrth3',xyzWork(:,ijkN1)
        cofgp(:) = cofg(:,ijkC4)
        cofgm(:) = cofg(:,ijkC4m)
        cfc(:)   = cfac(:,ijkN4)
!        cfc(:)   = 0.5_RFREAL*(cofgp(:)+cofgm(:))
        xyzc(:)  = xyzOrth(:,ijkN1)
        CALL RFLO_ProjectQuadCorner( ocell,eps,cofgp,cofgm,cfc,xyzc,xyzp )
        xyzWork(:,ijkN1) = xyzWork(:,ijkN1) + xyzp(:)
!write(*,*)'xyzOrth4',xyzp(:),xyzc(:)
!write(*,*)'xyzOrth4',xyzWork(:,ijkN1)
  END SUBROUTINE AccumulateXyzWork

END SUBROUTINE RFLO_LaplaceGridOrtho1

!******************************************************************************
!
! Purpose: Enforce orthogonality to Laplacian smoothing
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates, and 
!                 previous grid movements
!
! Output: region%levels%grid%xyz = new grid movements.
!
! Notes: on entry xyz holds the actual grid, on exit xyz holds the grid motion
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridOrtho2( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
                            RFLO_GetDimensPhysNodes, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, ibn, ien, iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijk, ijkN, ijkN1, ijkN2, ijkN3, ijkN4, errorFlag

  REAL(RFREAL) :: rnfac, ocell, eps, a11, a12, a21, a22, rh1, rh2, rja, c1, c2
  REAL(RFREAL) :: delt(XCOORD:ZCOORD), f1(XCOORD:ZCOORD), f2(XCOORD:ZCOORD)
  REAL(RFREAL) :: xyzp(XCOORD:ZCOORD)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOrig(:,:), xyzOld(:,:), xyzOrth(:,:)
  REAL(RFREAL), POINTER :: xyzTemp(:,:)
  REAL(RFREAL), ALLOCATABLE :: xyzWork(:,:)

  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_LaplaceGridOrtho2',&
       'RFLO_ModLaplaceSmoothing.F90' )

! get local parameters --------------------------------------------------------

  rnfac = 1._RFREAL/12._RFREAL
  ocell = global%moveGridOrthCell
  eps   = 100._RFREAL*EPSILON(1._RFREAL)

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  xyz     => region%levels(iLev)%grid%xyz
  xyzOrig => region%levels(iLev)%gridOld%xyz
  xyzOld  => region%levels(iLev)%grid%xyzOld
  xyzOrth => region%levels(iLev)%grid%xyzOrth
  xyzTemp => region%levels(iLev)%grid%xyzTemp

  ALLOCATE( xyzWork(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  xyzWork = 0._RFREAL

! reset motion vectors --------------------------------------------------------

  xyzTemp = xyzOld
  DO ijk=ibn,ien
    xyzOld(:,ijk)  = xyzOld(:,ijk) + xyzOrig(:,ijk)
    xyzOrth(:,ijk) = xyz(   :,ijk) + xyzOrig(:,ijk)
  ENDDO

!if (region%iRegionGlobal==1) write(*,*)'xyzOrth0',xyzOrth(:,IndIJK(2,2,2 ,iNOff,ijNOff))

! reset motion vectors --------------------------------------------------------

  DO k=kpnbeg+1,kpnend-1
    DO j=jpnbeg+1,jpnend-1
      DO i=ipnbeg+1,ipnend-1

! ----- I-Face

        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN1    = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        ijkN3    = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijkN4    = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

        CALL AccumulateXyzWork
!        if (i==2 .AND. j==4 .AND. k==2)  &
!           write(*,*)'xyzOrthI',xyzWork(:,IndIJK(2,4,2 ,iNOff,ijNOff))/4 &
!              -xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))

! ----- J-Face

        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN1    = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        ijkN2    = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijkN3    = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        ijkN4    = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)

        CALL AccumulateXyzWork
!        if (i==2 .AND. j==4 .AND. k==2) &
!           write(*,*)'xyzOrthJ',xyzWork(:,IndIJK(2,4,2 ,iNOff,ijNOff))/8 &
!               -xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))

! ----- K-Face

        ijkN     = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkN1    = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijkN2    = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkN3    = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ijkN4    = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)

        CALL AccumulateXyzWork
!        if (i==2 .AND. j==4 .AND. k==2) &
!           write(*,*)'xyzOrthK',xyzWork(:,IndIJK(2,4,2 ,iNOff,ijNOff))/12 &
!               -xyzOrig(:,IndIJK(2,4,2 ,iNOff,ijNOff))

      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! average accumulated grid coordinates

  DO k=kpnbeg+1,kpnend-1
    DO j=jpnbeg+1,jpnend-1
      DO i=ipnbeg+1,ipnend-1
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        xyzOrth(:,ijk)  = rnfac*xyzWork(:,ijk) - xyzOrig(:,ijk)
        xyz(    :,ijk)  = xyzOrth(:,ijk)
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

  xyzOld = xyzTemp

!if (region%iRegionGlobal==1) write(*,*)'xyzOrth1',xyz(:,IndIJK(2,2,2 ,iNOff,ijNOff)) + xyzOrig(:,IndIJK(2,2,2 ,iNOff,ijNOff))
!if (region%iRegionGlobal==1) write(*,*)'xyzOrth2',xyzOld(:,IndIJK(2,2,2 ,iNOff,ijNOff)) + xyzOrig(:,IndIJK(2,2,2 ,iNOff,ijNOff))

! deallocate temporary memory -------------------------------------------------

  DEALLOCATE( xyzWork,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

CONTAINS

! ==============================================================================
! subroutine for accumulation of projected cell corners coordinates
! ==============================================================================

  SUBROUTINE AccumulateXyzWork

    IMPLICIT NONE

    delt(:) = xyzOrth(:,ijkN)-xyzOld(:,ijkN)  

    f1(:)   = xyzOld(:,ijkN1)-xyzOld(:,ijkN)
    f2(:)   = xyzOld(:,ijkN2)-xyzOld(:,ijkN)
    a11     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN1) )
    a12     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN1) )
    a21     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN2) )
    a22     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN2) )
    rh1     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN1) )
    rh2     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN2) )
    rja     = 1._RFREAL/(a11*a22-a12*a21+eps)
    c1      = (a22*rh1 - a12*rh2)*rja
    c2      = (a11*rh2 - a21*rh1)*rja
    xyzp(:) = ocell*(xyzOld(:,ijkN) + c1*f1(:) + c2*f2(:)) + &
              (1._RFREAL-ocell)*xyzOrth(:,ijkN)
    xyzWork(:,ijkN) = xyzWork(:,ijkN) + xyzp(:)

    f1(:)   = xyzOld(:,ijkN2)-xyzOld(:,ijkN)
    f2(:)   = xyzOld(:,ijkN3)-xyzOld(:,ijkN)
    a11     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN2) )
    a12     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN2) )
    a21     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN3) )
    a22     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN3) )
    rh1     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN2) )
    rh2     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN3) )
    rja     = 1._RFREAL/(a11*a22-a12*a21+eps)
    c1      = (a22*rh1 - a12*rh2)*rja
    c2      = (a11*rh2 - a21*rh1)*rja
    xyzp(:) = ocell*(xyzOld(:,ijkN) + c1*f1(:) + c2*f2(:)) + &
              (1._RFREAL-ocell)*xyzOrth(:,ijkN)
    xyzWork(:,ijkN) = xyzWork(:,ijkN) + xyzp(:)

    f1(:)   = xyzOld(:,ijkN3)-xyzOld(:,ijkN)
    f2(:)   = xyzOld(:,ijkN4)-xyzOld(:,ijkN)
    a11     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN3) )
    a12     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN3) )
    a21     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN4) )
    a22     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN4) )
    rh1     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN3) )
    rh2     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN4) )
    rja     = 1._RFREAL/(a11*a22-a12*a21+eps)
    c1      = (a22*rh1 - a12*rh2)*rja
    c2      = (a11*rh2 - a21*rh1)*rja
    xyzp(:) = ocell*(xyzOld(:,ijkN) + c1*f1(:) + c2*f2(:)) + &
              (1._RFREAL-ocell)*xyzOrth(:,ijkN)
    xyzWork(:,ijkN) = xyzWork(:,ijkN) + xyzp(:)

    f1(:)   = xyzOld(:,ijkN4)-xyzOld(:,ijkN)
    f2(:)   = xyzOld(:,ijkN1)-xyzOld(:,ijkN)
    a11     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN4) )
    a12     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN4) )
    a21     = DOT_PRODUCT( f1(:),xyzOld(:,ijkN1) )
    a22     = DOT_PRODUCT( f2(:),xyzOld(:,ijkN1) )
    rh1     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN4) )
    rh2     = DOT_PRODUCT( delt(:),xyzOld(:,ijkN1) )
    rja     = 1._RFREAL/(a11*a22-a12*a21+eps)
    c1      = (a22*rh1 - a12*rh2)*rja
    c2      = (a11*rh2 - a21*rh1)*rja
    xyzp(:) = ocell*(xyzOld(:,ijkN) + c1*f1(:) + c2*f2(:)) + &
              (1._RFREAL-ocell)*xyzOrth(:,ijkN)
    xyzWork(:,ijkN) = xyzWork(:,ijkN) + xyzp(:)

  END SUBROUTINE AccumulateXyzWork

END SUBROUTINE RFLO_LaplaceGridOrtho2

!******************************************************************************
!
! Purpose:
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates.
!
! Output: region%levels%grid%xyz = grid movements.
!
! Notes: on entry and exit, xyz holds the grid motion
!
!******************************************************************************

SUBROUTINE RFLO_ProjectQuadCorner( ocell,eps,cofgp,cofgm,cfc,xyzc,xyzp )

  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: ocell, eps
  REAL(RFREAL) :: cofgp(:), cofgm(:), cfc(:), xyzc(:), xyzp(:)

! ... loop variables

! ... local variables
  REAL(RFREAL) :: cfcgp(XCOORD:ZCOORD), cfcgm(XCOORD:ZCOORD)
  REAL(RFREAL) :: sqrp, sqrm, cp, cm, fac

!******************************************************************************

  cfcgp(:) = cofgp(:) - cfc(:)
  cfcgm(:) = cofgm(:) - cfc(:)

  sqrp = cfcgp(XCOORD)*cfcgp(XCOORD) + &
         cfcgp(YCOORD)*cfcgp(YCOORD) + &
         cfcgp(ZCOORD)*cfcgp(ZCOORD)

  sqrm = cfcgm(XCOORD)*cfcgm(XCOORD) + &
         cfcgm(YCOORD)*cfcgm(YCOORD) + &
         cfcgm(ZCOORD)*cfcgm(ZCOORD)

  cp = DOT_PRODUCT( (cfc(:)-xyzc(:)),cfcgp(:) )/(sqrp + EPSILON( 1._RFREAL ))
  cm = DOT_PRODUCT( (cfc(:)-xyzc(:)),cfcgm(:) )/(sqrm + EPSILON( 1._RFREAL ))

  fac = 0.5_RFREAL
!  IF (cp < eps .OR. cm < eps ) fac = 1._RFREAL
  xyzp(:) = xyzc(:) + ocell*fac*(cp*cfcgp(:) + cm*cfcgm(:))
!  xyzp(:) = xyzc(:) + ocell*(ABS(cp)*cp*cfcgp(:) + ABS(cm)*cm*cfcgm(:))/ &
!                      (ABS(cp)+ABS(cm))

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_ProjectQuadCorner

! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModLaplaceSmoothing

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModLaplaceSmoothing.F90,v $
! Revision 1.13  2009/08/27 14:04:50  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.12  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/03/05 21:51:34  wasistho
! changed computational space coordinates to be based on initial grid
!
! Revision 1.9  2005/11/28 20:39:16  wasistho
! removed file-Id in rflo_projectQuadCorner
!
! Revision 1.8  2005/11/28 18:47:19  wasistho
! outcommented print statement in gridOrtho2
!
! Revision 1.7  2005/11/28 18:29:56  gzheng
! removed duplicated RFLO_CalcFaceVectors in USE statement, which broke gnu f95 compiler, also fixed 3 single line statement which is too long for some compilers.
!
! Revision 1.6  2005/11/17 00:46:26  wasistho
! dont touch degenerate points in moveGridSolve
!
! Revision 1.5  2005/11/16 20:45:30  wasistho
! update on MoveGridFoms (3), still less robust than (2)
!
! Revision 1.4  2005/11/01 09:01:58  wasistho
! added volume TFI for movement jump
!
! Revision 1.3  2005/10/28 07:36:33  wasistho
! gridOrtho is private
!
! Revision 1.2  2005/10/27 06:29:46  wasistho
! increased denominator by epsilon in ProjectQuadCorner
!
! Revision 1.1  2005/10/27 06:02:46  wasistho
! initial import
!
!
!
! ******************************************************************************













