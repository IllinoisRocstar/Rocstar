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
! Purpose: Suite for global grid-motion routines with elliptic PDE smoothings.
!
! Description: The procedure is similar to MOVEGRID_GLOBAL method, except
!              the Laplacian smoothing is replaced by elliptic PDE smoothing.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModMoveGridElliptGlo.F90,v 1.18 2009/08/27 14:04:50 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModMoveGridElliptGlo

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
  PUBLIC :: RFLO_MoveGridElliptGlo, &
            RFLO_MgElliptBndDeformation

! private : RFLO_MgElliptSurfacesGlo
!           RFLO_MgElliptInterfacesGlo
!           RFLO_MgElliptEdgesGlo
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModMoveGridElliptGlo.F90,v $ $Revision: 1.18 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


!******************************************************************************
!
! Purpose: redistribute grid nodes according to the movement of the
!          boundaries. This function smoothes the grid globally by
!          solving an elliptic PDE for node coordinates.
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions%levels%grid%xyz = new grid coordinates.
!
! Notes: grid%xyz temporarily stores nodal displacements. The deformation
!        is applied to the finest grid first.
!
!******************************************************************************

SUBROUTINE RFLO_MoveGridElliptGlo( regions )

  USE ModInterfaces, ONLY : RFLO_MoveGridSurfaces, RFLO_MoveGridInterfaces, &
        RFLO_ChangeInteriorGrid, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_CalcControlVolumes, &
        RFLO_CalcCellCentroids, RFLO_C2fAvgCoeffs, RFLO_C2eAvgCoeffs, &
        RFLO_CheckMetrics, RFLO_CalcGridSpeeds, RFLO_BoundaryDeformation
  USE RFLO_ModGridMetrics,      ONLY : RFLO_GridQualityGlobal
  USE RFLO_ModLaplaceSmoothing, ONLY : RFLO_LaplaceGridSmoo
  USE RFLO_ModElliptSmoothing,  ONLY : RFLO_ElliptGridSmoo, &
                                       RFLO_ElliptGridSmooRegion

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iter, iPatch, ijk

! ... local variables
  LOGICAL :: someMoved

  INTEGER :: bcType

  REAL(RFREAL)          :: resid, globalResid
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

  TYPE(t_grid), POINTER   :: grid, gridOld
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch
#ifdef GENX
  DOUBLE PRECISION :: dAlpha
#endif

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MoveGridElliptGlo',&
       'RFLO_ModMoveGridElliptGlo.F90' )

#ifdef GENX
! update geometry buffers -----------------------------------------------------

  dAlpha = global%dtMin/global%dTimeSystem
  CALL COM_call_function( global%genxHandleGm,1,dAlpha )
#endif

! receive and distribute deformations for each region -------------------------

  CALL RFLO_MgElliptSurfacesGlo( regions,someMoved )

! fix interfaces between regions ----------------------------------------------

  IF (someMoved .eqv. .true.) THEN
    CALL RFLO_MgElliptInterfacesGlo( regions )
  ENDIF

! update grid, dummy, corner and edge cells -----------------------------------

  IF (someMoved) THEN
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE .AND. &           ! on my processor
          regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! ----- change the interior grid

        grid    => regions(iReg)%levels(1)%grid
        gridOld => regions(iReg)%levels(1)%gridOld
        CALL RFLO_ChangeInteriorGrid( regions(iReg),grid%boundMoved, &
                                      grid%edgeMoved,grid%arcLen12, &
                                      grid%arcLen34,grid%arcLen56, &
                                      gridOld%xyzOld,grid%xyz )

! ----- update coarse grids and dummy cells

        CALL RFLO_GenerateCoarseGrids( regions(iReg) )   ! coarsen finest grid
        CALL RFLO_CopyGeometryDummy( regions(iReg) )     ! copy to dummy nodes
        CALL RFLO_ExtrapolateGeometry( regions(iReg) )   ! extrapolate
      ENDIF     ! region on this processor and active, grid moving
    ENDDO       ! iReg
    CALL RFLO_ExchangeGeometry( regions )                ! exchange geometry
  ENDIF

! smooth grid by solving Laplace equation -------------------------------------

  IF (someMoved) THEN
    resid = 0._RFREAL
    DO iter=1,global%moveGridNiter
      CALL RFLO_LaplaceGridSmoo( regions,resid )
    ENDDO
!    DO iter=1,global%moveGridNiter
!      CALL RFLO_ElliptGridSmoo( regions,resid )
!    ENDDO

    IF (global%verbLevel /= VERBOSE_NONE) THEN
#ifdef MPI
      CALL MPI_Reduce( resid,globalResid,1,MPI_RFREAL,MPI_SUM, &
                       MASTERPROC,global%mpiComm,global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
           __LINE__ )
#else
      globalResid = resid
#endif
      IF (global%myProcid == MASTERPROC .AND. &
          global%verbLevel>=VERBOSE_HIGH) THEN
        WRITE(STDOUT,2000) SOLVER_NAME,global%skewness,global%minVol
        WRITE(STDOUT,1000) SOLVER_NAME,global%moveGridNiter,SQRT(globalResid)
      ENDIF
    ENDIF    ! verbLevel
  ENDIF      ! someMoved

! update grid, dummy, corner and edge cells -----------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! --- change xyz from coordinates to deformations

      xyz    => regions(iReg)%levels(1)%grid%xyz
      xyzOld => regions(iReg)%levels(1)%gridOld%xyz

      DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
        xyz(XCOORD,ijk) = xyz(XCOORD,ijk) - xyzOld(XCOORD,ijk)
        xyz(YCOORD,ijk) = xyz(YCOORD,ijk) - xyzOld(YCOORD,ijk)
        xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) - xyzOld(ZCOORD,ijk)
      ENDDO

! --- redistribute deformations at boundaries

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld
      grid%boundMoved(:) = .true.
      grid%edgeMoved(:)  = .true.

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
!        IF ((bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE)) THEN
!          grid%boundMoved(patch%lbound) = .false.
!        ENDIF  ! bcType
        IF (bcType.eq.BC_SYMMETRY) THEN
          grid%boundMoved(patch%lbound) = .false.
        ENDIF  ! bcType
      ENDDO    ! iPatch

      CALL RFLO_BoundaryDeformation( regions(iReg),grid%boundMoved, &
                                     grid%edgeMoved,grid%arcLen12, &
                                     grid%arcLen34,grid%arcLen56, &
                                     gridOld%xyzOld,grid%xyz )

      CALL RFLO_MgElliptBndDeformation( regions(iReg) )

! --- change xyz from deformations to coordinates

      CALL RFLO_ChangeInteriorGrid( regions(iReg),grid%boundMoved, &
                                    grid%edgeMoved,grid%arcLen12, &
                                    grid%arcLen34,grid%arcLen56, &
                                    gridOld%xyzOld,grid%xyz )

! --- perform volume smoothing based on 3D Elliptic PDE

      IF (regions(iReg)%blockShape==REGION_SHAPE_NORMAL) THEN
        DO iter=1,global%moveGridViter
          CALL RFLO_ElliptGridSmooRegion( regions(iReg),resid )
        ENDDO
      ENDIF

! --- update coarse grids and dummy cells

      CALL RFLO_GenerateCoarseGrids( regions(iReg) )   ! coarsen finest grid
      CALL RFLO_CopyGeometryDummy( regions(iReg) )     ! copy to dummy nodes
      CALL RFLO_ExtrapolateGeometry( regions(iReg) )   ! extrapolate
    ENDIF     ! region on this processor and active, grid moving
  ENDDO       ! iReg

  CALL RFLO_ExchangeGeometry( regions )                ! exchange geometry

! print residual of 3D Elliptic PDE smoothing ---------------------------------

  IF (global%verbLevel /= VERBOSE_NONE) THEN
#ifdef MPI
    CALL MPI_Reduce( resid,globalResid,1,MPI_RFREAL,MPI_SUM, &
                     MASTERPROC,global%mpiComm,global%mpierr )
    IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
         __LINE__ )
#else
    globalResid = resid
#endif
    IF (global%myProcid == MASTERPROC .AND. &
          global%verbLevel>=VERBOSE_HIGH) THEN
      WRITE(STDOUT,1000) SOLVER_NAME,global%moveGridViter,SQRT(globalResid)
    ENDIF
  ENDIF    ! verbLevel

! calculate new metrics and grid speeds ---------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE .AND. &             ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN           ! and moving
      CALL RFLO_CalcFaceVectors( regions(iReg) )         ! faces
      CALL RFLO_CalcControlVolumes( regions(iReg) )      ! volumes
      CALL RFLO_CalcCellCentroids( regions(iReg) )       ! cell centroids
      IF (regions(iReg)%mixtInput%faceEdgeAvg==FE_AVG_LINEAR) &
        CALL RFLO_C2fAvgCoeffs( regions(iReg) )          ! cell2face averaging
      CALL RFLO_C2eAvgCoeffs( regions(iReg) )            ! cell2edge averaging
      CALL RFLO_CheckMetrics( iReg,regions(iReg) )       ! check metrics
      CALL RFLO_CalcGridSpeeds( regions(iReg) )          ! grid speeds
    ENDIF   ! region on this processor and active, grid moving
  ENDDO     ! iReg

! global grid quality measure -------------------------------------------------

  CALL RFLO_GridQualityGlobal( regions )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(A,1X,'Block-Elliptic-PDE grid motion:',I6,1PE13.4)
2000 FORMAT(A,1X,'global skewness, minvol:',2(1PE14.5))

END SUBROUTINE RFLO_MoveGridElliptGlo


!******************************************************************************
!
! Purpose: calculate node displacements on non-external flat patches
!          (finest grid only).
!
! Description: none.
!
! Input: region     = grid dimensions
!
! Output: xyz = updated deformations at boundaries.
!
!******************************************************************************

SUBROUTINE RFLO_MgElliptBndDeformation( region )

  USE RFLO_ModElliptSmoothing, ONLY : RFLO_ElliptGridJac2D, &
                                      RFLO_ElliptGridGauss2D, &
                                      RFLO_ElliptGridSOR2D
  USE RFLO_ModMoveGridUtil, ONLY : RFLO_MoveGridQFlatPatch, &
                                   RFLO_MoveGridCurvedPatch

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, ijk

! ... local variables
  REAL(RFREAL), POINTER   :: xyz(:,:), xyzOld(:,:)
  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_MgElliptBndDeformation',&
       'RFLO_ModMoveGridElliptGlo.F90' )

! transform deformations to coordinates ---------------------------------------

  xyz    => region%levels(1)%grid%xyz
  xyzOld => region%levels(1)%gridOld%xyz

  DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
    xyz(XCOORD,ijk) = xyz(XCOORD,ijk) + xyzOld(XCOORD,ijk)
    xyz(YCOORD,ijk) = xyz(YCOORD,ijk) + xyzOld(YCOORD,ijk)
    xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) + xyzOld(ZCOORD,ijk)
  ENDDO

! move nodes on boundaries with active edges ----------------------------------

  DO iPatch=1,region%nPatches
    patch  => region%levels(1)%patches(iPatch)

    IF (patch%bcMotion/=BC_EXTERNAL) THEN 

      IF (patch%bndFlat .EQV. .TRUE.) THEN
        CALL RFLO_ElliptGridJac2D( region,patch,iPatch )
!        CALL RFLO_ElliptGridGauss2D( region,patch,iPatch )
!        CALL RFLO_ElliptGridSOR2D( region,patch,iPatch )
      ELSEIF ((.NOT. patch%bndFlat) .AND. (patch%dirFlat > 0)) THEN
        CALL RFLO_MoveGridQFlatPatch( region,patch,iPatch )
      ELSEIF ((.NOT. patch%bndFlat) .AND. (patch%dirFlat < 0)) THEN
!need test        WRITE(STDOUT,*)'curved patch found, ireg, ipatch', &
!need test                        region%iRegionGlobal, iPatch
!need test        WRITE(STDOUT,*)'try activate RFLO_MoveGridCurvedPatch'
!need test        WRITE(STDOUT,*)'and turn off 2nd RFLO_BoundaryDeformation'
!need test        CALL RFLO_MoveGridCurvedPatch( region,patch,iPatch )
      ENDIF

    ENDIF      ! not.external
  ENDDO        ! iPatch

! transform coordinates to deformations ---------------------------------------

  DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
    xyz(XCOORD,ijk) = xyz(XCOORD,ijk) - xyzOld(XCOORD,ijk)
    xyz(YCOORD,ijk) = xyz(YCOORD,ijk) - xyzOld(YCOORD,ijk)
    xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) - xyzOld(ZCOORD,ijk)
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgElliptBndDeformation


!******************************************************************************
!
! Purpose: receive and distribute the deformations of surfaces
!          in block-wise manner.
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions%levels%grid%xyz = deformations at the boundaries
!         someMoved               = parts of grid moved.
!
! Notes: grid%xyz temporarily stores nodal displacements. The deformation
!        is applied to the finest grid first.
!
!******************************************************************************

SUBROUTINE RFLO_MgElliptSurfacesGlo( regions,someMoved )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY       : t_grid
  USE ModInterfaces, ONLY : RFLO_GetDeformation, RFLO_ArcLengthBounds, &
                          RFLO_EdgeDeformation, RFLO_EdgeDeformationStraight, &
                          RFLO_BoundaryDeformation
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  LOGICAL :: someMoved

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_grid), POINTER   :: grid, gridOld
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgElliptSurfacesGlo',&
       'RFLO_ModMoveGridElliptGlo.F90' )

! move grid separately for each region ----------------------------------------

  someMoved = .false.

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE .AND. &            ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN          ! and moving

      grid      => regions(iReg)%levels(1)%grid
      gridOld   => regions(iReg)%levels(1)%gridOld
      someMoved =  .true.

! --- store the old grid

      gridOld%indSvel  = grid%indSvel
      gridOld%ipc      = grid%ipc
      gridOld%jpc      = grid%jpc
      gridOld%kpc      = grid%kpc
      gridOld%xyz(:,:) = grid%xyz(:,:)
      gridOld%si(:,:)  = grid%si(:,:)
      gridOld%sj(:,:)  = grid%sj(:,:)
      gridOld%sk(:,:)  = grid%sk(:,:)
      gridOld%vol(:)   = grid%vol(:)

! --- calculate arclengths between boundaries

      CALL RFLO_ArcLengthBounds( regions(iReg),gridOld%xyzOld, &
                                 grid%arcLen12,grid%arcLen34,grid%arcLen56 )

! --- get the boundary deformations

      CALL RFLO_GetDeformation( regions(iReg),grid%boundMoved,grid%xyz )

! --- calculate deformations at remaining edges

      grid%edgeMoved(:) = .FALSE.
      CALL RFLO_MgElliptEdgesGlo( regions(iReg),grid%boundMoved, &
                                  grid%edgeMoved, grid%arcLen12, &
                                  grid%arcLen34, grid%arcLen56, &
                                  gridOld%xyzOld,gridOld%xyz,grid%xyz )

! --- correct deformations at straight edges

      IF (global%moveGridNiter < 1) THEN
        CALL RFLO_EdgeDeformationStraight( regions(iReg),grid%boundMoved, &
                                 grid%edgeStraight,grid%edgeMoved, &
                                 grid%arcLen12,grid%arcLen34,grid%arcLen56, &
                                 gridOld%xyzOld,gridOld%xyz,grid%xyz )
      ENDIF

! --- calculate deformations at remaining boundaries

      CALL RFLO_BoundaryDeformation( regions(iReg),grid%boundMoved, &
                                     grid%edgeMoved,grid%arcLen12, &
                                     grid%arcLen34,grid%arcLen56, &
                                     gridOld%xyzOld,grid%xyz )

    ENDIF   ! region on this processor and active, grid moving
  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgElliptSurfacesGlo


!******************************************************************************
!
! Purpose: exchange deformations between the regions as to ensure
!          matching grid nodes at the interfaces.
!
! Description: none.
!
! Input: regions = data of all grid regions, deformations.
!
! Output: regions%levels%grid%xyz = deformations at the boundaries.
!
! Notes: grid%xyz temporarily stores nodal displacements. The deformation
!        is applied to the finest grid first.
!
!******************************************************************************

SUBROUTINE RFLO_MgElliptInterfacesGlo( regions )

  USE ModInterfaces, ONLY : RFLO_ExchangeDnodeCopy, RFLO_EdgeDeformation, &
        RFLO_BoundaryDeformation, RFLO_ExchangeDnodeSend, &
        RFLO_ExchangeDnodeRecv, RFLO_ClearSendRequests
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, iPass

! ... local variables
  INTEGER :: bcType, iRegSrc, iPatchSrc

  TYPE(t_grid), POINTER   :: grid, gridOld, gridSrc
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch, patchSrc

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgElliptInterfacesGlo',&
       'RFLO_ModMoveGridElliptGlo.F90' )

! fix interfaces between regions ----------------------------------------------

  DO iPass=1,2

! - copy / send deformations

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
                                           patch,patchSrc,.false., &
                                           grid%xyz,gridSrc%xyz )
              CALL RFLO_MgElliptEdgesGlo( regions(iReg),grid%boundMoved, &
                                       grid%edgeMoved, grid%arcLen12, &
                                       grid%arcLen34, grid%arcLen56, &
                                       gridOld%xyzOld,gridOld%xyz,grid%xyz )
              CALL RFLO_BoundaryDeformation( regions(iReg),grid%boundMoved, &
                                       grid%edgeMoved,grid%arcLen12, &
                                       grid%arcLen34,grid%arcLen56, &
                                       gridOld%xyzOld,grid%xyz )
            ELSE
              CALL RFLO_ExchangeDnodeSend( regions(iReg),regions(iRegSrc), &
                                           patch,grid%xyz )
            ENDIF
          ENDIF  ! bcType
        ENDDO    ! iPatch

      ENDIF  ! region on this processor and active, grid moving
    ENDDO    ! iReg

! - receive deformations

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
                                           patch,patchSrc,.false.,grid%xyz )
              CALL RFLO_MgElliptEdgesGlo( regions(iReg),grid%boundMoved, &
                                       grid%edgeMoved, grid%arcLen12, &
                                       grid%arcLen34, grid%arcLen56, &
                                       gridOld%xyzOld,gridOld%xyz,grid%xyz )
              CALL RFLO_BoundaryDeformation( regions(iReg),grid%boundMoved, &
                                       grid%edgeMoved,grid%arcLen12, &
                                       grid%arcLen34,grid%arcLen56, &
                                       gridOld%xyzOld,grid%xyz )
            ENDIF
          ENDIF  ! bcType
        ENDDO    ! iPatch

      ENDIF  ! region on this processor and active, grid moving
    ENDDO    ! iReg

! - clear send requests

    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
          regions(iReg)%active==ACTIVE .AND. &           ! on my processor
          regions(iReg)%mixtInput%moveGrid) THEN         ! and moving
        CALL RFLO_ClearSendRequests( regions,iReg,.true. )
      ENDIF
    ENDDO

  ENDDO    ! iPass

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgElliptInterfacesGlo


!******************************************************************************
!
! Purpose: calculate node displacements on those edges whose end points have
!          moved, but the associated boundaries were not updated yet (finest
!          grid only).
!
! Description: points along an edge are shifted using 1-D linear transfinite
!              interpolation (TFI).
!
! Input: region     = grid dimensions
!        boundMoved = flag for boundaries of a region which have moved
!        arcLen12   = arclength between i=const. boundaries for each j, k
!        arcLen34   = arclength between j=const. boundaries for each k, i
!        arcLen56   = arclength between k=const. boundaries for each i, j
!        xyzOld     = grid from previous time step.
!
! Output: edgeMoved = flag if discretization at an edge was changed
!         dNode     = updated deformations at edges.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************

SUBROUTINE RFLO_MgElliptEdgesGlo( region,boundMoved,edgeMoved, &
                                  arcLen12,arcLen34,arcLen56,xyzRef,xyzOld, &
                                  dNode )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes, &
                            RFLO_Tfint1d
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6), edgeMoved(12)

  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:), xyzRef(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iEdge, ind

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, l1c, l2c
  INTEGER :: indBeg, indEnd, ijkN, ijkN1, ijkNBeg, ijkNEnd, iNOff, ijNOff
  INTEGER :: switch(12,9)

  REAL(RFREAL) :: arcLen, ds, s, dN(3), dNBeg(3), dNEnd(3)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_MgElliptEdgesGlo',&
       'RFLO_ModMoveGridElliptGlo.F90' )

! get dimensions --------------------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! set edge switch -------------------------------------------------------------
! switch(:,1) = begins at boundary
! switch(:,2) = ends on boundary
! switch(:,3) = right boundary
! switch(:,4) = left boundary
! switch(:,5) = direction (from-to boundary)
! switch(:,6) = start index
! switch(:,7) = end index
! switch(:,8) = constant index in 1st direction
! switch(:,9) = constant index in 2nd direction

  switch( 1,:) = (/5, 6, 1, 3, 56, kpnbeg, kpnend, ipnbeg, jpnbeg/)
  switch( 2,:) = (/3, 4, 1, 6, 34, jpnbeg, jpnend, kpnend, ipnbeg/)
  switch( 3,:) = (/5, 6, 1, 4, 56, kpnbeg, kpnend, ipnbeg, jpnend/)
  switch( 4,:) = (/3, 4, 1, 5, 34, jpnbeg, jpnend, kpnbeg, ipnbeg/)
  switch( 5,:) = (/5, 6, 2, 3, 56, kpnbeg, kpnend, ipnend, jpnbeg/)
  switch( 6,:) = (/3, 4, 2, 6, 34, jpnbeg, jpnend, kpnend, ipnend/)
  switch( 7,:) = (/5, 6, 2, 4, 56, kpnbeg, kpnend, ipnend, jpnend/)
  switch( 8,:) = (/3, 4, 2, 5, 34, jpnbeg, jpnend, kpnbeg, ipnend/)
  switch( 9,:) = (/1, 2, 3, 5, 12, ipnbeg, ipnend, jpnbeg, kpnbeg/)
  switch(10,:) = (/1, 2, 3, 6, 12, ipnbeg, ipnend, jpnbeg, kpnend/)
  switch(11,:) = (/1, 2, 4, 5, 12, ipnbeg, ipnend, jpnend, kpnbeg/)
  switch(12,:) = (/1, 2, 4, 6, 12, ipnbeg, ipnend, jpnend, kpnend/)

! edge movement flag ----------------------------------------------------------

  IF (boundMoved(1)) THEN
    edgeMoved( 1) = .true.; edgeMoved( 2) = .true.
    edgeMoved( 3) = .true.; edgeMoved( 4) = .true.
  ENDIF
  IF (boundMoved(2)) THEN
    edgeMoved( 5) = .true.; edgeMoved( 6) = .true.
    edgeMoved( 7) = .true.; edgeMoved( 8) = .true.
  ENDIF
  IF (boundMoved(3)) THEN
    edgeMoved( 1) = .true.; edgeMoved( 5) = .true.
    edgeMoved( 9) = .true.; edgeMoved(10) = .true.
  ENDIF
  IF (boundMoved(4)) THEN
    edgeMoved( 3) = .true.; edgeMoved( 7) = .true.
    edgeMoved(11) = .true.; edgeMoved(12) = .true.
  ENDIF
  IF (boundMoved(5)) THEN
    edgeMoved( 4) = .true.; edgeMoved( 8) = .true.
    edgeMoved( 9) = .true.; edgeMoved(11) = .true.
  ENDIF
  IF (boundMoved(6)) THEN
    edgeMoved( 2) = .true.; edgeMoved( 6) = .true.
    edgeMoved(10) = .true.; edgeMoved(12) = .true.
  ENDIF

! loop over all 12 edges ------------------------------------------------------

  DO iEdge=1,12
    IF ((boundMoved(switch(iEdge,1)) .OR. boundMoved(switch(iEdge,2))) .AND. &
        ((.NOT.boundMoved(switch(iEdge,3))) .OR. &
         (.NOT.boundMoved(switch(iEdge,4)))) .AND. &
         (.NOT.edgeMoved(iEdge))) THEN

      edgeMoved(iEdge) = .true.

      ds     = 0._RFREAL
      indBeg = switch(iEdge,6)
      indEnd = switch(iEdge,7)
      l1c    = switch(iEdge,8)
      l2c    = switch(iEdge,9)
      DO ind=indBeg+1,indEnd-1
        IF (switch(iEdge,5) == 12) THEN
          ijkN     = IndIJK(ind   ,l1c,l2c,iNOff,ijNOff)
          ijkN1    = IndIJK(ind-1 ,l1c,l2c,iNOff,ijNOff)
          ijkNBeg  = IndIJK(indBeg,l1c,l2c,iNOff,ijNOff)
          ijkNEnd  = IndIJK(indEnd,l1c,l2c,iNOff,ijNOff)
          arcLen   = arcLen12(l1c,l2c)
        ELSE IF (switch(iEdge,5) == 34) THEN
          ijkN     = IndIJK(l2c,ind   ,l1c,iNOff,ijNOff)
          ijkN1    = IndIJK(l2c,ind-1 ,l1c,iNOff,ijNOff)
          ijkNBeg  = IndIJK(l2c,indBeg,l1c,iNOff,ijNOff)
          ijkNEnd  = IndIJK(l2c,indEnd,l1c,iNOff,ijNOff)
          arcLen   = arcLen34(l1c,l2c)
        ELSE IF (switch(iEdge,5) == 56) THEN
          ijkN     = IndIJK(l1c,l2c,ind   ,iNOff,ijNOff)
          ijkN1    = IndIJK(l1c,l2c,ind-1 ,iNOff,ijNOff)
          ijkNBeg  = IndIJK(l1c,l2c,indBeg,iNOff,ijNOff)
          ijkNEnd  = IndIJK(l1c,l2c,indEnd,iNOff,ijNOff)
          arcLen   = arcLen56(l1c,l2c)
        ENDIF
        dNBeg(:) = dNode(:,ijkNBeg) + xyzOld(:,ijkNBeg) - xyzRef(:,ijkNBeg)
        dNEnd(:) = dNode(:,ijkNEnd) + xyzOld(:,ijkNEnd) - xyzRef(:,ijkNEnd)

        ds = ds + SQRT((xyzRef(XCOORD,ijkN)-xyzRef(XCOORD,ijkN1))**2 + &
                       (xyzRef(YCOORD,ijkN)-xyzRef(YCOORD,ijkN1))**2 + &
                       (xyzRef(ZCOORD,ijkN)-xyzRef(ZCOORD,ijkN1))**2)
        s  = ds/arcLen

        CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
        dNode(:,ijkN) = dN(:) + xyzRef(:,ijkN) - xyzOld(:,ijkN) 
      ENDDO   ! i

    ENDIF     ! boundMoved
  ENDDO       ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_MgElliptEdgesGlo


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModMoveGridElliptGlo


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModMoveGridElliptGlo.F90,v $
! Revision 1.18  2009/08/27 14:04:50  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.17  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.16  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.15  2006/05/19 03:14:49  wasistho
! commented MovedGridCurvedPatch call
!
! Revision 1.14  2006/03/18 11:49:28  wasistho
! called GridQualityGlobal
!
! Revision 1.13  2006/03/18 11:03:45  wasistho
! screen printed global skewness and minvol
!
! Revision 1.12  2006/03/18 09:25:01  wasistho
! modified mgElliptEdgesGlo
!
! Revision 1.11  2006/03/18 08:18:12  wasistho
! added mgElliptSurfacesGlo
!
! Revision 1.10  2006/03/17 06:38:17  wasistho
! added call to non-flat-surface smoother
!
! Revision 1.9  2006/03/16 08:26:49  wasistho
! added interfacesGlo and edgesGlo
!
! Revision 1.8  2006/03/08 06:41:01  wasistho
! added moveGridViter
!
! Revision 1.7  2006/03/05 19:08:24  wasistho
! set computational space coordinates from initial grid
!
! Revision 1.6  2006/03/04 04:35:36  wasistho
! call elliptGridSmooRegion if blockShape is normal
!
! Revision 1.5  2006/03/03 05:23:09  wasistho
! made public MgElliptBndDeformation
!
! Revision 1.4  2006/03/03 04:14:53  wasistho
! activated elliptGridSmooRegion
!
! Revision 1.3  2006/03/02 00:23:44  wasistho
! prepared elliptic pde grid motion
!
! Revision 1.2  2006/02/23 21:35:17  wasistho
! initialzed resid
!
! Revision 1.1  2006/02/11 03:42:14  wasistho
! added ModMoveGridElliptGlo/Fra
!
!
!
! ******************************************************************************











