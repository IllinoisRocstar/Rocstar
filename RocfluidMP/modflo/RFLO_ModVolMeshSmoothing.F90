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
! Purpose: Suite of volume mesh smoothing (VMS) routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModVolMeshSmoothing.F90,v 1.8 2009/08/27 14:04:51 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModVolMeshSmoothing

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
  PUBLIC :: RFLO_MoveGridVms

! private : RFLO_VmsInit
!           RFLO_VmsAverageVertices
!           RFLO_VmsLaplaceIterate
!           RFLO_VmsLaplacePerturb
!           RFLO_VmsLaplaceProcedure
!           RFLO_VmsProjectVertices
!           RFLO_VmsRestoreBoundDeform
 
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModVolMeshSmoothing.F90,v $ $Revision: 1.8 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: redistribute grid nodes according to the movement of the
!          boundaries. This function smoothes the grid globally by
!          volume mesh smoothing based on Laplacian propagation.
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

SUBROUTINE RFLO_MoveGridVms( regions )

  USE ModInterfaces, ONLY : RFLO_MoveGridSurfaces, RFLO_MoveGridInterfaces, &
        RFLO_ChangeInteriorGrid, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_CalcControlVolumes, &
        RFLO_CalcCellCentroids, RFLO_C2fAvgCoeffs, RFLO_C2eAvgCoeffs, &
        RFLO_CheckMetrics, RFLO_CalcGridSpeeds, RFLO_BoundaryDeformation
  USE RFLO_ModLaplaceSmoothing, ONLY : RFLO_LaplaceGridSmoo

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

  CALL RegisterFunction( global,'RFLO_MoveGridVms',&
       'RFLO_ModVolMeshSmoothing.F90' )

#ifdef GENX
! update geometry buffers -----------------------------------------------------

  dAlpha = global%dtMin/global%dTimeSystem
  CALL COM_call_function( global%genxHandleGm,1,dAlpha )
#endif

! receive and distribute deformations for each region -------------------------

  CALL RFLO_VmsInit( regions,someMoved )

! smooth grid by solving Laplace equation -------------------------------------

  IF (someMoved) THEN
    DO iter=1,global%moveGridNiter
      CALL RFLO_VmsLaplaceIterate( regions,iter,resid )
    ENDDO

    IF (global%verbLevel /= VERBOSE_NONE) THEN
#ifdef MPI
      CALL MPI_Reduce( resid,globalResid,1,MPI_RFREAL,MPI_SUM, &
                       MASTERPROC,global%mpiComm,global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
           __LINE__ )
#else
      globalResid = resid
#endif
      IF (global%myProcid == MASTERPROC) THEN
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
 !         grid%boundMoved(patch%lbound) = .false.
 !       ENDIF  ! bcType
        IF ((bcType.eq.BC_SYMMETRY)) THEN
          grid%boundMoved(patch%lbound) = .false.
        ENDIF  ! bcType
      ENDDO    ! iPatch
      CALL RFLO_BoundaryDeformation( regions(iReg),grid%boundMoved, &
                                     grid%edgeMoved,grid%arcLen12, &
                                     grid%arcLen34,grid%arcLen56, &
                                     gridOld%xyzOld,grid%xyz )

! --- change xyz from deformations to coordinates

      CALL RFLO_ChangeInteriorGrid( regions(iReg),grid%boundMoved, &
                                    grid%edgeMoved,grid%arcLen12, &
                                    grid%arcLen34,grid%arcLen56, &
                                    gridOld%xyzOld,grid%xyz )

! --- update coarse grids and dummy cells

      CALL RFLO_GenerateCoarseGrids( regions(iReg) )   ! coarsen finest grid
      CALL RFLO_CopyGeometryDummy( regions(iReg) )     ! copy to dummy nodes
      CALL RFLO_ExtrapolateGeometry( regions(iReg) )   ! extrapolate
    ENDIF     ! region on this processor and active, grid moving
  ENDDO       ! iReg

  CALL RFLO_ExchangeGeometry( regions )                ! exchange geometry

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

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(A,1X,'VMS grid motion: ',I6,1PE13.4)

END SUBROUTINE RFLO_MoveGridVms

!******************************************************************************
!
! Purpose: initial procedure for volume mesh smoothing.
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: someMoved = parts of grid moved.
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_VmsInit( regions,someMoved )

  USE ModInterfaces, ONLY : RFLO_GetDeformation
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

  CALL RegisterFunction( global,'RFLO_VmsInit',&
       'RFLO_ModVolMeshSmoothing.F90' )

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

! --- get the boundary deformations

      CALL RFLO_GetDeformation( regions(iReg),grid%boundMoved,grid%xyz )

    ENDIF   ! region on this processor and active, grid moving
  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_VmsInit

!******************************************************************************
!
! Purpose: smooth the distribution of grid points by iterating simplified
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

SUBROUTINE RFLO_VmsLaplaceIterate( regions,iter,resid )

  USE ModInterfaces, ONLY : RFLO_ChangeInteriorGrid, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_ExchangeDnodeCopy, &
        RFLO_ExchangeDnodeSend, RFLO_ExchangeDnodeRecv, RFLO_ClearSendRequests,&
        RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE RFLO_ModLaplaceSmoothing, ONLY : RFLO_LaplaceGridPatch

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iter
  REAL(RFREAL) :: resid

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, ijk, i, j, k

! ... local variables
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff
  INTEGER :: bcType, iRegSrc, iPatchSrc

  REAL(RFREAL) :: dx, dy, dz
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER   :: grid, gridOld, gridSrc
  TYPE(t_patch), POINTER  :: patch, patchSrc

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_VmsLaplaceIterate',&
       'RFLO_ModVolMeshSmoothing.F90' )

! smooth grid region-wise -----------------------------------------------------

  resid = 0._RFREAL

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! --- compute movements in the interior and along the boundaries

      CALL RFLO_VmsLaplaceProcedure( regions(iReg),iter )

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
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
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

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_VmsLaplaceIterate

!******************************************************************************
!
! Purpose: average vertices.
!
! Description: defined grid coordinates by averaging 12 projected vertices,
!              4 in each face-direction.
!
! Input: region = data of current region.
!
! Output: region%levels%grid%xyz = new coordinates.
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_VmsAverageVertices( region,vxyz )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  REAL(RFREAL), POINTER :: vxyz(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: iLev, ipnbeg,ipnend,jpnbeg,jpnend,kpnbeg,kpnend, iNOff,ijNOff, ijk
  REAL(RFREAL) :: rd
  REAL(RFREAL), POINTER :: xyz(:,:)

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_VmsAverageVertices',&
       'RFLO_ModVolMeshSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz => region%levels(iLev)%grid%xyz

! compute new coordinates -----------------------------------------------------

  rd  = 1._RFREAL/12._RFREAL

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)

        xyz(XCOORD,ijk) = vxyz(XCOORD,ijk)*rd
        xyz(YCOORD,ijk) = vxyz(YCOORD,ijk)*rd
        xyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)*rd
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_VmsAverageVertices

!******************************************************************************
!
! Purpose: conduct one Jacobi iteration to obtain initial grid movements
!          in the flow domain (boundaries included). This initial movements
!          function as perturbation.
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates.
!
! Output: region%levels%grid%xyz = grid movements.
!
! Notes: on entry, xyz holds node coordinates from a previous smoothing
!        step. On exit xyz contains the perturbed coordinates. 
!
!******************************************************************************

SUBROUTINE RFLO_VmsLaplacePerturb( region )

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
  INTEGER :: ijk, ip1, im1, jp1, jm1, kp1, km1

  REAL(RFREAL) :: rxi, ryi, rzi, rxj, ryj, rzj, rxk, ryk, rzk, rd
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:), xyzOrig(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_VmsLaplacePerturb',&
       'RFLO_ModVolMeshSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  xyz     => region%levels(iLev)%grid%xyz
  xyzOld  => region%levels(iLev)%grid%xyzOld
  xyzOrig => region%levels(iLev)%gridOld%xyz

  rd = 1._RFREAL/6._RFREAL

! reset motion vectors --------------------------------------------------------

  DO ijk=ibn,ien
    xyzOld(XCOORD,ijk) = xyz(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
    xyzOld(YCOORD,ijk) = xyz(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
    xyzOld(ZCOORD,ijk) = xyz(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
  ENDDO

! compute new coordinates -----------------------------------------------------

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ip1  = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        im1  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        jp1  = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        jm1  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        kp1  = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        km1  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

        rxi = xyzOld(XCOORD,im1) + xyzOld(XCOORD,ip1)
        ryi = xyzOld(YCOORD,im1) + xyzOld(YCOORD,ip1)
        rzi = xyzOld(ZCOORD,im1) + xyzOld(ZCOORD,ip1)

        rxj = xyzOld(XCOORD,jm1) + xyzOld(XCOORD,jp1)
        ryj = xyzOld(YCOORD,jm1) + xyzOld(YCOORD,jp1)
        rzj = xyzOld(ZCOORD,jm1) + xyzOld(ZCOORD,jp1)

        rxk = xyzOld(XCOORD,km1) + xyzOld(XCOORD,kp1)
        ryk = xyzOld(YCOORD,km1) + xyzOld(YCOORD,kp1)
        rzk = xyzOld(ZCOORD,km1) + xyzOld(ZCOORD,kp1)

        xyz(XCOORD,ijk) = xyz(XCOORD,ijk) + (rxi+rxj+rxk)*rd
        xyz(YCOORD,ijk) = xyz(YCOORD,ijk) + (ryi+ryj+ryk)*rd
        xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) + (rzi+rzj+rzk)*rd
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_VmsLaplacePerturb

!******************************************************************************
!
! Purpose: step by step procedure for Laplacian volume mesh smoothing.
!
! Description: none.
!
! Input: region = data of current grid region.
!
! Output: region%levels%grid%xyz = grid movements.
!
! Notes: on entry, xyz holds node coordinates from a previous smoothing
!        step. On exit however, xyz contains only the grid motion. 
!
!******************************************************************************

SUBROUTINE RFLO_VmsLaplaceProcedure( region,iter )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, &
                            RFLO_CalcCellCentroids        
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER :: iter

! ... loop variables
  INTEGER :: in

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: ibn,ien, idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER :: iNOff, ijNOff, errFl
  REAL(RFREAL), POINTER :: dxyz(:,:), vxyz(:,:)

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_VmsLaplaceProcedure',&
       'RFLO_ModVolMeshSmoothing.F90' )

! allocate buffers ------------------------------------------------------------

  CALL RFLO_GetDimensDummyNodes( region,1,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,1,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  ALLOCATE( dxyz(XCOORD:ZCOORD,ibn:ien), stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( vxyz(XCOORD:ZCOORD,ibn:ien), stat=errFl ); IF (errFl>0) GOTO 88

! restore grid coordinate, dxyz=original-deformation

  IF (iter==1) THEN
    dxyz = region%levels(1)%grid%xyz
    DO in = ibn,ien
      region%levels(1)%grid%xyz(XCOORD,in)= &
             region%levels(1)%gridOld%xyz(XCOORD,in) + dxyz(XCOORD,in)
      region%levels(1)%grid%xyz(YCOORD,in)= &
             region%levels(1)%gridOld%xyz(YCOORD,in) + dxyz(YCOORD,in)
      region%levels(1)%grid%xyz(ZCOORD,in)= &
             region%levels(1)%gridOld%xyz(ZCOORD,in) + dxyz(ZCOORD,in)
    ENDDO
  ENDIF

! initial cell centroid, xyz=actual-grid --------------------------------------

  CALL RFLO_CalcCellCentroids( region )

! perturb grid, xyz=actual-grid -----------------------------------------------

  CALL RFLO_VmsLaplacePerturb( region )

! restore interacting boundary grid, xyz=actual grid --------------------------

  CALL RFLO_VmsRestoreBoundDeform( region,ibn,ien,dxyz )

! face vector of perturbed grid, xyz=actual-grid ------------------------------

  CALL RFLO_CalcFaceVectors( region )
  
! project grid vertices to new planes, xyz=actual-grid ------------------------

  CALL RFLO_VmsProjectVertices( region,vxyz )

! restore interacting boundary grid, xyz=actual grid --------------------------

  CALL RFLO_VmsRestoreBoundDeform( region,ibn,ien,dxyz )

! finally average vertices, xyz=actual-grid -----------------------------------

  CALL RFLO_VmsAverageVertices( region,vxyz )

! transform from coordinate to deformation ------------------------------------

  DO in = ibn,ien
    region%levels(1)%grid%xyz(XCOORD,in)= &
           region%levels(1)%grid%xyz(XCOORD,in)- &
           region%levels(1)%gridOld%xyz(XCOORD,in) 
    region%levels(1)%grid%xyz(YCOORD,in)= &
           region%levels(1)%grid%xyz(YCOORD,in)- &
           region%levels(1)%gridOld%xyz(YCOORD,in) 
    region%levels(1)%grid%xyz(ZCOORD,in)= &
           region%levels(1)%grid%xyz(ZCOORD,in)- &
           region%levels(1)%gridOld%xyz(ZCOORD,in) 
  ENDDO

! deallocate temporary arrays -------------------------------------------------

  DEALLOCATE( dxyz, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( vxyz, stat=errFl ); IF (errFl>0) GOTO 99

  GOTO 999

! finalize --------------------------------------------------------------------

88   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

99   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_VmsLaplaceProcedure

!******************************************************************************
!
! Purpose: projects grid vertices to cell face plane of perturbed grid.
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates.
!
! Output: vxyz = sum of 12 projected cell vertices, 4 in each direction
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_VmsProjectVertices( region,vxyz )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetDimensPhysNodes, &
                            RFLO_GetNodeOffset, RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  REAL(RFREAL), POINTER :: vxyz(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: iLev, iNOff, ijNOff, iCOff, ijCOff
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, ibn, ien
  INTEGER :: ijk, m1, m2, m12, ic, errFl

  REAL(RFREAL) :: sn, c1, c2, c3, c4, wc, wm
  REAL(RFREAL), POINTER :: xyz(:,:), si(:,:), sj(:,:), sk(:,:), cofg(:,:)
  REAL(RFREAL), ALLOCATABLE :: fc(:,:), snx(:), sny(:), snz(:)

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_VmsProjectVertices',&
       'RFLO_ModVolMeshSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  xyz  => region%levels(iLev)%grid%xyz
  cofg => region%levels(iLev)%grid%cofg

! set pointers ----------------------------------------------------------------

  si => region%levels(iLev)%grid%si
  sj => region%levels(iLev)%grid%sj
  sk => region%levels(iLev)%grid%sk

  ALLOCATE( fc(XCOORD:ZCOORD,ibn:ien), stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( snx(ibn:ien),              stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( sny(ibn:ien),              stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( snz(ibn:ien),              stat=errFl ); IF (errFl>0) GOTO 88

! accumulate projected coordinates --------------------------------------------
! i-faces

  fc = 0._RFREAL

  DO k=kdnbeg+1,kdnend-1
    DO j=jdnbeg+1,jdnend-1
      DO i=idnbeg+1,idnend-1
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ic  = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        m1  = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)

        wm = 0.5_RFREAL
        wc = 0.5_RFREAL

        fc(XCOORD,ijk) = wc*cofg(XCOORD,ic) + wm*cofg(XCOORD,m1)
        fc(YCOORD,ijk) = wc*cofg(YCOORD,ic) + wm*cofg(YCOORD,m1)
        fc(ZCOORD,ijk) = wc*cofg(ZCOORD,ic) + wm*cofg(ZCOORD,m1)

        sn  = SQRT( si(XCOORD,ijk)**2+si(YCOORD,ijk)**2+si(ZCOORD,ijk)**2 )
        snx(ijk) = si(XCOORD,ijk)/(sn+1.E-12_RFREAL)
        sny(ijk) = si(YCOORD,ijk)/(sn+1.E-12_RFREAL)
        snz(ijk) = si(ZCOORD,ijk)/(sn+1.E-12_RFREAL)

      ENDDO
    ENDDO
  ENDDO

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        m1  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        m2  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        m12 = IndIJK(i  ,j-1,k-1,iNOff,ijNOff)

        c1 = snx(ijk)*(xyz(XCOORD,ijk)-fc(XCOORD,ijk))+ &
             sny(ijk)*(xyz(YCOORD,ijk)-fc(YCOORD,ijk))+ &
             snz(ijk)*(xyz(ZCOORD,ijk)-fc(ZCOORD,ijk))

        c2 = snx(m1)*(xyz(XCOORD,ijk)-fc(XCOORD,m1))+ &
             sny(m1)*(xyz(YCOORD,ijk)-fc(YCOORD,m1))+ &
             snz(m1)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m1))

        c3 = snx(m2)*(xyz(XCOORD,ijk)-fc(XCOORD,m2))+ &
             sny(m2)*(xyz(YCOORD,ijk)-fc(YCOORD,m2))+ &
             snz(m2)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m2))

        c4 = snx(m12)*(xyz(XCOORD,ijk)-fc(XCOORD,m12))+ &
             sny(m12)*(xyz(YCOORD,ijk)-fc(YCOORD,m12))+ &
             snz(m12)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m12))

        vxyz(XCOORD,ijk) = xyz(XCOORD,ijk)-c1*snx(ijk)
        vxyz(YCOORD,ijk) = xyz(YCOORD,ijk)-c1*sny(ijk)
        vxyz(ZCOORD,ijk) = xyz(ZCOORD,ijk)-c1*snz(ijk)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m1)-c2*snx(m1)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m1)-c2*sny(m1)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m1)-c2*snz(m1)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m2)-c3*snx(m2)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m2)-c3*sny(m2)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m2)-c3*snz(m2)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m12)-c4*snx(m12)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m12)-c4*sny(m12)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m12)-c4*snz(m12)

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! j-faces

  fc = 0._RFREAL

  DO k=kdnbeg,kdnend-1
    DO j=jdnbeg,jdnend-1
      DO i=idnbeg,idnend-1
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ic  = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        m1  = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)

        wm = 0.5_RFREAL
        wc = 0.5_RFREAL

        fc(XCOORD,ijk) = wc*cofg(XCOORD,ic) + wm*cofg(XCOORD,m1)
        fc(YCOORD,ijk) = wc*cofg(YCOORD,ic) + wm*cofg(YCOORD,m1)
        fc(ZCOORD,ijk) = wc*cofg(ZCOORD,ic) + wm*cofg(ZCOORD,m1)

        sn  = SQRT( sj(XCOORD,ijk)**2+sj(YCOORD,ijk)**2+sj(ZCOORD,ijk)**2 )
        snx(ijk) = sj(XCOORD,ijk)/(sn+1.E-12_RFREAL)
        sny(ijk) = sj(YCOORD,ijk)/(sn+1.E-12_RFREAL)
        snz(ijk) = sj(ZCOORD,ijk)/(sn+1.E-12_RFREAL)

      ENDDO
    ENDDO
  ENDDO

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        m1  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        m2  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        m12 = IndIJK(i-1,j  ,k-1,iNOff,ijNOff)

        c1 = snx(ijk)*(xyz(XCOORD,ijk)-fc(XCOORD,ijk))+ &
             sny(ijk)*(xyz(YCOORD,ijk)-fc(YCOORD,ijk))+ &
             snz(ijk)*(xyz(ZCOORD,ijk)-fc(ZCOORD,ijk))

        c2 = snx(m1)*(xyz(XCOORD,ijk)-fc(XCOORD,m1))+ &
             sny(m1)*(xyz(YCOORD,ijk)-fc(YCOORD,m1))+ &
             snz(m1)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m1))

        c3 = snx(m2)*(xyz(XCOORD,ijk)-fc(XCOORD,m2))+ &
             sny(m2)*(xyz(YCOORD,ijk)-fc(YCOORD,m2))+ &
             snz(m2)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m2))

        c4 = snx(m12)*(xyz(XCOORD,ijk)-fc(XCOORD,m12))+ &
             sny(m12)*(xyz(YCOORD,ijk)-fc(YCOORD,m12))+ &
             snz(m12)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m12))

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,ijk)-c1*snx(ijk)
        vxyz(YCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(YCOORD,ijk)-c1*sny(ijk)
        vxyz(ZCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(ZCOORD,ijk)-c1*snz(ijk)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m1)-c2*snx(m1)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m1)-c2*sny(m1)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m1)-c2*snz(m1)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m2)-c3*snx(m2)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m2)-c3*sny(m2)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m2)-c3*snz(m2)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m12)-c4*snx(m12)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m12)-c4*sny(m12)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m12)-c4*snz(m12)

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! k-faces

  fc = 0._RFREAL

  DO k=kdnbeg,kdnend-1
    DO j=jdnbeg,jdnend-1
      DO i=idnbeg,idnend-1
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ic  = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        m1  = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)

        wm = 0.5_RFREAL
        wc = 0.5_RFREAL

        fc(XCOORD,ijk) = wc*cofg(XCOORD,ic) + wm*cofg(XCOORD,m1)
        fc(YCOORD,ijk) = wc*cofg(YCOORD,ic) + wm*cofg(YCOORD,m1)
        fc(ZCOORD,ijk) = wc*cofg(ZCOORD,ic) + wm*cofg(ZCOORD,m1)

        sn  = SQRT( sk(XCOORD,ijk)**2+sk(YCOORD,ijk)**2+sk(ZCOORD,ijk)**2 )
        snx(ijk) = sk(XCOORD,ijk)/(sn+1.E-12_RFREAL)
        sny(ijk) = sk(YCOORD,ijk)/(sn+1.E-12_RFREAL)
        snz(ijk) = sk(ZCOORD,ijk)/(sn+1.E-12_RFREAL)

      ENDDO
    ENDDO
  ENDDO

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        m1  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        m2  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        m12 = IndIJK(i-1,j-1,k  ,iNOff,ijNOff)

        c1 = snx(ijk)*(xyz(XCOORD,ijk)-fc(XCOORD,ijk))+ &
             sny(ijk)*(xyz(YCOORD,ijk)-fc(YCOORD,ijk))+ &
             snz(ijk)*(xyz(ZCOORD,ijk)-fc(ZCOORD,ijk))

        c2 = snx(m1)*(xyz(XCOORD,ijk)-fc(XCOORD,m1))+ &
             sny(m1)*(xyz(YCOORD,ijk)-fc(YCOORD,m1))+ &
             snz(m1)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m1))

        c3 = snx(m2)*(xyz(XCOORD,ijk)-fc(XCOORD,m2))+ &
             sny(m2)*(xyz(YCOORD,ijk)-fc(YCOORD,m2))+ &
             snz(m2)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m2))

        c4 = snx(m12)*(xyz(XCOORD,ijk)-fc(XCOORD,m12))+ &
             sny(m12)*(xyz(YCOORD,ijk)-fc(YCOORD,m12))+ &
             snz(m12)*(xyz(ZCOORD,ijk)-fc(ZCOORD,m12))

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,ijk)-c1*snx(ijk)
        vxyz(YCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(YCOORD,ijk)-c1*sny(ijk)
        vxyz(ZCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(ZCOORD,ijk)-c1*snz(ijk)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m1)-c2*snx(m1)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m1)-c2*sny(m1)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m1)-c2*snz(m1)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m2)-c3*snx(m2)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m2)-c3*sny(m2)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m2)-c3*snz(m2)

        vxyz(XCOORD,ijk) = vxyz(XCOORD,ijk)+xyz(XCOORD,m12)-c4*snx(m12)
        vxyz(YCOORD,ijk) = vxyz(YCOORD,ijk)+xyz(YCOORD,m12)-c4*sny(m12)
        vxyz(ZCOORD,ijk) = vxyz(ZCOORD,ijk)+xyz(ZCOORD,m12)-c4*snz(m12)

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! deallocate temporary arrays -------------------------------------------------

  DEALLOCATE( fc , stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( snx, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( sny, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( snz, stat=errFl ); IF (errFl>0) GOTO 99

  GOTO 999

! finalize --------------------------------------------------------------------

88   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

99   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_VmsProjectVertices

!******************************************************************************
!
! Purpose: restore boundary deformation at interacting boundaries.
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates.
!        dxyz   = original movement
!
! Output: region%levels%grid%xyz = actual grid.
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_VmsRestoreBoundDeform( region,ibn,ien,dxyz )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER :: ibn, ien
  REAL(RFREAL), POINTER :: dxyz(:,:)

! ... loop variables
  INTEGER :: iPatch, i, j, k

! ... local variables
  TYPE(t_patch), POINTER :: patch

  INTEGER :: iLev, ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, ijkNB
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOrig(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_VmsRestoreBoundDeform',&
       'RFLO_ModVolMeshSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz     => region%levels(iLev)%grid%xyz
  xyzOrig => region%levels(iLev)%gridOld%xyz

! restore original boundary movements

  DO iPatch=1,region%nPatches
    patch => region%levels(iLev)%patches(iPatch)

    IF (patch%bcMotion == BC_EXTERNAL) THEN
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkNB  = IndIJK(i,j,k,iNOff,ijNOff)
            xyz(XCOORD,ijkNB) = xyzOrig(XCOORD,ijkNB)+ dxyz(XCOORD,ijkNB)
            xyz(YCOORD,ijkNB) = xyzOrig(YCOORD,ijkNB)+ dxyz(YCOORD,ijkNB)
            xyz(ZCOORD,ijkNB) = xyzOrig(ZCOORD,ijkNB)+ dxyz(ZCOORD,ijkNB) 
          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k
    ENDIF       ! bcCoupled
  ENDDO         ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_VmsRestoreBoundDeform

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLO_ModVolMeshSmoothing

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModVolMeshSmoothing.F90,v $
! Revision 1.8  2009/08/27 14:04:51  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.7  2008/12/06 08:44:17  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:28  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/03/05 21:52:57  wasistho
! changed computational space coordinates to be based on initial grid
!
! Revision 1.4  2005/10/27 05:59:11  wasistho
! added USE RFLO_ModLaplaceSmoothin
!
! Revision 1.3  2005/06/13 21:47:44  wasistho
! changed patch%bcCoupled to patch%bcMotion
!
! Revision 1.2  2005/05/28 06:11:51  wasistho
! cosmetics
!
! Revision 1.1  2005/05/21 00:16:46  wasistho
! added RFLO_ModVolMeshSmoothing
!
!
!
! ******************************************************************************














