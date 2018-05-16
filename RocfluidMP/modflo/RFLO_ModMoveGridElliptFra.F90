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
! Purpose: Suite for frame grid-motion routines with elliptic PDE smoothings.
!
! Description: The procedure is similar to MOVEGRID_FRAME method, except
!              the Laplacian smoothing is replaced by elliptic PDE smoothing.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModMoveGridElliptFra.F90,v 1.16 2009/08/27 14:04:50 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModMoveGridElliptFra

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
  PUBLIC :: RFLO_MoveGridElliptFra

! private : RFLO_MgElliptSurfacesFra
!           RFLO_MgElliptInterfacesFra
!           RFLO_MgElliptEdgesFra
!           RFLO_MgElliptBoundaries
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModMoveGridElliptFra.F90,v $ $Revision: 1.16 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


!******************************************************************************
!
! Purpose: redistribute grid nodes according to the movement of the
!          boundaries. This function smoothes the grid globally by
!          volume and surface mesh smoothing based on Elliptic PDE.
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

SUBROUTINE RFLO_MoveGridElliptFra( regions )

  USE ModInterfaces, ONLY : RFLO_MoveGridInterfaces, &
        RFLO_ChangeInteriorGrid, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_CalcControlVolumes, &
        RFLO_CalcCellCentroids, RFLO_CalcFaceCentroids, RFLO_C2fAvgCoeffs, &
        RFLO_C2eAvgCoeffs, RFLO_CheckMetrics, RFLO_CalcGridSpeeds, &
        RFLO_BoundaryDeformation
  USE RFLO_ModGridMetrics,       ONLY: RFLO_GridQualityGlobal
  USE RFLO_ModLaplaceSmoothing,  ONLY: RFLO_LaplaceGridSmoo
  USE RFLO_ModElliptSmoothing,   ONLY: RFLO_ElliptGridSmoo, &
                                       RFLO_ElliptGridSmooRegion
  USE RFLO_ModMoveGridElliptGlo, ONLY: RFLO_MgElliptBndDeformation 

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

  CALL RegisterFunction( global,'RFLO_MoveGridElliptFra',&
       'RFLO_ModMoveGridElliptFra.F90' )

#ifdef GENX
! update geometry buffers -----------------------------------------------------

  dAlpha = global%dtMin/global%dTimeSystem
  CALL COM_call_function( global%genxHandleGm,1,dAlpha )
#endif

! receive and distribute deformations for each region -------------------------

  CALL RFLO_MgElliptSurfacesFra( regions,someMoved )

! fix interfaces between regions ----------------------------------------------

  IF (someMoved) THEN
    CALL RFLO_MgElliptInterfacesFra( regions )
  ENDIF

! update grid, dummy, corner and edge cells -----------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! --- change the interior grid

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld
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

! smooth grid by Laplacian method ---------------------------------------------

  IF (global%moveGridNiter < 1 .AND. &
      global%moveGridViter < 1 .AND. &
      global%moveGridSiter < 1) THEN                   ! TFI only
    IF (global%verbLevel >= VERBOSE_HIGH) THEN
      IF (global%myProcid == MASTERPROC) THEN
        WRITE(STDOUT,4000) SOLVER_NAME,global%skewness,global%minVol
        WRITE(STDOUT,1000) SOLVER_NAME, &
                  global%moveGridNiter, global%moveGridAmplifX, &
                  global%moveGridAmplifY, global%moveGridAmplifZ, &
                  global%moveGridPower, global%moveGridOrthDir, &
                  global%moveGridOrthWghtX, global%moveGridOrthWghtY, &
                  global%moveGridOrthWghtZ
      ENDIF  ! masterproc
    ENDIF    ! verbLevel
    GOTO 888
  ENDIF      ! niter<1

  IF (someMoved) THEN
    resid = 0._RFREAL
    DO iter=1,global%moveGridNiter
      CALL RFLO_LaplaceGridSmoo( regions,resid )
    ENDDO
!    DO iter=1,global%moveGridNiter
!      CALL RFLO_ElliptGridSmoo( regions,resid )
!    ENDDO

    IF (global%verbLevel >= VERBOSE_HIGH) THEN
#ifdef MPI
      CALL MPI_Reduce( resid,globalResid,1,MPI_RFREAL,MPI_SUM, &
                       MASTERPROC,global%mpiComm,global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
           __LINE__ )
#else
      globalResid = resid
#endif
      IF (global%myProcid == MASTERPROC) THEN
        WRITE(STDOUT,4000) SOLVER_NAME,global%skewness,global%minVol
        WRITE(STDOUT,2000) SOLVER_NAME, &
                global%moveGridNiter, global%moveGridAmplifX, &
                global%moveGridAmplifY,global%moveGridAmplifZ, &
                global%moveGridPower,global%moveGridOrthDir, &
                global%moveGridOrthWghtX, global%moveGridOrthWghtY, &
                global%moveGridOrthWghtZ, SQRT(globalResid)
      ENDIF  ! masterproc
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
!   TEMPORARY - DEBUGGING OF SURFACE MOTION
      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType

        IF (bcType .eq. BC_SYMMETRY) THEN
          grid%boundMoved(patch%lbound) = .false.
        ENDIF  ! bcType
      ENDDO    ! iPatch
!   TEMPORARY - DEBUGGING OF SURFACE MOTION

!      WRITE(*,*) ' in here '
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

  IF (global%verbLevel >= VERBOSE_HIGH) THEN
#ifdef MPI
    CALL MPI_Reduce( resid,globalResid,1,MPI_RFREAL,MPI_SUM, &
                     MASTERPROC,global%mpiComm,global%mpierr )
    IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
         __LINE__ )
#else
    globalResid = resid
#endif
    IF (global%myProcid == MASTERPROC) THEN
      WRITE(STDOUT,3000) SOLVER_NAME,global%moveGridViter,SQRT(globalResid)
    ENDIF
  ENDIF    ! verbLevel

888 CONTINUE

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

1000 FORMAT(A,1X,'Global-Elliptic-PDE grid motion:', &
            I5,1X,4(1PE9.2),I4,3(1PE9.2))
2000 FORMAT(A,1X,'Global-Elliptic-PDE grid motion:', &
            I5,1X,4(1PE9.2),I4,3(1PE9.2),1PE12.4)
3000 FORMAT(A,1X,'Global-Elliptic-PDE grid motion:', &
            I5,1PE12.4)
4000 FORMAT(A,1X,'global skewness, minvol:',2(1PE14.5))

END SUBROUTINE RFLO_MoveGridElliptFra


!******************************************************************************
!
! Purpose: receive and distribute the deformations of surfaces
!          in block-wise manner.
!
! Description: surface smoothing is based on TFI
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

SUBROUTINE RFLO_MgElliptSurfacesFra( regions,someMoved )

  USE ModInterfaces, ONLY : RFLO_GetDeformation, RFLO_ArcLengthBounds, &
                            RFLO_EdgeDeformationStraight, &
                            RFLO_BoundaryDeformation
  USE RFLO_ModMoveGridFrame, ONLY : RFLO_MgFrameBroadCast, &
                                    RFLO_MgFrameCorrectNeighbors, &
                                    RFLO_MgFrameMoveCorners, &
                                    RFLO_MgFrameOrthoShift
  IMPLICIT NONE

! ... parameters
  LOGICAL :: someMoved

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iter

! ... local variables
  TYPE(t_grid), POINTER   :: grid, gridOld
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgElliptSurfacesFra',&
       'RFLO_ModMoveGridElliptFra.F90' )

! obtain external deformations ------------------------------------------------

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

      grid%xyzOld(:,:) = grid%xyz(:,:)

    ENDIF   ! region on this processor and active, grid moving
  ENDDO     ! iReg

! broadcast and compute block corners deformation -----------------------------

  iter = 1
  CALL RFLO_MgFrameBroadCast( regions,1,iter )
  CALL RFLO_MgFrameCorrectNeighbors( regions )
  CALL RFLO_MgFrameMoveCorners( regions )

  DO iter = 2,10
    CALL RFLO_MgFrameBroadCast( regions,1,iter )
    CALL RFLO_MgFrameMoveCorners( regions )
  ENDDO
  CALL RFLO_MgFrameBroadCast( regions,2,1 )             ! broadcast cBuff
  CALL RFLO_MgFrameOrthoShift( regions )                ! orth. to solid surf.

! compute edge and boundary deformations --------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE .AND. &            ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN          ! and moving

      grid      => regions(iReg)%levels(1)%grid
      gridOld   => regions(iReg)%levels(1)%gridOld

! --- calculate deformations at remaining edges

      grid%edgeMoved(:) = .FALSE.
      CALL RFLO_MgElliptEdgesFra( regions(iReg),grid%boundMoved, &
                               grid%allExternal,grid%edgeMoved,grid%arcLen12, &
                               grid%arcLen34,grid%arcLen56,gridOld%xyzOld, &
                               gridOld%xyz,grid%xyz )

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

END SUBROUTINE RFLO_MgElliptSurfacesFra


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
!        xyzOld     = initial grid.
!
! Output: edgeMoved = flag if discretization at an edge was changed
!         dNode     = updated deformations at edges.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************

SUBROUTINE RFLO_MgElliptEdgesFra( region,boundMoved,allExternal,edgeMoved, &
                                  arcLen12,arcLen34,arcLen56,xyzRef,xyzOld, &
                                  dNode )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes, &
                            RFLO_Tfint1d
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6), allExternal(6), edgeMoved(12)

  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:), xyzRef(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iEdge, ind

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, l1c, l2c
  INTEGER :: indBeg, indEnd, ijkN, ijkN1, ijkNBeg, ijkNEnd, iNOff, ijNOff
  INTEGER :: switch(12,11), interType, iEdgeGlo

  REAL(RFREAL) :: arcLen, ds, s, dN(3), dNBeg(3), dNEnd(3)
  LOGICAL :: interact

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_MgElliptEdgesFra',&
       'RFLO_ModMoveGridElliptFra.F90')

! get dimensions --------------------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! set edge switch -------------------------------------------------------------
! switch(:,1)  = begins at boundary
! switch(:,2)  = ends on boundary
! switch(:,3)  = right boundary
! switch(:,4)  = left boundary
! switch(:,5)  = direction (from-to boundary)
! switch(:,6)  = start index
! switch(:,7)  = end index
! switch(:,8)  = constant index in 1st direction
! switch(:,9)  = constant index in 2nd direction
! switch(:,10) = start corner number
! switch(:,11) = end corner number

  switch( 1,:) = (/5, 6, 1, 3, 56, kpnbeg, kpnend, ipnbeg, jpnbeg,  1,  2/)
  switch( 2,:) = (/3, 4, 1, 6, 34, jpnbeg, jpnend, kpnend, ipnbeg,  2,  3/)
  switch( 3,:) = (/5, 6, 1, 4, 56, kpnbeg, kpnend, ipnbeg, jpnend,  4,  3/)
  switch( 4,:) = (/3, 4, 1, 5, 34, jpnbeg, jpnend, kpnbeg, ipnbeg,  1,  4/)
  switch( 5,:) = (/5, 6, 2, 3, 56, kpnbeg, kpnend, ipnend, jpnbeg,  5,  6/)
  switch( 6,:) = (/3, 4, 2, 6, 34, jpnbeg, jpnend, kpnend, ipnend,  6,  7/)
  switch( 7,:) = (/5, 6, 2, 4, 56, kpnbeg, kpnend, ipnend, jpnend,  8,  7/)
  switch( 8,:) = (/3, 4, 2, 5, 34, jpnbeg, jpnend, kpnbeg, ipnend,  5,  8/)
  switch( 9,:) = (/1, 2, 3, 5, 12, ipnbeg, ipnend, jpnbeg, kpnbeg,  1,  5/)
  switch(10,:) = (/1, 2, 3, 6, 12, ipnbeg, ipnend, jpnbeg, kpnend,  2,  6/)
  switch(11,:) = (/1, 2, 4, 5, 12, ipnbeg, ipnend, jpnend, kpnbeg,  4,  8/)
  switch(12,:) = (/1, 2, 4, 6, 12, ipnbeg, ipnend, jpnend, kpnend,  3,  7/)

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
    IF (.NOT.edgeMoved(iEdge)) THEN

      edgeMoved(iEdge) = .true.

      ds     = 0._RFREAL
      indBeg = switch(iEdge,6)
      indEnd = switch(iEdge,7)
      l1c    = switch(iEdge,8)
      l2c    = switch(iEdge,9)

      iEdgeGlo = iEdge
      IF (iEdge==11) iEdgeGlo=12
      IF (iEdge==12) iEdgeGlo=11
      interact  = region%levels(iLev)%edgeCells(iEdgeGlo)%interact
      interType = region%levels(iLev)%edgeCells(iEdgeGlo)%interType

      IF (((region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,10))==1 .OR. &
            region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,11))==1) .AND. &
           ((interact .EQV. .true.) .AND. (interType==EDGE_INTERACT_FULL))) &
                                    .OR. &
          region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,10))==2 .OR. &
          region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,11))==2) THEN

!      IF (region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,10))==1 .OR. &
!          region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,11))==1 .OR. &
!          region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,10))==2 .OR. &
!          region%levels(iLev)%grid%nghbor(3,1,switch(iEdge,11))==2) THEN

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
      ENDIF   ! nghbor
    ENDIF     ! edgeMoved
  ENDDO       ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_MgElliptEdgesFra


!******************************************************************************
!
! Purpose: calculate node displacements on those boundaries whose edges
!          have moved but which were not marked as moving (finest grid only).
!
! Description: none.
!
! Input: region     = grid dimensions
!        boundMoved = flag for boundaries of a region which have moved
!        edgeMoved  = flag for edges whose nodes have moved
!        arcLen12   = arclength between i=const. boundaries for each j, k
!        arcLen34   = arclength between j=const. boundaries for each k, i
!        arcLen56   = arclength between k=const. boundaries for each i, j
!        xyzOld     = initial grid.
!
! Output: dNode = updated deformations at boundaries.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************

SUBROUTINE RFLO_MgElliptBoundaries( region,boundMoved,edgeMoved, &
                                    arcLen12,arcLen34,arcLen56,  &
                                    xyzOld,dNode )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes, &
                            RFLO_Tfint2d
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6), edgeMoved(12)

  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iBound, l1, l2

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: l1b, l1e, l2b, l2e, lc, ijkN, ijkE(4), ijkEm(4), iNOff, ijNOff
  INTEGER :: switch(6,9)

  LOGICAL :: sum12

  REAL(RFREAL) :: arcLen(4), ds(4), s(4)
  REAL(RFREAL) :: corner(3,8), e1(3), e2(3), e3(3), e4(3), &
                  p1(3), p2(3), p3(3), p4(3), dN(3)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_MgElliptBoundaries',&
       'RFLO_ModMoveGridElliptFra.F90' )

! get dimensions --------------------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! set boundary switch ---------------------------------------------------------
! switch(:,1-4) = numbers of the 4 edges of a boundary
! switch(:,5-6) = first/last index in l1-direction
! switch(:,7-8) = first/last index in l2-direction
! switch(:,  9) = constant index

  switch(1,:) = (/ 1,  2,  3,  4, jpnbeg, jpnend, kpnbeg, kpnend, ipnbeg/)
  switch(2,:) = (/ 5,  6,  7,  8, jpnbeg, jpnend, kpnbeg, kpnend, ipnend/)
  switch(3,:) = (/ 1,  5,  9, 10, kpnbeg, kpnend, ipnbeg, ipnend, jpnbeg/)
  switch(4,:) = (/ 3,  7, 11, 12, kpnbeg, kpnend, ipnbeg, ipnend, jpnend/)
  switch(5,:) = (/ 4,  8,  9, 11, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg/)
  switch(6,:) = (/ 2,  6, 10, 12, ipnbeg, ipnend, jpnbeg, jpnend, kpnend/)

! store displacements at corners ----------------------------------------------

  corner(:,1) = dNode(:,IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff))
  corner(:,2) = dNode(:,IndIJK(ipnbeg,jpnbeg,kpnend,iNOff,ijNOff))
  corner(:,3) = dNode(:,IndIJK(ipnbeg,jpnend,kpnend,iNOff,ijNOff))
  corner(:,4) = dNode(:,IndIJK(ipnbeg,jpnend,kpnbeg,iNOff,ijNOff))
  corner(:,5) = dNode(:,IndIJK(ipnend,jpnbeg,kpnbeg,iNOff,ijNOff))
  corner(:,6) = dNode(:,IndIJK(ipnend,jpnbeg,kpnend,iNOff,ijNOff))
  corner(:,7) = dNode(:,IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff))
  corner(:,8) = dNode(:,IndIJK(ipnend,jpnend,kpnbeg,iNOff,ijNOff))

! move nodes on boundaries with active edges ----------------------------------

  DO iBound=1,6
    IF ((.NOT.boundMoved(iBound)) .AND. &
        (edgeMoved(switch(iBound,1)) .OR. edgeMoved(switch(iBound,2)) .OR. &
         edgeMoved(switch(iBound,3)) .OR. edgeMoved(switch(iBound,4)))) THEN

!    IF ((edgeMoved(switch(iBound,1)) .OR. edgeMoved(switch(iBound,2)) .OR. &
!         edgeMoved(switch(iBound,3)) .OR. edgeMoved(switch(iBound,4)))) THEN

      l1b = switch(iBound,5)
      l1e = switch(iBound,6)
      l2b = switch(iBound,7)
      l2e = switch(iBound,8)
      lc  = switch(iBound,9)

      IF (iBound == 1) THEN
        p1(:) = corner(:,1)
        p2(:) = corner(:,4)
        p3(:) = corner(:,3)
        p4(:) = corner(:,2)
      ELSE IF (iBound == 2) THEN
        p1(:) = corner(:,5)
        p2(:) = corner(:,8)
        p3(:) = corner(:,7)
        p4(:) = corner(:,6)
      ELSE IF (iBound == 3) THEN
        p1(:) = corner(:,1)
        p2(:) = corner(:,2)
        p3(:) = corner(:,6)
        p4(:) = corner(:,5)
      ELSE IF (iBound == 4) THEN
        p1(:) = corner(:,4)
        p2(:) = corner(:,3)
        p3(:) = corner(:,7)
        p4(:) = corner(:,8)
      ELSE IF (iBound == 5) THEN
        p1(:) = corner(:,1)
        p2(:) = corner(:,5)
        p3(:) = corner(:,8)
        p4(:) = corner(:,4)
      ELSE IF (iBound == 6) THEN
        p1(:) = corner(:,2)
        p2(:) = corner(:,6)
        p3(:) = corner(:,7)
        p4(:) = corner(:,3)
      ENDIF

      ds(1:2) = 0._RFREAL
      DO l2=l2b+1,l2e-1

        sum12   = .true.
        ds(3:4) = 0._RFREAL
        DO l1=l1b+1,l1e-1
          IF (iBound==1 .OR. iBound==2) THEN
            ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
            ijkE(1)   = IndIJK(lc,jpnbeg,l2    ,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(lc,jpnbeg,l2-1  ,iNOff,ijNOff)
            ijkE(2)   = IndIJK(lc,jpnend,l2    ,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(lc,jpnend,l2-1  ,iNOff,ijNOff)
            ijkE(3)   = IndIJK(lc,l1    ,kpnbeg,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(lc,l1-1  ,kpnbeg,iNOff,ijNOff)
            ijkE(4)   = IndIJK(lc,l1    ,kpnend,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(lc,l1-1  ,kpnend,iNOff,ijNOff)
            arcLen(1) = arcLen56(lc,jpnbeg)
            arcLen(2) = arcLen56(lc,jpnend)
            arcLen(3) = arcLen34(kpnbeg,lc)
            arcLen(4) = arcLen34(kpnend,lc)
          ELSE IF (iBound==3 .OR. iBound==4) THEN
            ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
            ijkE(1)   = IndIJK(l2    ,lc,kpnbeg,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(l2-1  ,lc,kpnbeg,iNOff,ijNOff)
            ijkE(2)   = IndIJK(l2    ,lc,kpnend,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(l2-1  ,lc,kpnend,iNOff,ijNOff)
            ijkE(3)   = IndIJK(ipnbeg,lc,l1    ,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(ipnbeg,lc,l1-1  ,iNOff,ijNOff)
            ijkE(4)   = IndIJK(ipnend,lc,l1    ,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(ipnend,lc,l1-1  ,iNOff,ijNOff)
            arcLen(1) = arclen12(lc,kpnbeg)
            arcLen(2) = arclen12(lc,kpnend)
            arcLen(3) = arcLen56(ipnbeg,lc)
            arcLen(4) = arcLen56(ipnend,lc)
          ELSE IF (iBound==5 .OR. iBound==6) THEN
            ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
            ijkE(1)   = IndIJK(ipnbeg,l2    ,lc,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(ipnbeg,l2-1  ,lc,iNOff,ijNOff)
            ijkE(2)   = IndIJK(ipnend,l2    ,lc,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(ipnend,l2-1  ,lc,iNOff,ijNOff)
            ijkE(3)   = IndIJK(l1    ,jpnbeg,lc,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(l1-1  ,jpnbeg,lc,iNOff,ijNOff)
            ijkE(4)   = IndIJK(l1    ,jpnend,lc,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(l1-1  ,jpnend,lc,iNOff,ijNOff)
            arcLen(1) = arcLen34(lc,ipnbeg)
            arcLen(2) = arcLen34(lc,ipnend)
            arcLen(3) = arcLen12(jpnbeg,lc)
            arcLen(4) = arcLen12(jpnend,lc)
          ENDIF
          IF (sum12) THEN
            ds(1) = ds(1) + &
                    SQRT((xyzOld(XCOORD,ijkE(1))-xyzOld(XCOORD,ijkEm(1)))**2 + &
                         (xyzOld(YCOORD,ijkE(1))-xyzOld(YCOORD,ijkEm(1)))**2 + &
                         (xyzOld(ZCOORD,ijkE(1))-xyzOld(ZCOORD,ijkEm(1)))**2)
            ds(2) = ds(2) + &
                    SQRT((xyzOld(XCOORD,ijkE(2))-xyzOld(XCOORD,ijkEm(2)))**2 + &
                         (xyzOld(YCOORD,ijkE(2))-xyzOld(YCOORD,ijkEm(2)))**2 + &
                         (xyzOld(ZCOORD,ijkE(2))-xyzOld(ZCOORD,ijkEm(2)))**2)
            sum12 = .false.
          ENDIF
          ds(3) = ds(3) + &
                  SQRT((xyzOld(XCOORD,ijkE(3))-xyzOld(XCOORD,ijkEm(3)))**2 + &
                       (xyzOld(YCOORD,ijkE(3))-xyzOld(YCOORD,ijkEm(3)))**2 + &
                       (xyzOld(ZCOORD,ijkE(3))-xyzOld(ZCOORD,ijkEm(3)))**2)
          ds(4) = ds(4) + &
                  SQRT((xyzOld(XCOORD,ijkE(4))-xyzOld(XCOORD,ijkEm(4)))**2 + &
                       (xyzOld(YCOORD,ijkE(4))-xyzOld(YCOORD,ijkEm(4)))**2 + &
                       (xyzOld(ZCOORD,ijkE(4))-xyzOld(ZCOORD,ijkEm(4)))**2)
          s(:)  = ds(:)/arcLen(:)
          e1(:) = dNode(:,ijkE(1))
          e2(:) = dNode(:,ijkE(2))
          e3(:) = dNode(:,ijkE(3))
          e4(:) = dNode(:,ijkE(4))
          CALL RFLO_Tfint2d( s(1),s(2),s(3),s(4),e1,e2,e3,e4,p1,p2,p3,p4,dN )
          dNode(:,ijkN) = dN(:)
        ENDDO  ! l1
      ENDDO    ! l2

    ENDIF      ! edgeMoved
  ENDDO        ! iBound

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_MgElliptBoundaries


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

SUBROUTINE RFLO_MgElliptInterfacesFra( regions )

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

  CALL RegisterFunction( global,'RFLO_MgElliptInterfacesFra',&
       'RFLO_ModMoveGridElliptFra.F90' )

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
              CALL RFLO_MgElliptEdgesFra( regions(iReg),grid%boundMoved, &
                                       grid%allExternal,grid%edgeMoved, &
                                       grid%arcLen12,grid%arcLen34, &
                                       grid%arcLen56,gridOld%xyzOld, &
                                       gridOld%xyz,grid%xyz )
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
              CALL RFLO_MgElliptEdgesFra( regions(iReg),grid%boundMoved, &
                                       grid%allExternal,grid%edgeMoved, &
                                       grid%arcLen12,grid%arcLen34, &
                                       grid%arcLen56,gridOld%xyzOld, &
                                       gridOld%xyz,grid%xyz )
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

END SUBROUTINE RFLO_MgElliptInterfacesFra


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModMoveGridElliptFra

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModMoveGridElliptFra.F90,v $
! Revision 1.16  2009/08/27 14:04:50  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.15  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2006/03/18 13:26:21  wasistho
! added orthDir and orthWghtX,Y,Z
!
! Revision 1.12  2006/03/18 11:48:51  wasistho
! cosmetics
!
! Revision 1.11  2006/03/18 11:04:10  wasistho
! screen printed global skewness and minvol
!
! Revision 1.10  2006/03/18 09:24:47  wasistho
! modified mgElliptEdgesFra
!
! Revision 1.9  2006/03/18 08:17:48  wasistho
! added mgElliptSurfacesFra
!
! Revision 1.8  2006/03/16 08:27:09  wasistho
! added interfacesFra and edgesFra
!
! Revision 1.7  2006/03/15 06:38:19  wasistho
! added region and global skewness
!
! Revision 1.6  2006/03/14 04:41:13  wasistho
! added RFLO_EdgeDeformationStraight
!
! Revision 1.5  2006/03/09 19:14:25  wasistho
! prepared for global elliptic grid motion
!
! Revision 1.4  2006/03/05 19:08:31  wasistho
! set computational space coordinates from initial grid
!
! Revision 1.3  2006/03/04 04:35:41  wasistho
! call elliptGridSmooRegion if blockShape is normal
!
! Revision 1.2  2006/03/03 06:07:05  wasistho
! enabled global elliptic PDE grid motion
!
! Revision 1.1  2006/02/11 03:42:14  wasistho
! added ModMoveGridElliptGlo/Fra
!
!
!
! ******************************************************************************











