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
! Purpose: redistribute grid nodes according to the movement of the
!          boundaries. This function smoothes the grid globally by
!          solving the Laplace equation for node coordinates.
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
!
! $Id: RFLO_MoveGridGlobal.F90,v 1.9 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_MoveGridGlobal( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY       : t_grid
  USE ModInterfaces, ONLY : RFLO_MoveGridSurfaces, RFLO_MoveGridInterfaces, &
        RFLO_ChangeInteriorGrid, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_CalcControlVolumes, &
        RFLO_CalcCellCentroids, RFLO_C2fAvgCoeffs, RFLO_C2eAvgCoeffs, &
        RFLO_CheckMetrics, RFLO_CalcGridSpeeds, RFLO_BoundaryDeformation, &
        RFLO_GridRemesh
  USE RFLO_ModLaplaceSmoothing, ONLY : RFLO_LaplaceGridSmoo
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iter, iPatch, ijk

! ... local variables
  LOGICAL :: someMoved, someRemesh

  INTEGER :: bcType, iRemesh, jRemesh, nRemesh

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

  CALL RegisterFunction( global,'RFLO_MoveGridGlobal',&
  'RFLO_MoveGridGlobal.F90' )

#ifdef GENX
! update geometry buffers -----------------------------------------------------

  dAlpha = global%dtMin/global%dTimeSystem
  CALL COM_call_function( global%genxHandleGm,1,dAlpha )
#endif

! receive and distribute deformations for each region -------------------------

  CALL RFLO_MoveGridSurfaces( regions,someMoved )

! fix interfaces between regions ----------------------------------------------

  IF (someMoved) THEN
    CALL RFLO_MoveGridInterfaces( regions )
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
    DO iter=1,global%moveGridNiter
      CALL RFLO_LaplaceGridSmoo( regions,resid )
    ENDDO

    IF (global%verbLevel /= VERBOSE_NONE) THEN
#ifdef MPI
      CALL MPI_Reduce( resid,globalResid,1,MPI_RFREAL,MPI_SUM, &
                       MASTERPROC,global%mpiComm,global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
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
        IF ((bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE)) THEN
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

  someRemesh = .FALSE.
  iRemesh = 0
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
!      IF (regions(iReg)%levels(1)%grid%remesh==1) THEN
!        CALL RFLO_GridRemesh( regions(iReg) )            ! grid remeshing
!        iRemesh=1
!      ENDIF
      CALL RFLO_CalcGridSpeeds( regions(iReg) )          ! grid speeds
    ENDIF   ! region on this processor and active, grid moving
  ENDDO     ! iReg

#ifdef MPI
  CALL MPI_ALLREDUCE( iRemesh, nRemesh, 1, MPI_INTEGER, MPI_SUM, &
                      global%mpiComm, global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  IF (nRemesh > 0) someRemesh = .TRUE.
#endif

  IF (someRemesh) THEN
    CALL RFLO_ExchangeGeometry( regions )    ! exchange geometry
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
          regions(iReg)%active==ACTIVE .AND. &             ! on my processor
          iRemesh==1) THEN                                  ! and remeshing
        CALL RFLO_CalcFaceVectors( regions(iReg) )         ! faces
        CALL RFLO_CalcControlVolumes( regions(iReg) )      ! volumes
        CALL RFLO_CalcCellCentroids( regions(iReg) )       ! cell centroids
        IF (regions(iReg)%mixtInput%faceEdgeAvg==FE_AVG_LINEAR) &
          CALL RFLO_C2fAvgCoeffs( regions(iReg) )          ! cell2face averaging
        CALL RFLO_C2eAvgCoeffs( regions(iReg) )            ! cell2edge averaging
      ENDIF   ! region on this processor and active, grid moving
    ENDDO     ! iReg
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(A,1X,'Block Laplacian grid motion: ',I6,1PE13.4)

END SUBROUTINE RFLO_MoveGridGlobal

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_MoveGridGlobal.F90,v $
! Revision 1.9  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/03/05 19:13:57  wasistho
! set computational space coordinates from initial grid
!
! Revision 1.6  2006/03/03 00:48:28  wasistho
! renamed global grid motion to block laplacian motion
!
! Revision 1.5  2005/11/03 02:40:28  wasistho
! activate boundaryDeformation and chngInteriorGrid
!
! Revision 1.4  2005/10/27 05:51:48  wasistho
! added USE RFLO_ModLaplaceSmoothing
!
! Revision 1.3  2005/05/28 05:52:32  wasistho
! refrain from grid remeshing
!
! Revision 1.2  2005/05/27 01:51:17  wasistho
! added remeshing
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.10  2004/09/02 02:35:26  wasistho
! added face-edge averaging input-option parameter in Rocflo
!
! Revision 1.9  2004/08/03 22:46:18  wasistho
! added RFLO_c2eAvgCoeffs
!
! Revision 1.8  2004/07/30 17:26:12  wasistho
! provide cell2face averaging coefficients
!
! Revision 1.7  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/08/27 23:58:10  jblazek
! Removed 2nd interface to RFLO_ChangeInteriorGrid.
!
! Revision 1.2  2003/08/25 21:51:24  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.1  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
!******************************************************************************







